# Unit tests for src/geometry/domains/multiplyconnecteddomains.jl: the MultiplyConnectedDomain type, its constructor and its
#   domain queries (genus, hole_domain_ids, is_hole_curve, component_of_domain_id, validate_holes).
# Reference objects (fixtures from holes_fixtures.jl): the square [0,4]^2 as outer region with a circular hole (domain id 2) and a
# square hole [2,3]^2 (domain id 3); the composite outer region of StadiumBilliard; the quarter-square outer region of the D2
# rectangle. The curve `orientation` field is never read to decide holes: a hole stored with orientation -1 is still a hole, and a
# concave outer piece (inverted circle, orientation -1, region outside the circle) is not. No eigen solves. Helpers use the prefix `mc_`.
#
# Symmetry group G' of a billiard with holes (B4; _symmetry_group behind get_symmetry_group, its forwarding get_symmetries,
#   _first_incompatibility, _log_symmetry_group): G' is the subset of the billiard registry (original elements and order, same sym_ids)
#   that every hole is compatible with, so it is the intersection of compatible_symmetries over the holes. The fixtures sit in the
#   quarter-square billiard with the D2 registry (Y, XY, X = sym_ids 1, 2, 3): a D4 disc at the origin and four discs placed by D2 keep
#   all of D2, the D2 disc about (0.5, 0) and the half disc about (0.3, 0) keep {X}, a DomainHole keeps nothing. Group algebra only,
#   the holes of the intersection cases need not be disjoint because hf_HoleBilliard does not run validate_holes.

using StaticArrays: SVector
using BilliardGeometry: _first_incompatibility, _log_symmetry_group, _symmetry_group

mc_outer() = hf_rect_domain(0.0, 0.0, 4.0, 4.0, 1)
mc_circle_hole(id = 2) = DomainHole(hf_circle_domain(0.3, 1.0, 1.0, id))
mc_square_hole(id = 3) = DomainHole(hf_rect_domain(2.0, 2.0, 3.0, 3.0, id))

@testset "MultiplyConnectedDomain fields" begin
    outer = mc_outer()
    hole_a = mc_circle_hole()
    hole_b = mc_square_hole()
    dom = MultiplyConnectedDomain(outer, [hole_a, hole_b], 4)

    @test dom isa MultiplyConnectedDomain{Float64}
    @test dom isa AbsMultiplyConnectedDomain
    @test !(dom isa AbsSimpleDomain)
    @test fieldnames(MultiplyConnectedDomain) == (:subdomains, :holes, :id)
    @test dom.id == 4
    @test dom.subdomains isa Vector{AbsSimpleDomain}
    @test dom.holes isa Vector{AbsDomainHole}
    @test length(dom.subdomains) == 1
    @test dom.subdomains[1] === outer
    # the holes are stored as given, in order
    @test length(dom.holes) == 2
    @test dom.holes[1] === hole_a
    @test dom.holes[2] === hole_b
    # the old boundary / corners / genus fields are gone
    @test !hasfield(MultiplyConnectedDomain, :boundary)
    @test !hasfield(MultiplyConnectedDomain, :corners)
    @test !hasfield(MultiplyConnectedDomain, :genus)
end

@testset "MultiplyConnectedDomain outer region forms" begin
    # a single simple domain, a vector of simple domains and a composite domain are flattened into `subdomains`
    simple = MultiplyConnectedDomain(mc_outer(), [mc_circle_hole()], 1)
    @test length(simple.subdomains) == 1

    stadium = StadiumBilliard(1.0).fundamental_domain
    @test stadium isa CompositeDomain
    from_composite = MultiplyConnectedDomain(stadium, [mc_circle_hole(7)], 1)
    @test length(from_composite.subdomains) == 2
    @test [sd.id for sd in from_composite.subdomains] == [1, 2]
    from_vector = MultiplyConnectedDomain(stadium.subdomains, [mc_circle_hole(7)], 1)
    @test length(from_vector.subdomains) == 2
    @test from_vector.subdomains == from_composite.subdomains
    @test from_vector.subdomains isa Vector{AbsSimpleDomain}

    # a vector of the concrete outer subdomain type works too
    @test length(MultiplyConnectedDomain([mc_outer()], [mc_circle_hole()], 1).subdomains) == 1
end

@testset "MultiplyConnectedDomain hole forms" begin
    # a bare domain in `holes` is wrapped in a DomainHole, a DomainHole is kept as is
    bare = hf_circle_domain(0.3, 1.0, 1.0, 2)
    kept = mc_square_hole()
    dom = MultiplyConnectedDomain(mc_outer(), [bare, kept], 1)
    @test dom.holes[1] isa DomainHole
    @test dom.holes[1].domain === bare
    @test dom.holes[2] === kept

    # no holes at all is a valid (simply connected) multiply connected domain object
    none = MultiplyConnectedDomain(mc_outer(), [], 1)
    @test isempty(none.holes)
    @test none.holes isa Vector{AbsDomainHole}
    @test genus(none) == 0

    # a composite hole keeps its two loops in one hole object
    dom2 = MultiplyConnectedDomain(mc_outer(), [hf_two_loop_hole(5)], 1)
    @test length(dom2.holes) == 1
end

@testset "MultiplyConnectedDomain element type" begin
    outer32 = hf_rect_domain(0f0, 0f0, 2f0, 2f0, 1)
    hole32 = hf_circle_domain(0.2f0, 1f0, 1f0, 2)
    dom32 = MultiplyConnectedDomain(outer32, [hole32], 1)
    @test dom32 isa MultiplyConnectedDomain{Float32}
    @test genus(dom32) == 1
    # the id is stored as Int64
    @test dom32.id isa Int64
    @test MultiplyConnectedDomain(mc_outer(), [mc_circle_hole()], Int32(3)).id === Int64(3)
end

@testset "MultiplyConnectedDomain argument errors" begin
    outer = mc_outer()

    # an empty outer region
    @test_throws ArgumentError MultiplyConnectedDomain(AbsSimpleDomain[], [mc_circle_hole()], 1)
    @test_throws "at least one outer subdomain" MultiplyConnectedDomain(AbsSimpleDomain[], [mc_circle_hole()], 1)

    # the outer region may not itself be a multiply connected domain
    inner_mcd = MultiplyConnectedDomain(outer, [mc_circle_hole()], 1)
    @test_throws ArgumentError MultiplyConnectedDomain(inner_mcd, [mc_square_hole()], 1)

    # a hole may not contain a multiply connected domain, whether wrapped or bare
    nested = MultiplyConnectedDomain(hf_rect_domain(2.0, 2.0, 3.0, 3.0, 3), [], 9)
    @test_throws "MultiplyConnectedDomain" MultiplyConnectedDomain(outer, [DomainHole(nested)], 1)
    @test_throws ArgumentError MultiplyConnectedDomain(outer, [nested], 1)

    # curve domain_ids may not overlap between the outer region and a hole, or between two holes
    @test_throws "overlap" MultiplyConnectedDomain(outer, [mc_circle_hole(1)], 1)
    @test_throws ArgumentError MultiplyConnectedDomain(outer, [mc_circle_hole(2), mc_square_hole(2)], 1)
    @test_throws ArgumentError MultiplyConnectedDomain(outer, [mc_circle_hole(2), hf_two_loop_hole(1)], 1)

    # the exterior boundary of the outer region must be closed: three of the four edges of the square
    open_outer = SimpleDomain{Float64}(outer.boundary[1:3], outer.corners[1:3], 1)
    @test_throws ArgumentError MultiplyConnectedDomain(open_outer, [mc_circle_hole()], 1)
    @test_throws "outer boundary" MultiplyConnectedDomain(open_outer, [mc_circle_hole()], 1)

    # every hole loop must be closed
    sq = hf_rect_domain(2.0, 2.0, 3.0, 3.0, 3)
    open_hole = DomainHole(SimpleDomain{Float64}(sq.boundary[1:3], sq.corners[1:3], 3))
    @test_throws ArgumentError MultiplyConnectedDomain(outer, [open_hole], 1)
    @test_throws "hole 1" MultiplyConnectedDomain(outer, [open_hole], 1)
    @test_throws "hole 2" MultiplyConnectedDomain(outer, [mc_circle_hole(), open_hole], 1)

    # the id is typed as an integer
    @test_throws MethodError MultiplyConnectedDomain(outer, [mc_circle_hole()], 1.5)
end

@testset "genus" begin
    outer = mc_outer()
    # genus is the total number of hole loops, not of hole objects or hole curves
    @test genus(MultiplyConnectedDomain(outer, [], 1)) == 0
    @test genus(MultiplyConnectedDomain(outer, [mc_circle_hole()], 1)) == 1
    @test genus(MultiplyConnectedDomain(outer, [mc_circle_hole(), mc_square_hole()], 1)) == 2
    # one composite hole of two loops has genus 2
    @test genus(MultiplyConnectedDomain(outer, [hf_two_loop_hole(5)], 1)) == 2
    @test genus(MultiplyConnectedDomain(outer, [mc_circle_hole(), hf_two_loop_hole(5)], 1)) == 3
    # a square hole has four curves but one loop
    @test genus(MultiplyConnectedDomain(outer, [mc_square_hole()], 1)) == 1
    # billiards forward to their fundamental domain
    @test genus(hf_concentric_annulus()) == 1
    @test genus(hf_rect_with_hole()) == 1
    @test genus(hf_d2_with_hole()) == 1
    # simply connected domains and billiards have genus 0
    @test genus(mc_outer()) == 0
    @test genus(RectangleBilliard(1.0, 0.6)) == 0
end

@testset "hole_domain_ids" begin
    outer = mc_outer()
    dom = MultiplyConnectedDomain(outer, [mc_circle_hole(2), mc_square_hole(3)], 1)
    @test hole_domain_ids(dom) == [2, 3]
    # all ids of a composite hole are listed, sorted
    @test hole_domain_ids(MultiplyConnectedDomain(outer, [mc_circle_hole(9), hf_two_loop_hole(5)], 1)) == [5, 6, 9]
    @test hole_domain_ids(MultiplyConnectedDomain(outer, [], 1)) == Int[]
    # billiards forward to their fundamental domain
    @test hole_domain_ids(hf_shifted_annulus()) == [2]
    @test hole_domain_ids(hf_d2_with_hole()) == [5]
    # a symmetric hole lists the ids of all its pieces (id_stride 100): two half discs for the annulus, eight pieces for Sinai
    @test hole_domain_ids(hf_concentric_annulus()) == [2, 102]
    @test hole_domain_ids(SinaiBilliard(1.0, 0.3)) == [2, 102, 202, 302, 402, 502, 602, 702]
    # a domain without holes has none
    @test hole_domain_ids(outer) == Int[]
    @test hole_domain_ids(RectangleBilliard(1.0, 0.6)) == Int[]
end

@testset "is_hole_curve" begin
    outer = mc_outer()
    dom = MultiplyConnectedDomain(outer, [mc_circle_hole(2), mc_square_hole(3)], 1)
    @test all(c -> !is_hole_curve(dom, c), outer.boundary)
    @test is_hole_curve(dom, dom.holes[1].domain.boundary[1])
    @test all(c -> is_hole_curve(dom, c), dom.holes[2].domain.boundary)
    # decided by the domain_id of the curve, the orientation field is not read
    inverted = CircleSegment(0.3, 2pi, 0.0, [1.0, 1.0]; orientation = -1, domain_id = 2)
    @test is_hole_curve(dom, inverted)
    flipped_outer = LineSegment([0.0, 0.0], [4.0, 0.0]; orientation = -1, domain_id = 1)
    @test !is_hole_curve(dom, flipped_outer)

    # billiards forward to their fundamental domain
    ann = hf_shifted_annulus()
    outer_curve, hole_curve = get_all_curves(ann)
    @test !is_hole_curve(ann, outer_curve)
    @test is_hole_curve(ann, hole_curve)
    # symmetric annulus: the two outer curves (arc, wall) are not hole curves, the three curves of each half-disc piece are
    sym_ann = hf_concentric_annulus()
    @test [is_hole_curve(sym_ann, c) for c in get_all_curves(sym_ann)] == [false, false, true, true, true, true, true, true]
    # Sinai: the three outer curves of the fundamental triangle are not hole curves, the 8 pieces of 3 curves of the obstacle are
    sinai = SinaiBilliard(1.0, 0.3)
    @test [is_hole_curve(sinai, c) for c in get_all_curves(sinai)] == [fill(false, 3); fill(true, 24)]
    # unreduced (shifted) Sinai: the four sides are not hole curves, the full circle is
    sinai_shift = SinaiBilliard(1.0, 0.3; center = [0.5, -0.2])
    @test [is_hole_curve(sinai_shift, c) for c in get_all_curves(sinai_shift)] == [false, false, false, false, true]
    @test !is_hole_curve(mc_outer(), outer.boundary[1])
end

@testset "component_of_domain_id" begin
    outer = mc_outer()
    dom = MultiplyConnectedDomain(outer, [mc_circle_hole(2), hf_two_loop_hole(5)], 1)
    # 0 for the outer region, i for the i-th hole object (both loops of the composite hole share group 2)
    @test component_of_domain_id(dom) == Dict(1 => 0, 2 => 1, 5 => 2, 6 => 2)
    @test component_of_domain_id(hf_shifted_annulus()) == Dict(1 => 0, 2 => 1)
    @test component_of_domain_id(hf_d2_with_hole()) == Dict(1 => 0, 5 => 1)
    # all pieces of a symmetric hole belong to the one hole component
    @test component_of_domain_id(hf_concentric_annulus()) == Dict(1 => 0, 2 => 1, 102 => 1)
    @test component_of_domain_id(SinaiBilliard(1.0, 0.3)) == Dict(1 => 0, [2 + 100 * k => 1 for k in 0:7]...)
    # non-multiply-connected domains and billiards: everything is the outer component
    @test component_of_domain_id(outer) == Dict(1 => 0)
    @test component_of_domain_id(StadiumBilliard(1.0)) == Dict(1 => 0, 2 => 0)
end

@testset "concave outer piece is not a hole" begin
    # an inverted circle (orientation -1) bounding the outer region is a concave piece, not a hole
    bite = hf_bitten_square()
    @test count(c -> c.orientation == -1, bite.boundary) == 1
    dom = MultiplyConnectedDomain(bite, [], 1)
    @test genus(dom) == 0
    @test hole_domain_ids(dom) == Int[]
    @test all(c -> !is_hole_curve(dom, c), bite.boundary)
    @test boundary_component_is_hole(dom) == [false]
    @test length(boundary_components(dom)) == 1
    # the area follows from the Green integral over the closed loop: 1 - pi/16
    @test area(bite) ≈ 1 - pi / 16 atol = 1e-12
    b = hf_HoleBilliard{Float64}(dom, register_symmetries())
    @test area(b) ≈ 1 - pi / 16 atol = 1e-12
    # membership follows the domain function of the inverted arc: inside means outside the disk
    @test is_inside(dom, SVector(0.8, 0.8)) === true
    @test is_inside(dom, SVector(0.1, 0.1)) === false
    @test is_inside(dom, SVector(1.2, 0.5)) === false
end

@testset "validate_holes" begin
    outer_circle = hf_circle_domain(2.0, 0.0, 0.0, 1)
    # accepted: holes inside the outer region and disjoint from each other
    valid = MultiplyConnectedDomain(outer_circle, [hf_circle_domain(0.3, 1.0, 0.0, 2), hf_circle_domain(0.3, -1.0, 0.0, 3)], 1)
    @test validate_holes(valid, register_symmetries()) === nothing
    @test validate_holes(hf_concentric_annulus().fundamental_domain, register_symmetries()) === nothing
    @test validate_holes(hf_eccentric_annulus().fundamental_domain, register_symmetries()) === nothing
    @test validate_holes(hf_shifted_annulus().fundamental_domain, register_symmetries()) === nothing
    # the symmetric annulus and Sinai hole are valid with the registry of their billiard
    @test validate_holes(hf_concentric_annulus().fundamental_domain, hf_concentric_annulus().symmetries) === nothing
    sinai = SinaiBilliard(1.0, 0.3)
    @test validate_holes(sinai.fundamental_domain, sinai.symmetries) === nothing
    @test validate_holes(hf_rect_with_hole().fundamental_domain, register_symmetries()) === nothing
    # D2 billiard: the hole is in the fundamental quarter
    @test validate_holes(hf_d2_with_hole().fundamental_domain, D2_symmetry()) === nothing

    # a hole outside the outer region is rejected
    outside = MultiplyConnectedDomain(outer_circle, [hf_circle_domain(0.3, 5.0, 0.0, 2)], 1)
    @test_throws ArgumentError validate_holes(outside, register_symmetries())
    @test_throws "hole 1 is not inside the outer region" validate_holes(outside, register_symmetries())
    # the offending hole is named
    second_outside = MultiplyConnectedDomain(outer_circle, [hf_circle_domain(0.3, 1.0, 0.0, 2), hf_circle_domain(0.3, 5.0, 0.0, 3)], 1)
    @test_throws "hole 2" validate_holes(second_outside, register_symmetries())

    # a hole inside another hole is rejected (the sample point of hole 2 lies in hole 1)
    nested = MultiplyConnectedDomain(outer_circle, [hf_circle_domain(1.0, 0.0, 0.0, 2), hf_circle_domain(0.2, 0.1, 0.0, 3)], 1)
    @test_throws ArgumentError validate_holes(nested, register_symmetries())
    @test_throws "hole 2 overlaps hole 1" validate_holes(nested, register_symmetries())

    # the outer region is reconstructed by the registry: a hole in another quadrant is inside it only for D2
    other_quadrant = MultiplyConnectedDomain(hf_quarter_square_domain(), [hf_circle_domain(0.1, -0.5, -0.5, 5)], 1)
    @test validate_holes(other_quadrant, D2_symmetry()) === nothing
    @test_throws ArgumentError validate_holes(other_quadrant, register_symmetries())
    # outside the full square it is rejected even with the registry
    far = MultiplyConnectedDomain(hf_quarter_square_domain(), [hf_circle_domain(0.1, -3.0, -3.0, 5)], 1)
    @test_throws ArgumentError validate_holes(far, D2_symmetry())
end

# Quarter-square billiard (D2 registry unless given) holding the given holes.
mc_hole_billiard(holes...; reg = D2_symmetry()) = hf_HoleBilliard{Float64}(MultiplyConnectedDomain(hf_quarter_square_domain(), collect(AbsDomainHole, holes), 1), reg)

# sym_ids of a registry in order.
mc_ids(reg) = [g.sym_id for g in reg]

@testset "get_symmetry_group (billiard with holes)" begin
    # one hole: the compatible subgroup, original elements, registry order and sym_ids
    g4 = get_symmetry_group(mc_hole_billiard(hf_d4_symmetric_hole()))
    @test g4 isa SymmetryRegistry
    @test mc_ids(g4) == [1, 2, 3]
    @test [typeof(g) for g in g4] == [YAxisReflection, XYAxisReflection, XAxisReflection]
    @test mc_ids(get_symmetry_group(mc_hole_billiard(hf_symmetrized_disc()))) == [1, 2, 3]
    @test mc_ids(get_symmetry_group(mc_hole_billiard(hf_d2_symmetric_hole()))) == [3]
    @test mc_ids(get_symmetry_group(mc_hole_billiard(hf_half_disc_hole(); reg = register_symmetries(XAxisReflection())))) == [1]
    @test mc_ids(get_symmetry_group(mc_hole_billiard(hf_half_disc_hole()))) == [3]
    # a DomainHole has no symmetry: G' is trivial
    @test length(get_symmetry_group(mc_hole_billiard(DomainHole(hf_circle_domain(0.2, 0.3, 0.4, 5))))) == 0
    @test length(get_symmetry_group(hf_d2_with_hole())) == 0
    # an empty registry stays empty
    @test length(get_symmetry_group(hf_rect_with_hole())) == 0
    @test length(get_symmetry_group(hf_shifted_annulus())) == 0
    # origin-centered production billiards are fully symmetric: G' = G (the registry elements, identity not counted)
    @test mc_ids(get_symmetry_group(hf_concentric_annulus())) == [1]
    @test [typeof(g) for g in get_symmetry_group(hf_concentric_annulus())] == [XAxisReflection]
    @test length(get_symmetry_group(SinaiBilliard(1.0, 0.3))) == 7
    @test mc_ids(get_symmetry_group(SinaiBilliard(1.0, 0.3))) == mc_ids(SinaiBilliard(1.0, 0.3).symmetries)
    @test length(get_symmetry_group(SquareWithinSquareBilliard(1.0; positions = [0.0], sizes = [0.4]))) == 7
    # a shifted center (Sinai) and the default two-obstacle configuration keep the empty registry
    @test length(get_symmetry_group(SinaiBilliard(1.0, 0.3; center = [0.5, -0.2]))) == 0
    @test length(get_symmetry_group(SquareWithinSquareBilliard(1.0))) == 0

    # several holes: intersection over the holes
    other_disc = DomainHole(hf_circle_domain(0.05, 0.6, 0.6, 9))
    @test mc_ids(get_symmetry_group(mc_hole_billiard(hf_d4_symmetric_hole(), hf_half_disc_hole()))) == [3]
    @test mc_ids(get_symmetry_group(mc_hole_billiard(hf_d4_symmetric_hole(), hf_symmetrized_disc()))) == [1, 2, 3]
    @test length(get_symmetry_group(mc_hole_billiard(hf_d4_symmetric_hole(), other_disc))) == 0
    # G' is a subset of the registry, whatever the holes
    for b in (mc_hole_billiard(hf_d2_symmetric_hole()), mc_hole_billiard(hf_d4_symmetric_hole(), hf_half_disc_hole()))
        @test issubset(mc_ids(get_symmetry_group(b)), mc_ids(b.symmetries))
    end

    # D4 registry: the D4 disc keeps all seven elements, the D2 disc only X
    @test mc_ids(get_symmetry_group(mc_hole_billiard(hf_d4_symmetric_hole(); reg = D4_symmetry()))) == collect(1:7)
    @test mc_ids(get_symmetry_group(mc_hole_billiard(hf_d2_symmetric_hole(); reg = D4_symmetry()))) == [7]

    # billiards without holes keep returning their registry
    circ = CircleBilliard(1.0)
    @test get_symmetry_group(circ) === circ.symmetries
end

@testset "get_symmetries (billiard with holes)" begin
    # forwards to get_symmetry_group: generators of G', not of the outer registry
    gens4 = get_symmetries(mc_hole_billiard(hf_d4_symmetric_hole()))
    @test typeof.(gens4) == (YAxisReflection, XAxisReflection)
    @test [g.sym_id for g in gens4] == [1, 3]
    gens2 = get_symmetries(mc_hole_billiard(hf_d2_symmetric_hole()))
    @test typeof.(gens2) == (XAxisReflection,)
    @test only(gens2).sym_id == 3
    @test get_symmetries(mc_hole_billiard(DomainHole(hf_circle_domain(0.2, 0.3, 0.4, 5)))) == ()
    @test get_symmetries(hf_rect_with_hole()) == ()
    # the generators of G' generate G'
    b = mc_hole_billiard(hf_symmetrized_disc())
    @test get_symmetries(b) == get_symmetries(get_symmetry_group(b))
    # a registry is accepted directly as before, and the billiard without holes is unchanged
    @test typeof.(get_symmetries(RectangleBilliard(1.0, 0.6))) == (YAxisReflection, XAxisReflection)
end

@testset "_symmetry_group" begin
    # the domain-dispatched method behind get_symmetry_group: simply connected domains return the registry itself
    circ = CircleBilliard(1.0)
    @test _symmetry_group(circ.fundamental_domain, circ) === circ.symmetries
    rect = RectangleBilliard(1.0, 0.6)
    @test _symmetry_group(rect.fundamental_domain, rect) === rect.symmetries
    # a multiply connected domain returns the compatible subgroup, original elements and registry order
    for (hole, ids) in ((hf_d4_symmetric_hole(), [1, 2, 3]), (hf_symmetrized_disc(), [1, 2, 3]), (hf_d2_symmetric_hole(), [3]),
                        (hf_half_disc_hole(), [3]), (DomainHole(hf_circle_domain(0.2, 0.3, 0.4, 5)), Int[]))
        b = mc_hole_billiard(hole)
        g = _symmetry_group(b.fundamental_domain, b)
        @test g isa SymmetryRegistry
        @test mc_ids(g) == ids
        @test mc_ids(g) == mc_ids(get_symmetry_group(b))
    end
    # the elements are the objects of the billiard registry, not copies with new ids
    b = mc_hole_billiard(hf_d2_symmetric_hole())
    @test only(collect(_symmetry_group(b.fundamental_domain, b))) == b.symmetries[3]
    # an empty registry stays empty whatever the holes
    nb = hf_rect_with_hole()
    @test length(_symmetry_group(nb.fundamental_domain, nb)) == 0
    # no holes in a multiply connected domain: G' = G
    bare = hf_HoleBilliard{Float64}(MultiplyConnectedDomain(hf_quarter_square_domain(), [], 1), D2_symmetry())
    @test mc_ids(_symmetry_group(bare.fundamental_domain, bare)) == [1, 2, 3]
end

@testset "_first_incompatibility" begin
    # D2 disc about (0.5, 0) in a D2 billiard: Y is the first element that does not map the hole onto itself
    i, g = _first_incompatibility(mc_hole_billiard(hf_d2_symmetric_hole()).fundamental_domain, D2_symmetry())
    @test i == 1
    @test g isa YAxisReflection
    # the first failing hole is reported by index
    dom = mc_hole_billiard(hf_d4_symmetric_hole(), hf_half_disc_hole()).fundamental_domain
    i, g = _first_incompatibility(dom, D2_symmetry())
    @test i == 2
    @test g isa YAxisReflection
    # a DomainHole breaks the first registered element
    i, g = _first_incompatibility(hf_d2_with_hole().fundamental_domain, D2_symmetry())
    @test (i, typeof(g)) == (1, YAxisReflection)
    # every hole invariant: nothing
    @test _first_incompatibility(mc_hole_billiard(hf_d4_symmetric_hole()).fundamental_domain, D2_symmetry()) === nothing
    @test _first_incompatibility(mc_hole_billiard(hf_symmetrized_disc()).fundamental_domain, D2_symmetry()) === nothing
    # empty registry: nothing to violate
    @test _first_incompatibility(hf_d2_with_hole().fundamental_domain, register_symmetries()) === nothing
    # a domain without holes
    @test _first_incompatibility(RectangleBilliard(1.0, 0.6).fundamental_domain, D2_symmetry()) === nothing
    @test _first_incompatibility(hf_quarter_square_domain(), D2_symmetry()) === nothing
end

@testset "_log_symmetry_group" begin
    # G' != G: one @info line with both orders and the first failing hole and element
    b = mc_hole_billiard(hf_d2_symmetric_hole())
    @test_logs (:info, "hf: outer symmetry group G has order 4, symmetry group G' has order 2; hole 1 is not invariant under YAxisReflection(1)") _log_symmetry_group("hf", b)
    # G' == G: no failing-hole clause
    @test_logs (:info, "hf: outer symmetry group G has order 4, symmetry group G' has order 4") _log_symmetry_group("hf", mc_hole_billiard(hf_d4_symmetric_hole()))
    @test_logs (:info, "hf: outer symmetry group G has order 4, symmetry group G' has order 4") _log_symmetry_group("hf", mc_hole_billiard(hf_symmetrized_disc()))
    # DomainHole: G' trivial
    @test_logs (:info, r"G' has order 1; hole 1 is not invariant under YAxisReflection") _log_symmetry_group("hf", hf_d2_with_hole())
    # the function returns nothing
    @test (@test_logs (:info, r"G has order 4") _log_symmetry_group("hf", b)) === nothing
    # an empty registry logs nothing
    @test (@test_logs _log_symmetry_group("hf", hf_rect_with_hole())) === nothing
    @test (@test_logs _log_symmetry_group("hf", hf_shifted_annulus())) === nothing
    # the origin-centered production constructors log their groups on every construction (G' = G, so no failing-hole clause)
    @test_logs (:info, "SinaiBilliard: outer symmetry group G has order 8, symmetry group G' has order 8") SinaiBilliard(1.0, 0.3)
    @test_logs (:info, "AnnularBilliard: outer symmetry group G has order 2, symmetry group G' has order 2") AnnularBilliard(0.3, 0.2)
    @test_logs (:info, "SquareWithinSquareBilliard: outer symmetry group G has order 8, symmetry group G' has order 8") SquareWithinSquareBilliard(1.0; positions = [0.0], sizes = [0.4])
    # an explicit shifted center, or the default two-obstacle square, has an empty registry and stays silent
    @test_logs SinaiBilliard(1.0, 0.3; center = [0.5, -0.2])
    @test_logs AnnularBilliard(0.3, 0.2; center = [0.5, 0.1])
    @test_logs SquareWithinSquareBilliard(1.0)
end
