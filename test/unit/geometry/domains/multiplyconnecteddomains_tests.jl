# Unit tests for src/geometry/domains/multiplyconnecteddomains.jl: the MultiplyConnectedDomain type, its constructor and its
#   domain queries (genus, hole_domain_ids, is_hole_curve, component_of_domain_id, validate_holes).
# Reference objects (fixtures from holes_fixtures.jl): the square [0,4]^2 as outer region with a circular hole (domain id 2) and a
# square hole [2,3]^2 (domain id 3); the composite outer region of StadiumBilliard; the quarter-square outer region of the D2
# rectangle. The curve `orientation` field is never read to decide holes: a hole stored with orientation -1 is still a hole, and a
# concave outer piece (inverted circle, orientation -1, region outside the circle) is not. No eigen solves. Helpers use the prefix `mc_`.

using StaticArrays: SVector

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
    @test hole_domain_ids(hf_concentric_annulus()) == [2]
    @test hole_domain_ids(hf_d2_with_hole()) == [5]
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
    ann = hf_concentric_annulus()
    outer_curve, hole_curve = get_all_curves(ann)
    @test !is_hole_curve(ann, outer_curve)
    @test is_hole_curve(ann, hole_curve)
    # the Sinai circle is a hole curve, the four sides of its square are not
    sinai = SinaiBilliard(1.0, 0.3)
    @test [is_hole_curve(sinai, c) for c in get_all_curves(sinai)] == [false, false, false, false, true]
    @test !is_hole_curve(mc_outer(), outer.boundary[1])
end

@testset "component_of_domain_id" begin
    outer = mc_outer()
    dom = MultiplyConnectedDomain(outer, [mc_circle_hole(2), hf_two_loop_hole(5)], 1)
    # 0 for the outer region, i for the i-th hole object (both loops of the composite hole share group 2)
    @test component_of_domain_id(dom) == Dict(1 => 0, 2 => 1, 5 => 2, 6 => 2)
    @test component_of_domain_id(hf_concentric_annulus()) == Dict(1 => 0, 2 => 1)
    @test component_of_domain_id(hf_d2_with_hole()) == Dict(1 => 0, 5 => 1)
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
