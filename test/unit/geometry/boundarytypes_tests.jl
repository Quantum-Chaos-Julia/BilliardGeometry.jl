# Unit tests for boundary conditions and boundary/domain lookup on billiards:
#   src/geometry/boundarytypes.jl   (SpecularReflection, Transparent, PeriodicX, SymmetryWall, QuantumSolverIgnore,
#                                    get_boundary_curves, get_all_domains, get_domain, get_all_curves, get_curve,
#                                    update_boundary_condition, genus, boundary_components, geometric_boundary_components)
# Reference objects: CircleBilliard(1.0), the quarter-disk fundamental domain with the arc (SpecularReflection, segment 1)
# and the two symmetry walls (segments 2, 3); RectangleBilliard(1.0, 0.6) with two physical edges (segments 1, 2) and two
# symmetry walls (segments 3, 4); StadiumBilliard(1.0), a composite domain of two subdomains (ids 1, 2) joined by a
# Transparent curve; AnnularBilliard(0.3, 0.1), a multiply connected domain of genus 1 with the outer circle first.
# get_all_curves is called on billiards and on simple and composite domains (the CircleBilliard fundamental domain and the
# StadiumBilliard fundamental domain, whose curve count is the sum over its subdomains). update_boundary_condition is called
# on a simple-domain billiard (CircleBilliard) and on composite-domain billiards (StadiumBilliard with Float64 and Float32,
# MushroomBilliard with three subdomains), where the result must keep the billiard type and its element type, change only the
# targeted curve and leave the original billiard untouched. No eigen solves, no pinned wavefunction values. Helper functions
# and constants use the prefix `bt_`.

using StaticArrays: SVector

# Fundamental-domain fixtures; each is cheap to build.
bt_circle() = CircleBilliard(1.0)
bt_rect() = RectangleBilliard(1.0, 0.6)
bt_stadium() = StadiumBilliard(1.0)
bt_annulus() = AnnularBilliard(0.3, 0.1)

# True if `updated` has the same curves (by domain and segment id, in the same order) as `original` and every curve other than
# (domain_id, segment_id) keeps its boundary condition.
function bt_other_bcs_unchanged(original, updated, domain_id, segment_id)
    old = get_all_curves(original)
    new = get_all_curves(updated)
    length(old) == length(new) || return false
    for (c_old, c_new) in zip(old, new)
        (c_old.domain_id == c_new.domain_id && c_old.segment_id == c_new.segment_id) || return false
        (c_old.domain_id == domain_id && c_old.segment_id == segment_id) && continue
        c_old.bc == c_new.bc || return false
    end
    return true
end

@testset "boundary condition types" begin
    @test SpecularReflection() isa AbsBoundaryCondition
    @test QuantumSolverIgnore() isa AbsBoundaryCondition
    @test Transparent(2) isa AbsBoundaryCondition
    @test PeriodicX(3) isa AbsBoundaryCondition
    @test SymmetryWall(1, 2) isa AbsBoundaryCondition

    @test Transparent(2).next_id == 2
    @test PeriodicX(3).next_id == 3
    # next_id is stored as given, as Int64
    @test Transparent(5).next_id isa Int64
    @test PeriodicX(5).next_id isa Int64

    # parameter-free conditions compare equal, parametrized ones compare by their fields
    @test SpecularReflection() == SpecularReflection()
    @test Transparent(2) == Transparent(2)
    @test Transparent(2) != Transparent(3)
    @test PeriodicX(1) != PeriodicX(2)
    @test SymmetryWall(1, 2) == SymmetryWall(1, 2)
    @test SymmetryWall(1, 2) != SymmetryWall(2, 1)
end

@testset "SymmetryWall" begin
    # Round-trip: fields are read back unchanged, and a `SymmetryWall`-tagged
    # curve is excluded from `get_boundary_curves`'s physical-boundary filter
    # (`typeof(crv.bc) <: SpecularReflection`), since it is a fundamental-
    # domain cut, not part of the physical boundary. `RectangleBilliard`
    # tags its two symmetry walls with `SymmetryWall(1,2)`/`SymmetryWall(3,2)`
    # (see rectangle.jl), giving a real fixture rather than a synthetic curve.
    wall = SymmetryWall(2, 3)
    @test wall.sym_id == 2
    @test wall.sector_id == 3

    rect = RectangleBilliard(1.0, 0.6)
    bc = get_boundary_curves(rect)
    @test length(bc) == 2 # only the two SpecularReflection curves; both SymmetryWall curves excluded
    @test all(c -> typeof(c.bc) <: SpecularReflection, bc)

    walls = filter(c -> c.bc isa SymmetryWall, rect.fundamental_domain.boundary)
    @test length(walls) == 2
    @test Set((w.bc.sym_id, w.bc.sector_id) for w in walls) == Set([(1, 2), (3, 2)])
end

@testset "get_boundary_curves" begin
    # circle: only the arc is physical
    circle = bt_circle()
    bc = get_boundary_curves(circle)
    @test length(bc) == 1
    @test bc[1] isa CircleSegment
    @test bc[1].radius == 1.0
    @test bc[1].bc isa SpecularReflection

    # the billiard method forwards to its fundamental domain
    @test get_boundary_curves(circle.fundamental_domain) == bc

    # rectangle: the two physical edges are the right and top sides
    rect = bt_rect()
    bc = get_boundary_curves(rect)
    @test Set(c.segment_id for c in bc) == Set([1, 2])
    @test all(c -> c isa LineSegment, bc)

    # composite stadium: the Transparent joint and both symmetry walls are excluded
    bc = get_boundary_curves(bt_stadium())
    @test length(bc) == 2
    @test all(c -> c.bc isa SpecularReflection, bc)
end

@testset "get_all_domains" begin
    circle = bt_circle()
    doms = get_all_domains(circle)
    @test length(doms) == 1
    @test doms[1] === circle.fundamental_domain

    stadium = bt_stadium()
    doms = get_all_domains(stadium)
    @test length(doms) == 2
    @test [d.id for d in doms] == [1, 2]
    @test doms == stadium.fundamental_domain.subdomains
end

@testset "get_domain" begin
    circle = bt_circle()
    @test get_domain(circle, 1) === circle.fundamental_domain
    @test get_domain(circle, 5) === nothing

    stadium = bt_stadium()
    @test get_domain(stadium, 1) === stadium.fundamental_domain.subdomains[1]
    @test get_domain(stadium, 2) === stadium.fundamental_domain.subdomains[2]
    @test get_domain(stadium, 3) === nothing
end

@testset "get_all_curves" begin
    # a simple-domain billiard lists every curve of its boundary, whatever the boundary condition
    circle = bt_circle()
    crvs = get_all_curves(circle)
    @test length(crvs) == 3
    @test crvs == circle.fundamental_domain.boundary
    @test [c.segment_id for c in crvs] == [1, 2, 3]

    rect = bt_rect()
    @test length(get_all_curves(rect)) == 4
    @test count(c -> c.bc isa SymmetryWall, get_all_curves(rect)) == 2
    @test count(c -> c.bc isa SpecularReflection, get_all_curves(rect)) == 2
end

@testset "get_all_curves (::AbsDomain)" begin
    # G8: get_all_curves(domain) iterates the domain itself for a simple domain instead of wrapping it, MethodError on SimpleDomain
    circle = bt_circle()
    dom = circle.fundamental_domain
    crvs = get_all_curves(dom)
    @test crvs isa Vector{AbsCurve}
    @test length(crvs) == 3
    @test crvs == dom.boundary
    @test crvs == get_all_curves(circle)

    # G8: composite domain, the curves of all subdomains in subdomain order
    stadium = bt_stadium()
    cdom = stadium.fundamental_domain
    crvs = get_all_curves(cdom)
    @test crvs isa Vector{AbsCurve}
    @test length(crvs) == sum(length(sd.boundary) for sd in cdom.subdomains)
    @test length(crvs) == 7
    @test crvs == vcat((sd.boundary for sd in cdom.subdomains)...)
    @test crvs == get_all_curves(stadium)
    @test [c.domain_id for c in crvs] == [1, 1, 1, 2, 2, 2, 2]
end

@testset "get_curve" begin
    circle = bt_circle()
    arc = get_curve(circle, 1, 1)
    @test arc isa CircleSegment
    @test arc.bc isa SpecularReflection
    ywall = get_curve(circle, 1, 2)
    @test ywall isa LineSegment
    @test ywall.bc == SymmetryWall(1, 2)
    xwall = get_curve(circle, 1, 3)
    @test xwall.bc == SymmetryWall(3, 2)

    # no curve with the requested (domain_id, segment_id) pair
    @test get_curve(circle, 1, 9) === nothing
    @test get_curve(circle, 2, 1) === nothing

    # composite curves are searched by their subcurves' ids
    sub1 = LineSegment([0.0, 0.0], [1.0, 0.0]; domain_id = 4, segment_id = 1)
    sub2 = LineSegment([1.0, 0.0], [1.0, 1.0]; domain_id = 4, segment_id = 2)
    comp = CompositeCurve(AbsCurve[sub1, sub2])
    @test get_curve(comp, 4, 2) === sub2
    @test get_curve(comp, 4, 3) === nothing
end

@testset "update_boundary_condition" begin
    # turning the y-wall of the circle into a physical boundary adds a physical curve to the fundamental domain
    circle = bt_circle()
    updated = update_boundary_condition(circle, 1, 2, SpecularReflection())
    @test updated isa CircleBilliard{Float64}
    @test get_curve(updated, 1, 2).bc isa SpecularReflection
    @test length(get_boundary_curves(updated)) == 2
    # the other curves keep their conditions, and the original billiard is unchanged
    @test get_curve(updated, 1, 1).bc isa SpecularReflection
    @test get_curve(updated, 1, 3).bc == SymmetryWall(3, 2)
    @test get_curve(circle, 1, 2).bc == SymmetryWall(1, 2)
    @test length(get_boundary_curves(circle)) == 1
end

@testset "update_boundary_condition (StadiumBilliard)" begin
    # G7: @set on a composite-domain billiard rebuilds it through a constructor that does not exist,
    # MethodError: no constructor StadiumBilliard(::CompositeDomain, ::SymmetryRegistry)
    stadium = bt_stadium()
    updated = update_boundary_condition(stadium, 1, 1, QuantumSolverIgnore())
    @test updated isa StadiumBilliard{Float64}
    @test updated.symmetries == stadium.symmetries
    @test length(updated.fundamental_domain.subdomains) == 2
    # the target curve changed and kept its geometry
    target = get_curve(updated, 1, 1)
    @test target isa CircleSegment
    @test target.bc isa QuantumSolverIgnore
    @test target.radius == 1.0
    @test target.domain_id == 1 && target.segment_id == 1
    # all other curves keep their conditions
    @test bt_other_bcs_unchanged(stadium, updated, 1, 1)
    @test get_curve(updated, 1, 2).bc == Transparent(2)
    @test get_curve(updated, 2, 4).bc == Transparent(1)
    # the original billiard is untouched
    @test get_curve(stadium, 1, 1).bc isa SpecularReflection
    @test length(get_boundary_curves(stadium)) == 2
    @test length(get_boundary_curves(updated)) == 1

    # a curve of the second subdomain: the Transparent joint becomes a specular wall
    updated2 = update_boundary_condition(stadium, 2, 4, SpecularReflection())
    @test updated2 isa StadiumBilliard{Float64}
    @test get_curve(updated2, 2, 4).bc isa SpecularReflection
    @test bt_other_bcs_unchanged(stadium, updated2, 2, 4)
    @test get_curve(stadium, 2, 4).bc == Transparent(1)

    # G7: the element type is preserved
    stadium32 = StadiumBilliard(1.0f0)
    updated32 = update_boundary_condition(stadium32, 1, 1, QuantumSolverIgnore())
    @test updated32 isa StadiumBilliard{Float32}
    @test get_curve(updated32, 1, 1).bc isa QuantumSolverIgnore
    @test get_curve(updated32, 1, 1).radius isa Float32
    @test bt_other_bcs_unchanged(stadium32, updated32, 1, 1)
    @test get_curve(stadium32, 1, 1).bc isa SpecularReflection
end

@testset "update_boundary_condition (MushroomBilliard)" begin
    # G7: same MethodError as for the stadium, a composite domain with three subdomains (cap circle, triangle, stem)
    mushroom = MushroomBilliard(1.0)
    updated = update_boundary_condition(mushroom, 1, 1, QuantumSolverIgnore())
    @test updated isa MushroomBilliard{Float64}
    @test updated.symmetries == mushroom.symmetries
    @test length(updated.fundamental_domain.subdomains) == 3
    @test get_curve(updated, 1, 1) isa LineSegment
    @test get_curve(updated, 1, 1).bc isa QuantumSolverIgnore
    @test bt_other_bcs_unchanged(mushroom, updated, 1, 1)
    @test get_curve(updated, 1, 3).bc == Transparent(2)
    # the original billiard is untouched
    @test get_curve(mushroom, 1, 1).bc isa SpecularReflection

    # a curve of the third subdomain (stem bottom, segment 2)
    updated3 = update_boundary_condition(mushroom, 3, 2, QuantumSolverIgnore())
    @test updated3 isa MushroomBilliard{Float64}
    @test get_curve(updated3, 3, 2).bc isa QuantumSolverIgnore
    @test bt_other_bcs_unchanged(mushroom, updated3, 3, 2)
    @test get_curve(mushroom, 3, 2).bc isa SpecularReflection

    # G7: the element type is preserved
    mushroom32 = MushroomBilliard(1.0f0)
    updated32 = update_boundary_condition(mushroom32, 1, 1, QuantumSolverIgnore())
    @test updated32 isa MushroomBilliard{Float32}
    @test get_curve(updated32, 1, 1).bc isa QuantumSolverIgnore
end

@testset "genus" begin
    @test genus(bt_circle()) == 0
    @test genus(bt_circle().fundamental_domain) == 0
    @test genus(bt_rect()) == 0
    # composite domains are simply connected
    @test genus(bt_stadium()) == 0
    # the annulus has one hole
    annulus = bt_annulus()
    @test genus(annulus) == 1
    @test genus(annulus.fundamental_domain) == 1
end

@testset "boundary_components" begin
    # a simple domain has a single component made of its connected physical curves
    circle = bt_circle()
    comps = boundary_components(circle)
    @test length(comps) == 1
    @test comps[1] == get_boundary_curves(circle)
    @test boundary_components(circle.fundamental_domain) == comps

    rect = bt_rect()
    comps = boundary_components(rect)
    @test length(comps) == 1
    @test length(comps[1]) == 2

    # composite domain: still one component with the two physical curves
    comps = boundary_components(bt_stadium())
    @test length(comps) == 1
    @test length(comps[1]) == 2

    # annulus: the outer circle comes first, then the hole
    annulus = bt_annulus()
    comps = boundary_components(annulus)
    @test length(comps) == genus(annulus) + 1
    @test length(comps[1]) == 1 && length(comps[2]) == 1
    @test comps[1][1].radius == 1.0 && comps[1][1].domain_id == 1
    @test comps[2][1].radius == 0.3 && comps[2][1].domain_id == 2
end

@testset "geometric_boundary_components" begin
    # unlike boundary_components, the symmetry walls closing the fundamental domain are kept
    circle = bt_circle()
    comps = geometric_boundary_components(circle)
    @test length(comps) == 1
    @test length(comps[1]) == 3
    @test count(c -> c.bc isa SymmetryWall, comps[1]) == 2
    @test count(c -> c.bc isa SpecularReflection, comps[1]) == 1
    @test geometric_boundary_components(circle.fundamental_domain) == comps

    rect = bt_rect()
    comps = geometric_boundary_components(rect)
    @test length(comps) == 1
    @test length(comps[1]) == 4

    # the Transparent joint of the composite stadium is not part of the exterior
    comps = geometric_boundary_components(bt_stadium())
    @test length(comps) == 1
    @test !any(c -> c.bc isa Transparent, comps[1])

    # annulus: two components, outer first
    comps = geometric_boundary_components(bt_annulus())
    @test length(comps) == 2
    @test comps[1][1].radius == 1.0
    @test comps[2][1].radius == 0.3
end

# Billiards with holes (MultiplyConnectedDomain). References are the fixtures of holes_fixtures.jl: concentric and eccentric
# annulus, rectangle [0,3] x [0,2] with the rectangular hole [1,1.5] x [0.5,1.5] (no symmetry), and the D2 billiard whose outer
# region is the quarter square [0,1]^2 with the circular hole of radius 0.2 about (0.3, 0.4) in the full frame.

@testset "get_outer_boundary_curves" begin
    # without holes it equals get_boundary_curves, for domains and billiards
    circle = bt_circle()
    @test get_outer_boundary_curves(circle) == get_boundary_curves(circle)
    @test get_outer_boundary_curves(circle.fundamental_domain) == get_boundary_curves(circle.fundamental_domain)
    stadium = bt_stadium()
    @test get_outer_boundary_curves(stadium) == get_boundary_curves(stadium)
    @test length(get_outer_boundary_curves(stadium)) == 2

    # with holes only the outer physical curves of the fundamental frame are returned
    ann = hf_concentric_annulus()
    outer = get_outer_boundary_curves(ann)
    @test length(outer) == 1
    @test outer[1].radius == 1.0 && outer[1].domain_id == 1
    @test get_outer_boundary_curves(ann.fundamental_domain) == outer
    @test length(get_outer_boundary_curves(hf_rect_with_hole())) == 4
    # D2: the two physical edges of the quarter square, the symmetry walls are not physical
    d2 = get_outer_boundary_curves(hf_d2_with_hole())
    @test length(d2) == 2
    @test all(c -> c.domain_id == 1 && c.bc isa SpecularReflection, d2)
end

@testset "get_hole_boundary_curves" begin
    # no holes: empty
    @test isempty(get_hole_boundary_curves(bt_circle()))
    @test isempty(get_hole_boundary_curves(bt_stadium().fundamental_domain))
    @test get_hole_boundary_curves(bt_circle()) isa Vector{AbsCurve}

    ann = hf_concentric_annulus()
    holes = get_hole_boundary_curves(ann)
    @test length(holes) == 1
    @test holes[1].radius == 0.3 && holes[1].domain_id == 2
    @test get_hole_boundary_curves(ann.fundamental_domain) == holes
    # rectangular hole: four curves of the loop, counterclockwise
    rect = get_hole_boundary_curves(hf_rect_with_hole())
    @test length(rect) == 4
    @test all(c -> c.domain_id == 2, rect)
    # D2: the hole is in the full frame
    d2 = get_hole_boundary_curves(hf_d2_with_hole())
    @test length(d2) == 1
    @test d2[1].center == SVector(0.3, 0.4)
    # holes of several objects: loops one after the other
    dom = MultiplyConnectedDomain(hf_rect_domain(0.0, 0.0, 4.0, 4.0, 1), [hf_circle_domain(0.3, 3.0, 3.0, 2), hf_two_loop_hole(5)], 1)
    @test [c.domain_id for c in get_hole_boundary_curves(dom)] == [2, 5, 5, 5, 5, 6, 6, 6, 6]
end

@testset "get_boundary_curves (MultiplyConnectedDomain)" begin
    # domain level: the outer curves only
    ann = hf_concentric_annulus()
    @test get_boundary_curves(ann.fundamental_domain) == get_outer_boundary_curves(ann)
    @test length(get_boundary_curves(ann.fundamental_domain)) == 1
    @test length(get_boundary_curves(hf_rect_with_hole().fundamental_domain)) == 4

    # billiard level, empty registry: outer curves followed by the hole curves (the frames coincide)
    bc = get_boundary_curves(ann)
    @test length(bc) == 2
    @test [c.domain_id for c in bc] == [1, 2]
    @test bc[1] == get_outer_boundary_curves(ann)[1]
    @test bc[2] == get_hole_boundary_curves(ann)[1]
    rect = get_boundary_curves(hf_rect_with_hole())
    @test length(rect) == 8
    @test [c.domain_id for c in rect] == [1, 1, 1, 1, 2, 2, 2, 2]

    # nonempty registry with a hole: the frames differ, ArgumentError pointing to full_boundary_components
    d2 = hf_d2_with_hole()
    @test_throws ArgumentError get_boundary_curves(d2)
    @test_throws "full_boundary_components" get_boundary_curves(d2)
    # with a nonempty registry but no holes the outer curves are returned
    no_holes = hf_HoleBilliard{Float64}(MultiplyConnectedDomain(hf_quarter_square_domain(), [], 1), D2_symmetry())
    @test length(get_boundary_curves(no_holes)) == 2
end

@testset "get_all_domains (billiards with holes)" begin
    # outer subdomains followed by the simple domains of the holes
    ann = hf_concentric_annulus()
    doms = get_all_domains(ann)
    @test [d.id for d in doms] == [1, 2]
    @test doms[1] === ann.fundamental_domain.subdomains[1]
    @test doms[2] === ann.fundamental_domain.holes[1].domain
    @test get_all_domains(ann.fundamental_domain) == doms
    @test [d.id for d in get_all_domains(hf_d2_with_hole())] == [1, 5]
    # composite outer region and a composite hole
    dom = MultiplyConnectedDomain(StadiumBilliard(1.0).fundamental_domain, [hf_two_loop_hole(5)], 1)
    @test [d.id for d in get_all_domains(dom)] == [1, 2, 5, 6]
    # holes are found by get_domain and contribute to get_all_curves
    @test get_domain(ann, 2) === ann.fundamental_domain.holes[1].domain
    @test get_domain(ann, 3) === nothing
    @test length(get_all_curves(ann)) == 2
    @test length(get_all_curves(hf_rect_with_hole())) == 8
    # the symmetry walls of the D2 quarter are listed too: 4 outer curves and the hole circle
    @test length(get_all_curves(hf_d2_with_hole())) == 5
    @test get_curve(ann, 2, 1) isa CircleSegment
    @test get_curve(ann, 2, 1).radius == 0.3
end

@testset "update_boundary_condition (billiards with holes)" begin
    # a hole curve is found by (domain_id, segment_id) in the holes
    ann = hf_eccentric_annulus()
    updated = update_boundary_condition(ann, 2, 1, QuantumSolverIgnore())
    @test updated isa AnnularBilliard{Float64}
    @test get_curve(updated, 2, 1).bc isa QuantumSolverIgnore
    @test get_curve(updated, 2, 1).radius == 0.2
    @test updated.fundamental_domain.holes[1] isa DomainHole
    @test updated.symmetries == ann.symmetries
    # the other curve keeps its condition and the original billiard is unchanged
    @test get_curve(updated, 1, 1).bc isa SpecularReflection
    @test get_curve(ann, 2, 1).bc isa SpecularReflection
    @test bt_other_bcs_unchanged(ann, updated, 2, 1)

    # an outer curve of the multiply connected domain is found in the subdomains
    updated_outer = update_boundary_condition(ann, 1, 1, QuantumSolverIgnore())
    @test get_curve(updated_outer, 1, 1).bc isa QuantumSolverIgnore
    @test get_curve(updated_outer, 2, 1).bc isa SpecularReflection
    @test bt_other_bcs_unchanged(ann, updated_outer, 1, 1)

    # a curve of a multi-curve hole
    rect = hf_rect_with_hole()
    updated_rect = update_boundary_condition(rect, 2, 3, QuantumSolverIgnore())
    @test updated_rect isa hf_HoleBilliard{Float64}
    @test get_curve(updated_rect, 2, 3).bc isa QuantumSolverIgnore
    @test bt_other_bcs_unchanged(rect, updated_rect, 2, 3)
    @test get_curve(rect, 2, 3).bc isa SpecularReflection

    # unknown (domain_id, segment_id): ArgumentError
    @test_throws ArgumentError update_boundary_condition(ann, 3, 1, QuantumSolverIgnore())
    @test_throws ArgumentError update_boundary_condition(ann, 2, 9, QuantumSolverIgnore())
end

@testset "boundary_components (billiards with holes)" begin
    # outer component first, then one component per hole loop
    for (b, lens, ids) in ((hf_concentric_annulus(), [1, 1], [1, 2]), (hf_eccentric_annulus(), [1, 1], [1, 2]),
                           (hf_rect_with_hole(), [4, 4], [1, 2]), (hf_d2_with_hole(), [2, 1], [1, 5]))
        comps = boundary_components(b)
        @test [length(c) for c in comps] == lens
        @test [c[1].domain_id for c in comps] == ids
        @test boundary_components(b.fundamental_domain) == comps
        @test length(comps) == genus(b) + 1
    end
    # composite hole of two loops: two hole components
    dom = MultiplyConnectedDomain(hf_rect_domain(0.0, 0.0, 4.0, 4.0, 1), [hf_circle_domain(0.3, 3.0, 3.0, 2), hf_two_loop_hole(5)], 1)
    @test [length(c) for c in boundary_components(dom)] == [4, 1, 4, 4]
    @test length(boundary_components(dom)) == genus(dom) + 1
end

@testset "boundary_component_is_hole" begin
    @test boundary_component_is_hole(bt_circle()) == [false]
    @test boundary_component_is_hole(bt_stadium()) == [false]
    @test boundary_component_is_hole(bt_circle().fundamental_domain) == [false]
    @test boundary_component_is_hole(hf_concentric_annulus()) == [false, true]
    @test boundary_component_is_hole(hf_eccentric_annulus()) == [false, true]
    @test boundary_component_is_hole(hf_rect_with_hole()) == [false, true]
    @test boundary_component_is_hole(hf_d2_with_hole()) == [false, true]
    dom = MultiplyConnectedDomain(hf_rect_domain(0.0, 0.0, 4.0, 4.0, 1), [hf_circle_domain(0.3, 3.0, 3.0, 2), hf_two_loop_hole(5)], 1)
    @test boundary_component_is_hole(dom) == [false, true, true, true]
    # one flag per boundary component
    @test length(boundary_component_is_hole(dom)) == length(boundary_components(dom))
end

@testset "geometric_boundary_components (billiards with holes)" begin
    # the outer component keeps its symmetry walls (fundamental frame), then one closed loop per hole loop
    @test [length(c) for c in geometric_boundary_components(hf_concentric_annulus())] == [1, 1]
    @test [length(c) for c in geometric_boundary_components(hf_rect_with_hole())] == [4, 4]
    d2 = geometric_boundary_components(hf_d2_with_hole())
    @test [length(c) for c in d2] == [4, 1]
    @test count(c -> c.bc isa SymmetryWall, d2[1]) == 2
    @test length(geometric_boundary_components(hf_d2_with_hole().fundamental_domain)) == 2
end
