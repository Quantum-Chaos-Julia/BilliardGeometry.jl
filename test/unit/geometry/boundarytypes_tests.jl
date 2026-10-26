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
