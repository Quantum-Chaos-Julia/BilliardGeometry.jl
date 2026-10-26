# Unit tests for src/geometry/domains/holes.jl: the DomainHole type and its interface
#   (DomainHole, hole_domains, hole_leaves, hole_loop_count, get_all_domains, get_all_curves, get_boundary_curves,
#    boundary_components, is_inside, area, symmetry_compatible, _set_curve_bc).
# Reference objects (fixtures from holes_fixtures.jl): the circular hole of radius 0.2 about (0.3, 0.4), the rectangular hole
# [1, 1.5] x [0.5, 1.5], and a composite hole made of two disjoint rectangles (two loops). Areas, membership and loop counts are
# closed forms. No eigen solves. Helpers use the prefix `hl_`.

using StaticArrays: SVector
using BilliardGeometry: _set_curve_bc

# Absolute tolerance for closed-form areas of the holes (Green quadrature, worst observed 4e-16, margin 1000x).
const hl_ATOL = 1e-12

hl_circle_hole() = DomainHole(hf_circle_domain(0.2, 0.3, 0.4, 5))
hl_rect_hole() = DomainHole(hf_rect_domain(1.0, 0.5, 1.5, 1.5, 2))
# Two disjoint rectangles [0,1]x[0,1] (domain id 3) and [2,3]x[0,1] (domain id 4) as one composite hole object.
hl_two_loop_hole() = hf_two_loop_hole(3)

@testset "DomainHole constructor" begin
    h = hl_circle_hole()
    @test h isa DomainHole
    @test h isa AbsDomainHole
    @test !(h isa AbsDomain)
    @test h.domain isa SimpleDomain{Float64}
    @test fieldnames(DomainHole) == (:domain,)
    @test hl_two_loop_hole().domain isa CompositeDomain
    @test DomainHole(hl_two_loop_hole().domain) isa DomainHole{CompositeDomain}
end

@testset "hole_domains" begin
    h = hl_circle_hole()
    @test hole_domains(h) isa Vector{AbsDomain}
    @test length(hole_domains(h)) == 1
    @test hole_domains(h)[1] === h.domain
    comp = hl_two_loop_hole()
    @test hole_domains(comp)[1] === comp.domain
end

@testset "hole_leaves" begin
    h = hl_rect_hole()
    @test hole_leaves(h) isa Vector{AbsDomainHole}
    @test length(hole_leaves(h)) == 1
    @test hole_leaves(h)[1] === h
end

@testset "hole_loop_count" begin
    @test hole_loop_count(hl_circle_hole()) == 1
    @test hole_loop_count(hl_rect_hole()) == 1
    # a composite hole of two disjoint rectangles has two closed loops
    @test hole_loop_count(hl_two_loop_hole()) == 2
end

@testset "get_all_domains (DomainHole)" begin
    h = hl_circle_hole()
    @test get_all_domains(h) == AbsSimpleDomain[h.domain]
    comp = hl_two_loop_hole()
    @test [d.id for d in get_all_domains(comp)] == [3, 4]
end

@testset "get_all_curves (DomainHole)" begin
    @test length(get_all_curves(hl_circle_hole())) == 1
    @test get_all_curves(hl_circle_hole())[1] isa CircleSegment
    rect = hl_rect_hole()
    @test length(get_all_curves(rect)) == 4
    @test [c.segment_id for c in get_all_curves(rect)] == [1, 2, 3, 4]
    @test [c.domain_id for c in get_all_curves(hl_two_loop_hole())] == [3, 3, 3, 3, 4, 4, 4, 4]
end

@testset "get_boundary_curves (DomainHole)" begin
    h = hl_circle_hole()
    bc = get_boundary_curves(h)
    @test length(bc) == 1
    @test bc[1].radius == 0.2
    @test bc[1].center == SVector(0.3, 0.4)
    @test length(get_boundary_curves(hl_rect_hole())) == 4
    @test length(get_boundary_curves(hl_two_loop_hole())) == 8
    # a symmetry wall is not physical and is dropped from the physical boundary
    walled = DomainHole(RectangleBilliard(1.0, 1.0).fundamental_domain)
    @test length(get_boundary_curves(walled)) == 2
end

@testset "boundary_components (DomainHole)" begin
    # one closed loop per component, the hole region on the left (counterclockwise for a convex hole)
    comps = boundary_components(hl_circle_hole())
    @test length(comps) == 1
    @test length(comps[1]) == 1

    rect_comps = boundary_components(hl_rect_hole())
    @test length(rect_comps) == 1
    @test length(rect_comps[1]) == 4
    loop = rect_comps[1]
    # closed: every curve ends where the next starts, the last ends where the first starts
    for j in 1:4
        @test curve(loop[j], 1.0) ≈ curve(loop[mod1(j + 1, 4)], 0.0) atol = hl_ATOL
    end
    # counterclockwise: the first edge runs along +x at the bottom
    @test curve(loop[1], 0.0) ≈ [1.0, 0.5] atol = hl_ATOL
    @test curve(loop[1], 1.0) ≈ [1.5, 0.5] atol = hl_ATOL

    two = boundary_components(hl_two_loop_hole())
    @test length(two) == 2
    @test [length(c) for c in two] == [4, 4]
    @test [c[1].domain_id for c in two] == [3, 4]
end

@testset "is_inside (DomainHole)" begin
    circ = hl_circle_hole()
    @test is_inside(circ, SVector(0.3, 0.4)) === true
    @test is_inside(circ, SVector(0.45, 0.4)) === true
    @test is_inside(circ, SVector(0.6, 0.4)) === false
    @test is_inside(circ, SVector(0.0, 0.0)) === false
    # the hole is asymmetric: the mirror image of an inside point is outside
    @test is_inside(circ, SVector(-0.3, 0.4)) === false
    @test is_inside(circ, [SVector(0.3, 0.4), SVector(-0.3, 0.4)]) == [true, false]

    rect = hl_rect_hole()
    @test is_inside(rect, SVector(1.25, 1.0)) === true
    @test is_inside(rect, SVector(1.6, 1.0)) === false
    @test is_inside(rect, SVector(1.25, 0.4)) === false

    # a composite hole is the union of its pieces
    two = hl_two_loop_hole()
    @test is_inside(two, SVector(0.5, 0.5)) === true
    @test is_inside(two, SVector(2.5, 0.5)) === true
    @test is_inside(two, SVector(1.5, 0.5)) === false
    @test is_inside(two, [SVector(0.5, 0.5), SVector(1.5, 0.5), SVector(2.5, 0.5)]) == [true, false, true]
end

@testset "area (DomainHole)" begin
    @test area(hl_circle_hole()) ≈ pi * 0.2^2 atol = hl_ATOL
    @test area(hl_rect_hole()) ≈ 0.5 atol = hl_ATOL
    # sum over the loops
    @test area(hl_two_loop_hole()) ≈ 2.0 atol = hl_ATOL
    # the hole area is positive and independent of the curve orientation field (never read)
    inverted = DomainHole(SimpleDomain{Float64}(AbsCurve[CircleSegment(0.2, 2pi, 0.0, [0.3, 0.4]; orientation = -1, domain_id = 5)], [SVector(0.5, 0.4)], 5))
    @test area(inverted) ≈ pi * 0.2^2 atol = hl_ATOL
    # the quadrature tolerance is forwarded
    @test area(hl_circle_hole(); rtol = 1e-6) ≈ pi * 0.2^2 atol = 1e-8
end

@testset "symmetry_compatible (DomainHole)" begin
    h = hl_circle_hole()
    # a DomainHole has no symmetry: compatible only with the empty registry
    @test symmetry_compatible(h, register_symmetries()) === true
    @test symmetry_compatible(h, D2_symmetry()) === false
    @test symmetry_compatible(h, register_symmetries(XAxisReflection())) === false
    @test symmetry_compatible(hl_rect_hole(), register_symmetries()) === true
    @test symmetry_compatible(hl_two_loop_hole(), register_symmetries(XYAxisReflection())) === false
end

@testset "_set_curve_bc (DomainHole)" begin
    h = hl_rect_hole()
    new_h = _set_curve_bc(h, 2, 3, QuantumSolverIgnore())
    @test new_h isa DomainHole
    @test get_all_curves(new_h)[3].bc isa QuantumSolverIgnore
    # the other curves keep their condition, and the original hole is unchanged
    @test all(c -> c.bc isa SpecularReflection, get_all_curves(new_h)[[1, 2, 4]])
    @test all(c -> c.bc isa SpecularReflection, get_all_curves(h))
    # composite holes are searched by (domain_id, segment_id)
    two = hl_two_loop_hole()
    new_two = _set_curve_bc(two, 4, 2, Transparent(1))
    @test new_two isa DomainHole
    @test get_all_curves(new_two)[6].bc == Transparent(1)
    @test count(c -> !(c.bc isa SpecularReflection), get_all_curves(new_two)) == 1
    # no curve with that id pair: nothing
    @test _set_curve_bc(h, 9, 1, QuantumSolverIgnore()) === nothing
    @test _set_curve_bc(h, 2, 9, QuantumSolverIgnore()) === nothing
end
