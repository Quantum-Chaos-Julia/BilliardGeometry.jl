# Unit tests for the annular billiard constructor and structure:
#   src/geometry/billiards/annular.jl   (AnnularBilliard{T}, AnnularBilliard(R_inner, shift; R_outer, center))
# Reference object: a MultiplyConnectedDomain with id 1, genus 1 and an empty symmetry registry, made of an outer SimpleDomain
# (full CircleSegment, orientation 1, domain id 1) and a DomainHole around a SimpleDomain holding the full counterclockwise
# CircleSegment (orientation +1, domain id 2, segment id 1) of the hole, whose centre is displaced by `shift` along the x axis.
# Closed forms: area pi (R_outer^2 - R_inner^2), membership in the ring / hole / outside. Guards: R_outer > 0, R_inner > 0,
# R_inner < R_outer, and the hole must lie inside the outer circle (validate_holes). No eigen solves.

using StaticArrays: SVector

@testset "AnnularBilliard fields" begin
    ann = AnnularBilliard(0.3, 0.2)
    @test ann isa AnnularBilliard{Float64}
    @test ann isa AbsBilliard
    @test ann.fundamental_domain isa MultiplyConnectedDomain{Float64}
    @test ann.symmetries isa SymmetryRegistry
    @test fieldnames(AnnularBilliard) == (:fundamental_domain, :symmetries)
    @test @inferred(AnnularBilliard(0.3, 0.2)) isa AnnularBilliard{Float64}
end

@testset "AnnularBilliard symmetries" begin
    @test isempty(AnnularBilliard(0.3, 0.2).symmetries)
end

@testset "AnnularBilliard curves" begin
    dom = AnnularBilliard(0.3, 0.2).fundamental_domain
    @test dom.id == 1
    # one outer subdomain with the outer circle, one DomainHole around a simple domain with the hole circle
    @test length(dom.subdomains) == 1
    @test length(dom.holes) == 1
    @test dom.holes[1] isa DomainHole
    outer_dom = dom.subdomains[1]
    hole_dom = dom.holes[1].domain
    @test outer_dom isa SimpleDomain{Float64}
    @test hole_dom isa SimpleDomain{Float64}
    @test outer_dom.id == 1
    @test hole_dom.id == 2
    @test length(outer_dom.boundary) == 1
    @test length(hole_dom.boundary) == 1
    outer = outer_dom.boundary[1]
    inner = hole_dom.boundary[1]
    @test outer isa CircleSegment{Float64}
    @test inner isa CircleSegment{Float64}
    @test outer.bc isa SpecularReflection
    @test inner.bc isa SpecularReflection

    # both circles are counterclockwise (orientation +1, hole region on the left), with their own domain id
    @test [outer.orientation, inner.orientation] == [1, 1]
    @test [outer.domain_id, inner.domain_id] == [1, 2]
    @test [outer.segment_id, inner.segment_id] == [1, 1]
    @test outer.arc_angle == 2 * pi
    @test inner.arc_angle == 2 * pi

    # outer circle of radius one about the origin, inner circle of radius R_inner about (shift, 0)
    @test outer.radius == 1.0
    @test outer.center == SVector(0.0, 0.0)
    @test inner.radius == 0.3
    @test inner.center == SVector(0.2, 0.0)
end

@testset "AnnularBilliard topology" begin
    ann = AnnularBilliard(0.3, 0.2)
    @test genus(ann) == 1
    @test genus(ann.fundamental_domain) == 1
    @test !hasfield(MultiplyConnectedDomain, :genus)
    # one component for the outer circle and one for the hole, one curve each
    comps = boundary_components(ann)
    @test length(comps) == genus(ann) + 1
    @test [length(c) for c in comps] == [1, 1]
    @test [c[1].domain_id for c in comps] == [1, 2]
    @test boundary_component_is_hole(ann) == [false, true]
    @test length(geometric_boundary_components(ann)) == genus(ann) + 1
    @test hole_domain_ids(ann) == [2]
    @test component_of_domain_id(ann) == Dict(1 => 0, 2 => 1)
end

@testset "AnnularBilliard geometry" begin
    # area pi (R_outer^2 - R_inner^2), concentric and eccentric
    ann = AnnularBilliard(0.3, 0.1)
    @test area(ann) ≈ pi * (1 - 0.3^2) rtol = 1e-12
    @test area(AnnularBilliard(0.3, 0.0)) ≈ pi * (1 - 0.3^2) rtol = 1e-12
    @test area(AnnularBilliard(0.2, 0.3)) ≈ pi * (1 - 0.2^2) rtol = 1e-12
    @test area(AnnularBilliard(0.3, 0.1; R_outer = 2.0)) ≈ pi * (4 - 0.3^2) rtol = 1e-12
    # without symmetry the fundamental area is the area
    @test fundamental_area(ann) ≈ area(ann) rtol = 1e-14
    # membership: ring, hole, outside
    @test is_inside(ann, SVector(0.8, 0.0)) === true
    @test is_inside(ann, SVector(0.1, 0.0)) === false      # hole center (shift 0.1)
    @test is_inside(ann, SVector(1.2, 0.0)) === false
    # a ring point close to the hole is inside, a hole point is not
    @test is_inside(ann, SVector(0.1, 0.35)) === true
    @test is_inside(ann, SVector(0.1, 0.25)) === false
    # no true corners anywhere
    @test isempty(corner_angles(ann))
    @test isempty(corner_angles(ann; fundamental = false))
    # full boundary: outer loop then hole loop
    @test length(full_boundary(ann)) == 2
    @test full_boundary_components(ann).is_hole == [false, true]
end

@testset "AnnularBilliard keywords" begin
    ann = AnnularBilliard(0.3, 0.2; R_outer = 2.0, center = [0.1, 0.0])
    dom = ann.fundamental_domain
    outer = dom.subdomains[1].boundary[1]
    inner = dom.holes[1].domain.boundary[1]
    @test outer.radius == 2.0
    @test outer.center == SVector(0.1, 0.0)
    @test inner.radius == 0.3
    @test inner.center ≈ SVector(0.3, 0.0) atol = 1e-14
    @test area(ann) ≈ pi * (4 - 0.09) rtol = 1e-12
end

@testset "AnnularBilliard argument guards" begin
    # R_outer > 0
    @test_throws ArgumentError AnnularBilliard(0.3, 0.2; R_outer = 0.0)
    @test_throws ArgumentError AnnularBilliard(0.3, 0.2; R_outer = -1.0)
    # R_inner > 0
    @test_throws ArgumentError AnnularBilliard(0.0, 0.2)
    @test_throws ArgumentError AnnularBilliard(-0.3, 0.2)
    # R_inner < R_outer
    @test_throws ArgumentError AnnularBilliard(1.0, 0.0)
    @test_throws ArgumentError AnnularBilliard(1.5, 0.0)
    @test_throws ArgumentError AnnularBilliard(0.5, 0.0; R_outer = 0.5)
    # the sample point of the hole loop must lie inside the outer circle (validate_holes)
    @test_throws ArgumentError AnnularBilliard(0.3, 3.0)
end

@testset "AnnularBilliard eltype promotion" begin
    ann32 = AnnularBilliard(0.3f0, 0.2f0)
    @test ann32 isa AnnularBilliard{Float32}
    @test ann32.fundamental_domain isa MultiplyConnectedDomain{Float32}
    @test all(c -> c isa CircleSegment{Float32}, get_all_curves(ann32))
    @test genus(ann32) == 1
    @test area(ann32) isa Float32
end
