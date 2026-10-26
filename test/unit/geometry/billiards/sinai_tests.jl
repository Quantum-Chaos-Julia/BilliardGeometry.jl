# Unit tests for the Sinai billiard constructor and structure:
#   src/geometry/billiards/sinai.jl   (SinaiBilliard{T}, SinaiBilliard(a, R; center))
# Reference object: a MultiplyConnectedDomain with id 1, genus 1 and an empty symmetry registry, made of an outer SimpleDomain
# (the four counterclockwise LineSegments of the square [-a, a]^2 about `center`: bottom, right, top, left, orientation 1, domain
# id 1, segment ids 1..4) and a DomainHole around a SimpleDomain holding the full counterclockwise CircleSegment of radius R
# about `center` (orientation +1, domain id 2, segment id 1).
# Closed forms: area 4 a^2 - pi R^2, membership in the square / hole / outside, square corner angles pi/2. Guards: a > 0,
# R > 0, R < a. The centre is arbitrary. No eigen solves.

using StaticArrays: SVector

@testset "SinaiBilliard fields" begin
    sin_b = SinaiBilliard(1.0, 0.3)
    @test sin_b isa SinaiBilliard{Float64}
    @test sin_b isa AbsBilliard
    @test sin_b.fundamental_domain isa MultiplyConnectedDomain{Float64}
    @test sin_b.symmetries isa SymmetryRegistry
    @test fieldnames(SinaiBilliard) == (:fundamental_domain, :symmetries)
    @test @inferred(SinaiBilliard(1.0, 0.3)) isa SinaiBilliard{Float64}
end

@testset "SinaiBilliard symmetries" begin
    sin_b = SinaiBilliard(1.0, 0.3)
    @test isempty(sin_b.symmetries)
    @test get_symmetries(sin_b) == ()
end

@testset "SinaiBilliard curves" begin
    dom = SinaiBilliard(1.0, 0.3).fundamental_domain
    @test dom.id == 1
    # one outer subdomain with the square, one DomainHole around a simple domain with the circle
    @test length(dom.subdomains) == 1
    @test length(dom.holes) == 1
    @test dom.holes[1] isa DomainHole
    outer_dom = dom.subdomains[1]
    hole_dom = dom.holes[1].domain
    @test outer_dom isa SimpleDomain{Float64}
    @test hole_dom isa SimpleDomain{Float64}
    @test [outer_dom.id, hole_dom.id] == [1, 2]
    @test length(outer_dom.boundary) == 4
    @test length(hole_dom.boundary) == 1
    @test all(c -> c isa LineSegment{Float64,SpecularReflection}, outer_dom.boundary)
    inner = hole_dom.boundary[1]
    @test inner isa CircleSegment{Float64,SpecularReflection}

    # the square and the circle are counterclockwise (orientation +1, hole region on the left), each with its own domain id
    @test [c.orientation for c in outer_dom.boundary] == [1, 1, 1, 1]
    @test [c.domain_id for c in outer_dom.boundary] == [1, 1, 1, 1]
    @test [c.segment_id for c in outer_dom.boundary] == [1, 2, 3, 4]
    @test inner.orientation == 1
    @test inner.domain_id == 2
    @test inner.segment_id == 1
    @test inner.arc_angle == 2 * pi

    # sides in the order bottom, right, top, left, each starting where the previous one ends
    sides = outer_dom.boundary
    @test [s.pt0 for s in sides] == [SVector(-1.0, -1.0), SVector(1.0, -1.0), SVector(1.0, 1.0), SVector(-1.0, 1.0)]
    @test [s.pt1 for s in sides] == [SVector(1.0, -1.0), SVector(1.0, 1.0), SVector(-1.0, 1.0), SVector(-1.0, -1.0)]
    @test all(s -> s.length == 2.0, sides)
    # the vertices are the end points of the curves, in order
    @test outer_dom.corners == [SVector(1.0, -1.0), SVector(1.0, 1.0), SVector(-1.0, 1.0), SVector(-1.0, -1.0)]
    @test hole_dom.corners == [SVector(0.3, 0.0)]

    # circle of radius R about the origin
    @test inner.radius == 0.3
    @test inner.center == SVector(0.0, 0.0)
end

@testset "SinaiBilliard topology" begin
    sin_b = SinaiBilliard(1.0, 0.3)
    @test genus(sin_b) == 1
    @test genus(sin_b.fundamental_domain) == 1
    # one component for the square (four curves) and one for the circle
    comps = boundary_components(sin_b)
    @test length(comps) == genus(sin_b) + 1
    @test [length(c) for c in comps] == [4, 1]
    @test [[c.domain_id for c in comp] for comp in comps] == [[1, 1, 1, 1], [2]]
    @test boundary_component_is_hole(sin_b) == [false, true]
    @test length(geometric_boundary_components(sin_b)) == genus(sin_b) + 1
    @test hole_domain_ids(sin_b) == [2]
    @test component_of_domain_id(sin_b) == Dict(1 => 0, 2 => 1)
    # only the circle is a hole curve
    @test [is_hole_curve(sin_b, c) for c in get_all_curves(sin_b)] == [false, false, false, false, true]
end

@testset "SinaiBilliard geometry" begin
    sin_b = SinaiBilliard(1.0, 0.3)
    # area 4 a^2 - pi R^2; without symmetry the fundamental area is the area
    @test area(sin_b) ≈ 4 - pi * 0.3^2 rtol = 1e-12
    @test area(SinaiBilliard(2.0, 0.5)) ≈ 16 - pi * 0.25 rtol = 1e-12
    @test fundamental_area(sin_b) ≈ area(sin_b) rtol = 1e-14
    # membership: square minus the disk, in all four quadrants
    @test is_inside(sin_b, SVector(0.5, 0.5)) === true
    @test is_inside(sin_b, SVector(-0.5, -0.5)) === true
    @test is_inside(sin_b, SVector(-0.9, 0.1)) === true
    @test is_inside(sin_b, SVector(0.31, 0.0)) === true
    @test is_inside(sin_b, SVector(0.1, 0.1)) === false     # inside the hole
    @test is_inside(sin_b, SVector(-0.1, -0.2)) === false
    @test is_inside(sin_b, SVector(0.29, 0.0)) === false
    @test is_inside(sin_b, SVector(1.1, 0.5)) === false
    # corners: the four right angles of the square, with or without the symmetry reduction
    @test corner_angles(sin_b) ≈ fill(pi / 2, 4) atol = 1e-12
    @test corner_angles(sin_b; fundamental = false) ≈ fill(pi / 2, 4) atol = 1e-12
    # full boundary: square loop then circle loop
    @test length(full_boundary(sin_b)) == 5
    @test full_boundary_components(sin_b).is_hole == [false, true]
end

@testset "SinaiBilliard keywords" begin
    # the square and the obstacle share an arbitrary center
    sin_b = SinaiBilliard(1.0, 0.3; center = [0.5, -0.2])
    dom = sin_b.fundamental_domain
    sides = dom.subdomains[1].boundary
    @test sides[1].pt0 ≈ SVector(-0.5, -1.2) atol = 1e-14
    @test sides[3].pt0 ≈ SVector(1.5, 0.8) atol = 1e-14
    inner = dom.holes[1].domain.boundary[1]
    @test inner.center ≈ SVector(0.5, -0.2) atol = 1e-14
    @test inner.radius == 0.3
    @test area(sin_b) ≈ 4 - pi * 0.09 rtol = 1e-12
    @test is_inside(sin_b, SVector(0.5, -0.2)) === false
    @test is_inside(sin_b, SVector(1.2, -0.2)) === true
    @test is_inside(sin_b, SVector(0.0, 0.0)) === true
    @test is_inside(sin_b, SVector(2.0, 0.0)) === false
    # the explicit origin is the default
    @test SinaiBilliard(1.0, 0.3; center = [0.0, 0.0]) isa SinaiBilliard{Float64}
end

@testset "SinaiBilliard argument guards" begin
    # a > 0
    @test_throws ArgumentError SinaiBilliard(0.0, 0.3)
    @test_throws ArgumentError SinaiBilliard(-1.0, 0.3)
    # R > 0
    @test_throws ArgumentError SinaiBilliard(1.0, 0.0)
    @test_throws ArgumentError SinaiBilliard(1.0, -0.1)
    # R < a
    @test_throws ArgumentError SinaiBilliard(1.0, 1.0)
    @test_throws ArgumentError SinaiBilliard(1.0, 1.5)
end

@testset "SinaiBilliard eltype promotion" begin
    sin32 = SinaiBilliard(1.0f0, 0.3f0)
    @test sin32 isa SinaiBilliard{Float32}
    @test sin32.fundamental_domain isa MultiplyConnectedDomain{Float32}
    @test sin32.fundamental_domain.subdomains[1].boundary[1] isa LineSegment{Float32}
    @test sin32.fundamental_domain.holes[1].domain.boundary[1] isa CircleSegment{Float32}
    @test eltype(sin32.fundamental_domain.subdomains[1].corners) == SVector{2,Float32}
    @test all(c -> c isa Union{LineSegment{Float32},CircleSegment{Float32}}, get_all_curves(sin32))
    @test genus(sin32) == 1
    @test area(sin32) isa Float32
end
