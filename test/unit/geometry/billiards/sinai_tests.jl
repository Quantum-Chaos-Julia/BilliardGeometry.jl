# Unit tests for the Sinai billiard constructor (structural checks only):
#   src/geometry/billiards/sinai.jl   (SinaiBilliard{T}, SinaiBilliard(a, R; center))
# Reference object: the quadrant [0, a]^2 with a quarter-disc hole of radius R about the origin, a SimpleDomain of five curves
# (two physical walls, y-axis wall, hole arc with orientation -1 and domain id 2, x-axis wall); D2 symmetry registry.
# Only structure is checked: curve count and types, orientations, ids, genus and components. Guards: a > 0, R > 0, R < a, centre at origin.

using StaticArrays: SVector

@testset "SinaiBilliard fields" begin
    sin_b = SinaiBilliard(1.0, 0.3)
    @test sin_b isa SinaiBilliard{Float64}
    @test sin_b isa AbsBilliard
    @test sin_b.fundamental_domain isa SimpleDomain{Float64}
    @test sin_b.symmetries isa SymmetryRegistry
    @test fieldnames(SinaiBilliard) == (:fundamental_domain, :symmetries)
    @test @inferred(SinaiBilliard(1.0, 0.3)) isa SinaiBilliard{Float64}
end

@testset "SinaiBilliard symmetries" begin
    reg = SinaiBilliard(1.0, 0.3).symmetries
    @test length(reg) == 3
    @test reg[1] isa YAxisReflection
    @test reg[2] isa XYAxisReflection
    @test reg[3] isa XAxisReflection
    @test [s.sym_id for s in reg] == [1, 2, 3]
    @test reg == D2_symmetry()
end

@testset "SinaiBilliard curves" begin
    dom = SinaiBilliard(1.0, 0.3).fundamental_domain
    @test dom.id == 1
    @test length(dom.boundary) == 5
    @test length(dom.corners) == 5
    @test [typeof(c) for c in dom.boundary] == [
        LineSegment{Float64,SpecularReflection}, LineSegment{Float64,SpecularReflection}, LineSegment{Float64,SymmetryWall},
        CircleSegment{Float64,SpecularReflection}, LineSegment{Float64,SymmetryWall},
    ]
    @test [c.segment_id for c in dom.boundary] == [1, 2, 3, 1, 4]
    @test [c.bc for c in dom.boundary] == [SpecularReflection(), SpecularReflection(), SymmetryWall(1, 2), SpecularReflection(), SymmetryWall(3, 2)]

    # the hole arc is the only curve with orientation -1 and carries its own domain id
    @test [c.orientation for c in dom.boundary] == [1, 1, 1, -1, 1]
    @test [c.domain_id for c in dom.boundary] == [1, 1, 1, 2, 1]
    hole = dom.boundary[4]
    @test hole.radius == 0.3
    @test hole.arc_angle == pi / 2
    @test hole.center == SVector(0.0, 0.0)

    # outer corners of the quadrant
    @test dom.corners[1:3] == [SVector(1.0, 0.0), SVector(1.0, 1.0), SVector(0.0, 1.0)]
end

@testset "SinaiBilliard topology" begin
    sin_b = SinaiBilliard(1.0, 0.3)
    # the hole is part of the flat boundary of a SimpleDomain, so the domain itself has genus zero and one component
    @test genus(sin_b) == 0
    @test length(boundary_components(sin_b)) == 1
    @test length(geometric_boundary_components(sin_b)) == 1
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
    # the centre must be the origin
    @test_throws ArgumentError SinaiBilliard(1.0, 0.3; center=[0.1, 0.0])
    @test_throws ArgumentError SinaiBilliard(1.0, 0.3; center=[0.0, 0.1])
    @test SinaiBilliard(1.0, 0.3; center=[0.0, 0.0]) isa SinaiBilliard{Float64}
end

@testset "SinaiBilliard eltype promotion" begin
    sin32 = SinaiBilliard(1.0f0, 0.3f0)
    @test sin32 isa SinaiBilliard{Float32}
    @test sin32.fundamental_domain isa SimpleDomain{Float32}
    @test sin32.fundamental_domain.boundary[1] isa LineSegment{Float32}
    @test sin32.fundamental_domain.boundary[4] isa CircleSegment{Float32}
    @test eltype(sin32.fundamental_domain.corners) == SVector{2,Float32}
end
