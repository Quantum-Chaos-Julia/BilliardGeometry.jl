# Unit tests for the rectangle billiard constructor (structural checks only):
#   src/geometry/billiards/rectangle.jl   (RectangleBilliard{T}, RectangleBilliard(a, b; center))
# Reference object: the first quadrant [0, a] x [0, b] of the 2a x 2b rectangle, four LineSegments (two physical walls, a y-axis wall
# SymmetryWall(1, 2) and an x-axis wall SymmetryWall(3, 2)); D2 symmetry registry. Guard: centre at the origin.

using StaticArrays: SVector

@testset "RectangleBilliard fields" begin
    rect = RectangleBilliard(2.0, 1.0)
    @test rect isa RectangleBilliard{Float64}
    @test rect isa AbsBilliard
    @test rect.fundamental_domain isa SimpleDomain{Float64}
    @test rect.symmetries isa SymmetryRegistry
    @test fieldnames(RectangleBilliard) == (:fundamental_domain, :symmetries)
    @test @inferred(RectangleBilliard(2.0, 1.0)) isa RectangleBilliard{Float64}
    # both half-sides default to one
    @test RectangleBilliard().fundamental_domain.corners[2] == SVector(1.0, 1.0)
end

@testset "RectangleBilliard symmetries" begin
    reg = RectangleBilliard(2.0, 1.0).symmetries
    @test length(reg) == 3
    @test reg[1] isa YAxisReflection
    @test reg[2] isa XYAxisReflection
    @test reg[3] isa XAxisReflection
    @test [s.sym_id for s in reg] == [1, 2, 3]
    @test reg == D2_symmetry()
end

@testset "RectangleBilliard boundary" begin
    dom = RectangleBilliard(2.0, 1.0).fundamental_domain
    @test dom.id == 1
    @test length(dom.boundary) == 4
    @test all(c -> c isa LineSegment{Float64}, dom.boundary)
    @test [c.segment_id for c in dom.boundary] == [1, 2, 3, 4]
    @test all(c -> c.domain_id == 1, dom.boundary)
    @test [c.bc for c in dom.boundary] == [SpecularReflection(), SpecularReflection(), SymmetryWall(1, 2), SymmetryWall(3, 2)]

    # the walls chain counterclockwise through the corners (a, 0), (a, b), (0, b), (0, 0)
    @test dom.corners == [SVector(2.0, 0.0), SVector(2.0, 1.0), SVector(0.0, 1.0), SVector(0.0, 0.0)]
    @test [c.pt0 for c in dom.boundary] == dom.corners
    @test [c.pt1 for c in dom.boundary] == [dom.corners[2], dom.corners[3], dom.corners[4], dom.corners[1]]
end

@testset "RectangleBilliard center guard" begin
    @test_throws ArgumentError RectangleBilliard(2.0, 1.0; center=[0.1, 0.0])
    @test_throws ArgumentError RectangleBilliard(2.0, 1.0; center=[0.0, 0.1])
    @test RectangleBilliard(2.0, 1.0; center=[0.0, 0.0]) isa RectangleBilliard{Float64}
end

@testset "RectangleBilliard eltype promotion" begin
    rect32 = RectangleBilliard(2.0f0, 1.0f0)
    @test rect32 isa RectangleBilliard{Float32}
    @test rect32.fundamental_domain isa SimpleDomain{Float32}
    @test all(c -> c isa LineSegment{Float32}, rect32.fundamental_domain.boundary)
    @test eltype(rect32.fundamental_domain.corners) == SVector{2,Float32}
end
