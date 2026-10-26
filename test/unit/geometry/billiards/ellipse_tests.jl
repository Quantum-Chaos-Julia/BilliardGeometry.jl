# Unit tests for the ellipse billiard constructor (structural checks only):
#   src/geometry/billiards/ellipse.jl   (EllipseBilliard{T}, EllipseBilliard(a, b; center))
# Reference object: the quarter of x^2/a^2 + y^2/b^2 = 1 as a PolarSegment arc closed by a y-axis wall (SymmetryWall(1, 2)) and an
# x-axis wall (SymmetryWall(3, 2)); D2 symmetry registry. Guards: a > 0, b > 0, centre at the origin.

using StaticArrays: SVector

@testset "EllipseBilliard fields" begin
    ell = EllipseBilliard(2.0, 1.0)
    @test ell isa EllipseBilliard{Float64}
    @test ell isa AbsBilliard
    @test ell.fundamental_domain isa SimpleDomain{Float64}
    @test ell.symmetries isa SymmetryRegistry
    @test fieldnames(EllipseBilliard) == (:fundamental_domain, :symmetries)
    @test @inferred(EllipseBilliard(2.0, 1.0)) isa EllipseBilliard{Float64}
end

@testset "EllipseBilliard symmetries" begin
    reg = EllipseBilliard(2.0, 1.0).symmetries
    @test length(reg) == 3
    @test reg[1] isa YAxisReflection
    @test reg[2] isa XYAxisReflection
    @test reg[3] isa XAxisReflection
    @test [s.sym_id for s in reg] == [1, 2, 3]
    @test reg == D2_symmetry()
end

@testset "EllipseBilliard boundary" begin
    dom = EllipseBilliard(2.0, 1.0).fundamental_domain
    @test dom.id == 1
    @test length(dom.boundary) == 3
    arc, ywall, xwall = dom.boundary
    @test arc isa PolarSegment{Float64}
    @test ywall isa LineSegment{Float64}
    @test xwall isa LineSegment{Float64}
    @test [c.segment_id for c in dom.boundary] == [1, 2, 3]
    @test all(c -> c.domain_id == 1, dom.boundary)
    @test [c.bc for c in dom.boundary] == [SpecularReflection(), SymmetryWall(1, 2), SymmetryWall(3, 2)]

    # the quarter arc runs from the vertex (a, 0) to the co-vertex (0, b); R is the larger semi-axis
    @test arc.arc_angle == pi / 2
    @test arc.shift_angle == 0.0
    @test arc.R == 2.0
    @test arc.center == SVector(0.0, 0.0)
    @test dom.corners[1] ≈ SVector(2.0, 0.0) atol = 1e-14
    @test dom.corners[2] ≈ SVector(0.0, 1.0) atol = 1e-14
    @test dom.corners[3] == SVector(0.0, 0.0)
    @test ywall.pt1 == SVector(0.0, 0.0)
    @test xwall.pt0 == SVector(0.0, 0.0)
    # R follows the larger semi-axis
    @test EllipseBilliard(1.0, 3.0).fundamental_domain.boundary[1].R == 3.0
end

@testset "EllipseBilliard argument guards" begin
    @test_throws ArgumentError EllipseBilliard(0.0, 1.0)
    @test_throws ArgumentError EllipseBilliard(-1.0, 1.0)
    @test_throws ArgumentError EllipseBilliard(1.0, 0.0)
    @test_throws ArgumentError EllipseBilliard(1.0, -2.0)
    @test_throws ArgumentError EllipseBilliard(2.0, 1.0; center=[0.1, 0.0])
    @test_throws ArgumentError EllipseBilliard(2.0, 1.0; center=[0.0, 0.1])
end

@testset "EllipseBilliard eltype promotion" begin
    ell32 = EllipseBilliard(2.0f0, 1.0f0)
    @test ell32 isa EllipseBilliard{Float32}
    @test ell32.fundamental_domain isa SimpleDomain{Float32}
    @test ell32.fundamental_domain.boundary[1] isa PolarSegment{Float32}
    @test ell32.fundamental_domain.boundary[2] isa LineSegment{Float32}
    @test eltype(ell32.fundamental_domain.corners) == SVector{2,Float32}
end
