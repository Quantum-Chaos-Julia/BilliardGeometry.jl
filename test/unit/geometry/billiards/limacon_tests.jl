# Unit tests for the limacon billiard constructor (structural checks only):
#   src/geometry/billiards/limacon.jl   (LimaconBilliard{T}, LimaconBilliard(a; R))
# Reference object: the upper half of r(phi) = 1 + a cos(phi) (radius R = 1), a polar arc over [0, pi] with id 1 closed by a
# SymmetryWall on the x axis (sym_id 1); one XAxisReflection symmetry. Guard: |a| <= 1 keeps the curve simple.

using StaticArrays: SVector

@testset "LimaconBilliard fields" begin
    lim = LimaconBilliard(0.5)
    @test lim isa LimaconBilliard{Float64}
    @test lim isa AbsBilliard
    @test lim.fundamental_domain isa SimpleDomain{Float64}
    @test lim.symmetries isa SymmetryRegistry
    @test fieldnames(LimaconBilliard) == (:fundamental_domain, :symmetries)
    @test @inferred(LimaconBilliard(0.5)) isa LimaconBilliard{Float64}
end

@testset "LimaconBilliard symmetries" begin
    reg = LimaconBilliard(0.5).symmetries
    @test length(reg) == 1
    @test reg[1] isa XAxisReflection
    @test reg[1].sym_id == 1
end

@testset "LimaconBilliard boundary" begin
    dom = LimaconBilliard(0.5).fundamental_domain
    @test dom.id == 1
    @test length(dom.boundary) == 2
    arc, wall = dom.boundary
    @test arc isa FourierCoeffPolarSegment{Float64}
    @test wall isa LineSegment{Float64}
    @test [c.segment_id for c in dom.boundary] == [1, 2]
    @test all(c -> c.domain_id == 1, dom.boundary)
    @test arc.bc isa SpecularReflection
    @test wall.bc == SymmetryWall(1, 2)

    # the arc covers phi in [0, pi] of r = R + a cos(phi): coefficients [b_1, a_1] = [0, a]
    @test arc.R == 1.0
    @test arc.coef == [0.0, 0.5]
    @test arc.arc_angle == Float64(pi)
    @test arc.shift_angle == 0.0

    # the end points on the x axis are (1 + a, 0) and (-(1 - a), 0); the wall runs from the second to the first
    @test dom.corners[1] ≈ SVector(1.5, 0.0) atol = 1e-14
    @test dom.corners[2] ≈ SVector(-0.5, 0.0) atol = 1e-14
    @test wall.pt0 ≈ dom.corners[2] atol = 1e-14
    @test wall.pt1 ≈ dom.corners[1] atol = 1e-14
end

@testset "LimaconBilliard radius keyword" begin
    dom = LimaconBilliard(0.5; R=2.0).fundamental_domain
    @test dom.boundary[1].R == 2.0
    @test dom.boundary[1].coef == [0.0, 0.5]
end

@testset "LimaconBilliard argument guards" begin
    # |a| <= 1 is required; the limiting values are accepted
    @test_throws ArgumentError LimaconBilliard(1.1)
    @test_throws ArgumentError LimaconBilliard(-1.5)
    @test LimaconBilliard(1.0) isa LimaconBilliard{Float64}
    @test LimaconBilliard(-1.0) isa LimaconBilliard{Float64}
end

@testset "LimaconBilliard eltype promotion" begin
    lim32 = LimaconBilliard(0.5f0)
    @test lim32 isa LimaconBilliard{Float32}
    @test lim32.fundamental_domain isa SimpleDomain{Float32}
    @test lim32.fundamental_domain.boundary[1] isa FourierCoeffPolarSegment{Float32}
    @test lim32.fundamental_domain.boundary[2] isa LineSegment{Float32}
    @test eltype(lim32.fundamental_domain.corners) == SVector{2,Float32}
    @test LimaconBilliard(0.5f0; R=2.0f0) isa LimaconBilliard{Float32}
end
