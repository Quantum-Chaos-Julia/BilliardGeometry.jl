# Unit tests for the polar billiard constructor (structural checks only):
#   src/geometry/billiards/polar.jl   (PolarDomain{T}, PolarBilliard{T}, PolarBilliard(coef; center))
# Reference object: the closed curve r(phi) = 1 + sum_n (b_n sin(n phi) + a_n cos(n phi)) given by Fourier coefficients
# [b_1, a_1, b_2, a_2, ...], one FourierCoeffPolarSegment, no symmetries. The constructor has no argument guards.

using StaticArrays: SVector

@testset "PolarBilliard fields" begin
    pol = PolarBilliard([0.0, 0.1])
    @test pol isa PolarBilliard{Float64}
    @test pol isa AbsBilliard
    @test pol.fundamental_domain isa PolarDomain{Float64}
    @test pol.fundamental_domain isa AbsPolarDomain
    @test pol.symmetries isa SymmetryRegistry
    @test fieldnames(PolarBilliard) == (:fundamental_domain, :symmetries)
    @test fieldnames(PolarDomain) == (:boundary, :corners, :id)
    @test @inferred(PolarBilliard([0.0, 0.1])) isa PolarBilliard{Float64}
end

@testset "PolarBilliard symmetries" begin
    @test isempty(PolarBilliard([0.0, 0.1]).symmetries)
end

@testset "PolarBilliard boundary" begin
    dom = PolarBilliard([0.0, 0.1]).fundamental_domain
    @test dom.id == 1
    @test length(dom.boundary) == 1
    seg = dom.boundary[1]
    @test seg isa FourierCoeffPolarSegment{Float64}
    @test seg.coef == [0.0, 0.1]
    @test seg.arc_angle == 2 * pi
    @test seg.center == SVector(0.0, 0.0)
    @test seg.bc isa SpecularReflection
    # the single corner is the point of the curve at phi = 0, r(0) = R + a_1
    @test dom.corners == [SVector(1.1, 0.0)]
end

@testset "PolarBilliard center keyword" begin
    dom = PolarBilliard([0.0, 0.1]; center=[0.3, -0.2]).fundamental_domain
    @test dom.boundary[1].center == SVector(0.3, -0.2)
    @test dom.corners[1] ≈ SVector(1.4, -0.2) atol = 1e-14
end

@testset "PolarBilliard eltype promotion" begin
    pol32 = PolarBilliard(Float32[0.0, 0.1])
    @test pol32 isa PolarBilliard{Float32}
    @test pol32.fundamental_domain isa PolarDomain{Float32}
    @test eltype(pol32.fundamental_domain.corners) == SVector{2,Float32}
end
