# Unit tests for the C3 billiard constructor (structural checks only):
#   src/geometry/billiards/c3.jl   (C3Billiard{T}, C3Billiard(a; scale))
# Reference object: the 2 pi / 3 sector of r(phi) = scale/2 (1 + a (cos(3 phi) - sin(6 phi))), a FourierCoeffPolarSegment closed by
# two SymmetryWalls through the origin; C3 registry (rotations by 2 pi / 3 and 4 pi / 3, sym_ids 1, 2). Guard: scale > 0.

using StaticArrays: SVector

@testset "C3Billiard fields" begin
    c3 = C3Billiard(0.2)
    @test c3 isa C3Billiard{Float64}
    @test c3 isa AbsBilliard
    @test c3.fundamental_domain isa SimpleDomain{Float64}
    @test c3.symmetries isa SymmetryRegistry
    @test fieldnames(C3Billiard) == (:fundamental_domain, :symmetries)
    @test @inferred(C3Billiard(0.2)) isa C3Billiard{Float64}
end

@testset "C3Billiard symmetries" begin
    reg = C3Billiard(0.2).symmetries
    @test length(reg) == 2
    @test all(s -> s isa NFoldRotation{Float64}, reg)
    @test [s.sym_id for s in reg] == [1, 2]
    @test [s.m for s in reg] == [1, 2]
    @test all(s -> s.order == 3, reg)
end

@testset "C3Billiard boundary" begin
    a = 0.2
    dom = C3Billiard(a).fundamental_domain
    @test dom.id == 1
    @test length(dom.boundary) == 3
    arc, wall1, wall0 = dom.boundary
    @test arc isa FourierCoeffPolarSegment{Float64}
    @test wall1 isa LineSegment{Float64}
    @test wall0 isa LineSegment{Float64}
    @test [c.segment_id for c in dom.boundary] == [1, 2, 3]
    @test all(c -> c.domain_id == 1, dom.boundary)
    @test [c.bc for c in dom.boundary] == [SpecularReflection(), SymmetryWall(1, 2), SymmetryWall(2, 2)]

    # R = scale/2; a_3 = a scale/2 multiplies cos(3 phi) and b_6 = -a scale/2 multiplies sin(6 phi)
    @test arc.R == 0.5
    @test length(arc.coef) == 12
    @test arc.coef[6] == a / 2
    @test arc.coef[11] == -a / 2
    @test count(!iszero, arc.coef) == 2
    @test arc.arc_angle ≈ 2 * pi / 3 atol = 1e-14

    # cos(3 phi) = 1 and sin(6 phi) = 0 at phi = 0 and 2 pi / 3, so both end points have radius R (1 + a)
    r = 0.5 * (1 + a)
    @test dom.corners[1] ≈ SVector(r, 0.0) atol = 1e-14
    @test dom.corners[2] ≈ r * SVector(cos(2 * pi / 3), sin(2 * pi / 3)) atol = 1e-14
    @test dom.corners[3] == SVector(0.0, 0.0)
    @test wall1.pt1 == dom.corners[3]
    @test wall0.pt0 == dom.corners[3]
end

@testset "C3Billiard scale keyword" begin
    dom = C3Billiard(0.2; scale=2.0).fundamental_domain
    @test dom.boundary[1].R == 1.0
    @test dom.corners[1] ≈ SVector(1.2, 0.0) atol = 1e-14
end

@testset "C3Billiard argument guards" begin
    @test_throws ArgumentError C3Billiard(0.2; scale=0.0)
    @test_throws ArgumentError C3Billiard(0.2; scale=-1.0)
end

@testset "C3Billiard eltype promotion" begin
    c32 = C3Billiard(0.2f0)
    @test c32 isa C3Billiard{Float32}
    @test c32.fundamental_domain isa SimpleDomain{Float32}
    @test c32.fundamental_domain.boundary[1] isa FourierCoeffPolarSegment{Float32}
    @test c32.fundamental_domain.boundary[2] isa LineSegment{Float32}
    @test eltype(c32.fundamental_domain.corners) == SVector{2,Float32}
    @test all(s -> s isa NFoldRotation{Float32}, c32.symmetries)
    @test C3Billiard(0.2f0; scale=2.0f0) isa C3Billiard{Float32}
end
