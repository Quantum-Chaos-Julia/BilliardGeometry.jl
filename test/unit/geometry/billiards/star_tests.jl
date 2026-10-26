# Unit tests for the star billiard constructor (structural checks only):
#   src/geometry/billiards/star.jl   (StarBilliard{T}, StarBilliard(R, a, n; center))
# Reference object: the 2 pi / n sector of r(phi) = R + a cos(n phi), a FourierCoeffPolarSegment closed by two SymmetryWalls through
# the origin; Cn registry of the n - 1 rotations by 2 pi m / n (sym_ids 1..n-1). Guards: n >= 2, R > |a|, centre at the origin.

using StaticArrays: SVector

@testset "StarBilliard fields" begin
    star = StarBilliard(1.0, 0.2, 3)
    @test star isa StarBilliard{Float64}
    @test star isa AbsBilliard
    @test star.fundamental_domain isa SimpleDomain{Float64}
    @test star.symmetries isa SymmetryRegistry
    @test fieldnames(StarBilliard) == (:fundamental_domain, :symmetries)
    @test @inferred(StarBilliard(1.0, 0.2, 3)) isa StarBilliard{Float64}
end

@testset "StarBilliard symmetries" begin
    reg = StarBilliard(1.0, 0.2, 4).symmetries
    @test length(reg) == 3
    @test all(s -> s isa NFoldRotation{Float64}, reg)
    @test [s.sym_id for s in reg] == [1, 2, 3]
    @test [s.m for s in reg] == [1, 2, 3]
    @test all(s -> s.order == 4, reg)
    # one rotation fewer than the order of the symmetry
    @test length(StarBilliard(1.0, 0.2, 2).symmetries) == 1
    @test length(StarBilliard(1.0, 0.2, 5).symmetries) == 4
end

@testset "StarBilliard boundary" begin
    n = 3
    dom = StarBilliard(1.0, 0.2, n).fundamental_domain
    @test dom.id == 1
    @test length(dom.boundary) == 3
    arc, wall1, wall0 = dom.boundary
    @test arc isa FourierCoeffPolarSegment{Float64}
    @test wall1 isa LineSegment{Float64}
    @test wall0 isa LineSegment{Float64}
    @test [c.segment_id for c in dom.boundary] == [1, 2, 3]
    @test all(c -> c.domain_id == 1, dom.boundary)
    # the walls border the copies reached by the rotations by +2 pi / n (sym_id 1) and by +2 pi (n - 1) / n (sym_id n - 1)
    @test [c.bc for c in dom.boundary] == [SpecularReflection(), SymmetryWall(1, 2), SymmetryWall(n - 1, 2)]

    # the coefficient vector [b_1, a_1, ..., b_n, a_n] has only a_n = a nonzero
    @test arc.R == 1.0
    @test length(arc.coef) == 2 * n
    @test arc.coef == [0.0, 0.0, 0.0, 0.0, 0.0, 0.2]
    @test arc.arc_angle ≈ 2 * pi / n atol = 1e-14

    # cos(n phi) = 1 at both ends of the sector, so both end points have radius R + a
    @test dom.corners[1] ≈ SVector(1.2, 0.0) atol = 1e-14
    @test dom.corners[2] ≈ 1.2 * SVector(cos(2 * pi / n), sin(2 * pi / n)) atol = 1e-14
    @test dom.corners[3] == SVector(0.0, 0.0)
    @test wall1.pt0 ≈ dom.corners[2] atol = 1e-14
    @test wall1.pt1 == dom.corners[3]
    @test wall0.pt0 == dom.corners[3]
    @test wall0.pt1 ≈ dom.corners[1] atol = 1e-14
end

@testset "StarBilliard argument guards" begin
    # n >= 2
    @test_throws ArgumentError StarBilliard(1.0, 0.2, 1)
    @test_throws ArgumentError StarBilliard(1.0, 0.2, 0)
    # R > |a| keeps r(phi) positive
    @test_throws ArgumentError StarBilliard(1.0, 1.0, 3)
    @test_throws ArgumentError StarBilliard(1.0, -1.0, 3)
    @test_throws ArgumentError StarBilliard(1.0, -1.5, 3)
    # the centre must be the origin
    @test_throws ArgumentError StarBilliard(1.0, 0.2, 3; center=[0.1, 0.0])
    @test_throws ArgumentError StarBilliard(1.0, 0.2, 3; center=[0.0, 0.1])
    @test StarBilliard(1.0, 0.2, 3; center=[0.0, 0.0]) isa StarBilliard{Float64}
end

@testset "StarBilliard eltype promotion" begin
    star32 = StarBilliard(1.0f0, 0.2f0, 4)
    @test star32 isa StarBilliard{Float32}
    @test star32.fundamental_domain isa SimpleDomain{Float32}
    @test star32.fundamental_domain.boundary[1] isa FourierCoeffPolarSegment{Float32}
    @test star32.fundamental_domain.boundary[2] isa LineSegment{Float32}
    @test eltype(star32.fundamental_domain.corners) == SVector{2,Float32}
    @test all(s -> s isa NFoldRotation{Float32}, star32.symmetries)
end
