# Unit tests for the Prosen billiard constructor (structural checks only):
#   src/geometry/billiards/prosen.jl   (ProsenBilliard{T}, ProsenBilliard(a; center))
# Reference object: the quadrant of r(phi) = 1 + a cos(4 phi), a FourierCoeffPolarSegment over [0, pi/2] closed by a y-axis wall
# (SymmetryWall(1, 2)) and an x-axis wall (SymmetryWall(3, 2)); D2 symmetry registry. Guard: centre at the origin.

using StaticArrays: SVector

@testset "ProsenBilliard fields" begin
    pro = ProsenBilliard(0.2)
    @test pro isa ProsenBilliard{Float64}
    @test pro isa AbsBilliard
    @test pro.fundamental_domain isa SimpleDomain{Float64}
    @test pro.symmetries isa SymmetryRegistry
    @test fieldnames(ProsenBilliard) == (:fundamental_domain, :symmetries)
    @test @inferred(ProsenBilliard(0.2)) isa ProsenBilliard{Float64}
end

@testset "ProsenBilliard symmetries" begin
    reg = ProsenBilliard(0.2).symmetries
    @test length(reg) == 3
    @test reg[1] isa YAxisReflection
    @test reg[2] isa XYAxisReflection
    @test reg[3] isa XAxisReflection
    @test [s.sym_id for s in reg] == [1, 2, 3]
    @test reg == D2_symmetry()
end

@testset "ProsenBilliard boundary" begin
    a = 0.2
    dom = ProsenBilliard(a).fundamental_domain
    @test dom.id == 1
    @test length(dom.boundary) == 3
    arc, ywall, xwall = dom.boundary
    @test arc isa FourierCoeffPolarSegment{Float64}
    @test ywall isa LineSegment{Float64}
    @test xwall isa LineSegment{Float64}
    # the boundary is listed arc, y-axis wall, x-axis wall; the x-axis wall carries segment id 2 and the y-axis wall id 3
    @test [c.segment_id for c in dom.boundary] == [1, 3, 2]
    @test all(c -> c.domain_id == 1, dom.boundary)
    @test [c.bc for c in dom.boundary] == [SpecularReflection(), SymmetryWall(1, 2), SymmetryWall(3, 2)]

    # a_4 = a is the only nonzero coefficient, R = 1, the arc spans a quarter turn
    @test arc.R == 1.0
    @test length(arc.coef) == 8
    @test arc.coef == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, a]
    @test arc.arc_angle ≈ pi / 2 atol = 1e-14

    # cos(4 phi) = 1 at phi = 0 and pi/2, so both end points have radius 1 + a
    @test dom.corners[1] ≈ SVector(1 + a, 0.0) atol = 1e-14
    @test dom.corners[2] ≈ SVector(0.0, 1 + a) atol = 1e-14
    @test dom.corners[3] == SVector(0.0, 0.0)
    @test ywall.pt0 ≈ dom.corners[2] atol = 1e-14
    @test ywall.pt1 == dom.corners[3]
    @test xwall.pt0 == dom.corners[3]
    @test xwall.pt1 ≈ dom.corners[1] atol = 1e-14
end

@testset "ProsenBilliard center guard" begin
    @test_throws ArgumentError ProsenBilliard(0.2; center=[0.1, 0.0])
    @test_throws ArgumentError ProsenBilliard(0.2; center=[0.0, 0.1])
    @test ProsenBilliard(0.2; center=[0.0, 0.0]) isa ProsenBilliard{Float64}
end

@testset "ProsenBilliard eltype promotion" begin
    pro32 = ProsenBilliard(0.2f0)
    @test pro32 isa ProsenBilliard{Float32}
    @test pro32.fundamental_domain isa SimpleDomain{Float32}
    @test pro32.fundamental_domain.boundary[1] isa FourierCoeffPolarSegment{Float32}
    @test pro32.fundamental_domain.boundary[2] isa LineSegment{Float32}
    @test eltype(pro32.fundamental_domain.corners) == SVector{2,Float32}
end
