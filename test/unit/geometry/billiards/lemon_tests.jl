# Unit tests for the lemon billiard constructor (structural checks only):
#   src/geometry/billiards/lemon.jl   (LemonBilliard{T}, LemonBilliard(half_width))
# Reference object: intersection of two unit disks centred at (0, ±a). Fundamental domain is the first quadrant: a CircleSegment arc
# centred at (0, -a) from the corner (sqrt(1 - a^2), 0) to (0, 1 - a), closed by a y-axis wall (SymmetryWall(1, 2)) and an x-axis wall
# (SymmetryWall(3, 2)); D2 symmetry registry. Guard: 0 < half_width < 1.

using StaticArrays: SVector

@testset "LemonBilliard fields" begin
    b = LemonBilliard(0.5)
    @test b isa LemonBilliard{Float64}
    @test b isa AbsBilliard
    @test b.fundamental_domain isa SimpleDomain{Float64}
    @test b.symmetries isa SymmetryRegistry
    @test fieldnames(LemonBilliard) == (:fundamental_domain, :symmetries)
    @test @inferred(LemonBilliard(0.5)) isa LemonBilliard{Float64}
    @test @inferred(LemonBilliard(0.5f0)) isa LemonBilliard{Float32}
end

@testset "LemonBilliard symmetries" begin
    reg = LemonBilliard(0.5).symmetries
    @test length(reg) == 3
    @test reg[1] isa YAxisReflection
    @test reg[2] isa XYAxisReflection
    @test reg[3] isa XAxisReflection
    @test [s.sym_id for s in reg] == [1, 2, 3]
    @test reg == D2_symmetry()
end

@testset "LemonBilliard boundary" begin
    a = 0.5
    h = sqrt(1 - a^2)
    dom = LemonBilliard(a).fundamental_domain
    @test dom.id == 1
    @test length(dom.boundary) == 3
    arc, ywall, xwall = dom.boundary
    @test arc isa CircleSegment{Float64}
    @test ywall isa LineSegment{Float64}
    @test xwall isa LineSegment{Float64}
    @test [c.segment_id for c in dom.boundary] == [1, 2, 3]
    @test all(c -> c.domain_id == 1, dom.boundary)
    @test [c.bc for c in dom.boundary] == [SpecularReflection(), SymmetryWall(1, 2), SymmetryWall(3, 2)]

    # vertices
    @test length(dom.corners) == 3
    @test dom.corners[1] ≈ SVector(h, 0.0) atol = 1e-14
    @test dom.corners[2] ≈ SVector(0.0, 1 - a) atol = 1e-14
    @test dom.corners[3] == SVector(0.0, 0.0)

    # arc: unit circle about (0, -a), endpoints coincide with the vertices
    @test arc.radius == 1.0
    @test arc.center == SVector(0.0, -a)
    @test curve(arc, 0.0) ≈ dom.corners[1] atol = 1e-12
    @test curve(arc, 1.0) ≈ dom.corners[2] atol = 1e-12
    for t in range(0.0, 1.0; length=11)
        p = curve(arc, t)
        @test hypot(p[1], p[2] + a) ≈ 1.0 atol = 1e-12
        @test p[1] >= -1e-12 && p[2] >= -1e-12
    end

    # walls connect: corner -> origin -> vertex, and arc end meets the y-wall start
    @test ywall.pt0 ≈ dom.corners[2] atol = 1e-14
    @test ywall.pt1 == SVector(0.0, 0.0)
    @test xwall.pt0 == SVector(0.0, 0.0)
    @test xwall.pt1 ≈ dom.corners[1] atol = 1e-14
    @test curve(arc, 1.0) ≈ ywall.pt0 atol = 1e-12
    @test curve(arc, 0.0) ≈ xwall.pt1 atol = 1e-12
end

@testset "LemonBilliard argument guards" begin
    @test_throws ArgumentError LemonBilliard(0.0)
    @test_throws ArgumentError LemonBilliard(-0.5)
    @test_throws ArgumentError LemonBilliard(1.0)
    @test_throws ArgumentError LemonBilliard(1.5)
    @test_throws ArgumentError LemonBilliard(1.0f0)
end

@testset "LemonBilliard eltype promotion" begin
    b32 = LemonBilliard(0.5f0)
    @test b32.fundamental_domain isa SimpleDomain{Float32}
    @test b32.fundamental_domain.boundary[1] isa CircleSegment{Float32}
    @test b32.fundamental_domain.boundary[2] isa LineSegment{Float32}
    @test b32.fundamental_domain.boundary[3] isa LineSegment{Float32}
    @test eltype(b32.fundamental_domain.corners) == SVector{2,Float32}
end
