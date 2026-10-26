# Unit tests for the circle billiard constructor (structural checks only):
#   src/geometry/billiards/circle.jl   (CircleBilliard{T}, CircleBilliard(R; center))
# Reference object: the quarter disc of radius R with D2 symmetry (YAxisReflection, XYAxisReflection, XAxisReflection with sym_ids
# 1, 2, 3), a CircleSegment arc closed by a y-axis wall (SymmetryWall(1, 2)) and an x-axis wall (SymmetryWall(3, 2)).
# Guard: the centre must be the origin because the D2 registry reflects about the coordinate axes.

using StaticArrays: SVector

@testset "CircleBilliard fields" begin
    circ = CircleBilliard(2.0)
    @test circ isa CircleBilliard{Float64}
    @test circ isa AbsBilliard
    @test circ.fundamental_domain isa SimpleDomain{Float64}
    @test circ.symmetries isa SymmetryRegistry
    @test fieldnames(CircleBilliard) == (:fundamental_domain, :symmetries)
    @test @inferred(CircleBilliard(2.0)) isa CircleBilliard{Float64}
    # the radius defaults to one
    @test CircleBilliard().fundamental_domain.boundary[1].radius == 1.0
end

@testset "CircleBilliard symmetries" begin
    reg = CircleBilliard(2.0).symmetries
    @test length(reg) == 3
    @test reg[1] isa YAxisReflection
    @test reg[2] isa XYAxisReflection
    @test reg[3] isa XAxisReflection
    @test [s.sym_id for s in reg] == [1, 2, 3]
    @test reg == D2_symmetry()
end

@testset "CircleBilliard boundary" begin
    dom = CircleBilliard(2.0).fundamental_domain
    @test dom.id == 1
    @test length(dom.boundary) == 3
    arc, ywall, xwall = dom.boundary
    @test arc isa CircleSegment{Float64}
    @test ywall isa LineSegment{Float64}
    @test xwall isa LineSegment{Float64}
    @test [c.segment_id for c in dom.boundary] == [1, 2, 3]
    @test all(c -> c.domain_id == 1, dom.boundary)
    @test [c.bc for c in dom.boundary] == [SpecularReflection(), SymmetryWall(1, 2), SymmetryWall(3, 2)]

    # quarter arc of radius R about the origin from (R, 0) to (0, R), closed through the origin
    @test arc.radius == 2.0
    @test arc.arc_angle == pi / 2
    @test arc.shift_angle == 0.0
    @test arc.center == SVector(0.0, 0.0)
    @test dom.corners == [SVector(2.0, 0.0), SVector(0.0, 2.0), SVector(0.0, 0.0)]
    @test ywall.pt0 == SVector(0.0, 2.0)
    @test ywall.pt1 == SVector(0.0, 0.0)
    @test xwall.pt0 == SVector(0.0, 0.0)
    @test xwall.pt1 == SVector(2.0, 0.0)
end

@testset "CircleBilliard center guard" begin
    @test_throws ArgumentError CircleBilliard(1.0; center=[0.1, 0.0])
    @test_throws ArgumentError CircleBilliard(1.0; center=[0.0, 0.1])
    # an explicit origin is accepted
    @test CircleBilliard(1.0; center=[0.0, 0.0]) isa CircleBilliard{Float64}
end

@testset "CircleBilliard eltype promotion" begin
    circ32 = CircleBilliard(2.0f0)
    @test circ32 isa CircleBilliard{Float32}
    @test circ32.fundamental_domain isa SimpleDomain{Float32}
    @test circ32.fundamental_domain.boundary[1] isa CircleSegment{Float32}
    @test circ32.fundamental_domain.boundary[2] isa LineSegment{Float32}
    @test eltype(circ32.fundamental_domain.corners) == SVector{2,Float32}
    @test circ32.symmetries[1] isa YAxisReflection
    # a Float64 centre is converted to the element type of the radius
    @test CircleBilliard(2.0f0; center=[0.0, 0.0]) isa CircleBilliard{Float32}
end
