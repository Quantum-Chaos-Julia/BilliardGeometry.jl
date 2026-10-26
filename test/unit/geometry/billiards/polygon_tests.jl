# Unit tests for the polygon billiard constructor (structural checks only):
#   src/geometry/billiards/polygon.jl   (PolygonBilliard{T}, PolygonBilliard(vertices; symmetries))
# Reference object: the unit square given by four counterclockwise vertices, one LineSegment per edge with id 1..n, no symmetries
# unless a registry is passed. Guards: at least three vertices and counterclockwise orientation (positive signed area).

using StaticArrays: SVector

const pb_SQUARE = SVector{2,Float64}[SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(1.0, 1.0), SVector(0.0, 1.0)]

@testset "PolygonBilliard fields" begin
    poly = PolygonBilliard(pb_SQUARE)
    @test poly isa PolygonBilliard{Float64}
    @test poly isa AbsBilliard
    @test poly.fundamental_domain isa SimpleDomain{Float64}
    @test poly.symmetries isa SymmetryRegistry
    @test fieldnames(PolygonBilliard) == (:fundamental_domain, :symmetries)
end

@testset "PolygonBilliard symmetries" begin
    @test isempty(PolygonBilliard(pb_SQUARE).symmetries)
    # a registry passed as keyword is stored unchanged
    reg = register_symmetries(DiagonalReflection())
    @test PolygonBilliard(pb_SQUARE; symmetries=reg).symmetries === reg
    @test PolygonBilliard(pb_SQUARE; symmetries=reg).symmetries[1].sym_id == 1
end

@testset "PolygonBilliard boundary" begin
    dom = PolygonBilliard(pb_SQUARE).fundamental_domain
    @test dom.id == 1
    @test dom.corners == pb_SQUARE
    @test length(dom.boundary) == 4
    @test all(c -> c isa LineSegment{Float64}, dom.boundary)
    @test [c.segment_id for c in dom.boundary] == [1, 2, 3, 4]
    @test all(c -> c.domain_id == 1, dom.boundary)
    @test all(c -> c.bc isa SpecularReflection, dom.boundary)
    # edge i joins vertex i to vertex i + 1, the last edge closes the polygon
    @test [c.pt0 for c in dom.boundary] == pb_SQUARE
    @test [c.pt1 for c in dom.boundary] == pb_SQUARE[[2, 3, 4, 1]]

    # a triangle has three edges
    tri = PolygonBilliard(pb_SQUARE[1:3]).fundamental_domain
    @test length(tri.boundary) == 3
    @test [c.segment_id for c in tri.boundary] == [1, 2, 3]
end

@testset "PolygonBilliard argument guards" begin
    # at least three vertices
    @test_throws ArgumentError PolygonBilliard(SVector{2,Float64}[])
    @test_throws ArgumentError PolygonBilliard(pb_SQUARE[1:2])
    # clockwise order gives a negative signed area
    @test_throws ArgumentError PolygonBilliard(reverse(pb_SQUARE))
    # collinear vertices give zero signed area
    @test_throws ArgumentError PolygonBilliard(SVector{2,Float64}[SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(2.0, 0.0)])
end

@testset "PolygonBilliard eltype promotion" begin
    poly32 = PolygonBilliard(SVector{2,Float32}.(pb_SQUARE))
    @test poly32 isa PolygonBilliard{Float32}
    @test poly32.fundamental_domain isa SimpleDomain{Float32}
    @test all(c -> c isa LineSegment{Float32}, poly32.fundamental_domain.boundary)
    @test eltype(poly32.fundamental_domain.corners) == SVector{2,Float32}
end
