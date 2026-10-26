# Unit tests for src/geometry/domains/polygons.jl: the Polygon domain type and its constructor.
# Reference polygons are the unit square, right triangle, equilateral triangle and regular pentagon, whose interior angles are known; helpers use the prefix `pg_`.

using StaticArrays: SVector

# Absolute tolerance for interior angles in radians (worst observed 2.2e-16, margin 4500x).
const pg_ANGLE_ATOL = 1e-12
# Relative tolerance for Float32 interior angles (worst observed about 1 ulp = 6e-8, margin 16x).
const pg_RTOL_F32 = 1e-6

pg_square(id = 1; kwargs...) = Polygon([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]], id; kwargs...)
pg_triangle(id = 1; kwargs...) = Polygon([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], id; kwargs...)

@testset "Polygon fields" begin
    sq = pg_square(3)
    @test sq isa Polygon{Float64}
    @test sq isa AbsSimpleDomain
    @test fieldnames(Polygon) == (:boundary, :corners, :angles, :id)
    @test sq.id == 3
    @test sq.boundary isa Vector{AbsCurve}
    @test sq.corners isa Vector{SVector{2,Float64}}
    @test sq.angles isa Vector{Float64}
    @test sq.corners == [SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(1.0, 1.0), SVector(0.0, 1.0)]
    @test length(sq.boundary) == length(sq.corners) == length(sq.angles) == 4

    # edge i runs from corner i to the next corner (cyclically), is numbered i and carries the domain id
    for i in 1:4
        edge = sq.boundary[i]
        @test edge isa LineSegment{Float64}
        @test edge.pt0 == sq.corners[i]
        @test edge.pt1 == sq.corners[mod1(i + 1, 4)]
        @test edge.segment_id == i
        @test edge.domain_id == 3
        @test edge.orientation == 1
        @test edge.length ≈ 1.0
    end
    # consecutive edges close up into a loop
    @test all(is_connected(sq.boundary))
    @test is_closed(sq.boundary; check_periodic = true)
end

@testset "Polygon angles" begin
    # square: four right angles
    @test all(a -> isapprox(a, pi / 2; atol = pg_ANGLE_ATOL), pg_square().angles)

    # right triangle with legs along the axes: pi/2 at the origin and pi/4 at the other two corners, in corner order
    @test pg_triangle().angles ≈ [pi / 2, pi / 4, pi / 4] atol = pg_ANGLE_ATOL

    # equilateral triangle: pi/3 everywhere
    equilateral = Polygon([[0.0, 0.0], [1.0, 0.0], [0.5, sqrt(3) / 2]], 1)
    @test all(a -> isapprox(a, pi / 3; atol = pg_ANGLE_ATOL), equilateral.angles)

    # regular pentagon of circumradius 1: interior angle 3 pi / 5 and sum (n - 2) pi
    pentagon = Polygon([[cos(2pi * k / 5), sin(2pi * k / 5)] for k in 0:4], 1)
    @test all(a -> isapprox(a, 3pi / 5; atol = pg_ANGLE_ATOL), pentagon.angles)
    @test sum(pentagon.angles) ≈ 3pi atol = pg_ANGLE_ATOL
    @test sum(pg_square().angles) ≈ 2pi atol = pg_ANGLE_ATOL
    @test sum(pg_triangle().angles) ≈ pi atol = pg_ANGLE_ATOL

    # the angle does not depend on the traversal direction or on a translation of the polygon
    clockwise = Polygon([[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0]], 1)
    @test all(a -> isapprox(a, pi / 2; atol = pg_ANGLE_ATOL), clockwise.angles)
    moved = Polygon([[10.0, -3.0], [11.0, -3.0], [10.0, -2.0]], 1)
    @test moved.angles ≈ pg_triangle().angles atol = pg_ANGLE_ATOL
end

@testset "Polygon boundary conditions" begin
    # every edge gets the specular reflection by default
    @test all(e -> e.bc == SpecularReflection(), pg_square().boundary)

    # custom conditions are assigned to the edges in order: edge i runs from corner i to corner i + 1
    bcs = AbsBoundaryCondition[SpecularReflection(), Transparent(2), SymmetryWall(1, 2)]
    tri = pg_triangle(1; bcs = bcs)
    @test [e.bc for e in tri.boundary] == bcs
    @test tri.boundary[2] isa LineSegment{Float64,Transparent}

    # fewer conditions than corners cannot be matched to the edges
    @test_throws BoundsError pg_triangle(1; bcs = [SpecularReflection()])
end

@testset "Polygon element type" begin
    corners32 = [Float32[0, 0], Float32[1, 0], Float32[0, 1]]
    tri32 = Polygon(corners32, 1)
    @test tri32 isa Polygon{Float32}
    @test tri32.corners isa Vector{SVector{2,Float32}}
    @test tri32.angles isa Vector{Float32}
    @test all(e -> e isa LineSegment{Float32}, tri32.boundary)
    @test tri32.angles ≈ Float32[pi / 2, pi / 4, pi / 4] rtol = pg_RTOL_F32

    # corners of mixed precision are promoted to the widest type
    @test Polygon([Float32[0, 0], [1.0, 0.0], Float32[0, 1]], 1) isa Polygon{Float64}

    # tuples and static vectors are accepted as corners
    @test Polygon([(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)], 1).corners == pg_triangle().corners
    @test Polygon([SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(0.0, 1.0)], 1).angles == pg_triangle().angles
end

@testset "Polygon argument errors" begin
    # corners must be points of the plane
    @test_throws DimensionMismatch Polygon([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], 1)
end
