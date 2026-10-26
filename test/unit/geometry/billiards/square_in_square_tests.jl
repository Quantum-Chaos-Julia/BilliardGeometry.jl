# Unit tests for the square-within-square billiard constructor (structural checks only):
#   src/geometry/billiards/square_in_square.jl   (SquareWithinSquareBilliard{T}, SquareWithinSquareBilliard(a; positions, sizes))
# Reference object: the half of [-a, a]^2 below the diagonal y = x with square holes centred on the diagonal at `positions` with side
# lengths `sizes`; one DiagonalReflection symmetry. The fundamental domain is a SimpleDomain of 3 + 3 n LineSegments for n holes:
# two outer walls, the diagonal walls (SymmetryWall(1, 2)) and two hole edges (orientation -1, domain ids 2..n+1) per hole.
# Only structure is checked. Guards: a > 0, equal lengths of positions and sizes, positive sizes, holes strictly inside, no overlap.

using StaticArrays: SVector

@testset "SquareWithinSquareBilliard fields" begin
    sws = SquareWithinSquareBilliard(1.0)
    @test sws isa SquareWithinSquareBilliard{Float64}
    @test sws isa AbsBilliard
    @test sws.fundamental_domain isa SimpleDomain{Float64}
    @test sws.symmetries isa SymmetryRegistry
    @test fieldnames(SquareWithinSquareBilliard) == (:fundamental_domain, :symmetries)
    @test @inferred(SquareWithinSquareBilliard(1.0)) isa SquareWithinSquareBilliard{Float64}
end

@testset "SquareWithinSquareBilliard symmetries" begin
    reg = SquareWithinSquareBilliard(1.0).symmetries
    @test length(reg) == 1
    @test reg[1] isa DiagonalReflection
    @test reg[1].sym_id == 1
end

@testset "SquareWithinSquareBilliard curves" begin
    dom = SquareWithinSquareBilliard(1.0).fundamental_domain
    @test dom.id == 1
    # two holes by default: 3 + 3 * 2 curves, one corner per curve
    @test length(dom.boundary) == 9
    @test length(dom.corners) == 9
    @test all(c -> c isa LineSegment{Float64}, dom.boundary)
    @test [c.bc for c in dom.boundary] == [
        SpecularReflection(), SpecularReflection(), SymmetryWall(1, 2),
        SpecularReflection(), SpecularReflection(), SymmetryWall(1, 2),
        SpecularReflection(), SpecularReflection(), SymmetryWall(1, 2),
    ]
    # hole edges carry orientation -1 and the domain ids of their hole; the hole with the larger position is built first
    @test [c.orientation for c in dom.boundary] == [1, 1, 1, -1, -1, 1, -1, -1, 1]
    @test [c.domain_id for c in dom.boundary] == [1, 1, 1, 3, 3, 1, 2, 2, 1]
    @test [c.segment_id for c in dom.boundary] == [1, 2, 3, 1, 2, 4, 1, 2, 5]

    # the outer walls are the bottom and right sides of [-a, a]^2
    @test dom.corners[1:3] == [SVector(-1.0, -1.0), SVector(1.0, -1.0), SVector(1.0, 1.0)]
    @test dom.boundary[1].pt1 == SVector(1.0, -1.0)
    @test dom.boundary[2].pt1 == SVector(1.0, 1.0)
    # the boundary closes on the lower-left corner
    @test dom.boundary[end].pt1 == SVector(-1.0, -1.0)
end

@testset "SquareWithinSquareBilliard hole count" begin
    one_hole = SquareWithinSquareBilliard(1.0; positions=[0.0], sizes=[0.4]).fundamental_domain
    @test length(one_hole.boundary) == 6
    @test [c.orientation for c in one_hole.boundary] == [1, 1, 1, -1, -1, 1]
    @test [c.domain_id for c in one_hole.boundary] == [1, 1, 1, 2, 2, 1]

    three_holes = SquareWithinSquareBilliard(1.0; positions=[0.5, -0.5, 0.0], sizes=[0.2, 0.2, 0.2]).fundamental_domain
    @test length(three_holes.boundary) == 12
    @test sort(unique(c.domain_id for c in three_holes.boundary if c.orientation == -1)) == [2, 3, 4]
    @test count(c -> c.orientation == -1, three_holes.boundary) == 6

    # unsorted positions are sorted with their sizes: the hole at 0.5 has side 0.1 and its diagonal wall ends at 0.5 + 0.05
    unsorted = SquareWithinSquareBilliard(1.0; positions=[0.5, -0.5], sizes=[0.1, 0.2]).fundamental_domain
    @test unsorted.boundary[3].pt1 ≈ SVector(0.55, 0.55) atol = 1e-14
end

@testset "SquareWithinSquareBilliard topology" begin
    sws = SquareWithinSquareBilliard(1.0)
    # the holes are part of the flat boundary of a SimpleDomain, so the domain itself has genus zero and one component
    @test genus(sws) == 0
    @test length(boundary_components(sws)) == 1
    @test length(geometric_boundary_components(sws)) == 1
end

@testset "SquareWithinSquareBilliard argument guards" begin
    # a > 0
    @test_throws ArgumentError SquareWithinSquareBilliard(0.0)
    @test_throws ArgumentError SquareWithinSquareBilliard(-1.0)
    # positions and sizes have the same length
    @test_throws ArgumentError SquareWithinSquareBilliard(1.0; positions=[0.0], sizes=[0.2, 0.3])
    # sizes are positive
    @test_throws ArgumentError SquareWithinSquareBilliard(1.0; positions=[0.0], sizes=[0.0])
    @test_throws ArgumentError SquareWithinSquareBilliard(1.0; positions=[0.0], sizes=[-0.2])
    # holes lie strictly inside the outer square, on either side
    @test_throws ArgumentError SquareWithinSquareBilliard(1.0; positions=[0.9], sizes=[0.3])
    @test_throws ArgumentError SquareWithinSquareBilliard(1.0; positions=[-0.9], sizes=[0.3])
    @test_throws ArgumentError SquareWithinSquareBilliard(1.0; positions=[0.9], sizes=[0.2])
    # holes may neither overlap nor touch
    @test_throws ArgumentError SquareWithinSquareBilliard(1.0; positions=[-0.5, -0.4], sizes=[0.4, 0.4])
    @test_throws ArgumentError SquareWithinSquareBilliard(1.0; positions=[-0.5, -0.3], sizes=[0.2, 0.2])
end

@testset "SquareWithinSquareBilliard eltype promotion" begin
    sws32 = SquareWithinSquareBilliard(1.0f0)
    @test sws32 isa SquareWithinSquareBilliard{Float32}
    @test sws32.fundamental_domain isa SimpleDomain{Float32}
    @test all(c -> c isa LineSegment{Float32}, sws32.fundamental_domain.boundary)
    @test eltype(sws32.fundamental_domain.corners) == SVector{2,Float32}
    # Float64 positions and sizes are converted to the element type of a
    sws32 = SquareWithinSquareBilliard(1.0f0; positions=[-0.5, 0.5], sizes=[0.25, 0.25])
    @test sws32 isa SquareWithinSquareBilliard{Float32}
    @test all(c -> c isa LineSegment{Float32}, sws32.fundamental_domain.boundary)
end
