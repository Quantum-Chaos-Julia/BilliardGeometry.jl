# Unit tests for src/geometry/segments/compositecurves.jl: the CompositeCurve type and its constructor.
# Reference values are lengths and end points of explicit polygonal and circular pieces; helpers use the prefix `cc_`.

using StaticArrays: SVector

# Relative tolerance for lengths summed from closed-form piece lengths (worst observed 1 ulp = 2e-16, margin 50x).
const cc_RTOL = 1e-14
# Relative tolerance for Float32 lengths (worst observed 1 ulp = 6e-8, margin 16x).
const cc_RTOL_F32 = 1e-6
# Absolute tolerance for end points of circular arcs (worst observed 1.2e-16 from cos(pi), margin 8000x).
const cc_POINT_ATOL = 1e-12

# Right triangle with vertices (0,0), (1,0), (1,1): legs of length 1 and hypotenuse sqrt(2), pieces in boundary order.
cc_triangle() = [LineSegment([0.0, 0.0], [1.0, 0.0]), LineSegment([1.0, 0.0], [1.0, 1.0]), LineSegment([1.0, 1.0], [0.0, 0.0])]

@testset "CompositeCurve fields" begin
    legs = cc_triangle()
    cc = CompositeCurve(legs)
    @test cc isa CompositeCurve{Float64}
    @test cc isa AbsCompositeCurve
    @test fieldnames(CompositeCurve) == (:subcurves, :end_lengths, :corners, :length)
    @test cc.subcurves isa Vector{AbsCurve}
    @test cc.end_lengths isa Vector{Float64}
    @test cc.corners isa Vector{SVector{2,Float64}}
    @test cc.length isa Float64

    # the pieces are kept in the order given, since consecutive pieces already share their end points
    @test length(cc.subcurves) == 3
    @test all(i -> cc.subcurves[i] === legs[i], 1:3)

    # end_lengths holds the cumulative arc length at each junction, starting from zero: 0, 1, 2, 2 + sqrt(2)
    @test length(cc.end_lengths) == 4
    @test cc.end_lengths[1] == 0.0
    @test cc.end_lengths[2] ≈ 1.0 rtol = cc_RTOL
    @test cc.end_lengths[3] ≈ 2.0 rtol = cc_RTOL
    @test cc.end_lengths[4] ≈ 2.0 + sqrt(2.0) rtol = cc_RTOL
    @test cc.length ≈ 2.0 + sqrt(2.0) rtol = cc_RTOL
    @test cc.end_lengths[end] == cc.length

    # the increments are the lengths of the pieces
    @test diff(cc.end_lengths) ≈ [c.length for c in cc.subcurves] rtol = cc_RTOL

    # corners hold the start of the first piece and the end of every piece, so there is one more corner than pieces
    @test cc.corners == [SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(1.0, 1.0), SVector(0.0, 0.0)]
    @test length(cc.corners) == length(cc.subcurves) + 1
end

@testset "CompositeCurve (ordering of pieces)" begin
    # a piece order that is already a connected chain is kept, wherever the chain starts
    legs = cc_triangle()
    rotated = CompositeCurve([legs[3], legs[1], legs[2]])
    @test [c.pt0 for c in rotated.subcurves] == [SVector(1.0, 1.0), SVector(0.0, 0.0), SVector(1.0, 0.0)]
    @test rotated.corners == [SVector(1.0, 1.0), SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(1.0, 1.0)]
    @test rotated.length ≈ 2.0 + sqrt(2.0) rtol = cc_RTOL

    # pieces given out of order are chained so that every piece starts where the previous one ends
    shuffled = [legs[1], legs[3], legs[2]]
    cc = CompositeCurve(shuffled)
    @test [c.pt0 for c in cc.subcurves] == [SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(1.0, 1.0)]
    @test cc.corners == [SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(1.0, 1.0), SVector(0.0, 0.0)]
    @test cc.end_lengths[2:end] ≈ [1.0, 2.0, 2.0 + sqrt(2.0)] rtol = cc_RTOL
    @test cc.length ≈ 2.0 + sqrt(2.0) rtol = cc_RTOL
    # the caller's vector is not reordered
    @test shuffled[2] === legs[3]
    @test shuffled[3] === legs[2]
end

@testset "CompositeCurve (single piece)" begin
    line = LineSegment([0.0, 0.0], [3.0, 4.0])
    cc = CompositeCurve([line])
    @test cc.length ≈ 5.0 rtol = cc_RTOL
    @test cc.end_lengths ≈ [0.0, 5.0] rtol = cc_RTOL
    @test cc.corners == [SVector(0.0, 0.0), SVector(3.0, 4.0)]
    @test length(cc.subcurves) == 1
    @test cc.subcurves[1] === line
end

@testset "CompositeCurve (mixed piece types)" begin
    # half disk: the upper unit semicircle followed by the diameter; length pi + 2
    arc = CircleSegment(1.0, pi, 0.0, [0.0, 0.0])
    diameter = LineSegment([-1.0, 0.0], [1.0, 0.0])
    cc = CompositeCurve(AbsCurve[arc, diameter])
    @test cc isa CompositeCurve{Float64}
    @test cc.subcurves[1] === arc
    @test cc.subcurves[2] === diameter
    @test cc.end_lengths ≈ [0.0, pi, pi + 2.0] rtol = cc_RTOL
    @test cc.length ≈ pi + 2.0 rtol = cc_RTOL

    # the arc starts at (1, 0) and ends at (-1, 0), where the diameter starts; the loop closes at (1, 0)
    @test length(cc.corners) == 3
    @test cc.corners[1] ≈ SVector(1.0, 0.0) atol = cc_POINT_ATOL
    @test cc.corners[2] ≈ SVector(-1.0, 0.0) atol = cc_POINT_ATOL
    @test cc.corners[3] ≈ SVector(1.0, 0.0) atol = cc_POINT_ATOL

    # the same loop traversed from the diameter: it starts at (-1, 0), and the loop closes there
    from_diameter = CompositeCurve(AbsCurve[diameter, arc])
    @test from_diameter.length ≈ pi + 2.0 rtol = cc_RTOL
    @test from_diameter.corners[1] ≈ SVector(-1.0, 0.0) atol = cc_POINT_ATOL
    @test from_diameter.corners[3] ≈ SVector(-1.0, 0.0) atol = cc_POINT_ATOL
end

@testset "CompositeCurve (element type)" begin
    # the numeric type is taken from the first piece
    legs32 = [LineSegment(Float32[0, 0], Float32[1, 0]), LineSegment(Float32[1, 0], Float32[1, 1])]
    cc = CompositeCurve(legs32)
    @test cc isa CompositeCurve{Float32}
    @test cc.end_lengths isa Vector{Float32}
    @test cc.corners isa Vector{SVector{2,Float32}}
    @test cc.length isa Float32
    @test cc.length ≈ 2.0f0 rtol = cc_RTOL_F32
    @test cc.corners == SVector{2,Float32}[SVector(0.0f0, 0.0f0), SVector(1.0f0, 0.0f0), SVector(1.0f0, 1.0f0)]

    # abstractly typed piece vectors give the same type
    @test CompositeCurve(AbsCurve[legs32...]) isa CompositeCurve{Float32}
end

@testset "CompositeCurve (argument errors)" begin
    # the type and the start corner are read from the first piece, so an empty piece list fails with a BoundsError
    @test_throws BoundsError CompositeCurve(AbsCurve[])
end
