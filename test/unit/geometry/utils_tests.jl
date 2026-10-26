# Unit tests for curve-chain and small geometric helpers:
#   src/geometry/utils.jl   (is_overlaping, is_connected, is_closed, connect_curves, signed_angle, find_unique_elements,
#                            point_curve_parameter)
# Reference objects: the unit square traversed counter-clockwise as four LineSegments (0,0)->(1,0)->(1,1)->(0,1)->(0,0), whose
# consecutive end/start points coincide exactly; a gap piece (5,5)->(0,5) that breaks the chain; two separate L-shaped chains for
# connect_curves; the unit circle quadrants of full_boundary(CircleBilliard(1.0)) whose joints carry roundoff of order 1e-16 at
# coordinates equal to zero; the oriented angle between planar vectors, atan(a x b, a . b); the chord (0,0)->(2,1) and the
# circle arc x = cos(pi t), y = sin(pi t) for the inverse parametrization. Coordinates equal to zero need an absolute tolerance
# sqrt(eps(T)) in is_overlaping, a relative-only isapprox rejects (-1.8e-16, 6e-17). No eigen solves, no pinned wavefunction
# values. Helper functions use the prefix `ut_`.

using StaticArrays: SVector
using BilliardGeometry: signed_angle   # not exported
# Unit square sides in counter-clockwise order, segment_id marks the position.
function ut_square_sides()
    return [LineSegment([0.0, 0.0], [1.0, 0.0]; segment_id=1),
            LineSegment([1.0, 0.0], [1.0, 1.0]; segment_id=2),
            LineSegment([1.0, 1.0], [0.0, 1.0]; segment_id=3),
            LineSegment([0.0, 1.0], [0.0, 0.0]; segment_id=4)]
end

# Position of the curve `c` in `list` by identity, nothing if absent.
ut_position(list, c) = findfirst(x -> x === c, list)

@testset "is_overlaping" begin
    # identical and clearly separated points
    @test is_overlaping((0.3, 0.4), (0.3, 0.4))
    @test is_overlaping(SVector(0.3, 0.4), SVector(0.3, 0.4))
    @test !is_overlaping((0.0, 0.0), (1.0, 0.0))
    @test !is_overlaping((0.0, 0.0), (0.0, 1.0))

    # both coordinates must agree: a mismatch in one of them is enough to reject
    @test !is_overlaping((0.5, 1.0), (0.5, 1.0 + 1e-3))
    @test !is_overlaping((0.5, 1.0), (0.5 + 1e-3, 1.0))
    @test !is_overlaping((0.0, 0.0), (1e-3, 1e-3))

    # G3: isapprox without atol fails for coordinates that are roundoff-close to zero
    @test is_overlaping((-1.8e-16, -1.0), (6e-17, -1.0))
    @test is_overlaping((6e-17, -1.8e-16), (-1.0e-16, 3e-17))
    @test is_overlaping((0.0, 1.0), (1e-17, 1.0))

    # G3: points 1e-3 apart near zero stay distinct
    @test !is_overlaping((0.0, -1.0), (1e-3, -1.0))
    @test !is_overlaping((-1.0e-3, 0.0), (0.0, 0.0))

    # G3: Float32 points, tolerance sqrt(eps(Float32)) ~ 3.4e-4
    @test is_overlaping((-1.8f-8, -1.0f0), (6f-9, -1.0f0))
    @test is_overlaping(SVector(0.0f0, 1.0f0), SVector(1f-8, 1.0f0))
    @test !is_overlaping((0.0f0, -1.0f0), (1f-3, -1.0f0))
    @test !is_overlaping(SVector(0.0f0, 0.0f0), SVector(0.0f0, 1f-3))
end

@testset "is_connected" begin
    sq = ut_square_sides()

    # closed square: every end point meets the next start point, including the periodic link
    for start_initial in (true, false)
        conn = is_connected(sq; start_initial=start_initial)
        @test conn isa Vector{Bool}
        @test length(conn) == 4
        @test all(conn)
    end

    # open chain of three sides: only the periodic link (last end -> first start) is broken
    open3 = sq[1:3]
    conn = is_connected(open3; start_initial=true)
    @test conn isa Vector{Bool}
    @test conn == [false, true, true]
    # start_initial = false pairs end[i] with start[i+1], so the periodic link is the last entry
    conn = is_connected(open3; start_initial=false)
    @test conn isa Vector{Bool}
    @test conn == [true, true, false]

    # a gap piece in the middle breaks the chain at two links
    gapped = [sq[1], sq[2], LineSegment([5.0, 5.0], [0.0, 5.0])]
    @test is_connected(gapped; start_initial=true) == [false, true, false]
    @test is_connected(gapped; start_initial=false) == [true, false, false]
end

@testset "is_closed" begin
    sq = ut_square_sides()

    # a closed square for both pairing conventions, with and without the periodic link
    @test is_closed(sq)
    @test is_closed(sq; start_initial=false)
    @test is_closed(sq; check_periodic=true)
    @test is_closed(sq; check_periodic=true, start_initial=false)

    # the open chain of three sides is closed in the default non-periodic sense but not periodically
    open3 = sq[1:3]
    @test is_closed(open3)
    @test is_closed(open3; start_initial=false)
    @test !is_closed(open3; check_periodic=true)
    @test !is_closed(open3; check_periodic=true, start_initial=false)

    # an interior gap breaks the chain in every mode
    gapped = [sq[1], sq[2], LineSegment([5.0, 5.0], [0.0, 5.0])]
    @test !is_closed(gapped)
    @test !is_closed(gapped; start_initial=false)
    @test !is_closed(gapped; check_periodic=true)

    # G3: the quadrant arcs of the circle meet only up to roundoff at zero coordinates
    circle_boundary = full_boundary(CircleBilliard(1.0))
    @test is_closed(circle_boundary)
    @test all(is_connected(circle_boundary))
    @test is_closed(circle_boundary; check_periodic=true)
end

@testset "connect_curves" begin
    sq = ut_square_sides()

    # an already closed list is returned as it is
    @test connect_curves(sq) === sq

    # shuffled sides are chained into traversal order, starting with the first curve of the input
    shuffled = AbsCurve[sq[3], sq[1], sq[4], sq[2]]
    chained = connect_curves(shuffled)
    @test chained isa Vector{AbsCurve}
    @test length(chained) == 4
    @test chained[1] === shuffled[1]
    @test all(c -> ut_position(chained, c) !== nothing, sq)
    @test [c.segment_id for c in chained] == [3, 4, 1, 2]
    @test is_closed(chained; check_periodic=true)
    # the input list is not reordered
    @test [c.segment_id for c in shuffled] == [3, 1, 4, 2]

    # a concretely typed input vector is accepted too
    typed = [sq[2], sq[4], sq[1], sq[3]]
    @test typed isa Vector{<:LineSegment}
    chained = connect_curves(typed)
    @test [c.segment_id for c in chained] == [2, 3, 4, 1]

    # disconnected pieces stay in separate chains, each chain keeps its internal order
    a1 = LineSegment([0.0, 0.0], [1.0, 0.0])
    a2 = LineSegment([1.0, 0.0], [1.0, 1.0])
    b1 = LineSegment([5.0, 5.0], [6.0, 5.0])
    b2 = LineSegment([6.0, 5.0], [6.0, 6.0])
    pieces = AbsCurve[a2, b1, a1, b2]
    chained = connect_curves(pieces)
    @test length(chained) == 4
    @test all(c -> ut_position(chained, c) !== nothing, pieces)
    @test ut_position(chained, a2) == ut_position(chained, a1) + 1
    @test ut_position(chained, b2) == ut_position(chained, b1) + 1
    @test !is_closed(chained)

    # G3: joints that agree only up to roundoff at zero coordinates are still chained
    s1 = LineSegment([6e-17, 0.0], [1.0, 0.0]; segment_id=1)
    s2 = LineSegment([1.0, 0.0], [1.0, 1.0]; segment_id=2)
    s3 = LineSegment([1.0, 1.0], [0.0, 1.0]; segment_id=3)
    s4 = LineSegment([0.0, 1.0], [-1.8e-16, 0.0]; segment_id=4)
    chained = connect_curves(AbsCurve[s3, s1, s4, s2])
    @test length(chained) == 4
    @test [c.segment_id for c in chained] == [3, 4, 1, 2]
    @test is_closed(chained; check_periodic=true)
end

@testset "signed_angle" begin
    ex, ey = [1.0, 0.0], [0.0, 1.0]
    # quarter turns, counter-clockwise positive
    @test signed_angle(ex, ey) ≈ pi / 2 atol = 1e-14
    @test signed_angle(ey, ex) ≈ -pi / 2 atol = 1e-14
    # parallel vectors
    @test signed_angle(ex, 3 .* ex) == 0.0
    # antiparallel vectors give +-pi
    @test abs(signed_angle(ex, -ex)) ≈ pi atol = 1e-14
    # eighth turn and scale invariance
    @test signed_angle(ex, [1.0, 1.0]) ≈ pi / 4 atol = 1e-14
    @test signed_angle(2.5 .* ex, 0.1 .* [1.0, 1.0]) ≈ pi / 4 atol = 1e-14
    # antisymmetric in its arguments
    a, b = [0.3, -1.2], [2.0, 0.7]
    @test signed_angle(a, b) ≈ -signed_angle(b, a) atol = 1e-14
    # rotation by theta gives theta
    theta = 0.7
    @test signed_angle(a, [cos(theta) -sin(theta); sin(theta) cos(theta)] * a) ≈ theta atol = 1e-14
    # SVector and Float32 input
    @test signed_angle(SVector(1.0f0, 0.0f0), SVector(0.0f0, 1.0f0)) isa Float32
    @test signed_angle(SVector(1.0f0, 0.0f0), SVector(0.0f0, 1.0f0)) ≈ Float32(pi / 2)
end

@testset "find_unique_elements" begin
    # G9: the input must not be mutated
    v = [0, 0, 1, 1, 2]
    v0 = copy(v)
    result = find_unique_elements(v)
    @test v == v0
    @test result == [0, 1, 2]

    # sorted-unique result for runs of equal floats, input untouched
    w = [0.5, 0.5, 0.5, 1.0, 2.0, 2.0]
    w0 = copy(w)
    @test find_unique_elements(w) == [0.5, 1.0, 2.0]
    @test w == w0

    # G9: near-equal values merge with atol = sqrt(eps(T)); the last element of each run is kept
    x = [0.0, 1e-17, 0.5, 0.5, 1.0]
    x0 = copy(x)
    result = find_unique_elements(x)
    @test length(result) == 3
    @test result == [1e-17, 0.5, 1.0]
    @test x == x0

    # a run of three near-equal values keeps only the last one
    y = [1.0, 1.0 + 1e-12, 1.0 + 2e-12, 2.0]
    @test find_unique_elements(y) == [1.0 + 2e-12, 2.0]

    # G9: values 1e-3 apart are distinct
    z = [0.1, 0.101, 0.102]
    @test find_unique_elements(z) == z

    # G9: Float32 tolerance sqrt(eps(Float32)) ~ 3.4e-4
    f = Float32[0.0, 1f-8, 1.0]
    result = find_unique_elements(f)
    @test result isa Vector{Float32}
    @test result == Float32[1f-8, 1.0]

    # degenerate lengths
    @test find_unique_elements([1.0]) == [1.0]
    @test isempty(find_unique_elements(Float64[]))
end

@testset "point_curve_parameter" begin
    # LineSegment: y(t) = t has the single root t0, so the parameter is returned directly
    line = LineSegment([0.0, 0.0], [2.0, 1.0])
    for t0 in (0.1, 0.3, 0.9)
        @test point_curve_parameter(line, curve(line, t0)) ≈ t0 atol = 1e-8
    end
    # plain Vector points are accepted
    @test point_curve_parameter(line, [1.0, 0.5]) ≈ 0.5 atol = 1e-8

    # CircleSegment, quarter arc: y(t) = sin(pi t / 2) is monotone on [0, 1]
    quarter = CircleSegment(1.0, pi / 2, 0.0, 0.0, 0.0)
    for t0 in (0.2, 0.4, 0.85)
        @test point_curve_parameter(quarter, curve(quarter, t0)) ≈ t0 atol = 1e-8
    end

    # CircleSegment, half arc: y(t) = sin(pi t) has two roots, the x coordinate picks the matching one
    half = CircleSegment(1.0, pi, 0.0, 0.0, 0.0)
    for t0 in (0.25, 0.75)
        @test point_curve_parameter(half, curve(half, t0)) ≈ t0 atol = 1e-8
    end
end
