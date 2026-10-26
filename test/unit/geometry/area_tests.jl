# Unit tests for src/geometry/area.jl: area (curve, composite, vector, billiard), fundamental_area, corner_angles, symmetry_reduction_factor, maximal_symmetry.
# Reference values are closed forms (circle, rectangle, right triangle, regular polygon, squares by Green's theorem); helpers use the prefix `ar_`.

using StaticArrays: SVector
using BilliardGeometry: AbsSymmetry, symmetry_reduction_factor, maximal_symmetry

# Relative tolerance for closed-form areas (worst observed relative error 2e-16, margin 5000x).
const ar_RTOL = 1e-12
# Relative tolerance for Float32 areas (worst observed relative error 7.6e-8, margin 13x).
const ar_RTOL_F32 = 1e-6
# Absolute tolerance for corner angles in radians (worst observed 4.4e-16, margin 2000x).
const ar_ANGLE_ATOL = 1e-12

# Closed polygonal loop through the given vertices, one LineSegment per edge, in the vertex order given.
function ar_loop(vertices...)
    n = length(vertices)
    return [LineSegment(vertices[i], vertices[mod1(i + 1, n)]) for i in 1:n]
end

# Axis-aligned square of side s with lower left corner (x0, y0), counterclockwise (ccw = true) or clockwise.
function ar_square(x0, y0, s; ccw = true)
    ps = (SVector(x0, y0), SVector(x0 + s, y0), SVector(x0 + s, y0 + s), SVector(x0, y0 + s))
    return ccw ? ar_loop(ps...) : ar_loop(reverse(ps)...)
end

# Symmetry subtype without a method of its own, to reach the generic fallback of symmetry_reduction_factor.
struct ar_PlainSymmetry <: AbsSymmetry end

@testset "area (curve)" begin
    # a full circle of radius R has area pi R^2, wherever its center is (Green's theorem on a closed curve)
    circ = CircleSegment(2.0, 2pi, 0.0, [0.3, -0.2])
    @test area(circ) ≈ 4pi rtol = ar_RTOL
    @test area(CircleSegment(2.0, 2pi, 0.0, [0.0, 0.0])) ≈ 4pi rtol = ar_RTOL
    # reversing the traversal (negative arc_angle) does not change the area: the result is |.|/2
    @test area(CircleSegment(2.0, -2pi, 0.0, [0.3, -0.2])) ≈ 4pi rtol = ar_RTOL

    # a half circle arc about the origin: the closing chord passes through the origin and adds nothing, so the value is pi R^2 / 2
    @test area(CircleSegment(1.0, pi, 0.0, [0.0, 0.0])) ≈ pi / 2 rtol = ar_RTOL
    # a line through the origin contributes nothing
    @test area(LineSegment([0.0, 0.0], [1.0, 1.0])) == 0.0

    # the quadrature tolerance is forwarded
    @test area(circ; rtol = 1e-6) ≈ 4pi rtol = ar_RTOL
    @test area(circ) isa Float64
    @test area(CircleSegment(1f0, 2f0 * Float32(pi), 0f0, Float32[0, 0])) isa Float32
    @test area(CircleSegment(1f0, 2f0 * Float32(pi), 0f0, Float32[0, 0])) ≈ Float32(pi) rtol = ar_RTOL_F32
end

@testset "area (composite curve)" begin
    # right triangle with legs 2 and 1: area 1, for both traversal directions
    ccw = CompositeCurve(ar_loop(SVector(0.0, 0.0), SVector(2.0, 0.0), SVector(2.0, 1.0)))
    cw = CompositeCurve(ar_loop(SVector(0.0, 0.0), SVector(2.0, 1.0), SVector(2.0, 0.0)))
    @test area(ccw) ≈ 1.0 rtol = ar_RTOL
    @test area(cw) ≈ 1.0 rtol = ar_RTOL
    @test area(ccw) isa Float64

    # translation invariance of the closed loop: the same triangle moved by (10, -3)
    moved = CompositeCurve(ar_loop(SVector(10.0, -3.0), SVector(12.0, -3.0), SVector(12.0, -2.0)))
    @test area(moved) ≈ 1.0 rtol = ar_RTOL

    # half disk: arc of radius 1 followed by its diameter, area pi/2
    half = CompositeCurve(AbsCurve[CircleSegment(1.0, pi, 0.0, [0.0, 0.0]), LineSegment([-1.0, 0.0], [1.0, 0.0])])
    @test area(half) ≈ pi / 2 rtol = ar_RTOL

    # Float32 loop
    f32 = CompositeCurve(ar_loop(SVector(0f0, 0f0), SVector(1f0, 0f0), SVector(1f0, 1f0)))
    @test area(f32) isa Float32
    @test area(f32) ≈ 0.5f0 rtol = ar_RTOL_F32
end

@testset "area (vector of curves)" begin
    @test area(AbsCurve[]) == 0.0

    # unit right triangle, both traversal directions, concrete and abstract element types
    tri = ar_loop(SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(1.0, 1.0))
    @test area(tri) ≈ 0.5 rtol = ar_RTOL
    @test area(AbsCurve[tri...]) ≈ 0.5 rtol = ar_RTOL
    @test area(reverse(tri)) ≈ 0.5 rtol = ar_RTOL
    @test (@inferred area(tri)) isa Float64

    # contributions are summed with their orientation before the absolute value is taken:
    # a ccw outer square of side 2 with a cw inner square of side 1 gives 4 - 1, with a ccw inner square 4 + 1
    outer = ar_square(0.0, 0.0, 2.0)
    @test area(vcat(outer, ar_square(0.5, 0.5, 1.0; ccw = false))) ≈ 3.0 rtol = ar_RTOL
    @test area(vcat(outer, ar_square(0.5, 0.5, 1.0; ccw = true))) ≈ 5.0 rtol = ar_RTOL
    # the sign of the total follows the outer loop: a cw outer loop with a ccw inner loop gives the same 3
    @test area(vcat(ar_square(0.0, 0.0, 2.0; ccw = false), ar_square(0.5, 0.5, 1.0; ccw = true))) ≈ 3.0 rtol = ar_RTOL

    # hole boundaries as stored by the hole billiards: the full circle of radius R = 0.3 has area pi R^2,
    # the quarter arc about the origin (Sinai hole) has the sector area pi R^2 / 4
    @test area(AnnularBilliard(0.3, 0.2).fundamental_domain.boundary[2:2]) ≈ pi * 0.3^2 rtol = ar_RTOL
    @test area(SinaiBilliard(1.0, 0.3).fundamental_domain.boundary[4:4]) ≈ pi * 0.3^2 / 4 rtol = ar_RTOL
end

@testset "area (billiard)" begin
    # circle: pi R^2, independent of the symmetry reduction because the full boundary is used
    @test area(CircleBilliard(1.5)) ≈ pi * 1.5^2 rtol = ar_RTOL
    @test (@inferred area(CircleBilliard(1.5))) isa Float64
    @test area(CircleBilliard(1f0)) isa Float32
    @test area(CircleBilliard(1f0)) ≈ Float32(pi) rtol = ar_RTOL_F32
    # the keyword is forwarded to the quadrature
    @test area(CircleBilliard(1.5); rtol = 1e-6) ≈ pi * 1.5^2 rtol = ar_RTOL

    # RectangleBilliard(a, b) is the rectangle [-a, a] x [-b, b]: area 4 a b
    @test area(RectangleBilliard(1.0, 0.6)) ≈ 4 * 0.6 rtol = ar_RTOL
    @test (@inferred area(RectangleBilliard(1.0, 0.6))) isa Float64

    # right isosceles triangle with unit legs, no symmetry
    @test area(TriangleBilliard(pi / 2, 1.0)) ≈ 0.5 rtol = ar_RTOL

    # polygons: 2 x 1 rectangle and the regular pentagon of circumradius 1, area (5/2) sin(2 pi / 5)
    @test area(PolygonBilliard([SVector(0.0, 0.0), SVector(2.0, 0.0), SVector(2.0, 1.0), SVector(0.0, 1.0)])) ≈ 2.0 rtol = ar_RTOL
    pentagon = PolygonBilliard([SVector(cos(2pi * k / 5), sin(2pi * k / 5)) for k in 0:4])
    @test area(pentagon) ≈ 2.5 * sin(2pi / 5) rtol = ar_RTOL

    # stadium of radius 1 and half-width hw (two half disks and a 2hw x 2 rectangle): pi + 4 hw
    @test area(StadiumBilliard(0.7)) ≈ pi + 4 * 0.7 rtol = ar_RTOL
end

@testset "area (billiards with holes)" begin
    # G1: area(b) adds the area of the holes instead of subtracting it (hole curves are stored with orientation = -1 and the Green integral ignores it)
    @testset "AnnularBilliard" begin
        # AnnularBilliard(R_inner, shift): unit disk with a hole of radius R_inner centered at (shift, 0): pi (1 - R_inner^2)
        @test_broken area(AnnularBilliard(0.3, 0.2)) ≈ pi * (1 - 0.3^2) rtol = ar_RTOL
    end
    @testset "SinaiBilliard" begin
        # SinaiBilliard(a, R): square [-a, a]^2 with a disk of radius R removed: 4 a^2 - pi R^2
        @test_broken area(SinaiBilliard(1.0, 0.3)) ≈ 4 - pi * 0.3^2 rtol = ar_RTOL
    end
    @testset "SquareWithinSquareBilliard" begin
        # square [-1, 1]^2 with two holes of side 1/4: 4 - 2 (1/4)^2
        @test_broken area(SquareWithinSquareBilliard(1.0)) ≈ 4 - 2 * 0.25^2 rtol = ar_RTOL
    end
end

@testset "symmetry_reduction_factor (symmetry)" begin
    @test symmetry_reduction_factor(XAxisReflection()) == 2
    @test symmetry_reduction_factor(YAxisReflection()) == 2
    @test symmetry_reduction_factor(XYAxisReflection()) == 4
    @test symmetry_reduction_factor(DiagonalReflection()) == 2
    @test symmetry_reduction_factor(AntiDiagonalReflection()) == 2
    # a rotation of order N divides the area by N, whatever its power m
    @test symmetry_reduction_factor(NFoldRotation(3, 1)) == 3
    @test symmetry_reduction_factor(NFoldRotation(5, 2)) == 5
    # composite reflections are not reduced
    @test symmetry_reduction_factor(CompositeReflection(XAxisReflection(), YAxisReflection())) == 1
    # a symmetry without a method of its own falls back to 1
    @test symmetry_reduction_factor(ar_PlainSymmetry()) == 1
    @test (@inferred symmetry_reduction_factor(XYAxisReflection())) isa Int
end

@testset "symmetry_reduction_factor (billiard)" begin
    # the factor is that of the largest generator in the registry
    @test symmetry_reduction_factor(CircleBilliard(1.0)) == 4       # D2
    @test symmetry_reduction_factor(StadiumBilliard(0.7)) == 4      # D2
    @test symmetry_reduction_factor(SquareWithinSquareBilliard(1.0)) == 2   # diagonal reflection
    @test symmetry_reduction_factor(C3Billiard(0.1)) == 3           # C3
    # no registered symmetry: factor 1
    @test symmetry_reduction_factor(TriangleBilliard(pi / 2, 1.0)) == 1
    @test symmetry_reduction_factor(AnnularBilliard(0.3, 0.2)) == 1
    @test (@inferred symmetry_reduction_factor(CircleBilliard(1.0))) isa Int
end

@testset "maximal_symmetry" begin
    @test maximal_symmetry(CircleBilliard(1.0)) isa XYAxisReflection
    @test maximal_symmetry(RectangleBilliard(1.0, 0.6)) isa XYAxisReflection
    @test maximal_symmetry(SquareWithinSquareBilliard(1.0)) isa DiagonalReflection
    c3 = maximal_symmetry(C3Billiard(0.1))
    @test c3 isa NFoldRotation
    @test c3.order == 3
    # the returned generator is the registered one, with its sym_id
    @test maximal_symmetry(CircleBilliard(1.0)).sym_id == 2
    # empty registry: nothing
    @test maximal_symmetry(TriangleBilliard(pi / 2, 1.0)) === nothing
    @test maximal_symmetry(AnnularBilliard(0.3, 0.2)) === nothing
end

@testset "fundamental_area" begin
    # the physical area divided by the reduction factor
    @test fundamental_area(CircleBilliard(1.5)) ≈ pi * 1.5^2 / 4 rtol = ar_RTOL
    @test fundamental_area(RectangleBilliard(1.0, 0.6)) ≈ 0.6 rtol = ar_RTOL
    @test fundamental_area(StadiumBilliard(0.7)) ≈ (pi + 4 * 0.7) / 4 rtol = ar_RTOL
    # without symmetry the fundamental domain is the whole billiard
    @test fundamental_area(TriangleBilliard(pi / 2, 1.0)) ≈ 0.5 rtol = ar_RTOL
    # fundamental area times reduction factor gives the physical area
    b = CircleBilliard(1.5)
    @test fundamental_area(b) * symmetry_reduction_factor(b) ≈ area(b) rtol = ar_RTOL
    c = C3Billiard(0.1)
    @test fundamental_area(c) * symmetry_reduction_factor(c) ≈ area(c) rtol = ar_RTOL
    @test (@inferred fundamental_area(CircleBilliard(1.5))) isa Float64
end

@testset "fundamental_area (billiards with holes)" begin
    # G1: follows area(b), which adds the hole area
    @testset "SinaiBilliard" begin
        @test_broken fundamental_area(SinaiBilliard(1.0, 0.3)) ≈ (4 - pi * 0.3^2) / 4 rtol = ar_RTOL
    end
    @testset "SquareWithinSquareBilliard" begin
        @test_broken fundamental_area(SquareWithinSquareBilliard(1.0)) ≈ (4 - 2 * 0.25^2) / 2 rtol = ar_RTOL
    end
end

@testset "corner_angles (polygons)" begin
    # the interior angles of an n-gon sum to (n - 2) pi
    tri = TriangleBilliard(pi / 2, 1.0)
    @test sort(corner_angles(tri)) ≈ [pi / 4, pi / 4, pi / 2] atol = ar_ANGLE_ATOL
    @test sum(corner_angles(tri)) ≈ pi atol = ar_ANGLE_ATOL

    rect = PolygonBilliard([SVector(0.0, 0.0), SVector(2.0, 0.0), SVector(2.0, 1.0), SVector(0.0, 1.0)])
    @test corner_angles(rect) ≈ fill(pi / 2, 4) atol = ar_ANGLE_ATOL
    @test sum(corner_angles(rect)) ≈ 2pi atol = ar_ANGLE_ATOL

    # regular pentagon: five angles 3 pi / 5
    pentagon = PolygonBilliard([SVector(cos(2pi * k / 5), sin(2pi * k / 5)) for k in 0:4])
    @test corner_angles(pentagon) ≈ fill(3pi / 5, 5) atol = ar_ANGLE_ATOL
    @test sum(corner_angles(pentagon)) ≈ 3pi atol = ar_ANGLE_ATOL

    # a general triangle TriangleBilliard(gamma, chi): angles alpha + beta + gamma = pi
    gen = TriangleBilliard(1.2, 2.0)
    @test sum(corner_angles(gen)) ≈ pi atol = ar_ANGLE_ATOL
    @test any(a -> isapprox(a, 1.2; atol = ar_ANGLE_ATOL), corner_angles(gen))
end

@testset "corner_angles (fundamental domain)" begin
    # fundamental domain of RectangleBilliard: the quadrant rectangle, four right angles
    @test corner_angles(RectangleBilliard(1.0, 0.6)) ≈ fill(pi / 2, 4) atol = ar_ANGLE_ATOL
    # fundamental domain of the circle: the quarter disk, three right angles (arc/wall twice, wall/wall once)
    @test corner_angles(CircleBilliard(1.0)) ≈ fill(pi / 2, 3) atol = ar_ANGLE_ATOL
    # stadium fundamental domain: quarter stadium, three right angles
    @test corner_angles(StadiumBilliard(0.7)) ≈ fill(pi / 2, 3) atol = ar_ANGLE_ATOL
    @test eltype(corner_angles(CircleBilliard(1.0))) == Float64
    @test eltype(corner_angles(CircleBilliard(1f0))) == Float32
end

@testset "corner_angles (full boundary)" begin
    # fundamental = false uses the full physical boundary: a smooth circle has no corners, a rectangle keeps its four
    @test isempty(corner_angles(CircleBilliard(1.0); fundamental = false))
    @test eltype(corner_angles(CircleBilliard(1.0); fundamental = false)) == Float64
    @test corner_angles(RectangleBilliard(1.0, 0.6); fundamental = false) ≈ fill(pi / 2, 4) atol = ar_ANGLE_ATOL
end

@testset "corner_angles (angle_tol)" begin
    # turning angles below angle_tol are smooth joins, not corners: a right-angle corner turns by pi/2 < 2
    @test isempty(corner_angles(RectangleBilliard(1.0, 0.6); angle_tol = 2.0))
    @test length(corner_angles(RectangleBilliard(1.0, 0.6); angle_tol = 1.0)) == 4
end

@testset "corner_angles (billiards with holes)" begin
    # G2: the inner arc of SinaiBilliard is not chained into the fundamental boundary, so two corners are reported as 0 instead of pi/2
    sinai = SinaiBilliard(1.0, 0.3)
    @test length(corner_angles(sinai)) == 5
    @test_broken corner_angles(sinai) ≈ fill(pi / 2, 5) atol = ar_ANGLE_ATOL
end
