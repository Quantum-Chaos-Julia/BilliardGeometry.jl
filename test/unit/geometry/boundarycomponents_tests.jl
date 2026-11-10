# Unit tests for composite boundary components (a vector of connected curves parametrized by one global periodic parameter):
#   src/geometry/boundarycomponents.jl   (component_lengths, _junction_angle, _is_true_corner, _component_corner_locations,
#                                         _component_corner_angles, print_component_junctions, _global_t_to_segment_u,
#                                         _eval_composite_geom_global_t, and their methods taking the precomputed result
#                                         (lens, cum, Ltot) of component_lengths, which must equal the comp methods exactly)
#   src/geometry/boundarytypes.jl        (_group_curves_by_domain_id, used here for the grouping of flat curve lists)
# Reference objects: the right-angle pair of unit lines, the CCW unit square, the CCW right triangle (0,0),(1,0),(0,1), the
# half disk (arc of the unit circle followed by its diameter) and the full unit circle. For these the junction angles, the
# cumulative lengths, the global parameter t = 2 pi s / L of an arc length s and the interior corner angles are known in closed
# form. Corner angles are only tested for convex components traversed counterclockwise, where the unsigned junction angle
# gives the interior angle pi - tau; the signed turning angle and the hole convention (angle pi + tau on a hole loop traversed
# with the hole on the left) are tested separately. No eigen solves, no pinned
# wavefunction values. Helper functions and constants use the prefix `bc_`.

using StaticArrays: SVector
using LinearAlgebra: norm
using Test: collect_test_logs
using BilliardGeometry: _is_true_corner, _junction_angle, _component_corner_locations, _component_corner_angles,
                        _global_t_to_segment_u, _eval_composite_geom_global_t, _group_curves_by_domain_id,
                        _signed_junction_angle

# Absolute tolerance for closed-form comparisons of angles and trigonometric expressions (a few ulps of Float64).
const bc_ATOL = 1e-12
# Absolute tolerance for Float32 results (a few ulps of Float32).
const bc_ATOL_F32 = 1.0f-6

# Closed counterclockwise unit square, right triangle and half disk as composite components.
bc_square() = [LineSegment([0.0, 0.0], [1.0, 0.0]), LineSegment([1.0, 0.0], [1.0, 1.0]),
               LineSegment([1.0, 1.0], [0.0, 1.0]), LineSegment([0.0, 1.0], [0.0, 0.0])]
bc_triangle() = [LineSegment([0.0, 0.0], [1.0, 0.0]), LineSegment([1.0, 0.0], [0.0, 1.0]), LineSegment([0.0, 1.0], [0.0, 0.0])]
bc_half_disk() = AbsCurve[CircleSegment(1.0, pi, 0.0, [0.0, 0.0]), LineSegment([-1.0, 0.0], [1.0, 0.0])]

@testset "component_lengths" begin
    # two perpendicular unit line segments forming a right-angle corner
    line1 = LineSegment([0.0, 0.0], [1.0, 0.0])
    line2 = LineSegment([1.0, 0.0], [1.0, 1.0])
    comp = [line1, line2]
    lens, cum, Ltot = component_lengths(comp)
    @test lens == [1.0, 1.0]
    @test cum == [0.0, 1.0, 2.0]
    @test Ltot == 2.0

    # different lengths: cumulative sums, cum[1] = 0 and Ltot = cum[end]
    lens3, cum3, Ltot3 = component_lengths([LineSegment([0.0, 0.0], [3.0, 4.0]), LineSegment([3.0, 4.0], [3.0, 5.0])])
    @test lens3 == [5.0, 1.0]
    @test cum3 == [0.0, 5.0, 6.0]
    @test Ltot3 == cum3[end] == 6.0

    # the half disk has length pi + 2
    _, _, Lhd = component_lengths(bc_half_disk())
    @test Lhd ≈ pi + 2 atol = bc_ATOL

    # a single curve and the element type
    l1, c1, L1 = component_lengths([line1])
    @test l1 == [1.0] && c1 == [0.0, 1.0] && L1 == 1.0
    @test @inferred(component_lengths([LineSegment(Float32[0, 0], Float32[1, 0])]))[3] === 1.0f0
end

@testset "_junction_angle and _is_true_corner" begin
    line1 = LineSegment([0.0, 0.0], [1.0, 0.0])
    line2 = LineSegment([1.0, 0.0], [1.0, 1.0])
    line3 = LineSegment([1.0, 0.0], [2.0, 0.0])
    @test _is_true_corner(line1, line2, Float64)
    # two colinear segments have a smooth join (no corner)
    @test !_is_true_corner(line1, line3, Float64)

    # unsigned turning angle in [0, pi]
    @test _junction_angle(line1, line2, Float64) ≈ pi / 2 atol = bc_ATOL
    @test _junction_angle(line1, line3, Float64) == 0.0
    # a reversal of direction is a turning angle of pi, a left and a right turn are not distinguished
    back = LineSegment([1.0, 0.0], [0.0, 0.0])
    right = LineSegment([1.0, 0.0], [1.0, -1.0])
    @test _junction_angle(line1, back, Float64) ≈ pi atol = bc_ATOL
    @test _junction_angle(line1, right, Float64) ≈ pi / 2 atol = bc_ATOL

    # a tangent-continuous join of a line and an arc is smooth: the arc starts at (1,0) with direction (0,1)
    arc = CircleSegment(1.0, pi, 0.0, [0.0, 0.0])
    up = LineSegment([1.0, -1.0], [1.0, 0.0])
    @test !_is_true_corner(up, arc, Float64)
    # the angle_tol keyword turns a small turning angle into a smooth join
    @test _is_true_corner(line1, line2, Float64; angle_tol=1.0)
    @test !_is_true_corner(line1, line2, Float64; angle_tol=2.0)
end

@testset "_component_corner_locations" begin
    line1 = LineSegment([0.0, 0.0], [1.0, 0.0])
    line2 = LineSegment([1.0, 0.0], [1.0, 1.0])
    corners = _component_corner_locations(Float64, [line1, line2])
    @test length(corners) == 2
    @test 0.0 in corners
    @test isapprox(maximum(corners), pi; atol=1e-10)

    # two colinear segments have a smooth join (no corner)
    line3 = LineSegment([1.0, 0.0], [2.0, 0.0])
    @test isempty(_component_corner_locations(Float64, [line1, line3]))

    # the square has four corners at the cumulative quarter positions of the global parameter, seam corner first
    @test _component_corner_locations(Float64, bc_square()) ≈ [0.0, pi / 2, pi, 3 * pi / 2] atol = bc_ATOL
    # the half disk has corners at the seam and after the arc, t = 2 pi * pi / (pi + 2)
    @test _component_corner_locations(Float64, bc_half_disk()) ≈ [0.0, 2 * pi * pi / (pi + 2)] atol = bc_ATOL
    # a full circle is smooth everywhere, including its periodic seam
    @test isempty(_component_corner_locations(Float64, [CircleSegment(1.0, 2 * pi, 0.0, [0.0, 0.0])]))
end

@testset "_component_corner_angles" begin
    # interior angles of a convex counterclockwise component: pi minus the turning angle, seam corner first
    @test _component_corner_angles(Float64, bc_square()) ≈ fill(pi / 2, 4) atol = bc_ATOL
    @test _component_corner_angles(Float64, bc_triangle()) ≈ [pi / 2, pi / 4, pi / 4] atol = bc_ATOL
    @test sum(_component_corner_angles(Float64, bc_triangle())) ≈ pi atol = bc_ATOL
    @test sum(_component_corner_angles(Float64, bc_square())) ≈ 2 * pi atol = bc_ATOL
    @test _component_corner_angles(Float64, bc_half_disk()) ≈ [pi / 2, pi / 2] atol = bc_ATOL

    # the angles of a Polygon domain are in the same convention and vertex order as the corners of its boundary
    poly = Polygon([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], 1)
    @test _component_corner_angles(Float64, poly.boundary) ≈ poly.angles atol = bc_ATOL

    # smooth joins do not produce corners: the unit square split into collinear pieces has the same four angles
    split_square = [LineSegment([0.0, 0.0], [0.5, 0.0]), LineSegment([0.5, 0.0], [1.0, 0.0]),
                    LineSegment([1.0, 0.0], [1.0, 1.0]), LineSegment([1.0, 1.0], [0.0, 1.0]), LineSegment([0.0, 1.0], [0.0, 0.0])]
    @test _component_corner_angles(Float64, split_square) ≈ fill(pi / 2, 4) atol = bc_ATOL
    @test length(_component_corner_angles(Float64, split_square)) == length(_component_corner_locations(Float64, split_square))
    @test isempty(_component_corner_angles(Float64, [CircleSegment(1.0, 2 * pi, 0.0, [0.0, 0.0])]))
end

@testset "_signed_junction_angle" begin
    line1 = LineSegment([0.0, 0.0], [1.0, 0.0])
    up = LineSegment([1.0, 0.0], [1.0, 1.0])
    down = LineSegment([1.0, 0.0], [1.0, -1.0])
    straight = LineSegment([1.0, 0.0], [2.0, 0.0])
    back = LineSegment([1.0, 0.0], [0.0, 0.0])
    # signed turning angle: positive for a left turn, negative for a right turn, zero when straight, pi for a reversal
    @test _signed_junction_angle(line1, up, Float64) ≈ pi / 2 atol = bc_ATOL
    @test _signed_junction_angle(line1, down, Float64) ≈ -pi / 2 atol = bc_ATOL
    @test _signed_junction_angle(line1, straight, Float64) == 0.0
    @test _signed_junction_angle(line1, back, Float64) ≈ pi atol = bc_ATOL
    # its absolute value is the unsigned _junction_angle
    @test abs(_signed_junction_angle(line1, down, Float64)) ≈ _junction_angle(line1, down, Float64) atol = bc_ATOL
end

@testset "_component_corner_angles (hole)" begin
    # hole = true: the angle of the region outside the loop, pi + signed turning angle, seam corner first
    @test _component_corner_angles(Float64, bc_square(); hole = true) ≈ fill(3 * pi / 2, 4) atol = bc_ATOL
    # right triangle (0,0), (1,0), (0,1): turning angles pi/2, 3 pi / 4, 3 pi / 4 in seam-first order
    @test _component_corner_angles(Float64, bc_triangle(); hole = true) ≈ [3 * pi / 2, 7 * pi / 4, 7 * pi / 4] atol = bc_ATOL
    # for a convex loop traversed counterclockwise the hole angle is 2 pi minus the interior angle
    @test _component_corner_angles(Float64, bc_triangle(); hole = true) ≈ 2 * pi .- _component_corner_angles(Float64, bc_triangle()) atol = bc_ATOL
    # a smooth loop has no corners, smooth joins of collinear pieces do not count
    @test isempty(_component_corner_angles(Float64, [CircleSegment(1.0, 2 * pi, 0.0, [0.0, 0.0])]; hole = true))
    # the angle_tol keyword is honoured
    @test isempty(_component_corner_angles(Float64, bc_square(); angle_tol = 2.0, hole = true))
    # hole = false is the default
    @test _component_corner_angles(Float64, bc_square(); hole = false) == _component_corner_angles(Float64, bc_square())
end

@testset "print_component_junctions" begin
    # one @info record per join
    logs, _ = collect_test_logs() do
        print_component_junctions(bc_square())
    end
    @test length(logs) == 4
    @test all(l -> l.level == Base.CoreLogging.Info && l.message == "Junction", logs)
    kw(l, k) = Dict(l.kwargs)[k]
    @test count(l -> kw(l, :classification) == "true corner", logs) == 4
    @test count(l -> kw(l, :classification) == "smooth join", logs) == 0
    @test [(kw(l, :from), kw(l, :to)) for l in logs] == [(1, 2), (2, 3), (3, 4), (4, 1)]
    # the last join closes the loop at the seam sigma = 0
    @test kw(logs[4], :σ) == 0.0
    @test kw(logs[1], :angle) ≈ pi / 2 atol = bc_ATOL
    # the parameter of the first join is 2 pi * 1/4
    @test kw(logs[1], :σ) ≈ pi / 2 atol = bc_ATOL

    # smooth joins are labelled as such; the keyword T selects the number type of the printed values
    smooth = [LineSegment([0.0, 0.0], [1.0, 0.0]), LineSegment([1.0, 0.0], [2.0, 0.0])]
    smooth_logs, _ = collect_test_logs() do
        print_component_junctions(smooth; T=Float32)
    end
    @test length(smooth_logs) == 2
    @test count(l -> kw(l, :classification) == "smooth join", smooth_logs) == 2
    @test count(l -> kw(l, :classification) == "true corner", smooth_logs) == 0
    @test all(l -> kw(l, :angle) isa Float32, smooth_logs)
    @test all(l -> kw(l, :σ) isa Float32, smooth_logs)
    # a large angle_tol makes every join smooth
    tol_logs, _ = collect_test_logs() do
        print_component_junctions(bc_square(); angle_tol=2.0)
    end
    @test length(tol_logs) == 4
    @test count(l -> kw(l, :classification) == "smooth join", tol_logs) == 4
end

@testset "_global_t_to_segment_u" begin
    comp = bc_square()   # four unit sides, global parameter t = 2 pi s / 4
    @test _global_t_to_segment_u(Float64, comp, 0.0) == (1, 0.0)
    @test _global_t_to_segment_u(Float64, comp, pi / 4) == (1, 0.5)
    # t = pi / 2 is the start of the second side, t = 3 pi / 2 + pi / 4 is the middle of the fourth
    @test _global_t_to_segment_u(Float64, comp, pi / 2) == (2, 0.0)
    j, u = _global_t_to_segment_u(Float64, comp, 3 * pi / 2 + pi / 4)
    @test j == 4
    @test u ≈ 0.5 atol = bc_ATOL
    # the end of the period wraps to the start of the first segment
    @test _global_t_to_segment_u(Float64, comp, 2 * pi) == (1, 0.0)

    # unequal lengths: sides of length 1 and 3, so the second side starts at t = pi / 2 and t = pi is at u = 1/3
    uneven = [LineSegment([0.0, 0.0], [1.0, 0.0]), LineSegment([1.0, 0.0], [1.0, 3.0])]
    @test _global_t_to_segment_u(Float64, uneven, pi / 4) == (1, 0.5)
    j2, u2 = _global_t_to_segment_u(Float64, uneven, Float64(pi))
    @test j2 == 2
    @test u2 ≈ 1 / 3 atol = bc_ATOL

    # the mapping returns the point of the curve at the corresponding arc length
    for t in (0.1, 1.0, 2.5, 4.0, 6.0)
        jt, ut = _global_t_to_segment_u(Float64, comp, t)
        s = t / (2 * pi) * 4
        @test jt == clamp(floor(Int, s) + 1, 1, 4)
        @test ut ≈ s - floor(s) atol = bc_ATOL
    end
    @test @inferred(_global_t_to_segment_u(Float64, comp, 1.0)) isa Tuple{Int,Float64}
    # Float32 components: the middle of the second of two unit sides (total length 2, so t = 3 pi / 2 is s = 1.5)
    comp32 = [LineSegment(Float32[0, 0], Float32[1, 0]), LineSegment(Float32[1, 0], Float32[1, 1])]
    j32, u32 = _global_t_to_segment_u(Float32, comp32, Float32(3 * pi / 2))
    @test j32 == 2
    @test u32 isa Float32
    @test u32 ≈ 0.5f0 atol = bc_ATOL_F32
end

@testset "_global_t_to_segment_u (precomputed component_lengths)" begin
    # the method taking (lens, cum, Ltot) is the comp method with the lengths hoisted out of the loop: results are identical
    for comp in (bc_square(), bc_triangle(), bc_half_disk())
        lens, cum, Ltot = component_lengths(comp)
        tj = 2 * pi * cum[2] / Ltot                                      # end of the first segment (a junction)
        for t in (0.0, 1e-9, 0.3, tj, 2.0, Float64(pi), 4.4, 2 * pi - 1e-9, 2 * pi)
            @test _global_t_to_segment_u(Float64, lens, cum, Ltot, t) == _global_t_to_segment_u(Float64, comp, t)
        end
    end

    # a junction (t = pi / 2 of the square) starts the next segment, t = 0 and t = 2 pi the first one
    lens, cum, Ltot = component_lengths(bc_square())
    @test _global_t_to_segment_u(Float64, lens, cum, Ltot, 0.0) == (1, 0.0)
    @test _global_t_to_segment_u(Float64, lens, cum, Ltot, pi / 2) == (2, 0.0)
    @test _global_t_to_segment_u(Float64, lens, cum, Ltot, 2 * pi) == (1, 0.0)

    # Float32
    comp32 = [LineSegment(Float32[0, 0], Float32[1, 0]), LineSegment(Float32[1, 0], Float32[1, 1])]
    lens32, cum32, Ltot32 = component_lengths(comp32)
    for t in (0.0f0, 1.0f0, Float32(pi), Float32(3 * pi / 2), 2.0f0 * Float32(pi))
        @test _global_t_to_segment_u(Float32, lens32, cum32, Ltot32, t) == _global_t_to_segment_u(Float32, comp32, t)
    end
    @test @inferred(_global_t_to_segment_u(Float64, lens, cum, Ltot, 1.0)) isa Tuple{Int,Float64}
end

@testset "_eval_composite_geom_global_t" begin
    # polygon: points of the sides, derivatives with respect to t are the chord times du/dt = L / (2 pi l) = 2 / pi, no curvature
    comp = bc_square()
    xy, γt, γtt = _eval_composite_geom_global_t(Float64, comp, pi / 4)
    @test xy ≈ SVector(0.5, 0.0) atol = bc_ATOL
    @test γt ≈ SVector(2 / pi, 0.0) atol = bc_ATOL
    @test γtt == SVector(0.0, 0.0)
    xy2, γt2, _ = _eval_composite_geom_global_t(Float64, comp, 3 * pi / 4)
    @test xy2 ≈ SVector(1.0, 0.5) atol = bc_ATOL
    @test γt2 ≈ SVector(0.0, 2 / pi) atol = bc_ATOL

    # a unit-speed-in-t circle: the full circle of radius R is R (cos t, sin t), so gamma_t = R (-sin t, cos t) and
    # gamma_tt = -R (cos t, sin t)
    R = 1.7
    circle = [CircleSegment(R, 2 * pi, 0.0, [0.0, 0.0])]
    for t in (0.3, 1.9, 4.4)
        xyc, γtc, γttc = _eval_composite_geom_global_t(Float64, circle, t)
        @test xyc ≈ R * SVector(cos(t), sin(t)) atol = bc_ATOL
        @test γtc ≈ R * SVector(-sin(t), cos(t)) atol = bc_ATOL
        @test γttc ≈ -R * SVector(cos(t), sin(t)) atol = bc_ATOL
    end

    # the point agrees with `curve` of the active segment at the local parameter, also across the half-disk junction
    hd = bc_half_disk()
    tj = 2 * pi * pi / (pi + 2)   # global parameter of the junction between arc and diameter
    xya, _, _ = _eval_composite_geom_global_t(Float64, hd, tj - 1e-9)
    xyb, _, _ = _eval_composite_geom_global_t(Float64, hd, tj + 1e-9)
    @test xya ≈ SVector(-1.0, 0.0) atol = 1e-7
    @test xyb ≈ SVector(-1.0, 0.0) atol = 1e-7
    xymid, _, _ = _eval_composite_geom_global_t(Float64, hd, 2 * pi * (pi + 1) / (pi + 2))   # middle of the diameter
    @test xymid ≈ SVector(0.0, 0.0) atol = bc_ATOL
end

@testset "_eval_composite_geom_global_t (precomputed component_lengths)" begin
    # the method taking (lens, cum, Ltot) returns exactly what the comp method returns: point and both derivatives
    for comp in (bc_square(), bc_half_disk(), [CircleSegment(1.7, 2 * pi, 0.0, [0.0, 0.0])])
        lens, cum, Ltot = component_lengths(comp)
        tj = 2 * pi * cum[2] / Ltot                                        # end of the first segment (a junction; the seam for the circle)
        for t in (0.0, 1e-9, 0.3, tj, 2.0, Float64(pi), 4.4, 2 * pi - 1e-9, 2 * pi)
            @test _eval_composite_geom_global_t(Float64, comp, lens, cum, Ltot, t) == _eval_composite_geom_global_t(Float64, comp, t)
        end
    end

    # the analytic square values also hold for the hoisted method: t = pi / 4 is the middle of the first side
    sq = bc_square()
    lens, cum, Ltot = component_lengths(sq)
    xy, γt, γtt = _eval_composite_geom_global_t(Float64, sq, lens, cum, Ltot, pi / 4)
    @test xy ≈ SVector(0.5, 0.0) atol = bc_ATOL
    @test γt ≈ SVector(2 / pi, 0.0) atol = bc_ATOL
    @test γtt == SVector(0.0, 0.0)
    # the junction belongs to the next side, the end of the period wraps to the start of the first side
    @test _eval_composite_geom_global_t(Float64, sq, lens, cum, Ltot, pi / 2)[1] ≈ SVector(1.0, 0.0) atol = bc_ATOL
    @test _eval_composite_geom_global_t(Float64, sq, lens, cum, Ltot, 2 * pi)[1] ≈ SVector(0.0, 0.0) atol = bc_ATOL

    # Float32
    comp32 = [LineSegment(Float32[0, 0], Float32[1, 0]), LineSegment(Float32[1, 0], Float32[1, 1])]
    lens32, cum32, Ltot32 = component_lengths(comp32)
    for t in (0.0f0, 1.0f0, Float32(3 * pi / 2), 2.0f0 * Float32(pi))
        @test _eval_composite_geom_global_t(Float32, comp32, lens32, cum32, Ltot32, t) == _eval_composite_geom_global_t(Float32, comp32, t)
    end
end

@testset "_group_curves_by_domain_id" begin
    # groups a flat physical-boundary curve list into connected components by curve `domain_id`, preserving first-seen order:
    # curves with distinct `domain_id`s become distinct components, curves sharing one `domain_id` stay in a single component.
    line1b = LineSegment([0.0, 0.0], [1.0, 0.0]; domain_id=1)
    line2b = LineSegment([1.0, 0.0], [1.0, 1.0]; domain_id=2)
    comps = _group_curves_by_domain_id(AbsCurve[line1b, line2b])
    @test length(comps) == 2
    @test comps[1] == [line1b]
    @test comps[2] == [line2b]

    line1 = LineSegment([0.0, 0.0], [1.0, 0.0])
    line2 = LineSegment([1.0, 0.0], [1.0, 1.0])
    comps_same = _group_curves_by_domain_id(AbsCurve[line1, line2])
    @test length(comps_same) == 1
    @test comps_same[1] == [line1, line2]

    # first-seen order of the ids is kept, interleaved curves are collected in their group
    a = LineSegment([0.0, 0.0], [1.0, 0.0]; domain_id=5, segment_id=1)
    b = LineSegment([1.0, 0.0], [1.0, 1.0]; domain_id=3, segment_id=1)
    c = LineSegment([1.0, 1.0], [0.0, 1.0]; domain_id=5, segment_id=2)
    ordered = _group_curves_by_domain_id(AbsCurve[a, b, c])
    @test ordered == [AbsCurve[a, c], AbsCurve[b]]
    @test isempty(_group_curves_by_domain_id(AbsCurve[]))
end
