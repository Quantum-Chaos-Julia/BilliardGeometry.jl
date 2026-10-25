# Unit tests for the Poincare-Birkhoff coordinates (arc length s along the boundary, sine of the incidence angle p):
#   src/geometry/poincarebirkhoff.jl   (PoincareBirkhoff, fundamental_pb_coords, get_pb_curve, pb_coords, pb_sectors)
# Reference objects: one LineSegment and one half-circle CircleSegment for the one-curve coordinates, and the
# RectangleBilliard(1.0, 0.6) (fundamental boundary of two lines of length 0.6 and 1.0 with registry Y, XY, X). On a curve
# the coordinate s of a point is the previous boundary length plus the arc length to the point, and |p| is the sine of the angle
# between the velocity and the boundary normal: p = 0 along the normal and |p| = 1 along the tangent. Image sectors follow
# apply_symmetry_pb: reflection sectors run backwards (s -> (k + 1) L - s, p -> -p), XY continues the traversal (s -> k L + s).
# Billiards with holes are not used. No eigen solves, no pinned wavefunction values. Helper functions and constants
# use the prefix `pb_`.

using StaticArrays: SVector
using BilliardGeometry: fundamental_pb_coords, signed_angle

# Absolute tolerance for closed-form comparisons (arc length integration is adaptive, so slightly looser than ulps).
const pb_ATOL = 1e-8

# Bottom edge of the square [0, 2]^2 traversed counterclockwise, and the upper half of the unit circle.
pb_line() = LineSegment([0.0, 0.0], [2.0, 0.0])
pb_half_circle() = CircleSegment(1.0, pi, 0.0, [0.0, 0.0])

# Unit vector rotated by `alpha` from `u`.
pb_rotate(u, alpha) = SVector(cos(alpha) * u[1] - sin(alpha) * u[2], sin(alpha) * u[1] + cos(alpha) * u[2])

@testset "PoincareBirkhoff" begin
    pb = PoincareBirkhoff(0.25, -0.5)
    @test pb.s == 0.25 && pb.p == -0.5
    @test pb isa AbsCoords
    @test PoincareBirkhoff(1.0f0, 0.5f0) isa PoincareBirkhoff{Float32}
end

@testset "fundamental_pb_coords (LineSegment)" begin
    line = pb_line()
    pt = SVector(0.5, 0.0)
    n = domain_gradient_vector(line, pt)
    n = n ./ sqrt(sum(abs2, n))                              # unit normal of the line
    tang = SVector(1.0, 0.0)

    # s is the length of the previous curves plus the arc length from the start of this curve
    s, _ = fundamental_pb_coords(3.0, line, pt, n)
    @test s ≈ 3.5 atol = pb_ATOL
    @test fundamental_pb_coords(0.0, line, SVector(1.5, 0.0), n)[1] ≈ 1.5 atol = pb_ATOL

    # velocity along the normal: p = 0, independent of its length
    @test fundamental_pb_coords(0.0, line, pt, n)[2] ≈ 0.0 atol = pb_ATOL
    @test fundamental_pb_coords(0.0, line, pt, 4.0 * n)[2] ≈ 0.0 atol = pb_ATOL
    # velocity along the tangent: |p| = 1, and reversing the tangent flips the sign
    p_fwd = fundamental_pb_coords(0.0, line, pt, tang)[2]
    p_bwd = fundamental_pb_coords(0.0, line, pt, -tang)[2]
    @test abs(p_fwd) ≈ 1.0 atol = pb_ATOL
    @test p_bwd ≈ -p_fwd atol = pb_ATOL
    # |p| is the sine of the angle to the normal, and p is odd in that angle
    for alpha in (0.2, 0.7, 1.2)
        p_plus = fundamental_pb_coords(0.0, line, pt, pb_rotate(n, alpha))[2]
        p_minus = fundamental_pb_coords(0.0, line, pt, pb_rotate(n, -alpha))[2]
        @test abs(p_plus) ≈ sin(alpha) atol = pb_ATOL
        @test p_minus ≈ -p_plus atol = pb_ATOL
        @test p_plus ≈ sin(signed_angle(pb_rotate(n, alpha), n)) atol = pb_ATOL     # the defining formula
    end
    # only the direction of the velocity matters
    @test fundamental_pb_coords(0.0, line, pt, 7.0 * pb_rotate(n, 0.7))[2] ≈ fundamental_pb_coords(0.0, line, pt, pb_rotate(n, 0.7))[2] atol = pb_ATOL
end

@testset "fundamental_pb_coords (CircleSegment)" begin
    arc = pb_half_circle()
    pt = SVector(0.0, 1.0)                                    # a quarter of the way around, arc length pi/2
    n = domain_gradient_vector(arc, pt)
    n = n ./ sqrt(sum(abs2, n))
    @test abs(n[2]) ≈ 1.0 atol = pb_ATOL                      # the normal at the top of the circle is vertical (radial)

    @test fundamental_pb_coords(0.0, arc, pt, n)[1] ≈ pi / 2 atol = pb_ATOL
    @test fundamental_pb_coords(2.0, arc, pt, n)[1] ≈ 2.0 + pi / 2 atol = pb_ATOL
    @test fundamental_pb_coords(0.0, arc, pt, n)[2] ≈ 0.0 atol = pb_ATOL
    @test abs(fundamental_pb_coords(0.0, arc, pt, SVector(1.0, 0.0))[2]) ≈ 1.0 atol = pb_ATOL
    for alpha in (0.3, 0.9)
        @test abs(fundamental_pb_coords(0.0, arc, pt, pb_rotate(n, alpha))[2]) ≈ sin(alpha) atol = pb_ATOL
    end
end

@testset "get_pb_curve" begin
    rect = RectangleBilliard(1.0, 0.6)
    P = get_boundary_curves(rect)
    comp = CompositeCurve(P)
    # the returned offset is the length of all curves before the requested one
    l1, c1 = get_pb_curve(comp, P[1].segment_id, P[1].domain_id)
    l2, c2 = get_pb_curve(comp, P[2].segment_id, P[2].domain_id)
    @test l1 == 0.0
    @test l2 ≈ P[1].length atol = pb_ATOL
    @test c1.segment_id == P[1].segment_id && c2.segment_id == P[2].segment_id
    @test c2.length ≈ P[2].length atol = pb_ATOL
    # an unknown segment or domain id finds no curve
    @test get_pb_curve(comp, 99, P[1].domain_id) === nothing
    @test get_pb_curve(comp, P[1].segment_id, 99) === nothing
end

@testset "pb_coords" begin
    rect = RectangleBilliard(1.0, 0.6)
    nsym = length(rect.symmetries)
    crv = first(get_boundary_curves(rect))
    L = CompositeCurve(get_boundary_curves(rect)).length
    pt = curve(crv, 0.5)
    tang = tangent(crv, 0.5)
    vel = SVector(tang[2], -tang[1]) + tang                   # 45 degrees from the normal: |p| = sin(pi/4)

    base = pb_coords(rect, crv.segment_id, crv.domain_id, 1, pt, vel)
    @test base isa PoincareBirkhoff
    @test base.s ≈ 0.3 atol = pb_ATOL                         # midpoint of the first curve of length 0.6
    @test abs(base.p) ≈ sin(pi / 4) atol = pb_ATOL

    # the sector images follow apply_symmetry_pb, and each stays inside the window of its own copy
    for sym_sector in 1:(nsym + 1)
        coords = pb_coords(rect, crv.segment_id, crv.domain_id, sym_sector, pt, vel)
        lo, hi = (sym_sector - 1) * L, sym_sector * L
        @test lo - pb_ATOL <= coords.s <= hi + pb_ATOL
    end
    # sector 2 is the Y reflection (reversing), sector 3 is XY (continuing), sector 4 is the X reflection (reversing)
    c2 = pb_coords(rect, crv.segment_id, crv.domain_id, 2, pt, vel)
    c3 = pb_coords(rect, crv.segment_id, crv.domain_id, 3, pt, vel)
    c4 = pb_coords(rect, crv.segment_id, crv.domain_id, 4, pt, vel)
    @test c2.s ≈ 2 * L - base.s atol = pb_ATOL
    @test c2.p ≈ -base.p atol = pb_ATOL
    @test c3.s ≈ 2 * L + base.s atol = pb_ATOL
    @test c3.p ≈ base.p atol = pb_ATOL
    @test c4.s ≈ 4 * L - base.s atol = pb_ATOL
    @test c4.p ≈ -base.p atol = pb_ATOL
end

@testset "pb_sectors" begin
    # rectangle: the fundamental breakpoints 0, 0.6, 1.6 and the images under Y (reversed), XY and X (reversed)
    rect = RectangleBilliard(1.0, 0.6)
    L = 1.6
    sectors = pb_sectors(rect)
    @test issorted(sectors)
    @test length(sectors) == length(unique(sectors))
    @test isapprox(first(sectors), 0.0; atol = pb_ATOL)
    @test isapprox(last(sectors), (length(rect.symmetries) + 1) * L; atol = pb_ATOL)
    @test sectors ≈ [0.0, 0.6, 1.6, 2.6, 3.2, 3.8, 4.8, 5.8, 6.4] atol = pb_ATOL
    # the sectors are symmetric about the middle of the full boundary
    @test sectors ≈ (4 * L) .- reverse(sectors) atol = pb_ATOL

    # no symmetry: the breakpoints are the cumulative curve lengths
    tri = TriangleBilliard(1.0, 1.0)
    _, cum, Ltot = component_lengths(get_boundary_curves(tri))
    @test pb_sectors(tri) ≈ cum atol = pb_ATOL
    @test last(pb_sectors(tri)) ≈ Ltot atol = pb_ATOL
end
