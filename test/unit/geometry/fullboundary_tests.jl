# Unit tests for the reconstruction of the complete physical boundary from the fundamental domain and its symmetries:
#   src/geometry/fullboundary.jl   (full_boundary, _full_boundary_data, _apply_symmetry_to_curve, _polar_symmetry_image,
#                                   _reverse_curve)
# Reference objects: one LineSegment, one CircleSegment, one FourierCoeffPolarSegment and one PolarSegment, each placed off
# the origin and off the axes so that no symmetry maps the curve onto itself. The defining identity of an image is
#   curve(_apply_symmetry_to_curve(sym, c), t) == apply_symmetry(sym, curve(c, t)),
# and a reversed curve runs backwards over the same trace, curve(_reverse_curve(c), t) == curve(c, 1 - t). The provenance
# of the reconstructed boundary of the rectangle (sectors identity, Y, XY, X) and of the C3 billiard is known from the
# construction rules: reflections reverse the curve order and each curve, rotations and the pi rotation XY do neither.
# A polar image reuses the radial function of its source, so it is the exact image only under a symmetry that leaves the
# radial function invariant; the polar fixtures are therefore D4-symmetric (cos 4 phi) and are only mapped by the five
# reflections and the rotations of order 4, and the C3 fundamental curve only by rotations of order 3.
# The orientation field of line images keeps domain_fun consistent with the source curve (G13). Closure of the
# reconstructed loop and its total length on real billiards are covered by the integration tier. No eigen solves,
# no pinned wavefunction values. Helper functions and
# constants use the prefix `fb_`.

using StaticArrays: SVector
using BilliardGeometry: _apply_symmetry_to_curve, _reverse_curve, _full_boundary_data

# Absolute tolerance for pointwise comparisons of curve images (a few ulps of Float64 after one trigonometric evaluation).
const fb_ATOL = 1e-12

# Parameter values at which curve images are compared.
const fb_TS = (0.0, 0.2, 0.5, 0.8, 1.0)

# All reflection generators and two rotations (any symmetry maps a line or a circular arc to a line or arc).
const fb_SYMMETRIES = (XAxisReflection(), YAxisReflection(), XYAxisReflection(), DiagonalReflection(),
                       AntiDiagonalReflection(), NFoldRotation(3, 1), NFoldRotation(4, 3))
# The symmetries of a D4-symmetric radial function.
const fb_D4_SYMMETRIES = (XAxisReflection(), YAxisReflection(), XYAxisReflection(), DiagonalReflection(),
                          AntiDiagonalReflection(), NFoldRotation(4, 1), NFoldRotation(4, 3))
# The orientation-reversing subset, whose images are reversed by `full_boundary`.
const fb_REFLECTIONS = (XAxisReflection(), YAxisReflection(), DiagonalReflection(), AntiDiagonalReflection())

# Off-axis curves with non-default metadata, so that the transfer of domain_id and segment_id is visible.
# The Fourier coefficients are ordered (sin phi, cos phi, sin 2 phi, cos 2 phi, ...), so the eighth entry is cos 4 phi.
fb_line() = LineSegment([0.3, 0.1], [1.2, 0.7]; domain_id = 2, segment_id = 5)
fb_arc() = CircleSegment(0.8, 1.1, 0.4, [0.5, 0.2]; domain_id = 2, segment_id = 5)
fb_fourier() = FourierCoeffPolarSegment([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3]; arc_angle = 1.3, shift_angle = 0.4, center = [0.3, -0.2],
                                        domain_id = 2, segment_id = 5)
fb_polar() = PolarSegment(φ -> 1.0 + 0.2 * cos(4 * φ); arc_angle = 1.3, shift_angle = 0.4, center = [0.3, -0.2],
                          domain_id = 2, segment_id = 5)

# `curve(img, t) == apply_symmetry(sym, curve(c, t))` on the grid `fb_TS`.
function fb_image_matches(sym, c, img)
    return all(t -> isapprox(curve(img, t), apply_symmetry(sym, curve(c, t)); atol = fb_ATOL), fb_TS)
end

# `curve(rc, t) == curve(c, 1 - t)` on the grid `fb_TS`.
function fb_reversal_matches(c, rc)
    return all(t -> isapprox(curve(rc, t), curve(c, 1 - t); atol = fb_ATOL), fb_TS)
end

@testset "_apply_symmetry_to_curve (LineSegment)" begin
    c = fb_line()
    for sym in fb_SYMMETRIES
        img = _apply_symmetry_to_curve(sym, c)
        @test img isa LineSegment
        @test fb_image_matches(sym, c, img)
        @test isapprox(img.pt0, apply_symmetry(sym, c.pt0); atol = fb_ATOL)
        @test isapprox(img.pt1, apply_symmetry(sym, c.pt1); atol = fb_ATOL)
        @test isapprox(img.length, c.length; atol = fb_ATOL)
        @test img.domain_id == c.domain_id && img.segment_id == c.segment_id
        @test img.bc == c.bc
    end
    # the image under X is the mirror image about the x axis
    img = _apply_symmetry_to_curve(XAxisReflection(), c)
    @test img.pt0 == [0.3, -0.1] && img.pt1 == [1.2, -0.7]
end

@testset "_apply_symmetry_to_curve (LineSegment orientation)" begin
    # G13: orientation-preserving images (det > 0) keep the orientation field, reflections flip it
    for o in (1, -1)
        c = LineSegment([0.3, 0.1], [1.2, 0.7]; orientation = o)
        for sym in (XYAxisReflection(), NFoldRotation(3, 1), NFoldRotation(4, 3))
            @test _apply_symmetry_to_curve(sym, c).orientation == o
        end
        for sym in fb_REFLECTIONS
            @test _apply_symmetry_to_curve(sym, c).orientation == -o
        end
    end
    # domain_fun of the image at the image of a point equals domain_fun of the source at the point
    pts = (SVector(0.1, 0.9), SVector(0.8, -0.3), SVector(-0.4, 0.2))
    for o in (1, -1), sym in fb_SYMMETRIES
        c = LineSegment([0.3, 0.1], [1.2, 0.7]; orientation = o)
        img = _apply_symmetry_to_curve(sym, c)
        rimg = _reverse_curve(img)
        for p in pts
            @test isapprox(domain_fun(img, apply_symmetry(sym, p)), domain_fun(c, p); atol = fb_ATOL)
            @test isapprox(domain_fun(rimg, apply_symmetry(sym, p)), domain_fun(c, p); atol = fb_ATOL)
        end
    end
end

@testset "_apply_symmetry_to_curve (CircleSegment)" begin
    c = fb_arc()
    for sym in fb_SYMMETRIES
        img = _apply_symmetry_to_curve(sym, c)
        @test img isa CircleSegment
        @test fb_image_matches(sym, c, img)
        @test isapprox(img.radius, c.radius; atol = fb_ATOL)
        @test isapprox(img.length, c.length; atol = fb_ATOL)
        @test img.domain_id == c.domain_id && img.segment_id == c.segment_id
    end
    # a reflection reverses the sense of the arc, a rotation keeps it, and the center is mapped by the symmetry
    @test _apply_symmetry_to_curve(XAxisReflection(), c).arc_angle == -c.arc_angle
    @test _apply_symmetry_to_curve(NFoldRotation(3, 1), c).arc_angle == c.arc_angle
    @test isapprox(_apply_symmetry_to_curve(XAxisReflection(), c).center, [0.5, -0.2]; atol = fb_ATOL)
    @test isapprox(_apply_symmetry_to_curve(XYAxisReflection(), c).center, [-0.5, -0.2]; atol = fb_ATOL)
end

@testset "_apply_symmetry_to_curve (FourierCoeffPolarSegment)" begin
    # default full-circle curve of a polar billiard, as used for the D2 reconstruction
    orig = PolarBilliard([0.0, 0.0, 0.0, 0.3]).fundamental_domain.boundary[1]
    for sym in (XAxisReflection(), YAxisReflection(), XYAxisReflection())
        img = _apply_symmetry_to_curve(sym, orig)
        @test img isa FourierCoeffPolarSegment
        @test fb_image_matches(sym, orig, img)
    end

    # off-center, shifted arc under every symmetry of its radial function
    c = fb_fourier()
    for sym in fb_D4_SYMMETRIES
        img = _apply_symmetry_to_curve(sym, c)
        @test img isa FourierCoeffPolarSegment
        @test fb_image_matches(sym, c, img)
        @test img.coef == c.coef
        @test img.domain_id == c.domain_id && img.segment_id == c.segment_id
    end
    @test _apply_symmetry_to_curve(XAxisReflection(), c).arc_angle == -c.arc_angle
    @test _apply_symmetry_to_curve(NFoldRotation(4, 1), c).arc_angle == c.arc_angle

    # C3 fundamental curve under its two rotations
    c3_curve = first(get_boundary_curves(C3Billiard(0.2)))
    for sym in (NFoldRotation(3, 1), NFoldRotation(3, 2))
        img = _apply_symmetry_to_curve(sym, c3_curve)
        @test img isa FourierCoeffPolarSegment
        @test fb_image_matches(sym, c3_curve, img)
    end
end

@testset "_apply_symmetry_to_curve (PolarSegment)" begin
    c = fb_polar()
    for sym in fb_D4_SYMMETRIES
        img = _apply_symmetry_to_curve(sym, c)
        @test img isa PolarSegment
        @test fb_image_matches(sym, c, img)
        @test img.domain_id == c.domain_id && img.segment_id == c.segment_id
    end
    @test _apply_symmetry_to_curve(YAxisReflection(), c).arc_angle == -c.arc_angle
    @test _apply_symmetry_to_curve(NFoldRotation(4, 3), c).arc_angle == c.arc_angle
    @test isapprox(_apply_symmetry_to_curve(NFoldRotation(4, 3), c).center, apply_symmetry(NFoldRotation(4, 3), c.center); atol = fb_ATOL)
end

@testset "_reverse_curve (LineSegment)" begin
    c = fb_line()
    rc = _reverse_curve(c)
    @test rc isa LineSegment
    @test rc.pt0 == c.pt1 && rc.pt1 == c.pt0                  # endpoints swapped
    @test fb_reversal_matches(c, rc)
    @test isapprox(rc.length, c.length; atol = fb_ATOL)
    @test rc.domain_id == c.domain_id && rc.segment_id == c.segment_id
    @test fb_reversal_matches(rc, _reverse_curve(rc))          # reversing twice restores the parametrization
end

@testset "_reverse_curve (CircleSegment)" begin
    c = fb_arc()
    rc = _reverse_curve(c)
    @test rc isa CircleSegment
    @test rc.arc_angle == -c.arc_angle
    @test rc.shift_angle == c.shift_angle + c.arc_angle
    @test rc.radius == c.radius && rc.center == c.center
    @test isapprox(curve(rc, 0.0), curve(c, 1.0); atol = fb_ATOL)   # endpoints swapped
    @test isapprox(curve(rc, 1.0), curve(c, 0.0); atol = fb_ATOL)
    @test fb_reversal_matches(c, rc)
    @test isapprox(rc.length, c.length; atol = fb_ATOL)
    @test rc.domain_id == c.domain_id && rc.segment_id == c.segment_id
    @test fb_reversal_matches(rc, _reverse_curve(rc))
end

@testset "_reverse_curve (FourierCoeffPolarSegment)" begin
    c = fb_fourier()
    rc = _reverse_curve(c)
    @test rc isa FourierCoeffPolarSegment
    @test rc.arc_angle == -c.arc_angle
    @test rc.shift_angle == c.shift_angle + c.arc_angle
    @test rc.coef == c.coef && rc.center == c.center
    @test isapprox(curve(rc, 0.0), curve(c, 1.0); atol = fb_ATOL)
    @test isapprox(curve(rc, 1.0), curve(c, 0.0); atol = fb_ATOL)
    @test fb_reversal_matches(c, rc)
    @test rc.domain_id == c.domain_id && rc.segment_id == c.segment_id
    @test fb_reversal_matches(rc, _reverse_curve(rc))
end

@testset "_reverse_curve (PolarSegment)" begin
    c = fb_polar()
    rc = _reverse_curve(c)
    @test rc isa PolarSegment
    @test rc.arc_angle == -c.arc_angle
    @test rc.shift_angle == c.shift_angle + c.arc_angle
    @test rc.center == c.center
    @test isapprox(curve(rc, 0.0), curve(c, 1.0); atol = fb_ATOL)
    @test isapprox(curve(rc, 1.0), curve(c, 0.0); atol = fb_ATOL)
    @test fb_reversal_matches(c, rc)
    @test rc.domain_id == c.domain_id && rc.segment_id == c.segment_id
    @test fb_reversal_matches(rc, _reverse_curve(rc))
end

@testset "reversed reflection images (as built by full_boundary)" begin
    # the image of a reflection runs against the boundary traversal; reversing it restores the sense and
    # keeps the physical trace, so curve(rev(img), t) == sym(curve(c, 1 - t)) with a positive arc angle
    for c in (fb_line(), fb_arc(), fb_fourier(), fb_polar())
        for sym in fb_REFLECTIONS
            img = _reverse_curve(_apply_symmetry_to_curve(sym, c))
            @test typeof(img).name == typeof(c).name
            @test all(t -> isapprox(curve(img, t), apply_symmetry(sym, curve(c, 1 - t)); atol = fb_ATOL), fb_TS)
            if !(c isa LineSegment)
                @test img.arc_angle == c.arc_angle
            end
        end
    end
end

@testset "full_boundary" begin
    # trivial symmetry: the fundamental curves already are the complete boundary
    tri = TriangleBilliard(1.0, 1.0)
    @test full_boundary(tri) == get_boundary_curves(tri)

    # D2 stadium: four copies of the fundamental curves
    stad = StadiumBilliard(0.5)
    @test length(full_boundary(stad)) == 4 * length(get_boundary_curves(stad))
    # C3 billiard: three copies of the single fundamental arc
    c3 = C3Billiard(0.2)
    @test length(full_boundary(c3)) == 3 * length(get_boundary_curves(c3))
    @test all(c -> c isa FourierCoeffPolarSegment, full_boundary(c3))

    # D2-symmetric polar billiard: the fundamental curve plus one image per non-identity symmetry
    polar = PolarBilliard([0.0, 0.0, 0.0, 0.3])
    polar_d2 = PolarBilliard{Float64}(polar.fundamental_domain, D2_symmetry())
    fb = full_boundary(polar_d2)
    @test length(fb) == 4
    @test all(c -> c isa FourierCoeffPolarSegment, fb)

    # the first sector is the fundamental boundary itself, and full_boundary drops the provenance vectors
    rect = RectangleBilliard(1.0, 0.6)
    @test full_boundary(rect)[1:2] == get_boundary_curves(rect)
    @test full_boundary(rect) == first(_full_boundary_data(rect))

    # G13: every wall of the unfolded rectangle sees an interior point on its inner side
    for p in (SVector(0.1, 0.1), SVector(-0.5, 0.3), SVector(0.7, -0.2))
        @test all(c -> domain_fun(c, p) < 0, full_boundary(rect))
    end
end

@testset "_full_boundary_data" begin
    # D2 rectangle: registry (Y, XY, X) gives the sectors 2, 3, 4 behind the fundamental sector 1
    rect = RectangleBilliard(1.0, 0.6)
    P = get_boundary_curves(rect)
    curves, source, sector, reversed = _full_boundary_data(rect)
    @test length(curves) == length(source) == length(sector) == length(reversed) == 8
    @test reversed isa BitVector
    # reflections (Y, X) list their images in reverse order and reversed, the pi rotation XY keeps the order
    @test source == [1, 2, 2, 1, 1, 2, 2, 1]
    @test sector == [1, 1, 2, 2, 3, 3, 4, 4]
    @test reversed == Bool[0, 0, 1, 1, 0, 0, 1, 1]
    @test curves[1:2] == P

    # every reconstructed curve is the (possibly reversed) image of its source curve under its sector's symmetry
    syms = rect.symmetries
    for i in eachindex(curves)
        sym = i <= 2 ? nothing : syms[sector[i] - 1]
        for t in fb_TS
            tt = reversed[i] ? 1 - t : t
            base = curve(P[source[i]], tt)
            expected = isnothing(sym) ? base : apply_symmetry(sym, base)
            @test isapprox(curve(curves[i], t), expected; atol = fb_ATOL)
        end
    end

    # C3 billiard: rotations keep curve order and parametrization
    c3 = C3Billiard(0.2)
    _, source3, sector3, reversed3 = _full_boundary_data(c3)
    @test source3 == [1, 1, 1]
    @test sector3 == [1, 2, 3]
    @test !any(reversed3)

    # no symmetry: only sector 1, nothing reversed
    tri = TriangleBilliard(1.0, 1.0)
    _, source_t, sector_t, reversed_t = _full_boundary_data(tri)
    @test source_t == [1, 2, 3]
    @test sector_t == [1, 1, 1]
    @test !any(reversed_t)
end
