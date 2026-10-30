# Unit tests for the reconstruction of the complete physical boundary from the fundamental domain and its symmetries:
#   src/geometry/fullboundary.jl   (full_boundary, _full_boundary_data, _apply_symmetry_to_curve, _polar_symmetry_image,
#                                   _reverse_curve, _translate_curve, _translate_domain, _shift_bc, _shift_curve_ids,
#                                   _apply_symmetry_to_domain)
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
# reconstructed loop and its total length on real billiards are covered by the integration tier. Billiards with holes
# (hole fixtures of holes_fixtures.jl): full_boundary_components, BoundaryComponentMap, the MCD method of full_boundary and the
# split-out _reconstruct_curves; the outer loop is reconstructed by the registry, hole loops are taken in the full frame. The `blocks` field of
# BoundaryComponentMap (one vector of dihedral keys per group, identity first) is the registry keys for the outer group, the hole-group keys
# of a SymmetricHole or SymmetrizedHoles and the identity alone for a DomainHole. No eigen solves,
# no pinned wavefunction values. Helper functions and
# constants use the prefix `fb_`.

using StaticArrays: SVector
using BilliardGeometry: _apply_symmetry_to_curve, _reverse_curve, _full_boundary_data, _reconstruct_curves,
                        _block_keys, _IDENTITY_KEY, symmetry_key

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

@testset "_reconstruct_curves" begin
    # the billiard method is _reconstruct_curves of the fundamental curves and the billiard registry
    rect = RectangleBilliard(1.0, 0.6)
    P = get_boundary_curves(rect)
    @test all(_reconstruct_curves(P, rect.symmetries) .== _full_boundary_data(rect))
    # empty registry: the curves themselves, sector 1, nothing reversed
    curves, source, sector, reversed = _reconstruct_curves(P, register_symmetries())
    @test curves == P
    @test source == [1, 2]
    @test sector == [1, 1]
    @test !any(reversed)
    # D2 registry: the same sectors as for the rectangle billiard
    _, source_d2, sector_d2, _ = _reconstruct_curves(P, D2_symmetry())
    @test source_d2 == [1, 2, 2, 1, 1, 2, 2, 1]
    @test sector_d2 == [1, 1, 2, 2, 3, 3, 4, 4]
    # for a billiard with holes only the outer boundary is reconstructed
    d2 = hf_d2_with_hole()
    curves_h, source_h, sector_h, _ = _full_boundary_data(d2)
    @test length(curves_h) == 8
    @test sector_h == [1, 1, 2, 2, 3, 3, 4, 4]
    @test source_h == [1, 2, 2, 1, 1, 2, 2, 1]
end

@testset "BoundaryComponentMap" begin
    @test fieldnames(BoundaryComponentMap) == (:curves, :is_hole, :group, :blocks)
    map = BoundaryComponentMap([AbsCurve[LineSegment([0.0, 0.0], [1.0, 0.0])]], [false], [0], [[(false, 0 // 1)]])
    @test map.curves isa Vector{Vector{AbsCurve}}
    @test map.is_hole == [false]
    @test map.group == [0]
    @test map.blocks isa Vector{Vector{Tuple{Bool,Rational{Int}}}}
    @test map.blocks == [[(false, 0 // 1)]]
    # the constructor takes four arguments, the old three-argument form no longer exists
    @test_throws MethodError BoundaryComponentMap([AbsCurve[LineSegment([0.0, 0.0], [1.0, 0.0])]], [false], [0])
end

# True when the loop `crvs` is closed: every curve ends where the next one starts (cyclically).
fb_is_closed_loop(crvs) = all(isapprox(curve(crvs[j], 1.0), curve(crvs[mod1(j + 1, length(crvs))], 0.0); atol = 1e-10) for j in eachindex(crvs))

@testset "full_boundary_components" begin
    # without holes: one component equal to full_boundary
    rect = RectangleBilliard(1.0, 0.6)
    map = full_boundary_components(rect)
    @test map isa BoundaryComponentMap
    @test length(map.curves) == 1
    @test map.curves[1] == full_boundary(rect)
    @test map.is_hole == [false]
    @test map.group == [0]
    @test full_boundary_components(TriangleBilliard(1.0, 1.0)).group == [0]

    # concentric annulus: outer circle and hole circle
    ann = full_boundary_components(hf_concentric_annulus())
    @test [length(c) for c in ann.curves] == [1, 1]
    @test ann.is_hole == [false, true]
    @test ann.group == [0, 1]
    @test ann.curves[1][1].radius == 1.0
    @test ann.curves[2][1].radius == 0.3
    @test all(fb_is_closed_loop, ann.curves)

    # eccentric annulus
    ecc = full_boundary_components(hf_eccentric_annulus())
    @test [length(c) for c in ecc.curves] == [1, 1]
    @test ecc.is_hole == [false, true]
    @test ecc.group == [0, 1]
    @test ecc.curves[2][1].center == SVector(0.3, 0.0)

    # rectangle with a rectangular hole, no symmetry: two loops of four curves
    rh = full_boundary_components(hf_rect_with_hole())
    @test [length(c) for c in rh.curves] == [4, 4]
    @test rh.is_hole == [false, true]
    @test rh.group == [0, 1]
    @test all(fb_is_closed_loop, rh.curves)
    # outer loop counterclockwise, hole loop counterclockwise (hole on the left): Green areas 6 and 0.5
    @test area(rh.curves[1]) ≈ 6.0 atol = 1e-12
    @test area(rh.curves[2]) ≈ 0.5 atol = 1e-12

    # D2 billiard: the outer loop is reconstructed from the quarter (8 curves), the hole stays one curve in the full frame
    d2 = full_boundary_components(hf_d2_with_hole())
    @test [length(c) for c in d2.curves] == [8, 1]
    @test d2.is_hole == [false, true]
    @test d2.group == [0, 1]
    @test fb_is_closed_loop(d2.curves[1])
    @test area(d2.curves[1]) ≈ 4.0 atol = 1e-12
    @test d2.curves[2][1].center == SVector(0.3, 0.4)
    @test d2.curves[2][1].radius == 0.2
    # the outer loop is the first loop of the reconstruction, not a coordinate-matched reordering
    @test d2.curves[1] == first(_full_boundary_data(hf_d2_with_hole()))

    # several hole objects: groups count the hole objects, the loops of one object share its group
    dom = MultiplyConnectedDomain(hf_rect_domain(0.0, 0.0, 4.0, 4.0, 1), [hf_circle_domain(0.3, 3.0, 3.0, 2), hf_two_loop_hole(5)], 1)
    multi = full_boundary_components(hf_HoleBilliard{Float64}(dom, register_symmetries()))
    @test [length(c) for c in multi.curves] == [4, 1, 4, 4]
    @test multi.is_hole == [false, true, true, true]
    @test multi.group == [0, 1, 2, 2]
end

@testset "full_boundary_components (blocks)" begin
    ID = _IDENTITY_KEY
    # without holes: one group whose blocks are the identity followed by the registry, as _block_keys of the billiard
    for b in (RectangleBilliard(1.0, 0.6), CircleBilliard(1.0), C3Billiard(0.2))
        map = full_boundary_components(b)
        @test map.blocks == [_block_keys(b)]
        @test map.blocks[1][1] == ID
        @test length(map.blocks) == length(unique(map.group))
    end
    @test full_boundary_components(TriangleBilliard(1.0, 1.0)).blocks == [[ID]]

    # DomainHole groups carry the identity block alone, one entry per group (outer first)
    @test full_boundary_components(hf_concentric_annulus()).blocks == [[ID], [ID]]
    @test full_boundary_components(hf_rect_with_hole()).blocks == [[ID], [ID]]
    dom = MultiplyConnectedDomain(hf_rect_domain(0.0, 0.0, 4.0, 4.0, 1), [hf_circle_domain(0.3, 3.0, 3.0, 2), hf_two_loop_hole(5)], 1)
    multi = full_boundary_components(hf_HoleBilliard{Float64}(dom, register_symmetries()))
    @test multi.blocks == [[ID], [ID], [ID]]       # outer, hole 1, hole 2 (two loops share the group of hole 2)
    @test multi.group == [0, 1, 2, 2]
    d2dom = full_boundary_components(hf_d2_with_hole())
    @test d2dom.blocks == [_block_keys(hf_d2_with_hole()), [ID]]
    @test d2dom.blocks == [[ID, symmetry_key(YAxisReflection()), symmetry_key(XYAxisReflection()), symmetry_key(XAxisReflection())], [ID]]

    # symmetric holes: the keys of the hole group, in chain order, identity first
    mk(hole) = hf_HoleBilliard{Float64}(MultiplyConnectedDomain(hf_quarter_square_domain(), [hole], 1), D2_symmetry())
    d4 = full_boundary_components(mk(hf_d4_symmetric_hole()))
    @test d4.group == [0, 1]
    @test length.(d4.blocks) == [4, 8]
    @test d4.blocks[2] == _block_keys(D4_symmetry())
    @test d4.blocks[2][1] == ID
    @test d4.blocks[1] == _block_keys(D2_symmetry())
    half = full_boundary_components(mk(hf_half_disc_hole()))
    @test half.blocks[2] == [ID, symmetry_key(XAxisReflection())]

    # copies of a SymmetrizedHoles: one loop per block, all in group 1
    sym = full_boundary_components(mk(hf_symmetrized_disc()))
    @test sym.group == [0, 1, 1, 1, 1]
    @test length.(sym.blocks) == [4, 4]
    @test sym.blocks[2] == _block_keys(D2_symmetry())
    # the number of loops of a group is the number of its blocks for symmetric holes
    @test count(==(1), sym.group) == length(sym.blocks[2])
    @test count(==(1), d4.group) == 1               # one loop of 8 curves: the blocks are contiguous curve ranges of that loop
    @test length(d4.curves[2]) == length(d4.blocks[2])

    # DomainHole in a D2 billiard: its group has the identity block alone
    dh = full_boundary_components(hf_d2_with_hole())
    @test dh.blocks[2] == [ID]
    @test length(dh.blocks) == 2
end

@testset "full_boundary (billiards with holes)" begin
    # the concatenation of the loops of full_boundary_components, outer loop first
    for b in (hf_concentric_annulus(), hf_eccentric_annulus(), hf_rect_with_hole(), hf_d2_with_hole())
        loops = full_boundary_components(b).curves
        @test full_boundary(b) == reduce(vcat, loops)
        @test length(full_boundary(b)) == sum(length, loops)
    end
    @test length(full_boundary(hf_concentric_annulus())) == 2
    @test length(full_boundary(hf_d2_with_hole())) == 9
    @test full_boundary(hf_d2_with_hole())[end].radius == 0.2
    # the hole curves keep domain_fun pointing into the hole: an interior point of the hole is on the inner side of its boundary
    ann = hf_concentric_annulus()
    @test domain_fun(full_boundary(ann)[2], SVector(0.0, 0.0)) < 0
    @test domain_fun(full_boundary(ann)[1], SVector(0.0, 0.0)) < 0
end

################################################################################
##################### translation and domain images (B3) #######################
################################################################################
# _translate_curve, _translate_domain, _shift_bc, _shift_curve_ids, _apply_symmetry_to_domain (src/geometry/fullboundary.jl).
# Evidence: the translated curve is the old one plus the vector at every parameter, with unchanged length, tangent and metadata; the
# image of a domain under g is g applied pointwise (reflections: curves and corners in reverse order, each curve reversed, so that
# curve(img_k, t) = g(curve(c_{n+1-k}, 1-t))); the unit square keeps its signed Green area 1 (region on the left) under every
# element; a rotation by pi/2 about the center of a unit square maps the k-th side onto the (k+1)-th; an offset shifts domain ids
# and the neighbour ids of Transparent/PeriodicX only.

using BilliardGeometry: _translate_curve, _translate_domain, _shift_bc, _shift_curve_ids, _apply_symmetry_to_domain, _area_integral

const fb_V = SVector(1.0, 2.0)

# Curves with a non-default boundary condition and orientation, off the origin.
fb_line_bc() = LineSegment([0.3, 0.1], [1.2, 0.7]; bc = SymmetryWall(3, 4), orientation = -1, domain_id = 2, segment_id = 5)
fb_arc_bc() = CircleSegment(0.8, 1.1, 0.4, [0.5, 0.2]; bc = SymmetryWall(3, 4), orientation = -1, domain_id = 2, segment_id = 5)
fb_fourier_bc() = FourierCoeffPolarSegment([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3]; arc_angle = 1.3, shift_angle = 0.4, center = [0.3, -0.2],
                                           bc = SymmetryWall(3, 4), orientation = -1, domain_id = 2, segment_id = 5)
fb_polar_bc() = PolarSegment(φ -> 1.0 + 0.2 * cos(4 * φ); arc_angle = 1.3, shift_angle = 0.4, center = [0.3, -0.2],
                             bc = SymmetryWall(3, 4), orientation = -1, domain_id = 2, segment_id = 5)

# Unit square [1,2] x [0.5,1.5] counterclockwise (segment ids 1..4, domain id 3), signed Green area 1.
fb_unit_square() = hf_rect_domain(1.0, 0.5, 2.0, 1.5, 3)

# The same square with a different boundary condition on each side, domain id `id`.
function fb_bc_square(id)
    p = (SVector(1.0, 0.5), SVector(2.0, 0.5), SVector(2.0, 1.5), SVector(1.0, 1.5))
    bcs = AbsBoundaryCondition[SpecularReflection(), Transparent(7), SymmetryWall(3, 4), PeriodicX(9)]
    curves = AbsCurve[LineSegment(p[k], p[mod1(k + 1, 4)]; bc = bcs[k], domain_id = id, segment_id = k) for k in 1:4]
    return SimpleDomain{Float64}(curves, collect(p), id)
end

# Signed Green area of the whole boundary of a simple domain.
fb_signed_area(domain) = sum(_area_integral(c) for c in domain.boundary) / 2

# Pointwise approximate equality of two equally long point lists.
fb_points_match(a, b) = length(a) == length(b) && all(isapprox(x, y; atol = fb_ATOL) for (x, y) in zip(a, b))

# Every curve ends where the next one starts, the last one where the first one starts.
fb_chain_closed(curves) = all(k -> isapprox(curve(curves[k], 1.0), curve(curves[mod1(k + 1, length(curves))], 0.0); atol = fb_ATOL), eachindex(curves))

# Same trace and parametrization as the vector-shifted source at the grid `fb_TS`, and the same metadata.
function fb_translated_matches(c, tc, v)
    pts = all(t -> isapprox(curve(tc, t), curve(c, t) + v; atol = fb_ATOL), fb_TS)
    tans = all(t -> isapprox(tangent(tc, t), tangent(c, t); atol = fb_ATOL), fb_TS)
    meta = tc.bc == c.bc && tc.orientation == c.orientation && tc.domain_id == c.domain_id && tc.segment_id == c.segment_id
    return pts && tans && meta && isapprox(tc.length, c.length; atol = fb_ATOL)
end

@testset "_translate_curve (LineSegment)" begin
    c = fb_line_bc()
    tc = _translate_curve(c, fb_V)
    @test tc isa LineSegment{Float64}
    @test fb_translated_matches(c, tc, fb_V)
    @test tc.pt0 ≈ c.pt0 + fb_V atol = fb_ATOL
    @test tc.pt1 ≈ c.pt1 + fb_V atol = fb_ATOL
    # the source is not modified and the inverse translation gives the source back
    @test c.pt0 == [0.3, 0.1]
    @test fb_translated_matches(c, _translate_curve(tc, -fb_V), SVector(0.0, 0.0))
    @test fb_translated_matches(c, _translate_curve(c, SVector(0.0, 0.0)), SVector(0.0, 0.0))
end

@testset "_translate_curve (CircleSegment)" begin
    c = fb_arc_bc()
    tc = _translate_curve(c, fb_V)
    @test tc isa CircleSegment{Float64}
    @test fb_translated_matches(c, tc, fb_V)
    @test tc.center ≈ c.center + fb_V atol = fb_ATOL
    @test tc.radius == c.radius && tc.arc_angle == c.arc_angle && tc.shift_angle == c.shift_angle
    @test fb_translated_matches(c, _translate_curve(tc, -fb_V), SVector(0.0, 0.0))
end

@testset "_translate_curve (FourierCoeffPolarSegment)" begin
    c = fb_fourier_bc()
    tc = _translate_curve(c, fb_V)
    @test tc isa FourierCoeffPolarSegment
    @test fb_translated_matches(c, tc, fb_V)
    @test tc.center ≈ c.center + fb_V atol = fb_ATOL
    # only the center moves: coefficients, scale and angles are kept
    @test tc.coef == c.coef && tc.R == c.R && tc.arc_angle == c.arc_angle && tc.shift_angle == c.shift_angle
    @test fb_translated_matches(c, _translate_curve(tc, -fb_V), SVector(0.0, 0.0))
end

@testset "_translate_curve (PolarSegment)" begin
    c = fb_polar_bc()
    tc = _translate_curve(c, fb_V)
    @test tc isa PolarSegment
    @test fb_translated_matches(c, tc, fb_V)
    @test tc.center ≈ c.center + fb_V atol = fb_ATOL
    @test tc.r_func === c.r_func && tc.R == c.R && tc.arc_angle == c.arc_angle && tc.shift_angle == c.shift_angle
    @test fb_translated_matches(c, _translate_curve(tc, -fb_V), SVector(0.0, 0.0))
end

@testset "_translate_domain" begin
    dom = fb_unit_square()
    td = _translate_domain(dom, fb_V)
    @test td isa SimpleDomain{Float64}
    @test td.id == dom.id
    @test length(td.boundary) == 4
    @test [c.segment_id for c in td.boundary] == [1, 2, 3, 4]
    @test all(k -> fb_translated_matches(dom.boundary[k], td.boundary[k], fb_V), 1:4)
    @test fb_points_match(td.corners, [p + fb_V for p in dom.corners])
    @test fb_chain_closed(td.boundary)
    @test fb_signed_area(td) ≈ 1.0 atol = fb_ATOL
    # the source is untouched and the inverse translation returns the corners
    @test dom.corners[1] == SVector(1.0, 0.5)
    @test fb_points_match(_translate_domain(td, -fb_V).corners, dom.corners)
    # any simple domain type gives a SimpleDomain, here a CircleWedge with an arc
    wedge = hf_d2_quarter_disc()
    tw = _translate_domain(wedge, fb_V)
    @test tw isa SimpleDomain{Float64}
    @test tw.id == wedge.id
    @test tw.boundary[1].center ≈ wedge.center + fb_V atol = fb_ATOL
    @test [c.bc for c in tw.boundary] == [c.bc for c in wedge.boundary]
    @test fb_chain_closed(tw.boundary)
end

@testset "_shift_bc" begin
    @test _shift_bc(Transparent(7), 10) == Transparent(17)
    @test _shift_bc(PeriodicX(9), 10) == PeriodicX(19)
    @test _shift_bc(Transparent(7), 0) == Transparent(7)
    @test _shift_bc(Transparent(7), -3) == Transparent(4)
    # every other condition has no neighbour id and is returned unchanged
    @test _shift_bc(SpecularReflection(), 10) == SpecularReflection()
    @test _shift_bc(SymmetryWall(3, 4), 10) == SymmetryWall(3, 4)
    @test _shift_bc(QuantumSolverIgnore(), 10) == QuantumSolverIgnore()
end

@testset "_shift_curve_ids" begin
    for c in (fb_line_bc(), fb_arc_bc(), fb_fourier_bc(), fb_polar_bc())
        shifted = _shift_curve_ids(c, 10)
        @test typeof(shifted) == typeof(c)
        @test shifted.domain_id == c.domain_id + 10
        @test shifted.segment_id == c.segment_id
        @test shifted.bc == SymmetryWall(3, 4)          # a SymmetryWall keeps sym_id and sector_id
        @test c.domain_id == 2                           # the source is untouched
        # zero offset returns the curve itself
        @test _shift_curve_ids(c, 0) === c
    end
    # the neighbour ids of Transparent and PeriodicX shift together with the domain id
    tr = LineSegment([0.0, 0.0], [1.0, 0.0]; bc = Transparent(7), domain_id = 2, segment_id = 1)
    @test _shift_curve_ids(tr, 10).bc == Transparent(17)
    @test _shift_curve_ids(tr, 10).domain_id == 12
    pr = LineSegment([0.0, 0.0], [1.0, 0.0]; bc = PeriodicX(9), domain_id = 2, segment_id = 1)
    @test _shift_curve_ids(pr, 10).bc == PeriodicX(19)
    @test _shift_curve_ids(CircleSegment(0.5, pi, 0.0, [0.0, 0.0]; bc = Transparent(7), domain_id = 2), 10).bc == Transparent(17)
end

@testset "_apply_symmetry_to_domain (reflections)" begin
    dom = fb_unit_square()
    inside = SVector(1.5, 1.0)
    outside_pts = (SVector(0.5, 1.0), SVector(2.5, 1.0), SVector(1.5, 0.2), SVector(1.5, 1.8))
    for sym in fb_REFLECTIONS
        img = _apply_symmetry_to_domain(sym, dom)
        @test img isa SimpleDomain{Float64}
        @test img.id == dom.id
        # the curve list is reversed, every curve reversed and mapped: curve(img_k, t) = g(curve(c_{n+1-k}, 1-t))
        @test [c.segment_id for c in img.boundary] == [4, 3, 2, 1]
        @test all(c -> c.domain_id == 3, img.boundary)
        for k in 1:4, t in fb_TS
            @test isapprox(curve(img.boundary[k], t), apply_symmetry(sym, curve(dom.boundary[5 - k], 1 - t)); atol = fb_ATOL)
        end
        @test fb_points_match(img.corners, [SVector{2}(apply_symmetry(sym, p)) for p in Iterators.reverse(dom.corners)])
        # one closed loop with the region on the left, the same signed area as the source
        @test fb_chain_closed(img.boundary)
        @test fb_signed_area(img) ≈ 1.0 atol = fb_ATOL
        # region: p is inside the source iff g(p) is inside the image
        @test all(is_inside(img, SVector{2}(apply_symmetry(sym, inside))))
        @test all(p -> !all(is_inside(img, SVector{2}(apply_symmetry(sym, p)))), outside_pts)
        # the domain function keeps its sign (negative inside) on the reversed curves
        for k in 1:4
            @test domain_fun(img.boundary[k], SVector{2}(apply_symmetry(sym, inside))) ≈ domain_fun(dom.boundary[5 - k], inside) atol = fb_ATOL
        end
    end
    # the image under X is the mirror image about the x axis, the square [1,2] x [-1.5,-0.5], corners in reverse order
    img = _apply_symmetry_to_domain(XAxisReflection(), dom)
    @test fb_points_match(img.corners, [SVector(1.0, -1.5), SVector(2.0, -1.5), SVector(2.0, -0.5), SVector(1.0, -0.5)])
end

@testset "_apply_symmetry_to_domain (rotations and the pi rotation)" begin
    dom = fb_unit_square()
    inside = SVector(1.5, 1.0)
    outside_pts = (SVector(0.5, 1.0), SVector(2.5, 1.0), SVector(1.5, 0.2), SVector(1.5, 1.8))
    for sym in (NFoldRotation(4, 1), NFoldRotation(3, 1), NFoldRotation(4, 3), XYAxisReflection())
        img = _apply_symmetry_to_domain(sym, dom)
        @test img isa SimpleDomain{Float64}
        @test img.id == dom.id
        # orientation-preserving: the order of the curves and corners is kept
        @test [c.segment_id for c in img.boundary] == [1, 2, 3, 4]
        @test all(c -> c.domain_id == 3, img.boundary)
        for k in 1:4, t in fb_TS
            @test isapprox(curve(img.boundary[k], t), apply_symmetry(sym, curve(dom.boundary[k], t)); atol = fb_ATOL)
        end
        @test fb_points_match(img.corners, [SVector{2}(apply_symmetry(sym, p)) for p in dom.corners])
        @test fb_chain_closed(img.boundary)
        @test fb_signed_area(img) ≈ 1.0 atol = fb_ATOL
        @test all(is_inside(img, SVector{2}(apply_symmetry(sym, inside))))
        @test all(p -> !all(is_inside(img, SVector{2}(apply_symmetry(sym, p)))), outside_pts)
    end
    # the pi rotation about the origin: the square [-2,-1] x [-1.5,-0.5]
    img = _apply_symmetry_to_domain(XYAxisReflection(), dom)
    @test fb_points_match(img.corners, [SVector(-1.0, -0.5), SVector(-2.0, -0.5), SVector(-2.0, -1.5), SVector(-1.0, -1.5)])
    # a wedge with an arc and symmetry walls: rotation by pi/2 maps the wedge [0, pi/4] onto [pi/2, 3pi/4]
    wedge = hf_d4_eighth_disc()
    img = _apply_symmetry_to_domain(NFoldRotation(4, 1), wedge)
    @test fb_chain_closed(img.boundary)
    @test all(is_inside(img, SVector(-0.05, 0.1)))
    @test !all(is_inside(img, SVector(0.1, 0.05)))
end

@testset "_apply_symmetry_to_domain (id_offset)" begin
    dom = fb_bc_square(3)
    # a rotation keeps the order: ids and neighbour ids shift, segment ids and SymmetryWall conditions stay
    img = _apply_symmetry_to_domain(NFoldRotation(4, 1), dom, 10)
    @test img.id == 13
    @test all(c -> c.domain_id == 13, img.boundary)
    @test [c.segment_id for c in img.boundary] == [1, 2, 3, 4]
    @test [c.bc for c in img.boundary] == AbsBoundaryCondition[SpecularReflection(), Transparent(17), SymmetryWall(3, 4), PeriodicX(19)]
    # a reflection reverses the list, the conditions travel with their curves
    img = _apply_symmetry_to_domain(YAxisReflection(), dom, 10)
    @test img.id == 13
    @test all(c -> c.domain_id == 13, img.boundary)
    @test [c.segment_id for c in img.boundary] == [4, 3, 2, 1]
    @test [c.bc for c in img.boundary] == AbsBoundaryCondition[PeriodicX(19), SymmetryWall(3, 4), Transparent(17), SpecularReflection()]
    # offset zero (also the default) changes no id
    for img0 in (_apply_symmetry_to_domain(NFoldRotation(4, 1), dom, 0), _apply_symmetry_to_domain(NFoldRotation(4, 1), dom))
        @test img0.id == 3
        @test all(c -> c.domain_id == 3, img0.boundary)
        @test [c.bc for c in img0.boundary] == [c.bc for c in dom.boundary]
    end
    # the source is untouched
    @test all(c -> c.domain_id == 3, dom.boundary)
    @test dom.boundary[2].bc == Transparent(7)
    # the wedge piece of a symmetric hole keeps the sym_id and sector_id of its walls
    wedge = hf_d4_eighth_disc()
    img = _apply_symmetry_to_domain(NFoldRotation(4, 1), wedge, 10)
    @test img.id == 12
    @test all(c -> c.domain_id == 12, img.boundary)
    @test [c.bc for c in img.boundary] == [c.bc for c in wedge.boundary]
    @test [c.segment_id for c in img.boundary] == [1, 2, 3]
    @test all(is_inside(img, SVector(-0.05, 0.1)))
end

@testset "_apply_symmetry_to_domain (with center)" begin
    dom = fb_unit_square()
    center = SVector(0.3, 0.1)
    # reflection about the horizontal line y = 0.1 through `center`: p -> center + g(p - center)
    sym = XAxisReflection()
    img = _apply_symmetry_to_domain(sym, dom, center, 0)
    @test img isa SimpleDomain{Float64}
    @test [c.segment_id for c in img.boundary] == [4, 3, 2, 1]
    for k in 1:4, t in fb_TS
        @test isapprox(curve(img.boundary[k], t), center + apply_symmetry(sym, curve(dom.boundary[5 - k], 1 - t) - center); atol = fb_ATOL)
    end
    @test fb_points_match(img.corners, [SVector(1.0, -1.3), SVector(2.0, -1.3), SVector(2.0, -0.3), SVector(1.0, -0.3)])
    @test fb_chain_closed(img.boundary)
    @test fb_signed_area(img) ≈ 1.0 atol = fb_ATOL
    # applying the reflection twice about the same center gives the source back
    back = _apply_symmetry_to_domain(sym, img, center, 0)
    @test [c.segment_id for c in back.boundary] == [1, 2, 3, 4]
    @test fb_points_match(back.corners, dom.corners)
    @test all(k -> all(t -> isapprox(curve(back.boundary[k], t), curve(dom.boundary[k], t); atol = fb_ATOL), fb_TS), 1:4)
    # ids shift with the offset and the center is not an id
    shifted = _apply_symmetry_to_domain(sym, dom, center, 100)
    @test shifted.id == 103
    @test all(c -> c.domain_id == 103, shifted.boundary)

    # a rotation by pi/2 about the center of the unit square maps the k-th side onto the (k+1)-th
    sq_center = SVector(1.5, 1.0)
    rot = _apply_symmetry_to_domain(NFoldRotation(4, 1), dom, sq_center, 0)
    for k in 1:4, t in fb_TS
        @test isapprox(curve(rot.boundary[k], t), curve(dom.boundary[mod1(k + 1, 4)], t); atol = fb_ATOL)
    end
    @test fb_points_match(rot.corners, [dom.corners[mod1(k + 1, 4)] for k in 1:4])

    # the zero center is the same as the method about the origin
    zero_center = _apply_symmetry_to_domain(sym, dom, SVector(0.0, 0.0), 5)
    origin = _apply_symmetry_to_domain(sym, dom, 5)
    @test zero_center.id == origin.id
    @test fb_points_match(zero_center.corners, origin.corners)
    @test [c.domain_id for c in zero_center.boundary] == [c.domain_id for c in origin.boundary]

    # an arc about its own center: the Y reflection through the center (0.5, 0) of the quarter disc maps [0, pi/2] onto [pi/2, pi]
    quarter = hf_d2_quarter_disc()
    img = _apply_symmetry_to_domain(YAxisReflection(), quarter, SVector(0.5, 0.0), 0)
    @test all(is_inside(img, SVector(0.45, 0.05)))
    @test !all(is_inside(img, SVector(0.55, 0.05)))
    @test fb_chain_closed(img.boundary)
    @test curve(img.boundary[3], 0.5) ≈ [0.5 - 0.1 * cos(pi / 4), 0.1 * sin(pi / 4)] atol = fb_ATOL
end
