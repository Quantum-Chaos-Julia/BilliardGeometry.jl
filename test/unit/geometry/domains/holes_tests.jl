# Unit tests for src/geometry/domains/holes.jl: the DomainHole type and its interface
#   (DomainHole, hole_domains, hole_leaves, hole_loop_count, get_all_domains, get_all_curves, get_boundary_curves,
#    boundary_components, is_inside, area, symmetry_compatible, _set_curve_bc).
# Reference objects (fixtures from holes_fixtures.jl): the circular hole of radius 0.2 about (0.3, 0.4), the rectangular hole
# [1, 1.5] x [0.5, 1.5], and a composite hole made of two disjoint rectangles (two loops). Areas, membership and loop counts are
# closed forms. No eigen solves. Helpers use the prefix `hl_`.
#
# Symmetric holes (B3 of the symmetric-holes plan): SymmetricHole, SymmetrizedHoles, their interface (hole_domains, hole_leaves,
#   hole_loop_count, get_all_domains, get_all_curves, get_boundary_curves, boundary_components, is_inside, area, symmetry_compatible,
#   _set_curve_bc), the helpers _hole_parts, _part_contains, _apply_symmetry_to_region, and the copy-overlap check that validate_holes
#   (multiplyconnecteddomains.jl) performs for SymmetrizedHoles.
# Reference objects (fixtures from holes_fixtures.jl): the disc of radius 0.3 about the origin as an eighth-disc piece with D4_symmetry()
#   (area pi 0.09, one loop of 8 arcs whose k-th arc starts at the angle (k-1) pi/4), the disc of radius 0.1 about (0.5, 0) as a quarter-disc
#   piece with D2_symmetry() acting about (0.5, 0), the half disc of radius 0.2 about (0.3, 0) with a single mirror, and four discs of radius
#   0.1 placed by D2_symmetry(). Evidence: closed-form areas, loop lengths 2 pi R, membership by radius, chain order of the arcs, positive
#   Green signed area of every loop (hole region on the left).
#
# Compatibility with a registry (B4): compatible_symmetries and the generic symmetry_compatible per hole type, with the helpers
#   _fixes_point, _is_special_axis and _hole_block_keys. A hole is mapped onto itself by g when the dihedral key of g is a key of the
#   hole group and g fixes the hole center (SymmetricHole), or when the key of g permutes the placing keys (SymmetrizedHoles of a
#   DomainHole); a DomainHole and a SymmetrizedHoles of a SymmetricHole have no symmetry of their own. The expected subgroups follow from
#   the dihedral algebra (D2 = {Y, XY, X} has sym_ids 1, 2, 3 and D4_symmetry() has sym_ids 1..7 with Y = 3, XY = 4, X = 7).

using StaticArrays: SVector
using BilliardGeometry: _set_curve_bc, _hole_parts, _part_contains, _apply_symmetry_to_region, _area_integral,
                        _fixes_point, _is_special_axis, _hole_block_keys, _block_keys, _IDENTITY_KEY, symmetry_key

# Absolute tolerance for closed-form areas of the holes (Green quadrature, worst observed 4e-16, margin 1000x).
const hl_ATOL = 1e-12

hl_circle_hole() = DomainHole(hf_circle_domain(0.2, 0.3, 0.4, 5))
hl_rect_hole() = DomainHole(hf_rect_domain(1.0, 0.5, 1.5, 1.5, 2))
# Two disjoint rectangles [0,1]x[0,1] (domain id 3) and [2,3]x[0,1] (domain id 4) as one composite hole object.
hl_two_loop_hole() = hf_two_loop_hole(3)

# Green signed area of a loop of curves (positive when the region is on the left of the traversal).
hl_signed_area(loop) = sum(_area_integral(c) for c in loop) / 2
# Every curve of the loop ends where the next one starts, the last one where the first one starts.
hl_loop_closed(loop) = all(j -> isapprox(curve(loop[j], 1.0), curve(loop[mod1(j + 1, length(loop))], 0.0); atol = hl_ATOL), eachindex(loop))
hl_loop_length(loop) = sum(c.length for c in loop)
hl_polar_pt(center, R, θ) = SVector(center[1] + R * cos(θ), center[2] + R * sin(θ))

@testset "DomainHole constructor" begin
    h = hl_circle_hole()
    @test h isa DomainHole
    @test h isa AbsDomainHole
    @test !(h isa AbsDomain)
    @test h.domain isa SimpleDomain{Float64}
    @test fieldnames(DomainHole) == (:domain,)
    @test hl_two_loop_hole().domain isa CompositeDomain
    @test DomainHole(hl_two_loop_hole().domain) isa DomainHole{CompositeDomain}
end

@testset "hole_domains" begin
    h = hl_circle_hole()
    @test hole_domains(h) isa Vector{AbsDomain}
    @test length(hole_domains(h)) == 1
    @test hole_domains(h)[1] === h.domain
    comp = hl_two_loop_hole()
    @test hole_domains(comp)[1] === comp.domain
end

@testset "hole_leaves" begin
    h = hl_rect_hole()
    @test hole_leaves(h) isa Vector{AbsDomainHole}
    @test length(hole_leaves(h)) == 1
    @test hole_leaves(h)[1] === h
end

@testset "hole_loop_count" begin
    @test hole_loop_count(hl_circle_hole()) == 1
    @test hole_loop_count(hl_rect_hole()) == 1
    # a composite hole of two disjoint rectangles has two closed loops
    @test hole_loop_count(hl_two_loop_hole()) == 2
end

@testset "get_all_domains (DomainHole)" begin
    h = hl_circle_hole()
    @test get_all_domains(h) == AbsSimpleDomain[h.domain]
    comp = hl_two_loop_hole()
    @test [d.id for d in get_all_domains(comp)] == [3, 4]
end

@testset "get_all_curves (DomainHole)" begin
    @test length(get_all_curves(hl_circle_hole())) == 1
    @test get_all_curves(hl_circle_hole())[1] isa CircleSegment
    rect = hl_rect_hole()
    @test length(get_all_curves(rect)) == 4
    @test [c.segment_id for c in get_all_curves(rect)] == [1, 2, 3, 4]
    @test [c.domain_id for c in get_all_curves(hl_two_loop_hole())] == [3, 3, 3, 3, 4, 4, 4, 4]
end

@testset "get_boundary_curves (DomainHole)" begin
    h = hl_circle_hole()
    bc = get_boundary_curves(h)
    @test length(bc) == 1
    @test bc[1].radius == 0.2
    @test bc[1].center == SVector(0.3, 0.4)
    @test length(get_boundary_curves(hl_rect_hole())) == 4
    @test length(get_boundary_curves(hl_two_loop_hole())) == 8
    # a symmetry wall is not physical and is dropped from the physical boundary
    walled = DomainHole(RectangleBilliard(1.0, 1.0).fundamental_domain)
    @test length(get_boundary_curves(walled)) == 2
end

@testset "boundary_components (DomainHole)" begin
    # one closed loop per component, the hole region on the left (counterclockwise for a convex hole)
    comps = boundary_components(hl_circle_hole())
    @test length(comps) == 1
    @test length(comps[1]) == 1

    rect_comps = boundary_components(hl_rect_hole())
    @test length(rect_comps) == 1
    @test length(rect_comps[1]) == 4
    loop = rect_comps[1]
    # closed: every curve ends where the next starts, the last ends where the first starts
    for j in 1:4
        @test curve(loop[j], 1.0) ≈ curve(loop[mod1(j + 1, 4)], 0.0) atol = hl_ATOL
    end
    # counterclockwise: the first edge runs along +x at the bottom
    @test curve(loop[1], 0.0) ≈ [1.0, 0.5] atol = hl_ATOL
    @test curve(loop[1], 1.0) ≈ [1.5, 0.5] atol = hl_ATOL

    two = boundary_components(hl_two_loop_hole())
    @test length(two) == 2
    @test [length(c) for c in two] == [4, 4]
    @test [c[1].domain_id for c in two] == [3, 4]
end

@testset "is_inside (DomainHole)" begin
    circ = hl_circle_hole()
    @test is_inside(circ, SVector(0.3, 0.4)) === true
    @test is_inside(circ, SVector(0.45, 0.4)) === true
    @test is_inside(circ, SVector(0.6, 0.4)) === false
    @test is_inside(circ, SVector(0.0, 0.0)) === false
    # the hole is asymmetric: the mirror image of an inside point is outside
    @test is_inside(circ, SVector(-0.3, 0.4)) === false
    @test is_inside(circ, [SVector(0.3, 0.4), SVector(-0.3, 0.4)]) == [true, false]

    rect = hl_rect_hole()
    @test is_inside(rect, SVector(1.25, 1.0)) === true
    @test is_inside(rect, SVector(1.6, 1.0)) === false
    @test is_inside(rect, SVector(1.25, 0.4)) === false

    # a composite hole is the union of its pieces
    two = hl_two_loop_hole()
    @test is_inside(two, SVector(0.5, 0.5)) === true
    @test is_inside(two, SVector(2.5, 0.5)) === true
    @test is_inside(two, SVector(1.5, 0.5)) === false
    @test is_inside(two, [SVector(0.5, 0.5), SVector(1.5, 0.5), SVector(2.5, 0.5)]) == [true, false, true]
end

@testset "area (DomainHole)" begin
    @test area(hl_circle_hole()) ≈ pi * 0.2^2 atol = hl_ATOL
    @test area(hl_rect_hole()) ≈ 0.5 atol = hl_ATOL
    # sum over the loops
    @test area(hl_two_loop_hole()) ≈ 2.0 atol = hl_ATOL
    # the hole area is positive and independent of the curve orientation field (never read)
    inverted = DomainHole(SimpleDomain{Float64}(AbsCurve[CircleSegment(0.2, 2pi, 0.0, [0.3, 0.4]; orientation = -1, domain_id = 5)], [SVector(0.5, 0.4)], 5))
    @test area(inverted) ≈ pi * 0.2^2 atol = hl_ATOL
    # the quadrature tolerance is forwarded
    @test area(hl_circle_hole(); rtol = 1e-6) ≈ pi * 0.2^2 atol = 1e-8
end

@testset "symmetry_compatible (DomainHole)" begin
    h = hl_circle_hole()
    # a DomainHole has no symmetry: compatible only with the empty registry
    @test symmetry_compatible(h, register_symmetries()) === true
    @test symmetry_compatible(h, D2_symmetry()) === false
    @test symmetry_compatible(h, register_symmetries(XAxisReflection())) === false
    @test symmetry_compatible(hl_rect_hole(), register_symmetries()) === true
    @test symmetry_compatible(hl_two_loop_hole(), register_symmetries(XYAxisReflection())) === false
end

@testset "_set_curve_bc (DomainHole)" begin
    h = hl_rect_hole()
    new_h = _set_curve_bc(h, 2, 3, QuantumSolverIgnore())
    @test new_h isa DomainHole
    @test get_all_curves(new_h)[3].bc isa QuantumSolverIgnore
    # the other curves keep their condition, and the original hole is unchanged
    @test all(c -> c.bc isa SpecularReflection, get_all_curves(new_h)[[1, 2, 4]])
    @test all(c -> c.bc isa SpecularReflection, get_all_curves(h))
    # composite holes are searched by (domain_id, segment_id)
    two = hl_two_loop_hole()
    new_two = _set_curve_bc(two, 4, 2, Transparent(1))
    @test new_two isa DomainHole
    @test get_all_curves(new_two)[6].bc == Transparent(1)
    @test count(c -> !(c.bc isa SpecularReflection), get_all_curves(new_two)) == 1
    # no curve with that id pair: nothing
    @test _set_curve_bc(h, 9, 1, QuantumSolverIgnore()) === nothing
    @test _set_curve_bc(h, 2, 9, QuantumSolverIgnore()) === nothing
end

################################################################################
############################## SymmetricHole ###################################
################################################################################

@testset "SymmetricHole constructor (D4 eighth disc)" begin
    h = hf_d4_symmetric_hole()
    @test h isa SymmetricHole{Float64}
    @test h isa AbsDomainHole
    @test !(h isa AbsDomain)
    @test fieldnames(SymmetricHole) == (:piece, :symmetries, :center, :pieces, :id_stride)
    @test h.center == SVector(0.0, 0.0)
    @test h.id_stride == 100
    @test length(h.symmetries) == 7
    # the identity piece first, then the images in registry order, shifted by k * id_stride
    @test h.pieces isa Vector{AbsSimpleDomain}
    @test length(h.pieces) == 8
    @test all(p -> p isa AbsSimpleDomain, h.pieces)
    @test all(p -> p isa SimpleDomain{Float64}, h.pieces[2:end])     # the images; the first piece is the given piece
    @test [p.id for p in h.pieces] == [2, 102, 202, 302, 402, 502, 602, 702]
    for p in h.pieces
        @test all(c -> c.domain_id == p.id, p.boundary)
        @test length(p.boundary) == 3
        @test count(c -> c.bc isa SpecularReflection, p.boundary) == 1
        # the walls keep their sym_id (1 = Diagonal, 7 = X of the piece) whatever the image
        @test all(c -> c.bc isa SpecularReflection || (c.bc isa SymmetryWall && c.bc.sym_id in (1, 7)), p.boundary)
    end
    # keyword id_stride and a vector center
    h10 = SymmetricHole(hf_d4_eighth_disc(), D4_symmetry(), [0.0, 0.0]; id_stride = 10)
    @test h10.id_stride == 10
    @test [p.id for p in h10.pieces] == [2, 12, 22, 32, 42, 52, 62, 72]
    @test h10.center == SVector(0.0, 0.0)
    # a trivial group (empty registry) keeps the piece
    disc = SymmetricHole(hf_circle_domain(0.2, 0.3, 0.4, 5), register_symmetries())
    @test length(disc.pieces) == 1
    @test area(disc) ≈ pi * 0.04 atol = hl_ATOL
end

@testset "SymmetricHole constructor (off-center centers)" begin
    d2 = hf_d2_symmetric_hole()
    @test d2.center == SVector(0.5, 0.0)
    @test length(d2.pieces) == 4
    @test [p.id for p in d2.pieces] == [2, 102, 202, 302]
    half = hf_half_disc_hole()
    @test half.center == SVector(0.3, 0.0)
    @test length(half.pieces) == 2
    @test [p.id for p in half.pieces] == [3, 103]
end

@testset "SymmetricHole constructor errors" begin
    wedge = hf_d4_eighth_disc()
    # a wall refers to a sym_id that is not in the hole group (sym_id 7 does not exist in a one-element registry)
    @test_throws ArgumentError SymmetricHole(wedge, register_symmetries(DiagonalReflection()))
    @test_throws "sym_id=7" SymmetricHole(wedge, register_symmetries(DiagonalReflection()))
    # registry not in chain order: the images do not follow each other around the loop
    scrambled = register_symmetries(XAxisReflection(), DiagonalReflection(), YAxisReflection(), XYAxisReflection(),
                                    AntiDiagonalReflection(), NFoldRotation(4, 1), NFoldRotation(4, 3))
    @test_throws ArgumentError SymmetricHole(wedge, scrambled)
    @test_throws "do not chain into one closed loop in registry order" SymmetricHole(wedge, scrambled)
    # a registry that is not closed under composition is not a group (R, R^2, R^3 needed)
    @test_throws ArgumentError SymmetricHole(wedge, register_symmetries(NFoldRotation(4, 1)))
    @test_throws "not closed under composition" SymmetricHole(wedge, register_symmetries(NFoldRotation(4, 1)))
    # duplicate geometric actions
    @test_throws ArgumentError SymmetricHole(wedge, register_symmetries(XAxisReflection(), XAxisReflection()))
    # empty piece
    @test_throws ArgumentError SymmetricHole(CompositeDomain(AbsSimpleDomain[]), register_symmetries())
end

@testset "hole interface (SymmetricHole)" begin
    h = hf_d4_symmetric_hole()
    @test hole_loop_count(h) == 1
    @test hole_leaves(h) isa Vector{AbsDomainHole}
    @test length(hole_leaves(h)) == 1 && hole_leaves(h)[1] === h
    @test hole_domains(h) isa Vector{AbsDomain}
    @test length(hole_domains(h)) == 8
    @test [d.id for d in hole_domains(h)] == [2, 102, 202, 302, 402, 502, 602, 702]
    @test get_all_domains(h) == h.pieces
    @test get_all_domains(h) !== h.pieces     # a copy: editing the returned vector does not touch the hole
    @test length(get_all_curves(h)) == 24
    @test sort(unique(c.domain_id for c in get_all_curves(h))) == [2, 102, 202, 302, 402, 502, 602, 702]
    arcs = get_boundary_curves(h)
    @test length(arcs) == 8
    @test all(c -> c isa CircleSegment && c.bc isa SpecularReflection, arcs)
    @test all(c -> c.radius ≈ 0.3, arcs)
    comps = boundary_components(h)
    @test length(comps) == 1
    loop = comps[1]
    @test length(loop) == 8
    # chain order: the k-th arc starts at the angle (k-1) pi/4 about the center, one closed loop of length 2 pi R
    for k in 1:8
        @test curve(loop[k], 0.0) ≈ hl_polar_pt((0.0, 0.0), 0.3, (k - 1) * pi / 4) atol = hl_ATOL
        @test curve(loop[k], 1.0) ≈ hl_polar_pt((0.0, 0.0), 0.3, k * pi / 4) atol = hl_ATOL
    end
    @test hl_loop_closed(loop)
    @test hl_loop_length(loop) ≈ 2 * pi * 0.3 atol = hl_ATOL
    # the region is on the left: positive Green signed area pi R^2
    @test hl_signed_area(loop) ≈ pi * 0.09 atol = hl_ATOL
end

@testset "hole interface (SymmetricHole, off-center groups)" begin
    d2 = hf_d2_symmetric_hole()
    @test hole_loop_count(d2) == 1
    @test length(get_boundary_curves(d2)) == 4
    loop = boundary_components(d2)[1]
    for k in 1:4
        @test curve(loop[k], 0.0) ≈ hl_polar_pt((0.5, 0.0), 0.1, (k - 1) * pi / 2) atol = hl_ATOL
    end
    @test hl_loop_closed(loop)
    @test hl_loop_length(loop) ≈ 2 * pi * 0.1 atol = hl_ATOL
    @test hl_signed_area(loop) ≈ pi * 0.01 atol = hl_ATOL

    half = hf_half_disc_hole()
    @test length(get_all_curves(half)) == 6
    loop = boundary_components(half)[1]
    @test length(loop) == 2
    @test curve(loop[1], 0.0) ≈ hl_polar_pt((0.3, 0.0), 0.2, 0.0) atol = hl_ATOL
    @test curve(loop[2], 0.0) ≈ hl_polar_pt((0.3, 0.0), 0.2, pi) atol = hl_ATOL
    @test hl_loop_closed(loop)
    @test hl_loop_length(loop) ≈ 2 * pi * 0.2 atol = hl_ATOL
    @test hl_signed_area(loop) ≈ pi * 0.04 atol = hl_ATOL
end

@testset "is_inside (SymmetricHole)" begin
    h = hf_d4_symmetric_hole()
    # 16 directions away from the mirror lines: inside at radius 0.2, outside at radius 0.35
    angles = [(j + 0.5) * pi / 8 for j in 0:15]
    @test all(θ -> is_inside(h, hl_polar_pt((0.0, 0.0), 0.2, θ)) === true, angles)
    @test all(θ -> is_inside(h, hl_polar_pt((0.0, 0.0), 0.35, θ)) === false, angles)
    @test is_inside(h, SVector(0.5, 0.0)) === false
    # points on a mirror line (a wall between pieces) are boundary points of the pieces and are avoided here
    @test is_inside(h, [SVector(0.1, -0.05), SVector(0.5, 0.0)]) == [true, false]

    d2 = hf_d2_symmetric_hole()
    @test is_inside(d2, SVector(0.52, 0.05)) === true
    @test is_inside(d2, SVector(0.45, -0.05)) === true
    @test is_inside(d2, SVector(0.5, 0.15)) === false
    @test is_inside(d2, SVector(0.0, 0.0)) === false

    half = hf_half_disc_hole()
    @test is_inside(half, SVector(0.3, -0.1)) === true
    @test is_inside(half, SVector(0.3, 0.1)) === true
    @test is_inside(half, SVector(0.3, 0.25)) === false
    @test is_inside(half, SVector(0.6, 0.0)) === false
end

@testset "area (SymmetricHole)" begin
    @test area(hf_d4_symmetric_hole()) ≈ pi * 0.09 atol = hl_ATOL
    @test area(hf_d2_symmetric_hole()) ≈ pi * 0.01 atol = hl_ATOL
    @test area(hf_half_disc_hole()) ≈ pi * 0.04 atol = hl_ATOL
    @test area(hf_d4_symmetric_hole(); rtol = 1e-6) ≈ pi * 0.09 atol = 1e-8
end

@testset "symmetry_compatible (SymmetricHole, SymmetrizedHoles)" begin
    # the empty registry is compatible with every hole
    for h in (hf_d4_symmetric_hole(), hf_d2_symmetric_hole(), hf_half_disc_hole(), hf_symmetrized_disc())
        @test symmetry_compatible(h, register_symmetries()) === true
    end
    # the D4 disc at the origin is invariant under every element of D2 and D4
    d4 = hf_d4_symmetric_hole()
    @test symmetry_compatible(d4, D2_symmetry()) === true
    @test symmetry_compatible(d4, D4_symmetry()) === true
    # the D2 disc about (0.5, 0): only the x axis mirror fixes the center, so the full D2 and D4 registries are not compatible
    d2 = hf_d2_symmetric_hole()
    @test symmetry_compatible(d2, D2_symmetry()) === false
    @test symmetry_compatible(d2, D4_symmetry()) === false
    @test symmetry_compatible(d2, register_symmetries(XAxisReflection())) === true
    @test symmetry_compatible(d2, register_symmetries(YAxisReflection())) === false
    # the half disc keeps only the x axis mirror
    half = hf_half_disc_hole()
    @test symmetry_compatible(half, register_symmetries(XAxisReflection())) === true
    @test symmetry_compatible(half, D2_symmetry()) === false
    @test symmetry_compatible(half, D4_symmetry()) === false
    # four discs placed by D2 are invariant under D2 but not under D4 (the diagonal moves the pattern)
    sh = hf_symmetrized_disc()
    @test symmetry_compatible(sh, D2_symmetry()) === true
    @test symmetry_compatible(sh, D4_symmetry()) === false
    @test symmetry_compatible(sh, register_symmetries(DiagonalReflection())) === false
end

@testset "compatible_symmetries (DomainHole)" begin
    # a hole without symmetry is only mapped to itself by the identity: no registered element is compatible
    for h in (hl_circle_hole(), hl_rect_hole(), hl_two_loop_hole())
        for reg in (register_symmetries(), D2_symmetry(), D4_symmetry(), register_symmetries(XAxisReflection()))
            c = compatible_symmetries(h, reg)
            @test c isa SymmetryRegistry
            @test length(c) == 0
        end
    end
end

@testset "compatible_symmetries (SymmetricHole)" begin
    ids(c) = [g.sym_id for g in c]
    # D4 disc at the origin: every element of a registry inside D4 is compatible, with its original sym_id and registry order
    d4 = hf_d4_symmetric_hole()
    full = compatible_symmetries(d4, D4_symmetry())
    @test full isa SymmetryRegistry
    @test length(full) == 7
    @test ids(full) == collect(1:7)
    @test typeof.(collect(full)) == typeof.(collect(D4_symmetry()))
    @test ids(compatible_symmetries(d4, D2_symmetry())) == [1, 2, 3]
    @test length(compatible_symmetries(d4, register_symmetries())) == 0
    # a group element outside the hole group (rotation by 2 pi/3 against D4) is dropped
    @test length(compatible_symmetries(d4, register_symmetries(NFoldRotation(3, 1)))) == 0

    # D2 disc about (0.5, 0): the key of Y and XY are in the hole group but do not fix the center, only X does
    d2 = hf_d2_symmetric_hole()
    @test ids(compatible_symmetries(d2, D2_symmetry())) == [3]
    @test typeof(only(collect(compatible_symmetries(d2, D2_symmetry())))) == XAxisReflection
    @test ids(compatible_symmetries(d2, D4_symmetry())) == [7]       # the X element of D4
    @test length(compatible_symmetries(d2, register_symmetries(YAxisReflection()))) == 0
    @test length(compatible_symmetries(d2, register_symmetries(XYAxisReflection()))) == 0

    # half disc about (0.3, 0): the hole group is {X}, so the keys of Y and XY are not in it
    half = hf_half_disc_hole()
    @test ids(compatible_symmetries(half, D2_symmetry())) == [3]
    @test ids(compatible_symmetries(half, D4_symmetry())) == [7]
    @test length(compatible_symmetries(half, register_symmetries(YAxisReflection()))) == 0

    # the hole group contains X but the center lies off the x axis: X does not fix the center
    off_axis_piece = CircleWedge(0.2, pi, 0.0, 0.3, 0.1, 3; bcs = AbsBoundaryCondition[SpecularReflection(), SymmetryWall(1, 2), SymmetryWall(1, 2)])
    off_axis = SymmetricHole(off_axis_piece, register_symmetries(XAxisReflection()), SVector(0.3, 0.1))
    @test length(compatible_symmetries(off_axis, register_symmetries(XAxisReflection()))) == 0
    @test symmetry_compatible(off_axis, register_symmetries(XAxisReflection())) === false

    # a mirror at a generic angle is not in the key set of the {X} hole group
    generic = AxisReflection{1 // 6}()
    @test length(compatible_symmetries(half, register_symmetries(generic))) == 0
    @test symmetry_compatible(half, register_symmetries(generic)) === false
    @test ids(compatible_symmetries(half, register_symmetries(XAxisReflection()))) == [1]
end

@testset "compatible_symmetries (SymmetrizedHoles)" begin
    ids(c) = [g.sym_id for g in c]
    # base DomainHole placed by D2: the pattern is invariant under the elements of D2 and no others of D4
    sh = hf_symmetrized_disc()
    @test ids(compatible_symmetries(sh, D2_symmetry())) == [1, 2, 3]
    @test ids(compatible_symmetries(sh, D4_symmetry())) == [3, 4, 7]    # Y, XY, X of D4
    @test length(compatible_symmetries(sh, register_symmetries())) == 0
    @test length(compatible_symmetries(sh, register_symmetries(DiagonalReflection()))) == 0
    @test length(compatible_symmetries(sh, register_symmetries(NFoldRotation(4, 1)))) == 0
    # placed by the x mirror alone: a pair of discs, invariant under X and not under Y
    pair_x = SymmetrizedHoles(DomainHole(hf_circle_domain(0.1, 0.4, 0.2, 5)), register_symmetries(XAxisReflection()))
    @test ids(compatible_symmetries(pair_x, D2_symmetry())) == [3]
    # placed by the y mirror alone: invariant under Y only
    pair_y = SymmetrizedHoles(DomainHole(hf_circle_domain(0.1, 0.4, 0.2, 5)), register_symmetries(YAxisReflection()))
    @test ids(compatible_symmetries(pair_y, D2_symmetry())) == [1]
    # empty placing registry: one disc without symmetry of its own, nothing but the empty registry is compatible
    single = SymmetrizedHoles(DomainHole(hf_circle_domain(0.1, 0.4, 0.2, 5)), register_symmetries())
    @test length(compatible_symmetries(single, D2_symmetry())) == 0
    @test symmetry_compatible(single, register_symmetries()) === true

    # base SymmetricHole: the copies are placed only, the object has no symmetry of its own
    placed = SymmetrizedHoles(hf_d2_symmetric_hole(), register_symmetries(YAxisReflection()))
    @test length(compatible_symmetries(placed, D2_symmetry())) == 0
    @test length(compatible_symmetries(placed, register_symmetries(YAxisReflection()))) == 0
    @test symmetry_compatible(placed, register_symmetries()) === true
    @test symmetry_compatible(placed, register_symmetries(YAxisReflection())) === false
end

@testset "symmetry_compatible (generic via compatible_symmetries)" begin
    # compatible exactly when the compatible subgroup is the whole registry
    holes = (hl_circle_hole(), hl_two_loop_hole(), hf_d4_symmetric_hole(), hf_d2_symmetric_hole(), hf_half_disc_hole(), hf_symmetrized_disc())
    regs = (register_symmetries(), D2_symmetry(), D4_symmetry(), register_symmetries(XAxisReflection()), register_symmetries(DiagonalReflection()))
    for h in holes, reg in regs
        @test symmetry_compatible(h, reg) == (length(compatible_symmetries(h, reg)) == length(reg))
        @test symmetry_compatible(h, reg) isa Bool
    end
end

@testset "_is_special_axis" begin
    # the four special axes 0, 1/4, 1/2, 3/4 (as mirror angle pi a)
    @test _is_special_axis(XAxisReflection())
    @test _is_special_axis(DiagonalReflection())
    @test _is_special_axis(YAxisReflection())
    @test _is_special_axis(AntiDiagonalReflection())
    # the angle is taken mod 1
    @test _is_special_axis(AxisReflection{5 // 4}())
    @test _is_special_axis(AxisReflection{1 // 1}())
    @test _is_special_axis(AxisReflection{-1 // 4}())
    # other angles are generic
    @test !_is_special_axis(AxisReflection{1 // 6}())
    @test !_is_special_axis(AxisReflection{1 // 3}())
    @test !_is_special_axis(AxisReflection{1 // 8}())
    @test !_is_special_axis(AxisReflection{5 // 6}())
end

@testset "_fixes_point" begin
    # special axes are exact: the point must lie on the line, any offset (even 1e-12) breaks it
    @test _fixes_point(XAxisReflection(), SVector(0.3, 0.0))
    @test !_fixes_point(XAxisReflection(), SVector(0.3, 1e-12))
    @test _fixes_point(YAxisReflection(), SVector(0.0, 0.7))
    @test !_fixes_point(YAxisReflection(), SVector(1e-12, 0.7))
    @test _fixes_point(DiagonalReflection(), SVector(0.4, 0.4))
    @test !_fixes_point(DiagonalReflection(), SVector(0.4, 0.4 + 1e-12))
    @test _fixes_point(AntiDiagonalReflection(), SVector(0.4, -0.4))
    @test !_fixes_point(AntiDiagonalReflection(), SVector(0.4, 0.4))
    # the pi rotation and the quarter turn fix the origin only
    @test _fixes_point(XYAxisReflection(), SVector(0.0, 0.0))
    @test !_fixes_point(XYAxisReflection(), SVector(0.3, 0.0))
    @test _fixes_point(NFoldRotation(4, 1), SVector(0.0, 0.0))
    @test !_fixes_point(NFoldRotation(4, 1), SVector(0.5, 0.0))
    @test _fixes_point(NFoldRotation(3, 2), SVector(0.0, 0.0))
    # a generic axis uses the tolerance sqrt(eps): points on the line at angle pi/6 and tiny offsets are fixed, a visible offset is not
    gen = AxisReflection{1 // 6}()
    on_line = 0.5 * SVector(cos(pi / 6), sin(pi / 6))
    normal = SVector(-sin(pi / 6), cos(pi / 6))
    @test _fixes_point(gen, SVector(0.0, 0.0))
    @test _fixes_point(gen, on_line)
    @test _fixes_point(gen, on_line + 1e-10 * normal)
    @test !_fixes_point(gen, on_line + 1e-6 * normal)
    @test !_fixes_point(gen, SVector(0.3, 0.0))
    # Float32 points: exact for special axes
    @test _fixes_point(XAxisReflection(), SVector(0.3f0, 0.0f0))
    @test !_fixes_point(XAxisReflection(), SVector(0.3f0, 1f-6))
end

@testset "_hole_block_keys" begin
    # a DomainHole is one block, the identity
    @test _hole_block_keys(hl_circle_hole()) == [_IDENTITY_KEY]
    @test _hole_block_keys(hl_two_loop_hole()) == [_IDENTITY_KEY]
    # SymmetricHole: the keys of its group in chain order, identity first
    d4 = hf_d4_symmetric_hole()
    @test _hole_block_keys(d4) == _block_keys(D4_symmetry())
    @test length(_hole_block_keys(d4)) == 8
    @test _hole_block_keys(d4)[1] == _IDENTITY_KEY
    @test length(_hole_block_keys(hf_d2_symmetric_hole())) == 4
    @test _hole_block_keys(hf_half_disc_hole()) == [_IDENTITY_KEY, symmetry_key(XAxisReflection())]
    # SymmetrizedHoles: the identity and the placing operations, for both base kinds
    @test _hole_block_keys(hf_symmetrized_disc()) == _block_keys(D2_symmetry())
    @test length(_hole_block_keys(hf_symmetrized_disc())) == hole_loop_count(hf_symmetrized_disc())
    @test _hole_block_keys(SymmetrizedHoles(hf_d2_symmetric_hole(), register_symmetries(YAxisReflection()))) ==
          [_IDENTITY_KEY, symmetry_key(YAxisReflection())]
end

@testset "_set_curve_bc (SymmetricHole)" begin
    h = hf_d4_symmetric_hole()
    # the id of an image piece (302 = piece id 2 + 3 * 100) edits the stored piece and rebuilds, so every copy changes
    new_h = _set_curve_bc(h, 302, 2, SymmetryWall(1, 9))
    @test new_h isa SymmetricHole{Float64}
    @test new_h.center == h.center && new_h.id_stride == h.id_stride
    seg2 = [c for c in get_all_curves(new_h) if c.segment_id == 2]
    @test length(seg2) == 8
    @test all(c -> c.bc == SymmetryWall(1, 9), seg2)
    @test [p.id for p in new_h.pieces] == [p.id for p in h.pieces]
    # the identity piece id works as well; the original hole is unchanged
    @test _set_curve_bc(h, 2, 2, SymmetryWall(1, 9)) isa SymmetricHole
    @test all(c -> c.bc == SymmetryWall(1, 2), [c for c in get_all_curves(h) if c.segment_id == 2])
    # no curve with that id pair: nothing
    @test _set_curve_bc(h, 999, 2, SymmetryWall(1, 9)) === nothing
    @test _set_curve_bc(h, 2, 9, SymmetryWall(1, 9)) === nothing
    # an ignored physical arc breaks the loop, which the constructor rejects
    @test_throws ArgumentError _set_curve_bc(h, 2, 1, QuantumSolverIgnore())
end

################################################################################
############################# SymmetrizedHoles #################################
################################################################################

@testset "SymmetrizedHoles constructor" begin
    sh = hf_symmetrized_disc()
    @test sh isa SymmetrizedHoles{DomainHole{SimpleDomain{Float64}}}
    @test sh isa AbsDomainHole
    @test !(sh isa AbsDomain)
    @test fieldnames(SymmetrizedHoles) == (:base, :symmetries, :holes, :id_stride)
    @test sh.id_stride == 100
    @test length(sh.symmetries) == 3
    @test sh.holes isa Vector{AbsDomainHole}
    @test length(sh.holes) == 4
    @test sh.holes[1] === sh.base
    @test all(c -> c isa DomainHole, sh.holes)
    # copies by the elements of the registry about the origin (Y, XY, X): fresh domain ids k * 100
    centers = [only(get_boundary_curves(c)).center for c in sh.holes]
    @test centers == [SVector(0.4, 0.2), SVector(-0.4, 0.2), SVector(-0.4, -0.2), SVector(0.4, -0.2)]
    @test [only(get_all_domains(c)).id for c in sh.holes] == [5, 105, 205, 305]
    @test all(c -> only(get_boundary_curves(c)).radius == 0.1, sh.holes)
    # keyword id_stride
    @test SymmetrizedHoles(DomainHole(hf_circle_domain(0.1, 0.4, 0.2, 5)), D2_symmetry(); id_stride = 7).id_stride == 7
    @test [only(get_all_domains(c)).id for c in SymmetrizedHoles(DomainHole(hf_circle_domain(0.1, 0.4, 0.2, 5)), D2_symmetry(); id_stride = 7).holes] == [5, 12, 19, 26]
    # empty registry: the base alone
    @test length(SymmetrizedHoles(DomainHole(hf_circle_domain(0.1, 0.4, 0.2, 5)), register_symmetries()).holes) == 1
end

@testset "SymmetrizedHoles constructor (SymmetricHole base)" begin
    # the quarter-disc hole about (0.5, 0) mirrored in the y axis: two discs, default stride 100 * (3 + 1)
    base = hf_d2_symmetric_hole()
    sh = SymmetrizedHoles(base, register_symmetries(YAxisReflection()))
    @test sh isa SymmetrizedHoles{SymmetricHole{Float64}}
    @test sh.id_stride == 400
    @test length(sh.holes) == 2
    @test sh.holes[1] === base
    @test sh.holes[2] isa DomainHole       # a placed copy carries no symmetry of its own
    @test [d.id for d in get_all_domains(sh)] == [2, 102, 202, 302, 402, 502, 602, 702]
    @test hole_loop_count(sh) == 2
    @test length(only(boundary_components(sh.holes[2]))) == 4
    @test area(sh) ≈ 2 * pi * 0.01 atol = hl_ATOL
end

@testset "SymmetrizedHoles constructor errors" begin
    # more than one loop in the base
    @test_throws ArgumentError SymmetrizedHoles(hl_two_loop_hole(), D2_symmetry())
    @test_throws "exactly one boundary loop" SymmetrizedHoles(hl_two_loop_hole(), D2_symmetry())
    # duplicate geometric actions in the registry
    @test_throws ArgumentError SymmetrizedHoles(hl_circle_hole(), register_symmetries(XAxisReflection(), XAxisReflection()))
end

@testset "hole interface (SymmetrizedHoles)" begin
    sh = hf_symmetrized_disc()
    @test hole_loop_count(sh) == 4
    @test hole_leaves(sh) isa Vector{AbsDomainHole}
    @test length(hole_leaves(sh)) == 1 && hole_leaves(sh)[1] === sh
    @test hole_domains(sh) isa Vector{AbsDomain}
    @test length(hole_domains(sh)) == 4
    @test [d.id for d in get_all_domains(sh)] == [5, 105, 205, 305]
    @test length(get_all_curves(sh)) == 4
    @test [c.domain_id for c in get_all_curves(sh)] == [5, 105, 205, 305]
    @test length(get_boundary_curves(sh)) == 4
    comps = boundary_components(sh)
    @test length(comps) == 4
    @test all(loop -> length(loop) == 1, comps)
    @test all(hl_loop_closed, comps)
    @test all(loop -> isapprox(hl_loop_length(loop), 2 * pi * 0.1; atol = hl_ATOL), comps)
    # every copy has its region on the left
    @test all(loop -> isapprox(hl_signed_area(loop), pi * 0.01; atol = hl_ATOL), comps)
    # the loops are consecutive in the order of `holes`
    @test [only(loop).center for loop in comps] == [only(get_boundary_curves(c)).center for c in sh.holes]
end

@testset "is_inside (SymmetrizedHoles)" begin
    sh = hf_symmetrized_disc()
    for c in (SVector(0.4, 0.2), SVector(-0.4, 0.2), SVector(-0.4, -0.2), SVector(0.4, -0.2))
        @test is_inside(sh, c) === true
        @test is_inside(sh, c + SVector(0.05, 0.0)) === true
        @test is_inside(sh, c + SVector(0.2, 0.0)) === false
    end
    @test is_inside(sh, SVector(0.0, 0.0)) === false
    @test is_inside(sh, SVector(0.4, -0.01)) === false
    @test is_inside(sh, [SVector(0.4, 0.2), SVector(0.0, 0.0), SVector(-0.4, -0.2)]) == [true, false, true]
    # the SymmetricHole base: both discs
    sh2 = SymmetrizedHoles(hf_d2_symmetric_hole(), register_symmetries(YAxisReflection()))
    @test is_inside(sh2, SVector(0.5, 0.05)) === true
    @test is_inside(sh2, SVector(-0.5, 0.05)) === true
    @test is_inside(sh2, SVector(0.0, 0.05)) === false
end

@testset "area (SymmetrizedHoles)" begin
    @test area(hf_symmetrized_disc()) ≈ 4 * pi * 0.01 atol = hl_ATOL
    @test area(SymmetrizedHoles(hf_d2_symmetric_hole(), register_symmetries(YAxisReflection()))) ≈ 2 * pi * 0.01 atol = hl_ATOL
    @test area(hf_symmetrized_disc(); rtol = 1e-6) ≈ 4 * pi * 0.01 atol = 1e-8
end

@testset "_set_curve_bc (SymmetrizedHoles)" begin
    sh = SymmetrizedHoles(hf_d2_symmetric_hole(), register_symmetries(YAxisReflection()))
    # id 502 is the image (offset 400) of base piece id 102; the edit goes through the base and rebuilds all copies
    new_sh = _set_curve_bc(sh, 502, 2, SymmetryWall(1, 9))
    @test new_sh isa SymmetrizedHoles
    @test new_sh.id_stride == sh.id_stride
    @test new_sh.base.piece.boundary[2].bc == SymmetryWall(1, 9)
    @test all(c -> c.bc == SymmetryWall(1, 9), [c for c in get_all_curves(new_sh) if c.segment_id == 2])
    @test all(c -> c.bc == SymmetryWall(1, 2), [c for c in get_all_curves(sh) if c.segment_id == 2])
    @test _set_curve_bc(sh, 999, 2, SymmetryWall(1, 9)) === nothing
    @test _set_curve_bc(hf_symmetrized_disc(), 999, 1, QuantumSolverIgnore()) === nothing
    # a disc base: the id of a copy (205 = 5 + 2 * 100) is found through the base
    @test_throws ArgumentError _set_curve_bc(hf_symmetrized_disc(), 205, 1, QuantumSolverIgnore())
end

################################################################################
############################# helpers and validation ###########################
################################################################################

@testset "_hole_parts" begin
    d = hl_circle_hole()
    @test _hole_parts(d) isa Vector{AbsDomainHole}
    @test length(_hole_parts(d)) == 1 && _hole_parts(d)[1] === d
    sym = hf_d4_symmetric_hole()
    @test length(_hole_parts(sym)) == 1 && _hole_parts(sym)[1] === sym
    sh = hf_symmetrized_disc()
    @test _hole_parts(sh) === sh.holes
    @test length(_hole_parts(sh)) == 4
end

@testset "_part_contains" begin
    atol = 1e-8
    d = hl_circle_hole()
    @test _part_contains(d, SVector(0.3, 0.4), atol) === true
    @test _part_contains(d, SVector(0.5, 0.4), atol) === true        # on the boundary
    @test _part_contains(d, SVector(0.501, 0.4), atol) === false
    h = hf_d4_symmetric_hole()
    @test _part_contains(h, SVector(0.1, 0.05), atol) === true
    @test _part_contains(h, SVector(0.3, 0.0), atol) === true        # on the boundary, on a mirror line
    @test _part_contains(h, SVector(0.31, 0.0), atol) === false
end

@testset "_apply_symmetry_to_region" begin
    # a simple domain: SimpleDomain with shifted ids
    rect = hf_rect_domain(0.0, 0.0, 1.0, 1.0, 3)
    img = _apply_symmetry_to_region(XAxisReflection(), rect, 10)
    @test img isa SimpleDomain
    @test img.id == 13
    @test all(c -> c.domain_id == 13, img.boundary)
    @test all(is_inside(img, SVector(0.5, -0.5)))
    @test !all(is_inside(img, SVector(0.5, 0.5)))
    # a composite domain: every subdomain is mapped
    comp = _apply_symmetry_to_region(XAxisReflection(), hl_two_loop_hole().domain, 10)
    @test comp isa CompositeDomain
    @test [sd.id for sd in comp.subdomains] == [13, 14]
    @test all(is_inside(comp.subdomains[2], SVector(2.5, -0.5)))
    @test !all(is_inside(comp.subdomains[2], SVector(2.5, 0.5)))
end

@testset "validate_holes (symmetric holes)" begin
    outer = hf_rect_domain(-1.0, -1.0, 1.0, 1.0, 1)
    # four placed discs inside the square, a symmetric hole centered in it
    @test validate_holes(MultiplyConnectedDomain(outer, [hf_symmetrized_disc()], 1), register_symmetries()) === nothing
    @test validate_holes(MultiplyConnectedDomain(outer, [hf_d4_symmetric_hole()], 1), register_symmetries()) === nothing
    @test validate_holes(MultiplyConnectedDomain(outer, [SymmetrizedHoles(hf_d2_symmetric_hole(), register_symmetries(YAxisReflection()))], 1), register_symmetries()) === nothing

    # a disc on a mirror line: its two copies coincide
    on_line = SymmetrizedHoles(DomainHole(hf_circle_domain(0.1, 0.4, 0.0, 5)), register_symmetries(XAxisReflection()))
    bad = MultiplyConnectedDomain(outer, [on_line], 1)
    @test_throws ArgumentError validate_holes(bad, register_symmetries())
    @test_throws "copies 1 and 2 of hole 1 overlap" validate_holes(bad, register_symmetries())
    # the hole index is the position in the hole list
    other = DomainHole(hf_circle_domain(0.1, -0.6, -0.6, 3))
    @test_throws "copies 1 and 2 of hole 2 overlap" validate_holes(MultiplyConnectedDomain(outer, [other, on_line], 1), register_symmetries())
    @test validate_holes(MultiplyConnectedDomain(outer, [other, hf_symmetrized_disc()], 1), register_symmetries()) === nothing
end

@testset "SymmetricHole in a MultiplyConnectedDomain" begin
    outer = hf_rect_domain(-1.0, -1.0, 1.0, 1.0, 1)
    dom = MultiplyConnectedDomain(outer, [hf_d4_symmetric_hole()], 1)
    @test genus(dom) == 1
    @test hole_domain_ids(dom) == [2, 102, 202, 302, 402, 502, 602, 702]
    @test component_of_domain_id(dom) == Dict(1 => 0, (id => 1 for id in 2:100:702)...)
    @test all(c -> is_hole_curve(dom, c), get_all_curves(dom.holes[1]))
    @test is_inside(dom, SVector(0.8, 0.0)) === true
    @test is_inside(dom, SVector(0.1, 0.0)) === false
    @test area(hf_HoleBilliard{Float64}(dom, register_symmetries())) ≈ 4 - pi * 0.09 atol = hl_ATOL
    # copies placed by a SymmetrizedHoles: one hole object, four loops
    dom4 = MultiplyConnectedDomain(outer, [hf_symmetrized_disc()], 1)
    @test genus(dom4) == 4
    @test hole_domain_ids(dom4) == [5, 105, 205, 305]
    @test component_of_domain_id(dom4) == Dict(1 => 0, 5 => 1, 105 => 1, 205 => 1, 305 => 1)
    # a piece id equal to an outer id is rejected
    clash = hf_rect_domain(-1.0, -1.0, 1.0, 1.0, 102)
    @test_throws ArgumentError MultiplyConnectedDomain(clash, [hf_d4_symmetric_hole()], 1)
end
