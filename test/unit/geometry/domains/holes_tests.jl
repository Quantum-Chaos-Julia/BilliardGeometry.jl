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

using StaticArrays: SVector
using BilliardGeometry: _set_curve_bc, _hole_parts, _part_contains, _apply_symmetry_to_region, _area_integral

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
    # the hole carries its own group and is never reduced by the billiard registry: compatible only with the empty registry
    for h in (hf_d4_symmetric_hole(), hf_d2_symmetric_hole(), hf_half_disc_hole(), hf_symmetrized_disc())
        @test symmetry_compatible(h, register_symmetries()) === true
        @test symmetry_compatible(h, D2_symmetry()) === false
        @test symmetry_compatible(h, D4_symmetry()) === false
    end
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
