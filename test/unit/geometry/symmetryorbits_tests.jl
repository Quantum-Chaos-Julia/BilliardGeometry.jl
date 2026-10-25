# Unit tests for the exact node permutations induced by a billiard's registered symmetry group on a periodic midpoint grid:
#   src/geometry/symmetryorbits.jl   (SymmetryOrbitMap and its accessors fundamental_size, full_size, length, orbit_size,
#                                     symmetry_orbit; symmetry_node_multiple; the exact group algebra _symmetry_key,
#                                     _compose_symmetry_keys, _rotation_key, _compose_rotation_keys, _symmetry_determinant,
#                                     _close_under_composition, _minimal_reflection_generators, _reflection_sector_keys,
#                                     _reflection_sector_permutation, _rotation_sector_keys, _rotation_sector_permutation;
#                                     get_symmetries; _sector_node_permutation, _boundary_symmetry_permutation,
#                                     _compose_index_permutations; _build_symmetry_orbit_map; symmetry_index_orbits)
#   src/geometry/symmetryregistry.jl (symmetry_of only)
# Reference objects: CircleBilliard(1.0) (D2: Y, XY, X registered in this order, so the four boundary sectors of the node grid
# are id, Y, XY, X), C3Billiard(0.2) (rotations by 2 pi/3 and 4 pi/3), SquareWithinSquareBilliard(1.0) (diagonal reflection)
# and a circular PolarBilliard with a single registered reflection. The permutations are integer algebra on the sector
# index: on N nodes the X reflection is q -> N + 1 - q, the Y reflection is q -> N/2 + 1 - q (cyclically), the rotation by
# 2 pi/n is q -> q + N/n (cyclically), and the diagonal and anti-diagonal reflections are q -> N + 1 - q. The orbit maps are
# therefore checked against hand-derived integer tables on small grids; the geometric meaning of the permutations on real
# boundary samples is covered by the integration tier. No eigen solves, no pinned wavefunction values. Helper functions and
# constants use the prefix `so_`.

using BilliardGeometry: AbsReflection, _symmetry_key, _compose_symmetry_keys, _rotation_key, _compose_rotation_keys,
                        _symmetry_determinant, _close_under_composition, _minimal_reflection_generators,
                        _reflection_sector_keys, _reflection_sector_permutation, _rotation_sector_keys,
                        _rotation_sector_permutation, _sector_node_permutation, _boundary_symmetry_permutation,
                        _compose_index_permutations, _build_symmetry_orbit_map, _SYMMETRY_IDENTITY_KEY, _rotation_identity_key

# Number of midpoint nodes of the small grids: a multiple of 4 (D2 reductions) and 8 (diagonal reductions).
const so_N = 8
# Number of midpoint nodes for the C3 grid, a multiple of 3 and 2 nodes per sector.
const so_N3 = 6

# Node permutations of the reflections on `so_N` nodes, as q -> image.
const so_RX = [8, 7, 6, 5, 4, 3, 2, 1]
const so_RY = [4, 3, 2, 1, 8, 7, 6, 5]
const so_RXY = [5, 6, 7, 8, 1, 2, 3, 4]

# Circular PolarBilliard (unit radius, center at the origin) whose registry holds exactly the given symmetries, with the
# angular origin rotated by `shift`. The registry is closed for one reflection or for a full group, which is all the index
# algebra inspects; the shape is irrelevant for it.
function so_polar_billiard(syms::AbsSymmetry...; shift = -pi / 4)
    seg = FourierCoeffPolarSegment(Float64[]; shift_angle = shift)
    dom = PolarDomain{Float64}([seg], [curve(seg, 0.0)], 1)
    return PolarBilliard{Float64}(dom, register_symmetries(syms...))
end

# Fundamental-orbit member counts: every orbit of `orbits` has exactly `orbit_size` members.
function so_member_counts(orbits)
    return [count(==(b), orbits.orbit_of) for b in 1:fundamental_size(orbits)]
end

@testset "symmetry_node_multiple" begin
    @test symmetry_node_multiple(XAxisReflection()) == 4
    @test symmetry_node_multiple(YAxisReflection()) == 4
    @test symmetry_node_multiple(XYAxisReflection()) == 4
    @test symmetry_node_multiple(DiagonalReflection()) == 8
    @test symmetry_node_multiple(AntiDiagonalReflection()) == 8
    @test symmetry_node_multiple(NFoldRotation(6, 1)) == 6
    @test symmetry_node_multiple(NFoldRotation(3, 2)) == 3
    # a composite reflection needs the least common multiple of its generators
    @test symmetry_node_multiple(CompositeReflection(XAxisReflection(), YAxisReflection())) == 4
    @test symmetry_node_multiple(CompositeReflection(XAxisReflection(), DiagonalReflection())) == 8
end

@testset "symmetry keys and composition" begin
    @test _SYMMETRY_IDENTITY_KEY == (1, 0, 0, 1)
    @test _symmetry_key(XAxisReflection()) == (1, 0, 0, -1)
    @test _symmetry_key(YAxisReflection()) == (-1, 0, 0, 1)
    @test _symmetry_key(XYAxisReflection()) == (-1, 0, 0, -1)
    @test _symmetry_key(DiagonalReflection()) == (0, 1, 1, 0)
    @test _symmetry_key(AntiDiagonalReflection()) == (0, -1, -1, 0)

    # the integer keys multiply like the matrices of the generators
    kx, ky = _symmetry_key(XAxisReflection()), _symmetry_key(YAxisReflection())
    @test _compose_symmetry_keys(kx, ky) == _symmetry_key(XYAxisReflection())
    @test _compose_symmetry_keys(ky, kx) == _symmetry_key(XYAxisReflection())
    @test _compose_symmetry_keys(kx, kx) == _SYMMETRY_IDENTITY_KEY
    @test _compose_symmetry_keys(_symmetry_key(DiagonalReflection()), _symmetry_key(AntiDiagonalReflection())) ==
          _symmetry_key(XYAxisReflection())
    for sym in (XAxisReflection(), YAxisReflection(), XYAxisReflection(), DiagonalReflection(), AntiDiagonalReflection())
        @test _compose_symmetry_keys(_SYMMETRY_IDENTITY_KEY, _symmetry_key(sym)) == _symmetry_key(sym)
    end

    # determinant -1 for pure mirror lines, +1 for the pi rotation XY
    @test _symmetry_determinant(XAxisReflection()) == -1
    @test _symmetry_determinant(YAxisReflection()) == -1
    @test _symmetry_determinant(DiagonalReflection()) == -1
    @test _symmetry_determinant(AntiDiagonalReflection()) == -1
    @test _symmetry_determinant(XYAxisReflection()) == 1

    # closure of {identity} under X and Y is the full D2 group
    closure = _close_under_composition(Set([_SYMMETRY_IDENTITY_KEY]), [kx, ky])
    @test closure == Set([_SYMMETRY_IDENTITY_KEY, kx, ky, _symmetry_key(XYAxisReflection())])
    # closure under one reflection stays of order 2
    @test length(_close_under_composition(Set([_SYMMETRY_IDENTITY_KEY]), [kx])) == 2
end

@testset "rotation keys and composition" begin
    @test _rotation_key(NFoldRotation(3, 1)) == (3, 1)
    @test _rotation_key(NFoldRotation(3, 4)) == (3, 1)      # the power is reduced modulo the order
    @test _rotation_key(NFoldRotation(5, 5)) == (5, 0)
    @test _rotation_identity_key(7) == (7, 0)
    @test _compose_rotation_keys((3, 1), (3, 1)) == (3, 2)
    @test _compose_rotation_keys((3, 2), (3, 2)) == (3, 1)
    @test _compose_rotation_keys((3, 1), (3, 2)) == (3, 0)
    @test_throws ArgumentError _compose_rotation_keys((3, 1), (4, 1))
end

@testset "reflection sector keys and permutations" begin
    circ = CircleBilliard(1.0)
    # identity followed by the registry (Y, XY, X)
    @test _reflection_sector_keys(circ) == [_SYMMETRY_IDENTITY_KEY, _symmetry_key(YAxisReflection()),
                                            _symmetry_key(XYAxisReflection()), _symmetry_key(XAxisReflection())]
    # a billiard with no symmetry has only the identity sector
    @test _reflection_sector_keys(TriangleBilliard(1.0, 1.0)) == [_SYMMETRY_IDENTITY_KEY]

    # left multiplication permutes the sectors (id, Y, XY, X); every row is a permutation of 1:4
    @test _reflection_sector_permutation(circ, XAxisReflection()) == [4, 3, 2, 1]
    @test _reflection_sector_permutation(circ, YAxisReflection()) == [2, 1, 4, 3]
    @test _reflection_sector_permutation(circ, XYAxisReflection()) == [3, 4, 1, 2]
    # an element outside the registered group is rejected
    @test_throws ArgumentError _reflection_sector_permutation(circ, DiagonalReflection())
    @test_throws ArgumentError _reflection_sector_permutation(TriangleBilliard(1.0, 1.0), XAxisReflection())
    @test_throws ArgumentError _reflection_sector_permutation(circ, CompositeReflection(XAxisReflection(), YAxisReflection()))

    # registries that are not reflection groups, or that repeat an action, are rejected
    @test_throws ArgumentError _reflection_sector_keys(C3Billiard(0.2))
    @test_throws ArgumentError _reflection_sector_keys(so_polar_billiard(XAxisReflection(), XAxisReflection()))
end

@testset "rotation sector keys and permutations" begin
    c3 = C3Billiard(0.2)
    @test _rotation_sector_keys(c3) == [(3, 0), (3, 1), (3, 2)]
    @test _rotation_sector_permutation(c3, NFoldRotation(3, 1)) == [2, 3, 1]
    @test _rotation_sector_permutation(c3, NFoldRotation(3, 2)) == [3, 1, 2]
    # a rotation of another order is not in the registered group
    @test_throws ArgumentError _rotation_sector_permutation(c3, NFoldRotation(4, 1))

    # no registry, a reflection registry, mixed registries and repeated rotations are rejected
    @test_throws ArgumentError _rotation_sector_keys(TriangleBilliard(1.0, 1.0))
    @test_throws ArgumentError _rotation_sector_keys(CircleBilliard(1.0))
    @test_throws ArgumentError _rotation_sector_keys(so_polar_billiard(NFoldRotation(3, 1), XAxisReflection()))
    @test_throws ArgumentError _rotation_sector_keys(so_polar_billiard(NFoldRotation(3, 1), NFoldRotation(4, 1)))
    @test_throws ArgumentError _rotation_sector_keys(so_polar_billiard(NFoldRotation(3, 1), NFoldRotation(3, 1)))
end

@testset "_sector_node_permutation" begin
    # orientation-preserving swap of two sectors shifts the nodes by half the grid
    @test _sector_node_permutation(8, [2, 1], false) == [5, 6, 7, 8, 1, 2, 3, 4]
    # orientation-reversing swap reverses the whole grid
    @test _sector_node_permutation(8, [2, 1], true) == [8, 7, 6, 5, 4, 3, 2, 1]
    # identity sector permutation: the node index is kept, or reversed inside each block
    @test _sector_node_permutation(6, [1, 2, 3], false) == collect(1:6)
    @test _sector_node_permutation(6, [1, 2, 3], true) == [2, 1, 4, 3, 6, 5]
    # the node count must be divisible by the number of sectors
    @test_throws ArgumentError _sector_node_permutation(7, [2, 1], false)
end

@testset "_boundary_symmetry_permutation" begin
    circ = CircleBilliard(1.0)
    @test _boundary_symmetry_permutation(circ, so_N, XAxisReflection()) == so_RX
    @test _boundary_symmetry_permutation(circ, so_N, YAxisReflection()) == so_RY
    @test _boundary_symmetry_permutation(circ, so_N, XYAxisReflection()) == so_RXY
    # a rotation preserves the local node index, so on C3 it shifts the nodes by N/3
    c3 = C3Billiard(0.2)
    @test _boundary_symmetry_permutation(c3, so_N3, NFoldRotation(3, 1)) == [3, 4, 5, 6, 1, 2]
    @test _boundary_symmetry_permutation(c3, so_N3, NFoldRotation(3, 2)) == [5, 6, 1, 2, 3, 4]
    # a composite reflection is a set of geometric actions, not one
    @test_throws ArgumentError _boundary_symmetry_permutation(circ, so_N, CompositeReflection(XAxisReflection(), YAxisReflection()))
end

@testset "_compose_index_permutations" begin
    @test _compose_index_permutations([2, 3, 1], [3, 1, 2]) == [1, 2, 3]
    # the result applies the second argument first: (a o b)[q] = a[b[q]]
    @test _compose_index_permutations([2, 1, 3], [1, 3, 2]) == [2, 3, 1]
    @test _compose_index_permutations(so_RX, so_RY) == so_RXY
    @test _compose_index_permutations(so_RX, so_RX) == collect(1:so_N)
    @test_throws DimensionMismatch _compose_index_permutations([1, 2], [1, 2, 3])
end

@testset "_build_symmetry_orbit_map" begin
    # two permutations of four nodes: identity and the swap of neighbours
    orbits = _build_symmetry_orbit_map(Float64, 4, Vector{Vector{Int}}([[1, 2, 3, 4], [2, 1, 4, 3]]))
    @test orbits isa SymmetryOrbitMap{Float64}
    @test orbits.fundamental_indices == [1, 3]
    @test orbits.orbit_of == [1, 1, 2, 2]
    @test orbits.fund_to_full == [1 3; 2 4]
    @test all(==(one(ComplexF64)), orbits.phase)            # trivial irrep factors by default
    @test orbits.fund_to_scale == ones(ComplexF64, 2, 2)

    # explicit irrep factors are attached to each image
    signed = _build_symmetry_orbit_map(Float64, 4, Vector{Vector{Int}}([[1, 2, 3, 4], [2, 1, 4, 3]]), ComplexF64[1, -1])
    @test signed.phase == ComplexF64[1, -1, 1, -1]
    @test signed.fund_to_scale == ComplexF64[1 1; -1 -1]

    # the scalar type follows the requested type
    @test _build_symmetry_orbit_map(Float32, 4, Vector{Vector{Int}}([[1, 2, 3, 4], [2, 1, 4, 3]])) isa SymmetryOrbitMap{Float32}

    # input validation
    @test_throws ArgumentError _build_symmetry_orbit_map(Float64, 5, Vector{Vector{Int}}([[1, 2, 3, 4, 5], [2, 1, 4, 3, 5]]))
    @test_throws DimensionMismatch _build_symmetry_orbit_map(Float64, 4, Vector{Vector{Int}}([[1, 2, 3, 4], [2, 1, 4, 3]]), ComplexF64[1])
    @test_throws DimensionMismatch _build_symmetry_orbit_map(Float64, 4, Vector{Vector{Int}}([[1, 2, 3, 4], [2, 1, 4]]))
    @test_throws ArgumentError _build_symmetry_orbit_map(Float64, 4, Vector{Vector{Int}}([[1, 2, 3, 4], [2, 1, 9, 4]]))
    # two images that coincide do not give free orbits of the group order
    @test_throws ArgumentError _build_symmetry_orbit_map(Float64, 4, Vector{Vector{Int}}([[1, 2, 3, 4], [1, 2, 3, 4]]))
end

@testset "SymmetryOrbitMap accessors" begin
    circ = CircleBilliard(1.0)
    orbits = symmetry_index_orbits(Float64, circ, so_N, XYAxisReflection())
    @test orbits isa SymmetryOrbitMap{Float64}
    @test fundamental_size(orbits) == orbits.fundamental_size == so_N ÷ 4
    @test full_size(orbits) == orbits.full_size == so_N
    @test length(orbits) == so_N
    @test orbit_size(orbits) == 4 == size(orbits.fund_to_full, 1)
    @test size(orbits.fund_to_full) == (4, so_N ÷ 4)
    @test size(orbits.fund_to_scale) == (4, so_N ÷ 4)
    @test length(orbits.fundamental_indices) == so_N ÷ 4
    @test length(orbits.orbit_of) == length(orbits.phase) == so_N
    @test @inferred(fundamental_size(orbits)) isa Int
    @test @inferred(orbit_size(orbits)) isa Int

    # symmetry_orbit returns the columns of the orbit tables
    for b in 1:fundamental_size(orbits)
        qs, χs = symmetry_orbit(orbits, b)
        @test qs == orbits.fund_to_full[:, b]
        @test χs == orbits.fund_to_scale[:, b]
        @test qs[1] == orbits.fundamental_indices[b]       # the representative is the identity image
    end

    # the scalar type of the irrep factors follows the requested type
    @test symmetry_index_orbits(Float32, circ, so_N, XAxisReflection()) isa SymmetryOrbitMap{Float32}
    @test eltype(symmetry_index_orbits(Float32, circ, so_N, XAxisReflection()).phase) == ComplexF32
end

@testset "symmetry_index_orbits (XAxisReflection)" begin
    circ = CircleBilliard(1.0)
    orbits = symmetry_index_orbits(Float64, circ, so_N, XAxisReflection())
    @test fundamental_size(orbits) == so_N ÷ 2
    @test length(orbits) == so_N
    @test orbit_size(orbits) == 2
    @test orbits.fundamental_indices == [1, 2, 3, 4]
    @test orbits.fund_to_full == [1 2 3 4; 8 7 6 5]
    @test orbits.orbit_of == [1, 2, 3, 4, 4, 3, 2, 1]
    @test all(==(2), so_member_counts(orbits))
    @test all(==(one(ComplexF64)), orbits.phase)

    # explicit antisymmetric character: representatives keep phase 1, their partners get -1
    odd = symmetry_index_orbits(Float64, circ, so_N, XAxisReflection(), ComplexF64(-1))
    @test all(==(one(ComplexF64)), odd.phase[odd.fundamental_indices])
    partners = setdiff(1:so_N, odd.fundamental_indices)
    @test all(==(ComplexF64(-1)), odd.phase[partners])
    @test odd.orbit_of == orbits.orbit_of                  # the character does not change the folding
end

@testset "symmetry_index_orbits (YAxisReflection)" begin
    circ = CircleBilliard(1.0)
    orbits = symmetry_index_orbits(Float64, circ, so_N, YAxisReflection())
    @test fundamental_size(orbits) == so_N ÷ 2
    @test length(orbits) == so_N
    @test orbit_size(orbits) == 2
    @test orbits.fundamental_indices == [1, 2, 5, 6]
    @test orbits.fund_to_full == [1 2 5 6; 4 3 8 7]
    @test orbits.orbit_of == [1, 2, 2, 1, 3, 4, 4, 3]
    @test all(==(2), so_member_counts(orbits))

    odd = symmetry_index_orbits(Float64, circ, so_N, YAxisReflection(), ComplexF64(-1))
    @test all(==(one(ComplexF64)), odd.phase[odd.fundamental_indices])
    partners = setdiff(1:so_N, odd.fundamental_indices)
    @test all(==(ComplexF64(-1)), odd.phase[partners])
end

@testset "symmetry_index_orbits (XYAxisReflection)" begin
    circ = CircleBilliard(1.0)
    orbits = symmetry_index_orbits(Float64, circ, so_N, XYAxisReflection())
    @test fundamental_size(orbits) == so_N ÷ 4
    @test length(orbits) == so_N
    @test orbit_size(orbits) == 4
    @test orbits.fundamental_indices == [1, 2]
    # images in the order identity, X, Y, XY
    @test orbits.fund_to_full == [1 2; 8 7; 4 3; 5 6]
    @test all(==(4), so_member_counts(orbits))
    @test all(==(one(ComplexF64)), orbits.phase)

    # each orbit member count and each member's recorded orbit and phase agree
    for b in 1:fundamental_size(orbits)
        qs, χs = symmetry_orbit(orbits, b)
        @test orbit_size(orbits) == length(qs) == 4
        for (q, χ) in zip(qs, χs)
            @test orbits.orbit_of[q] == b
            @test orbits.phase[q] == χ
        end
    end

    # independent characters for X and Y: the XY image carries their product
    mixed = symmetry_index_orbits(Float64, circ, so_N, XYAxisReflection(), ComplexF64(-1), ComplexF64(1))
    @test mixed.fund_to_scale[:, 1] == ComplexF64[1, -1, 1, -1]
    both = symmetry_index_orbits(Float64, circ, so_N, XYAxisReflection(), ComplexF64(-1), ComplexF64(-1))
    @test both.fund_to_scale[:, 1] == ComplexF64[1, -1, -1, 1]
end

@testset "symmetry_index_orbits (CompositeReflection)" begin
    circ = CircleBilliard(1.0)
    comp = CompositeReflection(XAxisReflection(), YAxisReflection())
    orbits = symmetry_index_orbits(Float64, circ, so_N, comp)
    xy = symmetry_index_orbits(Float64, circ, so_N, XYAxisReflection())
    # the closure of X and Y is the same D2 folding as the native XYAxisReflection
    @test fundamental_size(orbits) == fundamental_size(xy) == so_N ÷ 4
    @test orbit_size(orbits) == orbit_size(xy) == 4
    @test orbits.fundamental_indices == xy.fundamental_indices
    @test orbits.orbit_of == xy.orbit_of
    @test orbits.phase == xy.phase
    @test orbits.fund_to_full == xy.fund_to_full
    @test orbits.fund_to_scale == xy.fund_to_scale
    @test all(==(4), so_member_counts(orbits))

    # an XYAxisReflection entry expands into its two axis generators
    xyc = symmetry_index_orbits(Float64, circ, so_N, CompositeReflection(XYAxisReflection()))
    @test xyc.fund_to_full == xy.fund_to_full

    # one character per generator; the group element X*Y carries the product
    chars = symmetry_index_orbits(Float64, circ, so_N, comp, ComplexF64[-1, 1])
    @test chars.fund_to_scale[:, 1] == ComplexF64[1, -1, 1, -1]

    # invalid input
    @test_throws ArgumentError symmetry_index_orbits(Float64, circ, so_N, CompositeReflection(AbsReflection[]))
    @test_throws DimensionMismatch symmetry_index_orbits(Float64, circ, so_N, comp, ComplexF64[1])
    @test_throws ArgumentError symmetry_index_orbits(Float64, circ, 6, comp)
    # the same generator with conflicting characters is an inconsistent parity assignment
    @test_throws ArgumentError symmetry_index_orbits(Float64, circ, so_N, CompositeReflection(XAxisReflection(), XAxisReflection()), ComplexF64[1, -1])
end

@testset "symmetry_index_orbits (NFoldRotation)" begin
    c3 = C3Billiard(0.2)
    rot = NFoldRotation(3, 1)
    orbits = symmetry_index_orbits(Float64, c3, so_N3, rot)
    @test fundamental_size(orbits) == so_N3 ÷ 3
    @test length(orbits) == so_N3
    @test orbit_size(orbits) == 3
    @test orbits.fundamental_indices == [1, 2]
    @test orbits.fund_to_full == [1 2; 3 4; 5 6]
    @test orbits.orbit_of == [1, 2, 1, 2, 1, 2]
    @test all(==(3), so_member_counts(orbits))
    @test all(==(one(ComplexF64)), orbits.phase)           # sector 0 is the trivial representation

    # rotational sector l: the image under rot^j carries the phase exp(2 pi i l j / 3)
    sector1 = symmetry_index_orbits(Float64, c3, so_N3, rot, 1)
    ω = cis(2pi / 3)
    @test sector1.fund_to_scale[:, 1] ≈ ComplexF64[1, ω, ω^2]
    @test sector1.fund_to_full == orbits.fund_to_full

    # the node count must be a multiple of the rotation order; the rotation must be registered
    @test_throws ArgumentError symmetry_index_orbits(Float64, c3, 4, rot)
    @test_throws ArgumentError symmetry_index_orbits(Float64, c3, 8, NFoldRotation(4, 1))
    @test_throws ArgumentError symmetry_index_orbits(Float64, CircleBilliard(1.0), so_N3, rot)
end

@testset "symmetry_index_orbits (DiagonalReflection)" begin
    sws = SquareWithinSquareBilliard(1.0)
    orbits = symmetry_index_orbits(Float64, sws, so_N, DiagonalReflection())
    @test fundamental_size(orbits) == so_N ÷ 2
    @test length(orbits) == so_N
    @test orbit_size(orbits) == 2
    @test orbits.fund_to_full == [1 2 3 4; 8 7 6 5]
    @test all(==(2), so_member_counts(orbits))
    # the diagonal reflection needs a multiple of 8 nodes
    @test_throws ArgumentError symmetry_index_orbits(Float64, sws, 12, DiagonalReflection())
    # and it has to be registered
    @test_throws ArgumentError symmetry_index_orbits(Float64, CircleBilliard(1.0), so_N, DiagonalReflection())
end

@testset "symmetry_index_orbits (AntiDiagonalReflection)" begin
    # no shipped billiard registers this reflection alone; a circle with a shifted angular origin does
    circ_ad = so_polar_billiard(AntiDiagonalReflection())
    orbits = symmetry_index_orbits(Float64, circ_ad, so_N, AntiDiagonalReflection())
    @test fundamental_size(orbits) == so_N ÷ 2
    @test length(orbits) == so_N
    @test orbit_size(orbits) == 2
    @test orbits.fund_to_full == [1 2 3 4; 8 7 6 5]
    @test all(==(2), so_member_counts(orbits))
    @test_throws ArgumentError symmetry_index_orbits(Float64, circ_ad, 12, AntiDiagonalReflection())
    @test_throws ArgumentError symmetry_index_orbits(Float64, circ_ad, so_N, XAxisReflection())
end

@testset "symmetry_index_orbits (partition and counting)" begin
    circ = CircleBilliard(1.0)
    orbits = symmetry_index_orbits(Float64, circ, so_N, XYAxisReflection())
    @test fundamental_size(orbits) * orbit_size(orbits) == full_size(orbits)
    @test sum(b -> count(==(b), orbits.orbit_of), 1:fundamental_size(orbits)) == full_size(orbits)
    # the orbits cover every node exactly once
    all_members = Int[]
    for b in 1:fundamental_size(orbits)
        qs, _ = symmetry_orbit(orbits, b)
        append!(all_members, qs)
    end
    @test sort(all_members) == collect(1:full_size(orbits))
    @test length(all_members) == length(unique(all_members))

    # the node count must be a multiple of the reduction's node multiple
    @test_throws ArgumentError symmetry_index_orbits(Float64, circ, 6, XAxisReflection())
    @test_throws ArgumentError symmetry_index_orbits(Float64, circ, 6, YAxisReflection())
    @test_throws ArgumentError symmetry_index_orbits(Float64, circ, 6, XYAxisReflection())
    # a reflection outside the registry (no symmetry at all) is rejected
    @test_throws ArgumentError symmetry_index_orbits(Float64, TriangleBilliard(1.0, 1.0), so_N, XAxisReflection())
end

@testset "_minimal_reflection_generators" begin
    # D2 registry (Y, XY, X): the two mirror lines generate the group, XY is not needed
    gens = _minimal_reflection_generators(CircleBilliard(1.0).symmetries)
    @test gens isa Tuple
    @test typeof.(gens) == (YAxisReflection, XAxisReflection)
    # a single element is its own generator, also when it is a pi rotation
    @test typeof.(_minimal_reflection_generators(so_polar_billiard(XYAxisReflection()).symmetries)) == (XYAxisReflection,)
    @test typeof.(_minimal_reflection_generators(so_polar_billiard(DiagonalReflection()).symmetries)) == (DiagonalReflection,)
    # two mirror lines without their product are not a group
    @test_throws ArgumentError _minimal_reflection_generators(so_polar_billiard(XAxisReflection(), YAxisReflection()).symmetries)
end

@testset "get_symmetries" begin
    # D2 billiard: three registered elements, two generators (the mirror lines, not XY)
    gens_d2 = get_symmetries(RectangleBilliard(1.0, 0.6))
    @test length(gens_d2) == 2
    @test all(g -> g isa YAxisReflection || g isa XAxisReflection, gens_d2)
    @test !any(g -> g isa XYAxisReflection, gens_d2)
    @test gens_d2[1].sym_id == 1 && gens_d2[2].sym_id == 3         # keeps the registered sym_id tags

    # cyclic C3 billiard: generated by the rotation m = 1 alone
    gens_c3 = get_symmetries(C3Billiard(0.2))
    @test length(gens_c3) == 1
    @test only(gens_c3) isa NFoldRotation
    @test only(gens_c3).m == 1

    # single registered reflection
    @test get_symmetries(SquareWithinSquareBilliard(1.0)) isa Tuple{DiagonalReflection}

    # no registered symmetry: empty generating set
    @test get_symmetries(TriangleBilliard(1.0, 1.0)) == ()
end

@testset "symmetry_of" begin
    reg = RectangleBilliard(1.0, 0.6).symmetries
    @test symmetry_of(reg, 1) isa YAxisReflection
    @test symmetry_of(reg, 2) isa XYAxisReflection
    @test symmetry_of(reg, 3) isa XAxisReflection
    @test all(i -> symmetry_of(reg, i).sym_id == i, 1:3)
    # an unknown sym_id raises instead of returning nothing or a wrong generator
    @test_throws ArgumentError symmetry_of(reg, 999)
    @test_throws ArgumentError symmetry_of(reg, 0)
end
