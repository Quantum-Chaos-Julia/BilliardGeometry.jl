# Unit tests for the exact node permutations induced by a billiard's registered symmetry group on a periodic midpoint grid:
#   src/geometry/symmetryorbits.jl   (SymmetryOrbitMap and its accessors fundamental_size, full_size, length, orbit_size,
#                                     symmetry_orbit; symmetry_node_multiple; the exact dihedral group algebra _IDENTITY_KEY,
#                                     _compose_keys, _keys_closure, _block_keys, _block_permutation, _index_combinations,
#                                     _minimal_generators; get_symmetries, get_symmetry_group; _sector_node_permutation,
#                                     _block_node_permutation, _boundary_symmetry_permutation, _compose_index_permutations;
#                                     _build_symmetry_orbit_map, _orbit_map_from_generators; symmetry_index_orbits)
#   src/geometry/symmetry.jl         (symmetry_key, used as the label of each registered element)
#   src/geometry/symmetryregistry.jl (symmetry_of only)
# Reference objects: CircleBilliard(1.0) (D2: Y, XY, X registered in this order, so the four boundary sectors of the node grid
# are id, Y, XY, X), C3Billiard(0.2) (rotations by 2 pi/3 and 4 pi/3) and circular PolarBilliards whose registry is chosen freely
# (single diagonal or anti-diagonal reflection, D4_symmetry(), Cn_symmetry(6), mixed rotation orders); the shape is irrelevant for
# the index algebra. The dihedral key (reflects, a) labels an element by its action: rot(a) is the rotation by 2 pi a, refl(a) the
# mirror at angle pi a, and keys compose as rot(a)rot(b) = rot(a+b), rot(a)refl(b) = refl(a+b), refl(a)rot(b) = refl(a-b),
# refl(a)refl(b) = rot(a-b) with the angle taken mod 1; the key algebra is cross-checked against the 2x2 matrices of the elements.
# The permutations are integer algebra on the sector index: on N nodes the X reflection is q -> N + 1 - q, the Y reflection is
# q -> N/2 + 1 - q (cyclically), the rotation by 2 pi/n is q -> q + N/n (cyclically), and the diagonal and anti-diagonal
# reflections are q -> N + 1 - q. The orbit maps are therefore checked against hand-derived integer tables on small grids
# (D4 composite values come from an observed run, marked as such); the geometric meaning of the permutations on real boundary
# samples is covered by the integration tier. No eigen solves, no pinned wavefunction values. Helper functions and constants use
# the prefix `so_`.

using BilliardGeometry: AbsReflection, symmetry_key, _IDENTITY_KEY, _compose_keys, _keys_closure, _block_keys, _block_permutation,
                        _index_combinations, _minimal_generators, _sector_node_permutation, _block_node_permutation,
                        _boundary_symmetry_permutation, _compose_index_permutations, _build_symmetry_orbit_map,
                        _orbit_map_from_generators, _sym_matrix, symmetry_reduction_factor

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

# PolarBilliard whose registry is the D4 group in chain order (Diagonal, rotation by pi/2, YAxis, XY, AntiDiagonal, rotation by
# 3 pi/2, XAxis).
so_d4_billiard() = so_polar_billiard(collect(D4_symmetry())...)

# Dihedral keys: rotation by 2 pi a and mirror at angle pi a.
so_rot(a) = (false, a)
so_refl(a) = (true, a)

# The eight dihedral keys of D4: identity followed by the registered elements.
const so_D4_KEYS = vcat([_IDENTITY_KEY], [symmetry_key(s) for s in D4_symmetry()])

# 2x2 matrix of any registered symmetry element.
so_matrix(sym::NFoldRotation) = sym.sym_map.linear
so_matrix(sym::AbsReflection) = _sym_matrix(sym)

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

@testset "_compose_keys" begin
    @test _IDENTITY_KEY == (false, 0 // 1)

    # the four composition rules, with angles mod 1
    @test _compose_keys(so_rot(1 // 3), so_rot(1 // 3)) == so_rot(2 // 3)
    @test _compose_keys(so_rot(1 // 3), so_rot(2 // 3)) == _IDENTITY_KEY
    @test _compose_keys(so_rot(3 // 4), so_rot(1 // 2)) == so_rot(1 // 4)                 # 5/4 mod 1
    @test _compose_keys(so_rot(1 // 4), so_refl(1 // 4)) == so_refl(1 // 2)
    @test _compose_keys(so_refl(1 // 4), so_rot(1 // 4)) == so_refl(0 // 1)
    @test _compose_keys(so_refl(0 // 1), so_rot(1 // 4)) == so_refl(3 // 4)               # -1/4 mod 1
    @test _compose_keys(so_refl(1 // 4), so_refl(3 // 4)) == so_rot(1 // 2)
    @test _compose_keys(so_refl(3 // 4), so_refl(1 // 4)) == so_rot(1 // 2)
    @test _compose_keys(so_refl(1 // 4), so_refl(1 // 4)) == _IDENTITY_KEY

    # named reflections: X and Y give the pi rotation XY in either order, so do the two diagonals
    kx, ky = symmetry_key(XAxisReflection()), symmetry_key(YAxisReflection())
    kxy = symmetry_key(XYAxisReflection())
    @test _compose_keys(kx, ky) == kxy == so_rot(1 // 2)
    @test _compose_keys(ky, kx) == kxy
    @test _compose_keys(kx, kx) == _IDENTITY_KEY
    @test _compose_keys(symmetry_key(DiagonalReflection()), symmetry_key(AntiDiagonalReflection())) == kxy
    # rot(1/4) after the diagonal mirror is the y axis mirror, the diagonal mirror after rot(1/4) the x axis mirror
    @test _compose_keys(symmetry_key(NFoldRotation(4, 1)), symmetry_key(DiagonalReflection())) == symmetry_key(YAxisReflection())
    @test _compose_keys(symmetry_key(DiagonalReflection()), symmetry_key(NFoldRotation(4, 1))) == kx
    @test (@inferred _compose_keys(kx, ky)) isa Tuple{Bool,Rational{Int}}

    # the identity is neutral on both sides and every element has an inverse inside D4
    for k in so_D4_KEYS
        @test _compose_keys(_IDENTITY_KEY, k) == k
        @test _compose_keys(k, _IDENTITY_KEY) == k
        @test any(h -> _compose_keys(k, h) == _IDENTITY_KEY, so_D4_KEYS)
    end

    # D4 is closed under composition (also under squaring) and the composition is associative
    @test all(g -> all(h -> _compose_keys(g, h) in so_D4_KEYS, so_D4_KEYS), so_D4_KEYS)
    @test all(g -> _compose_keys(g, g) in so_D4_KEYS, so_D4_KEYS)
    @test all(g -> all(h -> all(k -> _compose_keys(_compose_keys(g, h), k) == _compose_keys(g, _compose_keys(h, k)),
                               so_D4_KEYS), so_D4_KEYS), so_D4_KEYS)
    # mirrors square to the identity, the quarter turns to the half turn
    @test all(g -> _compose_keys(g, g) == _IDENTITY_KEY, filter(k -> k[1], so_D4_KEYS))
    @test _compose_keys(so_rot(1 // 4), so_rot(1 // 4)) == so_rot(1 // 2)
    @test _compose_keys(so_rot(3 // 4), so_rot(3 // 4)) == so_rot(1 // 2)
end

@testset "_compose_keys (matrix cross-check)" begin
    # the key of a product is the key of the product of the 2x2 matrices, for every pair of registered D4 elements
    elements = collect(D4_symmetry())
    for g in elements, h in elements
        key = _compose_keys(symmetry_key(g), symmetry_key(h))
        M = so_matrix(g) * so_matrix(h)
        idx = findfirst(s -> symmetry_key(s) == key, elements)
        if isnothing(idx)
            # only the product of an element with its inverse leaves the registry: the identity
            @test key == _IDENTITY_KEY
            @test M ≈ [1.0 0.0; 0.0 1.0] atol = 1e-12
        else
            @test M ≈ so_matrix(elements[idx]) atol = 1e-12
        end
    end
end

@testset "_keys_closure" begin
    kx, ky, kxy = symmetry_key(XAxisReflection()), symmetry_key(YAxisReflection()), symmetry_key(XYAxisReflection())
    # no generators: the identity alone
    @test _keys_closure(typeof(kx)[]) == [_IDENTITY_KEY]
    # X and Y generate D2; the identity comes first
    d2 = _keys_closure([kx, ky])
    @test d2 isa Vector
    @test length(d2) == 4
    @test d2[1] == _IDENTITY_KEY
    @test Set(d2) == Set([_IDENTITY_KEY, kx, ky, kxy])
    # one reflection has order 2, the pi rotation too
    @test Set(_keys_closure([kx])) == Set([_IDENTITY_KEY, kx])
    @test length(_keys_closure([kxy])) == 2
    # Y and the diagonal generate D4: eight elements, exactly the registered ones plus the identity
    d4 = _keys_closure([ky, symmetry_key(DiagonalReflection())])
    @test length(d4) == 8
    @test Set(d4) == Set(so_D4_KEYS)
    # one rotation by 2 pi/6 generates the cyclic group of order 6
    c6 = _keys_closure([so_rot(1 // 6)])
    @test length(c6) == 6
    @test Set(c6) == Set(so_rot(k // 6) for k in 0:5)
    # redundant generators do not enlarge the closure
    @test Set(_keys_closure([kx, ky, kxy])) == Set(d2)
    # the closure is closed under composition
    @test all(g -> all(h -> _compose_keys(g, h) in d4, d4), d4)
end

@testset "_block_keys" begin
    circ = CircleBilliard(1.0)
    # identity followed by the registry (Y, XY, X)
    @test _block_keys(circ) == [_IDENTITY_KEY, symmetry_key(YAxisReflection()), symmetry_key(XYAxisReflection()),
                                symmetry_key(XAxisReflection())]
    # a registry gives the same keys as the billiard that holds it
    @test _block_keys(circ.symmetries) == _block_keys(circ)
    # a billiard with no symmetry has only the identity block
    @test _block_keys(TriangleBilliard(1.0, 1.0)) == [_IDENTITY_KEY]
    # cyclic billiard: rotations by 0, 2 pi/3, 4 pi/3
    @test _block_keys(C3Billiard(0.2)) == [so_rot(0 // 1), so_rot(1 // 3), so_rot(2 // 3)]
    # D4: identity then the seven registered elements in order
    @test _block_keys(so_d4_billiard()) == so_D4_KEYS

    # mixed reflection/rotation registries and registries mixing rotation orders are accepted
    @test length(_block_keys(so_polar_billiard(XAxisReflection(), NFoldRotation(3, 1)))) == 3
    @test _block_keys(so_polar_billiard(NFoldRotation(4, 1), NFoldRotation(2, 1), NFoldRotation(4, 3))) ==
          [so_rot(0 // 1), so_rot(1 // 4), so_rot(1 // 2), so_rot(3 // 4)]

    # only duplicate geometric actions are rejected
    @test_throws ArgumentError _block_keys(so_polar_billiard(XAxisReflection(), XAxisReflection()))
    @test_throws ArgumentError _block_keys(so_polar_billiard(NFoldRotation(3, 1), NFoldRotation(3, 1)))
    # the same action under two names: rotation by pi as XY and as NFoldRotation(2, 1), the diagonal as AxisReflection{5//4}
    @test_throws ArgumentError _block_keys(so_polar_billiard(XYAxisReflection(), NFoldRotation(2, 1)))
    @test_throws ArgumentError _block_keys(so_polar_billiard(DiagonalReflection(), AxisReflection{5 // 4}()))
    # a composite reflection is not a reconstruction-sector element (it cannot be registered, so build the registry directly)
    comp_reg = SymmetryRegistry(AbsSymmetry[CompositeReflection(XAxisReflection(), YAxisReflection())])
    @test_throws ArgumentError _block_keys(comp_reg)
end

@testset "_block_permutation" begin
    circ = CircleBilliard(1.0)
    keys = _block_keys(circ)
    # left multiplication permutes the blocks (id, Y, XY, X)
    @test _block_permutation(keys, XAxisReflection()) == [4, 3, 2, 1]
    @test _block_permutation(keys, YAxisReflection()) == [2, 1, 4, 3]
    @test _block_permutation(keys, XYAxisReflection()) == [3, 4, 1, 2]
    # the permutation of XY is the product of those of X and Y
    @test _compose_index_permutations(_block_permutation(keys, XAxisReflection()), _block_permutation(keys, YAxisReflection())) ==
          _block_permutation(keys, XYAxisReflection())
    # an element outside the registered group is rejected
    @test_throws ArgumentError _block_permutation(keys, DiagonalReflection())
    @test_throws ArgumentError _block_permutation(_block_keys(TriangleBilliard(1.0, 1.0)), XAxisReflection())
    @test_throws ArgumentError _block_permutation(keys, NFoldRotation(3, 1))

    # rotations on C3
    c3keys = _block_keys(C3Billiard(0.2))
    @test _block_permutation(c3keys, NFoldRotation(3, 1)) == [2, 3, 1]
    @test _block_permutation(c3keys, NFoldRotation(3, 2)) == [3, 1, 2]
    # a rotation of another order is not in the registered group
    @test_throws ArgumentError _block_permutation(c3keys, NFoldRotation(4, 1))

    # D4: g maps the identity block to its own block, and every row is a permutation of 1:8
    d4keys = _block_keys(so_d4_billiard())
    for g in D4_symmetry()
        p = _block_permutation(d4keys, g)
        @test sort(p) == 1:8
        @test p[1] == g.sym_id + 1
    end
    # the identity rotation acts trivially
    @test _block_permutation(d4keys, NFoldRotation(4, 0)) == collect(1:8)
end

@testset "_index_combinations" begin
    @test _index_combinations(3, 2) == [[1, 2], [1, 3], [2, 3]]
    @test _index_combinations(3, 1) == [[1], [2], [3]]
    @test _index_combinations(3, 3) == [[1, 2, 3]]
    @test _index_combinations(4, 2) == [[1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]]
    # more elements requested than available: no subsets
    @test isempty(_index_combinations(2, 3))
    # C(5, 3) = 10 strictly increasing subsets, distinct and in lexicographic order
    combos = _index_combinations(5, 3)
    @test length(combos) == 10
    @test all(c -> issorted(c; lt = <=) && all(in(1:5), c), combos)
    @test issorted(combos)
    @test allunique(combos)
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
    # the block-level helper gives the same tables
    @test _block_node_permutation(so_N, _block_keys(circ), XAxisReflection()) == so_RX
    @test _block_node_permutation(so_N3, _block_keys(c3), NFoldRotation(3, 1)) == [3, 4, 5, 6, 1, 2]
    # elements outside the registered group are rejected on every billiard
    @test_throws ArgumentError _boundary_symmetry_permutation(circ, so_N, DiagonalReflection())
    @test_throws ArgumentError _boundary_symmetry_permutation(TriangleBilliard(1.0, 1.0), so_N, XAxisReflection())
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
    circ_d = so_polar_billiard(DiagonalReflection())
    orbits = symmetry_index_orbits(Float64, circ_d, so_N, DiagonalReflection())
    @test fundamental_size(orbits) == so_N ÷ 2
    @test length(orbits) == so_N
    @test orbit_size(orbits) == 2
    @test orbits.fund_to_full == [1 2 3 4; 8 7 6 5]
    @test all(==(2), so_member_counts(orbits))
    # the diagonal reflection needs a multiple of 8 nodes
    @test_throws ArgumentError symmetry_index_orbits(Float64, circ_d, 12, DiagonalReflection())
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

@testset "symmetry_index_orbits (AxisReflection)" begin
    # any registered special axis is accepted, on a registry that holds the whole D4 group
    d4 = so_d4_billiard()
    for sym in (XAxisReflection(), YAxisReflection(), DiagonalReflection(), AntiDiagonalReflection())
        orbits = symmetry_index_orbits(Float64, d4, 16, sym)
        @test orbit_size(orbits) == 2
        @test fundamental_size(orbits) == 8
        @test sort(vec(orbits.fund_to_full)) == 1:16        # the orbits partition the nodes
        @test orbits.fund_to_full[2, :] == [_boundary_symmetry_permutation(d4, 16, sym)[q] for q in orbits.fundamental_indices]
    end
    # a mirror at a generic angle has no node multiple, even when registered
    generic = so_polar_billiard(AxisReflection{1 // 6}())
    @test_throws ArgumentError symmetry_node_multiple(AxisReflection{1 // 6}())
    @test_throws ArgumentError symmetry_index_orbits(Float64, generic, 12, AxisReflection{1 // 6}())
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

@testset "_orbit_map_from_generators" begin
    circ = CircleBilliard(1.0)
    keys = _block_keys(circ)
    # generators X and Y with characters (-1, 1) reproduce the native XY orbit map: images in the order id, X, Y, XY
    orbits = _orbit_map_from_generators(Float64, so_N, keys, AbsSymmetry[XAxisReflection(), YAxisReflection()], ComplexF64[-1, 1])
    ref = symmetry_index_orbits(Float64, circ, so_N, XYAxisReflection(), ComplexF64(-1), ComplexF64(1))
    @test orbits isa SymmetryOrbitMap{Float64}
    @test orbits.fund_to_full == ref.fund_to_full == [1 2; 8 7; 4 3; 5 6]
    @test orbits.fund_to_scale == ref.fund_to_scale
    @test orbits.orbit_of == ref.orbit_of
    # a single generator gives the order-2 group
    single = _orbit_map_from_generators(Float64, so_N, keys, [XAxisReflection()], ComplexF64[1])
    @test orbit_size(single) == 2
    @test single.fund_to_full == [1 2 3 4; 8 7 6 5]
    # a repeated generator with the same character is consistent
    @test orbit_size(_orbit_map_from_generators(Float64, so_N, keys, [XAxisReflection(), XAxisReflection()], ComplexF64[1, 1])) == 2

    # one character per generator
    @test_throws DimensionMismatch _orbit_map_from_generators(Float64, so_N, keys, [XAxisReflection(), YAxisReflection()], ComplexF64[1])
    # the same group element reached with two different characters is inconsistent
    @test_throws ArgumentError _orbit_map_from_generators(Float64, so_N, keys, [XAxisReflection(), XAxisReflection()], ComplexF64[1, -1])
    # a generator outside the registered group is rejected
    @test_throws ArgumentError _orbit_map_from_generators(Float64, so_N, keys, [DiagonalReflection()], ComplexF64[1])
end

@testset "symmetry_index_orbits (CompositeReflection, D4)" begin
    d4 = so_d4_billiard()
    comp = CompositeReflection(YAxisReflection(), DiagonalReflection())
    # Y and the diagonal generate D4; the character is +1 for Y and -1 for the diagonal
    orbits = symmetry_index_orbits(Float64, d4, 16, comp, ComplexF64[1, -1])
    @test orbit_size(orbits) == 8
    @test fundamental_size(orbits) == 2
    @test length(orbits) == 16
    # observed values on the polar D4 registry (shift -pi/4), members listed in the closure order of the generators
    # (re-observed after D4_symmetry() moved to chain order Diag, R1, Y, XY, Anti, R3, X)
    @test orbits.fund_to_full[:, 1] == [1, 8, 4, 13, 5, 12, 16, 9]
    @test orbits.fund_to_scale[:, 1] == ComplexF64[1, 1, -1, -1, -1, -1, 1, 1]
    # the orbits cover every node exactly once
    @test sort(vec(orbits.fund_to_full)) == 1:16
    @test all(==(8), so_member_counts(orbits))
    # the first image is the identity, the second image is the first generator (Y), the third the second generator
    @test orbits.fund_to_full[2, :] == [_boundary_symmetry_permutation(d4, 16, YAxisReflection())[q] for q in orbits.fundamental_indices]
    @test orbits.fund_to_full[3, :] == [_boundary_symmetry_permutation(d4, 16, DiagonalReflection())[q] for q in orbits.fundamental_indices]
    @test orbits.fund_to_scale[2, :] == ComplexF64[1, 1]
    @test orbits.fund_to_scale[3, :] == ComplexF64[-1, -1]

    # trivial characters give a trivial irrep factor on every image
    trivial = symmetry_index_orbits(Float64, d4, 16, comp)
    @test all(==(one(ComplexF64)), trivial.fund_to_scale)
    @test trivial.fund_to_full == orbits.fund_to_full

    # the node count must be a multiple of lcm(4, 8) = 8
    @test_throws ArgumentError symmetry_index_orbits(Float64, d4, 12, comp)
    # unsupported entries of the composite
    @test_throws ArgumentError symmetry_index_orbits(Float64, d4, 16, CompositeReflection(AxisReflection{1 // 6}()))
    @test_throws ArgumentError symmetry_index_orbits(Float64, d4, 16, CompositeReflection(CompositeReflection(XAxisReflection())))
end

@testset "symmetry_index_orbits (NFoldRotation subgroup orbits)" begin
    c6 = so_polar_billiard((NFoldRotation(6, m) for m in 1:5)...)
    # the rotation by 2 pi/6 shifts the nodes by N/6 = 2 and generates the whole group
    full = symmetry_index_orbits(Float64, c6, 12, NFoldRotation(6, 1))
    @test orbit_size(full) == 6
    @test fundamental_size(full) == 2
    @test full.fund_to_full[:, 1] == [1, 3, 5, 7, 9, 11]
    # the rotation by 2 pi/3 (m = 2) generates the subgroup of order 6 / gcd(6, 2) = 3
    sub = symmetry_index_orbits(Float64, c6, 12, NFoldRotation(6, 2))
    @test orbit_size(sub) == 3
    @test fundamental_size(sub) == 4
    @test sub.fund_to_full[:, 1] == [1, 5, 9]
    @test sort(vec(sub.fund_to_full)) == 1:12
    @test all(==(3), so_member_counts(sub))
    @test all(==(one(ComplexF64)), sub.fund_to_scale)       # sector 0: trivial representation
    # the rotation by pi (m = 3) generates the subgroup of order 2
    half = symmetry_index_orbits(Float64, c6, 12, NFoldRotation(6, 3))
    @test orbit_size(half) == 2
    @test half.fund_to_full[:, 1] == [1, 7]

    # sector l carries the phase exp(2 pi i l j / d) on the image j of the subgroup of order d
    sector = symmetry_index_orbits(Float64, c6, 12, NFoldRotation(6, 2), 1)
    ω = cis(2pi / 3)
    @test sector.fund_to_scale[:, 1] ≈ ComplexF64[1, ω, ω^2]
    @test sector.fund_to_full == sub.fund_to_full
    sector2 = symmetry_index_orbits(Float64, c6, 12, NFoldRotation(6, 1), 2)
    @test sector2.fund_to_scale[:, 1] ≈ [cis(2pi * 2 * j / 6) for j in 0:5]

    # the node count must be a multiple of the rotation order
    @test_throws ArgumentError symmetry_index_orbits(Float64, c6, 9, NFoldRotation(6, 2))

    # on C3 the rotation by 4 pi/3 is its own generator: the images of node 1 are 1, 5, 3 (not the m = 1 order 1, 3, 5)
    c3 = C3Billiard(0.2)
    r2 = symmetry_index_orbits(Float64, c3, so_N3, NFoldRotation(3, 2))
    @test orbit_size(r2) == 3
    @test r2.fund_to_full[:, 1] == [1, 5, 3]
    @test sort(vec(r2.fund_to_full)) == 1:so_N3

    # registries mixing rotation orders are accepted: C4 registered as rotations of order 4, 2 and 4
    c4 = so_polar_billiard(NFoldRotation(4, 1), NFoldRotation(2, 1), NFoldRotation(4, 3))
    pi_rot = symmetry_index_orbits(Float64, c4, 8, NFoldRotation(2, 1))
    @test orbit_size(pi_rot) == 2
    @test pi_rot.fund_to_full[:, 1] == [1, 5]
    @test orbit_size(symmetry_index_orbits(Float64, c4, 8, NFoldRotation(4, 1))) == 4
    # and on D4 the quarter turn acts together with the mirrors
    d4 = so_d4_billiard()
    quarter = symmetry_index_orbits(Float64, d4, 16, NFoldRotation(4, 1), 1)
    @test orbit_size(quarter) == 4
    @test quarter.fund_to_scale[:, 1] ≈ ComplexF64[1, im, -1, -im]
    @test sort(vec(quarter.fund_to_full)) == 1:16
end

@testset "symmetry registry: mixed registries accepted, duplicates rejected" begin
    # reflections and rotations in one registry, here the D4 group: accepted, orbit maps are available for every element
    d4 = so_d4_billiard()
    @test length(_block_keys(d4)) == 8
    @test fundamental_size(symmetry_index_orbits(Float64, d4, 16, XYAxisReflection())) == 4
    @test fundamental_size(symmetry_index_orbits(Float64, d4, 16, NFoldRotation(4, 3))) == 4
    # an element whose block does not exist in a non-closed mixed registry is rejected when it acts
    mixed = so_polar_billiard(XAxisReflection(), NFoldRotation(3, 1))
    @test_throws ArgumentError symmetry_index_orbits(Float64, mixed, 12, NFoldRotation(3, 1))
    # duplicate actions are rejected, whatever the order they were registered in
    dup = so_polar_billiard(XAxisReflection(), XAxisReflection())
    @test_throws ArgumentError symmetry_index_orbits(Float64, dup, so_N, XAxisReflection())
    @test_throws ArgumentError _boundary_symmetry_permutation(dup, so_N, XAxisReflection())
end

@testset "_minimal_generators" begin
    # D2 registry (Y, XY, X): the two mirror lines generate the group, XY is not needed
    gens = _minimal_generators(CircleBilliard(1.0).symmetries)
    @test gens isa Tuple
    @test typeof.(gens) == (YAxisReflection, XAxisReflection)
    @test [g.sym_id for g in gens] == [1, 3]
    @test typeof.(_minimal_generators(D2_symmetry())) == typeof.(gens)
    @test [g.sym_id for g in _minimal_generators(D2_symmetry())] == [1, 3]

    # D4 registry (chain order Diag, R1, Y, ...): the diagonal mirror and the y axis mirror generate all eight elements
    d4gens = _minimal_generators(D4_symmetry())
    @test typeof.(d4gens) == (DiagonalReflection, YAxisReflection)
    @test [g.sym_id for g in d4gens] == [1, 3]
    @test length(_keys_closure([symmetry_key(g) for g in d4gens])) == 8

    # cyclic registry: the rotation m = 1 alone, a rotation is only used when no mirror is registered
    c6gens = _minimal_generators(Cn_symmetry(6))
    @test length(c6gens) == 1
    @test only(c6gens) isa NFoldRotation
    @test only(c6gens).m == 1
    # C4 registered with mixed rotation orders: the first rotation generates
    c4 = so_polar_billiard(NFoldRotation(4, 1), NFoldRotation(2, 1), NFoldRotation(4, 3))
    @test _minimal_generators(c4.symmetries) isa Tuple{NFoldRotation}
    @test (only(_minimal_generators(c4.symmetries)).order, only(_minimal_generators(c4.symmetries)).m) == (4, 1)

    # a single element is its own generator, also when it is a pi rotation
    @test typeof.(_minimal_generators(so_polar_billiard(XYAxisReflection()).symmetries)) == (XYAxisReflection,)
    @test typeof.(_minimal_generators(so_polar_billiard(DiagonalReflection()).symmetries)) == (DiagonalReflection,)

    # no registered symmetry: empty generating set
    @test _minimal_generators(TriangleBilliard(1.0, 1.0).symmetries) == ()
    # two mirror lines without their product are not a group
    @test_throws ArgumentError _minimal_generators(so_polar_billiard(XAxisReflection(), YAxisReflection()).symmetries)
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
    @test get_symmetries(so_polar_billiard(DiagonalReflection())) isa Tuple{DiagonalReflection}

    # no registered symmetry: empty generating set
    @test get_symmetries(TriangleBilliard(1.0, 1.0)) == ()

    # D4 billiard: the same generators as the registry
    d4 = so_d4_billiard()
    @test typeof.(get_symmetries(d4)) == (DiagonalReflection, YAxisReflection)

    # a registry is accepted directly and gives the same generators as the billiard that holds it
    @test typeof.(get_symmetries(D2_symmetry())) == typeof.(get_symmetries(RectangleBilliard(1.0, 0.6)))
    @test [g.sym_id for g in get_symmetries(D2_symmetry())] == [1, 3]
    @test [g.sym_id for g in get_symmetries(D4_symmetry())] == [1, 3]
    @test get_symmetries(Cn_symmetry(1)) == ()
    c3 = C3Billiard(0.2)
    @test [g.m for g in get_symmetries(c3.symmetries)] == [g.m for g in get_symmetries(c3)] == [1]
end

@testset "get_symmetry_group" begin
    @test :get_symmetry_group in names(BilliardGeometry)
    circ = CircleBilliard(1.0)
    @test get_symmetry_group(circ) === circ.symmetries
    @test get_symmetry_group(circ) isa SymmetryRegistry
    @test length(get_symmetry_group(circ)) == 3
    @test length(get_symmetry_group(C3Billiard(0.2))) == 2
    @test length(get_symmetry_group(so_d4_billiard())) == 7
    # no registered symmetry: empty registry
    @test length(get_symmetry_group(TriangleBilliard(1.0, 1.0))) == 0
    # the group order is the registry length plus one
    @test symmetry_reduction_factor(circ) == length(get_symmetry_group(circ)) + 1
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
