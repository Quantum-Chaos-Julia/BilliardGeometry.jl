################################################################################
# DISCRETE BOUNDARY SYMMETRY ORBITS
#
# A BIM discretization samples each complete physical boundary component on a
# periodic midpoint grid. The action of a discrete symmetry on those nodes is
# determined exactly from the same registered symmetry group used by
# `full_boundary`; it is never inferred from floating-point coordinates.
#
# Boundary loops are grouped by the source that reconstructs them (outer
# boundary, one group per hole object). A group is split into one contiguous
# block of nodes per reconstruction block; an element of the symmetry group
# permutes the blocks by exact key multiplication, so loops exchanged by a
# symmetry (copies of a `SymmetrizedHoles`) are handled as blocks of one group.
################################################################################

"""
    SymmetryOrbitMap{T}

Exact folding of a periodic boundary discretization under a discrete symmetry
group.

## Description
`SymmetryOrbitMap` records the complete discrete symmetry orbit of every
fundamental boundary node together with the irreducible-representation factor
relating each image to its representative. The node permutations are derived
algebraically from the billiard's registered physical-boundary reconstruction;
no floating-point coordinate matching is used.

## Attributes
* `fundamental_indices::Vector{Int}`: Full-boundary indices chosen as orbit representatives.
* `orbit_of::Vector{Int}`: Fundamental-orbit index associated with each full-boundary node.
* `phase::Vector{Complex{T}}`: Irreducible-representation factor associated with each full-boundary node.
* `full_size::Int`: Number of nodes on the complete boundary component.
* `fundamental_size::Int`: Number of symmetry-reduced boundary nodes.
* `fund_to_full::Matrix{Int}`: Full node index of each group image of each representative.
* `fund_to_scale::Matrix{Complex{T}}`: Irrep factor of each group image.
"""
struct SymmetryOrbitMap{T<:Real}
    fundamental_indices::Vector{Int}
    orbit_of::Vector{Int}
    phase::Vector{Complex{T}}
    full_size::Int
    fundamental_size::Int
    fund_to_full::Matrix{Int}
    fund_to_scale::Matrix{Complex{T}}
end

"""
    fundamental_size(orbits::SymmetryOrbitMap) → n::Int

Return the number of symmetry-reduced boundary nodes.

## Arguments
* `orbits::SymmetryOrbitMap`: Boundary symmetry-orbit map.

## Returns
* `n::Int`: Number of fundamental boundary nodes.
"""
fundamental_size(orbits::SymmetryOrbitMap) = orbits.fundamental_size

"""
    length(orbits::SymmetryOrbitMap) → n::Int

Return the number of nodes on the complete boundary component.

## Arguments
* `orbits::SymmetryOrbitMap`: Boundary symmetry-orbit map.

## Returns
* `n::Int`: Number of complete-boundary nodes.
"""
Base.length(orbits::SymmetryOrbitMap) = orbits.full_size

"""
    full_size(orbits::SymmetryOrbitMap) → n::Int

Return the number of nodes on the complete boundary component.

## Arguments
* `orbits::SymmetryOrbitMap`: Boundary symmetry-orbit map.

## Returns
* `n::Int`: Number of complete-boundary nodes.
"""
full_size(orbits::SymmetryOrbitMap) = orbits.full_size

"""
    orbit_size(orbits::SymmetryOrbitMap) → ng::Int

Return the number of symmetry-group images in each boundary orbit.

## Arguments
* `orbits::SymmetryOrbitMap`: Boundary symmetry-orbit map.

## Returns
* `ng::Int`: Symmetry-group order.
"""
orbit_size(orbits::SymmetryOrbitMap) = size(orbits.fund_to_full, 1)

"""
    symmetry_orbit(orbits::SymmetryOrbitMap, b::Int) → (qs, χs)

Return the complete node indices and irrep factors in fundamental orbit `b`.

## Arguments
* `orbits::SymmetryOrbitMap`: Boundary symmetry-orbit map.
* `b::Int`: Fundamental-orbit index.

## Returns
* `qs`: View of complete-boundary node indices.
* `χs`: View of corresponding irreducible-representation factors.
"""
@inline function symmetry_orbit(orbits::SymmetryOrbitMap, b::Int)
    return @view(orbits.fund_to_full[:, b]), @view(orbits.fund_to_scale[:, b])
end

"""
    symmetry_node_multiple(symmetry::AbsSymmetry) → n::Int

Return the node-count multiple required by a boundary symmetry reduction.

## Arguments
* `symmetry::AbsSymmetry`: Active discrete symmetry.

## Returns
* `n::Int`: Required node-count multiple.
"""
symmetry_node_multiple(::XAxisReflection) = 4
symmetry_node_multiple(::YAxisReflection) = 4
symmetry_node_multiple(::XYAxisReflection) = 4
symmetry_node_multiple(::DiagonalReflection) = 8
symmetry_node_multiple(::AntiDiagonalReflection) = 8
symmetry_node_multiple(::AxisReflection{A}) where {A} = throw(ArgumentError("symmetry_node_multiple is undefined for AxisReflection{$A}; only the axes 0, 1/4, 1/2, 3/4 are supported."))
symmetry_node_multiple(sym::NFoldRotation) = sym.order
symmetry_node_multiple(sym::CompositeReflection) = foldl(lcm, (symmetry_node_multiple(ref) for ref in sym.reflections); init = 1)

################################################################################
# EXACT REGISTERED-GROUP ALGEBRA
################################################################################

# Dihedral key `(reflects, a)` of a symmetry element, `a` taken mod 1: `rot(a)` is the
# rotation by `2πa`, `refl(a)` the mirror through the origin at angle `πa`.
@inline symmetry_key(::AxisReflection{A}) where {A} = (true, mod(A, 1//1))
@inline symmetry_key(::XYAxisReflection) = (false, 1//2)
@inline symmetry_key(sym::NFoldRotation) = (false, mod(sym.m // sym.order, 1//1))

const _DihedralKey = Tuple{Bool,Rational{Int}}

# Key of the identity element.
const _IDENTITY_KEY = (false, 0 // 1)::_DihedralKey

# Key of g ∘ h. With g = (rg, a), h = (rh, b): the product reflects iff exactly
# one factor does, and the angle is a - b if g reflects, a + b otherwise.
@inline function _compose_keys(g::_DihedralKey, h::_DihedralKey)
    rg, a = g
    rh, b = h
    return (rg ⊻ rh, mod(rg ? a - b : a + b, 1 // 1))::_DihedralKey
end

# Grow the set of keys containing the identity under composition with the generator keys.
function _keys_closure(gen_keys::AbstractVector{<:_DihedralKey})
    closure = _DihedralKey[_IDENTITY_KEY]
    head = 1
    while head <= length(closure)
        a = closure[head]
        for g in gen_keys
            c = _compose_keys(a, g)
            c in closure || push!(closure, c)
        end
        head += 1
    end
    return closure
end

# Keys of the reconstruction blocks of a registry in the order used by
# `full_boundary`: identity followed by the registry.
function _block_keys(reg::SymmetryRegistry)
    keys = _DihedralKey[_IDENTITY_KEY]
    for sym in reg
        sym isa CompositeReflection && throw(ArgumentError("CompositeReflection cannot be a physical reconstruction-sector element"))
        push!(keys, symmetry_key(sym))
    end
    length(unique(keys)) == length(keys) || throw(ArgumentError("SymmetryRegistry contains duplicate geometric actions"))
    return keys
end
_block_keys(billiard::AbsBilliard) = _block_keys(billiard.symmetries)

# Exact permutation of the blocks under left multiplication by `g`.
function _block_permutation(block_keys::Vector{_DihedralKey}, g::AbsSymmetry)
    gk = symmetry_key(g)
    p = Vector{Int}(undef, length(block_keys))
    @inbounds for s in eachindex(block_keys)
        j = findfirst(==(_compose_keys(gk, block_keys[s])), block_keys)
        isnothing(j) && throw(ArgumentError("$(typeof(g)) is not contained in the billiard's registered physical symmetry group"))
        p[s] = j
    end
    return p
end

# All size-`k` index subsets of `1:n` in lexicographic order.
function _index_combinations(n::Int, k::Int)
    combos = Vector{Vector{Int}}()
    k > n && return combos
    c = collect(1:k)
    while true
        push!(combos, copy(c))
        i = k
        while i >= 1 && c[i] == n - k + i
            i -= 1
        end
        i == 0 && break
        c[i] += 1
        for j in i + 1:k
            c[j] = c[j - 1] + 1
        end
    end
    return combos
end

# Smallest set of registry elements whose closure is the whole registry
# group. Candidates are tried by increasing size, in registration order with
# mirror reflections before orientation-preserving elements, so that e.g. `D₂`'s
# `XYAxisReflection` (== `XAxisReflection` ∘ `YAxisReflection`) is only
# picked if the mirrors alone are not enough.
function _minimal_generators(reg::SymmetryRegistry)
    isempty(reg) && return ()
    elements = collect(reg)
    keys = _DihedralKey[symmetry_key(sym) for sym in elements]
    ngoal = length(_keys_closure(keys))
    ngoal == length(keys) + 1 || throw(ArgumentError("SymmetryRegistry is not closed under composition by its own elements"))
    order = vcat(findall(sym -> orientation_reversing(sym), elements),
                 findall(sym -> !orientation_reversing(sym), elements))
    for k in 1:length(order)
        for combo in _index_combinations(length(order), k)
            idx = order[combo]
            length(_keys_closure(keys[idx])) == ngoal && return Tuple(elements[idx])
        end
    end
    throw(ArgumentError("SymmetryRegistry is not closed under composition by its own elements"))
end

"""
    get_symmetries(billiard::AbsBilliard) → gens::Tuple{Vararg{AbsSymmetry}}

Return a minimal set of fundamental symmetry generators for `billiard`,
picked from its registered [`SymmetryRegistry`](@ref) (`billiard.symmetries`),
in registration order.

## Description
`billiard.symmetries` lists *every* non-identity element `full_boundary`
needs to reconstruct the complete physical boundary (e.g. a `D₂`-symmetric
billiard registered via [`D2_symmetry`](@ref) lists all three of
`YAxisReflection`, `XYAxisReflection`, `XAxisReflection`). For assigning
per-generator irreducible-representation characters (see `SymmetrySector` in
`QuantumBilliards.jl`), only a minimal *generating* subset is needed, e.g.
`(YAxisReflection(), XAxisReflection())` for `D₂`, since
`XYAxisReflection = XAxisReflection ∘ YAxisReflection` is already implied by
composing the other two.

Generators are found by exact dihedral-key composition (the algebra used by
[`symmetry_index_orbits`](@ref)): the smallest set of registered elements whose
closure is the whole group is returned, trying sets of increasing size in
registration order with mirror reflections before orientation-preserving
elements. For `D₂` this is `(YAxisReflection, XAxisReflection)`, for `D₄`
`(YAxisReflection, DiagonalReflection)`, and for a cyclic registry its
`sym_id=1` element.

## Arguments
* `billiard::AbsBilliard`: Billiard whose registered symmetry group is reduced to generators.

## Returns
* `gens::Tuple{Vararg{AbsSymmetry}}`: Minimal generating subset of `billiard.symmetries`.
"""
get_symmetries(billiard::Bi) where {Bi<:AbsBilliard} = _minimal_generators(get_symmetry_group(billiard))
get_symmetries(reg::SymmetryRegistry) = _minimal_generators(reg)

"""
    get_symmetry_group(billiard::AbsBilliard) → reg::SymmetryRegistry

Return the registry of the symmetry elements of the physical billiard, the
group whose irreps label the eigenstates.

## Description
For a billiard without holes this is `billiard.symmetries`. For a billiard with
holes it is the stabilizer `G' = {g ∈ G : g(holes) = holes}` of the registered
outer group `G`, the intersection over the holes of
[`compatible_symmetries`](@ref); `G' = {e}` (an empty registry) is the
unreduced case. The returned elements are those of `billiard.symmetries`
(same `sym_id`s, registry order); `billiard.symmetries` itself keeps
describing `G`, from which [`full_boundary`](@ref) reconstructs the outer
boundary.

## Arguments
* `billiard::AbsBilliard`: Billiard whose symmetry group is requested.

## Returns
* `reg`: Registry of the non-identity elements of the symmetry group.
"""
get_symmetry_group(billiard::Bi) where {Bi<:AbsBilliard} = _symmetry_group(billiard.fundamental_domain, billiard)

_symmetry_group(::AbsDomain, billiard::AbsBilliard) = billiard.symmetries

################################################################################
# EXACT MIDPOINT-GRID PERMUTATIONS
################################################################################

"""
    _sector_node_permutation(N::Int, sector_perm::Vector{Int}, reversing::Bool) → p::Vector{Int}

Lift an exact reconstruction-sector permutation to a periodic midpoint-node
permutation.

## Description
The complete invariant boundary component is divided into one equal
parametrization block per registered physical symmetry sector. An
orientation-preserving action preserves the local midpoint index; an
orientation-reversing action reverses it. All operations are exact integer
operations.

## Arguments
* `N::Int`: Number of nodes on the complete invariant boundary component.
* `sector_perm::Vector{Int}`: Exact permutation of reconstruction sectors.
* `reversing::Bool`: Whether the geometric action reverses boundary orientation.

## Returns
* `p::Vector{Int}`: Exact full-component node permutation.
"""
function _sector_node_permutation(N::Int, sector_perm::Vector{Int}, reversing::Bool)
    ng = length(sector_perm)
    N % ng == 0 || throw(ArgumentError("Node count N = $N must be divisible by physical symmetry-group order $ng"))
    M = N ÷ ng
    p = Vector{Int}(undef, N)
    @inbounds for q in 1:N
        j = q - 1
        s = j ÷ M
        u = j - s * M
        sp = sector_perm[s + 1] - 1
        up = reversing ? M - 1 - u : u
        p[q] = sp * M + up + 1
    end
    return p
end

# Exact midpoint-index permutation induced by one physical symmetry element on
# an invariant boundary component split into the blocks of `block_keys`.
function _block_node_permutation(N::Int, block_keys::Vector{_DihedralKey}, symmetry::AbsSymmetry)
    sectors = _block_permutation(block_keys, symmetry)
    return _sector_node_permutation(N, sectors, orientation_reversing(symmetry))
end

function _boundary_symmetry_permutation(billiard::Bi, N::Int, symmetry::AbsSymmetry) where {Bi<:AbsBilliard}
    return _block_node_permutation(N, _block_keys(billiard), symmetry)
end

_boundary_symmetry_permutation(::Bi, ::Int, ::CompositeReflection) where {Bi<:AbsBilliard} =
    throw(ArgumentError("CompositeReflection contains several geometric generators"))

# Compose two index permutations as a ∘ b, so the resulting image of q is
# a[b[q]].
@inline function _compose_index_permutations(a::Vector{Int}, b::Vector{Int})
    length(a) == length(b) || throw(DimensionMismatch("Cannot compose permutations of lengths $(length(a)) and $(length(b))"))
    return [a[b[q]] for q in eachindex(b)]
end

################################################################################
# ORBIT-MAP CONSTRUCTION
################################################################################

# Build the final orbit lookup tables from exact full-boundary permutations
# and their corresponding irreducible-representation factors.
function _build_symmetry_orbit_map(::Type{T}, N::Int, perms::Vector{Vector{Int}}, scales::Vector{Complex{T}} = ones(Complex{T}, length(perms))) where {T<:Real}
    ng = length(perms)
    N % ng == 0 || throw(ArgumentError("Node count N = $N must be divisible by symmetry-group order $ng"))
    length(scales) == ng || throw(DimensionMismatch("Received $ng permutations but $(length(scales)) irrep factors"))
    all(length(p) == N for p in perms) || throw(DimensionMismatch("Every symmetry permutation must have length N = $N"))
    nf = N ÷ ng
    fundamental_indices = Vector{Int}(undef, nf)
    orbit_of = Vector{Int}(undef, N)
    phase = Vector{Complex{T}}(undef, N)
    fund_to_full = Matrix{Int}(undef, ng, nf)
    fund_to_scale = Matrix{Complex{T}}(undef, ng, nf)
    seen = falses(N)
    b = 0
    @inbounds for q in 1:N
        seen[q] && continue
        b += 1
        b <= nf || throw(ArgumentError("Symmetry action does not produce free orbits of size $ng"))
        fundamental_indices[b] = q
        for g in 1:ng
            qi = perms[g][q]
            χ = scales[g]
            1 <= qi <= N || throw(ArgumentError("Invalid symmetry image index $qi for N = $N"))
            orbit_of[qi] = b
            phase[qi] = χ
            fund_to_full[g, b] = qi
            fund_to_scale[g, b] = χ
            seen[qi] = true
        end
    end
    b == nf || throw(ArgumentError("Symmetry action produced $b boundary orbits; expected $nf"))
    return SymmetryOrbitMap{T}(fundamental_indices, orbit_of, phase, N, nf, fund_to_full, fund_to_scale)
end

################################################################################
# PUBLIC SINGLE-COMPONENT API
################################################################################

# Orbit map of the group generated by `generators` acting on the `block_keys`
# blocks of a component with `N` nodes. The group is built by closure under
# composition of the exact node permutations; `characters` assigns one
# character to each generator, and a group element reached by two different
# words must receive the same character.
function _orbit_map_from_generators(::Type{T}, N::Int, block_keys::Vector{_DihedralKey}, generators::AbstractVector{<:AbsSymmetry}, characters::AbstractVector{Complex{T}}) where {T<:Real}
    length(generators) == length(characters) || throw(DimensionMismatch("Received $(length(generators)) generators but $(length(characters)) characters"))
    genperms = Vector{Int}[_block_node_permutation(N, block_keys, g) for g in generators]
    perms = Vector{Vector{Int}}([collect(1:N)])
    scales = Complex{T}[one(Complex{T})]
    head = 1
    while head <= length(perms)
        p = perms[head]
        χ = scales[head]
        for g in eachindex(genperms)
            pref = _compose_index_permutations(genperms[g], p)
            χnew = characters[g] * χ
            j = findfirst(==(pref), perms)
            if isnothing(j)
                push!(perms, pref)
                push!(scales, χnew)
            elseif scales[j] != χnew
                throw(ArgumentError("Characters of the generators are inconsistent: one group element is reached by different words with different character products"))
            end
        end
        head += 1
    end
    return _build_symmetry_orbit_map(T, N, perms, scales)
end


# Orbit maps on a component split into the blocks of `block_keys`. Methods
# dispatch on the symmetry type; the node count must be a multiple of the
# number of blocks (checked by the block permutation).
function _symmetry_orbits(::Type{T}, block_keys::Vector{_DihedralKey}, N::Int, symmetry::AxisReflection, character::Complex{T} = one(Complex{T})) where {T<:Real}
    refl = _block_node_permutation(N, block_keys, symmetry)
    return _build_symmetry_orbit_map(T, N, Vector{Vector{Int}}([collect(1:N), refl]), Complex{T}[one(Complex{T}), character])
end

function _symmetry_orbits(::Type{T}, block_keys::Vector{_DihedralKey}, N::Int, ::XYAxisReflection, character_x::Complex{T} = one(Complex{T}), character_y::Complex{T} = one(Complex{T})) where {T<:Real}
    rx = _block_node_permutation(N, block_keys, XAxisReflection())
    ry = _block_node_permutation(N, block_keys, YAxisReflection())
    rxy = _compose_index_permutations(rx, ry)
    χxy = character_x * character_y
    return _build_symmetry_orbit_map(T, N, Vector{Vector{Int}}([collect(1:N), rx, ry, rxy]), Complex{T}[one(Complex{T}), character_x, character_y, χxy])
end

function _symmetry_orbits(::Type{T}, block_keys::Vector{_DihedralKey}, N::Int, symmetry::NFoldRotation, sector::Int = 0) where {T<:Real}
    d = symmetry.order ÷ gcd(symmetry.order, symmetry.m)
    g = _block_node_permutation(N, block_keys, symmetry)
    perms = Vector{Vector{Int}}(undef, d)
    perms[1] = collect(1:N)
    @inbounds for l in 2:d
        perms[l] = _compose_index_permutations(g, perms[l - 1])
    end
    scales = Complex{T}[cis(T(2 * pi) * T(sector * l) / T(d)) for l in 0:d - 1]
    return _build_symmetry_orbit_map(T, N, perms, scales)
end

function _symmetry_orbits(::Type{T}, block_keys::Vector{_DihedralKey}, N::Int, symmetry::CompositeReflection, characters::Vector{Complex{T}} = ones(Complex{T}, length(symmetry.reflections))) where {T<:Real}
    isempty(symmetry.reflections) && throw(ArgumentError("CompositeReflection requires at least one reflection"))
    length(characters) == length(symmetry.reflections) || throw(DimensionMismatch("Received $(length(symmetry.reflections)) reflections but $(length(characters)) characters"))
    generators = AbsSymmetry[]
    genscales = Complex{T}[]
    for (k, ref) in enumerate(symmetry.reflections)
        χ = characters[k]
        if ref isa XYAxisReflection
            push!(generators, XAxisReflection(), YAxisReflection())
            push!(genscales, χ, χ)
        elseif ref isa AxisReflection
            push!(generators, ref)
            push!(genscales, χ)
        else
            throw(ArgumentError("Unsupported reflection type $(typeof(ref)) in CompositeReflection"))
        end
    end
    return _orbit_map_from_generators(T, N, block_keys, generators, genscales)
end

"""
    symmetry_index_orbits(::Type{T}, billiard::Bi, N::Int, symmetry::AbsSymmetry, character...) where {T<:Real,Bi<:AbsBilliard} → orbits::SymmetryOrbitMap{T}

Construct the exact symmetry-orbit map for the complete outer physical
boundary component.

## Description
The node permutation is derived from the billiard's registered physical
symmetry reconstruction and the integer midpoint-grid index. Boundary
coordinates are never inspected and no floating-point matching is performed.
For a billiard with holes use [`group_symmetry_orbits`](@ref), which also
addresses the hole groups and checks the symmetry against the group `G'`.

## Arguments
* `T::Type{<:Real}`: Real scalar type of the boundary discretization.
* `billiard::Bi`: Billiard defining the registered physical symmetry group.
* `N::Int`: Number of nodes on the complete outer physical component.
* `symmetry::AbsSymmetry`: Active symmetry reduction.
* `character`: Requested irreducible-representation character or rotational sector.

## Returns
* `orbits::SymmetryOrbitMap{T}`: Exact boundary symmetry-orbit map.
"""
function symmetry_index_orbits(::Type{T}, billiard::Bi, N::Int, symmetry::AbsSymmetry, character...) where {T<:Real,Bi<:AbsBilliard}
    block_keys = _block_keys(billiard)
    multiple = lcm(_orbit_node_multiple(symmetry), length(block_keys))
    N % multiple == 0 || throw(ArgumentError("$(typeof(symmetry)) requires N divisible by $multiple; received N = $N"))
    return _symmetry_orbits(T, block_keys, N, symmetry, character...)
end

# Node multiple used by `symmetry_index_orbits`: reflections about axes other than the four special ones need no
# multiple of their own, the number of reconstruction blocks is enforced separately.
_orbit_node_multiple(symmetry::AbsSymmetry) = symmetry_node_multiple(symmetry)
_orbit_node_multiple(symmetry::AxisReflection{A}) where {A} = A in (0 // 1, 1 // 4, 1 // 2, 3 // 4) ? symmetry_node_multiple(symmetry) : 1
_orbit_node_multiple(symmetry::CompositeReflection) = foldl(lcm, (_orbit_node_multiple(ref) for ref in symmetry.reflections); init = 1)

################################################################################
# GROUP ORBITS (BILLIARDS WITH HOLES)
################################################################################

# Symmetry elements whose characters the orbit map of `symmetry` assigns.
_acting_generators(symmetry::AbsSymmetry) = AbsSymmetry[symmetry]
_acting_generators(::XYAxisReflection) = AbsSymmetry[XAxisReflection(), YAxisReflection()]
_acting_generators(symmetry::CompositeReflection) = reduce(vcat, (_acting_generators(ref) for ref in symmetry.reflections); init = AbsSymmetry[])

# Block keys of a solver group: outer boundary `0`, `i`-th hole object `i`.
function _group_block_keys(::AbsDomain, billiard::AbsBilliard, group::Int)
    group == 0 || throw(ArgumentError("billiard has no holes; the only boundary group is 0, received $group"))
    return _block_keys(billiard)
end

function _group_block_keys(domain::AbsMultiplyConnectedDomain, billiard::AbsBilliard, group::Int)
    group == 0 && return _block_keys(billiard)
    1 <= group <= length(domain.holes) || throw(ArgumentError("boundary group must be in 0:$(length(domain.holes)); received $group"))
    return _hole_block_keys(domain.holes[group])
end

# Index of the first hole not invariant under `g` although `g` belongs to the registry, else `nothing`.
_blocking_hole(::AbsDomain, ::SymmetryRegistry, ::AbsSymmetry) = nothing

function _blocking_hole(domain::AbsMultiplyConnectedDomain, registry::SymmetryRegistry, g::AbsSymmetry)
    k = symmetry_key(g)
    any(s -> symmetry_key(s) == k, registry) || return nothing
    for (i, h) in enumerate(domain.holes)
        any(s -> symmetry_key(s) == k, compatible_symmetries(h, registry)) || return i
    end
    return nothing
end

function _require_symmetry_in_group(billiard::AbsBilliard, symmetry::AbsSymmetry)
    group = get_symmetry_group(billiard)
    group_keys = _DihedralKey[symmetry_key(g) for g in group]
    for g in _acting_generators(symmetry)
        symmetry_key(g) in group_keys && continue
        description = isempty(group) ? "{e}" : "{e, " * join(string.(group), ", ") * "}"
        msg = "$g is not in the symmetry group G' = $description of the billiard"
        hole = _blocking_hole(billiard.fundamental_domain, billiard.symmetries, g)
        isnothing(hole) || (msg *= "; hole $hole is not invariant under it")
        throw(ArgumentError(msg * ". Use symmetry = nothing or a subgroup of G'."))
    end
    return nothing
end

"""
    group_symmetry_orbits(::Type{T}, billiard::Bi, group::Int, N::Int, symmetry::AbsSymmetry, character...) where {T<:Real,Bi<:AbsBilliard} → orbits::SymmetryOrbitMap{T}

Construct the exact symmetry-orbit map of one boundary group of a billiard,
acted on by elements of the symmetry group `G'`.

## Description
A group is the set of loops reconstructed from one source: `group = 0` is the
outer boundary (blocks of the registry `G`), `group = i` the `i`-th hole object
(blocks of the hole group, the placing operations of a [`SymmetrizedHoles`](@ref),
or the identity alone for a [`DomainHole`](@ref)). The nodes of a group are the
concatenation of its loops, one contiguous block of `N / n_blocks` nodes per
block. An element of `G'` permutes the blocks by left multiplication of its
dihedral key and reverses the index inside a block when it reverses
orientation, so every orbit has `|G'|` nodes.

## Arguments
* `T::Type{<:Real}`: Real scalar type of the boundary discretization.
* `billiard::Bi`: Billiard defining the symmetry groups.
* `group::Int`: Boundary group, `0` for the outer boundary, `i` for the `i`-th hole object.
* `N::Int`: Number of nodes of the group; a multiple of its number of blocks.
* `symmetry::AbsSymmetry`: Acting symmetry reduction; every generator must belong to `G'` (see [`get_symmetry_group`](@ref)).
* `character`: Requested irreducible-representation character or rotational sector.

## Returns
* `orbits::SymmetryOrbitMap{T}`: Exact boundary symmetry-orbit map of the group.

Throws `ArgumentError` naming `G'` (and the first hole that breaks the
symmetry) when a generator is outside `G'`.
"""
function group_symmetry_orbits(::Type{T}, billiard::Bi, group::Int, N::Int, symmetry::AbsSymmetry, character...) where {T<:Real,Bi<:AbsBilliard}
    _require_symmetry_in_group(billiard, symmetry)
    return _symmetry_orbits(T, _group_block_keys(billiard.fundamental_domain, billiard, group), N, symmetry, character...)
end
