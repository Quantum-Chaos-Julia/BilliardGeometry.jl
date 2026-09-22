################################################################################
# DISCRETE BOUNDARY SYMMETRY ORBITS
#
# A BIM discretization samples each complete physical boundary component on a
# periodic midpoint grid. The action of a discrete symmetry on those nodes is
# determined exactly from the same registered symmetry group used by
# `full_boundary`; it is never inferred from floating-point coordinates.
#
# Current restriction:
#
#   Every physical boundary component must be individually invariant under the
#   active symmetry. Symmetries that exchange distinct `domain_id` components
#   are not supported.
#
# Thus the same construction applies independently to the outer boundary and,
# for multiply connected billiards, to every invariant hole boundary.
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
symmetry_node_multiple(sym::NFoldRotation) = sym.order
symmetry_node_multiple(sym::CompositeReflection) = foldl(lcm, (symmetry_node_multiple(ref) for ref in sym.reflections); init = 1)

################################################################################
# EXACT REGISTERED-GROUP ALGEBRA
################################################################################

# Exact integer matrix key for the identity transformation.
const _SYMMETRY_IDENTITY_KEY = (1, 0, 0, 1)

# Return the exact integer matrix key of each implemented reflection.
@inline _symmetry_key(::XAxisReflection) = (1, 0, 0, -1)
@inline _symmetry_key(::YAxisReflection) = (-1, 0, 0, 1)
@inline _symmetry_key(::XYAxisReflection) = (-1, 0, 0, -1)
@inline _symmetry_key(::DiagonalReflection) = (0, 1, 1, 0)
@inline _symmetry_key(::AntiDiagonalReflection) = (0, -1, -1, 0)

# Compose two exact 2×2 integer symmetry matrices, returning the key of A ∘ B.
@inline function _compose_symmetry_keys(A::NTuple{4, Int}, B::NTuple{4, Int})
    a, b, c, d = A
    e, f, g, h = B
    return (a * e + b * g, a * f + b * h, c * e + d * g, c * f + d * h)
end

# Represent an N-fold rotation exactly by its order and reduced rotation power.
@inline _rotation_key(sym::NFoldRotation) = (sym.order, mod(sym.m, sym.order))

# Return the exact identity key for an N-fold rotation group.
@inline _rotation_identity_key(n::Int) = (n, 0)

# Compose two powers of the same cyclic rotation group.
@inline function _compose_rotation_keys(a::Tuple{Int, Int}, b::Tuple{Int, Int})
    a[1] == b[1] || throw(ArgumentError("Cannot compose rotations of orders $(a[1]) and $(b[1])"))
    return (a[1], mod(a[2] + b[2], a[1]))
end

# Return the exact geometric keys of the physical reconstruction sectors in
# the same order used by `full_boundary`: identity followed by the registry.
function _reflection_sector_keys(billiard::Bi) where {Bi<:AbsBilliard}
    keys = NTuple{4, Int}[_SYMMETRY_IDENTITY_KEY]
    @inbounds for sym in billiard.symmetries
        sym isa AbsReflection || throw(ArgumentError("Expected a reflection-group SymmetryRegistry; found $(typeof(sym))"))
        sym isa CompositeReflection && throw(ArgumentError("CompositeReflection cannot be a physical reconstruction-sector element"))
        push!(keys, _symmetry_key(sym))
    end
    length(unique(keys)) == length(keys) || throw(ArgumentError("SymmetryRegistry contains duplicate geometric actions"))
    return keys
end

# Compute how one reflection permutes the physical reconstruction sectors by
# exact group multiplication, without inspecting any boundary coordinates.
function _reflection_sector_permutation(billiard::Bi, action::AbsReflection) where {Bi<:AbsBilliard}
    action isa CompositeReflection && throw(ArgumentError("CompositeReflection is a generator set, not one geometric action"))
    keys = _reflection_sector_keys(billiard)
    g = _symmetry_key(action)
    p = Vector{Int}(undef, length(keys))
    @inbounds for s in eachindex(keys)
        target = _compose_symmetry_keys(g, keys[s])
        j = findfirst(==(target), keys)
        isnothing(j) && throw(ArgumentError("$(typeof(action)) is not contained in the billiard's registered physical symmetry group"))
        p[s] = j
    end
    return p
end

# Return the exact cyclic keys of the physical rotational reconstruction
# sectors in the same order used by `full_boundary`.
function _rotation_sector_keys(billiard::Bi) where {Bi<:AbsBilliard}
    isempty(billiard.symmetries) && throw(ArgumentError("Billiard has no registered rotational symmetry"))
    firstsym = first(billiard.symmetries)
    firstsym isa NFoldRotation || throw(ArgumentError("Expected a rotational SymmetryRegistry; found $(typeof(firstsym))"))
    n = firstsym.order
    keys = Tuple{Int, Int}[_rotation_identity_key(n)]
    @inbounds for sym in billiard.symmetries
        sym isa NFoldRotation || throw(ArgumentError("SymmetryRegistry mixes rotational and reflection actions"))
        sym.order == n || throw(ArgumentError("SymmetryRegistry mixes rotation orders"))
        push!(keys, _rotation_key(sym))
    end
    length(unique(keys)) == length(keys) || throw(ArgumentError("SymmetryRegistry contains duplicate rotational actions"))
    return keys
end

# Compute how one N-fold rotation permutes the physical reconstruction sectors
# by exact cyclic-group addition.
function _rotation_sector_permutation(billiard::Bi, action::NFoldRotation) where {Bi<:AbsBilliard}
    keys = _rotation_sector_keys(billiard)
    g = _rotation_key(action)
    p = Vector{Int}(undef, length(keys))
    @inbounds for s in eachindex(keys)
        target = _compose_rotation_keys(g, keys[s])
        j = findfirst(==(target), keys)
        isnothing(j) && throw(ArgumentError("Rotation m = $(action.m) is not contained in the billiard's registered physical symmetry group"))
        p[s] = j
    end
    return p
end

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

# Construct the exact midpoint-index permutation induced by one physical
# reflection on an invariant boundary component.
function _boundary_symmetry_permutation(billiard::Bi, N::Int, symmetry::AbsReflection) where {Bi<:AbsBilliard}
    symmetry isa CompositeReflection && throw(ArgumentError("CompositeReflection contains several geometric generators"))
    sectors = _reflection_sector_permutation(billiard, symmetry)
    return _sector_node_permutation(N, sectors, _orientation_reversing(symmetry))
end

# Construct the exact midpoint-index permutation induced by one physical
# N-fold rotation on an invariant boundary component.
function _boundary_symmetry_permutation(billiard::Bi, N::Int, symmetry::NFoldRotation) where {Bi<:AbsBilliard}
    sectors = _rotation_sector_permutation(billiard, symmetry)
    return _sector_node_permutation(N, sectors, false)
end

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

"""
    symmetry_index_orbits(::Type{T}, billiard::Bi, N::Int, symmetry::AbsSymmetry, character...) where {T<:Real,Bi<:AbsBilliard} → orbits::SymmetryOrbitMap{T}

Construct the exact symmetry-orbit map for one complete invariant physical
boundary component.

## Description
The node permutation is derived from the billiard's registered physical
symmetry reconstruction and the integer midpoint-grid index. Boundary
coordinates are never inspected and no floating-point matching is performed.

The physical boundary component must be individually invariant under the
requested symmetry. Symmetries that exchange distinct physical components
are currently unsupported.

## Arguments
* `T::Type{<:Real}`: Real scalar type of the boundary discretization.
* `billiard::Bi`: Billiard defining the registered physical symmetry group.
* `N::Int`: Number of nodes on the complete invariant physical component.
* `symmetry::AbsSymmetry`: Active symmetry reduction.
* `character`: Requested irreducible-representation character or rotational sector.

## Returns
* `orbits::SymmetryOrbitMap{T}`: Exact boundary symmetry-orbit map.
"""
function symmetry_index_orbits(::Type{T}, billiard::Bi, N::Int, symmetry::XAxisReflection, character::Complex{T} = one(Complex{T})) where {T<:Real,Bi<:AbsBilliard}
    N % symmetry_node_multiple(symmetry) == 0 || throw(ArgumentError("XAxisReflection requires N divisible by $(symmetry_node_multiple(symmetry)); received N = $N"))
    id = collect(1:N)
    refl = _boundary_symmetry_permutation(billiard, N, symmetry)
    return _build_symmetry_orbit_map(T, N, Vector{Vector{Int}}([id, refl]), Complex{T}[one(Complex{T}), character])
end

function symmetry_index_orbits(::Type{T}, billiard::Bi, N::Int, symmetry::YAxisReflection, character::Complex{T} = one(Complex{T})) where {T<:Real,Bi<:AbsBilliard}
    N % symmetry_node_multiple(symmetry) == 0 || throw(ArgumentError("YAxisReflection requires N divisible by $(symmetry_node_multiple(symmetry)); received N = $N"))
    id = collect(1:N)
    refl = _boundary_symmetry_permutation(billiard, N, symmetry)
    return _build_symmetry_orbit_map(T, N, Vector{Vector{Int}}([id, refl]), Complex{T}[one(Complex{T}), character])
end

function symmetry_index_orbits(::Type{T}, billiard::Bi, N::Int, symmetry::DiagonalReflection, character::Complex{T} = one(Complex{T})) where {T<:Real,Bi<:AbsBilliard}
    N % symmetry_node_multiple(symmetry) == 0 || throw(ArgumentError("DiagonalReflection requires N divisible by $(symmetry_node_multiple(symmetry)); received N = $N"))
    id = collect(1:N)
    refl = _boundary_symmetry_permutation(billiard, N, symmetry)
    return _build_symmetry_orbit_map(T, N, Vector{Vector{Int}}([id, refl]), Complex{T}[one(Complex{T}), character])
end

function symmetry_index_orbits(::Type{T}, billiard::Bi, N::Int, symmetry::AntiDiagonalReflection, character::Complex{T} = one(Complex{T})) where {T<:Real,Bi<:AbsBilliard}
    N % symmetry_node_multiple(symmetry) == 0 || throw(ArgumentError("AntiDiagonalReflection requires N divisible by $(symmetry_node_multiple(symmetry)); received N = $N"))
    id = collect(1:N)
    refl = _boundary_symmetry_permutation(billiard, N, symmetry)
    return _build_symmetry_orbit_map(T, N, Vector{Vector{Int}}([id, refl]), Complex{T}[one(Complex{T}), character])
end

function symmetry_index_orbits(::Type{T}, billiard::Bi, N::Int, symmetry::XYAxisReflection, character_x::Complex{T} = one(Complex{T}), character_y::Complex{T} = one(Complex{T})) where {T<:Real,Bi<:AbsBilliard}
    N % symmetry_node_multiple(symmetry) == 0 || throw(ArgumentError("XYAxisReflection requires N divisible by $(symmetry_node_multiple(symmetry)); received N = $N"))
    id = collect(1:N)
    rx = _boundary_symmetry_permutation(billiard, N, XAxisReflection())
    ry = _boundary_symmetry_permutation(billiard, N, YAxisReflection())
    rxy = _compose_index_permutations(rx, ry)
    χxy = character_x * character_y
    return _build_symmetry_orbit_map(T, N, Vector{Vector{Int}}([id, rx, ry, rxy]), Complex{T}[one(Complex{T}), character_x, character_y, χxy])
end

function symmetry_index_orbits(::Type{T}, billiard::Bi, N::Int, symmetry::NFoldRotation, sector::Int = 0) where {T<:Real,Bi<:AbsBilliard}
    n = symmetry.order
    N % symmetry_node_multiple(symmetry) == 0 || throw(ArgumentError("NFoldRotation of order $n requires N divisible by $n; received N = $N"))
    g = _boundary_symmetry_permutation(billiard, N, NFoldRotation(n, 1; T = typeof(symmetry.angle)))
    perms = Vector{Vector{Int}}(undef, n)
    perms[1] = collect(1:N)
    @inbounds for l in 2:n
        perms[l] = _compose_index_permutations(g, perms[l - 1])
    end
    scales = Complex{T}[cis(T(2 * pi) * T(sector * l) / T(n)) for l in 0:n - 1]
    return _build_symmetry_orbit_map(T, N, perms, scales)
end

function symmetry_index_orbits(::Type{T}, billiard::Bi, N::Int, symmetry::CompositeReflection, characters::Vector{Complex{T}} = ones(Complex{T}, length(symmetry.reflections))) where {T<:Real,Bi<:AbsBilliard}
    isempty(symmetry.reflections) && throw(ArgumentError("CompositeReflection requires at least one reflection"))
    length(characters) == length(symmetry.reflections) || throw(DimensionMismatch("Received $(length(symmetry.reflections)) reflections but $(length(characters)) characters"))
    N % symmetry_node_multiple(symmetry) == 0 || throw(ArgumentError("CompositeReflection requires N divisible by $(symmetry_node_multiple(symmetry)); received N = $N"))
    genperms = Vector{Vector{Int}}()
    genscales = Complex{T}[]
    @inbounds for (k, ref) in enumerate(symmetry.reflections)
        χ = characters[k]
        if ref isa XYAxisReflection
            push!(genperms, _boundary_symmetry_permutation(billiard, N, XAxisReflection()))
            push!(genscales, χ)
            push!(genperms, _boundary_symmetry_permutation(billiard, N, YAxisReflection()))
            push!(genscales, χ)
        elseif ref isa XAxisReflection || ref isa YAxisReflection || ref isa DiagonalReflection || ref isa AntiDiagonalReflection
            push!(genperms, _boundary_symmetry_permutation(billiard, N, ref))
            push!(genscales, χ)
        else
            throw(ArgumentError("Unsupported reflection type $(typeof(ref)) in CompositeReflection"))
        end
    end
    perms = Vector{Vector{Int}}([collect(1:N)])
    scales = Complex{T}[one(Complex{T})]
    head = 1
    while head <= length(perms)
        p = perms[head]
        χ = scales[head]
        @inbounds for g in eachindex(genperms)
            pref = _compose_index_permutations(genperms[g], p)
            χnew = genscales[g] * χ
            j = findfirst(==(pref), perms)
            if isnothing(j)
                push!(perms, pref)
                push!(scales, χnew)
            elseif scales[j] != χnew
                throw(ArgumentError("Composite reflection parities are inconsistent"))
            end
        end
        head += 1
    end
    return _build_symmetry_orbit_map(T, N, perms, scales)
end