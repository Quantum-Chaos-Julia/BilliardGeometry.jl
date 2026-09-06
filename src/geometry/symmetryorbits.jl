"""
SymmetryOrbitMap{T}

`SymmetryOrbitMap` is a concrete type representing the folding of a fully
discretized (periodic) boundary onto a fundamental domain under a discrete
[`AbsSymmetry`](@ref) group.

## Description
Boundary-integral method (BIM) solvers discretize unknown boundary densities
directly (rather than expanding them in an [`AbsBasis`](@ref)), so a discrete
symmetry cannot be baked into basis functions the way
[`CornerAdaptedFourierBessel`](@ref) does for basis solvers. Instead, the
*complete* physical boundary is discretized, and a `SymmetryOrbitMap` records,
for every full-boundary node, which fundamental-domain node it is the
symmetry image of, together with the irreducible-representation phase factor
relating the two. This lets a BIM solver assemble its Fredholm matrix only on
the fundamental-domain indices while still summing source contributions over
every symmetry image of the full boundary.

## Attributes
* `fundamental_indices::Vector{Int}`: Indices into the full-boundary node array kept as the fundamental-domain representatives.
* `orbit_of::Vector{Int}`: For each full-boundary node index, the fundamental-domain index (position in `fundamental_indices`) it is folded onto.
* `phase::Vector{Complex{T}}`: Per full-boundary-node irreducible-representation phase factor relating the node to its fundamental-domain representative (`one(Complex{T})` for the trivial representation).
* `full_size::Int`: Number of nodes on the complete (unfolded) boundary.
* `fundamental_size::Int`: Number of nodes on the fundamental domain, `length(fundamental_indices)`.

## API
The following functions can be evaluated for any `SymmetryOrbitMap`:
- [`fundamental_size`](@ref)
- `Base.length`
"""
struct SymmetryOrbitMap{T<:Real}
    fundamental_indices::Vector{Int}
    orbit_of::Vector{Int}
    phase::Vector{Complex{T}}
    full_size::Int
    fundamental_size::Int
end

"""
    fundamental_size(orbits::SymmetryOrbitMap) → n::Int

Returns the number of nodes on the fundamental domain, `n = orbits.fundamental_size`.
"""
fundamental_size(orbits::SymmetryOrbitMap) = orbits.fundamental_size

"""
    length(orbits::SymmetryOrbitMap) → n::Int

Returns the number of nodes on the complete (unfolded) boundary, `n = orbits.full_size`.
"""
Base.length(orbits::SymmetryOrbitMap) = orbits.full_size

"""
    symmetry_node_multiple(symmetry::AbsSymmetry) → n::Int

Returns the full-boundary node count multiple required for
[`symmetry_index_orbits`](@ref) to fold the boundary onto a fundamental
domain under `symmetry` using exact integer index permutations (no
floating-point point matching).

## Description
For [`XAxisReflection`](@ref), [`YAxisReflection`](@ref) and
[`XYAxisReflection`](@ref) the node count must be divisible by `4` (matching
the `-develop` reference solvers' node-count sizing, which requires this even
though the reflection group itself has order `2`, so that the boundary can
also accommodate the combined `D₂` action exactly). For
[`NFoldRotation`](@ref), the multiple is the rotation order `sym.order`.
"""
symmetry_node_multiple(::XAxisReflection) = 4
symmetry_node_multiple(::YAxisReflection) = 4
symmetry_node_multiple(::XYAxisReflection) = 4
symmetry_node_multiple(sym::NFoldRotation) = sym.order

# Canonical periodic boundary index actions (boundary assumed sampled by
# midpoint nodes `s_mid(k,N) = 2π(k-1/2)/N` in canonical orientation, exactly
# as ported to `QuantumBilliards.jl`/`BilliardGeometry.jl`'s BIM `evaluate_points`
# methods). All actions below are exact integer permutations of node indices,
# so no floating-point tolerance/point-matching is required.
@inline _idx_reflect_x(q::Int, N::Int) = mod1(N-q+1, N)
@inline _idx_reflect_y(q::Int, N::Int) = mod1(N÷2-q+1, N)
@inline _idx_rotate_pi(q::Int, N::Int) = mod1(q+N÷2, N)
@inline _idx_rotate(q::Int, N::Int, n::Int, l::Int) = mod1(q+l*(N÷n), N)

"""
    _build_symmetry_orbit_map(::Type{T}, N::Int, perms::Vector{Vector{Int}}) where {T<:Real} → orbits::SymmetryOrbitMap{T}

Builds a [`SymmetryOrbitMap`](@ref) from a group of exact full-boundary index
permutations `perms` (one per group element, `perms[1]` the identity), using
the trivial representation (`phase` is `one(Complex{T})` everywhere).
"""
function _build_symmetry_orbit_map(::Type{T}, N::Int, perms::Vector{Vector{Int}}) where {T<:Real}
    ng = length(perms)
    N%ng==0 || throw(ArgumentError("Node count N=$N must be divisible by symmetry-group order $ng"))
    nf = N÷ng
    fundamental_indices = Vector{Int}(undef, nf)
    orbit_of = Vector{Int}(undef, N)
    phase = ones(Complex{T}, N)
    seen = falses(N)
    b = 0
    @inbounds for q in 1:N
        seen[q] && continue
        b += 1
        fundamental_indices[b] = q
        for g in 1:ng
            qi = perms[g][q]
            orbit_of[qi] = b
            seen[qi] = true
        end
    end
    return SymmetryOrbitMap{T}(fundamental_indices, orbit_of, phase, N, nf)
end

"""
    symmetry_index_orbits(::Type{T}, xy::AbstractVector{SVector{2,T}}, symmetry::AbsSymmetry) where {T<:Real} → orbits::SymmetryOrbitMap{T}

Constructs the [`SymmetryOrbitMap`](@ref) folding the fully discretized
boundary points `xy` onto a fundamental domain under `symmetry`.

## Description
The reduction uses the exact integer index permutation the symmetry induces
on a canonically ordered periodic boundary sampling (the same convention used
by the BIM solvers' `evaluate_points` methods), not floating-point nearest-
neighbor matching on `xy`: `length(xy)` must already be a multiple of
[`symmetry_node_multiple`](@ref)`(symmetry)`. Only the trivial representation
is supported (`phase` is `one(Complex{T})` for every node); odd/anti-symmetric
BIM sectors are not yet implemented.

## Arguments
* `T`: Real scalar type used by the boundary discretization.
* `xy`: Complete (unfolded) boundary points in Cartesian coordinates, canonically ordered.
* `symmetry`: The discrete symmetry the boundary points are invariant under.

## Returns
* `orbits`: A [`SymmetryOrbitMap{T}`](@ref) instance.
"""
function symmetry_index_orbits(::Type{T}, xy::AbstractVector{SVector{2,T}}, symmetry::XAxisReflection) where {T<:Real}
    N = length(xy)
    N%symmetry_node_multiple(symmetry)==0 || throw(ArgumentError("XAxisReflection requires N divisible by 4; received N=$N"))
    id = collect(1:N)
    refl = [_idx_reflect_x(q,N) for q in 1:N]
    return _build_symmetry_orbit_map(T, N, [id,refl])
end

function symmetry_index_orbits(::Type{T}, xy::AbstractVector{SVector{2,T}}, symmetry::YAxisReflection) where {T<:Real}
    N = length(xy)
    N%symmetry_node_multiple(symmetry)==0 || throw(ArgumentError("YAxisReflection requires N divisible by 4; received N=$N"))
    id = collect(1:N)
    refl = [_idx_reflect_y(q,N) for q in 1:N]
    return _build_symmetry_orbit_map(T, N, [id,refl])
end

function symmetry_index_orbits(::Type{T}, xy::AbstractVector{SVector{2,T}}, symmetry::XYAxisReflection) where {T<:Real}
    N = length(xy)
    N%symmetry_node_multiple(symmetry)==0 || throw(ArgumentError("XYAxisReflection requires N divisible by 4; received N=$N"))
    id = collect(1:N)
    rx = [_idx_reflect_x(q,N) for q in 1:N]
    ry = [_idx_reflect_y(q,N) for q in 1:N]
    rxy = [_idx_rotate_pi(q,N) for q in 1:N]
    return _build_symmetry_orbit_map(T, N, [id,rx,ry,rxy])
end

function symmetry_index_orbits(::Type{T}, xy::AbstractVector{SVector{2,T}}, symmetry::NFoldRotation) where {T<:Real}
    N = length(xy)
    n = symmetry_node_multiple(symmetry)
    N%n==0 || throw(ArgumentError("NFoldRotation of order $n requires N divisible by $n; received N=$N"))
    perms = [[_idx_rotate(q,N,n,l) for q in 1:N] for l in 0:n-1]
    return _build_symmetry_orbit_map(T, N, perms)
end
