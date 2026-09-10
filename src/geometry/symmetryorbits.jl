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
* `fund_to_full::Matrix{Int}`: `fund_to_full[g,b]` is the full-boundary node index of the `g`-th group image of fundamental node `b`.
* `fund_to_scale::Matrix{Complex{T}}`: `fund_to_scale[g,b]` is the irrep factor associated with `fund_to_full[g,b]`.

## API
The following functions can be evaluated for any `SymmetryOrbitMap`:
- [`fundamental_size`](@ref)
- [`full_size`](@ref)
- [`orbit_size`](@ref)
- [`symmetry_orbit`](@ref)
- `Base.length`
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

Returns the number of nodes on the fundamental domain, `n = orbits.fundamental_size`.
"""
fundamental_size(orbits::SymmetryOrbitMap) = orbits.fundamental_size

"""
    length(orbits::SymmetryOrbitMap) → n::Int

Returns the number of nodes on the complete (unfolded) boundary, `n = orbits.full_size`.
"""
Base.length(orbits::SymmetryOrbitMap) = orbits.full_size

"""
    full_size(orbits::SymmetryOrbitMap) → n::Int

Returns the number of nodes on the complete (unfolded) boundary, `n = orbits.full_size`.
"""
full_size(orbits::SymmetryOrbitMap) = orbits.full_size

"""
    orbit_size(orbits::SymmetryOrbitMap) → ng::Int

Returns the symmetry-group order (number of group images per fundamental node), `size(orbits.fund_to_full,1)`.
"""
orbit_size(orbits::SymmetryOrbitMap) = size(orbits.fund_to_full, 1)

"""
    symmetry_orbit(orbits::SymmetryOrbitMap, b::Int) → (qs, χs)

Returns the full-boundary indices and irrep factors representing the complete
symmetry orbit of fundamental node `b`: if `qs,χs = symmetry_orbit(orbits,b)`,
then a reduced boundary value `u_b` generates the full orbit according to
`u[qs[l]] = χs[l]*u_b`. The returned arrays are views into `orbits` and
therefore allocate no copies.
"""
@inline function symmetry_orbit(orbits::SymmetryOrbitMap, b::Int)
    return @view(orbits.fund_to_full[:,b]), @view(orbits.fund_to_scale[:,b])
end

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
also accommodate the combined `D₂` action exactly). [`DiagonalReflection`](@ref)
and [`AntiDiagonalReflection`](@ref) require divisibility by `8`. For
[`NFoldRotation`](@ref), the multiple is the rotation order `sym.order`. For a
[`CompositeReflection`](@ref), the multiple is the least common multiple of
its constituent reflections' multiples.
"""
symmetry_node_multiple(::XAxisReflection) = 4
symmetry_node_multiple(::YAxisReflection) = 4
symmetry_node_multiple(::XYAxisReflection) = 4
symmetry_node_multiple(::DiagonalReflection) = 8
symmetry_node_multiple(::AntiDiagonalReflection) = 8
symmetry_node_multiple(sym::NFoldRotation) = sym.order
symmetry_node_multiple(sym::CompositeReflection) = foldl(lcm, (symmetry_node_multiple(ref) for ref in sym.reflections); init=1)

# Canonical periodic boundary index actions (boundary assumed sampled by
# midpoint nodes `s_mid(k,N) = 2π(k-1/2)/N` in canonical orientation, exactly
# as ported to `QuantumBilliards.jl`/`BilliardGeometry.jl`'s BIM `evaluate_points`
# methods). All actions below are exact integer permutations of node indices,
# so no floating-point tolerance/point-matching is required.
@inline _idx_reflect_x(q::Int, N::Int) = mod1(N-q+1, N)
@inline _idx_reflect_y(q::Int, N::Int) = mod1(N÷2-q+1, N)
@inline _idx_reflect_diag_plus(q::Int, N::Int) = mod1(N÷4-q+1, N)
@inline _idx_reflect_diag_minus(q::Int, N::Int) = mod1(3*N÷4-q+1, N)
@inline _idx_rotate_pi(q::Int, N::Int) = mod1(q+N÷2, N)
@inline _idx_rotate(q::Int, N::Int, n::Int, l::Int) = mod1(q+l*(N÷n), N)

"""
    _build_symmetry_orbit_map(::Type{T}, N::Int, perms::Vector{Vector{Int}}, scales::Vector{Complex{T}}=ones(Complex{T},length(perms))) where {T<:Real} → orbits::SymmetryOrbitMap{T}

Builds a [`SymmetryOrbitMap`](@ref) from a group of exact full-boundary index
permutations `perms` (one per group element, `perms[1]` the identity) and
their corresponding irrep factors `scales` (defaulting to the trivial
representation, `phase` equal to `one(Complex{T})` everywhere).
"""
function _build_symmetry_orbit_map(::Type{T}, N::Int, perms::Vector{Vector{Int}}, scales::Vector{Complex{T}}=ones(Complex{T}, length(perms))) where {T<:Real}
    ng = length(perms)
    N%ng==0 || throw(ArgumentError("Node count N=$N must be divisible by symmetry-group order $ng"))
    length(scales)==ng || throw(DimensionMismatch("Received $ng permutations but $(length(scales)) irrep factors"))
    nf = N÷ng
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
        fundamental_indices[b] = q
        for g in 1:ng
            qi = perms[g][q]
            χ = scales[g]
            orbit_of[qi] = b
            phase[qi] = χ
            fund_to_full[g,b] = qi
            fund_to_scale[g,b] = χ
            seen[qi] = true
        end
    end
    return SymmetryOrbitMap{T}(fundamental_indices, orbit_of, phase, N, nf, fund_to_full, fund_to_scale)
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
[`symmetry_node_multiple`](@ref)`(symmetry)`. The per-node `phase` factors are
the irreducible-representation characters returned by
[`symmetry_irrep_character`](@ref), so odd/anti-symmetric BIM sectors are
supported for every `symmetry` with a `symmetry_irrep_character` method.

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
    χ = symmetry_irrep_character(T, symmetry)
    return _build_symmetry_orbit_map(T, N, [id,refl], Complex{T}[one(Complex{T}), χ])
end

function symmetry_index_orbits(::Type{T}, xy::AbstractVector{SVector{2,T}}, symmetry::YAxisReflection) where {T<:Real}
    N = length(xy)
    N%symmetry_node_multiple(symmetry)==0 || throw(ArgumentError("YAxisReflection requires N divisible by 4; received N=$N"))
    id = collect(1:N)
    refl = [_idx_reflect_y(q,N) for q in 1:N]
    χ = symmetry_irrep_character(T, symmetry)
    return _build_symmetry_orbit_map(T, N, [id,refl], Complex{T}[one(Complex{T}), χ])
end

function symmetry_index_orbits(::Type{T}, xy::AbstractVector{SVector{2,T}}, symmetry::XYAxisReflection) where {T<:Real}
    N = length(xy)
    N%symmetry_node_multiple(symmetry)==0 || throw(ArgumentError("XYAxisReflection requires N divisible by 4; received N=$N"))
    id = collect(1:N)
    rx = [_idx_reflect_x(q,N) for q in 1:N]
    ry = [_idx_reflect_y(q,N) for q in 1:N]
    rxy = [_idx_rotate_pi(q,N) for q in 1:N]
    # χ_x/χ_y are the characters of the individual axis-reflection generators
    # (same convention as the `CompositeReflection` expansion below); χ_xy, the
    # character of the combined π-rotation, is `symmetry_irrep_character`.
    χ_x = Complex{T}(symmetry.parity_y)
    χ_y = Complex{T}(symmetry.parity_x)
    χ_xy = symmetry_irrep_character(T, symmetry)
    return _build_symmetry_orbit_map(T, N, [id,rx,ry,rxy], Complex{T}[one(Complex{T}), χ_x, χ_y, χ_xy])
end

function symmetry_index_orbits(::Type{T}, xy::AbstractVector{SVector{2,T}}, symmetry::NFoldRotation) where {T<:Real}
    N = length(xy)
    n = symmetry_node_multiple(symmetry)
    N%n==0 || throw(ArgumentError("NFoldRotation of order $n requires N divisible by $n; received N=$N"))
    perms = [[_idx_rotate(q,N,n,l) for q in 1:N] for l in 0:n-1]
    scales = Complex{T}[cis(T(2*pi)*T(symmetry.sector*l)/T(n)) for l in 0:n-1]
    return _build_symmetry_orbit_map(T, N, perms, scales)
end

"""
    symmetry_index_orbits(::Type{T}, xy::AbstractVector{SVector{2,T}}, symmetry::DiagonalReflection) where {T<:Real} → orbits::SymmetryOrbitMap{T}

Builds the exact two-element boundary orbits generated by reflection across
the `y=x` diagonal, using the irrep factor `symmetry_irrep_character(T,symmetry)`.
"""
function symmetry_index_orbits(::Type{T}, xy::AbstractVector{SVector{2,T}}, symmetry::DiagonalReflection) where {T<:Real}
    N = length(xy)
    N%symmetry_node_multiple(symmetry)==0 || throw(ArgumentError("DiagonalReflection requires N divisible by 8; received N=$N"))
    id = collect(1:N)
    refl = [_idx_reflect_diag_plus(q,N) for q in 1:N]
    χ = symmetry_irrep_character(T, symmetry)
    return _build_symmetry_orbit_map(T, N, [id,refl], Complex{T}[one(Complex{T}), χ])
end

"""
    symmetry_index_orbits(::Type{T}, xy::AbstractVector{SVector{2,T}}, symmetry::AntiDiagonalReflection) where {T<:Real} → orbits::SymmetryOrbitMap{T}

Builds the exact two-element boundary orbits generated by reflection across
the `y=-x` anti-diagonal, using the irrep factor `symmetry_irrep_character(T,symmetry)`.
"""
function symmetry_index_orbits(::Type{T}, xy::AbstractVector{SVector{2,T}}, symmetry::AntiDiagonalReflection) where {T<:Real}
    N = length(xy)
    N%symmetry_node_multiple(symmetry)==0 || throw(ArgumentError("AntiDiagonalReflection requires N divisible by 8; received N=$N"))
    id = collect(1:N)
    refl = [_idx_reflect_diag_minus(q,N) for q in 1:N]
    χ = symmetry_irrep_character(T, symmetry)
    return _build_symmetry_orbit_map(T, N, [id,refl], Complex{T}[one(Complex{T}), χ])
end

"""
    symmetry_index_orbits(::Type{T}, xy::AbstractVector{SVector{2,T}}, symmetry::CompositeReflection) where {T<:Real} → orbits::SymmetryOrbitMap{T}

Builds the complete boundary-orbit map generated by a [`CompositeReflection`](@ref).

The constituent reflections are first expanded into primitive exact index
generators (`XYAxisReflection` contributes both coordinate-axis reflections).
Starting from the identity, the complete finite symmetry group is generated
by closure under composition; the irrep factor of a generated action is the
product of the factors of its generators. If the same geometric action is
generated with two different factors, the requested reflection parities are
inconsistent and an `ArgumentError` is thrown.
"""
function symmetry_index_orbits(::Type{T}, xy::AbstractVector{SVector{2,T}}, symmetry::CompositeReflection) where {T<:Real}
    N = length(xy)
    isempty(symmetry.reflections) && throw(ArgumentError("CompositeReflection requires at least one reflection"))
    genperms = Vector{Vector{Int}}()
    genscales = Complex{T}[]
    for ref in symmetry.reflections
        if ref isa XAxisReflection
            push!(genperms, [_idx_reflect_x(q,N) for q in 1:N])
            push!(genscales, symmetry_irrep_character(T, ref))
        elseif ref isa YAxisReflection
            push!(genperms, [_idx_reflect_y(q,N) for q in 1:N])
            push!(genscales, symmetry_irrep_character(T, ref))
        elseif ref isa XYAxisReflection
            push!(genperms, [_idx_reflect_x(q,N) for q in 1:N])
            push!(genscales, Complex{T}(ref.parity_y))
            push!(genperms, [_idx_reflect_y(q,N) for q in 1:N])
            push!(genscales, Complex{T}(ref.parity_x))
        elseif ref isa DiagonalReflection
            push!(genperms, [_idx_reflect_diag_plus(q,N) for q in 1:N])
            push!(genscales, symmetry_irrep_character(T, ref))
        elseif ref isa AntiDiagonalReflection
            push!(genperms, [_idx_reflect_diag_minus(q,N) for q in 1:N])
            push!(genscales, symmetry_irrep_character(T, ref))
        else
            throw(ArgumentError("Unsupported reflection type $(typeof(ref)) in CompositeReflection"))
        end
    end
    perms = Vector{Int}[collect(1:N)]
    scales = Complex{T}[one(Complex{T})]
    head = 1
    while head <= length(perms)
        p = perms[head]
        χ = scales[head]
        for g in eachindex(genperms)
            gp = genperms[g]
            pref = [gp[p[q]] for q in 1:N]
            χnew = genscales[g]*χ
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
