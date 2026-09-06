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

Returns the smallest full-boundary node count multiple compatible with folding
onto a fundamental domain under `symmetry` (e.g. `2` for a single reflection,
so that the boundary can be split evenly across the reflection axis).

!!! note "Migration status"
    API scaffold only: this generic function currently has no methods. Concrete
    methods for [`XAxisReflection`](@ref), [`YAxisReflection`](@ref),
    [`XYAxisReflection`](@ref) and [`NFoldRotation`](@ref) are pending (see the
    `QuantumBilliardsTests` migration plan).
"""
function symmetry_node_multiple end

"""
    symmetry_index_orbits(::Type{T}, xy::AbstractVector{SVector{2,T}}, symmetry::AbsSymmetry) where {T<:Real} → orbits::SymmetryOrbitMap{T}

Constructs the [`SymmetryOrbitMap`](@ref) folding the fully discretized boundary
points `xy` onto a fundamental domain under `symmetry`.

## Arguments
* `T`: Real scalar type used by the boundary discretization.
* `xy`: Complete (unfolded) boundary points in Cartesian coordinates.
* `symmetry`: The discrete symmetry the boundary points are invariant under.

## Returns
* `orbits`: A [`SymmetryOrbitMap{T}`](@ref) instance.

!!! note "Migration status"
    API scaffold only: this generic function currently has no methods pending
    the boundary-integral-method solver migration (see the
    `QuantumBilliardsTests` migration plan).
"""
function symmetry_index_orbits end
