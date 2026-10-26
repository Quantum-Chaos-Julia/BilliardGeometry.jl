################################################################################
############################# SYMMETRY REGISTRY ################################
################################################################################
# `SymmetryRegistry` is stored on `billiard.symmetries`. Every concrete
# `AbsSymmetry` carries a stable `sym_id::Int`, assigned once (in declaration
# order) by `register_symmetries` at billiard-construction time, so a symmetry
# is looked up by id via `symmetry_of` instead of by positional index.
#
# The registry lists every non-identity symmetry-group element needed by
# `full_boundary`/`poincarebirkhoff.jl`, not just a minimal generating set.
################################################################################

"""
SymmetryRegistry <: AbstractVector{AbsSymmetry}

`SymmetryRegistry` stores a billiard's discrete symmetry-group generators
(every non-identity element `full_boundary`/`poincarebirkhoff.jl` need to
reconstruct the full physical boundary/Poincaré sections), each tagged with a
stable `sym_id` assigned by [`register_symmetries`](@ref).

## Description
Behaves like a plain `Vector{AbsSymmetry}` for iteration/indexing/`length`
(it is a thin `AbstractVector` wrapper), so existing code that iterates
`billiard.symmetries` is unaffected. The only new capability is
[`symmetry_of`](@ref), a `sym_id`-keyed lookup that does not depend on
storage order.

## Attributes
* `generators::Vector{AbsSymmetry}`: The registry's tagged symmetry generators, in declaration order.
"""
struct SymmetryRegistry <: AbstractVector{AbsSymmetry}
    generators::Vector{AbsSymmetry}
end

Base.size(reg::SymmetryRegistry) = size(reg.generators)
Base.getindex(reg::SymmetryRegistry, i::Int) = reg.generators[i]
Base.IndexStyle(::Type{<:SymmetryRegistry}) = IndexLinear()

"""
    register_symmetries(gens::AbsSymmetry...) → reg::SymmetryRegistry

Assigns each generator in `gens` a stable `sym_id` (`1`-based, in argument
order) and returns the resulting [`SymmetryRegistry`](@ref).

## Description
Called once per billiard, at construction time, with the same full list of
non-identity symmetry-group elements the billiard's `symmetries` field used
to hold directly (e.g. `register_symmetries(YAxisReflection(),
XYAxisReflection(), XAxisReflection())` for a `D₂`-symmetric billiard). The
assigned `sym_id`s are the stable identifiers curves tag themselves with via
[`SymmetryWall`](@ref), instead of relying on `billiard.symmetries`' storage
position.

## Arguments
* `gens`: The billiard's non-identity discrete symmetry-group elements, in the order `full_boundary` should apply them.

## Returns
* `reg`: A [`SymmetryRegistry`](@ref) with each generator tagged by its `sym_id`.
"""
function register_symmetries(gens::AbsSymmetry...)
    tagged = ntuple(i -> _with_sym_id(gens[i], i), length(gens))
    return SymmetryRegistry(AbsSymmetry[tagged...])
end

"""
    symmetry_of(reg::SymmetryRegistry, sym_id::Int) → sym::AbsSymmetry

Looks up the generator in `reg` tagged with the given `sym_id`, independent
of its storage position.
"""
function symmetry_of(reg::SymmetryRegistry, sym_id::Int)
    i = findfirst(s -> s.sym_id == sym_id, reg.generators)
    isnothing(i) && throw(ArgumentError("No symmetry registered with sym_id=$sym_id"))
    return reg.generators[i]
end
