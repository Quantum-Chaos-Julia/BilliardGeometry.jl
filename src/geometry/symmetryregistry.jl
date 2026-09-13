################################################################################
############################# SYMMETRY REGISTRY ################################
################################################################################
# `SymmetryRegistry` replaces the old plain `Vector{AbsSymmetry}` stored on
# `billiard.symmetries`: every concrete `AbsSymmetry` generator now carries a
# stable `sym_id::Int`, assigned once (in declaration order) by
# `register_symmetries` at billiard-construction time. This is Layer 1 of the
# symmetry-framework split (see Step 15 of the migration plan): a purely
# geometric, ordering-independent way to look up "the billiard's k-th
# registered symmetry generator" by a stable id instead of a fragile
# positional index (`symmetries[sym_sector-1]`).
#
# `SymmetryRegistry` still lists every non-identity symmetry-group element
# `full_boundary`/`poincarebirkhoff.jl` need directly (not just a minimal
# generating set, exactly as the old `Vector{AbsSymmetry}` did), so every
# consumer that used to iterate/index `billiard.symmetries` keeps working
# unchanged; only positional-index lookups are rewritten to use `symmetry_of`.
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
