"""
    CompositeCurve{T} <: AbsCompositeCurve

Consecutive curves forming a boundary, with cumulative arc lengths.

## Description
The constructor `CompositeCurve(subcurves)` orders the curves with
[`connect_curves`](@ref) and fills the remaining fields. The numeric type `T` is
taken from the first curve. Used to define the arc-length coordinate `s` of
[`PoincareBirkhoff`](@ref) sections.

## Attributes
* `subcurves::Vector{AbsCurve}`: Connected curves in traversal order.
* `end_lengths::Vector{T}`: Arc length at the start of each curve followed by the total length, so `end_lengths[i]` is the arc length at the start of `subcurves[i]` and `end_lengths[end]` equals `length`.
* `corners::Vector{SVector{2,T}}`: Start point of the first curve followed by the end point of every curve.
* `length::T`: Total arc length.
"""
struct CompositeCurve{T} <: AbsCompositeCurve where {T<:Real}
    subcurves::Vector{AbsCurve}
    end_lengths::Vector{T}
    corners::Vector{SVector{2,T}}
    length::T
end

function CompositeCurve(subcurves) 
    type = typeof(subcurves[1]).parameters[1] 
    end_lengths =[zero(type)]
    corners = Vector{SVector{2,type}}(undef,0)

    crvs = connect_curves(subcurves)
    push!(corners, curve(crvs[1],0.0))

    L = zero(type)
    for crv in crvs
        L += crv.length
        push!(end_lengths, L)
        push!(corners, curve(crv,1.0)) 
    end
    return CompositeCurve{type}(crvs,end_lengths,corners,L)
end
