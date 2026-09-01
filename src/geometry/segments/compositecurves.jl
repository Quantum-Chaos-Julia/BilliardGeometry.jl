"""
    CompositeCurve{T}<:AbsCompositeCurve

Connected sequence of boundary curves treated as a single composite curve.

The subcurves are ordered and connected, `end_lengths` stores the cumulative
arc length at each subcurve endpoint, and `corners` stores the corresponding
junction points.

## Fields
* `subcurves::Vector{AbsCurve}`: Connected component curves.
* `end_lengths::Vector{T}`: Cumulative arc lengths, starting at zero.
* `corners::Vector{SVector{2,T}}`: Subcurve junction points.
* `length::T`: Total composite-curve length.
"""
struct CompositeCurve{T<:Real}<:AbsCompositeCurve
    subcurves::Vector{AbsCurve}
    end_lengths::Vector{T}
    corners::Vector{SVector{2,T}}
    length::T
end

"""
    CompositeCurve(subcurves)

Construct a connected composite curve from an ordered collection of subcurves.

Each constituent curve is parametrized over `t∈[0,1]`. The curves are first
passed through `connect_curves`, after which cumulative arc lengths and
junction points are constructed.

## Arguments
* `subcurves`: Ordered boundary curves forming the composite curve.

## Returns
* `CompositeCurve`: Connected composite curve with cumulative lengths and corners.
"""
function CompositeCurve(subcurves)
    isempty(subcurves)&&throw(ArgumentError("subcurves must contain at least one curve"))
    crvs=connect_curves(subcurves)
    T=typeof(crvs[1].length)
    end_lengths=Vector{T}(undef,length(crvs)+1)
    corners=Vector{SVector{2,T}}(undef,length(crvs)+1)
    end_lengths[1]=zero(T)
    corners[1]=SVector{2,T}(curve(crvs[1],zero(T)))
    L=zero(T)
    @inbounds for i in eachindex(crvs)
        L+=crvs[i].length
        end_lengths[i+1]=L
        corners[i+1]=SVector{2,T}(curve(crvs[i],one(T)))
    end
    return CompositeCurve{T}(crvs,end_lengths,corners,L)
end