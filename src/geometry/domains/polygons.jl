"""
    Polygon{T} <: AbsSimpleDomain

Simple domain bounded by straight [`LineSegment`](@ref)s through a list of
corners, edge `i` running from `corners[i]` to `corners[i+1]` (periodically).

## Attributes
* `boundary::Vector{AbsCurve}`: The edges, `boundary[i]` has `segment_id = i`.
* `corners::Vector{SVector{2,T}}`: Vertices in traversal order.
* `angles::Vector{T}`: Interior angle at each corner, between the edges meeting there, in `[0,π]`.
* `id::Int64`: Id of the domain, also the `domain_id` of every edge.
"""
struct Polygon{T} <: AbsSimpleDomain where T<:Real
    boundary::Vector{AbsCurve}
    corners::Vector{SVector{2,T}}
    angles::Vector{T}
    id::Int64
end

"""
    Polygon(corners, id; bcs = [SpecularReflection() for c in corners]) → polygon::Polygon

Construct a [`Polygon`](@ref) from its `corners`, given counterclockwise so the
interior lies on the left of every edge. The numeric type is the promotion of
the corner entries.

## Arguments
* `corners`: Vector of 2-element vectors.
* `id`: Id of the domain.

## Keyword arguments
* `bcs` : One boundary condition per edge, `bcs[i]` for the edge from `corners[i]` to `corners[i+1]`.
"""
function Polygon(corners, id; bcs = [SpecularReflection() for i in corners] )
    M = length(corners)
    type = promote_type((eltype(c) for c in corners)...)
    corners  = [SVector{2,type}(c) for c in corners]
    boundary = Vector{AbsCurve}(undef,M)
    angles = Vector{type}(undef,M)
    corners_1 = CircularArray(corners)
    for i in 1:length(corners)
        line = LineSegment(corners_1[i],corners_1[i+1]; bc = bcs[i],  domain_id=id, segment_id=i)
        phi = angle(corners_1[i+1] .- corners_1[i], corners_1[i-1] .- corners_1[i])
        boundary[i] = line
        angles[i] = phi
    end
    return Polygon(boundary,corners,angles,id)
end
