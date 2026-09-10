"""
    Polygon{T}<:AbsSimpleDomain

Polygonal simple domain defined by an ordered collection of corners.

The boundary consists of line segments joining consecutive corners, including
the final corner back to the first. `angles` stores the angle associated with
each corner.

## Fields
* `boundary::Vector{AbsCurve}`: Ordered polygon edges.
* `corners::Vector{SVector{2,T}}`: Ordered polygon vertices.
* `angles::Vector{T}`: Angles associated with the polygon corners.
* `id::Int64`: Domain identifier.
"""
struct Polygon{T} <: AbsSimpleDomain where T<:Real
    boundary::Vector{AbsCurve}
    corners::Vector{SVector{2,T}}
    angles::Vector{T}
    id::Int64
end

"""
    Polygon(corners,id;bcs=[SpecularReflection() for i in corners])

Construct a polygonal domain from an ordered collection of corners.

Consecutive corners are connected by `LineSegment`s, with the final corner
connected back to the first.

## Arguments
* `corners`: Ordered polygon vertices.
* `id`: Domain identifier.
* `bcs`: Boundary conditions assigned to the polygon edges.

## Returns
* `Polygon`: Constructed polygonal domain.
"""
function Polygon(corners, id; bcs = [SpecularReflection() for i in corners] )
    M = length(corners)
    type = eltype(corners[1])
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