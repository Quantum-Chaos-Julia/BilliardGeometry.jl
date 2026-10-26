
struct Polygon{T} <: AbsSimpleDomain where T<:Real
    boundary::Vector{AbsCurve}
    corners::Vector{SVector{2,T}}
    angles::Vector{T}
    id::Int64
end

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
