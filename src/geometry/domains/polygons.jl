
struct Polygon{T} <: AbsSimpleDomain where T<:Real
    boundary::Vector{AbsCurve}
    corners::Vector{SVector{2,T}}
    angles::Vector{T}
    id::Int64
end

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

#=
function Polygon(corners, boundary_conditions, id)
    type = eltype(corners[1])
    corners  = [SVector{2,type}(c) for c in corners]
    boundary = []
    corners_1 = CircularArray(corners)
    for i in 1:length(corners)
        line = LineSegment(corners_1[i],corners_1[i+1];bc = boundary_conditions[i])
        push!(boundary,line)
    end
    return Polygon(boundary,corners,id)
end
=#