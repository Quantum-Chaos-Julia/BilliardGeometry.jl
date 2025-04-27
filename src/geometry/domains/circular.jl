
struct CircleCap{T} <: AbsSimpleDomain where T<:Real
    boundary::Vector{AbsCurve}
    corners::Vector{SVector{2,T}}
    id::Int64
end


struct CircleWedge{T} <: AbsSimpleDomain where T<:Real
    boundary::Vector{AbsCurve}
    corners::Vector{SVector{2,T}}
    center::SVector{2,T}
    id::Int64
end

function CircleWedge(R, arc_angle, shift_angle, x0, y0, id; bcs::Vector{AbsBoundaryCondition} = [SpecularReflection(), SpecularReflection(),SpecularReflection()])
    center = SVector(x0,y0)
    circle = CircleSegment(R, arc_angle, shift_angle, x0, y0; bc = bcs[1])
    pt0 = curve(circle, 0.0)
    pt1 = curve(circle, 1.0)
    line1 = LineSegment(pt1, center; bc = bcs[2]) 
    line2 = LineSegment(center, pt0; bc = bc = bcs[3]) 
    boundary = [circle, line1, line2]
    corners = [pt0, pt1, center]
    return CircleWedge{typeof(x0)}(boundary,corners,center,id)
end