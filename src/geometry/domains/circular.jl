"""
    CircleCap{T}<:AbsSimpleDomain

Simple domain representing a circular cap.

## Fields
* `boundary::Vector{AbsCurve}`: Ordered boundary curves enclosing the domain.
* `corners::Vector{SVector{2,T}}`: Corner or boundary-junction points.
* `id::Int64`: Domain identifier.
"""
struct CircleCap{T} <: AbsSimpleDomain where T<:Real
    boundary::Vector{AbsCurve}
    corners::Vector{SVector{2,T}}
    id::Int64
end

"""
    CircleWedge{T}<:AbsSimpleDomain

Simple circular wedge domain bounded by one circular segment and two line
segments connecting the circular arc to `center`.

## Fields
* `boundary::Vector{AbsCurve}`: Ordered boundary curves of the wedge.
* `corners::Vector{SVector{2,T}}`: Boundary junction points.
* `center::SVector{2,T}`: Center of the circular wedge.
* `id::Int64`: Domain identifier.
"""
struct CircleWedge{T} <: AbsSimpleDomain where T<:Real
    boundary::Vector{AbsCurve}
    corners::Vector{SVector{2,T}}
    center::SVector{2,T}
    id::Int64
end

"""
    CircleWedge(R,arc_angle,shift_angle,x0,y0,id;bcs::Vector{AbsBoundaryCondition}=[SpecularReflection(),SpecularReflection(),SpecularReflection()])

Construct a circular wedge centered at `(x0,y0)`.

The circular boundary is parametrized over `t∈[0,1]`, beginning at
`shift_angle` and spanning `arc_angle`. Its two endpoints are connected to the
center by line segments.

## Arguments
* `R`: Radius of the circular segment.
* `arc_angle`: Angular extent of the circular segment.
* `shift_angle`: Starting polar angle of the circular segment.
* `x0`: x-coordinate of the wedge center.
* `y0`: y-coordinate of the wedge center.
* `id`: Domain identifier.
* `bcs::Vector{AbsBoundaryCondition}`: Boundary conditions for the circular segment and the two radial line segments.

## Returns
* `CircleWedge`: Constructed circular wedge domain.
"""
function CircleWedge(R, arc_angle, shift_angle, x0, y0, id; bcs::Vector{AbsBoundaryCondition} = [SpecularReflection(), SpecularReflection(),SpecularReflection()])
    center = SVector(x0,y0)
    circle = CircleSegment(R, arc_angle, shift_angle, x0, y0; bc = bcs[1], domain_id=id, segment_id=1)
    pt0 = curve(circle, 0.0)
    pt1 = curve(circle, 1.0)
    line1 = LineSegment(pt1, center; bc = bcs[2], domain_id=id, segment_id=2) 
    line2 = LineSegment(center, pt0; bc = bc = bcs[3], domain_id=id, segment_id=3) 
    boundary = [circle, line1, line2]
    corners = [pt0, pt1, center]
    return CircleWedge{typeof(x0)}(boundary,corners,center,id)
end