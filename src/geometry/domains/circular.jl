"""
    CircleWedge{T} <: AbsSimpleDomain

Circular sector bounded by an arc and the two radii joining its end points to
the center.

## Attributes
* `boundary::Vector{AbsCurve}`: The arc (`segment_id = 1`), the radius from the arc end to the center (`2`) and the radius from the center to the arc start (`3`).
* `corners::Vector{SVector{2,T}}`: Arc start, arc end and center.
* `center::SVector{2,T}`: Center of the circle.
* `id::Int64`: Id of the domain, also the `domain_id` of every curve.
"""
struct CircleWedge{T} <: AbsSimpleDomain where T<:Real
    boundary::Vector{AbsCurve}
    corners::Vector{SVector{2,T}}
    center::SVector{2,T}
    id::Int64
end

"""
    CircleWedge(R, arc_angle, shift_angle, x0, y0, id; bcs = [SpecularReflection(), SpecularReflection(), SpecularReflection()]) → wedge::CircleWedge

Construct a [`CircleWedge`](@ref) centered at `(x0, y0)`.

## Arguments
* `R`: Radius.
* `arc_angle`: Angle swept by the arc.
* `shift_angle`: Polar angle of the arc start.
* `x0`, `y0`: Center of the circle.
* `id`: Id of the domain.

## Keyword arguments
* `bcs::Vector{AbsBoundaryCondition}` : Boundary conditions of the arc and of the two radii, in the order of the `boundary` field.
"""
function CircleWedge(R, arc_angle, shift_angle, x0, y0, id; bcs::Vector{AbsBoundaryCondition} = [SpecularReflection(), SpecularReflection(),SpecularReflection()])
    center = SVector(x0,y0)
    circle = CircleSegment(R, arc_angle, shift_angle, x0, y0; bc = bcs[1], domain_id=id, segment_id=1)
    pt0 = curve(circle, 0.0)
    pt1 = curve(circle, 1.0)
    line1 = LineSegment(pt1, center; bc = bcs[2], domain_id=id, segment_id=2) 
    line2 = LineSegment(center, pt0; bc = bcs[3], domain_id=id, segment_id=3) 
    boundary = [circle, line1, line2]
    corners = [pt0, pt1, center]
    return CircleWedge{typeof(x0)}(boundary,corners,center,id)
end