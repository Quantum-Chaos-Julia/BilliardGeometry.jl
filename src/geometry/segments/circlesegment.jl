#circle segment
function circle_eq(R::T, arc_angle::T, shift_angle::T, center::SVector{2,T}, t) where T<:Real
    return SVector(R*cos(arc_angle*t+shift_angle) + center[1], R*sin(arc_angle*t+shift_angle)+center[2])
end

circle_domain(R, center, x, y) = @. hypot(y-center[2],x-center[1]) - R

"""
    CircleSegment{T,BC} <: AbsCurve{BC}

Circular arc `r(t) = center + radius*(cos φ, sin φ)` with `φ = shift_angle +
arc_angle*t`, `t ∈ [0,1]`.

## Description
A negative `arc_angle` traverses the arc clockwise. [`domain_fun`](@ref) is
`|pt - center| - radius` times `orientation`, so with `orientation = 1` the
inside is the disk and with `orientation = -1` the exterior of the circle (an
obstacle).

## Attributes
* `radius::T`: Radius.
* `arc_angle::T`: Signed angle swept by the arc, `2π` for a full circle.
* `shift_angle::T`: Polar angle of the start point.
* `center::SVector{2,T}`: Center of the circle.
* `orientation::Int64`: `+1` or `-1`, selects the inside.
* `length::T`: Arc length `radius*|arc_angle|`.
* `bc::BC`: Boundary condition, see [`AbsBoundaryCondition`](@ref).
* `domain_id::Int64`: Id of the domain the segment belongs to.
* `segment_id::Int64`: Id of the segment inside its domain.
"""
struct CircleSegment{T,BC}  <: AbsCurve{BC} where {T<:Real}
    radius::T
    arc_angle::T
    shift_angle::T
    center::SVector{2,T}
    orientation::Int64
    length::T
    bc::BC
    domain_id::Int64
    segment_id::Int64
end

"""
    CircleSegment(R, arc_angle, shift_angle, center; bc = SpecularReflection(), orientation = 1, domain_id = 1, segment_id = 1) → arc::CircleSegment
    CircleSegment(R, arc_angle, shift_angle, x0, y0; bc = SpecularReflection(), orientation = 1, domain_id = 1, segment_id = 1) → arc::CircleSegment

Construct a [`CircleSegment`](@ref) with the center given as a 2-element vector
or as the coordinates `x0`, `y0`. The numeric type is the promotion of all
numeric arguments.

## Arguments
* `R`: Radius.
* `arc_angle`: Signed angle swept by the arc.
* `shift_angle`: Polar angle of the start point.
* `center`, `x0`, `y0`: Center of the circle.

## Keyword arguments
* `bc = SpecularReflection()` : Boundary condition.
* `orientation = 1` : `+1` or `-1`, selects the inside.
* `domain_id = 1` : Id of the domain the segment belongs to.
* `segment_id = 1` : Id of the segment inside its domain.
"""
function CircleSegment(R, arc_angle, shift_angle, center; bc = SpecularReflection(), orientation = 1, domain_id=1, segment_id=1)
    T = promote_type(typeof(R), typeof(arc_angle), typeof(shift_angle), eltype(center))
    center = SVector{2,T}(center)
    L = T(R)*abs(T(arc_angle))
    return CircleSegment(T(R),T(arc_angle),T(shift_angle),center,orientation,L,bc,domain_id,segment_id)
end

function CircleSegment(R, arc_angle, shift_angle, x0, y0; bc = SpecularReflection(), orientation = 1, domain_id=1, segment_id=1)
    T = promote_type(typeof(R), typeof(arc_angle), typeof(shift_angle), typeof(x0), typeof(y0))
    center = SVector{2,T}(x0,y0)
    L = T(R)*abs(T(arc_angle))
    return CircleSegment(T(R),T(arc_angle),T(shift_angle),center,orientation,L,bc,domain_id,segment_id)
end

# returns SVector(x,y)
function curve(circle::L, t) where {L<:CircleSegment}
    return circle_eq(circle.radius, circle.arc_angle, circle.shift_angle, circle.center, t)
end

function curve(circle::L, ts::AbstractArray) where {L<:CircleSegment}
    let R = circle.radius, c = circle.center, a=circle.arc_angle, s=circle.shift_angle 
        return collect(circle_eq(R, a, s, c, t) for t in ts)
    end
end

# returns negative value inside
function domain_fun(circle::L, pt::SVector{2,T}) where {L<:CircleSegment, T<:Real}
    return circle_domain(circle.radius, circle.center, pt[1], pt[2])*circle.orientation
end

function domain_fun(circle::L, pts::AbstractArray) where {L<:CircleSegment}
    let R = circle.radius, c = circle.center, orientation = circle.orientation
    return collect(circle_domain(R, c, pt[1], pt[2])*orientation for pt in pts)
    end
end

# arc length from the start of the arc to the point `pt` on it, measured along the sense of traversal; the angle is
# wrapped into [-ε, 2π-ε) so that points at the start map to 0 (up to rounding) and the end of a full turn maps to 0
function arc_length(circle::L, pt::SVector{2,T}) where {L<:CircleSegment, T<:Real}
    center = circle.center
    ε = sqrt(eps(T))
    direction = circle.arc_angle < 0 ? -one(T) : one(T)
    swept = direction*(atan(pt[2]-center[2], pt[1]-center[1]) - circle.shift_angle)
    return circle.radius*(mod(swept + ε, 2*T(π)) - ε)
end
