#circle segment
function circle_eq(R::T, arc_angle::T, shift_angle::T, center::SVector{2,T}, t) where T<:Real
    return SVector(R*cos(arc_angle*t+shift_angle) + center[1], R*sin(arc_angle*t+shift_angle)+center[2])
end

circle_domain(R, center, x, y) = @. hypot(y-center[2],x-center[1]) - R

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
