#circle segment
function circle_eq(R::T, arc_angle::T, shift_angle::T, center::SVector{2,T}, t) where T<:Real
    return SVector(R*cos(arc_angle*t+shift_angle) + center[1], R*sin(arc_angle*t+shift_angle)+center[2])
end

circle_domain(R, center, x, y) = @. hypot(y-center[2],x-center[1]) - R

struct CircleSegment{T}  <: AbsCurve where T<:Real
    radius::T
    arc_angle::T
    shift_angle::T
    center::SVector{2,T}
    orientation::Int64
    length::T
end

function CircleSegment(R, arc_angle, shift_angle, x0, y0; orientation = 1)
    center = SVector(x0,y0)
    L = R*arc_angle 
    return CircleSegment(R,arc_angle,shift_angle,center,orientation,L)
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

# arc length
function arc_length(circle::L, pt::SVector{2,T}) where {L<:CircleSegment, T<:Real}
    center = circle.center
    angle = atan(pt[2]-center[2], pt[1]-center[1]) - circle.shift_angle
    return circle.radius*angle
end

#=
function arc_length(circle::L, pts::AbstractArray) where {L<:CircleSegment}
    return collect(arc_length(circle, pt) for pt in pts)
end
=#