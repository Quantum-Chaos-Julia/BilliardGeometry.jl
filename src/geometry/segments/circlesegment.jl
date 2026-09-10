"""
    CircleSegment{T,BC}<:AbsCurve{BC}

Circular boundary segment of radius `radius` parametrized by

    φ(t)=shift_angle+arc_angle*t,    t∈[0,1].

A negative `arc_angle` reverses the traversal direction while `orientation`
controls the sign of the domain function.
"""
struct CircleSegment{T<:Real,BC}<:AbsCurve{BC}
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

@inline function circle_eq(R::T,arc_angle::T,shift_angle::T,center::SVector{2,T},t) where {T<:Real}
    φ=arc_angle*t+shift_angle
    return SVector(R*cos(φ)+center[1],R*sin(φ)+center[2])
end

@inline circle_domain(R,center,x,y)=@. hypot(y-center[2],x-center[1])-R

"""
    CircleSegment(R,arc_angle,shift_angle,center;bc=SpecularReflection(),orientation=1,domain_id=1,segment_id=1)

Construct a circular boundary segment parametrized over `t∈[0,1]`.

## Arguments
* `R`: Circle radius.
* `arc_angle`: Signed angular extent traversed as `t` goes from `0` to `1`.
* `shift_angle`: Polar angle corresponding to `t=0`.
* `center`: Circle center.
* `bc`: Boundary condition.
* `orientation`: Domain orientation, normally `1` or `-1`.
* `domain_id`: Domain identifier.
* `segment_id`: Segment identifier.

## Returns
* `CircleSegment`: Constructed circular segment.
"""
function CircleSegment(R,arc_angle,shift_angle,center;bc=SpecularReflection(),orientation=1,domain_id=1,segment_id=1)
    T=promote_type(typeof(R),typeof(arc_angle),typeof(shift_angle),eltype(center))
    c=SVector{2,T}(center)
    return CircleSegment(T(R),T(arc_angle),T(shift_angle),c,Int64(orientation),abs(T(R)*T(arc_angle)),bc,Int64(domain_id),Int64(segment_id))
end

"""
    CircleSegment(R,arc_angle,shift_angle,x0,y0;bc=SpecularReflection(),orientation=1,domain_id=1,segment_id=1)

Construct a circular boundary segment centered at `(x0,y0)` and parametrized
over `t∈[0,1]`.

## Arguments
* `R`: Circle radius.
* `arc_angle`: Signed angular extent traversed as `t` goes from `0` to `1`.
* `shift_angle`: Polar angle corresponding to `t=0`.
* `x0`: Center x-coordinate.
* `y0`: Center y-coordinate.
* `bc`: Boundary condition.
* `orientation`: Domain orientation, normally `1` or `-1`.
* `domain_id`: Domain identifier.
* `segment_id`: Segment identifier.

## Returns
* `CircleSegment`: Constructed circular segment.
"""
function CircleSegment(R,arc_angle,shift_angle,x0,y0;bc=SpecularReflection(),orientation=1,domain_id=1,segment_id=1)
    return CircleSegment(R,arc_angle,shift_angle,SVector(x0,y0);bc=bc,orientation=orientation,domain_id=domain_id,segment_id=segment_id)
end

"""
    curve(circle::L,t) where {L<:CircleSegment}

Evaluate a circular segment at parameter `t∈[0,1]`.

## Arguments
* `circle::L`: Circular segment.
* `t`: Curve parameter in `[0,1]`.

## Returns
* `SVector{2}`: Cartesian boundary point.
"""
@inline function curve(circle::L,t) where {L<:CircleSegment}
    return circle_eq(circle.radius,circle.arc_angle,circle.shift_angle,circle.center,t)
end

"""
    curve(circle::L,ts::AbstractArray) where {L<:CircleSegment}

Evaluate a circular segment at parameters `t∈[0,1]`.

## Arguments
* `circle::L`: Circular segment.
* `ts::AbstractArray`: Curve parameters with entries in `[0,1]`.

## Returns
* `Vector`: Cartesian boundary points.
"""
function curve(circle::L,ts::AbstractArray) where {L<:CircleSegment}
    return [curve(circle,t) for t in ts]
end

"""
    domain_fun(circle::L,pt::SVector{2,T}) where {L<:CircleSegment,T<:Real}

Evaluate the oriented circle domain function.

## Arguments
* `circle::L`: Circular segment.
* `pt::SVector{2,T}`: Cartesian point.

## Returns
* `Real`: Oriented domain-function value, negative on the domain side.
"""
@inline function domain_fun(circle::L,pt::SVector{2,T}) where {L<:CircleSegment,T<:Real}
    return circle_domain(circle.radius,circle.center,pt[1],pt[2])*circle.orientation
end

"""
    domain_fun(circle::L,pts::AbstractArray) where {L<:CircleSegment}

Evaluate the oriented circle domain function at multiple points.

## Arguments
* `circle::L`: Circular segment.
* `pts::AbstractArray`: Cartesian points.

## Returns
* `Vector`: Oriented domain-function values.
"""
function domain_fun(circle::L,pts::AbstractArray) where {L<:CircleSegment}
    return [domain_fun(circle,pt) for pt in pts]
end

"""
    arc_length(circle::L,t::T) where {L<:CircleSegment,T<:Real}

Compute the physical arc length from parameter `0` to `t∈[0,1]`.

## Arguments
* `circle::L`: Circular segment.
* `t::T`: Curve parameter in `[0,1]`.

## Returns
* `T`: Arc length from `t=0` to `t`.
"""
@inline function arc_length(circle::L,t::T) where {L<:CircleSegment,T<:Real}
    return abs(circle.radius*circle.arc_angle)*t
end

"""
    arc_length(circle::L,ts::AbstractArray) where {L<:CircleSegment}

Compute physical arc lengths at parameters `t∈[0,1]`.

## Arguments
* `circle::L`: Circular segment.
* `ts::AbstractArray`: Curve parameters with entries in `[0,1]`.

## Returns
* `Vector`: Arc lengths from `t=0` to each parameter.
"""
function arc_length(circle::L,ts::AbstractArray) where {L<:CircleSegment}
    return [arc_length(circle,t) for t in ts]
end

"""
    arc_length(circle::L,pt::SVector{2,T}) where {L<:CircleSegment,T<:Real}

Compute the physical arc length from the segment start to a point on the
circular segment.

## Arguments
* `circle::L`: Circular segment.
* `pt::SVector{2,T}`: Point on the segment corresponding to some `t∈[0,1]`.

## Returns
* `T`: Arc length measured in the traversal direction from `t=0`.
"""
function arc_length(circle::L,pt::SVector{2,T}) where {L<:CircleSegment,T<:Real}
    φ=atan(pt[2]-circle.center[2],pt[1]-circle.center[1])
    Δφ=circle.arc_angle>=zero(circle.arc_angle) ? mod2pi(φ-circle.shift_angle) : mod2pi(circle.shift_angle-φ)
    return circle.radius*Δφ
end