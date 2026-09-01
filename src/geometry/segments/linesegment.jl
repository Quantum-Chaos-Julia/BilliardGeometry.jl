@inline line_eq(pt0::SVector{2,T},pt1::SVector{2,T},t) where {T<:Number}=@. (pt1-pt0)*t+pt0
@inline line_domain(x0,y0,x1,y1,x,y)=(y1-y0)*x-(x1-x0)*y+x1*y0-y1*x0

abstract type AbsLine{BC}<:AbsCurve{BC} end

"""
    LineSegment{T,BC}<:AbsLine{BC}

Straight boundary segment from `pt0` to `pt1` parametrized by

    x(t)=pt0+t*(pt1-pt0),    t∈[0,1].

The traversal direction is from `pt0` at `t=0` to `pt1` at `t=1`, while
`orientation` controls the sign of the domain function.
"""
struct LineSegment{T<:Real,BC}<:AbsLine{BC}
    pt0::SVector{2,T}
    pt1::SVector{2,T}
    orientation::Int64
    length::T
    bc::BC
    domain_id::Int64
    segment_id::Int64
end

"""
    LineSegment(pt0,pt1;bc=SpecularReflection(),orientation=1,domain_id=1,segment_id=1)

Construct a straight boundary segment parametrized over `t∈[0,1]`.

## Arguments
* `pt0`: Initial endpoint corresponding to `t=0`.
* `pt1`: Final endpoint corresponding to `t=1`.
* `bc`: Boundary condition.
* `orientation`: Domain orientation, normally `1` or `-1`.
* `domain_id`: Domain identifier.
* `segment_id`: Segment identifier.

## Returns
* `LineSegment`: Constructed line segment.
"""
function LineSegment(pt0,pt1;bc=SpecularReflection(),orientation=1,domain_id=1,segment_id=1)
    T=promote_type(eltype(pt0),eltype(pt1))
    p0=SVector{2,T}(pt0)
    p1=SVector{2,T}(pt1)
    L=norm(p1-p0)
    return LineSegment(p0,p1,Int64(orientation),L,bc,Int64(domain_id),Int64(segment_id))
end

"""
    curve(line::L,t) where {L<:AbsLine}

Evaluate a line segment at parameter `t∈[0,1]`.

## Arguments
* `line::L`: Line segment.
* `t`: Curve parameter in `[0,1]`.

## Returns
* `SVector{2}`: Cartesian boundary point.
"""
@inline function curve(line::L,t) where {L<:AbsLine}
    return line_eq(line.pt0,line.pt1,t)
end

"""
    curve(line::L,ts::AbstractArray) where {L<:AbsLine}

Evaluate a line segment at parameters `t∈[0,1]`.

## Arguments
* `line::L`: Line segment.
* `ts::AbstractArray`: Curve parameters with entries in `[0,1]`.

## Returns
* `Vector`: Cartesian boundary points.
"""
function curve(line::L,ts::AbstractArray) where {L<:AbsLine}
    return [curve(line,t) for t in ts]
end

"""
    domain_fun(line::L,pt::SVector{2,T}) where {L<:AbsLine,T<:Real}

Evaluate the oriented line domain function.

## Arguments
* `line::L`: Line segment.
* `pt::SVector{2,T}`: Cartesian point.

## Returns
* `Real`: Oriented domain-function value, negative on the domain side.
"""
@inline function domain_fun(line::L,pt::SVector{2,T}) where {L<:AbsLine,T<:Real}
    return line_domain(line.pt0[1],line.pt0[2],line.pt1[1],line.pt1[2],pt[1],pt[2])*line.orientation
end

"""
    domain_fun(line::L,pts::AbstractArray) where {L<:AbsLine}

Evaluate the oriented line domain function at multiple points.

## Arguments
* `line::L`: Line segment.
* `pts::AbstractArray`: Cartesian points.

## Returns
* `Vector`: Oriented domain-function values.
"""
function domain_fun(line::L,pts::AbstractArray) where {L<:AbsLine}
    return [domain_fun(line,pt) for pt in pts]
end

"""
    arc_length(line::L,t::T) where {L<:AbsLine,T<:Real}

Compute the physical arc length from parameter `0` to `t∈[0,1]`.

## Arguments
* `line::L`: Line segment.
* `t::T`: Curve parameter in `[0,1]`.

## Returns
* `T`: Arc length from `t=0` to `t`.
"""
@inline function arc_length(line::L,t::T) where {L<:AbsLine,T<:Real}
    return line.length*t
end

"""
    arc_length(line::L,ts::AbstractArray) where {L<:AbsLine}

Compute physical arc lengths at parameters `t∈[0,1]`.

## Arguments
* `line::L`: Line segment.
* `ts::AbstractArray`: Curve parameters with entries in `[0,1]`.

## Returns
* `Vector`: Arc lengths from `t=0` to each parameter.
"""
function arc_length(line::L,ts::AbstractArray) where {L<:AbsLine}
    return [arc_length(line,t) for t in ts]
end

"""
    arc_length(line::L,pt::SVector{2,T}) where {L<:AbsLine,T<:Real}

Compute the physical arc length from `pt0` to a point on the segment.

## Arguments
* `line::L`: Line segment.
* `pt::SVector{2,T}`: Point on the segment corresponding to some `t∈[0,1]`.

## Returns
* `T`: Arc length measured from `pt0`.
"""
@inline function arc_length(line::L,pt::SVector{2,T}) where {L<:AbsLine,T<:Real}
    return norm(pt-line.pt0)
end