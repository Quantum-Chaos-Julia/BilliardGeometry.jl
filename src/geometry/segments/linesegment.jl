
#line segments
line_eq(pt0::SVector{2,T}, pt1::SVector{2,T}, t) where {T<:Number} = @. (pt1 - pt0) * t + pt0
line_domain(x0,y0,x1,y1,x,y) = ((y1-y0)*x-(x1-x0)*y+x1*y0-y1*x0)

abstract type AbsLine{BC} <: AbsCurve{BC} end 

"""
    LineSegment{T,BC} <: AbsCurve{BC}

Straight boundary segment `r(t) = pt0 + t*(pt1 - pt0)`, `t ∈ [0,1]`.

## Description
[`domain_fun`](@ref) is negative to the left of the direction `pt0 → pt1` when
`orientation = 1`, so a counterclockwise traversed outer boundary has its
interior inside. `orientation = -1` selects the other side.

## Attributes
* `pt0::SVector{2,T}`: Start point.
* `pt1::SVector{2,T}`: End point.
* `orientation::Int64`: `+1` or `-1`, selects the inside half-plane.
* `length::T`: Euclidean length `|pt1 - pt0|`.
* `bc::BC`: Boundary condition, see [`AbsBoundaryCondition`](@ref).
* `domain_id::Int64`: Id of the domain the segment belongs to.
* `segment_id::Int64`: Id of the segment inside its domain.
"""
struct LineSegment{T, BC}  <: AbsLine{BC} where {T<:Real}
    pt0::SVector{2,T}
    pt1::SVector{2,T}
    orientation::Int64
    length::T
    bc::BC
    domain_id::Int64
    segment_id::Int64
end

#constructors
"""
    LineSegment(pt0, pt1; bc = SpecularReflection(), orientation = 1, domain_id = 1, segment_id = 1) → line::LineSegment

Construct a [`LineSegment`](@ref) from `pt0` to `pt1`. The numeric type is the
promotion of the element types of the two points and the length is computed
automatically.

## Arguments
* `pt0`, `pt1`: Start and end point, any 2-element vectors.

## Keyword arguments
* `bc = SpecularReflection()` : Boundary condition.
* `orientation = 1` : `+1` or `-1`, selects the inside half-plane.
* `domain_id = 1` : Id of the domain the segment belongs to.
* `segment_id = 1` : Id of the segment inside its domain.
"""
function LineSegment(pt0, pt1; bc = SpecularReflection(), orientation = 1, domain_id=1, segment_id=1) 
    T = promote_type(eltype(pt0), eltype(pt1))
    pt0 = SVector{2,T}(pt0)
    pt1 = SVector{2,T}(pt1)
    x, y = pt1 .- pt0        
    L = hypot(x,y)
    return LineSegment(pt0,pt1,orientation,L,bc, domain_id, segment_id)
end


# returns SVector(x,y)
"""
    curve(crv::AbsCurve, t::Real) → pt::SVector{2}
    curve(crv::AbsCurve, ts::AbstractArray) → pts::Vector{SVector{2}}

Point(s) `r(t)` of the curve at the parameter value(s) `t ∈ [0,1]`. Every
concrete [`AbsCurve`](@ref) implements this method; the parametrization is
linear for a [`LineSegment`](@ref), by angle for a [`CircleSegment`](@ref) and
by the polar angle for polar curves.
"""
function curve(line::L, t) where {L<:AbsLine}
    return line_eq(line.pt0,line.pt1,t)
end
function curve(line::L, ts::AbstractArray) where {L<:AbsLine}
    return collect(line_eq(line.pt0,line.pt1,t) for t in ts)
end

# returns negative value inside
"""
    domain_fun(crv::AbsCurve, pt::SVector{2}) → d
    domain_fun(crv::AbsCurve, pts::AbstractArray) → ds::Vector

Signed defining function of the region bounded by `crv`, multiplied by the
curve's `orientation`: negative inside, positive outside, zero on the curve
(or on its extension). For a [`LineSegment`](@ref) it is the signed distance to
the supporting line times the segment length, for a [`CircleSegment`](@ref)
`|pt - center| - radius` and for polar curves `r - r(φ)`. Used by
[`is_inside`](@ref) and [`domain_gradient_vector`](@ref).
"""
function domain_fun(line::L, pt::SVector{2,T}) where {L<:AbsLine, T<:Real}
    let pt0 = line.pt0 
        pt1 = line.pt1
        orientation = line.orientation
    return line_domain(pt0[1],pt0[2],pt1[1],pt1[2],pt[1],pt[2])*orientation
    end
end

function domain_fun(line::L, pts::AbstractArray) where {L<:AbsLine}
    let pt0 = line.pt0 
        pt1 = line.pt1
        orientation = line.orientation
    return collect(line_domain(pt0[1],pt0[2],pt1[1],pt1[2],pt[1],pt[2])*orientation for pt in pts)
    end
end

# arc length
function arc_length(line::L, pt::SVector{2,T}) where {L<:AbsLine, T<:Real}
    r0 = line.pt0
    x, y = pt .- r0
    return hypot(x, y)
end


