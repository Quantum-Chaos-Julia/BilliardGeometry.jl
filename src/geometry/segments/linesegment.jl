
#line segments
line_eq(pt0::SVector{2,T}, pt1::SVector{2,T}, t) where {T<:Number} = @. (pt1 - pt0) * t + pt0
line_domain(x0,y0,x1,y1,x,y) = ((y1-y0)*x-(x1-x0)*y+x1*y0-y1*x0)

abstract type AbsLine{BC} <: AbsCurve{BC} end 

# used to define outer walls of billiard
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
function LineSegment(pt0, pt1; bc = SpecularReflection(), orientation = 1, domain_id=1, segment_id=1) 
    pt0 = SVector{2,eltype(pt0)}(pt0)
    pt1 = SVector{2,eltype(pt1)}(pt1)
    x, y = pt1 .- pt0        
    L = hypot(x,y)
    return LineSegment(pt0,pt1,orientation,L,bc, domain_id, segment_id)
end


# returns SVector(x,y)
function curve(line::L, t) where {L<:AbsLine}
    return line_eq(line.pt0,line.pt1,t)
end
function curve(line::L, ts::AbstractArray) where {L<:AbsLine}
    return collect(line_eq(line.pt0,line.pt1,t) for t in ts)
end

# returns negative value inside
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


