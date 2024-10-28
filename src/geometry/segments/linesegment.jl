
#line segments
line_eq(pt0::SVector{2,T}, pt1::SVector{2,T}, t) where {T<:Number} = @. (pt1 - pt0) * t + pt0
line_domain(x0,y0,x1,y1,x,y) = ((y1-y0)*x-(x1-x0)*y+x1*y0-y1*x0)

abstract type AbsLine <: AbsCurve end


# used to define outer walls of billiard
struct LineSegment{T}  <: AbsLine where T<:Real
    pt0::SVector{2,T}
    pt1::SVector{2,T}
    orientation::Int64
    length::T
end

abstract type AbsVirtLine <: AbsLine end

# used to define subdomain axis of fundamnetal domain_fun for composite findamental domains like in stadium
struct VirtLineSegment{T}  <: AbsVirtLine where T<:Real
    pt0::SVector{2,T}
    pt1::SVector{2,T}
    orientation::Int64
    length::T
    domain_exit_idx::Int64 #encodes the symetry sector the particle transitions into when crossing the line
end

# used to define symetry axis of fundamnetal domain_fun
struct SymLineSegment{T,S}  <: AbsVirtLine where {T<:Real, S<:AbsSymmetry}
    pt0::SVector{2,T}
    pt1::SVector{2,T}
    orientation::Int64
    length::T
    symmetry::S
end

#constructors
function LineSegment(pt0::SVector{2,T}, pt1::SVector{2,T}; orientation = 1) where T<:Real
    x, y = pt1 .- pt0        
    L = hypot(x,y)
    return LineSegment(pt0,pt1,orientation,L)
end

function VirtLineSegment(pt0::SVector{2,T}, pt1::SVector{2,T}, domain_exit_idx; orientation = 1,) where T<:Real
    x, y = pt1 .- pt0        
    L = hypot(x,y)
    return VirtLineSegment(pt0,pt1,orientation,L,domain_exit_idx)
end

function SymLineSegment(pt0::SVector{2,T}, pt1::SVector{2,T}, sym; orientation = 1) where T<:Real
    x, y = pt1 .- pt0        
    L = hypot(x,y)
    return SymLineSegment(pt0,pt1,orientation,L, sym)
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
#=
function arc_length(line::L, pts::AbstractArray) where {L<:AbsLine}
    return collect(arc_length(line, pt) for pt in pts)
end
=#