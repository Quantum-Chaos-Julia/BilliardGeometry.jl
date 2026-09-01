"""
    SinaiBilliard{T}<:BilliardGeometry.AbsBilliard

Sinai billiard formed by a centered rectangular outer wall and a centered
circular obstacle.

The rectangular component begins at its positive x-axis midpoint and is
traversed counterclockwise. Its right side is split into two smooth pieces so
that the complete periodic component obeys the canonical exact symmetry-index
convention. The circular obstacle uses zero angular shift and therefore has the
same periodic origin.

The fundamental domain is the first-quadrant Sinai sector bounded by the
physical outer wall, the physical quarter-circle obstacle, and the two
coordinate-axis symmetry walls.
"""
struct SinaiBilliard{T}<:BilliardGeometry.AbsBilliard
    fundamental_domain::BilliardGeometry.SimpleDomain{T}
    full_boundary::Vector{Vector{BilliardGeometry.AbsCurve}}
    symmetries::Vector{BilliardGeometry.AbsSymmetry}
end

"""
    SinaiBilliard(a::T,b::T,R_inner::T;center=SVector{2,T}(zero(T),zero(T))) where {T<:Real}

Construct a centered Sinai billiard with outer dimensions `2a×2b` and a
centered circular obstacle of radius `R_inner`.

Both connected components use the canonical symmetry-compatible periodic
origin on the positive x-axis.
"""
function SinaiBilliard(a::T,b::T,R_inner::T;center=SVector{2,T}(zero(T),zero(T))) where {T<:Real}
    c=SVector{2,T}(center)
    iszero(c[1])&&iszero(c[2])||throw(ArgumentError("D2 symmetry currently requires center == (0,0); received center=$c"))
    zero(T)<R_inner<min(a,b)||throw(ArgumentError("Require 0<R_inner<min(a,b); received R_inner=$R_inner"))
    p0=c+SVector{2,T}(a,zero(T))
    p1=c+SVector{2,T}(a,b)
    p2=c+SVector{2,T}(-a,b)
    p3=c+SVector{2,T}(-a,-b)
    p4=c+SVector{2,T}(a,-b)
    q0=c+SVector{2,T}(R_inner,zero(T))
    q1=c+SVector{2,T}(zero(T),R_inner)
    q2=c+SVector{2,T}(zero(T),b)
    bc=BilliardGeometry.SpecularReflection()
    right_upper=BilliardGeometry.LineSegment(p0,p1;bc=bc,domain_id=1,segment_id=1)
    top=BilliardGeometry.LineSegment(p1,p2;bc=bc,domain_id=1,segment_id=2)
    left=BilliardGeometry.LineSegment(p2,p3;bc=bc,domain_id=1,segment_id=3)
    bottom=BilliardGeometry.LineSegment(p3,p4;bc=bc,domain_id=1,segment_id=4)
    right_lower=BilliardGeometry.LineSegment(p4,p0;bc=bc,domain_id=1,segment_id=5)
    inner=BilliardGeometry.CircleSegment(R_inner,T(2pi),zero(T),c;bc=bc,domain_id=1,segment_id=6)
    full_boundary=[
        BilliardGeometry.AbsCurve[right_upper,top,left,bottom,right_lower],
        BilliardGeometry.AbsCurve[inner]
    ]
    right=BilliardGeometry.LineSegment(p0,p1;bc=bc,domain_id=1,segment_id=1)
    top_q=BilliardGeometry.LineSegment(p1,q2;bc=bc,domain_id=1,segment_id=2)
    ywall=BilliardGeometry.LineSegment(q2,q1;bc=BilliardGeometry.ReflectionSymmetry(BilliardGeometry.YAxisReflection(),6),domain_id=1,segment_id=3)
    inner_q=BilliardGeometry.CircleSegment(R_inner,T(pi/2),zero(T),c;bc=bc,orientation=-1,domain_id=1,segment_id=4) # clockwise so the obstacle-facing boundary has the correct fundamental-domain orientation
    xwall=BilliardGeometry.LineSegment(q0,p0;bc=BilliardGeometry.ReflectionSymmetry(BilliardGeometry.XAxisReflection(),6),domain_id=1,segment_id=5)
    fundamental_boundary=BilliardGeometry.AbsCurve[right,top_q,ywall,inner_q,xwall]
    vertices=SVector{2,T}[p0,p1,q2,q1,q0]
    fundamental_domain=BilliardGeometry.SimpleDomain{T}(fundamental_boundary,vertices,1)
    symmetries=BilliardGeometry.AbsSymmetry[BilliardGeometry.D2_symmetry...]
    return SinaiBilliard{T}(fundamental_domain,full_boundary,symmetries)
end