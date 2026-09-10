struct RectangleBilliard{T}<:BilliardGeometry.AbsBilliard
    fundamental_domain::BilliardGeometry.SimpleDomain{T}
    full_boundary::Vector{BilliardGeometry.AbsCurve}
    symmetries::Vector{BilliardGeometry.AbsSymmetry}
end

"""
    RectangleBilliard(a::T=one(T),b::T=one(T);center=SVector{2,T}(zero(T),zero(T))) where {T<:Real}

Construct an origin-centered rectangle of half-width `a` and half-height `b`
with D2 reflection symmetry. The fundamental domain is the first-quadrant rectangle `[0,a]×[0,b]`.

## Arguments
* `a::T=one(T)`: Half-width in the x direction.
* `b::T=one(T)`: Half-height in the y direction.

## Keyword Arguments
* `center::SVector{2,T}=SVector{2,T}(zero(T),zero(T))`: Rectangle center. The current origin-centered D2 symmetry implementation requires `center=(0,0)`.

## Returns
* `billiard::RectangleBilliard{T}`: Rectangle with canonical full boundary, first-quadrant fundamental domain, and D2 symmetry descriptors.
"""
function RectangleBilliard(a::T=one(T),b::T=one(T);center=SVector{2,T}(zero(T),zero(T))) where {T<:Real}
    c=SVector{2,T}(center)
    iszero(c[1])&&iszero(c[2])||throw(ArgumentError("D2 symmetry currently requires center == (0,0); received center=$c"))
    cx,cy=c
    p1=SVector{2,T}(cx-a,cy-b)
    p2=SVector{2,T}(cx+a,cy-b)
    p3=SVector{2,T}(cx+a,cy+b)
    p4=SVector{2,T}(cx-a,cy+b)
    q0=c
    q1=SVector{2,T}(cx+a,cy)
    q3=SVector{2,T}(cx,cy+b)
    # Full physical boundary: canonical +x start, CCW orientation.
    right_upper=BilliardGeometry.LineSegment(q1,p3;bc=BilliardGeometry.SpecularReflection(),domain_id=1,segment_id=1)
    top_full=BilliardGeometry.LineSegment(p3,p4;bc=BilliardGeometry.SpecularReflection(),domain_id=1,segment_id=2)
    left_full=BilliardGeometry.LineSegment(p4,p1;bc=BilliardGeometry.SpecularReflection(),domain_id=1,segment_id=3)
    bottom_full=BilliardGeometry.LineSegment(p1,p2;bc=BilliardGeometry.SpecularReflection(),domain_id=1,segment_id=4)
    right_lower=BilliardGeometry.LineSegment(p2,q1;bc=BilliardGeometry.SpecularReflection(),domain_id=1,segment_id=5)
    full_boundary=BilliardGeometry.AbsCurve[right_upper,top_full,left_full,bottom_full,right_lower]
    # First-quadrant fundamental domain: q1 -> p3 -> q3 -> q0 -> q1.
    right=BilliardGeometry.LineSegment(q1,p3;bc=BilliardGeometry.SpecularReflection(),domain_id=1,segment_id=1)
    top=BilliardGeometry.LineSegment(p3,q3;bc=BilliardGeometry.SpecularReflection(),domain_id=1,segment_id=2)
    ywall=BilliardGeometry.LineSegment(q3,q0;bc=BilliardGeometry.ReflectionSymmetry(BilliardGeometry.YAxisReflection(),4),domain_id=1,segment_id=3)
    xwall=BilliardGeometry.LineSegment(q0,q1;bc=BilliardGeometry.ReflectionSymmetry(BilliardGeometry.XAxisReflection(),4),domain_id=1,segment_id=4)
    fundamental_domain=BilliardGeometry.SimpleDomain{T}(BilliardGeometry.AbsCurve[right,top,ywall,xwall],SVector{2,T}[q1,p3,q3,q0],1)
    symmetries=BilliardGeometry.AbsSymmetry[BilliardGeometry.D2_symmetry...]
    return RectangleBilliard{T}(fundamental_domain,full_boundary,symmetries)
end