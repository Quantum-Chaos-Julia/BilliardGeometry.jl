"""
    CircleBilliard{T}<:BilliardGeometry.AbsBilliard

Circular billiard of radius `R` with D2 reflection symmetry.

The complete physical boundary is one full `CircleSegment`. The fundamental
domain is the first-quadrant quarter disk bounded by one physical quarter-circle
arc and the positive coordinate-axis reflection boundaries.
"""
struct CircleBilliard{T}<:BilliardGeometry.AbsBilliard
    fundamental_domain::BilliardGeometry.SimpleDomain{T}
    full_boundary::Vector{BilliardGeometry.CircleSegment}
    symmetries::Vector{BilliardGeometry.AbsSymmetry}
end

"""
    CircleBilliard(R::T=one(T);center=SVector{2,T}(zero(T),zero(T))) where {T<:Real}

Construct a circular billiard of radius `R` with D2 reflection symmetry.

## Arguments
* `R::T`: Circle radius.

## Keyword Arguments
* `center::SVector{2,T}`: Circle center.

## Returns
* `billiard::CircleBilliard{T}`: Circular billiard containing the fundamental
  domain, complete physical boundary and D2 symmetry descriptors.
"""
function CircleBilliard(R::T=one(T);center=SVector{2,T}(zero(T),zero(T))) where {T<:Real}
    c=SVector{2,T}(center)
    iszero(c[1])&&iszero(c[2])||throw(ArgumentError("D2 symmetry currently requires center == (0,0); received center=$c"))
    cx,cy=c
    bc=BilliardGeometry.SpecularReflection()
    pR=SVector{2,T}(cx+R,cy)
    pT=SVector{2,T}(cx,cy+R)
    arc=BilliardGeometry.CircleSegment(R,T(pi/2),zero(T),c;bc=bc,domain_id=1,segment_id=1)
    ywall=BilliardGeometry.LineSegment(pT,c;bc=BilliardGeometry.ReflectionSymmetry(BilliardGeometry.YAxisReflection(),4),domain_id=1,segment_id=2)
    xwall=BilliardGeometry.LineSegment(c,pR;bc=BilliardGeometry.ReflectionSymmetry(BilliardGeometry.XAxisReflection(),4),domain_id=1,segment_id=3)
    fundamental_boundary=BilliardGeometry.AbsCurve[arc,ywall,xwall]
    vertices=SVector{2,T}[pR,pT,c]
    fundamental_domain=BilliardGeometry.SimpleDomain{T}(fundamental_boundary,vertices,1)
    full_boundary=[BilliardGeometry.CircleSegment(R,T(2pi),zero(T),c;bc=bc,domain_id=1,segment_id=1)]
    symmetries=BilliardGeometry.AbsSymmetry[BilliardGeometry.D2_symmetry...]
    return CircleBilliard{T}(fundamental_domain,full_boundary,symmetries)
end