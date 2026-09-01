struct TriangleBilliard{T}<:BilliardGeometry.AbsBilliard
    fundamental_domain::BilliardGeometry.SimpleDomain{T}
    full_boundary::Vector{BilliardGeometry.AbsCurve}
    symmetries::Vector{BilliardGeometry.AbsSymmetry}
end

"""
    TriangleBilliard(p1::SVector{2,T},p2::SVector{2,T},p3::SVector{2,T}) where {T<:Real}

Construct a triangular billiard from three counterclockwise vertices.

No symmetry is assumed. The fundamental domain is therefore the complete triangle.
"""
function TriangleBilliard(p1::SVector{2,T},p2::SVector{2,T},p3::SVector{2,T}) where {T<:Real}
    cross=(p2[1]-p1[1])*(p3[2]-p1[2])-(p2[2]-p1[2])*(p3[1]-p1[1])
    abs(cross)>eps(T)||throw(ArgumentError("triangle vertices are collinear"))
    if cross<0;p2,p3=p3,p2 end
    e1=BilliardGeometry.LineSegment(p1,p2;bc=BilliardGeometry.SpecularReflection(),domain_id=1,segment_id=1)
    e2=BilliardGeometry.LineSegment(p2,p3;bc=BilliardGeometry.SpecularReflection(),domain_id=1,segment_id=2)
    e3=BilliardGeometry.LineSegment(p3,p1;bc=BilliardGeometry.SpecularReflection(),domain_id=1,segment_id=3)
    full_boundary=BilliardGeometry.AbsCurve[e1,e2,e3]
    fundamental_domain=BilliardGeometry.SimpleDomain{T}(copy(full_boundary),SVector{2,T}[p1,p2,p3],1)
    return TriangleBilliard{T}(fundamental_domain,full_boundary,BilliardGeometry.AbsSymmetry[])
end

"""
    IsoscelesTriangleBilliard(a::T=one(T),h::T=one(T);center=SVector{2,T}(zero(T),zero(T))) where {T<:Real}

Construct an isosceles triangular billiard with half-base `a`, height `h`, and
reflection symmetry about the y-axis.

## Arguments
* `a::T=one(T)`: Half-width of the base.
* `h::T=one(T)`: Height of the triangle.

## Keyword Arguments
* `center::SVector{2,T}=SVector{2,T}(zero(T),zero(T))`: Translation of the triangle. Its x-coordinate must be zero for compatibility with the origin-centered y-axis reflection.

## Returns
* `billiard::TriangleBilliard{T}`: Isosceles triangle with full physical boundary, right-half fundamental domain, and `YAxisReflection` symmetry.
"""
function IsoscelesTriangleBilliard(a::T=one(T),h::T=one(T);center=SVector{2,T}(zero(T),zero(T))) where {T<:Real}
    c=SVector{2,T}(center)
    cx,cy=c
    iszero(cx)||throw(ArgumentError("YAxisReflection requires center[1] == 0; received center[1]=$cx"))
    pl=SVector{2,T}(cx-a,cy)
    pr=SVector{2,T}(cx+a,cy)
    pt=SVector{2,T}(cx,cy+h)
    pm=SVector{2,T}(cx,cy)
    # Full physical boundary: canonical start at +x and CCW orientation.
    right=BilliardGeometry.LineSegment(pr,pt;bc=BilliardGeometry.SpecularReflection(),domain_id=1,segment_id=1)
    left=BilliardGeometry.LineSegment(pt,pl;bc=BilliardGeometry.SpecularReflection(),domain_id=1,segment_id=2)
    base=BilliardGeometry.LineSegment(pl,pr;bc=BilliardGeometry.SpecularReflection(),domain_id=1,segment_id=3)
    full_boundary=BilliardGeometry.AbsCurve[right,left,base]
    # Right-half fundamental domain: pr -> pt -> pm -> pr.
    physical=BilliardGeometry.LineSegment(pr,pt;bc=BilliardGeometry.SpecularReflection(),domain_id=1,segment_id=1)
    symwall=BilliardGeometry.LineSegment(pt,pm;bc=BilliardGeometry.ReflectionSymmetry(BilliardGeometry.YAxisReflection(),2),domain_id=1,segment_id=2)
    halfbase=BilliardGeometry.LineSegment(pm,pr;bc=BilliardGeometry.SpecularReflection(),domain_id=1,segment_id=3)
    fundamental_domain=BilliardGeometry.SimpleDomain{T}(BilliardGeometry.AbsCurve[physical,symwall,halfbase],SVector{2,T}[pr,pt,pm],1)
    return TriangleBilliard{T}(fundamental_domain,full_boundary,BilliardGeometry.AbsSymmetry[BilliardGeometry.YAxisReflection()])
end