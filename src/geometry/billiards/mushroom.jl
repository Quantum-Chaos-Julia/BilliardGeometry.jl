"""
    MushroomBilliard{T}<:BilliardGeometry.AbsBilliard

Mushroom billiard with semicircular cap radius `R`, centered stem half-width `w`
and stem height `h`. The vertical y-axis is the reflection-symmetry axis.

The full boundary is split into exact Y-reflection partners:
    1 ↔ 2, 3 ↔ 8, 4 ↔ 7, 5 ↔ 6.
"""
struct MushroomBilliard{T}<:BilliardGeometry.AbsBilliard
    fundamental_domain::BilliardGeometry.SimpleDomain{T}
    full_boundary::Vector{Vector{BilliardGeometry.AbsCurve}}
    symmetries::Vector{BilliardGeometry.AbsSymmetry}
end

function MushroomBilliard(R::T,w::T,h::T;center=SVector{2,T}(zero(T),zero(T))) where {T<:Real}
    R>0||throw(ArgumentError("Require R>0; received R=$R"))
    0<w<R||throw(ArgumentError("Require 0<w<R; received w=$w, R=$R"))
    h>0||throw(ArgumentError("Require h>0; received h=$h"))
    c=SVector{2,T}(center)
    bc=BilliardGeometry.SpecularReflection()
    pR=c+SVector{2,T}(R,0);pT=c+SVector{2,T}(0,R);pL=c+SVector{2,T}(-R,0)
    sL=c+SVector{2,T}(-w,0);sR=c+SVector{2,T}(w,0)
    bL=c+SVector{2,T}(-w,-h);bM=c+SVector{2,T}(0,-h);bR=c+SVector{2,T}(w,-h)

    capR=BilliardGeometry.CircleSegment(R,T(pi/2),zero(T),c;bc=bc,domain_id=1,segment_id=1)
    capL=BilliardGeometry.CircleSegment(R,T(pi/2),T(pi/2),c;bc=bc,domain_id=1,segment_id=2)
    shelfL=BilliardGeometry.LineSegment(pL,sL;bc=bc,domain_id=1,segment_id=3)
    stemL=BilliardGeometry.LineSegment(sL,bL;bc=bc,domain_id=1,segment_id=4)
    bottomL=BilliardGeometry.LineSegment(bL,bM;bc=bc,domain_id=1,segment_id=5)
    bottomR=BilliardGeometry.LineSegment(bM,bR;bc=bc,domain_id=1,segment_id=6)
    stemR=BilliardGeometry.LineSegment(bR,sR;bc=bc,domain_id=1,segment_id=7)
    shelfR=BilliardGeometry.LineSegment(sR,pR;bc=bc,domain_id=1,segment_id=8)
    full_boundary=[BilliardGeometry.AbsCurve[capR,capL,shelfL,stemL,bottomL,bottomR,stemR,shelfR]]

    symbc=BilliardGeometry.ReflectionSymmetry(BilliardGeometry.YAxisReflection(),8)
    ywall=BilliardGeometry.LineSegment(pT,bM;bc=symbc,domain_id=1,segment_id=2)
    fundamental_boundary=BilliardGeometry.AbsCurve[capR,ywall,bottomR,stemR,shelfR]
    vertices=SVector{2,T}[pR,pT,bM,bR,sR]
    fundamental_domain=BilliardGeometry.SimpleDomain{T}(fundamental_boundary,vertices,1)
    symmetries=BilliardGeometry.AbsSymmetry[BilliardGeometry.YAxisReflection()]
    return MushroomBilliard{T}(fundamental_domain,full_boundary,symmetries)
end