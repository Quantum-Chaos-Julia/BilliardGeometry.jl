"""
    MushroomBilliard{T}<:BilliardGeometry.AbsBilliard

Mushroom billiard with semicircular cap radius `R`, centered stem half-width
`w` and stem height `h`.

The vertical y-axis is the reflection-symmetry axis.
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
    c=SVector{2,T}(center);bc=BilliardGeometry.SpecularReflection()
    pR=c+SVector{2,T}(R,0);pT=c+SVector{2,T}(0,R);pL=c+SVector{2,T}(-R,0)
    sL=c+SVector{2,T}(-w,0);sR=c+SVector{2,T}(w,0)
    bL=c+SVector{2,T}(-w,-h);bM=c+SVector{2,T}(0,-h);bR=c+SVector{2,T}(w,-h)
    cap=BilliardGeometry.CircleSegment(R,T(pi),zero(T),c;bc=bc,domain_id=1,segment_id=1)
    shelfL=BilliardGeometry.LineSegment(pL,sL;bc=bc,domain_id=1,segment_id=2)
    stemL=BilliardGeometry.LineSegment(sL,bL;bc=bc,domain_id=1,segment_id=3)
    bottom=BilliardGeometry.LineSegment(bL,bR;bc=bc,domain_id=1,segment_id=4)
    stemR=BilliardGeometry.LineSegment(bR,sR;bc=bc,domain_id=1,segment_id=5)
    shelfR=BilliardGeometry.LineSegment(sR,pR;bc=bc,domain_id=1,segment_id=6)
    full_boundary=[BilliardGeometry.AbsCurve[cap,shelfL,stemL,bottom,stemR,shelfR]]
    capR=BilliardGeometry.CircleSegment(R,T(pi/2),zero(T),c;bc=bc,domain_id=1,segment_id=1)
    ywall=BilliardGeometry.LineSegment(pT,bM;bc=BilliardGeometry.ReflectionSymmetry(BilliardGeometry.YAxisReflection(),2),domain_id=1,segment_id=2)
    bottomR=BilliardGeometry.LineSegment(bM,bR;bc=bc,domain_id=1,segment_id=3)
    stemRq=BilliardGeometry.LineSegment(bR,sR;bc=bc,domain_id=1,segment_id=4)
    shelfRq=BilliardGeometry.LineSegment(sR,pR;bc=bc,domain_id=1,segment_id=5)
    fundamental_boundary=BilliardGeometry.AbsCurve[capR,ywall,bottomR,stemRq,shelfRq]
    vertices=SVector{2,T}[pR,pT,bM,bR,sR]
    fundamental_domain=BilliardGeometry.SimpleDomain{T}(fundamental_boundary,vertices,1)
    symmetries=BilliardGeometry.AbsSymmetry[BilliardGeometry.YAxisReflection()]
    return MushroomBilliard{T}(fundamental_domain,full_boundary,symmetries)
end