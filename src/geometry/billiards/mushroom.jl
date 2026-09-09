"""
    MushroomBilliard{T}<:BilliardGeometry.AbsBilliard

Mushroom billiard with cap radius `R`, centered stem half-width `w` and
stem height `h`.

The periodic full-boundary origin is chosen on the positive x-axis so that
reflection across the y-axis acts exactly by the canonical integer map

    q -> mod1(N÷2-q+1,N).

The vertical y-axis is the reflection-symmetry axis.
"""
struct MushroomBilliard{T}<:BilliardGeometry.AbsBilliard
    fundamental_domain::BilliardGeometry.SimpleDomain{T}
    full_boundary::Vector{Vector{BilliardGeometry.AbsCurve}}
    symmetries::Vector{BilliardGeometry.AbsSymmetry}
end

function MushroomBilliard(R::T,w::T,h::T) where {T<:Real}
    R>0||throw(ArgumentError("Require R>0; received R=$R"))
    0<w<R||throw(ArgumentError("Require 0<w<R; received w=$w, R=$R"))
    h>0||throw(ArgumentError("Require h>0; received h=$h"))

    bc=SpecularReflection();ign=QuantumSolverIgnore()
    z=zero(T);cfull=SVector{2,T}(z,z)
    x0=((T(pi)+2)*R-2h)/4
    w<x0<R||throw(ArgumentError("Canonical Y-symmetry origin x0=$x0 must satisfy w<x0<R"))

    p0=cfull+SVector{2,T}(x0,z);pm=cfull+SVector{2,T}(-x0,z)
    pR=cfull+SVector{2,T}(R,z);pL=cfull+SVector{2,T}(-R,z)
    sL=cfull+SVector{2,T}(-w,z);sR=cfull+SVector{2,T}(w,z)
    bL=cfull+SVector{2,T}(-w,-h);bR=cfull+SVector{2,T}(w,-h)

    shelfRO=LineSegment(p0,pR;bc=bc,domain_id=1,segment_id=1)
    cap=CircleSegment(R,T(pi),z,cfull;bc=bc,domain_id=1,segment_id=2)
    shelfLO=LineSegment(pL,pm;bc=bc,domain_id=1,segment_id=3)
    shelfLI=LineSegment(pm,sL;bc=bc,domain_id=1,segment_id=4)
    stemL=LineSegment(sL,bL;bc=bc,domain_id=1,segment_id=5)
    bottom=LineSegment(bL,bR;bc=bc,domain_id=1,segment_id=6)
    stemR=LineSegment(bR,sR;bc=bc,domain_id=1,segment_id=7)
    shelfRI=LineSegment(sR,p0;bc=bc,domain_id=1,segment_id=8)
    full_boundary=[AbsCurve[shelfRO,cap,shelfLO,shelfLI,stemL,bottom,stemR,shelfRI]]

    cfd=SVector{2,T}(-w,z)
    pRf=cfd+SVector{2,T}(R,z);pTf=cfd+SVector{2,T}(z,R)
    sRf=cfd+SVector{2,T}(w,z)
    bMf=cfd+SVector{2,T}(z,-h);bRf=cfd+SVector{2,T}(w,-h)

    bottomR=LineSegment(bMf,bRf;bc=bc,domain_id=1,segment_id=1)
    stemRq=LineSegment(bRf,sRf;bc=ign,domain_id=1,segment_id=2)
    shelfRq=LineSegment(sRf,pRf;bc=ign,domain_id=1,segment_id=3)
    capR=CircleSegment(R,T(pi/2),z,cfd;bc=bc,domain_id=1,segment_id=4)
    ywall=LineSegment(pTf,bMf;bc=bc,domain_id=1,segment_id=5)

    fundamental_boundary=AbsCurve[bottomR,stemRq,shelfRq,capR,ywall]
    vertices=SVector{2,T}[bMf,bRf,sRf,pRf,pTf]
    fundamental_domain=SimpleDomain{T}(fundamental_boundary,vertices,1)
    return MushroomBilliard{T}(fundamental_domain,full_boundary,AbsSymmetry[])
end