struct StadiumBilliard{T}<:BilliardGeometry.AbsBilliard
    fundamental_domain::BilliardGeometry.SimpleDomain{T}
    full_boundary::Vector{Vector{BilliardGeometry.AbsCurve}}
    symmetries::Vector{BilliardGeometry.AbsSymmetry}
end

function StadiumBilliard(R::T,a::T;center=SVector{2,T}(zero(T),zero(T))) where {T<:Real}
    R>0||throw(ArgumentError("Require R>0; received R=$R"))
    a>0||throw(ArgumentError("Require a>0; received a=$a"))
    c=SVector{2,T}(center);bc=SpecularReflection();ign=QuantumSolverIgnore()
    cL=c+SVector{2,T}(-a,0);cR=c+SVector{2,T}(a,0)
    pR=c+SVector{2,T}(a+R,0);tR=c+SVector{2,T}(a,R);tM=c+SVector{2,T}(0,R)
    tL=c+SVector{2,T}(-a,R);pL=c+SVector{2,T}(-a-R,0)
    bL=c+SVector{2,T}(-a,-R);bM=c+SVector{2,T}(0,-R);bR=c+SVector{2,T}(a,-R)
    rightT=CircleSegment(R,T(pi/2),zero(T),cR;bc=bc,domain_id=1,segment_id=1)
    topR=LineSegment(tR,tM;bc=bc,domain_id=1,segment_id=2)
    topL=LineSegment(tM,tL;bc=bc,domain_id=1,segment_id=3)
    leftT=CircleSegment(R,T(pi/2),T(pi/2),cL;bc=bc,domain_id=1,segment_id=4)
    leftB=CircleSegment(R,T(pi/2),T(pi),cL;bc=bc,domain_id=1,segment_id=5)
    bottomL=LineSegment(bL,bM;bc=bc,domain_id=1,segment_id=6)
    bottomR=LineSegment(bM,bR;bc=bc,domain_id=1,segment_id=7)
    rightB=CircleSegment(R,T(pi/2),T(3pi/2),cR;bc=bc,domain_id=1,segment_id=8)
    full_boundary=[AbsCurve[rightT,topR,topL,leftT,leftB,bottomL,bottomR,rightB]]
    xwall=LineSegment(c,pR;bc=ign,domain_id=1,segment_id=1)
    arc=CircleSegment(R,T(pi/2),zero(T),cR;bc=bc,domain_id=1,segment_id=2)
    top=LineSegment(tR,tM;bc=bc,domain_id=1,segment_id=3)
    ywall=LineSegment(tM,c;bc=ign,domain_id=1,segment_id=4)
    fundamental_boundary=AbsCurve[xwall,arc,top,ywall]
    vertices=SVector{2,T}[c,pR,tR,tM]
    fundamental_domain=SimpleDomain{T}(fundamental_boundary,vertices,1)
    symmetries=AbsSymmetry[XAxisReflection(),YAxisReflection(),XYAxisReflection()]
    return StadiumBilliard{T}(fundamental_domain,full_boundary,symmetries)
end