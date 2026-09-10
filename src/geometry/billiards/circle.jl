#circular billiard, D2 symmetry (quadrant fundamental domain)
struct CircleBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::SimpleDomain{T}
    symmetries::Vector{AbsSymmetry}
end

function CircleBilliard(R::T=1.0; center=SVector{2,T}(zero(T),zero(T))) where T<:Real
    c = SVector{2,T}(center)
    (iszero(c[1]) && iszero(c[2])) || throw(ArgumentError("D2 symmetry requires center == (0,0); received center=$c"))
    cx, cy = c
    bc = SpecularReflection()
    pR = SVector{2,T}(cx+R,cy)
    pT = SVector{2,T}(cx,cy+R)
    arc = CircleSegment(R, T(pi/2), zero(T), c; bc=bc, domain_id=1, segment_id=1)
    ywall = LineSegment(pT, c; bc=ReflectionSymmetry(YAxisReflection(),4), domain_id=1, segment_id=2)
    xwall = LineSegment(c, pR; bc=ReflectionSymmetry(XAxisReflection(),4), domain_id=1, segment_id=3)
    fundamental_boundary = AbsCurve[arc, ywall, xwall]
    vertices = SVector{2,T}[pR, pT, c]
    fundamental_domain = SimpleDomain{T}(fundamental_boundary, vertices, 1)
    symmetries = AbsSymmetry[D2_symmetry...]
    return CircleBilliard{T}(fundamental_domain, symmetries)
end
