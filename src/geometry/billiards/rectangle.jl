#rectangle billiard 2a x 2b, D2 symmetry (first-quadrant fundamental domain)
struct RectangleBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::SimpleDomain{T}
    symmetries::Vector{AbsSymmetry}
end

function RectangleBilliard(a::T=1.0, b::T=1.0; center=SVector{2,T}(zero(T),zero(T))) where T<:Real
    c = SVector{2,T}(center)
    (iszero(c[1]) && iszero(c[2])) || throw(ArgumentError("D2 symmetry requires center == (0,0); received center=$c"))
    cx, cy = c
    bc = SpecularReflection()
    q0 = c
    q1 = SVector{2,T}(cx+a,cy)
    p3 = SVector{2,T}(cx+a,cy+b)
    q3 = SVector{2,T}(cx,cy+b)
    right = LineSegment(q1, p3; bc=bc, domain_id=1, segment_id=1)
    top = LineSegment(p3, q3; bc=bc, domain_id=1, segment_id=2)
    ywall = LineSegment(q3, q0; bc=ReflectionSymmetry(YAxisReflection(),4), domain_id=1, segment_id=3)
    xwall = LineSegment(q0, q1; bc=ReflectionSymmetry(XAxisReflection(),4), domain_id=1, segment_id=4)
    fundamental_boundary = AbsCurve[right, top, ywall, xwall]
    vertices = SVector{2,T}[q1, p3, q3, q0]
    fundamental_domain = SimpleDomain{T}(fundamental_boundary, vertices, 1)
    symmetries = AbsSymmetry[D2_symmetry...]
    return RectangleBilliard{T}(fundamental_domain, symmetries)
end
