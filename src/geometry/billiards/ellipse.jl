#ellipse billiard x^2/a^2 + y^2/b^2 = 1, D2 symmetry (quadrant fundamental domain)
struct EllipseBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::SimpleDomain{T}
    symmetries::Vector{AbsSymmetry}
end

function EllipseBilliard(a::T, b::T; center=SVector{2,T}(zero(T),zero(T))) where T<:Real
    a > zero(T) || throw(ArgumentError("a must be positive; received a=$a"))
    b > zero(T) || throw(ArgumentError("b must be positive; received b=$b"))
    c = SVector{2,T}(center)
    (iszero(c[1]) && iszero(c[2])) || throw(ArgumentError("D2 symmetry requires center == (0,0); received center=$c"))
    bc = SpecularReflection()
    r(phi) = a*b/sqrt((b*cos(phi))^2 + (a*sin(phi))^2)
    arc = PolarSegment(r; R=max(a,b), arc_angle=T(pi/2), shift_angle=zero(T), center=c, bc=bc, domain_id=1, segment_id=1)
    p0 = curve(arc, zero(T))
    p1 = curve(arc, one(T))
    ywall = LineSegment(p1, c; bc=ReflectionSymmetry(YAxisReflection(),4), domain_id=1, segment_id=2)
    xwall = LineSegment(c, p0; bc=ReflectionSymmetry(XAxisReflection(),4), domain_id=1, segment_id=3)
    fundamental_boundary = AbsCurve[arc, ywall, xwall]
    vertices = SVector{2,T}[p0, p1, c]
    fundamental_domain = SimpleDomain{T}(fundamental_boundary, vertices, 1)
    symmetries = AbsSymmetry[D2_symmetry...]
    return EllipseBilliard{T}(fundamental_domain, symmetries)
end
