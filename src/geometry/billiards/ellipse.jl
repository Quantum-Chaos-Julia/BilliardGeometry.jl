"""
    EllipseBilliard{T}<:BilliardGeometry.AbsBilliard

Ellipse with semi-axes `a` and `b`,

    x²/a²+y²/b²=1.

The full physical boundary is represented by a function-defined
`PolarSegment`. The fundamental domain is the first-quadrant sector.
"""
struct EllipseBilliard{T}<:BilliardGeometry.AbsBilliard
    fundamental_domain::BilliardGeometry.SimpleDomain{T}
    full_boundary::Vector{BilliardGeometry.AbsCurve}
    symmetries::Vector{BilliardGeometry.AbsSymmetry}
end

"""
    EllipseBilliard(a::T,b::T;center=SVector{2,T}(zero(T),zero(T))) where {T<:Real}

Construct an axis-aligned elliptical billiard with semi-axes `a` and `b`.

## Arguments
* `a::T`: Semi-axis along the x direction.
* `b::T`: Semi-axis along the y direction.
* `center::SVector{2,T}`: Center of the ellipse.

## Returns
* `EllipseBilliard{T}`: Constructed elliptical billiard.
"""
function EllipseBilliard(a::T,b::T;center=SVector{2,T}(zero(T),zero(T))) where {T<:Real}
    a>zero(T)||throw(ArgumentError("a must be positive; received a=$a"))
    b>zero(T)||throw(ArgumentError("b must be positive; received b=$b"))
    c=SVector{2,T}(center)
    iszero(c[1])&&iszero(c[2])||throw(ArgumentError("D2 symmetry currently requires center == (0,0); received center=$c"))
    bc=BilliardGeometry.SpecularReflection()
    r(φ)=a*b/sqrt((b*cos(φ))^2+(a*sin(φ))^2)
    full=BilliardGeometry.PolarSegment(r;R=max(a,b),arc_angle=T(2pi),shift_angle=zero(T),center=c,bc=bc,domain_id=1,segment_id=1)
    full_boundary=BilliardGeometry.AbsCurve[full]
    arc=BilliardGeometry.PolarSegment(r;R=max(a,b),arc_angle=T(pi/2),shift_angle=zero(T),center=c,bc=bc,domain_id=1,segment_id=1)
    p0=BilliardGeometry.curve(arc,zero(T))
    p1=BilliardGeometry.curve(arc,one(T))
    wall1=BilliardGeometry.LineSegment(p1,c;bc=BilliardGeometry.ReflectionSymmetry(BilliardGeometry.YAxisReflection(),4),domain_id=1,segment_id=2)
    wall0=BilliardGeometry.LineSegment(c,p0;bc=BilliardGeometry.ReflectionSymmetry(BilliardGeometry.XAxisReflection(),4),domain_id=1,segment_id=3)
    fundamental_boundary=BilliardGeometry.AbsCurve[arc,wall1,wall0]
    vertices=SVector{2,T}[p0,p1,c]
    fundamental_domain=BilliardGeometry.SimpleDomain{T}(fundamental_boundary,vertices,1)
    symmetries=BilliardGeometry.AbsSymmetry[BilliardGeometry.YAxisReflection(),BilliardGeometry.XYAxisReflection(),BilliardGeometry.XAxisReflection()]
    return EllipseBilliard{T}(fundamental_domain,full_boundary,symmetries)
end