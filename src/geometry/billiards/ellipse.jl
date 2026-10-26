"""
    EllipseBilliard{T} <: AbsBilliard

Elliptic billiard `x²/a² + y²/b² = 1`. The fundamental domain is the first
quadrant, with the `D₂` symmetry group of [`D2_symmetry`](@ref). The boundary is
a [`PolarSegment`](@ref).

## Attributes
* `fundamental_domain::SimpleDomain{T}`: The quarter ellipse.
* `symmetries::SymmetryRegistry`: Reflections `YAxisReflection`, `XYAxisReflection`, `XAxisReflection` with `sym_id` `1, 2, 3`.
"""
struct EllipseBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::SimpleDomain{T}
    symmetries::SymmetryRegistry
end

"""
    EllipseBilliard(a, b; center = (0,0)) → billiard::EllipseBilliard

Construct an [`EllipseBilliard`](@ref) with semi-axes `a` (along `x`) and `b`
(along `y`), both positive. The `center` must be the origin.

## Throws
`ArgumentError` for non-positive semi-axes or a nonzero `center`.
"""
function EllipseBilliard(a::T, b::T; center=SVector{2,T}(zero(T),zero(T))) where T<:Real
    a > zero(T) || throw(ArgumentError("a must be positive; received a=$a"))
    b > zero(T) || throw(ArgumentError("b must be positive; received b=$b"))
    c = SVector{2,T}(center)
    (iszero(c[1]) && iszero(c[2])) || throw(ArgumentError("D2 symmetry requires center == (0,0); received center=$c"))
    bc = SpecularReflection()
    symmetries = D2_symmetry() # sym_id 1=YAxisReflection, 2=XYAxisReflection, 3=XAxisReflection
    r(phi) = a*b/sqrt((b*cos(phi))^2 + (a*sin(phi))^2)
    arc = PolarSegment(r; R=max(a,b), arc_angle=T(pi/2), shift_angle=zero(T), center=c, bc=bc, domain_id=1, segment_id=1)
    p0 = curve(arc, zero(T))
    p1 = curve(arc, one(T))
    ywall = LineSegment(p1, c; bc=SymmetryWall(1,2), domain_id=1, segment_id=2)
    xwall = LineSegment(c, p0; bc=SymmetryWall(3,2), domain_id=1, segment_id=3)
    fundamental_boundary = AbsCurve[arc, ywall, xwall]
    vertices = SVector{2,T}[p0, p1, c]
    fundamental_domain = SimpleDomain{T}(fundamental_boundary, vertices, 1)
    return EllipseBilliard{T}(fundamental_domain, symmetries)
end
