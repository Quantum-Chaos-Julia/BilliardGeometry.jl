"""
    LimaconBilliard{T} <: AbsBilliard

Limaçon billiard `r(φ) = R + a*cos(φ)`, symmetric under reflection about the
`x` axis. The fundamental domain is the upper half, bounded by the arc
`0 ≤ φ ≤ π` and a [`SymmetryWall`](@ref) on the `x` axis.

## Attributes
* `fundamental_domain::SimpleDomain{T}`: The upper half.
* `symmetries::SymmetryRegistry`: The single [`XAxisReflection`](@ref).
"""
struct LimaconBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::SimpleDomain{T}
    symmetries::SymmetryRegistry
end

"""
    LimaconBilliard(a; R = 1) → billiard::LimaconBilliard

Construct a [`LimaconBilliard`](@ref) with `r(φ) = R + a*cos(φ)`.

## Arguments
* `a`: Deformation parameter, `|a| ≤ 1` as required for a simple curve.

## Keyword arguments
* `R = 1` : Constant part of the radius.

## Throws
`ArgumentError` if `|a| > 1`.
"""
function LimaconBilliard(a::T; R::T = one(T)) where T<:Real
    abs(a) <= one(T) || throw(ArgumentError("Limaçon requires |a|<=1 to remain a simple (non-self-intersecting) curve; received a=$a"))
    center = SVector{2,T}(zero(T), zero(T)); coef = SVector{2,T}(zero(T), a)
    bc = SpecularReflection(); symmetries = register_symmetries(XAxisReflection())

    arc = FourierCoeffPolarSegment(coef; R = R, arc_angle = T(pi), shift_angle = zero(T),
        center = center, bc = bc, domain_id = 1, segment_id = 1)

    p0 = curve(arc, zero(T)); p1 = curve(arc, one(T))
    wall = LineSegment(p1, p0; bc = SymmetryWall(1, 2), domain_id = 1, segment_id = 2)

    boundary = AbsCurve[arc, wall]; vertices = SVector{2,T}[p0, p1]
    domain = SimpleDomain{T}(boundary, vertices, 1)
    return LimaconBilliard{T}(domain, symmetries)
end
