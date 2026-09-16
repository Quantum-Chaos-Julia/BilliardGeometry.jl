# Limaçon billiard r(φ) = R(1 + a*cos(φ)); x-axis reflection symmetry.

struct LimaconBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::SimpleDomain{T}
    symmetries::SymmetryRegistry
end

function LimaconBilliard(a::T; R::T = one(T)) where T<:Real
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