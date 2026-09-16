# Sinai billiard: square [-a,a] × [-a,a] with a circular obstacle of radius R.
# Multiply connected; no symmetry reduction is applied. The square is the
# positively oriented outer boundary and the circle is a negatively oriented
# inner boundary (hole). Distinct domain_ids allow CompositeBIMSolver to group
# the outer square and inner circle as separate connected boundary components.

struct SinaiBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::MultiplyConnectedDomain{T}
    symmetries::SymmetryRegistry
end

function SinaiBilliard(a::T, R::T; center=SVector{2,T}(zero(T), zero(T))) where T<:Real
    a > zero(T) || throw(ArgumentError("a must be positive; received $a"))
    R > zero(T) || throw(ArgumentError("R must be positive; received $R"))
    c = SVector{2,T}(center)
    abs(c[1]) + R < a && abs(c[2]) + R < a || throw(ArgumentError("Circular obstacle must lie strictly inside the square; received a=$a, R=$R, center=$c"))

    bc = SpecularReflection()
    p1 = SVector{2,T}(-a, -a); p2 = SVector{2,T}(a, -a)
    p3 = SVector{2,T}(a, a); p4 = SVector{2,T}(-a, a)

    e1 = LineSegment(p1, p2; bc=bc, domain_id=1, segment_id=1)
    e2 = LineSegment(p2, p3; bc=bc, domain_id=1, segment_id=2)
    e3 = LineSegment(p3, p4; bc=bc, domain_id=1, segment_id=3)
    e4 = LineSegment(p4, p1; bc=bc, domain_id=1, segment_id=4)
    inner = CircleSegment(R, T(2*pi), zero(T), c; bc=bc, orientation=-1, domain_id=2, segment_id=1)

    outer = AbsCurve[e1, e2, e3, e4]
    holes = Vector{AbsCurve}[[inner]]
    vertices = SVector{2,T}[p1, p2, p3, p4, c + SVector{2,T}(R, zero(T))]

    fundamental_domain = MultiplyConnectedDomain(outer, holes, vertices, 1)
    symmetries = register_symmetries()
    return SinaiBilliard{T}(fundamental_domain, symmetries)
end