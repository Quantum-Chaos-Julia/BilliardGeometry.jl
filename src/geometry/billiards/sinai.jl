# Sinai billiard: square [-a,a] × [-a,a] with a circular obstacle of radius R.
#
# Both connected boundary components use the canonical symmetry-compatible
# parametrization starting on the positive x-axis. The outer square starts at
# (a,0) and is traversed counterclockwise; the inner circle starts at (R,0).
# This convention makes the exact index maps used by symmetry_index_orbits
# applicable independently to both components.
#
# The square is the positively oriented outer boundary and the circle is marked
# with orientation=-1 as a hole. CompositeBIMSolver samples both components in
# their canonical ordering and subsequently flips only the oriented differential
# data of the hole, leaving its symmetry indexing unchanged.

struct SinaiBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::MultiplyConnectedDomain{T}
    symmetries::SymmetryRegistry
end

function SinaiBilliard(a::T, R::T; center=SVector{2,T}(zero(T),zero(T))) where T<:Real
    a > zero(T) || throw(ArgumentError("a must be positive; received $a"))
    R > zero(T) || throw(ArgumentError("R must be positive; received $R"))
    c = SVector{2,T}(center)
    abs(c[1])+R < a && abs(c[2])+R < a || throw(ArgumentError("Circular obstacle must lie strictly inside the square; received a=$a, R=$R, center=$c"))
    bc = SpecularReflection(); z = zero(T)
    p0 = SVector{2,T}(a,z)
    p1 = SVector{2,T}(a,a); p2 = SVector{2,T}(-a,a)
    p3 = SVector{2,T}(-a,-a); p4 = SVector{2,T}(a,-a)
    e1 = LineSegment(p0,p1; bc=bc, domain_id=1, segment_id=1)
    e2 = LineSegment(p1,p2; bc=bc, domain_id=1, segment_id=2)
    e3 = LineSegment(p2,p3; bc=bc, domain_id=1, segment_id=3)
    e4 = LineSegment(p3,p4; bc=bc, domain_id=1, segment_id=4)
    e5 = LineSegment(p4,p0; bc=bc, domain_id=1, segment_id=5)
    inner = CircleSegment(R,T(2*pi),z,c; bc=bc, orientation=-1, domain_id=2, segment_id=1)
    outer = AbsCurve[e1,e2,e3,e4,e5]
    holes = Vector{AbsCurve}[[inner]]
    vertices = SVector{2,T}[p0,p1,p2,p3,p4,c+SVector{2,T}(R,z)]
    fundamental_domain = MultiplyConnectedDomain(outer,holes,vertices,1)
    symmetries = register_symmetries()
    return SinaiBilliard{T}(fundamental_domain,symmetries)
end