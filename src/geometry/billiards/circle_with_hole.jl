#concentric annular billiard (outer circle, radius R_outer) with a circular
#hole (obstacle) of radius R_inner. Multiply connected; no symmetry reduction
#is applied, so the fundamental domain is the full physical geometry. Both
#curves live in one SimpleDomain (not a CompositeDomain: is_inside for a
#CompositeDomain is a union over subdomains, but a domain-with-a-hole needs
#the intersection "inside outer AND outside inner" that a single SimpleDomain's
#all(...)-over-curves already gives): an outer full CircleSegment
#(orientation=1, physical wall) and an inner full CircleSegment
#(orientation=-1, so is_inside requires being outside the hole). Distinct
#`domain_id`s (1 outer, 2 inner) let `CompositeBIMSolver` group the
#discretized boundary into its two connected components.
struct AnnularBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::SimpleDomain{T}
    symmetries::Vector{AbsSymmetry}
end

function AnnularBilliard(R_outer::T, R_inner::T; center=SVector{2,T}(zero(T),zero(T))) where T<:Real
    R_outer > zero(T) || throw(ArgumentError("R_outer must be positive; received $R_outer"))
    R_inner > zero(T) || throw(ArgumentError("R_inner must be positive; received $R_inner"))
    R_inner < R_outer || throw(ArgumentError("R_inner must be smaller than R_outer; received R_inner=$R_inner, R_outer=$R_outer"))
    c = SVector{2,T}(center)
    bc = SpecularReflection()
    outer = CircleSegment(R_outer, T(2*pi), zero(T), c; bc=bc, orientation=1, domain_id=1, segment_id=1)
    inner = CircleSegment(R_inner, T(2*pi), zero(T), c; bc=bc, orientation=-1, domain_id=2, segment_id=1)
    vertices = SVector{2,T}[c+SVector{2,T}(R_outer,zero(T)), c+SVector{2,T}(R_inner,zero(T))]
    fundamental_domain = SimpleDomain{T}(AbsCurve[outer,inner], vertices, 1)
    symmetries = Vector{AbsSymmetry}(undef,0)
    return AnnularBilliard{T}(fundamental_domain, symmetries)
end
