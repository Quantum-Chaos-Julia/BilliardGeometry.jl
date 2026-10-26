#annular billiard (outer circle, radius R_outer) with a circular hole
#(obstacle) of radius R_inner whose center is shifted along the x axis by
#`shift`. Multiply connected; no symmetry reduction is applied, so the
#fundamental domain is the full physical geometry. The outer region is one
#SimpleDomain holding a full CircleSegment (domain_id=1); the hole is a
#DomainHole around a SimpleDomain holding a full counterclockwise
#CircleSegment (orientation=1, domain_id=2), i.e. traversed with the hole
#region on its left.
struct AnnularBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::MultiplyConnectedDomain{T}
    symmetries::SymmetryRegistry
end

function AnnularBilliard(R_inner::T, shift::T; R_outer = one(T), center=SVector{2,T}(zero(T),zero(T))) where T<:Real
    R_outer > zero(T) || throw(ArgumentError("R_outer must be positive; received $R_outer"))
    R_inner > zero(T) || throw(ArgumentError("R_inner must be positive; received $R_inner"))
    R_inner < R_outer || throw(ArgumentError("R_inner must be smaller than R_outer; received R_inner=$R_inner, R_outer=$R_outer"))
    c = SVector{2,T}(center)
    bc = SpecularReflection()
    outer = CircleSegment(R_outer, T(2*pi), zero(T), c; bc=bc, orientation=1, domain_id=1, segment_id=1)
    hole_center = c+SVector{2,T}(shift,zero(T))
    inner = CircleSegment(R_inner, T(2*pi), zero(T), hole_center; bc=bc, orientation=1, domain_id=2, segment_id=1)
    outer_domain = SimpleDomain{T}(AbsCurve[outer], SVector{2,T}[c+SVector{2,T}(R_outer,zero(T))], 1)
    hole_domain = SimpleDomain{T}(AbsCurve[inner], SVector{2,T}[hole_center+SVector{2,T}(R_inner,zero(T))], 2)
    symmetries = register_symmetries()
    fundamental_domain = MultiplyConnectedDomain(outer_domain, [DomainHole(hole_domain)], 1)
    validate_holes(fundamental_domain, symmetries)
    return AnnularBilliard{T}(fundamental_domain, symmetries)
end

