################################################################################
################################## DOMAIN HOLES ################################
################################################################################
# A hole is a region subtracted from the outer region of a
# `MultiplyConnectedDomain`. Holes live in the full physical frame and are
# self-contained: they are never replicated or reduced by the billiard's
# symmetry registry.
#
# Traversal convention: every closed hole loop is traversed with the hole
# region on its left (counterclockwise for a convex hole). The curve
# `orientation` field only selects which side of one curve is "inside" and is
# never read to decide holes, components or areas.
################################################################################

"""
    DomainHole{D<:AbsDomain} <: AbsDomainHole

`DomainHole` is the concrete [`AbsDomainHole`](@ref) wrapping one simple or composite domain that is subtracted from the outer region of a [`MultiplyConnectedDomain`](@ref).

## Description
The hole lives in the full physical frame and has no symmetry: it is never
replicated or reduced by the billiard symmetry registry, so it is
[`symmetry_compatible`](@ref) only with an empty registry. Every closed loop of
its boundary is traversed with the hole region on the left (counterclockwise for
a convex hole); the curve `orientation` field is never read to decide hole
status.

## Attributes
* `domain::D`: Full hole region, a simple or composite domain whose specular boundary consists of closed loops in the full physical frame, each traversed with the hole region on the left.

## API
The following functions can be evaluated for this type:
- [`hole_domains`](@ref), [`hole_loop_count`](@ref), [`hole_leaves`](@ref)
- [`is_inside`](@ref), [`area`](@ref), [`symmetry_compatible`](@ref)
- [`get_boundary_curves`](@ref), [`get_all_domains`](@ref), [`get_all_curves`](@ref), [`boundary_components`](@ref)
"""
struct DomainHole{D<:AbsDomain} <: AbsDomainHole
    domain::D
end

_domain_eltype(h::DomainHole) = _domain_eltype(h.domain)

"""
    hole_domains(hole::AbsDomainHole) → domains::Vector{AbsDomain}

Domains making up the region of `hole`; for a [`DomainHole`](@ref) the single wrapped domain.
"""
hole_domains(h::DomainHole) = AbsDomain[h.domain]

"""
    hole_leaves(hole::AbsDomainHole) → leaves::Vector{AbsDomainHole}

Elementary holes `hole` is composed of; a [`DomainHole`](@ref) is its own single leaf.
"""
hole_leaves(h::DomainHole) = AbsDomainHole[h]
get_all_domains(h::DomainHole) = get_all_domains(h.domain)
get_all_curves(h::DomainHole) = get_all_curves(h.domain)
get_boundary_curves(h::DomainHole) = get_boundary_curves(h.domain)
boundary_components(h::DomainHole) = _split_closed_loops(get_boundary_curves(h.domain))

"""
    hole_loop_count(hole::AbsDomainHole) → n::Int

Number of closed boundary loops of `hole`, each traversed with the hole region
on its left. Their sum over all holes is the [`genus`](@ref) of the domain.
"""
hole_loop_count(h::DomainHole) = length(boundary_components(h))

function _set_curve_bc(h::DomainHole, domain_id, segment_id, bc::AbsBoundaryCondition)
    new_domain = _set_curve_bc(h.domain, domain_id, segment_id, bc)
    isnothing(new_domain) && return nothing
    return DomainHole(new_domain)
end

# Region membership of a point for a simple or composite domain.
_region_contains(domain::AbsSimpleDomain, pt) = all(is_inside(domain, pt))
_region_contains(domain::AbsCompositeDomain, pt) = any(_region_contains(sd, pt) for sd in domain.subdomains)

is_inside(h::DomainHole, pt::SVector{2,T}) where {T<:Real} = _region_contains(h.domain, pt)
is_inside(h::DomainHole, pts::AbstractArray) = [is_inside(h, pt) for pt in pts]

"""
    symmetry_compatible(hole::AbsDomainHole, registry::SymmetryRegistry) → flag::Bool

`true` when the billiard symmetry group `registry` acts on `hole` through exact
group algebra, so that a symmetry-reduced solver may treat the hole through
its orbit structure. A [`DomainHole`](@ref) has no symmetry and is compatible
only with an empty registry.
"""
symmetry_compatible(::DomainHole, registry::SymmetryRegistry) = isempty(registry)
