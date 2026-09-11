struct SpecularReflection <:AbsBoundaryCondition end

struct Transparent <:AbsBoundaryCondition
    next_id::Int64
end

struct PeriodicX <:AbsBoundaryCondition
    next_id::Int64
end

struct ReflectionSymmetry{S} <:AbsBoundaryCondition where S<:AbsReflection
    symmetry::S
    N_sectors::Int64
end

struct QuantumSolverIgnore <:AbsBoundaryCondition end

# Shared body of get_boundary_curves(::AbsSimpleDomain)/(::AbsMultiplyConnectedDomain):
# both store their curves flat in domain.boundary, filter to physical
# (SpecularReflection) curves, and connect them into contiguous runs.
function _connected_physical_curves(boundary::Vector{AbsCurve})
    is_outer(crv) = typeof(crv.bc) <: SpecularReflection
    physical = filter(is_outer, boundary)
    return connect_curves(physical)
end

function get_boundary_curves(domain::D) where D<:AbsSimpleDomain
    return _connected_physical_curves(domain.boundary)
end

function get_boundary_curves(domain::D) where D<:AbsMultiplyConnectedDomain
    return _connected_physical_curves(domain.boundary)
end

function get_boundary_curves(composite_domain::D) where D<:AbsCompositeDomain
    boundary = Vector{AbsCurve}()
    for domain in composite_domain.subdomains
        subboundary = get_boundary_curves(domain)
        append!(boundary,subboundary)
    end
    return connect_curves(boundary)
end

function get_boundary_curves(billiard::B) where B<:AbsBilliard
    return get_boundary_curves(billiard.fundamental_domain)
end


function get_all_domains(billiard::B) where B<:AbsBilliard
    domain = billiard.fundamental_domain
    if typeof(domain) <: AbsCompositeDomain
        subdomains = domain.subdomains
    else
        subdomains = [domain]
    end
    return subdomains
end

function get_domain(billiard::B, id) where B<:AbsBilliard
    domains = get_all_domains(billiard)
    for dom in domains
        if dom.id == id
            return dom
        end
    end
end

function get_all_curves(billiard::B) where B<:AbsBilliard
    subdomains = get_all_domains(billiard)
    curves = Vector{AbsCurve}()
    for dom in subdomains
        append!(curves, dom.boundary)
    end
    return curves
end

function get_all_curves(domain::D) where D<:AbsDomain
    if typeof(domain) <: AbsCompositeDomain
        subdomains = domain.subdomains
    else
        subdomains = domain
    end
    curves = Vector{AbsCurve}()
    for dom in subdomains
        append!(curves, dom.boundary)
    end
    return curves
end

function get_curve(billiard::B, domain_id, segment_id) where B<:AbsBilliard
    curves = get_all_curves(billiard)
    for crv in curves
        if (domain_id == crv.domain_id && segment_id == crv.segment_id)
            return crv
        end
    end
end

function get_curve(composite_curve::C, domain_id, segment_id) where C<:AbsCompositeCurve
    curves = composite_curve.subcurves
    for crv in curves
        if (domain_id == crv.domain_id && segment_id == crv.segment_id)
            return crv
        end
    end
end


function update_boundary_condition(billiard::B, domain_id, segment_id, bc::BC) where {B<:AbsBilliard, BC<:AbsBoundaryCondition}
    fundamental_domain = billiard.fundamental_domain
    
    if typeof(fundamental_domain) <: AbsCompositeDomain
        subdomains = fundamental_domain.subdomains
        dom_idx = findfirst(d -> d.id == domain_id, subdomains)
        curves = subdomains[dom_idx].boundary
        seg_idx = findfirst(crv -> crv.segment_id == segment_id, curves)
        return @set billiard.fundamental_domain.subdomains[dom_idx].boundary[seg_idx].bc = bc
    else
        curves =fundamental_domain.boundary
        seg_idx = findfirst(crv -> crv.segment_id == segment_id, curves) 
        return @set billiard.fundamental_domain.boundary[seg_idx].bc = bc
    end
end

"""
    genus(domain::AbsDomain) → g::Int64
    genus(billiard::AbsBilliard) → g::Int64

Returns the genus (number of holes) of a domain/billiard: `0` for
[`AbsSimpleDomain`](@ref)/[`AbsCompositeDomain`](@ref), and the stored
`genus` field for [`AbsMultiplyConnectedDomain`](@ref) (e.g. [`AnnularBilliard`](@ref)).
"""
genus(domain::AbsDomain) = 0
genus(domain::AbsMultiplyConnectedDomain) = domain.genus
genus(billiard::B) where B<:AbsBilliard = genus(billiard.fundamental_domain)

# Groups a flat physical-boundary curve list into connected components by
# curve domain_id, preserving first-seen domain_id order and within-group
# curve order.
function _group_curves_by_domain_id(boundary::Vector{AbsCurve})
    ids = Int[]
    groups = Vector{Vector{AbsCurve}}()
    @inbounds for c in boundary
        idx = findfirst(==(c.domain_id), ids)
        if idx === nothing
            push!(ids, c.domain_id)
            push!(groups, AbsCurve[c])
        else
            push!(groups[idx], c)
        end
    end
    return groups
end

"""
    boundary_components(domain::AbsDomain) → components::Vector{Vector{AbsCurve}}
    boundary_components(billiard::AbsBilliard) → components::Vector{Vector{AbsCurve}}

Returns the connected physical boundary curves of `domain`/`billiard`, grouped
one vector per connected component (outer boundary first, then each hole for
[`AbsMultiplyConnectedDomain`](@ref)).

## Description
For [`AbsSimpleDomain`](@ref)/[`AbsCompositeDomain`](@ref) (genus `0`), this is
a single-element vector containing [`get_boundary_curves`](@ref)'s result. For
[`AbsMultiplyConnectedDomain`](@ref), the connected physical curves are
grouped by each curve's `domain_id`, preserving first-seen order (outer
boundary's `domain_id` listed first in `domain.boundary` by construction, see
[`MultiplyConnectedDomain`](@ref)).
"""
boundary_components(domain::AbsDomain) = [get_boundary_curves(domain)]
boundary_components(domain::AbsMultiplyConnectedDomain) = _group_curves_by_domain_id(get_boundary_curves(domain))
boundary_components(billiard::B) where B<:AbsBilliard = boundary_components(billiard.fundamental_domain)