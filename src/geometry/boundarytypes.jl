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

function get_boundary_curves(domain::D) where D<:AbsSimpleDomain
    is_outer(crv) = typeof(crv.bc) <: SpecularReflection
    boundary = filter(is_outer, domain.boundary)
    return connect_curves(boundary)
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