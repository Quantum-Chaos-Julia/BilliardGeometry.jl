"""
    SpecularReflection<:AbsBoundaryCondition

Specular reflecting boundary condition.

A trajectory incident on a boundary carrying this condition is reflected
specularly.
"""
struct SpecularReflection<:AbsBoundaryCondition end

"""
    Transparent<:AbsBoundaryCondition

Transparent boundary condition connecting the current domain to another domain.

## Fields
* `next_id::Int64`: Identifier of the domain reached across the boundary.
"""
struct Transparent<:AbsBoundaryCondition
    next_id::Int64
end

"""
    PeriodicX<:AbsBoundaryCondition

Periodic boundary condition in the x-direction.

## Fields
* `next_id::Int64`: Identifier of the domain connected through the periodic boundary.
"""
struct PeriodicX<:AbsBoundaryCondition
    next_id::Int64
end

"""
    ReflectionSymmetry{S}<:AbsBoundaryCondition where S<:AbsReflection

Boundary condition representing a reflection-symmetry boundary.

## Fields
* `symmetry::S`: Reflection symmetry associated with the boundary.
* `N_sectors::Int64`: Number of symmetry sectors associated with the construction.
"""
struct ReflectionSymmetry{S}<:AbsBoundaryCondition where S<:AbsReflection
    symmetry::S
    N_sectors::Int64
end

"""
    QuantumSolverIgnore<:AbsBoundaryCondition

Boundary condition marking a boundary segment that should be ignored by the
quantum solver.
"""
struct QuantumSolverIgnore<:AbsBoundaryCondition end

"""
    get_boundary_curves(domain::D) where D<:AbsSimpleDomain

Return the connected boundary formed from the specularly reflecting curves of a
simple domain.

Only curves whose boundary condition is `SpecularReflection` are retained.

## Arguments
* `domain::D`: Simple domain.

## Returns
* Connected boundary formed from the selected curves.
"""
function get_boundary_curves(domain::D) where D<:AbsSimpleDomain
    is_outer(crv)=typeof(crv.bc)<:SpecularReflection
    boundary=filter(is_outer,domain.boundary)
    return connect_curves(boundary)
end

"""
    get_boundary_curves(composite_domain::D) where D<:AbsCompositeDomain

Return the connected specularly reflecting boundary of a composite domain.

The boundary curves are collected from all simple subdomains and then passed to
`connect_curves`.

## Arguments
* `composite_domain::D`: Composite domain.

## Returns
* Connected boundary formed from the selected curves of all subdomains.
"""
function get_boundary_curves(composite_domain::D) where D<:AbsCompositeDomain
    boundary=Vector{AbsCurve}()
    for domain in composite_domain.subdomains
        subboundary=get_boundary_curves(domain)
        append!(boundary,subboundary)
    end
    return connect_curves(boundary)
end

"""
    get_boundary_curves(billiard::B) where B<:AbsBilliard

Return the boundary curves of a billiard's fundamental domain.

## Arguments
* `billiard::B`: Billiard.

## Returns
* Result of `get_boundary_curves` applied to `billiard.fundamental_domain`.
"""
function get_boundary_curves(billiard::B) where B<:AbsBilliard
    return get_boundary_curves(billiard.fundamental_domain)
end

"""
    get_all_domains(billiard::B) where B<:AbsBilliard

Return all simple domains contained in a billiard's fundamental domain.

For a composite fundamental domain, its `subdomains` are returned. Otherwise
the fundamental domain is returned as a one-element collection.

## Arguments
* `billiard::B`: Billiard.

## Returns
* Collection of simple domains.
"""
function get_all_domains(billiard::B) where B<:AbsBilliard
    domain=billiard.fundamental_domain
    if typeof(domain)<:AbsCompositeDomain
        subdomains=domain.subdomains
    else
        subdomains=[domain]
    end
    return subdomains
end

"""
    get_domain(billiard::B,id) where B<:AbsBilliard

Return the domain of a billiard whose identifier equals `id`.

## Arguments
* `billiard::B`: Billiard.
* `id`: Domain identifier.

## Returns
* Matching domain, if one is found.
"""
function get_domain(billiard::B,id) where B<:AbsBilliard
    domains=get_all_domains(billiard)
    for dom in domains
        if dom.id==id
            return dom
        end
    end
end

"""
    get_all_curves(billiard::B) where B<:AbsBilliard

Return all boundary curves belonging to the domains of a billiard's fundamental
domain.

## Arguments
* `billiard::B`: Billiard.

## Returns
* `Vector{AbsCurve}` containing all domain boundary curves.
"""
function get_all_curves(billiard::B) where B<:AbsBilliard
    subdomains=get_all_domains(billiard)
    curves=Vector{AbsCurve}()
    for dom in subdomains
        append!(curves,dom.boundary)
    end
    return curves
end

"""
    get_all_curves(domain::D) where D<:AbsDomain

Return all boundary curves belonging to a domain.

For a composite domain, curves are collected from all subdomains. Otherwise the
function uses the supplied domain directly.

## Arguments
* `domain::D`: Domain.

## Returns
* `Vector{AbsCurve}` containing the collected boundary curves.
"""
function get_all_curves(domain::D) where D<:AbsDomain
    if typeof(domain)<:AbsCompositeDomain
        subdomains=domain.subdomains
    else
        subdomains=domain
    end
    curves=Vector{AbsCurve}()
    for dom in subdomains
        append!(curves,dom.boundary)
    end
    return curves
end

"""
    get_curve(billiard::B,domain_id,segment_id) where B<:AbsBilliard

Return the boundary curve of a billiard matching `domain_id` and `segment_id`.

## Arguments
* `billiard::B`: Billiard.
* `domain_id`: Domain identifier.
* `segment_id`: Segment identifier.

## Returns
* Matching curve, if one is found.
"""
function get_curve(billiard::B,domain_id,segment_id) where B<:AbsBilliard
    curves=get_all_curves(billiard)
    for crv in curves
        if domain_id==crv.domain_id&&segment_id==crv.segment_id
            return crv
        end
    end
end

"""
    get_curve(composite_curve::C,domain_id,segment_id) where C<:AbsCompositeCurve

Return the subcurve of a composite curve matching `domain_id` and `segment_id`.

## Arguments
* `composite_curve::C`: Composite curve.
* `domain_id`: Domain identifier.
* `segment_id`: Segment identifier.

## Returns
* Matching subcurve, if one is found.
"""
function get_curve(composite_curve::C,domain_id,segment_id) where C<:AbsCompositeCurve
    curves=composite_curve.subcurves
    for crv in curves
        if domain_id==crv.domain_id&&segment_id==crv.segment_id
            return crv
        end
    end
end

"""
    update_boundary_condition(billiard::B,domain_id,segment_id,bc::BC) where {B<:AbsBilliard,BC<:AbsBoundaryCondition}

Return a billiard with the boundary condition of the selected fundamental-domain
segment replaced by `bc`.

For composite fundamental domains, `domain_id` selects the subdomain and
`segment_id` selects the curve within that subdomain. For simple fundamental
domains, `segment_id` selects the curve directly.

## Arguments
* `billiard::B`: Billiard.
* `domain_id`: Domain identifier.
* `segment_id`: Segment identifier.
* `bc::BC`: New boundary condition.

## Returns
* Updated billiard produced by `@set`.
"""
function update_boundary_condition(billiard::B,domain_id,segment_id,bc::BC) where {B<:AbsBilliard,BC<:AbsBoundaryCondition}
    fundamental_domain=billiard.fundamental_domain
    if typeof(fundamental_domain)<:AbsCompositeDomain
        subdomains=fundamental_domain.subdomains
        dom_idx=findfirst(d->d.id==domain_id,subdomains)
        curves=subdomains[dom_idx].boundary
        seg_idx=findfirst(crv->crv.segment_id==segment_id,curves)
        return @set billiard.fundamental_domain.subdomains[dom_idx].boundary[seg_idx].bc=bc
    else
        curves=fundamental_domain.boundary
        seg_idx=findfirst(crv->crv.segment_id==segment_id,curves)
        return @set billiard.fundamental_domain.boundary[seg_idx].bc=bc
    end
end