struct SpecularReflection <:AbsBoundaryCondition end

struct Transparent <:AbsBoundaryCondition
    next_id::Int64
end

struct PeriodicX <:AbsBoundaryCondition
    next_id::Int64
end

"""
SymmetryWall <: AbsBoundaryCondition

Boundary condition marking a curve as a symmetry wall: crossing it
corresponds to applying the billiard's registered symmetry generator
identified by `sym_id` (see [`SymmetryRegistry`](@ref)/
[`register_symmetries`](@ref)), landing in the neighboring fundamental-domain
copy `sector_id`. Covers both reflection walls and rotation wedge cuts
(`AbsRotation` generators), which previously had no dedicated boundary
condition at all.

## Attributes
* `sym_id::Int64`: `sym_id` of the billiard's registered symmetry generator this wall corresponds to.
* `sector_id::Int64`: Index of the neighboring fundamental-domain copy this wall borders.
"""
struct SymmetryWall <:AbsBoundaryCondition
    sym_id::Int64
    sector_id::Int64
end

struct QuantumSolverIgnore <:AbsBoundaryCondition end

# Shared body of get_boundary_curves(::AbsSimpleDomain): filters the flat
# domain.boundary to physical (SpecularReflection) curves and connects them
# into contiguous runs.
function _connected_physical_curves(boundary::Vector{AbsCurve})
    is_outer(crv) = typeof(crv.bc) <: SpecularReflection
    physical = filter(is_outer, boundary)
    return connect_curves(physical)
end

function get_boundary_curves(domain::D) where D<:AbsSimpleDomain
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

"""
    get_outer_boundary_curves(domain::AbsDomain) → curves::Vector{AbsCurve}
    get_outer_boundary_curves(billiard::AbsBilliard) → curves::Vector{AbsCurve}

Physical (specular) curves of the outer boundary in the fundamental frame of
the billiard. For domains without holes this equals [`get_boundary_curves`](@ref).
"""
get_outer_boundary_curves(domain::AbsDomain) = get_boundary_curves(domain)
function get_outer_boundary_curves(domain::AbsMultiplyConnectedDomain)
    boundary = Vector{AbsCurve}()
    for sd in domain.subdomains
        append!(boundary, get_boundary_curves(sd))
    end
    return connect_curves(boundary)
end
get_outer_boundary_curves(billiard::B) where B<:AbsBilliard = get_outer_boundary_curves(billiard.fundamental_domain)

"""
    get_hole_boundary_curves(domain::AbsDomain) → curves::Vector{AbsCurve}
    get_hole_boundary_curves(billiard::AbsBilliard) → curves::Vector{AbsCurve}

Physical curves of every hole loop in the full physical frame, one loop after
the other, each traversed with the hole region on the left. Empty when there
are no holes.
"""
get_hole_boundary_curves(::AbsDomain) = AbsCurve[]
function get_hole_boundary_curves(domain::AbsMultiplyConnectedDomain)
    curves = AbsCurve[]
    for h in domain.holes, loop in boundary_components(h)
        append!(curves, loop)
    end
    return curves
end
get_hole_boundary_curves(billiard::B) where B<:AbsBilliard = get_hole_boundary_curves(billiard.fundamental_domain)

"""
    get_boundary_curves(domain::AbsMultiplyConnectedDomain) → curves::Vector{AbsCurve}

Returns the outer physical curves only, see [`get_outer_boundary_curves`](@ref).
Hole curves live in a different frame, see [`get_hole_boundary_curves`](@ref).
"""
get_boundary_curves(domain::D) where D<:AbsMultiplyConnectedDomain = get_outer_boundary_curves(domain)

"""
    get_boundary_curves(billiard::AbsBilliard) → curves::Vector{AbsCurve}

Physical curves of the fundamental domain. For a billiard with a
[`MultiplyConnectedDomain`](@ref) and holes: an empty symmetry registry returns the
outer curves followed by the hole curves (the frames coincide); a nonempty
registry throws `ArgumentError`, use [`full_boundary_components`](@ref).
"""
function get_boundary_curves(billiard::B) where B<:AbsBilliard
    return _billiard_boundary_curves(billiard.fundamental_domain, billiard)
end
_billiard_boundary_curves(domain::AbsDomain, ::AbsBilliard) = get_boundary_curves(domain)
function _billiard_boundary_curves(domain::AbsMultiplyConnectedDomain, billiard::AbsBilliard)
    outer = get_outer_boundary_curves(domain)
    holes = get_hole_boundary_curves(domain)
    isempty(holes) && return outer
    isempty(billiard.symmetries) || throw(ArgumentError("get_boundary_curves is ambiguous for a billiard with holes and a nonempty symmetry registry (outer curves are in the fundamental frame, hole curves in the full frame); use full_boundary_components"))
    return AbsCurve[outer; holes]
end

"""
    get_all_domains(domain::AbsDomain) → domains::Vector{AbsSimpleDomain}
    get_all_domains(billiard::AbsBilliard) → domains::Vector{AbsSimpleDomain}

Simple domains making up `domain`/`billiard`. For a [`MultiplyConnectedDomain`](@ref)
the outer subdomains are followed by the simple domains of every hole, so
holes appear here but are not filled regions.
"""
get_all_domains(domain::AbsSimpleDomain) = AbsSimpleDomain[domain]
get_all_domains(domain::AbsCompositeDomain) = domain.subdomains
function get_all_domains(domain::AbsMultiplyConnectedDomain)
    domains = AbsSimpleDomain[domain.subdomains...]
    for h in domain.holes
        append!(domains, get_all_domains(h))
    end
    return domains
end
get_all_domains(billiard::B) where B<:AbsBilliard = get_all_domains(billiard.fundamental_domain)

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
    curves = Vector{AbsCurve}()
    for dom in get_all_domains(domain)
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


# Returns `domain` with the boundary condition of the curve identified by
# (domain_id, segment_id) replaced, or `nothing` when `domain` has no such curve.
function _set_curve_bc(domain::AbsSimpleDomain, domain_id, segment_id, bc::AbsBoundaryCondition)
    idx = findfirst(crv -> crv.domain_id == domain_id && crv.segment_id == segment_id, domain.boundary)
    isnothing(idx) && return nothing
    return @set domain.boundary[idx].bc = bc
end

function _set_curve_bc(domain::AbsCompositeDomain, domain_id, segment_id, bc::AbsBoundaryCondition)
    for (i, sd) in enumerate(domain.subdomains)
        new_sd = _set_curve_bc(sd, domain_id, segment_id, bc)
        isnothing(new_sd) && continue
        return @set domain.subdomains[i] = new_sd
    end
    return nothing
end

function _set_curve_bc(domain::AbsMultiplyConnectedDomain, domain_id, segment_id, bc::AbsBoundaryCondition)
    for (i, sd) in enumerate(domain.subdomains)
        new_sd = _set_curve_bc(sd, domain_id, segment_id, bc)
        isnothing(new_sd) && continue
        return @set domain.subdomains[i] = new_sd
    end
    for (i, h) in enumerate(domain.holes)
        new_h = _set_curve_bc(h, domain_id, segment_id, bc)
        isnothing(new_h) && continue
        return @set domain.holes[i] = new_h
    end
    return nothing
end

function update_boundary_condition(billiard::B, domain_id, segment_id, bc::BC) where {B<:AbsBilliard, BC<:AbsBoundaryCondition}
    fundamental_domain = billiard.fundamental_domain

    if typeof(fundamental_domain) <: AbsMultiplyConnectedDomain
        new_domain = _set_curve_bc(fundamental_domain, domain_id, segment_id, bc)
        isnothing(new_domain) && throw(ArgumentError("No curve with domain_id=$domain_id and segment_id=$segment_id"))
        return @set billiard.fundamental_domain = new_domain
    elseif typeof(fundamental_domain) <: AbsCompositeDomain
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
For [`AbsMultiplyConnectedDomain`](@ref), the first component is the outer
boundary (all subdomain physical curves, connected, fundamental frame) and is
followed by one component per hole loop (full frame, hole region on the left),
see [`boundary_component_is_hole`](@ref).
"""
boundary_components(domain::AbsDomain) = [get_boundary_curves(domain)]
function boundary_components(domain::AbsMultiplyConnectedDomain)
    components = Vector{Vector{AbsCurve}}()
    push!(components, get_outer_boundary_curves(domain))
    for h in domain.holes
        append!(components, boundary_components(h))
    end
    return components
end
boundary_components(billiard::B) where B<:AbsBilliard = boundary_components(billiard.fundamental_domain)

"""
    boundary_component_is_hole(domain) → flags::Vector{Bool}

`true` for every component of [`boundary_components`](@ref) that is a hole loop,
for a domain or billiard. All `false` when there are no holes.
"""
boundary_component_is_hole(domain::AbsDomain) = Bool[false]
boundary_component_is_hole(domain::AbsMultiplyConnectedDomain) = Bool[false; fill(true, genus(domain))]
boundary_component_is_hole(billiard::B) where B<:AbsBilliard = boundary_component_is_hole(billiard.fundamental_domain)

# Connects the complete exterior boundary of a simple/composite domain --
# every curve regardless of boundary-condition type EXCEPT Transparent
# (which marks an internal join stitching two AbsCompositeDomain subdomains
# together, not a real exterior edge). Unlike get_boundary_curves/
# boundary_components, which retain only SpecularReflection curves for
# solver-matrix construction, this keeps SymmetryWall/QuantumSolverIgnore/
# PeriodicX curves too, since they are real exterior edges of the
# fundamental domain.
_is_exterior_curve(crv) = !(typeof(crv.bc) <: Transparent)
_connected_geometric_curves(boundary::Vector{AbsCurve}) = connect_curves(filter(_is_exterior_curve, boundary))

_geometric_boundary(domain::AbsDomain) = _connected_geometric_curves(domain.boundary)
function _geometric_boundary(composite_domain::AbsCompositeDomain)
    boundary = Vector{AbsCurve}()
    for domain in composite_domain.subdomains
        append!(boundary, domain.boundary)
    end
    return _connected_geometric_curves(boundary)
end

"""
    geometric_boundary_components(domain::AbsDomain) → components::Vector{Vector{AbsCurve}}
    geometric_boundary_components(billiard::AbsBilliard) → components::Vector{Vector{AbsCurve}}

Returns the complete exterior boundary curves of `domain`/`billiard`,
grouped one vector per connected component (outer boundary first, then each
hole for [`AbsMultiplyConnectedDomain`](@ref)).

## Description
Unlike [`boundary_components`](@ref)/[`get_boundary_curves`](@ref), which
retain only `SpecularReflection` curves (the ones a solver actually places
quadrature points on), this keeps every real exterior curve regardless of
boundary condition -- including `SymmetryWall` cuts closing off a
symmetry-reduced fundamental domain -- while still excluding `Transparent`
curves (internal joins stitching together an [`AbsCompositeDomain`](@ref)'s
subdomains). This is needed by [`corner_angles`](@ref) to correctly locate
every true corner of the fundamental domain, including corners where a
physical curve meets a symmetry-wall cut (e.g. [`LimaconBilliard`](@ref)'s
reflection axis).
"""
geometric_boundary_components(domain::AbsDomain) = [_geometric_boundary(domain)]
function geometric_boundary_components(domain::AbsMultiplyConnectedDomain)
    outer = AbsCurve[crv for sd in domain.subdomains for crv in sd.boundary]
    components = Vector{Vector{AbsCurve}}()
    push!(components, _connected_geometric_curves(outer))
    for h in domain.holes
        hole_curves = AbsCurve[crv for dom in get_all_domains(h) for crv in dom.boundary]
        append!(components, _split_closed_loops(_connected_geometric_curves(hole_curves)))
    end
    return components
end
geometric_boundary_components(billiard::B) where B<:AbsBilliard = geometric_boundary_components(billiard.fundamental_domain)