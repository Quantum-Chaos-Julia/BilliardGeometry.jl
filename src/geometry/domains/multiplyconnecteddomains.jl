################################################################################
########################## MULTIPLY CONNECTED DOMAIN ###########################
################################################################################
# Ω = (⋃ subdomains) \ (⋃ holes). `subdomains` describe only the OUTER region
# in the billiard's fundamental frame (union semantics, as in
# `CompositeDomain`). Every hole is a self-contained `AbsDomainHole` in the
# full physical frame, never reconstructed by the billiard symmetry registry.
################################################################################

"""
    MultiplyConnectedDomain{T<:Real} <: AbsMultiplyConnectedDomain

Planar domain `Ω = (⋃ subdomains) \\ (⋃ holes)`.

## Attributes
* `subdomains::Vector{AbsSimpleDomain}`: Outer region (union), in the fundamental frame of the billiard.
* `holes::Vector{AbsDomainHole}`: Regions removed from the outer region, each in the full physical frame.
* `id::Int64`: Integer id of the domain.

## Description
The curve `orientation` field is never used to identify holes. Every subdomain
and hole boundary is a closed loop traversed with that domain's own region on
the left.
"""
struct MultiplyConnectedDomain{T<:Real} <: AbsMultiplyConnectedDomain
    subdomains::Vector{AbsSimpleDomain}
    holes::Vector{AbsDomainHole}
    id::Int64
end

_domain_eltype(::MultiplyConnectedDomain{T}) where {T} = T

_outer_subdomains(d::AbsSimpleDomain) = AbsSimpleDomain[d]
_outer_subdomains(d::AbsCompositeDomain) = AbsSimpleDomain[d.subdomains...]
_outer_subdomains(::AbsMultiplyConnectedDomain) = throw(ArgumentError("the outer region of a MultiplyConnectedDomain cannot itself be a MultiplyConnectedDomain"))
function _outer_subdomains(domains::AbstractVector)
    subdomains = AbsSimpleDomain[]
    for d in domains
        append!(subdomains, _outer_subdomains(d))
    end
    return subdomains
end

_as_hole(h::AbsDomainHole) = h
_as_hole(d::AbsDomain) = DomainHole(d)

"""
    MultiplyConnectedDomain(outer, holes, id) → domain::MultiplyConnectedDomain

Constructs `Ω = (⋃ subdomains) \\ (⋃ holes)`.

## Arguments
* `outer`: Simple domain, composite domain, or vector of simple domains making up the outer region (flattened into `subdomains`).
* `holes`: Vector of [`AbsDomainHole`](@ref); a bare `AbsDomain` entry is wrapped in a [`DomainHole`](@ref).
* `id`: Integer id of the domain.

## Throws
`ArgumentError` when a hole domain is itself a `MultiplyConnectedDomain`, when
the curve `domain_id`s of different outer/hole objects overlap, or when the
outer exterior boundary or a hole boundary is not a union of closed loops.
"""
function MultiplyConnectedDomain(outer, holes::AbstractVector, id::Integer)
    subdomains = _outer_subdomains(outer)
    isempty(subdomains) && throw(ArgumentError("MultiplyConnectedDomain requires at least one outer subdomain"))
    hs = AbsDomainHole[_as_hole(h) for h in holes]
    _validate_hole_structure(subdomains, hs)
    T = _domain_eltype(first(subdomains))
    return MultiplyConnectedDomain{T}(subdomains, hs, Int64(id))
end

_curve_domain_ids(domains) = Set{Int}(crv.domain_id for d in domains for crv in d.boundary)

function _validate_hole_structure(subdomains::Vector{AbsSimpleDomain}, holes::Vector{AbsDomainHole})
    outer_curves = AbsCurve[crv for sd in subdomains for crv in sd.boundary]
    try
        _split_closed_loops(_connected_geometric_curves(outer_curves))
    catch err
        err isa ArgumentError || rethrow()
        throw(ArgumentError("outer boundary of MultiplyConnectedDomain is not a union of closed loops: " * err.msg))
    end
    seen = _curve_domain_ids(subdomains)
    for (i, h) in enumerate(holes)
        for d in hole_domains(h)
            d isa AbsMultiplyConnectedDomain && throw(ArgumentError("hole $i contains a MultiplyConnectedDomain; holes must be simple or composite domains"))
        end
        try
            boundary_components(h)
        catch err
            err isa ArgumentError || rethrow()
            throw(ArgumentError("boundary of hole $i is not a union of closed loops: " * err.msg))
        end
        ids = _curve_domain_ids(get_all_domains(h))
        isdisjoint(seen, ids) || throw(ArgumentError("curve domain_ids of hole $i overlap those of the outer region or of another hole: $(sort!(collect(intersect(seen, ids))))"))
        union!(seen, ids)
    end
    return nothing
end

"""
    hole_domain_ids(domain::AbsDomain) → ids::Vector{Int}
    hole_domain_ids(billiard::AbsBilliard) → ids::Vector{Int}

Curve `domain_id`s contained in a hole of `domain`/`billiard`; empty when it
has no holes.
"""
hole_domain_ids(::AbsDomain) = Int[]
hole_domain_ids(domain::MultiplyConnectedDomain) = sort!(collect(Set{Int}(crv.domain_id for h in domain.holes for crv in get_all_curves(h))))
hole_domain_ids(billiard::B) where {B<:AbsBilliard} = hole_domain_ids(billiard.fundamental_domain)

"""
    is_hole_curve(domain, crv::AbsCurve) → flag::Bool

`true` when `crv` belongs to a hole of `domain` (a domain or billiard), decided
by the curve's `domain_id`.
"""
is_hole_curve(::AbsDomain, ::AbsCurve) = false
is_hole_curve(domain::MultiplyConnectedDomain, crv::AbsCurve) = crv.domain_id in hole_domain_ids(domain)
is_hole_curve(billiard::B, crv::AbsCurve) where {B<:AbsBilliard} = is_hole_curve(billiard.fundamental_domain, crv)

"""
    component_of_domain_id(domain) → map::Dict{Int,Int}

Maps each curve `domain_id` of `domain` (a domain or billiard) to its solver
group: `0` for the outer region, `i` for the `i`-th hole object.
"""
function component_of_domain_id(domain::MultiplyConnectedDomain)
    map = Dict{Int,Int}()
    for sd in domain.subdomains, crv in sd.boundary
        map[crv.domain_id] = 0
    end
    for (i, h) in enumerate(domain.holes), crv in get_all_curves(h)
        map[crv.domain_id] = i
    end
    return map
end
component_of_domain_id(domain::AbsDomain) = Dict{Int,Int}(crv.domain_id => 0 for crv in get_all_curves(domain))
component_of_domain_id(billiard::B) where {B<:AbsBilliard} = component_of_domain_id(billiard.fundamental_domain)

"""
    genus(domain::AbsDomain) → g::Int64
    genus(billiard::AbsBilliard) → g::Int64

Returns the genus (number of hole loops) of a domain/billiard: `0` for
[`AbsSimpleDomain`](@ref)/[`AbsCompositeDomain`](@ref), and the total loop
count of all holes for [`MultiplyConnectedDomain`](@ref).
"""
genus(domain::AbsDomain) = 0
genus(domain::MultiplyConnectedDomain) = sum(hole_loop_count, domain.holes; init=0)
genus(billiard::B) where {B<:AbsBilliard} = genus(billiard.fundamental_domain)

"""
    validate_holes(domain::AbsMultiplyConnectedDomain, registry::SymmetryRegistry) → nothing

Sample-point validation of the holes of `domain`. For one point on each hole
loop it checks that the point lies inside the outer region reconstructed by
`registry` (an image of the point under some element of the symmetry group,
including the identity, lies in a fundamental outer subdomain) and outside
every other hole.

This is not a proof that holes neither touch nor overlap each other or the
outer boundary.

## Throws
`ArgumentError` naming the offending hole.
"""
function validate_holes(domain::AbsMultiplyConnectedDomain, registry::SymmetryRegistry)
    outer_contains(pt) = any(all(is_inside(sd, pt)) for sd in domain.subdomains)
    for (i, h) in enumerate(domain.holes)
        for loop in boundary_components(h)
            pt = SVector{2}(curve(first(loop), 0.5))
            in_outer = outer_contains(pt) || any(outer_contains(SVector{2}(apply_symmetry(sym, pt))) for sym in registry)
            in_outer || throw(ArgumentError("hole $i is not inside the outer region"))
            for (j, other) in enumerate(domain.holes)
                j == i && continue
                is_inside(other, pt) && throw(ArgumentError("hole $i overlaps hole $j"))
            end
        end
    end
    return nothing
end
