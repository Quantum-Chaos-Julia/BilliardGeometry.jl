
"""
    SimpleDomain{T}<:AbsSimpleDomain

General simple domain represented by an ordered collection of boundary curves
and corner points.

## Fields
* `boundary::Vector{AbsCurve}`: Ordered boundary curves enclosing the domain.
* `corners::Vector{SVector{2,T}}`: Corner or boundary-junction points.
* `id::Int64`: Domain identifier.
"""
struct SimpleDomain{T} <: AbsSimpleDomain where T<:Real
    boundary::Vector{AbsCurve}
    corners::Vector{SVector{2,T}}
    id::Int64
end

"""
    CompositeDomain<:AbsCompositeDomain

Composite domain formed from a collection of simple subdomains.

## Fields
* `subdomains::Vector{AbsSimpleDomain}`: Simple domains comprising the composite domain.
"""
struct CompositeDomain <: AbsCompositeDomain
    subdomains::Vector{AbsSimpleDomain}
end

"""
    reset_ids!(domain::AbsCompositeDomain)

Generate consecutive identifiers for the simple domains contained in a
composite domain.

## Arguments
* `domain::AbsCompositeDomain`: Composite domain whose subdomain identifiers are reset.

## Returns
No explicit return value is specified by this function.
"""
function reset_ids!(domain::AbsCompositeDomain)
    for (i, sd)  in enumerate(domain.subdomains)
        @set sd.id = i
    end
end