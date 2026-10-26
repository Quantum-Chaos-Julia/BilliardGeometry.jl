"""
    SimpleDomain{T} <: AbsSimpleDomain

Region bounded by a list of curves.

## Attributes
* `boundary::Vector{AbsCurve}`: Boundary curves, including symmetry walls.
* `corners::Vector{SVector{2,T}}`: Corner points of the boundary.
* `id::Int64`: Id of the domain.
"""
struct SimpleDomain{T} <: AbsSimpleDomain where T<:Real
    boundary::Vector{AbsCurve}
    corners::Vector{SVector{2,T}}
    id::Int64
end

"""
    CompositeDomain <: AbsCompositeDomain

Union of several simple domains glued along [`Transparent`](@ref) curves, for
example the stadium of [`StadiumBilliard`](@ref).

## Attributes
* `subdomains::Vector{AbsSimpleDomain}`: The non-overlapping pieces.
"""
struct CompositeDomain <: AbsCompositeDomain
    subdomains::Vector{AbsSimpleDomain}
end

_domain_eltype(::SimpleDomain{T}) where T = T
_domain_eltype(::Polygon{T}) where T = T
_domain_eltype(::CircleWedge{T}) where T = T
_domain_eltype(domain::CompositeDomain) = _domain_eltype(first(domain.subdomains))

"""
    reset_ids!(domain::AbsCompositeDomain) → domain

Renumber the `id` field of the subdomains in place as `1, 2, ...` in list order.
The `domain_id` of the curves is not changed.
"""
function reset_ids!(domain::AbsCompositeDomain)
    for (i, sd)  in enumerate(domain.subdomains)
        domain.subdomains[i] = @set sd.id = i
    end
    return domain
end