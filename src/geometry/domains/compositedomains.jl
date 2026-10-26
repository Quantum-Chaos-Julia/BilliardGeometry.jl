
struct SimpleDomain{T} <: AbsSimpleDomain where T<:Real
    boundary::Vector{AbsCurve}
    corners::Vector{SVector{2,T}}
    id::Int64
end

struct CompositeDomain <: AbsCompositeDomain
    subdomains::Vector{AbsSimpleDomain}
end

_domain_eltype(::SimpleDomain{T}) where T = T
_domain_eltype(::Polygon{T}) where T = T
_domain_eltype(::CircleWedge{T}) where T = T
_domain_eltype(domain::CompositeDomain) = _domain_eltype(first(domain.subdomains))

function reset_ids!(domain::AbsCompositeDomain)
    for (i, sd)  in enumerate(domain.subdomains)
        domain.subdomains[i] = @set sd.id = i
    end
    return domain
end