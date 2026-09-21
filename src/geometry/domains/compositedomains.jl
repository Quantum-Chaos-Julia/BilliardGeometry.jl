
struct SimpleDomain{T} <: AbsSimpleDomain where T<:Real
    boundary::Vector{AbsCurve}
    corners::Vector{SVector{2,T}}
    id::Int64
end

struct CompositeDomain <: AbsCompositeDomain
    subdomains::Vector{AbsSimpleDomain}
end

function reset_ids!(domain::AbsCompositeDomain)
    for (i, sd)  in enumerate(domain.subdomains)
        domain.subdomains[i] = @set sd.id = i
    end
    return domain
end