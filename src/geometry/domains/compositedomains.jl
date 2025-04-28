
struct SimpleDomain{T} <: AbsSimpleDomain where T<:Real
    boundary::Vector{AbsCurve}
    corners::Vector{SVector{2,T}}
    id::Int64
end

struct CompositeDomain <: AbsCompositeDomain
    subdomains::Vector{AbsSimpleDomain}
end

#generate new ids for the simple domains that make up the complex domain
function reset_ids!(domain::AbsCompositeDomain)
    for (i, sd)  in enumerate(domain.subdomains)
        @set sd.id = i
    end
end