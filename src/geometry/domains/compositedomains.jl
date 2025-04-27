
struct CompositeDomain <: AbsCompositeDomain
    subdomains::Vector{AbsSimpleDomain}
end

#generate new ids for the simple domains that make up the complex domain
function reset_ids!(domain::AbsCompositeDomain)
    for (i, sd)  in enumerate(domain.subdomains)
        @set sd.id = i
    end
end