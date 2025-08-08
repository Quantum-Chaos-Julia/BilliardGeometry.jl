struct CompositeCurve{T} <: AbsCompositeCurve where {T<:Real}
    subcurves::Vector{AbsCurve}
    end_lengths::Vector{T}
    corners::Vector{SVector{2,T}}
    length::T
end

function CompositeCurve(subcurves) 
    type = typeof(subcurves[1]).parameters[1] 
    end_lengths = Vector{type}(undef,0)
    corners = Vector{SVector{2,type}}(undef,0)

    crvs = connect_curves(subcurves)
    push!(corners, curve(crvs[1],0.0))

    L = zero(type)
    for crv in crvs
        L += crv.length
        push!(end_lengths, L)
        push!(corners, curve(crv,1.0)) 
    end
    return CompositeCurve{type}(crvs,end_lengths,corners,L)
end
