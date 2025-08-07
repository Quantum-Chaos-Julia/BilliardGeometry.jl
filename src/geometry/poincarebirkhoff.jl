struct PoincareBirkhoff{T} <: AbsCoords where {T<:Real}
    s::T
    p::T
end 

function pb_coords(curve::C, pt::SVector{2,T}, velocity::SVector{2,T}) where {C<:AbsCurve, T<:Real}
    s = arc_length(curve, pt)
    g = domain_gradient_vector(curve, pt)
    n =  g./norm(g)
    v = velocity ./ norm(velocity) 
    p = sin(angle(v, n))
    return PoincareBirkhoff(s,p)
end