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

function get_pb_curve(composite_curve::C, domain_id, segment_id) where C<:AbsCompositeCurve
    curves = composite_curve.subcurves
    for (i,crv) in enumerate(curves)
        if (domain_id == crv.domain_id && segment_id == crv.segment_id)
            return composite_curve.end_lengths[i], crv
        end
    end
    return nothing, nothing
end