"""
    PoincareBirkhoff{T}<:AbsCoords where {T<:Real}

Poincaré-Birkhoff boundary coordinates.

The coordinate `s` gives the arclength position along the boundary and `p`
stores the momentum-like angular coordinate.

## Fields
* `s::T`: Arclength coordinate along the boundary.
* `p::T`: Poincaré-Birkhoff momentum coordinate.
"""
struct PoincareBirkhoff{T}<:AbsCoords where {T<:Real}
    s::T
    p::T
end

"""
    fundamental_pb_coords(L::T,curve::C,pt::SVector{2,T},velocity::SVector{2,T}) where {C<:AbsCurve,T<:Real}

Compute Poincaré-Birkhoff coordinates on a curve of the fundamental boundary.

## Arguments
* `L::T`: Arclength accumulated before `curve`.
* `curve::C`: Boundary curve containing `pt`.
* `pt::SVector{2,T}`: Boundary point.
* `velocity::SVector{2,T}`: Velocity at the boundary point.

## Returns
* Tuple `(s,p)` containing the boundary arclength coordinate and momentum coordinate.
"""
function fundamental_pb_coords(L::T,curve::C,pt::SVector{2,T},velocity::SVector{2,T}) where {C<:AbsCurve,T<:Real}
    s=arc_length(curve,pt)
    g=domain_gradient_vector(curve,pt)
    n=g./norm(g)
    v=velocity./norm(velocity)
    p=sin(signed_angle(v,n))
    return L+s,p
end

"""
    get_pb_curve(composite_curve::C,segment_id,domain_id) where C<:AbsCompositeCurve

Find the subcurve of a composite boundary identified by `domain_id` and `segment_id`.

## Arguments
* `composite_curve::C`: Composite boundary curve.
* `segment_id`: Segment identifier.
* `domain_id`: Domain identifier.

## Returns
* Tuple containing the corresponding entry of `composite_curve.end_lengths` and
  the matching subcurve, if one is found.
"""
function get_pb_curve(composite_curve::C,segment_id,domain_id) where C<:AbsCompositeCurve
    curves=composite_curve.subcurves
    for (i,crv) in enumerate(curves)
        if domain_id==crv.domain_id&&segment_id==crv.segment_id
            return composite_curve.end_lengths[i],crv
        end
    end
end

function pb_coords(billiard::B,subsegment::Int64,subdomain::Int64,sym_sector::Int64,pt::SVector{2,T},velocity::SVector{2,T}) where {B<:AbsBilliard,T<:Real}
    symmetries=billiard.symmetries
    fundamental_boundary=CompositeCurve(get_boundary_curves(billiard))
    L=fundamental_boundary.length
    l,crv=get_pb_curve(fundamental_boundary,subsegment,subdomain)
    s,p=fundamental_pb_coords(l,crv,pt,velocity) #l is length of prevoius curves
    if sym_sector>1
        s,p=apply_symmetry_pb(symmetries[sym_sector-1],sym_sector,s,p,L)
    end
    return PoincareBirkhoff(s,p)
end

function pb_sectors(billiard)
    boundary=CompositeCurve(get_boundary_curves(billiard))
    L=boundary.length
    ends=boundary.end_lengths
    all_sectors=copy(boundary.end_lengths)
    for sym_sector in 2:(length(billiard.symmetries)+1)
        println(sym_sector)
        sym_ends=[apply_symmetry_pb(billiard.symmetries[sym_sector-1],sym_sector,s,0.0,L)[1] for s in ends]
        append!(all_sectors,sym_ends)
    end
    return find_unique_elements(sort(all_sectors))
end