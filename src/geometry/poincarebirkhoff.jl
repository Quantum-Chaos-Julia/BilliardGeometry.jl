"""
    PoincareBirkhoff{T} <: AbsCoords

Poincaré–Birkhoff phase space coordinates of a boundary collision.

## Attributes
* `s::T`: Arc length along the boundary, measured from the start of the boundary traversal.
* `p::T`: Sine of the signed angle between the velocity and the outward normal, `p ∈ [-1,1]`.
"""
struct PoincareBirkhoff{T} <: AbsCoords where {T<:Real}
    s::T
    p::T
end 

function fundamental_pb_coords(L::T, curve::C, pt::SVector{2,T}, velocity::SVector{2,T}) where {C<:AbsCurve, T<:Real}
    s = arc_length(curve, pt)
    g = domain_gradient_vector(curve, pt)
    n =  g./norm(g)
    v = velocity ./ norm(velocity) 
    p = sin(signed_angle(v, n))
    return L+s, p
end

"""
    get_pb_curve(composite_curve::AbsCompositeCurve, segment_id, domain_id) → (l, crv)

Find the curve of `composite_curve` with the given `segment_id` and `domain_id`.

## Returns
* `l`: Arc length of `composite_curve` at the start of `crv`.
* `crv`: The curve. The function returns `nothing` if there is no such curve.
"""
function get_pb_curve(composite_curve::C, segment_id, domain_id) where C<:AbsCompositeCurve
    curves = composite_curve.subcurves
    for (i,crv) in enumerate(curves)
        if (domain_id == crv.domain_id && segment_id == crv.segment_id)
            return composite_curve.end_lengths[i], crv
        end
    end
end


"""
    pb_coords(billiard::AbsBilliard, subsegment::Int64, subdomain::Int64, sym_sector::Int64, pt::SVector{2,T}, velocity::SVector{2,T}) → pb::PoincareBirkhoff{T}

Poincaré–Birkhoff coordinates `(s,p)` of a collision with the full physical
boundary of `billiard`.

## Description
The collision point `pt` lies on a curve of the fundamental domain, identified by
`subdomain` and `subsegment`, and belongs to the `sym_sector`-th copy of the
fundamental boundary. The arc length `s` is measured along the connected
fundamental physical boundary and `p` is the sine of the signed angle between
`velocity` and the outward normal, obtained from the gradient of
[`domain_fun`](@ref). For `sym_sector > 1` the coordinates are mapped to the
image sector with [`apply_symmetry_pb`](@ref), consistently with
[`full_boundary`](@ref).

## Arguments
* `billiard`: Billiard.
* `subsegment`: `segment_id` of the curve containing `pt`.
* `subdomain`: `domain_id` of the curve containing `pt`.
* `sym_sector`: Index of the boundary copy, `1` for the fundamental domain.
* `pt`: Collision point in the fundamental frame, on the curve.
* `velocity`: Velocity at the collision, in the fundamental frame.

## Returns
* `pb`: The coordinates as a [`PoincareBirkhoff`](@ref).
"""
function pb_coords(billiard::B, subsegment::Int64, subdomain::Int64, sym_sector::Int64, pt::SVector{2,T}, velocity::SVector{2,T}) where {B<:AbsBilliard, T<:Real}
    symmetries = billiard.symmetries
    fundamental_boundary = CompositeCurve(get_boundary_curves(billiard))
    L = fundamental_boundary.length
    l, crv = get_pb_curve(fundamental_boundary, subsegment, subdomain)
    s, p = fundamental_pb_coords(l, crv, pt, velocity) #l is length of prevoius curves
    if sym_sector > 1 
        sym = symmetry_of(symmetries, sym_sector-1)
        s, p = apply_symmetry_pb(sym, sym_sector, s, p, L)
    end
    return PoincareBirkhoff(s,p)
end

"""
    pb_sectors(billiard::AbsBilliard) → s::Vector

Sorted distinct arc-length positions, along the full outer physical boundary,
at which consecutive boundary curves join, including `0`. The positions are
obtained from the fundamental boundary and its images under every registered
symmetry, see [`apply_symmetry_pb`](@ref).
"""
function pb_sectors(billiard)
    boundary = CompositeCurve(get_outer_boundary_curves(billiard))
    L = boundary.length
    ends = boundary.end_lengths
    all_sectors = copy(boundary.end_lengths)
    for sym_sector in 2:(length(billiard.symmetries)+1)
        sym = symmetry_of(billiard.symmetries, sym_sector-1)
        sym_ends = [apply_symmetry_pb(sym, sym_sector, s, 0.0, L)[1] for s in ends]
        append!(all_sectors, sym_ends)
    end
    return find_unique_elements(sort(all_sectors))
end