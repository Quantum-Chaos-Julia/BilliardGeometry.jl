################################################################################
############################# FULL PHYSICAL BOUNDARY ##########################
################################################################################
# `full_boundary(billiard)` reconstructs the complete physical boundary of a
# billiard from its fundamental domain plus its discrete symmetry generators
# (`billiard.symmetries`), instead of every billiard hand-maintaining its own
# `full_boundary` field. This is needed by boundary-integral (BIM) solvers,
# which discretize the *entire* physical boundary (unlike basis solvers,
# which only ever need the fundamental domain).
#
# `billiard.symmetries` lists every non-identity element of the symmetry group
# directly (e.g. `D2_symmetry = [YAxisReflection(),XYAxisReflection(),XAxisReflection()]`),
# not a minimal generating set, so reconstructing the full boundary does not
# require group closure: for each `sym in billiard.symmetries`, the image of
# the fundamental-domain's physical curves under `sym` is appended directly.
#
# Orientation-reversing symmetries (pure reflections: `XAxisReflection`,
# `YAxisReflection`, `DiagonalReflection`, `AntiDiagonalReflection`) reverse
# the sense of a continuous CCW boundary traversal, so their image curves must
# also be reversed (both individually and in list order) to continue the
# closed loop. Orientation-preserving symmetries (`XYAxisReflection` — the
# product of two reflections, i.e. a π rotation — and `NFoldRotation`)
# preserve traversal direction, so their images are appended in the same
# order with no per-curve reversal.
################################################################################

"""
    orientation_reversing(sym::AbsSymmetry) → flag::Bool

`true` for the pure reflections, which reverse the sense of boundary traversal,
and `false` for the symmetries that preserve it: [`XYAxisReflection`](@ref) (a
π rotation) and [`NFoldRotation`](@ref).
"""
orientation_reversing(::XAxisReflection) = true
orientation_reversing(::YAxisReflection) = true
orientation_reversing(::DiagonalReflection) = true
orientation_reversing(::AntiDiagonalReflection) = true
orientation_reversing(::XYAxisReflection) = false
orientation_reversing(::NFoldRotation) = false

# 2x2 linear-map matrix implementing `apply_symmetry(sym, ·)` for a rotation,
# reused to transform a `CircleSegment`'s center and angular parametrization
# (the reflection cases are defined alongside `apply_symmetry` in
# `symmetry.jl`, included earlier).
_sym_matrix(sym::NFoldRotation, ::Type{T}=Float64) where {T<:Real} = SMatrix{2,2,T}(sym.sym_map.linear)

# Pointwise image of a `LineSegment` under `sym` (endpoints only). `line_domain` is the cross product of the line
# direction with the offset of the point, which a linear map `M` scales by `det(M)`; the `orientation` field is therefore
# kept for rotations and flipped for reflections, so that `domain_fun` of the image at `sym(p)` equals `domain_fun` of
# the source at `p`.
function _apply_symmetry_to_curve(sym::AbsSymmetry, c::LineSegment{T}) where {T<:Real}
    pt0 = SVector{2,T}(apply_symmetry(sym, c.pt0))
    pt1 = SVector{2,T}(apply_symmetry(sym, c.pt1))
    orientation = orientation_reversing(sym) ? -c.orientation : c.orientation
    return LineSegment(pt0, pt1; bc=c.bc, orientation=orientation, domain_id=c.domain_id, segment_id=c.segment_id)
end

# Pointwise image of a `CircleSegment` under `sym`. Writing the 2x2 matrix
# `M=_sym_matrix(sym)` as a rotation (det M=+1, rotation angle φ) or a
# reflection (det M=-1, reflection axis angle ψ with M[2,1]=sin(2ψ),
# M[1,1]=cos(2ψ)), the image of `circle_eq(R,arc_angle,shift_angle,center,t)`
# under `M` is again a circular arc of the same radius, with
#   center' = M*center,
#   (shift_angle', arc_angle') = (shift_angle+φ, arc_angle)              if det M=+1,
#                                (2ψ-shift_angle, -arc_angle)            if det M=-1.
# A negative `arc_angle'` (reflection case) is resolved by `_reverse_curve`,
# which reverses the parametrization direction back to a positive arc angle.
function _apply_symmetry_to_curve(sym::AbsSymmetry, c::CircleSegment{T}) where {T<:Real}
    M = _sym_matrix(sym, T)
    center2 = SVector{2,T}(M*c.center)
    if det(M) > 0
        φ = atan(M[2,1], M[1,1])
        new_shift = c.shift_angle + φ
        new_arc = c.arc_angle
    else
        two_ψ = atan(M[2,1], M[1,1])
        new_shift = two_ψ - c.shift_angle
        new_arc = -c.arc_angle
    end
    return CircleSegment(c.radius, new_arc, new_shift, center2; bc=c.bc, orientation=c.orientation, domain_id=c.domain_id, segment_id=c.segment_id)
end

# Center/shift-angle/arc-angle image of a polar curve under `sym`, shared by
# both `FourierCoeffPolarSegment` and `PolarSegment`. Same rotation/reflection
# formula as `CircleSegment` above, since it only depends on the angular
# parametrization convention (rotating a polar curve about a possibly
# off-origin center rotates its local origin too, not just its angular
# phase); the radial function (`coef`/`r_func`) is invariant under rotation
# of the ambient frame and is carried through unchanged by the caller.
function _polar_symmetry_image(sym::AbsSymmetry, c::L) where {L<:AbsPolarCurve}
    T = eltype(c.center)
    M = _sym_matrix(sym, T)
    center2 = SVector{2,T}(M*c.center)
    if det(M) > 0
        φ = atan(M[2,1], M[1,1])
        new_shift = c.shift_angle + φ
        new_arc = c.arc_angle
    else
        two_ψ = atan(M[2,1], M[1,1])
        new_shift = two_ψ - c.shift_angle
        new_arc = -c.arc_angle
    end
    return center2, new_shift, new_arc
end

# Pointwise image of a `FourierCoeffPolarSegment`/`PolarSegment` under `sym`.
# A negative `new_arc` (reflection case) is resolved by `_reverse_curve`,
# which reverses the parametrization direction back to a positive arc angle.
function _apply_symmetry_to_curve(sym::AbsSymmetry, c::L) where {L<:FourierCoeffPolarSegment}
    center2, new_shift, new_arc = _polar_symmetry_image(sym, c)
    return FourierCoeffPolarSegment(c.coef; R=c.R, arc_angle=new_arc, shift_angle=new_shift, center=center2, orientation=c.orientation, bc=c.bc, domain_id=c.domain_id, segment_id=c.segment_id)
end
function _apply_symmetry_to_curve(sym::AbsSymmetry, c::L) where {L<:PolarSegment}
    center2, new_shift, new_arc = _polar_symmetry_image(sym, c)
    return PolarSegment(c.r_func; R=c.R, arc_angle=new_arc, shift_angle=new_shift, center=center2, orientation=c.orientation, bc=c.bc, domain_id=c.domain_id, segment_id=c.segment_id)
end

# Reverses a curve's parametrization direction (`t -> 1-t`), leaving its
# physical trace unchanged.
_reverse_curve(c::LineSegment) = LineSegment(c.pt1, c.pt0; bc=c.bc, orientation=-c.orientation, domain_id=c.domain_id, segment_id=c.segment_id)
function _reverse_curve(c::CircleSegment{T}) where {T<:Real}
    new_shift = c.shift_angle + c.arc_angle
    new_arc = -c.arc_angle
    return CircleSegment(c.radius, new_arc, new_shift, c.center; bc=c.bc, orientation=c.orientation, domain_id=c.domain_id, segment_id=c.segment_id)
end
function _reverse_curve(c::L) where {L<:FourierCoeffPolarSegment}
    new_shift = c.shift_angle + c.arc_angle
    new_arc = -c.arc_angle
    return FourierCoeffPolarSegment(c.coef; R=c.R, arc_angle=new_arc, shift_angle=new_shift, center=c.center, orientation=c.orientation, bc=c.bc, domain_id=c.domain_id, segment_id=c.segment_id)
end
function _reverse_curve(c::L) where {L<:PolarSegment}
    new_shift = c.shift_angle + c.arc_angle
    new_arc = -c.arc_angle
    return PolarSegment(c.r_func; R=c.R, arc_angle=new_arc, shift_angle=new_shift, center=c.center, orientation=c.orientation, bc=c.bc, domain_id=c.domain_id, segment_id=c.segment_id)
end

"""
    _reconstruct_curves(P, symmetries)

Reconstruct a closed curve set from its fundamental curves `P` and a
symmetry registry, together with its exact discrete construction provenance.

## Description
The fundamental curves `P` form reconstruction sector `1`. For every
non-identity symmetry `symmetries[k]`, sector `k + 1` is constructed directly
from the same fundamental curves.

For every reconstructed curve, `source` stores the corresponding fundamental
curve index, `sector` stores the reconstruction-sector index, and `reversed`
records whether the fundamental curve parametrization was reversed when
constructing that image.

The provenance is derived directly from the reconstruction procedure. No
geometric matching or floating-point inference is performed.

## Arguments
* `P`: Fundamental curves in traversal order.
* `symmetries`: Registry of non-identity symmetry elements.

## Returns
* `curves::Vector{AbsCurve}`: Reconstructed curves.
* `source::Vector{Int}`: Fundamental-curve source index of each curve.
* `sector::Vector{Int}`: Reconstruction-sector index of each curve.
* `reversed::BitVector`: Whether each reconstructed curve reverses its source parametrization.
"""
function _reconstruct_curves(P::AbstractVector{<:AbsCurve}, symmetries::SymmetryRegistry)
    np = length(P)
    curves = Vector{AbsCurve}(P)
    source = collect(1:np)
    sector = ones(Int, np)
    reversed = falses(np)
    @inbounds for (k, sym) in enumerate(symmetries)
        if orientation_reversing(sym)
            images = [_reverse_curve(_apply_symmetry_to_curve(sym, c)) for c in P]
            append!(curves, reverse(images))
            append!(source, np:-1:1)
            append!(sector, fill(k + 1, np))
            append!(reversed, trues(np))
        else
            images = [_apply_symmetry_to_curve(sym, c) for c in P]
            append!(curves, images)
            append!(source, 1:np)
            append!(sector, fill(k + 1, np))
            append!(reversed, falses(np))
        end
    end
    return curves, source, sector, reversed
end

"""
    _full_boundary_data(billiard::Bi) where {Bi<:AbsBilliard}

Reconstruct the complete outer physical boundary of `billiard` together with
its discrete construction provenance, see [`_reconstruct_curves`](@ref).

For a billiard with a [`MultiplyConnectedDomain`](@ref) only the outer
boundary is reconstructed; hole loops are handled by
[`full_boundary_components`](@ref).
"""
function _full_boundary_data(billiard::Bi) where {Bi<:AbsBilliard}
    return _reconstruct_curves(get_outer_boundary_curves(billiard), billiard.symmetries)
end

"""
    full_boundary(billiard::Bi) where {Bi<:AbsBilliard}

Reconstruct the complete physical boundary from the billiard's fundamental
physical curves and registered discrete symmetries.

## Description
Sector `1` is the fundamental physical boundary returned by
`get_boundary_curves(billiard)`. Each registered non-identity symmetry adds
one additional physical reconstruction sector. Reflection images are reversed
as required to preserve continuous boundary traversal orientation.

## Arguments
* `billiard::Bi`: Billiard whose complete physical boundary is requested.

## Returns
* `curves::Vector{AbsCurve}`: Complete reconstructed physical boundary.

For a billiard with holes this is the concatenation of the loops of
[`full_boundary_components`](@ref).
"""
function full_boundary(billiard::Bi) where {Bi<:AbsBilliard}
    return _full_boundary(billiard.fundamental_domain, billiard)
end

_full_boundary(::AbsDomain, billiard::AbsBilliard) = first(_full_boundary_data(billiard))
_full_boundary(::AbsMultiplyConnectedDomain, billiard::AbsBilliard) = reduce(vcat, full_boundary_components(billiard).curves)

"""
    BoundaryComponentMap

Closed boundary loops of a billiard together with their role.

## Attributes
* `curves::Vector{Vector{AbsCurve}}`: One closed loop per entry; the reconstructed outer boundary first, then every hole loop.
* `is_hole::Vector{Bool}`: Whether the loop is a hole loop.
* `group::Vector{Int}`: Solver group of the loop: `0` for the outer boundary, `i` for loops of the `i`-th hole object.
"""
struct BoundaryComponentMap
    curves::Vector{Vector{AbsCurve}}
    is_hole::Vector{Bool}
    group::Vector{Int}
end

"""
    full_boundary_components(billiard::AbsBilliard) → map::BoundaryComponentMap

Closed physical boundary loops of the complete billiard, derived from the
construction provenance (no coordinate matching). The outer loop is the outer
boundary reconstructed by the billiard symmetry registry; hole loops are taken
from the holes in the full frame, each traversed with the hole region on the
left. For billiards without holes the single result of [`full_boundary`](@ref)
is returned as one component.
"""
function full_boundary_components(billiard::B) where {B<:AbsBilliard}
    return _full_boundary_components(billiard.fundamental_domain, billiard)
end

_full_boundary_components(::AbsDomain, billiard::AbsBilliard) = BoundaryComponentMap([full_boundary(billiard)], [false], [0])

function _full_boundary_components(domain::AbsMultiplyConnectedDomain, billiard::AbsBilliard)
    curves = Vector{Vector{AbsCurve}}()
    push!(curves, first(_full_boundary_data(billiard)))
    is_hole = Bool[false]
    group = Int[0]
    for (i, h) in enumerate(domain.holes), loop in boundary_components(h)
        push!(curves, loop)
        push!(is_hole, true)
        push!(group, i)
    end
    return BoundaryComponentMap(curves, is_hole, group)
end
