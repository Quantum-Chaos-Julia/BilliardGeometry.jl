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

# `true` for the pure (orientation-reversing) reflections, `false` for
# symmetries that preserve boundary traversal direction (π-rotations,
# N-fold rotations).
_orientation_reversing(::XAxisReflection) = true
_orientation_reversing(::YAxisReflection) = true
_orientation_reversing(::DiagonalReflection) = true
_orientation_reversing(::AntiDiagonalReflection) = true
_orientation_reversing(::XYAxisReflection) = false
_orientation_reversing(::NFoldRotation) = false

# 2x2 linear-map matrix implementing `apply_symmetry(sym, ·)`, reused to
# transform a `CircleSegment`'s center and angular parametrization.
_sym_matrix(::XAxisReflection) = reflect_y.linear
_sym_matrix(::YAxisReflection) = reflect_x.linear
_sym_matrix(::XYAxisReflection) = reflect_xy.linear
_sym_matrix(::DiagonalReflection) = reflect_diag.linear
_sym_matrix(::AntiDiagonalReflection) = reflect_antidiag.linear
_sym_matrix(sym::NFoldRotation) = sym.sym_map.linear

# Pointwise image of a `LineSegment` under `sym` (endpoints only). The
# `orientation` field is negated since swapping the roles of "inside"/
# "outside" is exactly what a reflected/rotated line implies for the raw
# `line_domain` formula; `full_boundary` curves are not used for `is_inside`
# checks, so this is a minor correctness nicety rather than a requirement.
function _apply_symmetry_to_curve(sym::AbsSymmetry, c::LineSegment{T}) where {T<:Real}
    pt0 = SVector{2,T}(apply_symmetry(sym, c.pt0))
    pt1 = SVector{2,T}(apply_symmetry(sym, c.pt1))
    return LineSegment(pt0, pt1; bc=c.bc, orientation=-c.orientation, domain_id=c.domain_id, segment_id=c.segment_id)
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
    M = _sym_matrix(sym)
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

# Reverses a curve's parametrization direction (`t -> 1-t`), leaving its
# physical trace unchanged.
_reverse_curve(c::LineSegment) = LineSegment(c.pt1, c.pt0; bc=c.bc, orientation=-c.orientation, domain_id=c.domain_id, segment_id=c.segment_id)
function _reverse_curve(c::CircleSegment{T}) where {T<:Real}
    new_shift = c.shift_angle + c.arc_angle
    new_arc = -c.arc_angle
    return CircleSegment(c.radius, new_arc, new_shift, c.center; bc=c.bc, orientation=c.orientation, domain_id=c.domain_id, segment_id=c.segment_id)
end

"""
    full_boundary(billiard::Bi) where {Bi<:AbsBilliard} → curves::Vector{AbsCurve}

Reconstructs the complete, closed, canonically CCW-oriented physical boundary
of `billiard` from its fundamental domain and its discrete `symmetries`.

## Description
Starts from `get_boundary_curves(billiard)` — the connected `SpecularReflection`
curves of the fundamental domain, i.e. the physical part of the discretization
basis solvers already use — and appends, for every `sym in billiard.symmetries`,
the image of those same fundamental-domain curves under `sym`. Symmetry walls
(curves whose boundary condition is [`ReflectionSymmetry`](@ref)) and internal
subdomain seams (e.g. `Transparent`) are never included, matching
[`get_boundary_curves`](@ref)'s existing filtering. For a billiard with no
symmetries, `full_boundary(billiard) == get_boundary_curves(billiard)`.

Currently supports [`LineSegment`](@ref) and [`CircleSegment`](@ref) physical
curves under [`XAxisReflection`](@ref), [`YAxisReflection`](@ref),
[`XYAxisReflection`](@ref), [`DiagonalReflection`](@ref),
[`AntiDiagonalReflection`](@ref) and [`NFoldRotation`](@ref) symmetries;
generalizing to every curve/segment type is deferred until more billiards are
ported.

## Arguments
* `billiard`: The billiard whose complete physical boundary is reconstructed.

## Returns
* `curves`: The complete, closed physical boundary as a `Vector{AbsCurve}`.
"""
function full_boundary(billiard::Bi) where {Bi<:AbsBilliard}
    P = get_boundary_curves(billiard)
    curves = Vector{AbsCurve}(P)
    for sym in billiard.symmetries
        if _orientation_reversing(sym)
            images = [_reverse_curve(_apply_symmetry_to_curve(sym, c)) for c in P]
            append!(curves, reverse(images))
        else
            images = [_apply_symmetry_to_curve(sym, c) for c in P]
            append!(curves, images)
        end
    end
    return curves
end
