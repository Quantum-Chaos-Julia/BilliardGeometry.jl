"""
    JuliaBilliard{T} <: AbsBilliard

Disc of radius `R` with three congruent circular obstacles of radius `r`, centered at distance `d` from the origin at the polar angles `π/2`, `7π/6` and `11π/6`, as in the Julia logo.

## Description
The shape has the dihedral symmetry `D₃`: the rotations by `2π/3` and the three mirror lines
at the angles `π/6`, `π/2` and `5π/6`, each of which passes through one obstacle. The outer
region is the sector `π/6 ≤ φ ≤ π/2` of the disc, bounded by one physical arc and two
[`SymmetryWall`](@ref) radii (`domain_id = 1`), see [`Dn_symmetry`](@ref)`(3; shift = 1//6)`.
The domain is multiply connected with genus `3`. The obstacles form one
[`SymmetricHole`](@ref) acted on by the same group about the origin: its piece is the
half of the disc at `(0, d)` with `x ≥ 0` (`domain_id = 2`), closed by a [`SymmetryWall`](@ref) on
the mirror line, and its six images are the halves of the three obstacles.

## Attributes
* `fundamental_domain::MultiplyConnectedDomain{T}`: The sector of the disc with the three circular holes.
* `symmetries::SymmetryRegistry`: The five non-identity elements of [`Dn_symmetry`](@ref)`(3; shift = 1//6)`.
"""
struct JuliaBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::MultiplyConnectedDomain{T}
    symmetries::SymmetryRegistry
end

"""
    JuliaBilliard(r, d; R = 1) → billiard::JuliaBilliard

Construct a [`JuliaBilliard`](@ref).

## Arguments
* `r`: Radius of each obstacle, must be positive.
* `d`: Distance of the obstacle centers from the origin, must be positive.

## Keyword arguments
* `R = 1` : Radius of the outer disc.

## Throws
`ArgumentError` for non-positive `r`, `d` or `R`, when an obstacle does not lie strictly inside the disc (`d + r ≥ R`), or when two obstacles overlap or touch (`√3 d ≤ 2r`).
"""
function JuliaBilliard(r::T, d::T; R=one(T)) where T<:Real
    R = T(R)
    r > zero(T) || throw(ArgumentError("r must be positive; received $r"))
    d > zero(T) || throw(ArgumentError("d must be positive; received $d"))
    R > zero(T) || throw(ArgumentError("R must be positive; received $R"))
    d + r < R || throw(ArgumentError("obstacles must lie strictly inside the disc, need d + r < R; received d=$d, r=$r, R=$R"))
    sqrt(T(3))*d > 2*r || throw(ArgumentError("obstacles must not overlap or touch, need sqrt(3)*d > 2r; received d=$d, r=$r"))
    bc = SpecularReflection()
    c = SVector{2,T}(zero(T), zero(T))
    symmetries = Dn_symmetry(3; shift=1//6, T=T) # sym_id 1 = mirror at π/2, sym_id 5 = mirror at π/6
    arc = CircleSegment(R, T(pi)/3, T(pi)/6, c; bc=bc, orientation=1, domain_id=1, segment_id=1)
    p0 = curve(arc, zero(T))
    p1 = curve(arc, one(T))
    wall1 = LineSegment(p1, c; bc=SymmetryWall(1,2), orientation=1, domain_id=1, segment_id=2)
    wall0 = LineSegment(c, p0; bc=SymmetryWall(5,2), orientation=1, domain_id=1, segment_id=3)
    outer_domain = SimpleDomain{T}(AbsCurve[arc, wall1, wall0], SVector{2,T}[p0, p1, c], 1)
    hole_center = SVector{2,T}(zero(T), d)
    top = hole_center + SVector{2,T}(zero(T), r)
    bottom = hole_center - SVector{2,T}(zero(T), r)
    half_arc = CircleSegment(r, T(pi), -T(pi)/2, hole_center; bc=bc, orientation=1, domain_id=2, segment_id=1)
    mirror = LineSegment(top, bottom; bc=SymmetryWall(1,2), orientation=1, domain_id=2, segment_id=2)
    piece = SimpleDomain{T}(AbsCurve[half_arc, mirror], SVector{2,T}[top, bottom], 2)
    holes = SymmetricHole(piece, Dn_symmetry(3; shift=1//6, T=T), c)
    fundamental_domain = MultiplyConnectedDomain(outer_domain, [holes], 1)
    validate_holes(fundamental_domain, symmetries)
    billiard = JuliaBilliard{T}(fundamental_domain, symmetries)
    _log_symmetry_group("JuliaBilliard", billiard)
    return billiard
end
