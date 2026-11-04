"""
    SinaiBilliard{T} <: AbsBilliard

Sinai billiard: the square `[-a,a] × [-a,a]` with a concentric circular
obstacle of radius `R`.

## Description
The domain is multiply connected with genus `1`. With the default `center = (0,0)` the
billiard has the `D₄` symmetry of [`D4_symmetry`](@ref) and the outer region is
the fundamental triangle `0 ≤ y ≤ x ≤ a` (one physical side and two
[`SymmetryWall`](@ref) edges, `domain_id = 1`); the obstacle is a
[`SymmetricHole`](@ref) whose piece is the eighth of the disc of radius `R`
between the angles `0` and `π/4` (`domain_id = 2`), so the symmetry group of
the billiard is the full `D₄`. With another `center` no symmetry is used: the
outer region is the full square, a [`SimpleDomain`](@ref) of the four sides
traversed counterclockwise with `domain_id = 1`, and the obstacle is a
[`DomainHole`](@ref) with a full counterclockwise [`CircleSegment`](@ref) with
`domain_id = 2`, i.e. traversed with the hole region on its left.

## Attributes
* `fundamental_domain::MultiplyConnectedDomain{T}`: Fundamental triangle (or full square) with the circular hole.
* `symmetries::SymmetryRegistry`: [`D4_symmetry`](@ref) for `center = (0,0)`, empty otherwise.
"""
struct SinaiBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::MultiplyConnectedDomain{T}
    symmetries::SymmetryRegistry
end

"""
    SinaiBilliard(a, R; center = (0,0)) → billiard::SinaiBilliard

Construct a [`SinaiBilliard`](@ref).

## Arguments
* `a`: Half side length of the square.
* `R`: Radius of the obstacle, `0 < R < a`.

## Keyword arguments
* `center` : Center of the square and of the obstacle. The `D₄` symmetry reduction is used only for the origin.

## Throws
`ArgumentError` for non-positive `a` or `R`, or `R ≥ a`.
"""
function SinaiBilliard(a::T, R::T; center=SVector{2,T}(zero(T),zero(T))) where T<:Real
    a > zero(T) || throw(ArgumentError("a must be positive; received $a"))
    R > zero(T) || throw(ArgumentError("R must be positive; received $R"))
    R < a || throw(ArgumentError("R must be smaller than a; received a=$a, R=$R"))
    c = SVector{2,T}(center)
    bc = SpecularReflection()
    if iszero(c)
        q0 = SVector{2,T}(zero(T),zero(T)); q1 = SVector{2,T}(a,zero(T)); q2 = SVector{2,T}(a,a)
        side = LineSegment(q1,q2; bc=bc, orientation=1, domain_id=1, segment_id=1)
        diagonal_wall = LineSegment(q2,q0; bc=SymmetryWall(1,2), orientation=1, domain_id=1, segment_id=2)
        axis_wall = LineSegment(q0,q1; bc=SymmetryWall(7,8), orientation=1, domain_id=1, segment_id=3)
        outer_domain = SimpleDomain{T}(AbsCurve[side,diagonal_wall,axis_wall], SVector{2,T}[q1,q2,q0], 1)
        piece = CircleWedge(R, T(pi/4), zero(T), zero(T), zero(T), 2; bcs=AbsBoundaryCondition[bc, SymmetryWall(1,2), SymmetryWall(7,8)])
        symmetries = D4_symmetry()
        fundamental_domain = MultiplyConnectedDomain(outer_domain, [SymmetricHole(piece, D4_symmetry(), SVector{2,T}(zero(T),zero(T)))], 1)
    else
        p0 = c+SVector{2,T}(-a,-a); p1 = c+SVector{2,T}(a,-a); p2 = c+SVector{2,T}(a,a); p3 = c+SVector{2,T}(-a,a)
        bottom = LineSegment(p0,p1; bc=bc, orientation=1, domain_id=1, segment_id=1)
        right = LineSegment(p1,p2; bc=bc, orientation=1, domain_id=1, segment_id=2)
        top = LineSegment(p2,p3; bc=bc, orientation=1, domain_id=1, segment_id=3)
        left = LineSegment(p3,p0; bc=bc, orientation=1, domain_id=1, segment_id=4)
        inner = CircleSegment(R, T(2*pi), zero(T), c; bc=bc, orientation=1, domain_id=2, segment_id=1)
        outer_domain = SimpleDomain{T}(AbsCurve[bottom,right,top,left], SVector{2,T}[p1,p2,p3,p0], 1)
        hole_domain = SimpleDomain{T}(AbsCurve[inner], SVector{2,T}[c+SVector{2,T}(R,zero(T))], 2)
        symmetries = register_symmetries()
        fundamental_domain = MultiplyConnectedDomain(outer_domain, [DomainHole(hole_domain)], 1)
    end
    validate_holes(fundamental_domain, symmetries)
    billiard = SinaiBilliard{T}(fundamental_domain, symmetries)
    return billiard
end