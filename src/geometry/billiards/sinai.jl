"""
    SinaiBilliard{T} <: AbsBilliard

Sinai billiard: the square `[-a,a] × [-a,a]` with a concentric circular
obstacle of radius `R`.

## Description
The domain is multiply connected with genus `1`. No symmetry reduction is
applied, so the fundamental domain is the full physical geometry. The outer
region is a [`SimpleDomain`](@ref) holding the four sides of the square
traversed counterclockwise with `domain_id = 1`; the obstacle is a
[`DomainHole`](@ref) with a full counterclockwise [`CircleSegment`](@ref) with
`domain_id = 2`, i.e. traversed with the hole region on its left.

## Attributes
* `fundamental_domain::MultiplyConnectedDomain{T}`: Square with the circular hole.
* `symmetries::SymmetryRegistry`: Empty registry.
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
* `center` : Center of the square and of the obstacle.

## Throws
`ArgumentError` for non-positive `a` or `R`, or `R ≥ a`.
"""
function SinaiBilliard(a::T, R::T; center=SVector{2,T}(zero(T),zero(T))) where T<:Real
    a > zero(T) || throw(ArgumentError("a must be positive; received $a"))
    R > zero(T) || throw(ArgumentError("R must be positive; received $R"))
    R < a || throw(ArgumentError("R must be smaller than a; received a=$a, R=$R"))
    c = SVector{2,T}(center)
    bc = SpecularReflection()
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
    validate_holes(fundamental_domain, symmetries)
    return SinaiBilliard{T}(fundamental_domain, symmetries)
end