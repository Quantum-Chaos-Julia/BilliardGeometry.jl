"""
    AnnularBilliard{T} <: AbsBilliard

Circular billiard of radius `R_outer` with a circular hole (obstacle) of radius
`R_inner` whose center is displaced along the `x` axis.

## Description
The domain is multiply connected with genus `1`. With the default `center = (0,0)`
the billiard has the mirror symmetry `y → -y` of [`XAxisReflection`](@ref) (the
hole center lies on the `x` axis): the outer region is the upper half of the
disk (a [`CircleSegment`](@ref) arc and a [`SymmetryWall`](@ref) edge,
`domain_id = 1`) and the hole is a [`SymmetricHole`](@ref) whose piece is the
upper half of the hole disk (`domain_id = 2`). With another `center` no symmetry
is used: the outer region is a [`SimpleDomain`](@ref) holding a full
counterclockwise [`CircleSegment`](@ref) with `domain_id = 1` and the hole is a
[`DomainHole`](@ref) with a full counterclockwise circle with `domain_id = 2`,
i.e. traversed with the hole region on its left. The case `shift = 0` is the
concentric annulus.

## Attributes
* `fundamental_domain::MultiplyConnectedDomain{T}`: Outer disk (upper half for `center = (0,0)`) with the hole.
* `symmetries::SymmetryRegistry`: Registry of [`XAxisReflection`](@ref) for `center = (0,0)`, empty otherwise.
"""
struct AnnularBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::MultiplyConnectedDomain{T}
    symmetries::SymmetryRegistry
end

"""
    AnnularBilliard(R_inner, shift; R_outer = 1, center = (0,0)) → billiard::AnnularBilliard

Construct an [`AnnularBilliard`](@ref).

## Arguments
* `R_inner`: Radius of the hole, must satisfy `0 < R_inner < R_outer`.
* `shift`: Displacement of the hole center from the outer center along the `x` axis.

## Keyword arguments
* `R_outer = 1` : Radius of the outer circle.
* `center` : Center of the outer circle. The `x` axis mirror symmetry is used only for the origin.

## Throws
`ArgumentError` for non-positive radii, `R_inner ≥ R_outer`, or a hole that is not inside the outer disk (see [`validate_holes`](@ref)).
"""
function AnnularBilliard(R_inner::T, shift::T; R_outer = one(T), center=SVector{2,T}(zero(T),zero(T))) where T<:Real
    R_outer > zero(T) || throw(ArgumentError("R_outer must be positive; received $R_outer"))
    R_inner > zero(T) || throw(ArgumentError("R_inner must be positive; received $R_inner"))
    R_inner < R_outer || throw(ArgumentError("R_inner must be smaller than R_outer; received R_inner=$R_inner, R_outer=$R_outer"))
    c = SVector{2,T}(center)
    bc = SpecularReflection()
    hole_center = c+SVector{2,T}(shift,zero(T))
    if iszero(c)
        arc = CircleSegment(R_outer, T(pi), zero(T), c; bc=bc, orientation=1, domain_id=1, segment_id=1)
        left = SVector{2,T}(-R_outer,zero(T)); right = SVector{2,T}(R_outer,zero(T))
        wall = LineSegment(left, right; bc=SymmetryWall(1,2), orientation=1, domain_id=1, segment_id=2)
        outer_domain = SimpleDomain{T}(AbsCurve[arc,wall], SVector{2,T}[left,right], 1)
        piece = CircleWedge(R_inner, T(pi), zero(T), shift, zero(T), 2; bcs=AbsBoundaryCondition[bc, SymmetryWall(1,2), SymmetryWall(1,2)])
        symmetries = register_symmetries(XAxisReflection())
        fundamental_domain = MultiplyConnectedDomain(outer_domain, [SymmetricHole(piece, register_symmetries(XAxisReflection()), hole_center)], 1)
    else
        outer = CircleSegment(R_outer, T(2*pi), zero(T), c; bc=bc, orientation=1, domain_id=1, segment_id=1)
        inner = CircleSegment(R_inner, T(2*pi), zero(T), hole_center; bc=bc, orientation=1, domain_id=2, segment_id=1)
        outer_domain = SimpleDomain{T}(AbsCurve[outer], SVector{2,T}[c+SVector{2,T}(R_outer,zero(T))], 1)
        hole_domain = SimpleDomain{T}(AbsCurve[inner], SVector{2,T}[hole_center+SVector{2,T}(R_inner,zero(T))], 2)
        symmetries = register_symmetries()
        fundamental_domain = MultiplyConnectedDomain(outer_domain, [DomainHole(hole_domain)], 1)
    end
    validate_holes(fundamental_domain, symmetries)
    billiard = AnnularBilliard{T}(fundamental_domain, symmetries)
    
    return billiard
end

