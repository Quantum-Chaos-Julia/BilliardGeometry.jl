"""
    SquareWithinSquareBilliard{T} <: AbsBilliard

Square `[-a,a] × [-a,a]` containing axis-aligned square obstacles centered on
the diagonal `y = x`.

## Description
The domain is multiply connected with genus equal to the number of obstacles.
No symmetry reduction is applied, so the fundamental domain is the full
physical geometry. The outer region is a [`SimpleDomain`](@ref) holding the four
sides of the outer square traversed counterclockwise with `domain_id = 1`; each
obstacle is a [`DomainHole`](@ref) with four counterclockwise
[`LineSegment`](@ref)s and `domain_id = i + 1`, `i` counting the obstacles from
the lower left, i.e. traversed with the hole region on its left.

## Attributes
* `fundamental_domain::MultiplyConnectedDomain{T}`: Outer square with the square holes.
* `symmetries::SymmetryRegistry`: Empty registry.
"""
struct SquareWithinSquareBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::MultiplyConnectedDomain{T}
    symmetries::SymmetryRegistry
end

"""
    SquareWithinSquareBilliard(a = 1; positions = [-a/2, a/2], sizes = fill(a/4, 2)) → billiard::SquareWithinSquareBilliard

Construct a [`SquareWithinSquareBilliard`](@ref).

## Arguments
* `a`: Half side length of the outer square, must be positive.

## Keyword arguments
* `positions` : Coordinate `d` of the center `(d, d)` of each obstacle on the diagonal.
* `sizes` : Side length of each obstacle, positive and of the same length as `positions`.

## Throws
`ArgumentError` if the lengths differ, a size is not positive, an obstacle does not lie strictly inside the outer square, or two obstacles overlap or touch.
"""
function SquareWithinSquareBilliard(a::T=one(T); positions::AbstractVector{<:Real}=T[-a/2,a/2], sizes::AbstractVector{<:Real}=fill(a/4,2)) where {T<:Real}
    a > zero(T) || throw(ArgumentError("a must be positive; received $a"))
    length(positions) == length(sizes) || throw(ArgumentError("positions and sizes must have the same length"))
    pos = T.(positions); sz = T.(sizes); n = length(pos)
    all(>(zero(T)), sz) || throw(ArgumentError("All square side lengths must be positive"))
    p = sortperm(pos); pos = pos[p]; sz = sz[p]
    @inbounds for i in 1:n
        h = sz[i]/2
        -a < pos[i]-h && pos[i]+h < a || throw(ArgumentError("Square hole $i must lie strictly inside the outer square; received position=$(pos[i]), size=$(sz[i])"))
    end
    @inbounds for i in 1:n-1
        pos[i]+sz[i]/2 < pos[i+1]-sz[i+1]/2 || throw(ArgumentError("Square holes $i and $(i+1) overlap or touch"))
    end
    bc = SpecularReflection()
    outer_domain = _square_domain(SVector{2,T}(zero(T),zero(T)), a, bc, 1)
    holes = DomainHole[DomainHole(_square_domain(SVector{2,T}(pos[i],pos[i]), sz[i]/2, bc, i+1)) for i in 1:n]
    symmetries = register_symmetries()
    fundamental_domain = MultiplyConnectedDomain(outer_domain, holes, 1)
    validate_holes(fundamental_domain, symmetries)
    return SquareWithinSquareBilliard{T}(fundamental_domain, symmetries)
end

# Counterclockwise square of half side `h` about `c` as a closed SimpleDomain.
function _square_domain(c::SVector{2,T}, h::T, bc::AbsBoundaryCondition, domain_id::Integer) where {T<:Real}
    p0 = c+SVector{2,T}(-h,-h); p1 = c+SVector{2,T}(h,-h); p2 = c+SVector{2,T}(h,h); p3 = c+SVector{2,T}(-h,h)
    sides = AbsCurve[
        LineSegment(p0,p1; bc=bc, orientation=1, domain_id=domain_id, segment_id=1),
        LineSegment(p1,p2; bc=bc, orientation=1, domain_id=domain_id, segment_id=2),
        LineSegment(p2,p3; bc=bc, orientation=1, domain_id=domain_id, segment_id=3),
        LineSegment(p3,p0; bc=bc, orientation=1, domain_id=domain_id, segment_id=4),
    ]
    return SimpleDomain{T}(sides, SVector{2,T}[p1,p2,p3,p0], domain_id)
end