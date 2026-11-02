"""
    CutCircleBilliard{T} <: AbsBilliard

Union of two unit disks with centers `(±half_width, 0)`, displaced left-right. The two
circular outlines meet in re-entrant corners at `(0, ±sqrt(1 - half_width²))`. The
fundamental domain is the first quadrant with the `D₂` symmetry group of
[`D2_symmetry`](@ref).

## Attributes
* `fundamental_domain::SimpleDomain{T}`: The quadrant, a circular arc plus two [`SymmetryWall`](@ref) segments.
* `symmetries::SymmetryRegistry`: Reflections `YAxisReflection`, `XYAxisReflection`, `XAxisReflection` with `sym_id` `1, 2, 3`.
"""
struct CutCircleBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::SimpleDomain{T}
    symmetries::SymmetryRegistry
end

"""
    CutCircleBilliard(half_width) → billiard::CutCircleBilliard

Construct a [`CutCircleBilliard`](@ref) whose two unit disks have centers at
`(±half_width, 0)`.

## Arguments
* `half_width`: Half the distance between the centers, `0 < half_width < 1` so that the disks overlap.

## Throws
`ArgumentError` if `half_width` is not in `(0, 1)`.
"""
function CutCircleBilliard(half_width::T) where T<:Real
    (zero(T) < half_width < one(T)) || throw(ArgumentError("half_width must satisfy 0 < half_width < 1; received half_width=$half_width"))
    a = half_width
    h = sqrt(one(T) - a^2)
    symmetries = D2_symmetry() # sym_id 1=YAxisReflection, 2=XYAxisReflection, 3=XAxisReflection
    c = SVector{2,T}(zero(T), zero(T))
    pR = SVector{2,T}(one(T) + a, zero(T))
    pT = SVector{2,T}(zero(T), h)
    arc_angle = T(pi) - atan(h, a)
    arc = CircleSegment(one(T), arc_angle, zero(T), SVector{2,T}(a, zero(T)); bc=SpecularReflection(), domain_id=1, segment_id=1)
    ywall = LineSegment(pT, c; bc=SymmetryWall(1,2), domain_id=1, segment_id=2)
    xwall = LineSegment(c, pR; bc=SymmetryWall(3,2), domain_id=1, segment_id=3)
    fundamental_boundary = AbsCurve[arc, ywall, xwall]
    vertices = SVector{2,T}[pR, pT, c]
    fundamental_domain = SimpleDomain{T}(fundamental_boundary, vertices, 1)
    return CutCircleBilliard{T}(fundamental_domain, symmetries)
end
