"""
    GeneralizedSinaiBilliard{T} <: AbsBilliard

Generalized Sinai billiard of Barnett (Asymptotic rate of quantum ergodicity in
chaotic Euclidean billiards): the full domain is symmetric under reflections in
both axes, and its first-quadrant fundamental domain is bounded by the two
coordinate axes and two concave circular arcs that meet at the corner `(1,1)`.

## Description
The first arc crosses the x-axis at right angles at `(a,0)`, so its center lies
on the x-axis at `(cₐ,0)` with `cₐ = (2 - a²)/(2(1 - a))`, and it passes through
`(1,1)`. The second arc crosses the y-axis at right angles at `(0,b)` with
center `(0,c_b)`, `c_b = (2 - b²)/(2(1 - b))`. Both arcs are concave: the
billiard lies outside the two circles. The fundamental domain is a
[`SimpleDomain`](@ref) with two [`SymmetryWall`](@ref) segments on the axes and
the `D₂` symmetry group of [`D2_symmetry`](@ref). For `a == b` the billiard
also has the diagonal reflection: the fundamental domain is then the eighth
`0 ≤ y ≤ x` (x-axis, arc, diagonal) and the group is [`D4_symmetry`](@ref).

## Attributes
* `fundamental_domain::SimpleDomain{T}`: The quadrant, with `domain_id = 1`.
* `symmetries::SymmetryRegistry`: Reflections `YAxisReflection`, `XYAxisReflection`, `XAxisReflection` with `sym_id` `1, 2, 3`.
"""
struct GeneralizedSinaiBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::SimpleDomain{T}
    symmetries::SymmetryRegistry
end

"""
    GeneralizedSinaiBilliard(a, b) → billiard::GeneralizedSinaiBilliard

Construct a [`GeneralizedSinaiBilliard`](@ref).

## Arguments
* `a`: Distance from the origin to the arc on the x-axis, `0 < a < 1`.
* `b`: Distance from the origin to the arc on the y-axis, `0 < b < 1`.

## Throws
`ArgumentError` if `a` or `b` is not in `(0, 1)`.
"""
function GeneralizedSinaiBilliard(a::T, b::T) where T<:AbstractFloat
    (zero(T) < a < one(T)) || throw(ArgumentError("a must satisfy 0 < a < 1; received a=$a"))
    (zero(T) < b < one(T)) || throw(ArgumentError("b must satisfy 0 < b < 1; received b=$b"))
    cx = (2 - a^2)/(2*(1 - a))
    cy = (2 - b^2)/(2*(1 - b))
    o = SVector{2,T}(zero(T), zero(T))
    pa = SVector{2,T}(a, zero(T))
    pc = SVector{2,T}(one(T), one(T))
    pb = SVector{2,T}(zero(T), b)
    if a == b
        symmetries = D4_symmetry()
        xwall = LineSegment(o, pa; bc=SymmetryWall(7,8), domain_id=1, segment_id=3)
        arc_x = CircleSegment(cx - a, -atan(one(T), cx - one(T)), T(pi), SVector{2,T}(cx, zero(T)); bc=SpecularReflection(), orientation=-1, domain_id=1, segment_id=1)
        diagonal = LineSegment(pc, o; bc=SymmetryWall(1,2), domain_id=1, segment_id=2)
        fundamental_domain = SimpleDomain{T}(AbsCurve[arc_x, diagonal, xwall], SVector{2,T}[pa, pc, o], 1)
        return GeneralizedSinaiBilliard{T}(fundamental_domain, symmetries)
    end
    symmetries = D2_symmetry() # sym_id 1=YAxisReflection, 2=XYAxisReflection, 3=XAxisReflection
    xwall = LineSegment(o, pa; bc=SymmetryWall(3,2), domain_id=1, segment_id=4)
    arc_x = CircleSegment(cx - a, -atan(one(T), cx - one(T)), T(pi), SVector{2,T}(cx, zero(T)); bc=SpecularReflection(), orientation=-1, domain_id=1, segment_id=1)
    arc_y = CircleSegment(cy - b, -atan(one(T), cy - one(T)), -atan(cy - one(T), one(T)), SVector{2,T}(zero(T), cy); bc=SpecularReflection(), orientation=-1, domain_id=1, segment_id=2)
    ywall = LineSegment(pb, o; bc=SymmetryWall(1,2), domain_id=1, segment_id=3)
    fundamental_domain = SimpleDomain{T}(AbsCurve[arc_x, arc_y, ywall, xwall], SVector{2,T}[pa, pc, pb, o], 1)
    return GeneralizedSinaiBilliard{T}(fundamental_domain, symmetries)
end

GeneralizedSinaiBilliard(a::Real, b::Real) = GeneralizedSinaiBilliard(promote(float(a), float(b))...)
