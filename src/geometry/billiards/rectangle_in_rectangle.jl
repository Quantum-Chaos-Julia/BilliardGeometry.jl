using BilliardGeometry
using StaticArrays

################################################################################
# RECTANGLE WITHIN RECTANGLE BILLIARD
################################################################################

"""
    RectangleWithinRectangleBilliard{T} <: AbsBilliard

Multiply-connected rectangular billiard consisting of a rectangular outer wall
and a smaller centered rectangular obstacle.

The dimensions are

    outer: 2a_outer × 2b_outer
    inner: 2a_inner × 2b_inner

The full physical geometry is used; no symmetry reduction is applied.

Both connected boundary components begin at their positive x-axis midpoint.
The outer boundary is traversed counterclockwise and the inner boundary
clockwise, as required for a hole.

The two physical boundary components use distinct `domain_id`s:

    domain_id = 1  → outer rectangle
    domain_id = 2  → inner rectangle

so that they can be treated independently by `CompositeBIMSolver`.
"""
struct RectangleWithinRectangleBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::MultiplyConnectedDomain{T}
    symmetries::SymmetryRegistry
end

"""
    RectangleWithinRectangleBilliard(
        a_outer::T,
        b_outer::T,
        a_inner::T,
        b_inner::T;
        center=SVector{2,T}(zero(T),zero(T))
    ) where T<:Real

Construct a rectangle-within-rectangle billiard.

The outer rectangle has half-width `a_outer` and half-height `b_outer`.
The inner rectangular obstacle has half-width `a_inner` and half-height
`b_inner`. Both rectangles have the same center.

The outer boundary is oriented counterclockwise and the inner boundary
clockwise. Each component begins at the midpoint of its positive-x side.
"""
function RectangleWithinRectangleBilliard(a_outer::T, b_outer::T, a_inner::T, b_inner::T;
    center=SVector{2,T}(zero(T), zero(T))) where T<:Real

    a_outer > zero(T) || throw(ArgumentError("a_outer must be positive; received $a_outer"))
    b_outer > zero(T) || throw(ArgumentError("b_outer must be positive; received $b_outer"))
    a_inner > zero(T) || throw(ArgumentError("a_inner must be positive; received $a_inner"))
    b_inner > zero(T) || throw(ArgumentError("b_inner must be positive; received $b_inner"))
    a_inner < a_outer || throw(ArgumentError("a_inner must be smaller than a_outer; received a_inner=$a_inner, a_outer=$a_outer"))
    b_inner < b_outer || throw(ArgumentError("b_inner must be smaller than b_outer; received b_inner=$b_inner, b_outer=$b_outer"))

    c = SVector{2,T}(center); bc = SpecularReflection()

    ############################################################################
    # OUTER RECTANGLE — COUNTERCLOCKWISE
    ############################################################################

    o0 = c + SVector{2,T}( a_outer, zero(T))
    o1 = c + SVector{2,T}( a_outer,  b_outer)
    o2 = c + SVector{2,T}(-a_outer,  b_outer)
    o3 = c + SVector{2,T}(-a_outer, -b_outer)
    o4 = c + SVector{2,T}( a_outer, -b_outer)

    outer_right_upper = LineSegment(o0, o1; bc=bc, domain_id=1, segment_id=1)
    outer_top = LineSegment(o1, o2; bc=bc, domain_id=1, segment_id=2)
    outer_left = LineSegment(o2, o3; bc=bc, domain_id=1, segment_id=3)
    outer_bottom = LineSegment(o3, o4; bc=bc, domain_id=1, segment_id=4)
    outer_right_lower = LineSegment(o4, o0; bc=bc, domain_id=1, segment_id=5)

    outer = AbsCurve[
        outer_right_upper,
        outer_top,
        outer_left,
        outer_bottom,
        outer_right_lower
    ]

    ############################################################################
    # INNER RECTANGLE — CLOCKWISE
    ############################################################################

    i0 = c + SVector{2,T}( a_inner, zero(T))
    i1 = c + SVector{2,T}( a_inner, -b_inner)
    i2 = c + SVector{2,T}(-a_inner, -b_inner)
    i3 = c + SVector{2,T}(-a_inner,  b_inner)
    i4 = c + SVector{2,T}( a_inner,  b_inner)

    inner_right_lower = LineSegment(i0, i1; bc=bc, domain_id=2, segment_id=1)
    inner_bottom = LineSegment(i1, i2; bc=bc, domain_id=2, segment_id=2)
    inner_left = LineSegment(i2, i3; bc=bc, domain_id=2, segment_id=3)
    inner_top = LineSegment(i3, i4; bc=bc, domain_id=2, segment_id=4)
    inner_right_upper = LineSegment(i4, i0; bc=bc, domain_id=2, segment_id=5)

    inner = AbsCurve[
        inner_right_lower,
        inner_bottom,
        inner_left,
        inner_top,
        inner_right_upper
    ]

    ############################################################################
    # MULTIPLY-CONNECTED DOMAIN
    ############################################################################

    vertices = SVector{2,T}[
        o0, o1, o2, o3, o4,
        i0, i1, i2, i3, i4
    ]

    fundamental_domain = MultiplyConnectedDomain(
        outer,
        Vector{AbsCurve}[inner],
        vertices,
        1
    )

    symmetries = register_symmetries()
    return RectangleWithinRectangleBilliard{T}(fundamental_domain, symmetries)
end
