
################################################################################ RECTANGLE WITHIN RECTANGLE BILLIARD
################################################################################

# Centered rectangle-within-rectangle billiard.
#
# Both physical boundary components:
#   - start at their positive x-axis midpoint,
#   - are parametrized counterclockwise,
#   - have identical normalized periodic parametrizations.
#
# This ensures that the D2 symmetry orbit map acts identically on both
# connected components.
#
# The outer rectangle has orientation=1 and the inner rectangle orientation=-1,
# so the latter is interpreted as a hole without reversing its parametrization.
#
# domain_id=1 -> outer rectangle
# domain_id=2 -> inner rectangle

struct RectangleWithinRectangleBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::MultiplyConnectedDomain{T}
    symmetries::SymmetryRegistry
end

function RectangleWithinRectangleBilliard(a_outer::T, b_outer::T, a_inner::T, b_inner::T;
    center=SVector{2,T}(zero(T), zero(T))) where T<:Real

    a_outer > zero(T) || throw(ArgumentError("a_outer must be positive; received $a_outer"))
    b_outer > zero(T) || throw(ArgumentError("b_outer must be positive; received $b_outer"))
    a_inner > zero(T) || throw(ArgumentError("a_inner must be positive; received $a_inner"))
    b_inner > zero(T) || throw(ArgumentError("b_inner must be positive; received $b_inner"))
    a_inner < a_outer || throw(ArgumentError("a_inner must be smaller than a_outer"))
    b_inner < b_outer || throw(ArgumentError("b_inner must be smaller than b_outer"))

    c = SVector{2,T}(center)
    iszero(c[1]) && iszero(c[2]) || throw(ArgumentError("D2 symmetry requires center == (0,0); received center=$c"))
    bc = SpecularReflection()

    ############################################################################
    # OUTER RECTANGLE — CCW, ORIENTATION +1
    ############################################################################

    o0 = c + SVector{2,T}( a_outer, zero(T))
    o1 = c + SVector{2,T}( a_outer,  b_outer)
    o2 = c + SVector{2,T}(-a_outer,  b_outer)
    o3 = c + SVector{2,T}(-a_outer, -b_outer)
    o4 = c + SVector{2,T}( a_outer, -b_outer)

    outer_right_upper = LineSegment(o0, o1; bc=bc, orientation=1, domain_id=1, segment_id=1)
    outer_top = LineSegment(o1, o2; bc=bc, orientation=1, domain_id=1, segment_id=2)
    outer_left = LineSegment(o2, o3; bc=bc, orientation=1, domain_id=1, segment_id=3)
    outer_bottom = LineSegment(o3, o4; bc=bc, orientation=1, domain_id=1, segment_id=4)
    outer_right_lower = LineSegment(o4, o0; bc=bc, orientation=1, domain_id=1, segment_id=5)

    outer = AbsCurve[
        outer_right_upper, outer_top, outer_left,
        outer_bottom, outer_right_lower
    ]

    ############################################################################
    # INNER RECTANGLE — SAME CCW PARAMETRIZATION, ORIENTATION -1
    ############################################################################

    i0 = c + SVector{2,T}( a_inner, zero(T))
    i1 = c + SVector{2,T}( a_inner,  b_inner)
    i2 = c + SVector{2,T}(-a_inner,  b_inner)
    i3 = c + SVector{2,T}(-a_inner, -b_inner)
    i4 = c + SVector{2,T}( a_inner, -b_inner)

    inner_right_upper = LineSegment(i0, i1; bc=bc, orientation=-1, domain_id=2, segment_id=1)
    inner_top = LineSegment(i1, i2; bc=bc, orientation=-1, domain_id=2, segment_id=2)
    inner_left = LineSegment(i2, i3; bc=bc, orientation=-1, domain_id=2, segment_id=3)
    inner_bottom = LineSegment(i3, i4; bc=bc, orientation=-1, domain_id=2, segment_id=4)
    inner_right_lower = LineSegment(i4, i0; bc=bc, orientation=-1, domain_id=2, segment_id=5)

    inner = AbsCurve[
        inner_right_upper, inner_top, inner_left,
        inner_bottom, inner_right_lower
    ]

    ############################################################################
    # MULTIPLY-CONNECTED DOMAIN + D2 SYMMETRY
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

    symmetries = D2_symmetry()
    return RectangleWithinRectangleBilliard{T}(fundamental_domain, symmetries)
end