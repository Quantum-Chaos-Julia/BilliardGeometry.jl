"""
    RectangleWithinRectangleBilliard{T}<:BilliardGeometry.AbsBilliard

Multiply connected billiard formed by a rectangular outer wall and a smaller
centered rectangular interior obstacle.

Each rectangle is one connected physical boundary component. For symmetry-
compatible periodic discretization, each component begins at its positive
x-axis midpoint and is traversed counterclockwise. The right side is therefore
split into two smooth curve pieces; this introduces no additional geometric
corner.

The physical boundary is

    full_boundary = [outer_boundary,inner_boundary],

with the outer boundary first and the rectangular hole second. The fundamental
domain is the first-quadrant L-shaped region bounded by the physical outer and
inner walls and the two coordinate-axis symmetry walls.
"""
struct RectangleWithinRectangleBilliard{T}<:BilliardGeometry.AbsBilliard
    fundamental_domain::BilliardGeometry.SimpleDomain{T}
    full_boundary::Vector{Vector{BilliardGeometry.AbsCurve}}
    symmetries::Vector{BilliardGeometry.AbsSymmetry}
end

"""
    RectangleWithinRectangleBilliard(
        a_outer::T,
        b_outer::T,
        a_inner::T,
        b_inner::T;
        center=SVector{2,T}(zero(T),zero(T)),
    ) where {T<:Real} → RectangleWithinRectangleBilliard{T}

Construct a centered rectangular billiard containing a smaller centered
rectangular obstacle.

The dimensions are

    outer: 2a_outer × 2b_outer,
    inner: 2a_inner × 2b_inner.

Each rectangular component starts at its positive x-axis midpoint and follows
the canonical counterclockwise periodic orientation required by the exact
integer symmetry maps.

## Arguments
* `a_outer::T`: Outer half-width.
* `b_outer::T`: Outer half-height.
* `a_inner::T`: Inner half-width.
* `b_inner::T`: Inner half-height.

## Keyword Arguments
* `center::SVector{2,T}`: Common center.
"""
function RectangleWithinRectangleBilliard(a_outer::T,b_outer::T,a_inner::T,b_inner::T;center=SVector{2,T}(zero(T),zero(T))) where {T<:Real}
    c=SVector{2,T}(center)
    o0=c+SVector{2,T}(a_outer,zero(T))
    o1=c+SVector{2,T}(a_outer,b_outer)
    o2=c+SVector{2,T}(-a_outer,b_outer)
    o3=c+SVector{2,T}(-a_outer,-b_outer)
    o4=c+SVector{2,T}(a_outer,-b_outer)
    i0=c+SVector{2,T}(a_inner,zero(T))
    i1=c+SVector{2,T}(a_inner,b_inner)
    i2=c+SVector{2,T}(-a_inner,b_inner)
    i3=c+SVector{2,T}(-a_inner,-b_inner)
    i4=c+SVector{2,T}(a_inner,-b_inner)
    oy=c+SVector{2,T}(zero(T),b_outer)
    iy=c+SVector{2,T}(zero(T),b_inner)
    bc=BilliardGeometry.SpecularReflection()
    outer_right_upper=BilliardGeometry.LineSegment(o0,o1;bc=bc,domain_id=1,segment_id=1)
    outer_top=BilliardGeometry.LineSegment(o1,o2;bc=bc,domain_id=1,segment_id=2)
    outer_left=BilliardGeometry.LineSegment(o2,o3;bc=bc,domain_id=1,segment_id=3)
    outer_bottom=BilliardGeometry.LineSegment(o3,o4;bc=bc,domain_id=1,segment_id=4)
    outer_right_lower=BilliardGeometry.LineSegment(o4,o0;bc=bc,domain_id=1,segment_id=5)
    inner_right_upper=BilliardGeometry.LineSegment(i0,i1;bc=bc,domain_id=1,segment_id=6)
    inner_top=BilliardGeometry.LineSegment(i1,i2;bc=bc,domain_id=1,segment_id=7)
    inner_left=BilliardGeometry.LineSegment(i2,i3;bc=bc,domain_id=1,segment_id=8)
    inner_bottom=BilliardGeometry.LineSegment(i3,i4;bc=bc,domain_id=1,segment_id=9)
    inner_right_lower=BilliardGeometry.LineSegment(i4,i0;bc=bc,domain_id=1,segment_id=10)
    full_boundary=[
        BilliardGeometry.AbsCurve[outer_right_upper,outer_top,outer_left,outer_bottom,outer_right_lower],
        BilliardGeometry.AbsCurve[inner_right_upper,inner_top,inner_left,inner_bottom,inner_right_lower]
    ]
    outer_right=BilliardGeometry.LineSegment(o0,o1;bc=bc,domain_id=1,segment_id=1)
    outer_top_q=BilliardGeometry.LineSegment(o1,oy;bc=bc,domain_id=1,segment_id=2)
    ywall=BilliardGeometry.LineSegment(oy,iy;bc=BilliardGeometry.ReflectionSymmetry(BilliardGeometry.YAxisReflection(),4),domain_id=1,segment_id=3)
    inner_top_q=BilliardGeometry.LineSegment(iy,i1;bc=bc,domain_id=1,segment_id=4) # clockwise obstacle orientation in the fundamental domain
    inner_right_q=BilliardGeometry.LineSegment(i1,i0;bc=bc,domain_id=1,segment_id=5) # clockwise obstacle orientation in the fundamental domain
    xwall=BilliardGeometry.LineSegment(i0,o0;bc=BilliardGeometry.ReflectionSymmetry(BilliardGeometry.XAxisReflection(),4),domain_id=1,segment_id=6)
    fundamental_boundary=BilliardGeometry.AbsCurve[outer_right,outer_top_q,ywall,inner_top_q,inner_right_q,xwall]
    vertices=SVector{2,T}[o0,o1,oy,iy,i1,i0]
    fundamental_domain=BilliardGeometry.SimpleDomain{T}(fundamental_boundary,vertices,1)
    symmetries=BilliardGeometry.AbsSymmetry[BilliardGeometry.D2_symmetry...]
    return RectangleWithinRectangleBilliard{T}(fundamental_domain,full_boundary,symmetries)
end

"""
    DiagonalRectangleWithinRectangleBilliard(a_outer::T,a_inner::T,shift::T) where {T<:Real}

Construct the diagonal fundamental domain of a square-within-square billiard
whose inner square is shifted by `(shift,shift)`.

The geometry is translated so that the re-entrant corner is exactly at the
origin, as required by the corner-adapted Fourier-Bessel (CAFB) basis.

Before translation the vertices are

    A=(-a_outer,-a_outer)
    B=( a_outer,-a_outer)
    C=( a_outer, a_outer)
    D=( shift+a_inner, shift+a_inner)
    E=( shift+a_inner, shift-a_inner)
    F=( shift-a_inner, shift-a_inner)

with counterclockwise ordering

    A -> B -> C -> D -> E -> F -> A.

The whole geometry is translated by `-E`, so that

    E=(0,0)
    D=(0,2a_inner)
    F=(-2a_inner,0).

The two segments meeting at the re-entrant corner,

    D -> E
    E -> F,

form the `3π/2` CAFB corner. Since the corner-adapted basis satisfies the
boundary condition analytically on these two rays, these segments are marked
with `QuantumSolverIgnore()`.

The diagonal segments

    C -> D
    F -> A

are boundaries of the diagonal fundamental domain and are therefore retained
by the quantum solver.

The missing inner-square sector is the upper-left quadrant relative to the
CAFB origin, while the billiard occupies the complementary `3π/2` sector.
"""
struct DiagonalRectangleWithinRectangleBilliard{T}<:BilliardGeometry.AbsBilliard
    full_boundary::Vector{BilliardGeometry.AbsCurve}
    fundamental_domain::BilliardGeometry.SimpleDomain{T}
    symmetries::Vector{BilliardGeometry.AbsSymmetry}
end

function DiagonalRectangleWithinRectangleBilliard(a_outer::T,a_inner::T,shift::T) where {T<:Real}
    a_outer>zero(T)||throw(ArgumentError("a_outer must be positive"))
    a_inner>zero(T)||throw(ArgumentError("a_inner must be positive"))
    a_inner<a_outer||throw(ArgumentError("a_inner must be smaller than a_outer"))
    lo=shift-a_inner
    hi=shift+a_inner
    -a_outer<lo<hi<a_outer||throw(ArgumentError("shifted inner square must lie strictly inside the outer square"))
    A0=SVector{2,T}(-a_outer,-a_outer)
    B0=SVector{2,T}(a_outer,-a_outer)
    C0=SVector{2,T}(a_outer,a_outer)
    D0=SVector{2,T}(hi,hi)
    E0=SVector{2,T}(hi,lo)
    F0=SVector{2,T}(lo,lo)
    A=A0-E0
    B=B0-E0
    C=C0-E0
    D=D0-E0
    E=SVector{2,T}(zero(T),zero(T))
    F=F0-E0
    corner_vertical=BilliardGeometry.LineSegment(D,E;bc=BilliardGeometry.QuantumSolverIgnore(),domain_id=1,segment_id=1)
    corner_horizontal=BilliardGeometry.LineSegment(E,F;bc=BilliardGeometry.QuantumSolverIgnore(),domain_id=1,segment_id=2)
    diagonal_lower=BilliardGeometry.LineSegment(F,A;bc=BilliardGeometry.SpecularReflection(),domain_id=1,segment_id=3)
    outer_bottom=BilliardGeometry.LineSegment(A,B;bc=BilliardGeometry.SpecularReflection(),domain_id=1,segment_id=4)
    outer_right=BilliardGeometry.LineSegment(B,C;bc=BilliardGeometry.SpecularReflection(),domain_id=1,segment_id=5)
    diagonal_upper=BilliardGeometry.LineSegment(C,D;bc=BilliardGeometry.SpecularReflection(),domain_id=1,segment_id=6)
    boundary=BilliardGeometry.AbsCurve[corner_vertical,corner_horizontal,diagonal_lower,outer_bottom,outer_right,diagonal_upper]
    vertices=SVector{2,T}[D,E,F,A,B,C]
    fundamental_domain=BilliardGeometry.SimpleDomain(boundary,vertices,1)
    symmetries=BilliardGeometry.AbsSymmetry[BilliardGeometry.DiagonalReflection()]
    return DiagonalRectangleWithinRectangleBilliard{T}(boundary,fundamental_domain,symmetries)
end