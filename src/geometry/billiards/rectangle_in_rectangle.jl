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
    DiagonalRectangleWithinRectangleBilliard{T} <: BilliardGeometry.AbsBilliard

Diagonal fundamental-domain representation of a square-within-square billiard
whose inner square is shifted diagonally by `(shift,shift)`.

The outer square is centered at the origin and the inner square is centered at
`(shift,shift)`. Reflection across the diagonal `y=x` remains an exact symmetry,
so the billiard may be reduced to one diagonal half.

The fundamental-domain boundary is the polygon

    A -> B -> C -> D -> E -> F -> A,

with vertices

    A=(-a_outer,-a_outer)
    B=( a_outer,-a_outer)
    C=( a_outer, a_outer)
    D=( shift+a_inner, shift+a_inner)
    E=( shift+a_inner, shift-a_inner)
    F=( shift-a_inner, shift-a_inner).

The physical boundary pieces are `A -> B`, `B -> C`, `D -> E`, and `E -> F`.
The segments `C -> D` and `F -> A` lie on the diagonal symmetry line.

## Attributes
* `full_boundary`: Ordered boundary curves of the diagonal fundamental domain.
* `fundamental_domain`: Polygonal diagonal fundamental domain.
* `symmetries`: Symmetries of the billiard; here diagonal reflection.
* `cafb_corner`: Reentrant corner used by a single-corner CAFB construction.
* `cafb_corner_angle`: Opening angle at `cafb_corner`.
* `cafb_rotation_angle`: Local CAFB rotation angle at `cafb_corner`.
"""
struct DiagonalRectangleWithinRectangleBilliard{T}<:BilliardGeometry.AbsBilliard
    full_boundary::Vector{BilliardGeometry.AbsCurve}
    fundamental_domain::BilliardGeometry.SimpleDomain{T}
    symmetries::Vector{BilliardGeometry.AbsSymmetry}
    cafb_corner::SVector{2,T}
    cafb_corner_angle::T
    cafb_rotation_angle::T
end

"""
    DiagonalRectangleWithinRectangleBilliard(a_outer::T,a_inner::T,shift::T) where {T<:Real}

Construct the diagonal fundamental domain of a square-within-square billiard
whose inner square is shifted diagonally by `(shift,shift)`.

The outer square has half-width `a_outer` and is centered at the origin. The
inner square has half-width `a_inner` and is centered at `(shift,shift)`.

The resulting fundamental-domain vertices are

    A=(-a_outer,-a_outer)
    B=( a_outer,-a_outer)
    C=( a_outer, a_outer)
    D=( shift+a_inner, shift+a_inner)
    E=( shift+a_inner, shift-a_inner)
    F=( shift-a_inner, shift-a_inner),

ordered counterclockwise as

    A -> B -> C -> D -> E -> F -> A.

The diagonal segments `C -> D` and `F -> A` are symmetry boundaries, while the
remaining segments are physical boundary pieces.

## Arguments
* `a_outer::T`: Half-width of the outer square.
* `a_inner::T`: Half-width of the inner square.
* `shift::T`: Diagonal displacement of the inner-square center, so that its
  center is `(shift,shift)`.

## Returns
* `billiard::DiagonalRectangleWithinRectangleBilliard{T}`: Diagonal
  fundamental-domain billiard with diagonal-reflection symmetry.
"""
function DiagonalRectangleWithinRectangleBilliard(a_outer::T,a_inner::T,shift::T) where {T<:Real}
    a_outer>zero(T)||throw(ArgumentError("a_outer must be positive"))
    a_inner>zero(T)||throw(ArgumentError("a_inner must be positive"))
    a_inner<a_outer||throw(ArgumentError("a_inner must be smaller than a_outer"))
    lo=shift-a_inner
    hi=shift+a_inner
    -a_outer<lo<hi<a_outer||throw(ArgumentError("shifted inner square must lie strictly inside the outer square"))
    A=SVector{2,T}(-a_outer,-a_outer)
    B=SVector{2,T}(a_outer,-a_outer)
    C=SVector{2,T}(a_outer,a_outer)
    D=SVector{2,T}(hi,hi)
    E=SVector{2,T}(hi,lo)
    F=SVector{2,T}(lo,lo)
    bc=BilliardGeometry.SpecularReflection()
    c1=BilliardGeometry.LineSegment(A,B;bc=bc,domain_id=1,segment_id=1)
    c2=BilliardGeometry.LineSegment(B,C;bc=bc,domain_id=1,segment_id=2)
    c3=BilliardGeometry.LineSegment(C,D;bc=bc,domain_id=1,segment_id=3)
    c4=BilliardGeometry.LineSegment(D,E;bc=bc,domain_id=1,segment_id=4)
    c5=BilliardGeometry.LineSegment(E,F;bc=bc,domain_id=1,segment_id=5)
    c6=BilliardGeometry.LineSegment(F,A;bc=bc,domain_id=1,segment_id=6)
    boundary=BilliardGeometry.AbsCurve[c1,c2,c3,c4,c5,c6]
    vertices=SVector{2,T}[A,B,C,D,E,F]
    fundamental_domain=BilliardGeometry.SimpleDomain(boundary,vertices,1)
    symmetries=BilliardGeometry.AbsSymmetry[BilliardGeometry.DiagonalReflection()]
    return DiagonalRectangleWithinRectangleBilliard{T}(boundary,fundamental_domain,symmetries,E,T(3pi/2),T(pi))
end

"""
    DiagonallyShiftedRectangleWithinRectangleBilliard{T} <: BilliardGeometry.AbsBilliard

Full multiply-connected square-within-square billiard whose inner square is
shifted diagonally relative to the outer square.

The outer square is centered at the origin. The inner square is centered at
`(shift,shift)`, so reflection across the diagonal `y=x` remains an exact
symmetry, while reflection across the individual coordinate axes is generally
broken for nonzero `shift`.

The physical boundary consists of two disconnected components,

    full_boundary[1] = outer square,
    full_boundary[2] = inner square.

Each component begins at the midpoint of its positive-x side and is traversed
counterclockwise. The right side is split into two curve pieces so that this
starting convention is preserved without introducing an additional geometric
corner.

The associated `fundamental_domain` is the same diagonal half-domain represented
by [`DiagonalRectangleWithinRectangleBilliard`](@ref).

## Attributes
* `fundamental_domain`: Diagonal fundamental domain bounded by physical and
  diagonal-symmetry segments.
* `full_boundary`: Two disconnected physical boundary components, with the
  outer square first and the shifted inner square second.
* `symmetries`: Symmetries of the geometry; here diagonal reflection.
"""
struct DiagonallyShiftedRectangleWithinRectangleBilliard{T}<:BilliardGeometry.AbsBilliard
    fundamental_domain::BilliardGeometry.SimpleDomain{T}
    full_boundary::Vector{Vector{BilliardGeometry.AbsCurve}}
    symmetries::Vector{BilliardGeometry.AbsSymmetry}
end

"""
    DiagonallyShiftedRectangleWithinRectangleBilliard(a_outer::T,a_inner::T,shift::T) where {T<:Real}

Construct a full multiply-connected square-within-square billiard whose inner
square is shifted diagonally by `(shift,shift)`.
The outer square has half-width `a_outer` and is centered at `(0,0)`. The inner
square has half-width `a_inner` and is centered at `(shift,shift)`.
Because the displacement lies along `y=x`, diagonal reflection remains an exact
symmetry. The full physical boundary contains separate outer and inner
components, while the stored fundamental domain is the diagonal half-domain

    A -> B -> C -> D -> E -> F -> A,

with

    A=(-a_outer,-a_outer)
    B=( a_outer,-a_outer)
    C=( a_outer, a_outer)
    D=( shift+a_inner, shift+a_inner)
    E=( shift+a_inner, shift-a_inner)
    F=( shift-a_inner, shift-a_inner).

For periodic boundary discretization, each full square component begins at the
midpoint of its positive-x side and is traversed counterclockwise.

## Arguments
* `a_outer::T`: Half-width of the outer square.
* `a_inner::T`: Half-width of the inner square.
* `shift::T`: Diagonal displacement of the inner-square center, so that its
  center is `(shift,shift)`.

## Returns
* `billiard::DiagonallyShiftedRectangleWithinRectangleBilliard{T}`: Full
  multiply-connected billiard with diagonal-reflection symmetry.
"""
function DiagonallyShiftedRectangleWithinRectangleBilliard(a_outer::T,a_inner::T,shift::T) where {T<:Real}
    a_outer>zero(T)||throw(ArgumentError("a_outer must be positive"))
    a_inner>zero(T)||throw(ArgumentError("a_inner must be positive"))
    a_inner<a_outer||throw(ArgumentError("a_inner must be smaller than a_outer"))
    lo=shift-a_inner
    hi=shift+a_inner
    -a_outer<lo<hi<a_outer||throw(ArgumentError("shifted inner square must lie strictly inside the outer square"))
    bc=BilliardGeometry.SpecularReflection()
    o0=SVector{2,T}(a_outer,zero(T))
    o1=SVector{2,T}(a_outer,a_outer)
    o2=SVector{2,T}(-a_outer,a_outer)
    o3=SVector{2,T}(-a_outer,-a_outer)
    o4=SVector{2,T}(a_outer,-a_outer)
    i0=SVector{2,T}(hi,shift)
    i1=SVector{2,T}(hi,hi)
    i2=SVector{2,T}(lo,hi)
    i3=SVector{2,T}(lo,lo)
    i4=SVector{2,T}(hi,lo)
    outer_right_upper=BilliardGeometry.LineSegment(o0,o1;bc=bc,domain_id=1,segment_id=1)
    outer_top=BilliardGeometry.LineSegment(o1,o2;bc=bc,domain_id=1,segment_id=2)
    outer_left=BilliardGeometry.LineSegment(o2,o3;bc=bc,domain_id=1,segment_id=3)
    outer_bottom=BilliardGeometry.LineSegment(o3,o4;bc=bc,domain_id=1,segment_id=4)
    outer_right_lower=BilliardGeometry.LineSegment(o4,o0;bc=bc,domain_id=1,segment_id=5)
    inner_right_upper=BilliardGeometry.LineSegment(i0,i1;bc=bc,domain_id=2,segment_id=1)
    inner_top=BilliardGeometry.LineSegment(i1,i2;bc=bc,domain_id=2,segment_id=2)
    inner_left=BilliardGeometry.LineSegment(i2,i3;bc=bc,domain_id=2,segment_id=3)
    inner_bottom=BilliardGeometry.LineSegment(i3,i4;bc=bc,domain_id=2,segment_id=4)
    inner_right_lower=BilliardGeometry.LineSegment(i4,i0;bc=bc,domain_id=2,segment_id=5)
    outer=BilliardGeometry.AbsCurve[outer_right_upper,outer_top,outer_left,outer_bottom,outer_right_lower]
    inner=BilliardGeometry.AbsCurve[inner_right_upper,inner_top,inner_left,inner_bottom,inner_right_lower]
    A=SVector{2,T}(-a_outer,-a_outer)
    B=SVector{2,T}(a_outer,-a_outer)
    C=SVector{2,T}(a_outer,a_outer)
    D=SVector{2,T}(hi,hi)
    E=SVector{2,T}(hi,lo)
    F=SVector{2,T}(lo,lo)
    f1=BilliardGeometry.LineSegment(A,B;bc=bc,domain_id=1,segment_id=11)
    f2=BilliardGeometry.LineSegment(B,C;bc=bc,domain_id=1,segment_id=12)
    f3=BilliardGeometry.LineSegment(C,D;bc=bc,domain_id=1,segment_id=13)
    f4=BilliardGeometry.LineSegment(D,E;bc=bc,domain_id=1,segment_id=14)
    f5=BilliardGeometry.LineSegment(E,F;bc=bc,domain_id=1,segment_id=15)
    f6=BilliardGeometry.LineSegment(F,A;bc=bc,domain_id=1,segment_id=16)
    fundamental_boundary=BilliardGeometry.AbsCurve[f1,f2,f3,f4,f5,f6]
    vertices=SVector{2,T}[A,B,C,D,E,F]
    fundamental_domain=BilliardGeometry.SimpleDomain(fundamental_boundary,vertices,1)
    symmetries=BilliardGeometry.AbsSymmetry[BilliardGeometry.DiagonalReflection()]
    return DiagonallyShiftedRectangleWithinRectangleBilliard{T}(fundamental_domain,[outer,inner],symmetries)
end