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
Construct the diagonal fundamental domain of a square-within-square billiard whose inner square is shifted by `(shift,shift)`.
The returned simply connected fundamental domain is the half `y<=x`, with vertices
    A=(-a_outer,-a_outer)
    B=( a_outer,-a_outer)
    C=( a_outer, a_outer)
    D=( shift+a_inner, shift+a_inner)
    E=( shift+a_inner, shift-a_inner)
    F=( shift-a_inner, shift-a_inner)
ordered counterclockwise as
    A -> B -> C -> D -> E -> F -> A.
The diagonal segments `C -> D` and `F -> A` are artificial symmetry cuts and are marked with `QuantumSolverIgnore()`, since the corresponding condition is satisfied analytically by the symmetry-adapted basis.
The re-entrant corner `E` has opening angle `3π/2` and is the natural origin for a corner-adapted Fourier-Bessel basis.
"""
struct DiagonalRectangleWithinRectangleBilliard{T}<:BilliardGeometry.AbsBilliard
    full_boundary::Vector{BilliardGeometry.AbsCurve}
    fundamental_domain::BilliardGeometry.SimpleDomain{T}
    symmetries::Vector{BilliardGeometry.AbsSymmetry}
    cafb_corner::SVector{2,T}
    cafb_corner_angle::T
    cafb_rotation_angle::T
end
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
    D=SVector{2,T}(shift+a_inner,shift+a_inner)
    E=SVector{2,T}(shift+a_inner,shift-a_inner)
    F=SVector{2,T}(shift-a_inner,shift-a_inner)
    outer_bottom=BilliardGeometry.LineSegment(A,B;bc=BilliardGeometry.SpecularReflection(),domain_id=1,segment_id=1)
    outer_right=BilliardGeometry.LineSegment(B,C;bc=BilliardGeometry.SpecularReflection(),domain_id=1,segment_id=2)
    diagonal_upper=BilliardGeometry.LineSegment(C,D;bc=BilliardGeometry.QuantumSolverIgnore(),domain_id=1,segment_id=3)
    inner_right=BilliardGeometry.LineSegment(D,E;bc=BilliardGeometry.SpecularReflection(),domain_id=1,segment_id=4)
    inner_bottom=BilliardGeometry.LineSegment(E,F;bc=BilliardGeometry.SpecularReflection(),domain_id=1,segment_id=5)
    diagonal_lower=BilliardGeometry.LineSegment(F,A;bc=BilliardGeometry.QuantumSolverIgnore(),domain_id=1,segment_id=6)
    boundary=BilliardGeometry.AbsCurve[outer_bottom,outer_right,diagonal_upper,inner_right,inner_bottom,diagonal_lower]
    vertices=SVector{2,T}[A,B,C,D,E,F]
    fundamental_domain=BilliardGeometry.SimpleDomain(boundary,vertices,1)
    symmetries=BilliardGeometry.AbsSymmetry[BilliardGeometry.DiagonalReflection()]
    cafb_corner=E
    cafb_corner_angle=T(3pi/2)
    cafb_rotation_angle=T(pi)
    return DiagonalRectangleWithinRectangleBilliard{T}(boundary,fundamental_domain,symmetries,cafb_corner,cafb_corner_angle,cafb_rotation_angle)
end