"""
    LShapeBilliard{T} <: BilliardGeometry.AbsBilliard

MATLAB-logo L-shaped billiard with diagonal reflection symmetry.

The full physical domain is the square

    [-a,a]×[-a,a]

with the upper-right quadrant

    (0,a]×(0,a]

removed. Its counterclockwise boundary is

    A -> B -> C -> D -> E -> F -> A,

with

    A=(0,0)
    B=(a,0)
    C=(a,-a)
    D=(-a,-a)
    E=(-a,a)
    F=(0,a).

The reentrant corner `A=(0,0)` has opening angle `3π/2`.

The domain is invariant under diagonal reflection

    (x,y) -> (y,x).

The stored fundamental domain is one diagonal half of the L. For the odd
diagonal sector, parity `-1` corresponds to a Dirichlet condition on the
diagonal symmetry boundary.

## Attributes
* `fundamental_domain`: One diagonal half of the L-shaped domain.
* `full_boundary`: Complete six-segment physical L-shaped boundary.
* `symmetries`: Exact diagonal reflection symmetry.
* `cafb_corner`: Reentrant corner at the origin.
* `cafb_corner_angle`: Reentrant opening angle `3π/2`.
* `cafb_rotation_angle`: Local CAFB rotation angle at the reentrant corner.
"""
struct LShapeBilliard{T}<:BilliardGeometry.AbsBilliard
    fundamental_domain::BilliardGeometry.SimpleDomain{T}
    full_boundary::Vector{BilliardGeometry.AbsCurve}
    symmetries::Vector{BilliardGeometry.AbsSymmetry}
    cafb_corner::SVector{2,T}
    cafb_corner_angle::T
    cafb_rotation_angle::T
end

function LShapeBilliard(a::T) where {T<:Real}
    a>zero(T)||throw(ArgumentError("a must be positive"))
    z=zero(T)
    A=SVector{2,T}(z,z)
    B=SVector{2,T}(a,z)
    C=SVector{2,T}(a,-a)
    D=SVector{2,T}(-a,-a)
    E=SVector{2,T}(-a,a)
    F=SVector{2,T}(z,a)
    bc=BilliardGeometry.SpecularReflection()
    # Full physical L boundary, traversed counterclockwise:
    # A -> F -> E -> D -> C -> B -> A.
    c1=BilliardGeometry.LineSegment(A,F;bc=bc,domain_id=1,segment_id=1)
    c2=BilliardGeometry.LineSegment(F,E;bc=bc,domain_id=1,segment_id=2)
    c3=BilliardGeometry.LineSegment(E,D;bc=bc,domain_id=1,segment_id=3)
    c4=BilliardGeometry.LineSegment(D,C;bc=bc,domain_id=1,segment_id=4)
    c5=BilliardGeometry.LineSegment(C,B;bc=bc,domain_id=1,segment_id=5)
    c6=BilliardGeometry.LineSegment(B,A;bc=bc,domain_id=1,segment_id=6)
    full_boundary=BilliardGeometry.AbsCurve[c1,c2,c3,c4,c5,c6]
    # Diagonal fundamental domain y<=x, also traversed counterclockwise:
    # A -> D is the diagonal symmetry wall;
    # D -> C -> B -> A are physical boundary pieces.
    f1=BilliardGeometry.LineSegment(A,D;bc=BilliardGeometry.ReflectionSymmetry(BilliardGeometry.DiagonalReflection(),-1),domain_id=1,segment_id=1)
    f2=BilliardGeometry.LineSegment(D,C;bc=bc,domain_id=1,segment_id=2)
    f3=BilliardGeometry.LineSegment(C,B;bc=bc,domain_id=1,segment_id=3)
    f4=BilliardGeometry.LineSegment(B,A;bc=bc,domain_id=1,segment_id=4)
    fundamental_boundary=BilliardGeometry.AbsCurve[f1,f2,f3,f4]
    vertices=SVector{2,T}[A,D,C,B]
    fundamental_domain=BilliardGeometry.SimpleDomain(fundamental_boundary,vertices,1)
    symmetries=BilliardGeometry.AbsSymmetry[BilliardGeometry.DiagonalReflection(-1)]
    return LShapeBilliard{T}(fundamental_domain,full_boundary,symmetries,A,T(3pi/2),T(pi))
end