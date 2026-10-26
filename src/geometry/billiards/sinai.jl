"""
    SinaiBilliard{T} <: AbsBilliard

Sinai billiard: the square `[-a,a] × [-a,a]` with a circular obstacle of radius
`R` at the origin. The fundamental domain is the first quadrant of the square
with the quarter disk removed, with the `D₂` symmetry group of
[`D2_symmetry`](@ref). The obstacle arc is a [`CircleSegment`](@ref) with
`orientation = -1` and `domain_id = 2`.

## Attributes
* `fundamental_domain::SimpleDomain{T}`: The quadrant with the obstacle cut out.
* `symmetries::SymmetryRegistry`: Reflections `YAxisReflection`, `XYAxisReflection`, `XAxisReflection` with `sym_id` `1, 2, 3`.
"""
struct SinaiBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::SimpleDomain{T}
    symmetries::SymmetryRegistry
end

"""
    SinaiBilliard(a, R; center = (0,0)) → billiard::SinaiBilliard

Construct a [`SinaiBilliard`](@ref).

## Arguments
* `a`: Half side length of the square.
* `R`: Radius of the obstacle, `0 < R < a`.

## Keyword arguments
* `center` : Center of the square and obstacle; must be the origin.

## Throws
`ArgumentError` for non-positive `a` or `R`, `R ≥ a`, or a nonzero `center`.
"""
function SinaiBilliard(a::T, R::T; center=SVector{2,T}(zero(T),zero(T))) where T<:Real
    a > zero(T) || throw(ArgumentError("a must be positive; received $a"))
    R > zero(T) || throw(ArgumentError("R must be positive; received $R"))
    R < a || throw(ArgumentError("R must be smaller than a; received a=$a, R=$R"))
    c = SVector{2,T}(center)
    (iszero(c[1]) && iszero(c[2])) || throw(ArgumentError("D2 symmetry requires center == (0,0); received center=$c"))
    z = zero(T); bc = SpecularReflection()
    q1 = SVector{2,T}(a,z); p3 = SVector{2,T}(a,a); q3 = SVector{2,T}(z,a)
    p0 = SVector{2,T}(R,z); p1 = SVector{2,T}(z,R)
    symmetries = D2_symmetry()
    right = LineSegment(q1,p3; bc=bc, domain_id=1, segment_id=1)
    top = LineSegment(p3,q3; bc=bc, domain_id=1, segment_id=2)
    ywall = LineSegment(q3,p1; bc=SymmetryWall(1,2), domain_id=1, segment_id=3)
    inner = CircleSegment(R,T(pi/2),z,c; bc=bc, orientation=-1, domain_id=2, segment_id=1)
    xwall = LineSegment(p0,q1; bc=SymmetryWall(3,2), domain_id=1, segment_id=4)
    fundamental_boundary = AbsCurve[right,top,ywall,inner,xwall]
    vertices = SVector{2,T}[q1,p3,q3,p1,p0]
    fundamental_domain = SimpleDomain{T}(fundamental_boundary,vertices,1)
    return SinaiBilliard{T}(fundamental_domain,symmetries)
end