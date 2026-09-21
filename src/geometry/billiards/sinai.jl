struct SinaiBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::SimpleDomain{T}
    symmetries::SymmetryRegistry
end

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