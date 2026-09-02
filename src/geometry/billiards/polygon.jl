"""
    PolygonBilliard{T}<:BilliardGeometry.AbsBilliard

Simply connected polygonal billiard defined by an ordered collection of
vertices.

The vertices are stored implicitly through the boundary segments and must be
ordered counterclockwise around the physical domain. Consecutive vertices are
joined by straight [`BilliardGeometry.LineSegment`](@ref) boundary pieces, with
the final vertex connected back to the first.

For vertices

    v₁,v₂,...,vₙ,

the physical boundary is

    v₁ -> v₂ -> ... -> vₙ -> v₁.

Because the billiard is simply connected and no symmetry reduction is applied
by this type, the full physical domain is also its fundamental domain.

## Attributes
* `full_boundary`: Ordered physical boundary segments of the polygon.
* `fundamental_domain`: Full polygonal domain.
* `symmetries`: Symmetries attached to the billiard.
"""
struct PolygonBilliard{T}<:BilliardGeometry.AbsBilliard
    full_boundary::Vector{BilliardGeometry.AbsCurve}
    fundamental_domain::BilliardGeometry.SimpleDomain{T}
    symmetries::Vector{BilliardGeometry.AbsSymmetry}
end

"""
    PolygonBilliard(vertices::AbstractVector{<:SVector{2,T}};symmetries=BilliardGeometry.AbsSymmetry[]) where {T<:Real}

Construct a simply connected polygonal billiard from counterclockwise ordered
vertices.

Each pair of consecutive vertices is connected by a specularly reflecting
straight boundary segment, with the last vertex connected back to the first.
The polygon itself is used as the fundamental domain.

The constructor requires at least three vertices and checks that the signed
polygon area is positive, corresponding to counterclockwise orientation.
Self-intersection of the polygon is not checked.

## Arguments
* `vertices`: Polygon vertices in counterclockwise boundary order.
* `symmetries`: Optional symmetries to attach to the billiard.

## Returns
* `billiard::PolygonBilliard{T}`: Polygonal billiard with the supplied vertices
  and symmetries.
"""
function PolygonBilliard(vertices::AbstractVector{<:SVector{2,T}};symmetries=BilliardGeometry.AbsSymmetry[]) where {T<:Real}
    n=length(vertices)
    n>=3||throw(ArgumentError("a polygon requires at least three vertices"))
    verts=SVector{2,T}.(vertices)
    area2=zero(T)
    @inbounds for i=1:n
        p=verts[i]
        q=verts[mod1(i+1,n)]
        area2+=p[1]*q[2]-q[1]*p[2]
    end
    area2>zero(T)||throw(ArgumentError("polygon vertices must be ordered counterclockwise"))
    bc=BilliardGeometry.SpecularReflection()
    boundary=Vector{BilliardGeometry.AbsCurve}(undef,n)
    @inbounds for i=1:n
        boundary[i]=BilliardGeometry.LineSegment(verts[i],verts[mod1(i+1,n)];bc=bc,domain_id=1,segment_id=i)
    end
    fundamental_domain=BilliardGeometry.SimpleDomain{T}(boundary,verts,1)
    syms=BilliardGeometry.AbsSymmetry[symmetries...]
    return PolygonBilliard{T}(boundary,fundamental_domain,syms)
end

"""
    PentagonBilliard{T}<:BilliardGeometry.AbsBilliard

Regular pentagonal billiard with fivefold rotational symmetry.

The pentagon is centered at `center`, has circumradius `radius`, and is rotated
by `rotation_angle`. Its five physical vertices are ordered counterclockwise and
joined by specularly reflecting straight boundary segments.

The billiard has cyclic `C₅` rotational symmetry. Its fundamental domain is the
triangular sector

    O -> V₁ -> V₂ -> O,

where `O` is the center of the pentagon and `V₁,V₂` are two adjacent vertices.
The edge `V₁ -> V₂` is a physical boundary segment, while `O -> V₁` and
`V₂ -> O` are artificial rotational cuts and are therefore marked with
[`BilliardGeometry.QuantumSolverIgnore`](@ref).

The full physical boundary remains the complete five-sided polygon.

## Attributes
* `full_boundary`: Five ordered physical boundary segments of the pentagon.
* `fundamental_domain`: One-fifth triangular rotational fundamental domain.
* `symmetries`: Four nontrivial images of the fivefold rotational symmetry.
"""
struct PentagonBilliard{T}<:BilliardGeometry.AbsBilliard
    full_boundary::Vector{BilliardGeometry.AbsCurve}
    fundamental_domain::BilliardGeometry.SimpleDomain{T}
    symmetries::Vector{BilliardGeometry.AbsSymmetry}
end

"""
    PentagonBilliard(radius::T;center=SVector{2,T}(zero(T),zero(T)),rotation_angle=zero(T),sector::Int=0) where {T<:Real}

Construct a regular pentagonal billiard with fivefold rotational symmetry.

The five vertices lie on a circle of radius `radius` centered at `center`.
`rotation_angle` gives the polar angle of the first vertex, and successive
vertices are separated by `2π/5`.

The complete physical boundary is

    V₁ -> V₂ -> V₃ -> V₄ -> V₅ -> V₁.

The rotational fundamental domain is the triangle

    O -> V₁ -> V₂ -> O,

where `O=center`. The two radial edges are rotational identification cuts and
are marked with `QuantumSolverIgnore()`, while `V₁ -> V₂` is a physical
specular boundary.

The rotational symmetry images are constructed as

    Cn_symmetry(5,sector).

## Arguments
* `radius::T`: Circumradius of the regular pentagon.
* `center::SVector{2,T}`: Center of the pentagon.
* `rotation_angle`: Polar angle of the first vertex.
* `sector::Int`: Irreducible-representation sector of the `C₅` symmetry,
  interpreted modulo `5`.

## Returns
* `billiard::PentagonBilliard{T}`: Regular pentagonal billiard with a
  one-fifth triangular fundamental domain and `C₅` rotational symmetry.
"""
function PentagonBilliard(radius::T;center=SVector{2,T}(zero(T),zero(T)),rotation_angle=zero(T),sector::Int=0) where {T<:Real}
    radius>zero(T)||throw(ArgumentError("radius must be positive"))
    c=SVector{2,T}(center)
    θ0=T(rotation_angle)
    vertices=Vector{SVector{2,T}}(undef,5)
    @inbounds for i=1:5
        θ=θ0+T(2pi)*(i-1)/T(5)
        s,cθ=sincos(θ)
        vertices[i]=c+radius*SVector{2,T}(cθ,s)
    end
    bc=BilliardGeometry.SpecularReflection()
    full_boundary=Vector{BilliardGeometry.AbsCurve}(undef,5)
    @inbounds for i=1:5
        full_boundary[i]=BilliardGeometry.LineSegment(vertices[i],vertices[mod1(i+1,5)];bc=bc,domain_id=1,segment_id=i)
    end
    V1=vertices[1]
    V2=vertices[2]
    cut1=BilliardGeometry.LineSegment(c,V1;bc=BilliardGeometry.QuantumSolverIgnore(),domain_id=1,segment_id=1)
    physical=BilliardGeometry.LineSegment(V1,V2;bc=bc,domain_id=1,segment_id=2)
    cut2=BilliardGeometry.LineSegment(V2,c;bc=BilliardGeometry.QuantumSolverIgnore(),domain_id=1,segment_id=3)
    fundamental_boundary=BilliardGeometry.AbsCurve[cut1,physical,cut2]
    fundamental_vertices=SVector{2,T}[c,V1,V2]
    fundamental_domain=BilliardGeometry.SimpleDomain{T}(fundamental_boundary,fundamental_vertices,1)
    symmetries=BilliardGeometry.AbsSymmetry[BilliardGeometry.Cn_symmetry(5,sector)...]
    return PentagonBilliard{T}(full_boundary,fundamental_domain,symmetries)
end