"""
    PolygonBilliard{T} <: AbsBilliard

Simply connected polygonal billiard with specular edges, built from vertices
ordered counterclockwise.

## Attributes
* `fundamental_domain::SimpleDomain{T}`: The polygon, with edge `i` from vertex `i` to vertex `i+1` and `segment_id = i`.
* `symmetries::SymmetryRegistry`: Registered symmetries, empty by default.
"""
struct PolygonBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::SimpleDomain{T}
    symmetries::SymmetryRegistry
end

"""
    PolygonBilliard(vertices; symmetries = register_symmetries()) → billiard::PolygonBilliard

Construct a [`PolygonBilliard`](@ref).

## Arguments
* `vertices`: Vector of `SVector{2,T}`, at least three, in counterclockwise order.

## Keyword arguments
* `symmetries::SymmetryRegistry` : Symmetries of the polygon. The vertices are not reduced: they must already describe the fundamental domain, and no `SymmetryWall` edges are inserted.

## Throws
`ArgumentError` for fewer than three vertices or a non counterclockwise ordering (nonpositive signed area).
"""
function PolygonBilliard(vertices::AbstractVector{<:SVector{2,T}}; symmetries::SymmetryRegistry=register_symmetries()) where T<:Real
    n = length(vertices)
    n >= 3 || throw(ArgumentError("a polygon requires at least three vertices"))
    verts = SVector{2,T}.(vertices)
    area2 = zero(T)
    @inbounds for i in 1:n
        p = verts[i]
        q = verts[mod1(i+1,n)]
        area2 += p[1]*q[2] - q[1]*p[2]
    end
    area2 > zero(T) || throw(ArgumentError("polygon vertices must be ordered counterclockwise"))
    bc = SpecularReflection()
    boundary = Vector{AbsCurve}(undef,n)
    @inbounds for i in 1:n
        boundary[i] = LineSegment(verts[i], verts[mod1(i+1,n)]; bc=bc, domain_id=1, segment_id=i)
    end
    fundamental_domain = SimpleDomain{T}(boundary, verts, 1)
    return PolygonBilliard{T}(fundamental_domain, symmetries)
end
