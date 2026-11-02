"""
    PolygonBilliard{T} <: AbsBilliard

Simply connected polygonal billiard with specular edges, built from vertices
ordered counterclockwise.

## Attributes
* `fundamental_domain::Union{SimpleDomain{T},CompositeDomain}`: The polygon, with edge `i` from vertex `i` to vertex `i+1` and `segment_id = i`; a [`CompositeDomain`](@ref) of convex pieces when built from pieces (needed for non-convex polygons).
* `symmetries::SymmetryRegistry`: Registered symmetries, empty by default.
"""
struct PolygonBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::Union{SimpleDomain{T},CompositeDomain}
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

"""
    PolygonBilliard(pieces; symmetries = register_symmetries(), atol = 1e-10) → billiard::PolygonBilliard

Construct a possibly non-convex [`PolygonBilliard`](@ref) as a
[`CompositeDomain`](@ref) of convex `pieces`, each a [`Polygon`](@ref) with
`id` equal to its index in `pieces`.

## Arguments
* `pieces`: Vector of vertex vectors (`SVector{2,T}`), one per convex piece, each ordered either way. Pieces must not overlap.

## Keyword arguments
* `symmetries::SymmetryRegistry` : Symmetries of the polygon. The pieces must already describe the fundamental domain.
* `atol` : Tolerance for matching the shared edges of neighbouring pieces.

## Description
An edge of a piece that coincides (as a full edge) with an edge of another piece
becomes a [`Transparent`](@ref) join; every other edge is a [`SpecularReflection`](@ref) wall.

## Throws
`ArgumentError` for a piece with fewer than three vertices, a degenerate or non convex piece.
"""
function PolygonBilliard(pieces::AbstractVector{<:AbstractVector{<:SVector{2,T}}}; symmetries::SymmetryRegistry=register_symmetries(), atol::Real=1e-10) where T<:Real
    ccw = map(pieces) do piece
        length(piece) >= 3 || throw(ArgumentError("every piece requires at least three vertices"))
        verts = SVector{2,T}.(piece)
        n = length(verts)
        cross(i) = let p = verts[i], q = verts[mod1(i+1,n)], r = verts[mod1(i+2,n)]
            (q[1]-p[1])*(r[2]-q[2]) - (q[2]-p[2])*(r[1]-q[1])
        end
        turns = [cross(i) for i in 1:n]
        all(>=(zero(T)), turns) || all(<=(zero(T)), turns) || throw(ArgumentError("pieces must be convex"))
        sum(abs, turns) > zero(T) || throw(ArgumentError("degenerate piece"))
        sum(turns) > zero(T) ? verts : reverse(verts)
    end
    samepoint(p, q) = norm(p - q) <= atol
    domains = AbsSimpleDomain[]
    for (i, verts) in enumerate(ccw)
        n = length(verts)
        bcs = AbsBoundaryCondition[]
        for j in 1:n
            p, q = verts[j], verts[mod1(j+1,n)]
            neighbour = 0
            for (k, other) in enumerate(ccw)
                k == i && continue
                m = length(other)
                if any(l -> samepoint(other[l], q) && samepoint(other[mod1(l+1,m)], p), 1:m)
                    neighbour = k
                    break
                end
            end
            push!(bcs, neighbour == 0 ? SpecularReflection() : Transparent(neighbour))
        end
        push!(domains, Polygon(verts, i; bcs=bcs))
    end
    return PolygonBilliard{T}(CompositeDomain(domains), symmetries)
end
