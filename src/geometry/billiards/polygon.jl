#simply connected polygonal billiard from ordered (counterclockwise) vertices
struct PolygonBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::SimpleDomain{T}
    symmetries::Vector{AbsSymmetry}
end

function PolygonBilliard(vertices::AbstractVector{<:SVector{2,T}}; symmetries=AbsSymmetry[]) where T<:Real
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
    syms = AbsSymmetry[symmetries...]
    return PolygonBilliard{T}(fundamental_domain, syms)
end
