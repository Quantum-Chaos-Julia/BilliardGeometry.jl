struct SquareWithinSquareBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::SimpleDomain{T}
    symmetries::SymmetryRegistry
end

function SquareWithinSquareBilliard(a::T=one(T); positions::AbstractVector{<:Real}=T[-a/2,a/2], sizes::AbstractVector{<:Real}=fill(a/4,2)) where {T<:Real}
    a > zero(T) || throw(ArgumentError("a must be positive; received $a"))
    length(positions) == length(sizes) || throw(ArgumentError("positions and sizes must have the same length"))
    pos = T.(positions); sz = T.(sizes); n = length(pos)
    all(>(zero(T)), sz) || throw(ArgumentError("All square side lengths must be positive"))
    p = sortperm(pos); pos = pos[p]; sz = sz[p]
    @inbounds for i in 1:n
        h = sz[i]/2
        -a < pos[i]-h && pos[i]+h < a || throw(ArgumentError("Square hole $i must lie strictly inside the outer square; received position=$(pos[i]), size=$(sz[i])"))
    end
    @inbounds for i in 1:n-1
        pos[i]+sz[i]/2 < pos[i+1]-sz[i+1]/2 || throw(ArgumentError("Square holes $i and $(i+1) overlap or touch"))
    end
    z = zero(T); bc = SpecularReflection()
    p0 = SVector{2,T}(-a,-a); p1 = SVector{2,T}(a,-a); p2 = SVector{2,T}(a,a)
    bottom = LineSegment(p0,p1; bc=bc, domain_id=1, segment_id=1)
    right = LineSegment(p1,p2; bc=bc, domain_id=1, segment_id=2)
    boundary = AbsCurve[bottom,right]
    current = p2; wall_id = 3
    @inbounds for i in n:-1:1
        d = pos[i]; h = sz[i]/2
        qhi = SVector{2,T}(d+h,d+h)
        qx = SVector{2,T}(d+h,d-h)
        qlo = SVector{2,T}(d-h,d-h)
        push!(boundary,LineSegment(current,qhi; bc=SymmetryWall(1,2), domain_id=1, segment_id=wall_id)); wall_id += 1
        push!(boundary,LineSegment(qlo,qx; bc=bc, orientation=-1, domain_id=i+1, segment_id=1))
        push!(boundary,LineSegment(qx,qhi; bc=bc, orientation=-1, domain_id=i+1, segment_id=2))
        current = qlo
    end
    push!(boundary,LineSegment(current,p0; bc=SymmetryWall(1,2), domain_id=1, segment_id=wall_id))
    vertices = SVector{2,T}[curve(c,z) for c in boundary]
    fundamental_domain = SimpleDomain{T}(boundary,vertices,1)
    symmetries = register_symmetries(DiagonalReflection())
    return SquareWithinSquareBilliard{T}(fundamental_domain,symmetries)
end