"""
    _boundary_components(boundary) → comps::Vector

Normalizes a boundary representation into a vector of connected components.
If `boundary` is already a vector of component vectors (each component being
a vector of curve segments), it is returned unchanged. Otherwise every
supplied curve is treated as its own single-segment component.
"""
function _boundary_components(boundary)
    isempty(boundary) && throw(ArgumentError("Boundary cannot be empty"))
    return boundary[1] isa AbstractVector ? boundary : [[crv] for crv in boundary]
end

"""
    component_lengths(comp::Vector) → (lens, cum, Ltot)

Computes the length of each curve segment in `comp`, the cumulative lengths
`cum` (with `cum[1] = 0`), and the total component length `Ltot = cum[end]`.
"""
function component_lengths(comp::Vector)
    lens = [crv.length for crv in comp]
    cum = Vector{eltype(lens)}(undef, length(lens)+1)
    cum[1] = zero(eltype(lens))
    @inbounds for j in eachindex(lens)
        cum[j+1] = cum[j] + lens[j]
    end
    return lens, cum, cum[end]
end

@inline function _unit_tangent_at_start(crv, ::Type{T}) where {T<:Real}
    v = SVector{2,T}(tangent(crv, zero(T)))
    return v/hypot(v[1], v[2])
end

@inline function _unit_tangent_at_end(crv, ::Type{T}) where {T<:Real}
    v = SVector{2,T}(tangent(crv, one(T)))
    return v/hypot(v[1], v[2])
end

"""
    _junction_angle(cleft, cright, ::Type{T}) where {T<:Real} → angle::T

Unsigned turning angle in `[0,π]` between the tangent leaving `cleft` and the
tangent entering `cright`.
"""
@inline function _junction_angle(cleft, cright, ::Type{T}) where {T<:Real}
    tL = _unit_tangent_at_end(cleft, T)
    tR = _unit_tangent_at_start(cright, T)
    cr = tL[1]*tR[2] - tL[2]*tR[1]
    dt = clamp(tL[1]*tR[1] + tL[2]*tR[2], -one(T), one(T))
    return atan(abs(cr), dt)
end

"""
    _is_true_corner(cleft, cright, ::Type{T}; angle_tol=T(1e-8)) where {T<:Real} → flag::Bool

`true` when the junction between `cleft` and `cright` has a tangent
discontinuity larger than `angle_tol`.
"""
@inline function _is_true_corner(cleft, cright, ::Type{T}; angle_tol=T(1e-8)) where {T<:Real}
    return _junction_angle(cleft, cright, T) > angle_tol
end

"""
    _component_corner_locations(::Type{T}, comp::Vector; angle_tol=T(1e-8)) where {T<:Real} → corners::Vector{T}

Locations of true corners of a composite boundary component in the global
periodic parameter `σ ∈ [0,2π)`. Smooth joins are ignored. If the periodic
seam between the final and first segments is a true corner, it is
represented by `σ = 0`.
"""
function _component_corner_locations(::Type{T}, comp::Vector; angle_tol=T(1e-8)) where {T<:Real}
    _, cum, Ltot = component_lengths(comp)
    corners = T[]
    m = length(comp)
    _is_true_corner(comp[end], comp[1], T; angle_tol=angle_tol) && push!(corners, zero(T))
    @inbounds for j in 1:m-1
        _is_true_corner(comp[j], comp[j+1], T; angle_tol=angle_tol) && push!(corners, T(2*pi)*cum[j+1]/Ltot)
    end
    return corners
end

"""
    print_component_junctions(comp::Vector; T=Float64, angle_tol=1e-8)

Prints diagnostic information (global periodic parameter, tangent
discontinuity angle, corner classification) for every join of a composite
boundary component.
"""
function print_component_junctions(comp::Vector; T=Float64, angle_tol=1e-8)
    _, cum, Ltot = component_lengths(comp)
    m = length(comp)
    println("junction diagnostics:")
    @inbounds for j in 1:m
        jr = j == m ? 1 : j+1
        σ = j == m ? zero(T) : T(2*pi)*cum[j+1]/Ltot
        a = _junction_angle(comp[j], comp[jr], T)
        flag = a > T(angle_tol) ? "TRUE CORNER" : "smooth join"
        println("  $j -> $jr : σ = $σ, angle = $a, $flag")
    end
end
