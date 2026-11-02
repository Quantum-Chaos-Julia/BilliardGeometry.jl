"""
    SquashBilliard{T} <: AbsBilliard

Bunimovich squash: a unit circle centered at `(-half_width, 0)` and a circle of
radius `ratio` centered at `(half_width, 0)`, joined by the two common external
tangent lines. The tangent lines make an angle `α = asin((ratio - 1)/(2 half_width))`
with the `x` axis, so the circles join the straight part smoothly. For
`ratio = 1` the billiard is the [`StadiumBilliard`](@ref). The fundamental
domain is the upper half, a [`CompositeDomain`](@ref) of a left
[`CircleWedge`](@ref) (`id = 1`), a trapezoidal [`Polygon`](@ref) (`id = 2`) and
a right [`CircleWedge`](@ref) (`id = 3`), with the `x` axis a
[`SymmetryWall`](@ref).

## Attributes
* `fundamental_domain::CompositeDomain`: The upper half.
* `symmetries::SymmetryRegistry`: The single [`XAxisReflection`](@ref).
"""
struct SquashBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::CompositeDomain
    symmetries::SymmetryRegistry
end

"""
    SquashBilliard(half_width, ratio = 1) → billiard::SquashBilliard
    SquashBilliard(fundamental_domain::CompositeDomain, symmetries::SymmetryRegistry) → billiard::SquashBilliard

Construct a [`SquashBilliard`](@ref), or wrap a given composite fundamental
domain and symmetry registry.

## Arguments
* `half_width`: Half the distance between the circle centers, `half_width > 0`.
* `ratio = 1`: Radius of the right circle relative to the unit left circle, `ratio ≥ 1` and `ratio - 1 < 2 half_width`.

## Throws
`ArgumentError` if the parameters violate these conditions.
"""
function SquashBilliard(half_width::T1, ratio::T2 = one(T1)) where {T1<:Real,T2<:Real}
    T = float(promote_type(T1, T2))
    L = T(half_width); r = T(ratio)
    L > zero(T) || throw(ArgumentError("half_width must be positive; received half_width=$half_width"))
    r >= one(T) || throw(ArgumentError("ratio must satisfy ratio >= 1; received ratio=$ratio"))
    (r - one(T) < 2L) || throw(ArgumentError("ratio - 1 must be smaller than 2*half_width; received half_width=$half_width, ratio=$ratio"))
    symmetries = register_symmetries(XAxisReflection())
    α = asin((r - one(T)) / (2L))
    sα, cα = sincos(α)
    left_dom = CircleWedge(one(T), T(pi)/2 - α, T(pi)/2 + α, -L, zero(T), 1;
        bcs = [SpecularReflection(), SymmetryWall(1,2), Transparent(2)])
    trapezoid = Polygon([[-L, zero(T)], [L, zero(T)], [L - r*sα, r*cα], [-L - sα, cα]], 2;
        bcs = [SymmetryWall(1,2), Transparent(3), SpecularReflection(), Transparent(1)])
    right_dom = CircleWedge(r, T(pi)/2 + α, zero(T), L, zero(T), 3;
        bcs = [SpecularReflection(), Transparent(2), SymmetryWall(1,2)])
    return SquashBilliard{T}(CompositeDomain([left_dom, trapezoid, right_dom]), symmetries)
end

SquashBilliard(fundamental_domain::CompositeDomain, symmetries::SymmetryRegistry) =
    SquashBilliard{_domain_eltype(fundamental_domain)}(fundamental_domain, symmetries)
