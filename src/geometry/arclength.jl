function _arc_length_integrand(crv::C, t::T) where {C<:AbsCurve, T<:Real}
    return norm(tangent(crv, t))
end

"""
    arc_length(crv::AbsCurve, t::Real; rtol = sqrt(eps(T)), atol = 0) → s
    arc_length(crv::AbsCurve, ts::AbstractArray; rtol, atol) → ss::Vector
    arc_length(crv::AbsCurve, pt::SVector{2}; rtol, atol) → s

Arc length of `crv` from its start `t = 0` to the parameter `t`, to each
parameter in `ts`, or to the point `pt` on the curve.

## Description
The length is the integral of the speed,

```math
s(t) = \\int_0^t \\left\\| \\mathbf{r}'(\\tau) \\right\\| d\\tau,
```

evaluated with adaptive Gauss–Kronrod quadrature using [`tangent`](@ref). For a
point, the parameter is first found with [`invert_curve`](@ref). The point
methods are specialized for [`LineSegment`](@ref) (Euclidean distance from the
start) and [`CircleSegment`](@ref) (radius times the swept angle, measured along
the sense of traversal and wrapped so that the start maps to `0`), and then
take no tolerance keywords.

## Arguments
* `crv`: Curve.
* `t`, `ts`: Parameter value(s) in `[0,1]`.
* `pt`: Point on the curve.

## Keyword arguments
* `rtol = sqrt(eps(T))` : Relative tolerance of the quadrature.
* `atol = 0` : Absolute tolerance of the quadrature.

## Returns
* `s`, `ss`: Arc length(s).
"""
function arc_length(crv::C, t1::T; rtol=sqrt(eps(T)), atol=zero(T)) where {C<:AbsCurve, T<:Real}
    f(t) = _arc_length_integrand(crv,t)
    return quadgk(f, zero(T), t1; rtol, atol)[1]
end

function arc_length(crv::C, ts::AbstractArray; rtol=sqrt(eps(eltype(ts))), atol=zero(eltype(ts))) where {C<:AbsCurve}
    return [arc_length(crv,t; rtol, atol) for t in ts]
end

function arc_length(crv::C, pt::SVector{2,T}; rtol=sqrt(eps(T)), atol=zero(T)) where {C<:AbsCurve, T<:Real}
    t1 = invert_curve(crv, pt)
    f(t) = _arc_length_integrand(crv,t)
    return quadgk(f, zero(T), t1; rtol, atol)[1]
end

"""
    construct_arc_length_interpolation(crv::AbsCurve; rtol = 1e-10, n_samples = 100, interp_method = CubicSpline) → (s_of_t, t_of_s)

Build interpolants of the arc length as a function of the parameter and of its
inverse.

## Description
The arc length is integrated by quadrature at `n_samples` equally spaced
parameter values and interpolated with `interp_method` from
DataInterpolations.jl, once in each direction.

## Arguments
* `crv`: Curve.

## Keyword arguments
* `rtol = 1e-10` : Relative tolerance of the quadrature.
* `n_samples = 100` : Number of sample points in `t ∈ [0,1]`.
* `interp_method = CubicSpline` : Interpolation type, constructed as `interp_method(values, nodes)`.

## Returns
* `s_of_t`: Interpolant `t ↦ s`.
* `t_of_s`: Interpolant `s ↦ t`.
"""
function construct_arc_length_interpolation(crv::C; rtol=1e-10, n_samples=100, interp_method=CubicSpline) where {C<:AbsCurve}
    t_samples = collect(range(0.0, 1.0, length=n_samples))
    s_samples = zeros(n_samples)
    integrand(t) = _arc_length_integrand(crv, t)
    for i in 2:n_samples
        s_samples[i], _ = quadgk(integrand, 0.0, t_samples[i], rtol=rtol)
    end
    s_of_t = interp_method(s_samples, t_samples)
    t_of_s = interp_method(t_samples, s_samples)
    return s_of_t, t_of_s
end