function _arc_length_integrand(crv::C, t::T) where {C<:AbsCurve, T<:Real}
    f(t) = curve(crv,t)
    return norm(ForwardDiff.derivative(f, t))
end

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