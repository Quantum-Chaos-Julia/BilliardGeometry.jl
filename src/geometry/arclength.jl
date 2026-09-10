"""
    _arc_length_integrand(crv::C,t::T) where {C<:AbsCurve,T<:Real}

Evaluate the arclength integrand `|dr/dt|` of a parametrized curve.

## Arguments
* `crv::C`: Boundary curve.
* `t::T`: Curve parameter.

## Returns
* Speed `|dr/dt|` at `t`.
"""
function _arc_length_integrand(crv::C,t::T) where {C<:AbsCurve,T<:Real}
    f(t)=curve(crv,t)
    return norm(ForwardDiff.derivative(f,t))
end

"""
    arc_length(crv::C,t1::T;rtol=sqrt(eps(T)),atol=zero(T)) where {C<:AbsCurve,T<:Real}

Compute the arclength of `crv` from parameter `t=0` to `t=t1`.

The arclength is evaluated by adaptive quadrature of `|dr/dt|`.

## Arguments
* `crv::C`: Boundary curve.
* `t1::T`: Final curve parameter.
* `rtol`: Relative quadrature tolerance.
* `atol`: Absolute quadrature tolerance.

## Returns
* Arclength from `0` to `t1`.
"""
function arc_length(crv::C,t1::T;rtol=sqrt(eps(T)),atol=zero(T)) where {C<:AbsCurve,T<:Real}
    f(t)=_arc_length_integrand(crv,t)
    return quadgk(f,zero(T),t1;rtol,atol)[1]
end

"""
    arc_length(crv::C,ts::AbstractArray;rtol=sqrt(eps(eltype(ts))),atol=zero(eltype(ts))) where {C<:AbsCurve}

Compute the arclength of `crv` from `t=0` to each parameter value in `ts`.

## Arguments
* `crv::C`: Boundary curve.
* `ts::AbstractArray`: Curve parameter values.
* `rtol`: Relative quadrature tolerance.
* `atol`: Absolute quadrature tolerance.

## Returns
* Array of arclength values corresponding to `ts`.
"""
function arc_length(crv::C,ts::AbstractArray;rtol=sqrt(eps(eltype(ts))),atol=zero(eltype(ts))) where {C<:AbsCurve}
    return [arc_length(crv,t;rtol,atol) for t in ts]
end

"""
    arc_length(crv::C,pt::SVector{2,T};rtol=sqrt(eps(T)),atol=zero(T)) where {C<:AbsCurve,T<:Real}

Compute the arclength of `crv` from its parameter origin to the point `pt`.

The point is first mapped back to its curve parameter using `invert_curve`, after
which the arclength is evaluated by adaptive quadrature of `|dr/dt|`.

## Arguments
* `crv::C`: Boundary curve.
* `pt::SVector{2,T}`: Point on the curve.
* `rtol`: Relative quadrature tolerance.
* `atol`: Absolute quadrature tolerance.

## Returns
* Arclength from the curve origin to `pt`.
"""
function arc_length(crv::C,pt::SVector{2,T};rtol=sqrt(eps(T)),atol=zero(T)) where {C<:AbsCurve,T<:Real}
    t1=invert_curve(crv,pt)
    f(t)=_arc_length_integrand(crv,t)
    return quadgk(f,zero(T),t1;rtol,atol)[1]
end

"""
    construct_arc_length_interpolation(crv::C;rtol=1e-10,n_samples=100,interp_method=CubicSpline) where {C<:AbsCurve}

Construct interpolation maps between curve parameter `t` and arclength `s`.

The arclength is sampled at `n_samples` uniformly spaced parameter values on
`[0,1]`, and interpolation objects are constructed for both `s(t)` and its
inverse `t(s)`.

## Arguments
* `crv::C`: Boundary curve.
* `rtol`: Relative tolerance used in the arclength quadrature.
* `n_samples`: Number of parameter samples.
* `interp_method`: Interpolation constructor used for both mappings.

## Returns
* `s_of_t`: Interpolation mapping parameter `t` to arclength `s`.
* `t_of_s`: Interpolation mapping arclength `s` to parameter `t`.
"""
function construct_arc_length_interpolation(crv::C;rtol=1e-10,n_samples=100,interp_method=CubicSpline) where {C<:AbsCurve}
    t_samples=collect(range(0.0,1.0,length=n_samples))
    s_samples=zeros(n_samples)
    integrand(t)=_arc_length_integrand(crv,t)
    for i in 2:n_samples
        s_samples[i],_=quadgk(integrand,0.0,t_samples[i],rtol=rtol)
    end
    s_of_t=interp_method(s_samples,t_samples)
    t_of_s=interp_method(t_samples,s_samples)
    return s_of_t,t_of_s
end