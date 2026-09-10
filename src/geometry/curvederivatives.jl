"""
    tangent(crv::AbsCurve, ts::AbstractArray) → dr::Vector

Evaluate the first parameter derivative `dr/dt` of `crv` at every parameter
value in `ts`. This is a derivative with respect to the curve parameter `t`,
not a normalized unit tangent and not a derivative with respect to arclength.

Falls back to broadcasting the scalar-argument method; concrete curve types
implement `tangent(crv, t::Real)`.
"""
function tangent(crv::AbsCurve, ts::AbstractArray{<:Real})
    return [tangent(crv, t) for t in ts]
end

"""
    tangent_2(crv::AbsCurve, ts::AbstractArray) → ddr::Vector

Evaluate the second parameter derivative `d²r/dt²` of `crv` at every parameter
value in `ts`. Falls back to broadcasting the scalar-argument method.
"""
function tangent_2(crv::AbsCurve, ts::AbstractArray{<:Real})
    return [tangent_2(crv, t) for t in ts]
end

# Line segments: r(t) = pt0 + t*(pt1-pt0), affine in t.
@inline function tangent(line::L, t::Real) where {L<:AbsLine}
    return line.pt1 - line.pt0
end

@inline function tangent_2(line::L, t::Real) where {L<:AbsLine}
    return zero(line.pt0)
end

# Circular segments: phi(t) = shift_angle + arc_angle*t, r(t) = center + radius*(cos(phi),sin(phi)).
@inline function tangent(circle::L, t::Real) where {L<:CircleSegment}
    phi = circle.arc_angle*t + circle.shift_angle
    Ra = circle.radius*circle.arc_angle
    return SVector(-Ra*sin(phi), Ra*cos(phi))
end

@inline function tangent_2(circle::L, t::Real) where {L<:CircleSegment}
    phi = circle.arc_angle*t + circle.shift_angle
    Ra2 = circle.radius*circle.arc_angle^2
    return SVector(-Ra2*cos(phi), -Ra2*sin(phi))
end

# FourierCoeffPolarSegment: r(phi) = R + Σ aₙcos(nφ) + Σ bₙsin(nφ), coef = [b₁,a₁,b₂,a₂,...]
# (same coefficient convention as `polar_radius`). Both angular derivatives are
# evaluated analytically from the Fourier coefficients.
@inline function _polar_radius_derivative(polar::L, phi::T) where {L<:FourierCoeffPolarSegment,T<:Real}
    dr = zero(T)
    sin_coef = polar.coef[1:2:end]
    cos_coef = polar.coef[2:2:end]
    @inbounds for (n,a) in enumerate(cos_coef)
        dr -= n*a*sin(n*phi)
    end
    @inbounds for (n,b) in enumerate(sin_coef)
        dr += n*b*cos(n*phi)
    end
    return dr
end

@inline function _polar_radius_derivative_2(polar::L, phi::T) where {L<:FourierCoeffPolarSegment,T<:Real}
    ddr = zero(T)
    sin_coef = polar.coef[1:2:end]
    cos_coef = polar.coef[2:2:end]
    @inbounds for (n,a) in enumerate(cos_coef)
        ddr -= n^2*a*cos(n*phi)
    end
    @inbounds for (n,b) in enumerate(sin_coef)
        ddr -= n^2*b*sin(n*phi)
    end
    return ddr
end

@inline function tangent(polar::L, t::T) where {L<:FourierCoeffPolarSegment,T<:Real}
    phi = polar.shift_angle + t*polar.arc_angle
    r = polar_radius(polar, phi)
    dr = _polar_radius_derivative(polar, phi)
    dphi = polar.arc_angle
    return dphi*SVector(dr*cos(phi)-r*sin(phi), dr*sin(phi)+r*cos(phi))
end

@inline function tangent_2(polar::L, t::T) where {L<:FourierCoeffPolarSegment,T<:Real}
    phi = polar.shift_angle + t*polar.arc_angle
    r = polar_radius(polar, phi)
    dr = _polar_radius_derivative(polar, phi)
    ddr = _polar_radius_derivative_2(polar, phi)
    dphi2 = polar.arc_angle^2
    return dphi2*SVector(ddr*cos(phi)-2*dr*sin(phi)-r*cos(phi), ddr*sin(phi)+2*dr*cos(phi)-r*sin(phi))
end

# PolarSegment (arbitrary r_func(φ)): no analytic derivative is available, so
# `curve(polar_curve, t)` is differentiated via nested `ForwardDiff` calls.
function tangent(polar_curve::L, t::T) where {L<:PolarSegment,T<:Real}
    dx = ForwardDiff.derivative(u->curve(polar_curve,u)[1], t)
    dy = ForwardDiff.derivative(u->curve(polar_curve,u)[2], t)
    return SVector(dx,dy)
end

function tangent_2(polar_curve::L, t::T) where {L<:PolarSegment,T<:Real}
    ddx = ForwardDiff.derivative(u->ForwardDiff.derivative(v->curve(polar_curve,v)[1], u), t)
    ddy = ForwardDiff.derivative(u->ForwardDiff.derivative(v->curve(polar_curve,v)[2], u), t)
    return SVector(ddx,ddy)
end

"""
    tangent_vec(crv::AbsCurve, ts::AbstractArray) → t̂::Vector

Unit tangent vectors, `t̂(t) = v(t)/norm(v(t))` where `v(t) = tangent(crv,t)`.
"""
function tangent_vec(crv::AbsCurve, ts::AbstractArray{<:Real})
    ta = tangent(crv, ts)
    return [ti/norm(ti) for ti in ta]
end

"""
    normal_vec(crv::AbsCurve, ts::AbstractArray) → n̂::Vector

Unit outward normal vectors, obtained from the unit tangent by a clockwise
90° rotation, `(tx,ty) -> (ty,-tx)`. Valid for a CCW-oriented outer boundary;
hole boundaries need the opposite curve orientation for the same rule to give
the outward normal of the billiard domain.
"""
function normal_vec(crv::AbsCurve, ts::AbstractArray{<:Real})
    ta = tangent_vec(crv, ts)
    return [SVector(ti[2], -ti[1]) for ti in ta]
end

"""
    curvature(crv::AbsCurve, ts::AbstractArray) → κ::Vector
    curvature(crv::AbsCurve, t::Real) → κ::Real

Signed curvature `κ(t) = (x'y''-y'x'')/(x'^2+y'^2)^{3/2}`, reusing
`tangent`/`tangent_2` generically for any `<:AbsCurve`.
"""
function curvature(crv::AbsCurve, ts::AbstractArray{<:Real})
    dr = tangent(crv, ts)
    ddr = tangent_2(crv, ts)
    kappa = similar(ts)
    @inbounds for i in eachindex(ts)
        v = dr[i]
        a = ddr[i]
        den = hypot(v[1], v[2])^3
        kappa[i] = (v[1]*a[2]-v[2]*a[1])/den
    end
    return kappa
end

function curvature(crv::AbsCurve, t::Real)
    dr = tangent(crv, t)
    ddr = tangent_2(crv, t)
    den = hypot(dr[1], dr[2])^3
    return (dr[1]*ddr[2]-dr[2]*ddr[1])/den
end
