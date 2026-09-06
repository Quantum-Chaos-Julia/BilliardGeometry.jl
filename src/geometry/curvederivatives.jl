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

# PolarSegment: r(phi) = R + Σ aₙcos(nφ) + Σ bₙsin(nφ), coef = [b₁,a₁,b₂,a₂,...]
# (same coefficient convention as `polar_radius`). Both angular derivatives are
# evaluated analytically from the Fourier coefficients.
@inline function _polar_radius_derivative(polar::L, phi::T) where {L<:PolarSegment,T<:Real}
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

@inline function _polar_radius_derivative_2(polar::L, phi::T) where {L<:PolarSegment,T<:Real}
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

@inline function tangent(polar::L, t::T) where {L<:PolarSegment,T<:Real}
    phi = polar.shift_angle + t*polar.arc_angle
    r = polar_radius(polar, phi)
    dr = _polar_radius_derivative(polar, phi)
    dphi = polar.arc_angle
    return dphi*SVector(dr*cos(phi)-r*sin(phi), dr*sin(phi)+r*cos(phi))
end

@inline function tangent_2(polar::L, t::T) where {L<:PolarSegment,T<:Real}
    phi = polar.shift_angle + t*polar.arc_angle
    r = polar_radius(polar, phi)
    dr = _polar_radius_derivative(polar, phi)
    ddr = _polar_radius_derivative_2(polar, phi)
    dphi2 = polar.arc_angle^2
    return dphi2*SVector(ddr*cos(phi)-2*dr*sin(phi)-r*cos(phi), ddr*sin(phi)+2*dr*cos(phi)-r*sin(phi))
end
