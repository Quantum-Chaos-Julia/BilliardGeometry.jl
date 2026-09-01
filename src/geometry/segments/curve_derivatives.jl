"""
    tangent(crv::AbsCurve,ts::AbstractArray)

Evaluate the first parameter derivative of a curve at all parameter values in
`ts`.

For a parametrized boundary `r(t)=(x(t),y(t))`, `tangent(crv,t)` denotes the
derivative `dr/dt`.

This is the derivative with respect to the curve parameter `t`, not a normalized
unit tangent and not a derivative with respect to arclength.

## Arguments
* `crv::AbsCurve`: Boundary curve.
* `ts::AbstractArray`: Parameter values at which to evaluate the derivative.

## Returns
* `tangents`: First parameter derivative at each value in `ts`.
"""
function tangent(crv::AbsCurve,ts::AbstractArray{<:Real})
    return [tangent(crv,t) for t in ts]
end

"""
    tangent_2(crv::AbsCurve,ts::AbstractArray)

Evaluate the second parameter derivative of a curve at all parameter values in
`ts`.

For a parametrized boundary `r(t)`, `tangent_2(crv,t)` denotes the derivative
`d²r/dt²`.

## Arguments
* `crv::AbsCurve`: Boundary curve.
* `ts::AbstractArray`: Parameter values at which to evaluate the derivative.

## Returns
* `tangents_2`: Second parameter derivative at each value in `ts`.
"""
function tangent_2(crv::AbsCurve,ts::AbstractArray{<:Real})
    return [tangent_2(crv,t) for t in ts]
end

"""
    tangent(line::L,t) where {L<:AbsLine}

Return the first parameter derivative of a line segment.

For `r(t)=p₀+t(p₁-p₀)`, the derivative is constant,

    dr/dt = p₁-p₀.

## Arguments
* `line::L`: Line segment.
* `t`: Curve parameter. The result is independent of `t`.

## Returns
* `dr`: Constant first parameter derivative of the line.
"""
@inline function tangent(line::L,t::Real) where {L<:AbsLine}
    return line.pt1-line.pt0
end

"""
    tangent_2(line::L,t) where {L<:AbsLine}

Return the second parameter derivative of a line segment.

Because a line segment is affine in its curve parameter,

    d²r/dt² = 0.

## Arguments
* `line::L`: Line segment.
* `t`: Curve parameter. The result is independent of `t`.

## Returns
* `ddr`: Zero vector with the same coordinate type as the line.
"""
@inline function tangent_2(line::L,t::Real) where {L<:AbsLine}
    return zero(line.pt0)
end

"""
    tangent(circle::L,t) where {L<:CircleSegment}

Return the first parameter derivative of a circular segment.

For

    phi(t) = shift_angle + arc_angle*t,

and

    r(t) = center + radius*(cos(phi),sin(phi)),

the derivative is

    dr/dt = radius*arc_angle*(-sin(phi),cos(phi)).

## Arguments
* `circle::L`: Circular boundary segment.
* `t`: Curve parameter.

## Returns
* `dr`: First parameter derivative `dr/dt`.
"""
@inline function tangent(circle::L,t::Real) where {L<:CircleSegment}
    phi=circle.arc_angle*t+circle.shift_angle
    Ra=circle.radius*circle.arc_angle
    return SVector(-Ra*sin(phi),Ra*cos(phi))
end

"""
    tangent_2(circle::L,t) where {L<:CircleSegment}

Return the second parameter derivative of a circular segment.

For `phi(t)=shift_angle+arc_angle*t`,

    d²r/dt² = radius*arc_angle²*(-cos(phi),-sin(phi)).

## Arguments
* `circle::L`: Circular boundary segment.
* `t`: Curve parameter.

## Returns
* `ddr`: Second parameter derivative `d²r/dt²`.
"""
@inline function tangent_2(circle::L,t::Real) where {L<:CircleSegment}
    phi=circle.arc_angle*t+circle.shift_angle
    Ra2=circle.radius*circle.arc_angle^2
    return SVector(-Ra2*cos(phi),-Ra2*sin(phi))
end

"""
    _polar_radius_derivative(polar::L,phi::T) where {L<:PolarSegment,T<:Real} → T

Return the first angular derivative `dr/dphi` of the Fourier radius stored by a
`PolarSegment`.

For the coefficient convention

    r(phi) = R + Σ aₙ cos(n*phi) + Σ bₙ sin(n*phi),

with

    polar.coef = [b₁,a₁,b₂,a₂,...],

the derivative is evaluated analytically.
"""
@inline function _polar_radius_derivative(polar::L,phi::T) where {L<:PolarSegment,T<:Real}
    dr=zero(T)
    sin_coef=polar.coef[1:2:end]
    cos_coef=polar.coef[2:2:end]
    @inbounds for (n,a) in enumerate(cos_coef)
        dr-=n*a*sin(n*phi)
    end
    @inbounds for (n,b) in enumerate(sin_coef)
        dr+=n*b*cos(n*phi)
    end
    return dr
end

"""
    _polar_radius_derivative_2(polar::L,phi::T) where {L<:PolarSegment,T<:Real} → T

Return the second angular derivative `d²r/dphi²` of the Fourier radius stored by
a `PolarSegment`.

The derivative is evaluated analytically from the Fourier coefficients.
"""
@inline function _polar_radius_derivative_2(polar::L,phi::T) where {L<:PolarSegment,T<:Real}
    ddr=zero(T)
    sin_coef=polar.coef[1:2:end]
    cos_coef=polar.coef[2:2:end]
    @inbounds for (n,a) in enumerate(cos_coef)
        ddr-=n^2*a*cos(n*phi)
    end
    @inbounds for (n,b) in enumerate(sin_coef)
        ddr-=n^2*b*sin(n*phi)
    end
    return ddr
end

"""
    tangent(polar::L,t::T) where {L<:PolarSegment,T<:Real}

Return the first parameter derivative of a Fourier polar segment.

For

    phi(t) = shift_angle + arc_angle*t,

and

    r_vec(t) = center + r(phi)*(cos(phi),sin(phi)),

the derivative is

    dr_vec/dt =
        arc_angle*(
            r'(phi)*cos(phi)-r(phi)*sin(phi),
            r'(phi)*sin(phi)+r(phi)*cos(phi)
        ).

Both `r` and `r'` are evaluated directly from the Fourier representation.

## Arguments
* `polar::L`: Fourier polar boundary segment.
* `t::T`: Curve parameter.

## Returns
* `dr`: First parameter derivative `dr/dt`.
"""
@inline function tangent(polar::L,t::T) where {L<:PolarSegment,T<:Real}
    phi=polar.shift_angle+t*polar.arc_angle
    r=polar_radius(polar,phi)
    dr=_polar_radius_derivative(polar,phi)
    dphi=polar.arc_angle
    return dphi*SVector(dr*cos(phi)-r*sin(phi),dr*sin(phi)+r*cos(phi))
end

"""
    tangent_2(polar::L,t::T) where {L<:PolarSegment,T<:Real}

Return the second parameter derivative of a Fourier polar segment.

For `phi(t)=shift_angle+arc_angle*t`,

    d²r_vec/dt² =
        arc_angle²*(
            r''(phi)*cos(phi)-2*r'(phi)*sin(phi)-r(phi)*cos(phi),
            r''(phi)*sin(phi)+2*r'(phi)*cos(phi)-r(phi)*sin(phi)
        ).

The radius and both angular derivatives are evaluated analytically from the
Fourier coefficients.

## Arguments
* `polar::L`: Fourier polar boundary segment.
* `t::T`: Curve parameter.

## Returns
* `ddr`: Second parameter derivative `d²r/dt²`.
"""
@inline function tangent_2(polar::L,t::T) where {L<:PolarSegment,T<:Real}
    phi=polar.shift_angle+t*polar.arc_angle
    r=polar_radius(polar,phi)
    dr=_polar_radius_derivative(polar,phi)
    ddr=_polar_radius_derivative_2(polar,phi)
    dphi2=polar.arc_angle^2
    return dphi2*SVector(ddr*cos(phi)-2*dr*sin(phi)-r*cos(phi),ddr*sin(phi)+2*dr*cos(phi)-r*sin(phi))
end

"""
    tangent_vec(crv::L,ts::AbstractVector{T}) where {T<:Real,L<:AbsCurve}

Compute the unit tangent vectors of a curve at the parameter values `ts`.

If `v(t)=dr/dt`, the unit tangent is

    t_hat(t) = v(t)/norm(v(t)).

## Arguments
* `crv::L`: Boundary curve implementing [`tangent`](@ref).
* `ts::AbstractVector{T}`: Curve parameter values.

## Returns
* `tangents::Vector{SVector{2,T}}`: Unit tangent vector at each parameter value.
"""
function tangent_vec(crv::L,ts::AbstractVector{T}) where {T<:Real,L<:AbsCurve}
    ta=tangent(crv,ts)
    return [ti/norm(ti) for ti in ta]
end

"""
    normal_vec(crv::L,ts::AbstractVector{T}) where {T<:Real,L<:AbsCurve}

Compute unit normal vectors from the unit tangents by a clockwise rotation
through 90 degrees,

    (t_x,t_y) -> (t_y,-t_x).

For a counterclockwise-oriented outer boundary this gives the outward normal.
Hole boundaries must therefore have the opposite orientation for the same rule
to produce the outward normal of the billiard domain.

## Arguments
* `crv::L`: Boundary curve implementing [`tangent_vec`](@ref).
* `ts::AbstractVector{T}`: Curve parameter values.

## Returns
* `normals::Vector{SVector{2,T}}`: Unit normal vector at each parameter value.
"""
function normal_vec(crv::L,ts::AbstractVector{T}) where {T<:Real,L<:AbsCurve}
    ta=tangent_vec(crv,ts)
    return [SVector(ti[2],-ti[1]) for ti in ta]
end

"""
    curvature(crv::L,ts::AbstractVector{T}) where {T<:Real,L<:AbsCurve} → Vector{T}

Compute the signed curvature of a smooth planar curve at all parameter values
in `ts`.

Using

    r'(t)  = (x'(t),y'(t)),
    r''(t) = (x''(t),y''(t)),

the curvature is

    kappa(t) =
        (x'(t)*y''(t)-y'(t)*x''(t)) /
        (x'(t)^2+y'(t)^2)^(3/2).

The derivatives are obtained from [`tangent`](@ref) and [`tangent_2`](@ref), so
the analytic line, circle and Fourier-polar implementations are reused directly.

The sign follows the orientation of the curve parametrization.

## Arguments
* `crv::L`: Smooth boundary curve implementing `tangent` and `tangent_2`.
* `ts::AbstractVector{T}`: Parameter values at which to evaluate the curvature.

## Returns
* `kappa::Vector{T}`: Signed curvature at each parameter value.
"""
function curvature(crv::L,ts::AbstractVector{T}) where {T<:Real,L<:AbsCurve}
    dr=tangent(crv,ts)
    ddr=tangent_2(crv,ts)
    kappa=similar(ts)
    @inbounds for i in eachindex(ts)
        v=dr[i]
        a=ddr[i]
        den=hypot(v[1],v[2])^3
        kappa[i]=(v[1]*a[2]-v[2]*a[1])/den
    end
    return kappa
end

"""
    curvature(crv::L,t::T) where {T<:Real,L<:AbsCurve} → T

Compute the signed curvature of a smooth planar curve at a single parameter
value `t`.

The curvature is evaluated from the first and second parameter derivatives as

    kappa(t) =
        (x'(t)*y''(t)-y'(t)*x''(t)) /
        (x'(t)^2+y'(t)^2)^(3/2).

The sign follows the orientation of the curve parametrization.

## Arguments
* `crv::L`: Smooth boundary curve implementing `tangent` and `tangent_2`.
* `t::T`: Curve parameter.

## Returns
* `kappa::T`: Signed curvature at `t`.
"""
function curvature(crv::L,t::T) where {T<:Real,L<:AbsCurve}
    dr=tangent(crv,t)
    ddr=tangent_2(crv,t)
    den=hypot(dr[1],dr[2])^3
    return (dr[1]*ddr[2]-dr[2]*ddr[1])/den
end