"""
    FourierCoeffPolarSegment{T,BC,N} <: AbsPolarCurve{BC}

Polar boundary segment with Fourier radial function

```math
r(φ) = R + \\sum_{n} \\left[a_n \\cos(nφ) + b_n \\sin(nφ)\\right],
```

where `coef = [b₁,a₁,b₂,a₂,...]`, `φ = shift_angle + t*arc_angle` and
`t ∈ [0,1]`. First and second derivatives are evaluated analytically, see
[`tangent`](@ref).

## Attributes
* `R::T`: Mean radius.
* `coef::SVector{N,T}`: Interleaved sine and cosine coefficients `[b₁,a₁,b₂,a₂,...]`.
* `arc_angle::T`: Angle swept by the segment.
* `shift_angle::T`: Polar angle of the start point.
* `center::SVector{2,T}`: Origin of the polar coordinates.
* `orientation::Int64`: `+1` or `-1`, selects the inside.
* `length::T`: Arc length, computed by the constructor.
* `bc::BC`: Boundary condition, see [`AbsBoundaryCondition`](@ref).
* `domain_id::Int64`: Id of the domain the segment belongs to.
* `segment_id::Int64`: Id of the segment inside its domain.
"""
struct FourierCoeffPolarSegment{T,BC,N}  <: AbsPolarCurve{BC} where {N<:Int,T<:Real}
    R::T
    coef::SVector{N,T}
    arc_angle::T
    shift_angle::T
    center::SVector{2,T}
    orientation::Int64
    length::T
    bc::BC
    domain_id::Int64
    segment_id::Int64
end

"""
    FourierCoeffPolarSegment(coef; R = 1.0, arc_angle = 2π, shift_angle = 0.0, center = [0.0,0.0], orientation = 1, bc = SpecularReflection(), domain_id = 1, segment_id = 1) → seg::FourierCoeffPolarSegment

Construct a [`FourierCoeffPolarSegment`](@ref). The numeric type is the
promotion of the element types of all numeric arguments and the arc length is
computed by quadrature.

## Arguments
* `coef`: Interleaved coefficients `[b₁,a₁,b₂,a₂,...]`.

## Keyword arguments
* `R = 1.0` : Mean radius.
* `arc_angle = 2π` : Angle swept by the segment.
* `shift_angle = 0.0` : Polar angle of the start point.
* `center = [0.0,0.0]` : Origin of the polar coordinates.
* `orientation = 1` : `+1` or `-1`, selects the inside.
* `bc = SpecularReflection()` : Boundary condition.
* `domain_id = 1` : Id of the domain the segment belongs to.
* `segment_id = 1` : Id of the segment inside its domain.
"""
function FourierCoeffPolarSegment(coef; R=1.0, arc_angle=2.0*pi, shift_angle=0.0, center=[0.0,0.0], orientation=1, bc=SpecularReflection(), domain_id=1, segment_id=1)
    N = length(coef)
    T = promote_type(eltype(coef), typeof(R), typeof(arc_angle), typeof(shift_angle), eltype(center))
    coefs = SVector{N,T}(coef)
    polar_curve = FourierCoeffPolarSegment{T,typeof(bc),N}(T(R), coefs, T(arc_angle), T(shift_angle), SVector{2,T}(center), Int64(orientation), zero(T), bc, Int64(domain_id), Int64(segment_id))
    L = arc_length(polar_curve, one(T))
    return @set polar_curve.length = L
end

"""
    polar_radius(crv::AbsPolarCurve, φ) → r
    polar_radius(crv::AbsPolarCurve, φs::AbstractArray) → rs::Vector

Radial function `r(φ)` of a polar curve at the polar angle(s) `φ`, measured
about the curve `center`. For a [`FourierCoeffPolarSegment`](@ref) it is the
Fourier series, for a [`PolarSegment`](@ref) its `r_func`.
"""
function polar_radius(polar_segment::L, phi::T) where {L<:FourierCoeffPolarSegment, T<:Real}
    let radius = polar_segment.R, sin_coef = polar_segment.coef[1:2:end], cos_coef = polar_segment.coef[2:2:end]
        for (n,a) in enumerate(cos_coef)
            radius = radius + a*cos(phi*n)
        end
        for (n,b) in enumerate(sin_coef)
            radius = radius + b*sin(phi*n)
        end
        return radius
    end
end

"""
    PolarSegment{T,BC,F} <: AbsPolarCurve{BC}

Polar boundary segment defined by an arbitrary radial function `r_func(φ)`
(e.g. an interpolation/spline or any closure), not restricted to a Fourier
series. The parametrization is `φ = shift_angle + t*arc_angle`, `t ∈ [0,1]`.
Derivatives ([`tangent`](@ref), [`tangent_2`](@ref)) are obtained via
`ForwardDiff` rather than an analytic formula.

## Attributes
* `R::T`: Reference radius.
* `r_func::F`: Radial function `r_func(φ)`.
* `arc_angle::T`: Angle swept by the segment.
* `shift_angle::T`: Polar angle of the start point.
* `center::SVector{2,T}`: Origin of the polar coordinates.
* `orientation::Int64`: `+1` or `-1`, selects the inside.
* `length::T`: Arc length, computed by the constructor.
* `bc::BC`: Boundary condition, see [`AbsBoundaryCondition`](@ref).
* `domain_id::Int64`: Id of the domain the segment belongs to.
* `segment_id::Int64`: Id of the segment inside its domain.
"""
struct PolarSegment{T,BC,F}  <: AbsPolarCurve{BC} where {T<:Real}
    R::T
    r_func::F
    arc_angle::T
    shift_angle::T
    center::SVector{2,T}
    orientation::Int64
    length::T
    bc::BC
    domain_id::Int64
    segment_id::Int64
end

"""
    PolarSegment(r_func; R = 1.0, arc_angle = 2π, shift_angle = 0.0, center = [0.0,0.0], orientation = 1, bc = SpecularReflection(), domain_id = 1, segment_id = 1) → seg::PolarSegment

Construct a [`PolarSegment`](@ref) from the radial function `r_func(φ)`. The
numeric type is the promotion of the numeric keyword arguments and the arc
length is computed by quadrature. Keyword arguments are as in
[`FourierCoeffPolarSegment`](@ref).

## Arguments
* `r_func`: Callable returning the radius for a polar angle `φ`.
"""
function PolarSegment(r_func::F; R=1.0, arc_angle=2.0*pi, shift_angle=0.0, center=[0.0,0.0], orientation=1, bc=SpecularReflection(), domain_id=1, segment_id=1) where {F}
    T = promote_type(typeof(R), typeof(arc_angle), typeof(shift_angle), eltype(center))
    polar_curve = PolarSegment{T,typeof(bc),F}(T(R), r_func, T(arc_angle), T(shift_angle), SVector{2,T}(center), Int64(orientation), zero(T), bc, Int64(domain_id), Int64(segment_id))
    L = arc_length(polar_curve, one(T))
    return @set polar_curve.length = L
end

function polar_radius(polar_segment::L, phi::T) where {L<:PolarSegment, T<:Real}
    return polar_segment.r_func(phi)
end

function curve(polar_curve::L, t::T) where {L<:AbsPolarCurve, T<:Real}
    phi = polar_curve.shift_angle + t*polar_curve.arc_angle
    radius = polar_radius(polar_curve, phi)
    pt = Polar(radius, phi) #in polar coordinates
    return Translation(polar_curve.center)(CartesianFromPolar()(pt))
end

function curve(polar_curve::L, ts::AbstractArray) where {L<:AbsPolarCurve}
    phi = @. polar_curve.shift_angle + ts*polar_curve.arc_angle
    radius = polar_radius(polar_curve::L, phi)
    pts_polar = [Polar(r,th) for (r,th) in zip(radius,phi)]
    return [Translation(polar_curve.center)(CartesianFromPolar()(pt)) for pt in pts_polar]
end

#generic functions
function polar_radius(polar_curve::L, phi::AbstractArray) where {L<:AbsPolarCurve}
    return [polar_radius(polar_curve, t) for t in phi]
end

function polar_domain(polar_curve::L, pt) where {L<:AbsPolarCurve}
    pt_polar = PolarFromCartesian()(Translation(-polar_curve.center)(pt))
    R = polar_radius(polar_curve, pt_polar.θ)
    return @. (pt_polar.r - R)
end

# returns negative value inside
function domain_fun(polar_curve::L, pt::SVector{2,T}) where {L<:AbsPolarCurve, T<:Real}
    return polar_domain(polar_curve, pt)*polar_curve.orientation

end

function domain_fun(polar_curve::L, pts::AbstractArray) where {L<:AbsPolarCurve}
    return collect(polar_domain(polar_curve, pt)*polar_curve.orientation for pt in pts)
end



