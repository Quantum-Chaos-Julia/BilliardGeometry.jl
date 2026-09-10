"""
    PolarSegment{T,BC,F}<:AbsPolarCurve{BC}

Polar boundary segment defined by an arbitrary radial function `r_func(φ)`.

The angular parametrization is

    φ(t)=shift_angle+t*arc_angle,    t∈[0,1].

The function `r_func` returns the complete polar radius `r(φ)`.
"""
struct PolarSegment{T,BC,F}<:AbsPolarCurve{BC} where {T<:Real}
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
    FourierCoeffPolarSegment{T,BC,N}<:AbsPolarCurve{BC}

Polar boundary segment with Fourier radial function

    r(φ)=R+Σₙ[aₙ*cos(nφ)+bₙ*sin(nφ)],

where `coef=[b₁,a₁,b₂,a₂,...]`.
"""
struct FourierCoeffPolarSegment{T,BC,N}<:AbsPolarCurve{BC} where {N<:Int,T<:Real}
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
    PolarSegment(r_func::F;R=1.0,arc_angle=2pi,shift_angle=0.0,center=[0.0,0.0],orientation=1,bc=SpecularReflection(),domain_id=1,segment_id=1) where {F}

Construct a polar segment from the radial function `r_func(φ)`.

## Arguments
* `r_func::F`: Function returning the complete radius `r(φ)`.
* `R`: Reference radial scale.
* `arc_angle`: Angular extent of the segment.
* `shift_angle`: Initial polar angle.
* `center`: Polar-coordinate center.
* `orientation`: Domain orientation, `1` or `-1`.
* `bc`: Boundary condition.
* `domain_id`: Domain identifier.
* `segment_id`: Segment identifier.

## Returns
* `PolarSegment`: Constructed function-defined polar segment.
"""
function PolarSegment(r_func::F;R=1.0,arc_angle=2.0*pi,shift_angle=0.0,center=[0.0,0.0],orientation=1,bc=SpecularReflection(),domain_id=1,segment_id=1) where {F}
    T=promote_type(typeof(R),typeof(arc_angle),typeof(shift_angle),eltype(center))
    polar_curve=PolarSegment{T,typeof(bc),F}(T(R),r_func,T(arc_angle),T(shift_angle),SVector{2,T}(center),Int64(orientation),zero(T),bc,Int64(domain_id),Int64(segment_id))
    L=arc_length(polar_curve,one(T))
    return @set polar_curve.length=L
end

"""
    FourierCoeffPolarSegment(coef;R=1.0,arc_angle=2pi,shift_angle=0.0,center=[0.0,0.0],orientation=1,bc=SpecularReflection(),domain_id=1,segment_id=1)

Construct a Fourier-coefficient polar segment.

## Arguments
* `coef`: Fourier coefficients `[b₁,a₁,b₂,a₂,...]`.
* `R`: Constant radial term.
* `arc_angle`: Angular extent of the segment.
* `shift_angle`: Initial polar angle.
* `center`: Polar-coordinate center.
* `orientation`: Domain orientation, `1` or `-1`.
* `bc`: Boundary condition.
* `domain_id`: Domain identifier.
* `segment_id`: Segment identifier.

## Returns
* `FourierCoeffPolarSegment`: Constructed Fourier polar segment.
"""
function FourierCoeffPolarSegment(coef;R=1.0,arc_angle=2.0*pi,shift_angle=0.0,center=[0.0,0.0],orientation=1,bc=SpecularReflection(),domain_id=1,segment_id=1)
    N=length(coef)
    T=promote_type(eltype(coef),typeof(R),typeof(arc_angle),typeof(shift_angle),eltype(center))
    coefs=SVector{N,T}(coef)
    polar_curve=FourierCoeffPolarSegment{T,typeof(bc),N}(T(R),coefs,T(arc_angle),T(shift_angle),SVector{2,T}(center),Int64(orientation),zero(T),bc,Int64(domain_id),Int64(segment_id))
    L=arc_length(polar_curve,one(T))
    return @set polar_curve.length=L
end

"""
    polar_radius(polar_segment::PolarSegment,φ::T) where {T<:Real}

Evaluate the radial function of a function-defined polar segment.

## Arguments
* `polar_segment::PolarSegment`: Polar segment.
* `φ::T`: Polar angle.

## Returns
* `Real`: Radius `r(φ)`.
"""
@inline function polar_radius(polar_segment::PolarSegment,φ::T) where {T<:Real}
    return polar_segment.r_func(φ)
end

"""
    polar_radius(polar_segment::FourierCoeffPolarSegment,φ::T) where {T<:Real}

Evaluate the Fourier radial function.

## Arguments
* `polar_segment::FourierCoeffPolarSegment`: Fourier polar segment.
* `φ::T`: Polar angle.

## Returns
* `Real`: Radius `r(φ)`.
"""
function polar_radius(polar_segment::L,φ::T) where {L<:FourierCoeffPolarSegment,T<:Real}
    radius=polar_segment.R
    sin_coef=polar_segment.coef[1:2:end]
    cos_coef=polar_segment.coef[2:2:end]
    @inbounds for (n,a) in enumerate(cos_coef)
        radius+=a*cos(n*φ)
    end
    @inbounds for (n,b) in enumerate(sin_coef)
        radius+=b*sin(n*φ)
    end
    return radius
end

"""
    polar_radius(polar_curve::L,φs::AbstractArray) where {L<:AbsPolarCurve}

Evaluate the radial function at multiple polar angles.

## Arguments
* `polar_curve::L`: Polar segment.
* `φs::AbstractArray`: Polar angles.

## Returns
* `Vector`: Radii at the supplied angles.
"""
function polar_radius(polar_curve::L,φs::AbstractArray) where {L<:AbsPolarCurve}
    return [polar_radius(polar_curve,φ) for φ in φs]
end

"""
    curve(polar_curve::L,t::T) where {L<:AbsPolarCurve,T<:Real}

Evaluate a polar segment at parameter `t`.

## Arguments
* `polar_curve::L`: Polar segment.
* `t::T`: Curve parameter.

## Returns
* `SVector{2}`: Cartesian boundary point.
"""
function curve(polar_curve::L,t::T) where {L<:AbsPolarCurve,T<:Real}
    φ=polar_curve.shift_angle+t*polar_curve.arc_angle
    radius=polar_radius(polar_curve,φ)
    pt=Polar(radius,φ)
    return Translation(polar_curve.center)(CartesianFromPolar()(pt))
end

"""
    curve(polar_curve::L,ts::AbstractArray) where {L<:AbsPolarCurve}

Evaluate a polar segment at multiple parameters.

## Arguments
* `polar_curve::L`: Polar segment.
* `ts::AbstractArray`: Curve parameters.

## Returns
* `Vector`: Cartesian boundary points.
"""
function curve(polar_curve::L,ts::AbstractArray) where {L<:AbsPolarCurve}
    return [curve(polar_curve,t) for t in ts]
end

"""
    arc_length(polar_curve::L,t::T) where {L<:AbsPolarCurve,T<:Real}

Compute the arc length from parameter `0` to `t`.

## Arguments
* `polar_curve::L`: Polar segment.
* `t::T`: Upper curve parameter.

## Returns
* `Real`: Arc length.
"""
function arc_length(polar_curve::L,t::T) where {L<:AbsPolarCurve,T<:Real}
    integrand(u)=norm(tangent(polar_curve,u))
    length,_=quadgk(integrand,zero(T),t)
    return length
end

"""
    arc_length(polar_curve::L,ts::AbstractArray) where {L<:AbsPolarCurve}

Compute arc lengths at multiple curve parameters.

## Arguments
* `polar_curve::L`: Polar segment.
* `ts::AbstractArray`: Curve parameters.

## Returns
* `Vector`: Arc lengths from `0` to each parameter.
"""
function arc_length(polar_curve::L,ts::AbstractArray) where {L<:AbsPolarCurve}
    return [arc_length(polar_curve,t) for t in ts]
end

"""
    polar_domain(polar_curve::L,pt) where {L<:AbsPolarCurve}

Evaluate the unsigned polar domain function `ρ-r(φ)`.

## Arguments
* `polar_curve::L`: Polar segment.
* `pt`: Cartesian point.

## Returns
* `Real`: Value `ρ-r(φ)`.
"""
function polar_domain(polar_curve::L,pt) where {L<:AbsPolarCurve}
    pt_polar=PolarFromCartesian()(Translation(-polar_curve.center)(pt))
    R=polar_radius(polar_curve,pt_polar.θ)
    return pt_polar.r-R
end

"""
    domain_fun(polar_curve::L,pt::SVector{2,T}) where {L<:AbsPolarCurve,T<:Real}

Evaluate the oriented polar domain function.

## Arguments
* `polar_curve::L`: Polar segment.
* `pt::SVector{2,T}`: Cartesian point.

## Returns
* `Real`: Oriented domain-function value.
"""
function domain_fun(polar_curve::L,pt::SVector{2,T}) where {L<:AbsPolarCurve,T<:Real}
    return polar_domain(polar_curve,pt)*polar_curve.orientation
end

"""
    domain_fun(polar_curve::L,pts::AbstractArray) where {L<:AbsPolarCurve}

Evaluate the oriented polar domain function at multiple points.

## Arguments
* `polar_curve::L`: Polar segment.
* `pts::AbstractArray`: Cartesian points.

## Returns
* `Vector`: Oriented domain-function values.
"""
function domain_fun(polar_curve::L,pts::AbstractArray) where {L<:AbsPolarCurve}
    return [polar_domain(polar_curve,pt)*polar_curve.orientation for pt in pts]
end