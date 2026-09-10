"""
    FourierCoeffPolarSegment{T,BC,N} <: AbsPolarCurve{BC}

Polar boundary segment with Fourier radial function
`r(φ) = R + Σₙ[aₙ*cos(nφ)+bₙ*sin(nφ)]`, where `coef = [b₁,a₁,b₂,a₂,...]`.
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

function FourierCoeffPolarSegment(coef; R=1.0, arc_angle=2.0*pi, shift_angle=0.0, center=[0.0,0.0], orientation=1, bc=SpecularReflection(), domain_id=1, segment_id=1)
    N = length(coef)
    T = promote_type(eltype(coef), typeof(R), typeof(arc_angle), typeof(shift_angle), eltype(center))
    coefs = SVector{N,T}(coef)
    polar_curve = FourierCoeffPolarSegment{T,typeof(bc),N}(T(R), coefs, T(arc_angle), T(shift_angle), SVector{2,T}(center), Int64(orientation), zero(T), bc, Int64(domain_id), Int64(segment_id))
    L = arc_length(polar_curve, one(T))
    return @set polar_curve.length = L
end

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
series. Derivatives (`tangent`/`tangent_2`) are obtained via `ForwardDiff`
rather than an analytic formula (see `curvederivatives.jl`).
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



