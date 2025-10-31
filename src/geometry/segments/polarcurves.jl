struct PolarSegment{T,BC,N}  <: AbsPolarCurve{BC} where {N<:Int,T<:Real}
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

function PolarSegment(coef; R=1.0, arc_angle =2.0*pi, shift_angle=0.0, center = [0.0,0.0], orientation = 1, bc = SpecularReflection(), domain_id=1, segment_id=1)
    N = length(coef)
    type = eltype(coef)
    polar_curve = PolarSegment{type,typeof(bc),N}(R,coef,arc_angle,shift_angle,center,orientation,0.0,bc,domain_id,segment_id)
    L = arc_length(polar_curve, 1.0)
    return @set polar_curve.length = L
end

function polar_radius(polar_segment::L, phi::T) where {L<:PolarSegment, T<:Real}
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



