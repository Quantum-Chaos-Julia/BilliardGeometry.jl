struct PolarSegment{N,T,BC}  <: AbsPolarCurve{BC} where {N<:Int,T<:Real}
    R::T
    coef::SVector{N,T}
    arc_angle::T
    shift_angle::T
    center::SVector{2,T}
    orientation::Int64
    length::T
    bc::BC
end

function PolarSegment(coef; R=1.0, arc_angle =2.0*pi, shift_angle=0.0, center = [0.0,0.0], orientation = 1, bc = SpecularReflection())
    N = length(coef)
    type = eltype(coef)
    L = 1.0
    return PolarSegment{N,type,typeof(bc)}(R,coef,arc_angle,shift_angle,center,orientation,L,bc)
end

function polar_radius(R, coef, phi)
    let radius = R, sin_coef =coef[1:2:end], cos_coef = coef[2:2:end]
        for (n,a) in enumerate(cos_coef)
            radius = radius + a*cos(phi*n)
        end
        for (n,b) in enumerate(sin_coef)
            radius = radius + b*sin(phi*n)
        end
        return radius
    end
end

function polar_radius(R, coef, phi::AbstractArray)
    return [polar_radius(R, coef, t) for t in phi]
end

function polar_domain(R, coef, x, y, center)
    angle = atan(y-center[2], x-center[1]) 
    r = polar_radius(R, coef, angle)
    return @. (hypot(y-center[2],x-center[1]) - r)
end



function curve(polar_curve::L, t) where {L<:PolarSegment}
    phi = polar_curve.shift_angle + t*polar_curve.arc_angle
    radius = polar_radius(polar_curve.R, polar_curve.coef, phi)
    return @SVector [radius*cos(phi)+polar_curve.center[1], radius*sin(phi)+polar_curve.center[2]]
end


function curve(polar_curve::L, ts::AbstractArray) where {L<:PolarSegment}
    type = eltype(ts)
    phi = @. polar_curve.shift_angle + ts*polar_curve.arc_angle
    radius = polar_radius(polar_curve.R, polar_curve.coef, phi)
    x0 = polar_curve.center[1]
    y0 = polar_curve.center[2]
    rx = @. radius*cos(phi)+x0
    ry = @. radius*sin(phi)+y0
    return [SVector{2,type}(x,y) for (x,y) in zip(rx,ry)]
end

# returns negative value inside
function domain_fun(polar_curve::L, pt::SVector{2,T}) where {L<:PolarSegment, T<:Real}
    let R = polar_curve.R, coef = polar_curve.coef, center = polar_curve.center 
        return polar_domain(R, coef, pt[1], pt[2], center)*polar_curve.orientation
    end
end

function domain_fun(polar_curve::L, pts::AbstractArray) where {L<:PolarSegment}
    let R = polar_curve.R, coef = polar_curve.coef, center = polar_curve.center 
        return collect(polar_domain(R, coef, pt[1], pt[2], center)*polar_curve.orientation for pt in pts)
    end
end