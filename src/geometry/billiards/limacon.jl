

function limacon_eq(a::T, R::T, arc_angle::T, shift_angle::T, center::SVector{2,T}, t) where T<:Real
    #only shift angle = 0 is supported
    s, c = sincos(arc_angle*t+shift_angle)
    return SVector(R*(one(a)+a*c)*c - center[1], R*(one(a)+a*c)*s - center[2])
end

function limacon_domain(a, R, center, x, y, cusp; m=10.0)
    angle = atan(y+center[2], x+center[1]) 
    
    c = cos(angle)
    r = R*(one(a)+a*c)
    if x < cusp[1] && y<(cusp[2]-1e-12)
        return -m*y
    end
   
    return @. m*(hypot(y+center[2],x+center[1]) - r)
end

function limacon_arc_length(a, R, phi)
    m = 4*a/(one(a)+a)^2;
    res = Elliptic.E(phi, m)
    return 2.0 *(1.0 + a) * res
end

struct LimaconSegment{T} <: AbsCurve where T<:Real
    parameter::T
    radius::T
    arc_angle::T
    shift_angle::T
    center::SVector{2,T}
    orientation::Int64
    length::T
    cusp::SVector{2,T}
    start::SVector{2,T}
end

function LimaconSegment(a; orientation=1)
    type = typeof(a)
    R = one(type)
    arc = type(pi)
    shift = zero(type)
    center = SVector(zero(type),zero(type))
    L = limacon_arc_length(a, R, arc)
    
    r1 = limacon_eq(a, R, arc, shift, center, 0.0) 
    r0 = limacon_eq(a, R, arc, shift, center, 1.0)
    #limacon = LimaconSegment(a, R, arc, shift, center, orientation, L, r1)
    #recenter
    center = (r0 .+ r1)/2
    r1 = limacon_eq(a, R, arc, shift, center, 0.0) 
    r0 = limacon_eq(a, R, arc, shift, center, 1.0)
    return LimaconSegment(a, R, arc, shift, center, orientation, L, r0, r1)
end

function curve(limacon::L, t) where {L<:LimaconSegment}
    let a = limacon.parameter , R = limacon.radius, c = limacon.center, arc=limacon.arc_angle, s=limacon.shift_angle 
        return limacon_eq(a, R, arc, s, c, t)
    end
end
function curve(limacon::L, ts::AbstractArray) where {L<:LimaconSegment}
    let a = limacon.parameter , R = limacon.radius, c = limacon.center, arc=limacon.arc_angle, s=limacon.shift_angle 
        return collect(limacon_eq(a, R, arc, s, c, t) for t in ts)
    end
end

# returns negative value inside
function domain_fun(limacon::L, pt::SVector{2,T}) where {L<:LimaconSegment, T<:Real}
    let a = limacon.parameter , R = limacon.radius, c = limacon.center, arc=limacon.arc_angle, s=limacon.shift_angle 
        return limacon_domain(a, R, c, pt[1], pt[2], limacon.cusp)*limacon.orientation
    end
end

function domain_fun(limacon::L, pts::AbstractArray) where {L<:LimaconSegment}
    let a = limacon.parameter , R = limacon.radius, c = limacon.center, arc=limacon.arc_angle, s=limacon.shift_angle 
        return collect(limacon_domain(a, R, c, pt[1], pt[2], limacon.cusp)*limacon.orientation for pt in pts)
    end
end

# arc length
function arc_length(limacon::L, pt::SVector{2,T}) where {L<:LimaconSegment, T<:Real}
    let center = limacon.center, a=limacon.a, R=limacon.R
        angle = atan(pt[2]-center[2], pt[1]-center[1]) - limacon.shift_angle
        return limacon_arc_length(a,R,angle)
    end
end


####################################################################################

struct Limacon{T} <: AbsBilliard where T<:Real
    parameter::T
    subdomains::Vector{AbsDomain}
    symmetries::Vector{CoordinateTransformations.Transformation}
end


function Limacon(a)
    y_ref = YReflection(2)
    limacon = LimaconSegment(a)
    r0 = limacon.cusp
    r1 = limacon.start
    type = typeof(a)
    x_segment= SymLineSegment(r0, r1, y_ref)
    limacon_dom =  Domain{type}([limacon,x_segment],1)
    symmetries = [ident, reflect_y]
    return Limacon{type}(a, [limacon_dom], symmetries)
end
