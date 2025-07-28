

function limacon_eq(a::T, R::T, arc_angle::T, shift_angle::T, center::SVector{2,T}, t) where T<:Real
    #only shift angle = 0 is supported
    s, c = sincos(arc_angle*t+shift_angle)
    return SVector(R*(one(a)+a*c)*c - center[1], R*(one(a)+a*c)*s - center[2])
end


function limacon_arc_length(a, R, phi)
    m = 4*a/(one(a)+a)^2;
    res = Elliptic.E(phi, m)
    return 2.0 *(1.0 + a) * res
end

struct LimaconSegment{T,BC} <: AbsPolarCurve{BC} where T<:Real
    parameter::T
    R::T
    arc_angle::T
    shift_angle::T
    center::SVector{2,T}
    orientation::Int64
    length::T
    cusp::SVector{2,T}
    start::SVector{2,T}
    bc::BC
end

function LimaconSegment(a; orientation=1)
    type = typeof(a)
    R = one(type)
    arc = type(pi)
    shift = zero(type)
    center = SVector(zero(type),zero(type))
    L = limacon_arc_length(a, R, arc)
    bc = SpecularReflection()
    limacon = LimaconSegment(a, R, arc, shift, center, orientation, L, center, center, bc)
    r0, r1 = curve(limacon, [0.0,1.0])
    center = (r0 .+ r1)/2
    #recenter
    limacon = LimaconSegment(a, R, arc, shift, -center, orientation, L, r1, r0, bc)
    r0, r1 = curve(limacon, [0.0,1.0]) 
    limacon = LimaconSegment(a, R, arc, shift, -center, orientation, L, r1, r0, bc)   

    return limacon
end

function polar_radius(limacon::L, phi::T) where {L<:LimaconSegment, T<:Real}
   return limacon.R*(one(phi)+limacon.parameter *cos(phi))
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
    fundamental_domain::PolarDomain
    symmetries::Vector{AbsSymmetry}
end


function Limacon(a)
    limacon = LimaconSegment(a)
    r0 = limacon.cusp
    r1 = limacon.start
    type = typeof(a)
    x_segment = LineSegment(r0, r1; bc=ReflectionSymmetry(XAxisReflection(),2))
    limacon_dom =  PolarDomain{type}([limacon,x_segment],[r1,r0],1)
    symmetries = [XAxisReflection()]
    return Limacon{type}(limacon_dom, symmetries)
end
