

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

struct LimaconSegment{T<:Real,BC} <: AbsPolarCurve{BC}
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
    domain_id::Int64
    segment_id::Int64
end

function LimaconSegment(a; orientation=1, bc = SpecularReflection(), domain_id=1, segment_id=1 )
    type = typeof(a)
    R = one(type)
    arc = type(pi)
    shift = zero(type)
    center = SVector(zero(type),zero(type))
    L = 0.0
    
    limacon = LimaconSegment(a, R, arc, shift, center, orientation, L, center, center, bc, domain_id, segment_id)
    r0, r1 = curve(limacon, [0.0,1.0])
    center = (r0 .+ r1)/2
    #recenter
    limacon = LimaconSegment(a, R, arc, shift, -center, orientation, L, r1, r0, bc, domain_id, segment_id)
    r0, r1 = curve(limacon, [0.0,1.0]) 
    limacon = LimaconSegment(a, R, arc, shift, -center, orientation, L, r1, r0, bc, domain_id, segment_id)   
    L = arc_length(limacon, 1.0)
    return @set limacon.length = L
end

function polar_radius(limacon::L, phi::T) where {L<:LimaconSegment, T<:Real}
   return limacon.R*(one(phi)+limacon.parameter *cos(phi))
end

# arc length
#=
function arc_length(limacon::L, pt::SVector{2,T}) where {L<:LimaconSegment, T<:Real}
    let center = limacon.center, a=limacon.parameter, R=limacon.R
        pt_polar = PolarFromCartesian()(Translation(-center)(pt))
        return limacon_arc_length(a,R,pt_polar.θ)
    end
end
=#

####################################################################################

struct LimaconBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::PolarDomain{T}
    symmetries::Vector{AbsSymmetry}
end


function LimaconBilliard(a)
    limacon = LimaconSegment(a)
    r0 = limacon.cusp
    r1 = limacon.start
    type = typeof(a)
    x_segment = LineSegment(r0, r1; bc=ReflectionSymmetry(XAxisReflection(),2),segment_id=2)
    limacon_dom =  PolarDomain{type}([limacon,x_segment],[r1,r0],1)
    symmetries = [XAxisReflection()]
    return LimaconBilliard{type}(limacon_dom, symmetries)
end
