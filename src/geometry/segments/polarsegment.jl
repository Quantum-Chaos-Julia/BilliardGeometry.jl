
struct PolarSegment{T<:Real,BC,F}<:AbsPolarCurve{BC}
    rfun::F
    arc_angle::T
    shift_angle::T
    center::SVector{2,T}
    orientation::Int64
    length::T
    bc::BC
    domain_id::Int64
    segment_id::Int64
end

function PolarSegment(::Type{T},rfun::F;arc_angle=T(2*π),shift_angle=zero(T),center=SVector{2,T}(zero(T),zero(T)),orientation=1,bc=SpecularReflection(),domain_id=1,segment_id=1) where {T<:Real,F}
    seg=PolarSegment{T,typeof(bc),F}(rfun,T(arc_angle),T(shift_angle),center,orientation,zero(T),bc,domain_id,segment_id)
    L=arc_length(seg,one(T))
    return @set seg.length=L
end

@inline polar_radius(seg::PolarSegment{T},φ) where {T<:Real}=seg.rfun(φ)
@inline polar_radius_derivative(crv::AbsPolarCurve,φ)=ForwardDiff.derivative(Base.Fix1(polar_radius,crv),φ)

function PolarSegment(coef::AbstractVector{T};R=one(T),arc_angle=T(2*π),shift_angle=zero(T),center=SVector{2,T}(zero(T),zero(T)),orientation=1,bc=SpecularReflection(),domain_id=1,segment_id=1) where {T<:Real}
    @assert iseven(length(coef)) "coef must be [sin1,cos1,sin2,cos2,...] (even length)"
    c=SVector{length(coef),T}(coef)
    R=T(R)
    rfun=φ->begin
        r=R
        @inbounds for n in 1:length(c)÷2
            r+=c[2*n]*cos(n*φ)
            r+=c[2*n-1]*sin(n*φ)
        end
        return r
    end
    return PolarSegment(T,rfun;arc_angle=arc_angle,shift_angle=shift_angle,center=center,orientation=orientation,bc=bc,domain_id=domain_id,segment_id=segment_id)
end

@inline function curve(polar_curve::L,t) where {L<:AbsPolarCurve}
    phi=polar_curve.shift_angle+t*polar_curve.arc_angle
    radius=polar_radius(polar_curve,phi)
    pt=Polar(radius,phi) #in polar coordinates
    return Translation(polar_curve.center)(CartesianFromPolar()(pt))
end

function curve(crv::AbsPolarCurve,ts::AbstractArray)
    return collect(curve(crv,t) for t in ts)
end

function polar_radius(crv::AbsPolarCurve,phis::AbstractArray)
    T=eltype(phis)
    return collect(polar_radius(crv,T(φ)) for φ in phis)
end

@inline function polar_domain(polar_curve::L,pt) where {L<:AbsPolarCurve}
    pt_polar=PolarFromCartesian()(Translation(-polar_curve.center)(pt))
    R=polar_radius(polar_curve,pt_polar.θ)
    return @. (pt_polar.r-R)
end

# returns negative value inside
@inline function domain_fun(polar_curve::L,pt::SVector{2,T}) where {L<:AbsPolarCurve,T<:Real}
    return polar_domain(polar_curve,pt)*polar_curve.orientation
end

function domain_fun(polar_curve::L,pts::AbstractArray) where {L<:AbsPolarCurve}
    return collect(polar_domain(polar_curve,pt)*polar_curve.orientation for pt in pts)
end

@inline function tangent(crv::AbsPolarCurve,t)
    φ=crv.shift_angle+t*crv.arc_angle
    r=polar_radius(crv,φ)
    rp=polar_radius_derivative(crv,φ)
    s,c=sincos(φ)
    return crv.arc_angle*SVector(rp*c-r*s,rp*s+r*c)
end

function tangent(crv::AbsPolarCurve,ts::AbstractArray)
    T=eltype(ts)
    return collect(tangent(crv,T(t)) for t in ts)
end




