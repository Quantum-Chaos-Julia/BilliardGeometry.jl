"""
    C3Billiard{T}<:BilliardGeometry.AbsBilliard

Smooth threefold-rotationally symmetric billiard

    r(φ)=scale/2*(1+a*(cos(3φ)-sin(6φ))).

The full physical boundary is one `FourierCoeffPolarSegment`. The fundamental domain is
the sector `0≤φ≤2π/3`, with its radial edges excluded from the physical
boundary discretization.
"""
struct C3Billiard{T}<:BilliardGeometry.AbsBilliard
    fundamental_domain::BilliardGeometry.SimpleDomain{T}
    full_boundary::Vector{BilliardGeometry.AbsCurve}
    symmetries::Vector{BilliardGeometry.AbsSymmetry}
end

"""
    C3Billiard(a::T;scale::T=one(T)) where {T<:Real}

Construct the smooth C₃-symmetric polar billiard

    r(φ)=scale/2*(1+a*(cos(3φ)-sin(6φ))).

## Arguments
* `a::T`: Boundary-deformation amplitude.

## Keyword Arguments
* `scale::T`: Overall radial scale.

## Returns
* `C3Billiard{T}`: Constructed billiard with a `2π/3` fundamental sector.
"""
function C3Billiard(a::T;scale::T=one(T)) where {T<:Real}
    scale>zero(T)||throw(ArgumentError("scale must be positive; received scale=$scale"))
    c=SVector{2,T}(zero(T),zero(T))
    R=scale/T(2)
    amp=a*scale/T(2)
    coef=SVector{12,T}(ntuple(i->i==6 ? amp : i==11 ? -amp : zero(T),12))
    bc=BilliardGeometry.SpecularReflection()
    full=BilliardGeometry.FourierCoeffPolarSegment(coef;R=R,arc_angle=T(2pi),shift_angle=zero(T),center=c,bc=bc,domain_id=1,segment_id=1)
    full_boundary=BilliardGeometry.AbsCurve[full]
    arc=BilliardGeometry.FourierCoeffPolarSegment(coef;R=R,arc_angle=T(2pi/3),shift_angle=zero(T),center=c,bc=bc,domain_id=1,segment_id=1)
    p0=BilliardGeometry.curve(arc,zero(T))
    p1=BilliardGeometry.curve(arc,one(T))
    wall1=BilliardGeometry.LineSegment(p1,c;bc=BilliardGeometry.QuantumSolverIgnore(),domain_id=1,segment_id=2)
    wall0=BilliardGeometry.LineSegment(c,p0;bc=BilliardGeometry.QuantumSolverIgnore(),domain_id=1,segment_id=3)
    fundamental_boundary=BilliardGeometry.AbsCurve[arc,wall1,wall0]
    vertices=SVector{2,T}[p0,p1,c]
    fundamental_domain=BilliardGeometry.SimpleDomain{T}(fundamental_boundary,vertices,1)
    symmetries=BilliardGeometry.AbsSymmetry[]
    append!(symmetries,BilliardGeometry.Cn_symmetry(3))
    return C3Billiard{T}(fundamental_domain,full_boundary,symmetries)
end