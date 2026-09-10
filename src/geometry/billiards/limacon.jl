"""
    LimaconBilliard{T}<:BilliardGeometry.AbsBilliard

Limacon billiard with polar boundary

    r(φ)=R*(1+ε*cosφ).

The full physical boundary is represented by one closed
[`BilliardGeometry.FourierCoeffPolarSegment`](@ref). The billiard is symmetric under
reflection across the x-axis, so the fundamental domain is the upper half.

## Attributes
* `fundamental_domain::BilliardGeometry.SimpleDomain{T}`: Upper-half fundamental domain.
* `full_boundary::Vector{BilliardGeometry.AbsCurve}`: Complete physical boundary.
* `symmetries::Vector{BilliardGeometry.AbsSymmetry}`: Geometric symmetries of the billiard.
"""
struct LimaconBilliard{T}<:BilliardGeometry.AbsBilliard
    fundamental_domain::BilliardGeometry.SimpleDomain{T}
    full_boundary::Vector{BilliardGeometry.AbsCurve}
    symmetries::Vector{BilliardGeometry.AbsSymmetry}
end

"""
    LimaconBilliard(ε::T;R::T=one(T),center=SVector{2,T}(zero(T),zero(T))) where {T<:Real} → LimaconBilliard{T}

Construct the Limacon billiard

    r(φ)=R*(1+ε*cosφ).

## Arguments
* `ε::T`: Limacon deformation parameter.

## Keyword Arguments
* `R::T=one(T)`: Mean radius.
* `center::SVector{2,T}`: Billiard center.

## Returns
* `billiard::LimaconBilliard{T}`: Constructed Limacon billiard.
"""
function LimaconBilliard(ε::T;R::T=one(T),center=SVector{2,T}(zero(T),zero(T))) where {T<:Real}
    c=SVector{2,T}(center)
    iszero(c[2])||throw(ArgumentError("x-axis symmetry requires center[2] == 0; received center=$c"))
    coef=SVector{2,T}(zero(T),R*ε)
    bc=BilliardGeometry.SpecularReflection()
    arc=BilliardGeometry.FourierCoeffPolarSegment(coef;R=R,arc_angle=T(pi),shift_angle=zero(T),center=c,bc=bc,domain_id=1,segment_id=1)
    pR=BilliardGeometry.curve(arc,zero(T))
    pL=BilliardGeometry.curve(arc,one(T))
    wall=BilliardGeometry.LineSegment(pL,pR;bc=BilliardGeometry.ReflectionSymmetry(BilliardGeometry.XAxisReflection(),2),domain_id=1,segment_id=2)
    fundamental_boundary=BilliardGeometry.AbsCurve[arc,wall]
    vertices=SVector{2,T}[pR,pL]
    fundamental_domain=BilliardGeometry.SimpleDomain{T}(fundamental_boundary,vertices,1)
    full=BilliardGeometry.FourierCoeffPolarSegment(coef;R=R,arc_angle=T(2pi),shift_angle=zero(T),center=c,bc=bc,domain_id=1,segment_id=1)
    full_boundary=BilliardGeometry.AbsCurve[full]
    symmetries=BilliardGeometry.AbsSymmetry[BilliardGeometry.XAxisReflection()]
    return LimaconBilliard{T}(fundamental_domain,full_boundary,symmetries)
end