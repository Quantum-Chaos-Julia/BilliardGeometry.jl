"""
    StarBilliard{T}<:BilliardGeometry.AbsBilliard

Simply connected star-shaped billiard

    r(φ)=R+a*cos(nφ),

with `n`-fold rotational symmetry. The full boundary is one smooth closed
`FourierCoeffPolarSegment`, while the fundamental domain is the sector `0≤φ≤2π/n`.
"""
struct StarBilliard{T}<:BilliardGeometry.AbsBilliard
    fundamental_domain::BilliardGeometry.SimpleDomain{T}
    full_boundary::Vector{BilliardGeometry.AbsCurve}
    symmetries::Vector{BilliardGeometry.AbsSymmetry}
end

"""
    StarBilliard(R::T,a::T,n::Int;center=SVector{2,T}(zero(T),zero(T))) where {T<:Real}

Construct the smooth star billiard `r(φ)=R+a*cos(nφ)`. Requires `n≥2` and
`R>|a|`, so the radial function remains positive.

## Arguments
* `R::T`: Mean radius.
* `a::T`: Radial deformation amplitude.
* `n::Int`: Rotational symmetry order.

## Keyword Arguments
* `center::SVector{2,T}`: Billiard center.

## Returns
* `StarBilliard{T}`: Constructed star billiard.
"""
function StarBilliard(R::T,a::T,n::Int;center=SVector{2,T}(zero(T),zero(T))) where {T<:Real}
    n>=2||throw(ArgumentError("n must satisfy n≥2; received n=$n"))
    R>abs(a)||throw(ArgumentError("Require R>|a| so r(φ)>0; received R=$R, a=$a"))
    c=SVector{2,T}(center)
    coef=SVector{2*n,T}(ntuple(i->i==2*n ? a : zero(T),2*n))
    bc=BilliardGeometry.SpecularReflection()
    full=BilliardGeometry.FourierCoeffPolarSegment(coef;R=R,arc_angle=T(2pi),shift_angle=zero(T),center=c,bc=bc,domain_id=1,segment_id=1)
    full_boundary=BilliardGeometry.AbsCurve[full]
    arc=BilliardGeometry.FourierCoeffPolarSegment(coef;R=R,arc_angle=T(2pi/n),shift_angle=zero(T),center=c,bc=bc,domain_id=1,segment_id=1)
    p0=BilliardGeometry.curve(arc,zero(T))
    p1=BilliardGeometry.curve(arc,one(T))
    wall1=BilliardGeometry.LineSegment(p1,c;bc=BilliardGeometry.QuantumSolverIgnore(),domain_id=1,segment_id=2)
    wall0=BilliardGeometry.LineSegment(c,p0;bc=BilliardGeometry.QuantumSolverIgnore(),domain_id=1,segment_id=3)
    fundamental_boundary=BilliardGeometry.AbsCurve[arc,wall1,wall0]
    vertices=SVector{2,T}[p0,p1,c]
    fundamental_domain=BilliardGeometry.SimpleDomain{T}(fundamental_boundary,vertices,1)
    symmetries=BilliardGeometry.AbsSymmetry[BilliardGeometry.XAxisReflection()]
    if iseven(n)
        push!(symmetries,BilliardGeometry.YAxisReflection())
        push!(symmetries,BilliardGeometry.XYAxisReflection())
    end
    append!(symmetries,BilliardGeometry.Cn_symmetry(n))
    return StarBilliard{T}(fundamental_domain,full_boundary,symmetries)
end