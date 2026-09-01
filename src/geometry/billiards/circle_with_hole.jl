"""
    AnnularBilliard{T}<:BilliardGeometry.AbsBilliard

Concentric annular billiard bounded by two circles with radii `R_outer` and
`R_inner`.

The physical boundary consists of two connected components,

    full_boundary = [outer_boundary,inner_boundary],

where the first component is the outer wall and the second component is the
interior circular obstacle. Boundary orientation of the inner full component is
handled during boundary-point generation.

The fundamental domain is the first-quadrant annular sector bounded by the
outer and inner quarter circles and the two coordinate-axis reflection walls.
"""
struct AnnularBilliard{T}<:BilliardGeometry.AbsBilliard
    fundamental_domain::BilliardGeometry.SimpleDomain{T}
    full_boundary::Vector{Vector{BilliardGeometry.AbsCurve}}
    symmetries::Vector{BilliardGeometry.AbsSymmetry}
end

"""
    AnnularBilliard(R_outer::T,R_inner::T;center=SVector{2,T}(zero(T),zero(T))) where {T<:Real} → AnnularBilliard{T}

Construct a concentric annular billiard with outer radius `R_outer` and inner
radius `R_inner`.

## Arguments
* `R_outer::T`: Radius of the outer circular boundary.
* `R_inner::T`: Radius of the inner circular obstacle.

## Keyword Arguments
* `center::SVector{2,T}`: Common center of the two circular boundaries.

## Returns
* `billiard::AnnularBilliard{T}`: Constructed annular billiard.

## Throws
* `ArgumentError`: If `R_outer≤0`, `R_inner≤0`, or `R_inner≥R_outer`.
"""
function AnnularBilliard(R_outer::T,R_inner::T;center=SVector{2,T}(zero(T),zero(T))) where {T<:Real}
    R_outer>zero(T)||throw(ArgumentError("R_outer must be positive; received $R_outer"))
    R_inner>zero(T)||throw(ArgumentError("R_inner must be positive; received $R_inner"))
    R_inner<R_outer||throw(ArgumentError("R_inner must be smaller than R_outer; received R_inner=$R_inner, R_outer=$R_outer"))
    c=SVector{2,T}(center)
    iszero(c[1])&&iszero(c[2])||throw(ArgumentError("D2 symmetry currently requires center == (0,0); received center=$c"))
    bc=BilliardGeometry.SpecularReflection()
    outer=BilliardGeometry.CircleSegment(R_outer,T(2pi),zero(T),c;bc=bc,domain_id=1,segment_id=1)
    inner=BilliardGeometry.CircleSegment(R_inner,T(2pi),zero(T),c;bc=bc,domain_id=1,segment_id=2)
    full_boundary=[
        BilliardGeometry.AbsCurve[outer],
        BilliardGeometry.AbsCurve[inner]
    ]
    o0=c+SVector{2,T}(R_outer,zero(T))
    o1=c+SVector{2,T}(zero(T),R_outer)
    i0=c+SVector{2,T}(R_inner,zero(T))
    i1=c+SVector{2,T}(zero(T),R_inner)
    outer_q=BilliardGeometry.CircleSegment(R_outer,T(pi/2),zero(T),c;bc=bc,domain_id=1,segment_id=1)
    ywall=BilliardGeometry.LineSegment(o1,i1;bc=BilliardGeometry.ReflectionSymmetry(BilliardGeometry.YAxisReflection(),4),domain_id=1,segment_id=2)
    inner_q=BilliardGeometry.CircleSegment(R_inner,-T(pi/2),T(pi/2),c;orientation=-1,bc=bc,domain_id=1,segment_id=3)
    xwall=BilliardGeometry.LineSegment(i0,o0;bc=BilliardGeometry.ReflectionSymmetry(BilliardGeometry.XAxisReflection(),4),domain_id=1,segment_id=4)
    fundamental_boundary=BilliardGeometry.AbsCurve[outer_q,ywall,inner_q,xwall]
    vertices=SVector{2,T}[o0,o1,i1,i0]
    fundamental_domain=BilliardGeometry.SimpleDomain{T}(fundamental_boundary,vertices,1)
    symmetries=BilliardGeometry.AbsSymmetry[BilliardGeometry.D2_symmetry...]
    return AnnularBilliard{T}(fundamental_domain,full_boundary,symmetries)
end

"""
    CircleStarBilliard{T}<:BilliardGeometry.AbsBilliard

Multiply connected billiard formed by an outer circular wall and a centered
star-shaped polar obstacle

    r(φ)=R_inner+a*cos(nφ).

The physical boundary consists of two connected components,

    full_boundary = [outer_boundary,inner_boundary],

with the star-shaped obstacle represented exactly by a
[`BilliardGeometry.FourierCoeffPolarSegment`](@ref).

The fundamental domain is one rotational sector `0≤φ≤2π/n`, bounded by the
outer circular arc, the inner star-shaped arc, and two radial symmetry walls.

## Attributes
* `fundamental_domain::BilliardGeometry.SimpleDomain{T}`: Rotational fundamental sector.
* `full_boundary::Vector{Vector{BilliardGeometry.AbsCurve}}`: Physical boundary components.
* `symmetries::Vector{BilliardGeometry.AbsSymmetry}`: Geometric symmetries of the billiard.
"""
struct CircleStarBilliard{T}<:BilliardGeometry.AbsBilliard
    fundamental_domain::BilliardGeometry.SimpleDomain{T}
    full_boundary::Vector{Vector{BilliardGeometry.AbsCurve}}
    symmetries::Vector{BilliardGeometry.AbsSymmetry}
end

"""
    CircleStarBilliard(R_outer::T,R_inner::T,a::T,n::Int;center=SVector{2,T}(zero(T),zero(T))) where {T<:Real} → CircleStarBilliard{T}

Construct a circular billiard containing the centered star-shaped obstacle

    r(φ)=R_inner+a*cos(nφ).

## Arguments
* `R_outer::T`: Radius of the outer circular boundary.
* `R_inner::T`: Mean radius of the inner star-shaped obstacle.
* `a::T`: Amplitude of the radial deformation.
* `n::Int`: Number of angular lobes.

## Keyword Arguments
* `center::SVector{2,T}`: Common center of the outer boundary and inner obstacle.

## Returns
* `billiard::CircleStarBilliard{T}`: Constructed circle-star billiard.

## Throws
* `ArgumentError`: If the radii, deformation, or lobe number do not define a
  positive inner obstacle lying strictly inside the outer circle.
"""
function CircleStarBilliard(R_outer::T,R_inner::T,a::T,n::Int;center=SVector{2,T}(zero(T),zero(T))) where {T<:Real}
    R_outer>zero(T)||throw(ArgumentError("R_outer must be positive; received $R_outer"))
    R_inner>zero(T)||throw(ArgumentError("R_inner must be positive; received $R_inner"))
    n>=2||throw(ArgumentError("n must be at least 2; received $n"))
    R_inner-abs(a)>zero(T)||throw(ArgumentError("R_inner-|a| must be positive; received R_inner=$R_inner, a=$a"))
    R_inner+abs(a)<R_outer||throw(ArgumentError("The inner obstacle must lie strictly inside the outer circle; received R_inner+|a|=$(R_inner+abs(a)), R_outer=$R_outer"))
    c=SVector{2,T}(center)
    iszero(c[1])&&iszero(c[2])||throw(ArgumentError("rotational symmetry currently requires center == (0,0); received center=$c"))
    coef=SVector{2*n,T}(ntuple(i->i==2*n ? a : zero(T),2*n))
    bc=BilliardGeometry.SpecularReflection()
    outer=BilliardGeometry.CircleSegment(R_outer,T(2pi),zero(T),c;bc=bc,domain_id=1,segment_id=1)
    inner=BilliardGeometry.FourierCoeffPolarSegment(coef;R=R_inner,arc_angle=T(2pi),shift_angle=zero(T),center=c,bc=bc,domain_id=1,segment_id=2)
    full_boundary=[
        BilliardGeometry.AbsCurve[outer],
        BilliardGeometry.AbsCurve[inner]
    ]
    θ=T(2pi/n)
    o0=c+SVector{2,T}(R_outer,zero(T))
    o1=c+SVector{2,T}(R_outer*cos(θ),R_outer*sin(θ))
    inner_outer=BilliardGeometry.FourierCoeffPolarSegment(coef;R=R_inner,arc_angle=θ,shift_angle=zero(T),center=c,bc=bc,domain_id=1,segment_id=1)
    i0=BilliardGeometry.curve(inner_outer,zero(T))
    i1=BilliardGeometry.curve(inner_outer,one(T))
    outer_q=BilliardGeometry.CircleSegment(R_outer,θ,zero(T),c;bc=bc,domain_id=1,segment_id=1)
    wall1=BilliardGeometry.LineSegment(o1,i1;bc=BilliardGeometry.QuantumSolverIgnore(),domain_id=1,segment_id=2)
    inner_q=BilliardGeometry.FourierCoeffPolarSegment(coef;R=R_inner,arc_angle=-θ,shift_angle=θ,center=c,orientation=-1,bc=bc,domain_id=1,segment_id=3) # negative arc reverses traversal; orientation=-1 marks the obstacle side
    wall0=BilliardGeometry.LineSegment(i0,o0;bc=BilliardGeometry.QuantumSolverIgnore(),domain_id=1,segment_id=4)
    fundamental_boundary=BilliardGeometry.AbsCurve[outer_q,wall1,inner_q,wall0]
    vertices=SVector{2,T}[o0,o1,i1,i0]
    fundamental_domain=BilliardGeometry.SimpleDomain{T}(fundamental_boundary,vertices,1)
    symmetries=BilliardGeometry.AbsSymmetry[BilliardGeometry.XAxisReflection()]
    if iseven(n)
        push!(symmetries,BilliardGeometry.YAxisReflection())
        push!(symmetries,BilliardGeometry.XYAxisReflection())
    end
    append!(symmetries,BilliardGeometry.Cn_symmetry(n))
    return CircleStarBilliard{T}(fundamental_domain,full_boundary,symmetries)
end