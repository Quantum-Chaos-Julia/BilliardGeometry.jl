
"""
    PolarDomain{T} <: AbsPolarDomain

Simple domain bounded by polar curves.

## Attributes
* `boundary::Vector{AbsCurve}`: Boundary curves.
* `corners::Vector{SVector{2,T}}`: Corner points of the boundary.
* `id::Int64`: Id of the domain.
"""
struct PolarDomain{T} <: AbsPolarDomain where T<:Real
    boundary::Vector{AbsCurve}
    corners::Vector{SVector{2,T}}
    id::Int64
end

"""
    PolarBilliard{T} <: AbsBilliard

Billiard bounded by one closed [`FourierCoeffPolarSegment`](@ref), with no
symmetry reduction.

## Attributes
* `fundamental_domain::PolarDomain{T}`: The full domain.
* `symmetries::SymmetryRegistry`: Empty registry.
"""
struct PolarBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::PolarDomain{T}
    symmetries::SymmetryRegistry
end


"""
    PolarBilliard(coef; center = [0.0,0.0]) → billiard::PolarBilliard

Construct a [`PolarBilliard`](@ref) with radial function
`r(φ) = 1 + Σₙ [aₙ cos(nφ) + bₙ sin(nφ)]` over a full turn.

## Arguments
* `coef`: Interleaved Fourier coefficients `[b₁,a₁,b₂,a₂,...]`, see [`FourierCoeffPolarSegment`](@ref). The element type of `coef` is the numeric type of the billiard.

## Keyword arguments
* `center` : Origin of the polar coordinates.
"""
function PolarBilliard(coef; center=[0.0,0.0])
    type = eltype(coef)
    segment = FourierCoeffPolarSegment(coef;center=center )
    r0 = curve(segment, 0.0)
    dom =  PolarDomain{type}([segment],[r0],1)
    symmetries = register_symmetries()
    return PolarBilliard{type}(dom, symmetries)
end