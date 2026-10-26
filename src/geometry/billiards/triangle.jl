function triangle_corners(angles, x0, y0, h) #x0, y0 position of gamma corner
    alpha, beta, gamma = angles
    B = SVector((-h/tan(beta+alpha))-x0, h-y0)
    A = SVector(h/tan(alpha)+B[1], -y0)
    C =  SVector(-x0, -y0)
    return [A,B,C]
end

"""
    TriangleBilliard{T} <: AbsBilliard

Triangle with interior angles `α`, `β`, `γ` at the corners `A`, `B`, `C`
(`C` is the corner with angle `γ`), parametrized by `γ` and the ratio
`χ = β/α`, so that `α = (π - γ)/(1 + χ)` and `β = χα`. The fundamental domain is
the full triangle, a [`Polygon`](@ref), and there are no symmetries.

## Attributes
* `fundamental_domain::Polygon{T}`: The triangle with edges `A→B`, `B→C`, `C→A`.
* `symmetries::SymmetryRegistry`: Empty registry.
"""
struct TriangleBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::Polygon{T}
    symmetries::SymmetryRegistry
end

"""
    TriangleBilliard(gamma, chi; bcs = [SpecularReflection(), SpecularReflection(), SpecularReflection()], x0 = 0, y0 = 0, h = 1) → billiard::TriangleBilliard

Construct a [`TriangleBilliard`](@ref).

## Arguments
* `gamma`: Angle at corner `C`, `0 < gamma < π`.
* `chi`: Ratio `β/α` of the other two angles, positive.

## Keyword arguments
* `bcs` : Boundary conditions of the three edges. Use [`QuantumSolverIgnore`](@ref) to turn edges into virtual ones.
* `x0`, `y0` : The corner `C` is placed at `(-x0, -y0)`.
* `h = 1` : Height of corner `B` above the edge `C→A`, which fixes the scale.

## Throws
`ArgumentError` if the angles are not admissible or the arguments are not real.
"""
function TriangleBilliard(gamma, chi; bcs = [SpecularReflection(),SpecularReflection(),SpecularReflection()], x0=zero(gamma), y0=zero(gamma), h = one(gamma))
    T = promote_type(typeof(gamma), typeof(chi))
    T <: Real || throw(ArgumentError("gamma and chi must be Real; received types $(typeof(gamma)), $(typeof(chi))"))
    gamma, chi = T(gamma), T(chi)
    zero(T) < gamma < T(pi) || throw(ArgumentError("gamma must satisfy 0<gamma<pi; received gamma=$gamma"))
    chi > zero(T) || throw(ArgumentError("chi must be positive (so alpha,beta stay positive); received chi=$chi"))
    alpha = (T(pi)-gamma)/(1+chi)
    beta = alpha*chi
    angles = SVector(alpha, beta, gamma)
    #println("α=$alpha, β=$beta, γ=$gamma")
    corners = triangle_corners(angles, x0, y0, h)
    domain = Polygon(corners, 1; bcs)
    symmetries = register_symmetries()
    return TriangleBilliard{T}(domain, symmetries)
end