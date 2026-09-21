function triangle_corners(angles, x0, y0, h) #x0, y0 position of gamma corner
    alpha, beta, gamma = angles
    B = SVector((-h/tan(beta+alpha))-x0, h-y0)
    A = SVector(h/tan(alpha)+B[1], -y0)
    C =  SVector(-x0, -y0)
    return [A,B,C]
end

struct TriangleBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::Polygon{T}
    symmetries::SymmetryRegistry
end

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