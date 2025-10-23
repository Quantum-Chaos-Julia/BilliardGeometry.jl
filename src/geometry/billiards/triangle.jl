function triangle_corners(angles, x0, y0, h) #x0, y0 position of gamma corner
    alpha, beta, gamma = angles
    B = SVector((-h/tan(beta+alpha))-x0, h-y0)
    A = SVector(h/tan(alpha)+B[1], -y0)
    C =  SVector(-x0, -y0)
    return SVector(A,B,C)
end

struct TriangleBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::SimpleDomain
    symmetries::Vector{AbsSymmetry}
end

function TriangleBilliard(gamma, chi; bcs = [SpecularReflection() for i in corners], x0=zero(gamma), y0=zero(gamma), h = one(gamma))
    alpha = (pi-gamma)/(1+chi)
    beta = alpha*chi
    angles = SVector(alpha, beta, gamma)
    #println("α=$alpha, β=$beta, γ=$gamma")
    corners = triangle_corners(angles, x0, y0, h)
    domain = Polygon(corners, id; bcs)
    symmetries = []
    return TriaTriangleBilliardngle(domain, symmetries)
end