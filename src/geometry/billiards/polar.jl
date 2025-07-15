struct PolarBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::SimpleDomain
    symmetries::Vector{AbsSymmetry}
end


function PolarBilliard(coef; center=[0.0,0.0])
    type = eltype(coef)
    segment = PolarSegment(coef;center=center )
    r0 = curve(segment, 0.0)
    dom =  SimpleDomain{type}([segment],[r0],1)
    symmetries = []
    return PolarBilliard{type}(dom, symmetries)
end