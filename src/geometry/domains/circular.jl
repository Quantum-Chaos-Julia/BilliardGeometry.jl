
struct CircleCap{T} <: AbsSimpleDomain where T<:Real
    boundary::Vector{AbsCurve}
    corners::Vector{SVector{2,T}}
    id::Int64
end

function CircleCap(pt0,pt1,R)
    midpoint = (pt0 .+ pt1)./2.0
end
