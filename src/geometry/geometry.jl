include("symmetry.jl")
export XReflection, YReflection, XYReflection, new_sector, symmetry_action
include("linesegment.jl")
export LineSegment, VirtLineSegment, SymLineSegment
include("circlesegment.jl")
export CircleSegment
include("stadium.jl")
export Stadium
include("mushroom.jl")
export Mushroom
include("limacon.jl")
export Limacon, LimaconSegment


export Domain, Billiard, is_inside, curve, domain_fun, domain_gradient_vector
struct Domain{T} <: AbsDomain where T<:Real
    boundary::Vector{AbsCurve}
    id::Int64
end

function is_inside(domain::D, pt::SVector{2,T}) where {D<:AbsDomain, T<:Real}
    d = [is_inside(crv, pt) for crv in domain.boundary]
    return d
end

function is_inside(domain::D, pts::AbstractArray) where {D<:AbsDomain}
    d = [is_inside(crv, pts) for crv in domain.boundary]
    return  reduce(hcat,d)
end


#check if points inside for general curves
function is_inside(curve::C, pt::SVector{2,T}) where {C<:AbsCurve, T<:Real}
    return domain_fun(curve, pt) .< zero(eltype(pt)) 
end

function is_inside(curve::C, pts::AbstractArray) where {C<:AbsCurve}
    let
    d = domain_fun(curve, pts)
    return d .< zero(eltype(pts[1])) 
    end
end

#gradient of domain_function gives normal direcrion
function domain_gradient_vector(curve::C, pt::SVector{2,T}) where {C<:AbsCurve, T<:Real}
    f(r) = domain_fun(curve, r)
    g = ForwardDiff.gradient(f, pt)
    return g
end

function domain_gradient_vector(curve::C, pts::AbstractArray) where {C<:AbsCurve}
    f(r) = domain_fun(curve, r)
    gs = [ForwardDiff.gradient(f, pt) for pt in pts]
    return gs
end

#=
struct Billiard{T} <: AbsBilliard where T<:Real
    subdomains::Vector{AbsDomain}
    #symmetries::
end
=#