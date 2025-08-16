include("symmetry.jl")
export XAxisReflection, YAxisReflection, XYAxisReflection, NFoldRotation, apply_symmetry, apply_symmetry_pb, D2_symmetry, Cn_symmetry
include("boundarytypes.jl")
export SpecularReflection, Transparent, PeriodicX, ReflectionSymmetry, get_boundary_curves, get_all_curves, get_curve, get_all_domains, get_domain
include("segments/linesegment.jl")
export LineSegment
include("segments/polarcurves.jl")
export PolarSegment, polar_radius
include("segments/circlesegment.jl")
export CircleSegment
include("segments/compositecurves.jl")
export CompositeCurve
include("domains/polygons.jl")
export Polygon
include("domains/circular.jl")
export CircleWedge
include("domains/compositedomains.jl")
export SimpleDomain, CompositeDomain, reset_ids!
include("billiards/stadium.jl")
export Stadium
include("billiards/mushroom.jl")
export Mushroom
include("billiards/polar.jl")
export PolarBilliard, PolarDomain
include("billiards/limacon.jl")
export Limacon, LimaconSegment
include("poincarebirkhoff.jl")
export PoincareBirkhoff, pb_coords, get_pb_curve, pb_sectors

export is_inside, curve, domain_fun, domain_gradient_vector, arc_length

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

function arc_length(circle::L, pts::AbstractArray) where {L<:AbsCurve}
    return collect(arc_length(circle, pt) for pt in pts)
end


#=
struct Billiard{T} <: AbsBilliard where T<:Real
    subdomains::Vector{AbsDomain}
    #symmetries::
end
=#