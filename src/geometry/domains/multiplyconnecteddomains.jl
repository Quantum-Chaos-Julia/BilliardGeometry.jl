#Multiply connected planar domain (a region with one or more holes). Unlike
#AbsCompositeDomain (whose is_inside/get_boundary_curves semantics are a
#*union* over subdomains, appropriate for overlapping-piece shapes like
#StadiumBilliard), a domain with holes needs an *intersection* semantics
#("inside the outer curve AND outside every hole"), which the already
#AbsDomain-generic is_inside(domain, pt) = all(is_inside(crv, pt) for crv in
#domain.boundary) gives for free as long as hole curves are stored with
#reversed (orientation=-1) curves in the same flat `boundary` vector as the
#outer curve(s) - exactly like SimpleDomain, just under a distinct type name
#so the hole topology is explicit rather than an orientation-trick a reader
#has to already know to reproduce correctly.
struct MultiplyConnectedDomain{T} <: AbsMultiplyConnectedDomain where T<:Real
    boundary::Vector{AbsCurve}
    corners::Vector{SVector{2,T}}
    id::Int64
    genus::Int64
end

"""
    MultiplyConnectedDomain(outer_curves, hole_curve_groups, corners, id) → domain::MultiplyConnectedDomain

Constructs a [`MultiplyConnectedDomain`](@ref) from the outer boundary's
curves, one curve group per hole, and the domain's corner vertices.

## Arguments
* `outer_curves`: Vector of `AbsCurve` making up the outer boundary, listed first so that [`get_boundary_curves`](@ref)/`full_boundary` return the outer boundary before any hole (outermost-first, by construction).
* `hole_curve_groups`: `Vector{<:Vector{<:AbsCurve}}`, one curve group per hole.
* `corners`: Corner vertices of the domain.
* `id`: Integer id of the domain.

## Returns
* `domain`: A [`MultiplyConnectedDomain`](@ref) whose `genus` is `length(hole_curve_groups)`.
"""
function MultiplyConnectedDomain(outer_curves::Vector{<:AbsCurve}, hole_curve_groups::Vector{<:Vector{<:AbsCurve}}, corners::Vector{SVector{2,T}}, id::Int64) where T<:Real
    boundary = AbsCurve[outer_curves...; reduce(vcat, hole_curve_groups; init=AbsCurve[])...]
    genus = length(hole_curve_groups)
    return MultiplyConnectedDomain{T}(boundary, corners, id, genus)
end
