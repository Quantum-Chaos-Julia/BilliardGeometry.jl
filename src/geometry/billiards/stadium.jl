"""
    StadiumBilliard{T} <: AbsBilliard

Bunimovich stadium: a rectangle of length `2*half_width` and height `2` capped
by two semicircles of radius `1`. The fundamental domain is the first quadrant,
a [`CompositeDomain`](@ref) of a [`CircleWedge`](@ref) (`id = 1`) and a
[`Polygon`](@ref) (`id = 2`), with the `D₂` symmetry group of
[`D2_symmetry`](@ref).

## Attributes
* `fundamental_domain::CompositeDomain`: The quadrant.
* `symmetries::SymmetryRegistry`: Reflections `YAxisReflection`, `XYAxisReflection`, `XAxisReflection` with `sym_id` `1, 2, 3`.
"""
struct StadiumBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::CompositeDomain
    symmetries::SymmetryRegistry
end

"""
    StadiumBilliard(half_width) → billiard::StadiumBilliard
    StadiumBilliard(fundamental_domain::CompositeDomain, symmetries::SymmetryRegistry) → billiard::StadiumBilliard

Construct a [`StadiumBilliard`](@ref), or wrap a given composite fundamental
domain and symmetry registry.

## Arguments
* `half_width`: Half the length of the straight sections.
"""
function StadiumBilliard(half_width::T) where T<:Real
    symmetries = D2_symmetry() # sym_id 1=YAxisReflection, 2=XYAxisReflection, 3=XAxisReflection
    circle_dom = CircleWedge(one(T), T(pi/2), zero(T), half_width, zero(T), 1; 
    bcs = [SpecularReflection(),Transparent(2),SymmetryWall(3,4)])

    rectangle_dom = Polygon([[half_width,one(T)],[zero(T),one(T)],[zero(T),zero(T)],[half_width,zero(T)]],2;
    bcs = [SpecularReflection(),SymmetryWall(1,4),SymmetryWall(3,4),Transparent(1)])

    return StadiumBilliard{T}(CompositeDomain([circle_dom,rectangle_dom]), symmetries)
end

StadiumBilliard(fundamental_domain::CompositeDomain, symmetries::SymmetryRegistry) =
    StadiumBilliard{_domain_eltype(fundamental_domain)}(fundamental_domain, symmetries)


