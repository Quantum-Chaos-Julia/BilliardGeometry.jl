"""
    MushroomBilliard{T} <: AbsBilliard

Bunimovich mushroom: a semicircular cap of radius `R` on top of a rectangular
stem of width `2*half_width` and height `stem_heigth`. Reflection about the `y`
axis is the only symmetry, so the fundamental domain is the right half. It is a
[`CompositeDomain`](@ref) of three simple domains: the cap region between the
arc and the chord from `(0,R)` to `(half_width,0)` (`id = 1`), the triangle
below that chord (`id = 2`) and the stem (`id = 3`), glued along
[`Transparent`](@ref) curves.

## Attributes
* `fundamental_domain::CompositeDomain`: Right half of the mushroom.
* `symmetries::SymmetryRegistry`: The single [`YAxisReflection`](@ref).
"""
struct MushroomBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::CompositeDomain
    symmetries::SymmetryRegistry
end

"""
    MushroomBilliard(half_width, stem_heigth = 1; R = 1, origin = zeros(2)) → billiard::MushroomBilliard
    MushroomBilliard(fundamental_domain::CompositeDomain, symmetries::SymmetryRegistry) → billiard::MushroomBilliard

Construct a [`MushroomBilliard`](@ref), or wrap a given composite fundamental
domain and symmetry registry.

## Arguments
* `half_width`: Half the width of the stem.
* `stem_heigth`: Height of the stem.

## Keyword arguments
* `R = 1` : Radius of the cap.
* `origin` : Center of the cap circle; must be the origin because the symmetry is the reflection about the `y` axis.

## Throws
`ArgumentError` if `origin` is not `(0,0)`.
"""
function MushroomBilliard(half_width::T, stem_heigth::T=one(T); R=one(T), origin=zeros(T,2)) where T<:Real
    cx, cy = T(origin[1]), T(origin[2]) #center of circle segment
    (iszero(cx) && iszero(cy)) || throw(ArgumentError("YAxisReflection symmetry requires origin == (0,0); received origin=$origin"))
    Rv = T(R)
    symmetries = register_symmetries(YAxisReflection())
    #mushroom cap consists of two domains
    x_seg = LineSegment(SVector{2,T}(half_width,cy), SVector{2,T}(Rv,cy); bc = SpecularReflection(), domain_id=1, segment_id=1)
    circle = CircleSegment(Rv, T(pi/2), zero(T), cx, cy; bc = SpecularReflection(), domain_id=1, segment_id=2)
    chord1 = LineSegment(SVector{2,T}(cx,Rv), SVector{2,T}(half_width,cy);bc = Transparent(2), domain_id=1, segment_id=3)
    
    corners = SVector{2,T}[SVector{2,T}(Rv,cy),SVector{2,T}(cx,Rv),SVector{2,T}(half_width,cy)]
    circle_dom =  SimpleDomain{T}([x_seg,circle,chord1],corners,1)

    triangle_dom = Polygon([[cx,Rv],[cx,cy],[half_width,cy]],2; 
    bcs=[SymmetryWall(1,2),Transparent(3),Transparent(1)])

    #mushroom stem consists of one domain   
    stem_dom = Polygon([[cx,cy],[cx,-stem_heigth],[half_width,-stem_heigth],[half_width,cy]],3; 
    bcs=[SymmetryWall(1,2),SpecularReflection(),SpecularReflection(),Transparent(2)])

    return MushroomBilliard{T}(CompositeDomain([circle_dom, triangle_dom, stem_dom]), symmetries)
end

MushroomBilliard(fundamental_domain::CompositeDomain, symmetries::SymmetryRegistry) =
    MushroomBilliard{_domain_eltype(fundamental_domain)}(fundamental_domain, symmetries)


