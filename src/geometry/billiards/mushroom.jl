struct MushroomBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::CompositeDomain
    symmetries::SymmetryRegistry
end

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


