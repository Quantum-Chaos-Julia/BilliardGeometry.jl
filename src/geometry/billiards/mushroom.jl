struct Mushroom{T} <: AbsBilliard where T<:Real
    fundamental_domain::CompositeDomain
    symmetries::Vector{AbsSymmetry}
end

function Mushroom(half_width,stem_heigth=1.0;R=1.0,origin=[0.0,0.0])
    type = typeof(half_width)     
    cx,cy = origin #center of circle segment
    #mushroom cap consists of two domains
    circle = CircleSegment(R, pi/2, 0.0, cx, cy; bc = SpecularReflection())
    chord1 = LineSegment([cx,R], [half_width,cy];bc = Transparent(2))
    x_seg = LineSegment([half_width,cy], [R,cy]; bc = SpecularReflection())
    corners = [SVector{2,type}([R,cy]),SVector{2,type}([cx,R]),SVector{2,type}([half_width,cy])]
    circle_dom =  SimpleDomain{Float64}([x_seg,circle,chord1],corners,1)

    triangle_dom = Polygon([[cx,R],origin,[half_width,cy]],2; 
    bcs=[ReflectionSymmetry(YAxisReflection()),Transparent(3),Transparent(1)])

    #mushroom stem consists of one domain   
    stem_dom = Polygon([origin,[cx,-stem_heigth],[half_width,-stem_heigth],[half_width,cy]],3; 
    bcs=[ReflectionSymmetry(YAxisReflection()),SpecularReflection(),SpecularReflection(),Transparent(2)])

    symmetries = [YAxisReflection()] #order coresponds to symmetry sectors
    return Mushroom{typeof(half_width)}(CompositeDomain([circle_dom, triangle_dom, stem_dom]), symmetries)
end


