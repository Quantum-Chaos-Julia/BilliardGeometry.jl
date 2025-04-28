struct Stadium{T} <: AbsBilliard where T<:Real
    fundamental_domain::CompositeDomain
    symmetries::Vector{AbsSymmetry}
end

function Stadium(half_width)
    circle_dom = CircleWedge(1.0, pi/2, 0.0, half_width, 0.0, 1; 
    bcs = [SpecularReflection(),Transparent(2),ReflectionSymmetry(XAxisReflection(),4)])

    rectangle_dom = Polygon([[half_width,1.0],[0.0,1.0],[0.0,0.0],[half_width,0.0]],2;
    bcs = [SpecularReflection(),ReflectionSymmetry(YAxisReflection(),4),ReflectionSymmetry(XAxisReflection(),4),Transparent(1)])

    symmetries = D2_symmetry #order coresponds to symmetry sectors
    return Stadium{typeof(half_width)}(CompositeDomain([circle_dom,rectangle_dom]), symmetries)
end


