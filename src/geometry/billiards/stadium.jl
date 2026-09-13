struct StadiumBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::CompositeDomain
    symmetries::SymmetryRegistry
end

function StadiumBilliard(half_width)
    symmetries = D2_symmetry() # sym_id 1=YAxisReflection, 2=XYAxisReflection, 3=XAxisReflection
    circle_dom = CircleWedge(1.0, pi/2, 0.0, half_width, 0.0, 1; 
    bcs = [SpecularReflection(),Transparent(2),SymmetryWall(3,4)])

    rectangle_dom = Polygon([[half_width,1.0],[0.0,1.0],[0.0,0.0],[half_width,0.0]],2;
    bcs = [SpecularReflection(),SymmetryWall(1,4),SymmetryWall(3,4),Transparent(1)])

    return StadiumBilliard{typeof(half_width)}(CompositeDomain([circle_dom,rectangle_dom]), symmetries)
end


