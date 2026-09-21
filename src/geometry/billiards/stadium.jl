struct StadiumBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::CompositeDomain
    symmetries::SymmetryRegistry
end

function StadiumBilliard(half_width::T) where T<:Real
    symmetries = D2_symmetry() # sym_id 1=YAxisReflection, 2=XYAxisReflection, 3=XAxisReflection
    circle_dom = CircleWedge(one(T), T(pi/2), zero(T), half_width, zero(T), 1; 
    bcs = [SpecularReflection(),Transparent(2),SymmetryWall(3,4)])

    rectangle_dom = Polygon([[half_width,one(T)],[zero(T),one(T)],[zero(T),zero(T)],[half_width,zero(T)]],2;
    bcs = [SpecularReflection(),SymmetryWall(1,4),SymmetryWall(3,4),Transparent(1)])

    return StadiumBilliard{T}(CompositeDomain([circle_dom,rectangle_dom]), symmetries)
end


