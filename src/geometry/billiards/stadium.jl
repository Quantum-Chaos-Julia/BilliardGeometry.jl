struct Stadium{T} <: AbsBilliard where T<:Real
    subdomains::Vector{AbsDomain}
    symmetries::Vector{CoordinateTransformations.Transformation}
end

function Stadium(half_width)
    x_ref = XReflection(4)
    y_ref = YReflection(4)

    circle = CircleSegment(1.0, pi/2, 0.0, half_width, 0.0)
    x_segment = SymLineSegment(SVector(half_width,0.0),SVector(half_width+1.0,0.0), y_ref)
    y_segment = VirtLineSegment(SVector(half_width, 1.0), SVector(half_width,0.0),2)
    
    t_seg = LineSegment(SVector(half_width,1.0),SVector(0.0,1.0))
    l_seg = SymLineSegment(SVector(0.0,1.0),SVector(0.0,0.0), x_ref)
    b_seg = SymLineSegment(SVector(0.0,0.0),SVector(half_width,0.0) , y_ref)
    r_seg = VirtLineSegment(SVector(half_width,0.0),SVector(half_width,1.0),1)
    
    circle_dom =  Domain{Float64}([circle,x_segment,y_segment],1)
    rectangle_dom = Domain{Float64}([t_seg, l_seg, b_seg, r_seg],2)
    symmetries = [ident, reflect_x, reflect_xy, reflect_y] #order coresponds to symmetry sectors
    return Stadium{typeof(half_width)}([circle_dom,rectangle_dom], symmetries)
end


