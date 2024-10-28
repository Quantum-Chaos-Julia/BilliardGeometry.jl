struct Mushroom{T} <: AbsBilliard where T<:Real
    subdomains::Vector{AbsDomain}
    symmetries::Vector{CoordinateTransformations.Transformation}
end

function Mushroom(half_width,stem_heigth=1.0;R=1.0,origin=SVector(0.0,0.0))
    cx,cy = origin #center of circle segment
    x_ref = XReflection(2)
    
    #mushroom cap consists of two domains
    circle = CircleSegment(R, pi/2, 0.0, cx, cy)
    chord1 = VirtLineSegment(SVector(cx,R), SVector(half_width,cy),2)
    x_seg = LineSegment(SVector(half_width,cy), SVector(R,cy))
    circle_dom =  Domain{Float64}([x_seg,circle,chord1],1)


    y_seg = SymLineSegment(SVector(cx,R), origin, x_ref)
    x_virt_seg = VirtLineSegment(origin, SVector(half_width,cy),3)
    chord2 = VirtLineSegment( SVector(half_width,cy),SVector(cx,R),1)
    triangle_dom =  Domain{Float64}([y_seg,x_virt_seg,chord2],2)

    #mushroom stem consists of one domain
    left_seg = SymLineSegment(origin, SVector(cx,-stem_heigth), x_ref)
    bottom_seg = LineSegment(SVector(cx,-stem_heigth),SVector(half_width,-stem_heigth))
    right_seg = LineSegment(SVector(half_width,-stem_heigth), SVector(half_width,cy))
    top_seg = VirtLineSegment( SVector(half_width,cy), origin, 2)
   
    stem_dom = Domain{Float64}([left_seg, bottom_seg,right_seg, top_seg],3)
    symmetries = [ident, reflect_x] #order coresponds to symmetry sectors
    return Mushroom{typeof(half_width)}([circle_dom, triangle_dom, stem_dom], symmetries)
end


