
ident = IdentityTransformation()
reflect_x = LinearMap(SMatrix{2,2}([-1.0 0.0;0.0 1.0]))
reflect_y = LinearMap(SMatrix{2,2}([1.0 0.0;0.0 -1.0]))
reflect_xy = reflect_x ∘ reflect_y

struct XAxisReflection <: AbsReflection
end
struct YAxisReflection <: AbsReflection
end
struct XYAxisReflection <: AbsReflection 
end

function apply_symmetry(sym::XAxisReflection, pt::SVector{2,T})
    return reflect_y(pt)
end
function apply_symmetry(sym::XAxisReflection, pts)
    return [reflect_y(pt) for pt in pts]
end

function apply_symmetry(sym::YAxisReflection, pt::SVector{2,T})
    return reflect_x(pt)
end
function apply_symmetry(sym::YAxisReflection, pts)
    return [reflect_x(pt) for pt in pts]
end

function apply_symmetry(sym::XYAxisReflection, pt::SVector{2,T})
    return reflect_xy(pt)
end
function apply_symmetry(sym::XYAxisReflection, pts)
    return [reflect_xy(pt) for pt in pts]
end

D2_symmetry = [XAxisReflection(), XYAxisReflection(), YAxisReflection()]

abstract type AbsRotation <: AbsSymmetry end

struct NFoldRotation <: AbsRotation
    order::Int64
    m::Int64
    angle::Float64
    sym_map::LinearMap{SMatrix{2, 2, Float64, 4}}
end

function NFoldRotation(N,m)
    angle = 2*pi/N
    sym_map = LinearMap(RotZ(angle*m))
    return NFoldRotation(N,m,angle,sym_map) 
end

function apply_symmetry(sym::NFoldRotation, pts)
    return [sym.sym_map(pt) for pt in pts]
end


Cn_symmetry(n) = [NFoldRotation(n,i) for i in 1:(n-1)]