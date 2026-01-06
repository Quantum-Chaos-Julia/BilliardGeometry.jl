
@inline _x_coord_reflect(pt::SVector{2,T},sx::T) where {T<:Real}=SVector{2,T}((sx+sx)-pt[1],pt[2])
@inline _y_coord_reflect(pt::SVector{2,T},sy::T) where {T<:Real}=SVector{2,T}(pt[1],(sy+sy)-pt[2])
@inline _xy_coord_reflect(pt::SVector{2,T},sx::T,sy::T) where {T<:Real}=SVector{2,T}((sx+sx)-pt[1],(sy+sy)-pt[2])

@inline _x_normal_reflect(n::SVector{2,T}) where {T<:Real}=SVector{2,T}(-n[1],n[2])
@inline _y_normal_reflect(n::SVector{2,T}) where {T<:Real}=SVector{2,T}(n[1],-n[2])
@inline _xy_normal_reflect(n::SVector{2,T}) where {T<:Real}=-n

@inline function _rotate(pt::SVector{2,T},angle::T,c0::SVector{2,T}) where {T<:Real}
    dx=pt[1]-c0[1];dy=pt[2]-c0[2]
    sl,cl=sincos(angle)
    return SVector(c0[1]+cl*dx-sl*dy,c0[2]+sl*dx+cl*dy)
end
function _rotate(pts::AbstractVector{SVector{2,T}},angle::T,c0::SVector{2,T}) where {T<:Real}
    sl,cl=sincos(angle)
    return [SVector(c0[1]+cl*(pt[1]-c0[1])-sl*(pt[2]-c0[2]),c0[2]+sl*(pt[1]-c0[1])+cl*(pt[2]-c0[2])) for pt in pts]
    
end
@inline function _rotate_normal(pt::SVector{2,T},angle::T) where {T<:Real}
    sl,cl=sincos(angle)
    return SVector(cl*pt[1]-sl*pt[2],sl*pt[1]+cl*pt[2])
end
function _rotate_normal(pts::AbstractVector{SVector{2,T}},angle::T) where {T<:Real}
    sl,cl=sincos(angle)
    return [SVector(cl*pt[1]-sl*pt[2],sl*pt[1]+cl*pt[2]) for pt in pts]
end

struct XAxisReflection{T}<:AbsReflection
    y0::T
end

struct YAxisReflection{T}<:AbsReflection
    x0::T
end

struct XYAxisReflection{T}<:AbsReflection
    x0::T
    y0::T
end

struct Rotation{T}<:AbsRotation
    order::Int
    m::Int
    c0::SVector{2,T}
end

XAxisReflection()=XAxisReflection{Float64}(0.0)
XAxisReflection(::Type{T}) where {T<:Real}=XAxisReflection{T}(zero(T))
XAxisReflection(y0::T) where {T<:Real}=XAxisReflection{T}(y0)

YAxisReflection()=YAxisReflection{Float64}(0.0)
YAxisReflection(::Type{T}) where {T<:Real}=YAxisReflection{T}(zero(T))
YAxisReflection(x0::T) where {T<:Real}=YAxisReflection{T}(x0)

XYAxisReflection()=XYAxisReflection{Float64}(0.0,0.0)
XYAxisReflection(::Type{T}) where {T<:Real}=XYAxisReflection{T}(zero(T),zero(T))
XYAxisReflection(x0::T,y0::T) where {T<:Real} = XYAxisReflection{T}(x0,y0)

Rotation(order::Int,m::Int)=Rotation{Float64}(order,m,SVector{2,Float64}(0.0,0.0))
Rotation(::Type{T},order::Int,m::Int) where {T<:Real}=Rotation{T}(order,m,SVector{2,T}(zero(T),zero(T)))
Rotation(order::Int,m::Int,c0::SVector{2,T}) where {T<:Real}=Rotation{T}(order,m,c0)

#### XAxisReflection methods ####

@inline function apply_symmetry(sym::XAxisReflection,pt::SVector{2,T}) where {T<:Real}
    return _y_coord_reflect(pt,sym.y0)
end
function apply_symmetry(sym::XAxisReflection,pts::AbstractVector{SVector{2,T}}) where {T<:Real}
    return [_y_coord_reflect(pt,sym.y0) for pt in pts]
end
@inline function apply_symmetry_normals(sym::XAxisReflection,n::SVector{2,T}) where {T<:Real}
    return _y_normal_reflect(n)
end
function apply_symmetry_normals(sym::XAxisReflection,ns::AbstractVector{SVector{2,T}}) where {T<:Real}
    return [_y_normal_reflect(n) for n in ns]
end

#### YAxisReflection methods ####

@inline function apply_symmetry(sym::YAxisReflection,pt::SVector{2,T})  where {T<:Real}
    return _x_coord_reflect(pt,sym.x0)
end
function apply_symmetry(sym::YAxisReflection,pts::AbstractVector{SVector{2,T}}) where {T<:Real}
    return [_x_coord_reflect(pt,sym.x0) for pt in pts]
end
@inline function apply_symmetry_normals(sym::YAxisReflection,n::SVector{2,T}) where {T<:Real}
    return _x_normal_reflect(n)
end
function apply_symmetry_normals(sym::YAxisReflection,ns::AbstractVector{SVector{2,T}}) where {T<:Real}
    return [_x_normal_reflect(n) for n in ns]
end

#### XYAxisReflection methods ####

@inline function apply_symmetry(sym::XYAxisReflection,pt::SVector{2,T})  where {T<:Real}
    return _xy_coord_reflect(pt,sym.x0,sym.y0)
end
function apply_symmetry(sym::XYAxisReflection,pts::AbstractVector{SVector{2,T}}) where {T<:Real}
    return [_xy_coord_reflect(pt,sym.x0,sym.y0) for pt in pts]
end
@inline function apply_symmetry_normals(sym::XYAxisReflection,n::SVector{2,T}) where {T<:Real}
    return _xy_normal_reflect(n)
end
function apply_symmetry_normals(sym::XYAxisReflection,ns::AbstractVector{SVector{2,T}}) where {T<:Real}
    return [_xy_normal_reflect(n) for n in ns]
end

#### Rotation methods ####

@inline function apply_symmetry(sym::Rotation,pt::SVector{2,T}) where {T<:Real}
    return _rotate(pt,T(TWO_PI*sym.m/sym.order),sym.c0)
end
function apply_symmetry(sym::Rotation,pts::AbstractVector{SVector{2,T}}) where {T<:Real}
    angle=T(TWO_PI*sym.m/sym.order)
    return _rotate(pts,angle,sym.c0)
end
@inline function apply_symmetry_normals(sym::Rotation,n::SVector{2,T}) where {T<:Real}
    return _rotate_normal(n,T(TWO_PI*sym.m/sym.order))
end
function apply_symmetry_normals(sym::Rotation,ns::AbstractVector{SVector{2,T}}) where {T<:Real}
    angle=T(TWO_PI*sym.m/sym.order)
    return _rotate_normal(ns,angle)
end

#### Phase space methods for Reflections ###

function apply_symmetry_pb(sym::AbsReflection,sym_sector::Int,s::T,p::T,L::T) where {T<:Real}
    if sym_sector==1
        return s,p
    elseif sym_sector==2
        return 2*L-s,-p
    elseif sym_sector==3
        return 2*L+s,p
    elseif sym_sector==4    
        return 4*L-s,-p
    else
        error("Invalid symmetry sector: $sym_sector")
    end
end

#### Composite symmetry generators ####

function D2_symmetry(x0::T,y0::T) where {T<:Real}
    return [YAxisReflection(x0),XYAxisReflection(x0,y0),XAxisReflection(y0)]
end
function D2_symmetry()
    return D2_symmetry(0.0,0.0)
end
function Cn_symmetry(n::Int,c0::SVector{2,T}) where {T<:Real}
    return [Rotation(n,i,c0) for i in 0:(n-1)]
end
function Cn_symmetry(n::Int)
    return Cn_symmetry(n,SVector(0.0,0.0))
end