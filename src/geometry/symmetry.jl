using CoordinateTransformations
using Rotations
using StaticArrays

#struct Ident{T} <: AbsSymmetry end #Identity transformation


ident = IdentityTransformation()
reflect_x = LinearMap(SMatrix{2,2}([-1.0 0.0;0.0 1.0]))
reflect_y = LinearMap(SMatrix{2,2}([1.0 0.0;0.0 -1.0]))
reflect_xy = reflect_x ∘ reflect_y

abstract type AbsReflection <: AbsSymmetry end



struct XReflection{T} <: AbsReflection where T<:Real
    N_sectors::Int64
    phase::T #phase factor is only used in quantum billiards
end

XReflection(N_sectors) = XReflection(N_sectors, 1.0)

struct YReflection{T} <: AbsReflection where T<:Real
    N_sectors::Int64
    phase::T
end

YReflection(N_sectors) = YReflection(N_sectors, 1.0)

struct XYReflection{T} <: AbsReflection where T<:Real
    N_sectors::Int64
    phase::T
end

XYReflection(N_sectors) = XYReflection(N_sectors,1.0)

#reflections
#sector change when crossing y axis means x coordinate is reflected
# N_sectors is number of symmetry sectors 2 - for only x symmetry, 4 - for x and y symmetry 
function new_sector(sym_x::R, id) where  R <: XReflection
    N_sectors = sym_x.N_sectors
    
    if N_sectors == 2
        if id == 1
            new_id = 2
        else
            new_id = 1
        end 
    else # N_sectors = 4 case
        if id == 1
            new_id = 2
        end
        if id == 2
            new_id = 1
        end
        if id == 3
            new_id = 4
        end
        if id == 4
            new_id = 3
        end
    end
    return new_id
end

#sector change when crossing x axis means y coordinate is reflected
function new_sector(sym_y::R, id) where  R <: YReflection
    N_sectors = sym_y.N_sectors
    
    if N_sectors == 2
        if id == 1
            new_id = 2
        else
            new_id = 1
        end 
    else # N_sectors = 4 case
        if id == 1
            new_id = 4
        end
        if id == 2
            new_id = 3
        end
        if id == 3
            new_id = 2
        end
        if id == 4
            new_id = 1
        end
    end
    return new_id
end


abstract type AbsRotation <: AbsSymmetry end

#=
function symmetry_action(billiard::B, pt) where  B <: AbsBilliard
    let id = pt.sym_sector
        return billiard.symmetry[id](pt)
    end
end
=#

#Reflection symetries