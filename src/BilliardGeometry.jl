module BilliardGeometry
using StaticArrays 
using LinearAlgebra
using CoordinateTransformations, Rotations
using ForwardDiff
using CircularArrays
using Accessors

abstract type AbsBoundaryCondition end
abstract type AbsCurve{BC} end 
abstract type AbsSymmetry end
abstract type AbsDomain end
abstract type AbsSimpleDomain <: AbsDomain end
abstract type AbsComplexDomain <: AbsDomain end
abstract type AbsBilliard end


abstract type AbsReflection <: AbsSymmetry end

export AbsCurve, AbsDomain, AbsComplexDomain, AbsSimpleDomain, AbsBilliard, AbsBoundaryCondition

include("geometry/utils.jl")
include("geometry/geometry.jl")
include("geometry/boundarytypes.jl")


end
