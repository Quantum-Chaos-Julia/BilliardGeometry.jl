module BilliardGeometry
using StaticArrays 
using LinearAlgebra
using CoordinateTransformations, Rotations
using ForwardDiff
using CircularArrays
using Accessors
using FastGaussQuadrature
using StatsBase
using Elliptic

abstract type AbsBoundaryCondition end
abstract type AbsCurve{BC} end 
abstract type AbsPolarCurve{BC} <: AbsCurve{BC} end 
abstract type AbsSymmetry end
abstract type AbsDomain end
abstract type AbsSimpleDomain <: AbsDomain end
abstract type AbsPolarDomain <: AbsDomain end
abstract type AbsCompositeDomain <: AbsDomain end
abstract type AbsBilliard end

abstract type AbsSampler end
abstract type AbsReflection <: AbsSymmetry end

export AbsCurve,AbsPolarCurve, AbsDomain, AbsPolarDomain, AbsCompositeDomain, AbsSimpleDomain, AbsBilliard, AbsBoundaryCondition, AbsSampler

include("geometry/utils.jl")
include("geometry/geometry.jl")


include("quadrature/samplers.jl")
export LinearNodes, GaussLegendreNodes, FourierNodes, sample_points
end
