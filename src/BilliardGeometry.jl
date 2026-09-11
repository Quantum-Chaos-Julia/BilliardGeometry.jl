module BilliardGeometry
using StaticArrays 
using LinearAlgebra
using CoordinateTransformations, Rotations
using ForwardDiff
using CircularArrays
using Accessors
using FastGaussQuadrature, QuadGK
using StatsBase
using Elliptic
using Roots, Optim
using DataInterpolations
using FFTW

abstract type AbsBoundaryCondition end
abstract type AbsCurve{BC} end 
abstract type AbsPolarCurve{BC} <: AbsCurve{BC} end 
abstract type AbsCompositeCurve end
abstract type AbsSymmetry end
abstract type AbsCoords end
abstract type AbsDomain end
abstract type AbsSimpleDomain <: AbsDomain end
abstract type AbsPolarDomain <: AbsSimpleDomain end
abstract type AbsCompositeDomain <: AbsDomain end
abstract type AbsMultiplyConnectedDomain <: AbsDomain end
abstract type AbsBilliard end

abstract type AbsSampler end
abstract type AbsReflection <: AbsSymmetry end

export AbsCurve, AbsPolarCurve, AbsCompositeCurve, AbsCoords, AbsDomain, AbsPolarDomain, AbsCompositeDomain, AbsSimpleDomain, AbsMultiplyConnectedDomain, AbsBilliard, AbsBoundaryCondition, AbsSampler, AbsSymmetry

include("geometry/utils.jl")
include("geometry/geometry.jl")


include("quadrature/samplers.jl")
export LinearNodes, GaussLegendreNodes, FourierNodes, sample_points
include("quadrature/kressgrading.jl")
export kress_R!, kress_R_even!, kress_R_odd!, s_mid, kress_graded_nodes_data, multi_kress_graded_nodes_data
end
