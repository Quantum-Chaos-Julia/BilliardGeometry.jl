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
#using DataInterpolations

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
abstract type AbsBilliard end

abstract type AbsSampler end
abstract type AbsReflection <: AbsSymmetry end
abstract type AbsRotation <: AbsSymmetry end

export AbsCurve, AbsPolarCurve, AbsCompositeCurve, AbsCoords, AbsDomain, AbsPolarDomain, AbsCompositeDomain, AbsSimpleDomain, AbsBilliard, AbsBoundaryCondition, AbsSampler, AbsSymmetry

include("geometry/utils.jl")
include("geometry/geometry.jl")

include("geometry/segments/linesegment.jl")
include("geometry/segments/circlesegment.jl")
include("geometry/segments/compositecurves.jl")
include("geometry/segments/polarsegment.jl")

include("interpolations/cheb_interp.jl")

include("geometry/domains/circular.jl")
include("geometry/domains/compositedomains.jl")
include("geometry/domains/polygons.jl")

include("geometry/arclength.jl")
include("geometry/boundarytypes.jl")
include("geometry/inversions.jl")
include("geometry/poincarebirkhoff.jl")

include("quadrature/samplers.jl")
export LinearNodes, GaussLegendreNodes, FourierNodes, sample_points
end
