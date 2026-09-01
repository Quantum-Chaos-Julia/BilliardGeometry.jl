module BilliardGeometry

using StaticArrays
using LinearAlgebra
using CoordinateTransformations,Rotations
using ForwardDiff
using CircularArrays
using Accessors
using FastGaussQuadrature,QuadGK
using StatsBase
using Elliptic
using Roots,Optim
using DataInterpolations

abstract type AbsBoundaryCondition end
abstract type AbsCurve{BC} end
abstract type AbsPolarCurve{BC}<:AbsCurve{BC} end
abstract type AbsCompositeCurve end
abstract type AbsSymmetry end
abstract type AbsReflection<:AbsSymmetry end
abstract type AbsRotation<:AbsSymmetry end
abstract type AbsCoords end
abstract type AbsDomain end
abstract type AbsSimpleDomain<:AbsDomain end
abstract type AbsPolarDomain<:AbsSimpleDomain end
abstract type AbsCompositeDomain<:AbsDomain end
abstract type AbsBilliard end
abstract type AbsSampler end

include("geometry/symmetry.jl")
include("geometry/boundarytypes.jl")
include("geometry/segments/linesegment.jl")
include("geometry/segments/polarcurves.jl")
include("geometry/segments/circlesegment.jl")
include("geometry/utils.jl")
include("geometry/segments/compositecurves.jl")
include("geometry/domains/polygons.jl")
include("geometry/domains/circular.jl")
include("geometry/domains/compositedomains.jl")
include("geometry/segments/curve_derivatives.jl")
include("geometry/billiards/stadium.jl")
include("geometry/billiards/polar.jl")
include("geometry/geometry.jl")
include("geometry/inversions.jl")
include("geometry/arclength.jl")
include("geometry/poincarebirkhoff.jl")
include("quadrature/samplers.jl")

export AbsBoundaryCondition
export AbsCurve,AbsPolarCurve,AbsCompositeCurve
export AbsSymmetry,AbsReflection,AbsRotation
export AbsCoords
export AbsDomain,AbsSimpleDomain,AbsPolarDomain,AbsCompositeDomain
export AbsBilliard,AbsSampler

export XAxisReflection,YAxisReflection,XYAxisReflection
export DiagonalReflection,AntiDiagonalReflection
export NFoldRotation
export apply_symmetry,apply_symmetry_pb
export D2_symmetry,Cn_symmetry,symmetry_irrep_character

export SpecularReflection,QuantumSolverIgnore
export Transparent,PeriodicX,ReflectionSymmetry
export get_boundary_curves,get_all_curves,get_curve
export get_all_domains,get_domain,update_boundary_condition

export LineSegment
export PolarSegment,polar_radius
export CircleSegment
export CompositeCurve
export curve,domain_fun
export tangent,tangent_2,curvature

export is_overlaping,is_connected,is_closed
export angle,signed_angle
export connect_curves,find_unique_elements,point_curve_parameter

export Polygon
export CircleCap,CircleWedge
export SimpleDomain,CompositeDomain,PolarDomain
export reset_ids!

export TriangleBilliard
export StadiumBilliard
export MushroomBilliard
export PolarBilliard
export LimaconBilliard,LimaconSegment

export is_inside,curve,domain_fun,domain_gradient_vector
export invert_curve
export arc_length,construct_arc_length_interpolation

export PoincareBirkhoff
export pb_coords,get_pb_curve,pb_sectors

export LinearNodes,GaussLegendreNodes,FourierNodes
export sample_points,random_interior_points

end