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
using Roots,Optim
#using DataInterpolations

abstract type AbsBoundaryCondition end
abstract type AbsCurve{BC} end 
abstract type AbsPolarCurve{BC}<:AbsCurve{BC} end 
abstract type AbsCompositeCurve end
abstract type AbsSymmetry end
abstract type AbsCoords end
abstract type AbsDomain end
abstract type AbsSimpleDomain<:AbsDomain end
abstract type AbsPolarDomain<:AbsSimpleDomain end
abstract type AbsCompositeDomain<:AbsDomain end
abstract type AbsBilliard end
abstract type AbsSampler end
abstract type AbsReflection<:AbsSymmetry end
abstract type AbsRotation<:AbsSymmetry end

export AbsCurve,AbsPolarCurve,AbsCompositeCurve,AbsCoords,AbsDomain,AbsPolarDomain,AbsCompositeDomain,AbsSimpleDomain,AbsBilliard,AbsBoundaryCondition,AbsSampler,AbsSymmetry

include("geometry/utils.jl")
export is_overlaping,connect_curves,is_connected,is_closed,angle,signed_angle,find_unique_elements,point_curve_parameter
include("geometry/geometry.jl")
export curvature,is_inside,domain_gradient_vector
include("interpolations/cheb_interp.jl")
include("geometry/boundarytypes.jl")
export Transparent,PeriodicX,ReflectionSymmetry,QuantumSolverIgnore,get_boundary_curves,get_all_domains,get_domain,get_all_curves,get_curve,update_boundary_condition,SpecularReflection

include("geometry/segments/linesegment.jl")
export LineSegment,curve,domain_fun,arc_length,tangent,arc_length
include("geometry/segments/circlesegment.jl")
export CircleSegment
include("geometry/segments/compositecurves.jl")
export CompositeCurve
include("geometry/segments/polarsegment.jl")
export PolarSegment,polar_domain,polar_radius
include("geometry/domains/circular.jl")
export CircleCap,CircleWedge
include("geometry/domains/compositedomains.jl")
export SimpleDomain,CompositeDomain,reset_ids!
include("geometry/domains/polygons.jl")
export Polygon
include("geometry/billiards/limacon.jl")
export LimaconSegment
include("geometry/billiards/mushroom.jl")
export MushroomBilliard
include("geometry/billiards/stadium.jl")
export StadiumBilliard
include("geometry/billiards/triangle.jl")
export TriangleBilliard

include("geometry/arclength.jl")
export _arc_length_integrand,construct_arc_length_interpolation
include("geometry/inversions.jl")
export invert_curve,invert_curve_roots,invert_curve_optim,invert_curve_grid_search,invert_curve_nelder_mead
include("geometry/poincarebirkhoff.jl")
export fundamental_pb_coords,get_pb_curve,pb_coords,pb_sectors
include("quadrature/samplers.jl")
export LinearNodes,GaussLegendreNodes,FourierNodes,sample_points
end
