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
export is_overlaping, is_connected, is_closed, angle, connect_curves, find_unique_elements, point_curve_parameter

include("geometry/symmetry.jl")
include("geometry/symmetryregistry.jl")
include("geometry/symmetryorbits.jl")
include("geometry/boundarytypes.jl")
include("geometry/segments/linesegment.jl")
include("geometry/segments/polarcurves.jl")
include("geometry/segments/circlesegment.jl")
include("geometry/segments/compositecurves.jl")
include("geometry/domains/polygons.jl")
include("geometry/domains/circular.jl")
include("geometry/domains/compositedomains.jl")
include("geometry/domains/multiplyconnecteddomains.jl")
include("geometry/billiards/triangle.jl")
include("geometry/billiards/stadium.jl")
include("geometry/billiards/mushroom.jl")
include("geometry/billiards/polar.jl")
include("geometry/billiards/limacon.jl")
include("geometry/billiards/circle.jl")
include("geometry/billiards/ellipse.jl")
include("geometry/billiards/rectangle.jl")
include("geometry/billiards/polygon.jl")
include("geometry/billiards/star.jl")
include("geometry/billiards/c3.jl")
include("geometry/billiards/prosen.jl")
include("geometry/billiards/annular.jl")
include("geometry/fullboundary.jl")
include("geometry/inversions.jl")
include("geometry/arclength.jl")
include("geometry/curvederivatives.jl")
include("geometry/boundarycomponents.jl")
include("geometry/poincarebirkhoff.jl")
include("geometry/geometry.jl")
export XAxisReflection, YAxisReflection, XYAxisReflection, DiagonalReflection, AntiDiagonalReflection, CompositeReflection, NFoldRotation, apply_symmetry, apply_symmetry_pb, D2_symmetry, Cn_symmetry
export SymmetryRegistry, register_symmetries, symmetry_of
export SymmetryOrbitMap, fundamental_size, full_size, orbit_size, symmetry_orbit, symmetry_node_multiple, symmetry_index_orbits
export SpecularReflection, QuantumSolverIgnore, Transparent, PeriodicX, SymmetryWall, get_boundary_curves, get_all_curves, get_curve, get_all_domains, get_domain, update_boundary_condition, genus, boundary_components
export LineSegment
export PolarSegment, FourierCoeffPolarSegment, polar_radius
export CircleSegment
export CompositeCurve
export Polygon
export CircleWedge
export SimpleDomain, CompositeDomain, reset_ids!
export MultiplyConnectedDomain
export TriangleBilliard
export StadiumBilliard
export MushroomBilliard
export PolarBilliard, PolarDomain
export LimaconBilliard, LimaconSegment
export CircleBilliard
export EllipseBilliard
export RectangleBilliard
export PolygonBilliard
export StarBilliard
export C3Billiard
export ProsenBilliard
export AnnularBilliard
export full_boundary
export invert_curve
export arc_length, construct_arc_length_interpolation
export tangent, tangent_2, tangent_vec, normal_vec, curvature
export component_lengths, print_component_junctions
export PoincareBirkhoff, pb_coords, get_pb_curve, pb_sectors
export is_inside, curve, domain_fun, domain_gradient_vector

include("geometry/area.jl")
export area, fundamental_area, corner_angles

include("quadrature/samplers.jl")
export LinearNodes, GaussLegendreNodes, FourierNodes, sample_points
include("quadrature/kressgrading.jl")
export kress_R!, kress_R_even!, kress_R_odd!, s_mid, kress_graded_nodes_data, multi_kress_graded_nodes_data

include("geometry/billiards/sinai.jl")
export SinaiBilliard



end
