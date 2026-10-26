"""
    BilliardGeometry

Planar billiard geometry: boundary curves, domains, billiards, discrete
symmetries and boundary quadrature nodes.

## Description
A [`AbsBilliard`](@ref) stores a `fundamental_domain` (an [`AbsDomain`](@ref)
made of [`AbsCurve`](@ref) segments carrying an [`AbsBoundaryCondition`](@ref))
and a [`SymmetryRegistry`](@ref) of the non-identity elements of its discrete
symmetry group. The complete physical boundary is reconstructed from the two by
[`full_boundary`](@ref). All curves are parametrized by `t ∈ [0,1]`, and the
sign of [`domain_fun`](@ref) tells on which side of a curve a point lies.

## Exports
* Abstract interfaces: [`AbsCurve`](@ref), [`AbsPolarCurve`](@ref), [`AbsCompositeCurve`](@ref), [`AbsDomain`](@ref), [`AbsSimpleDomain`](@ref), [`AbsPolarDomain`](@ref), [`AbsCompositeDomain`](@ref), [`AbsMultiplyConnectedDomain`](@ref), [`AbsDomainHole`](@ref), [`AbsBilliard`](@ref), [`AbsBoundaryCondition`](@ref), [`AbsSymmetry`](@ref), [`AbsCoords`](@ref), [`AbsSampler`](@ref).
* Boundary conditions: [`SpecularReflection`](@ref), [`Transparent`](@ref), [`PeriodicX`](@ref), [`SymmetryWall`](@ref), [`QuantumSolverIgnore`](@ref).
* Curves: [`LineSegment`](@ref), [`CircleSegment`](@ref), [`FourierCoeffPolarSegment`](@ref), [`PolarSegment`](@ref), [`CompositeCurve`](@ref).
* Domains: [`Polygon`](@ref), [`CircleWedge`](@ref), [`SimpleDomain`](@ref), [`CompositeDomain`](@ref), [`PolarDomain`](@ref), [`MultiplyConnectedDomain`](@ref), [`DomainHole`](@ref).
* Billiards: [`TriangleBilliard`](@ref), [`StadiumBilliard`](@ref), [`MushroomBilliard`](@ref), [`PolarBilliard`](@ref), [`LimaconBilliard`](@ref), [`CircleBilliard`](@ref), [`EllipseBilliard`](@ref), [`RectangleBilliard`](@ref), [`PolygonBilliard`](@ref), [`StarBilliard`](@ref), [`C3Billiard`](@ref), [`ProsenBilliard`](@ref), [`AnnularBilliard`](@ref), [`SinaiBilliard`](@ref), [`SquareWithinSquareBilliard`](@ref).
* Symmetries: [`XAxisReflection`](@ref), [`YAxisReflection`](@ref), [`XYAxisReflection`](@ref), [`DiagonalReflection`](@ref), [`AntiDiagonalReflection`](@ref), [`CompositeReflection`](@ref), [`NFoldRotation`](@ref), [`SymmetryRegistry`](@ref), [`SymmetryOrbitMap`](@ref).
* Boundary samplers and Kress grading: [`LinearNodes`](@ref), [`GaussLegendreNodes`](@ref), [`FourierNodes`](@ref), [`sample_points`](@ref), [`kress_graded_nodes_data`](@ref), [`multi_kress_graded_nodes_data`](@ref).
"""
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

"""
    AbsBoundaryCondition

Supertype of the boundary conditions a curve carries in its `bc` field.

## Description
The condition decides how a curve is treated by solvers and by the geometry
routines: [`SpecularReflection`](@ref) curves are physical walls,
[`SymmetryWall`](@ref) curves are cuts closing a symmetry-reduced fundamental
domain, [`Transparent`](@ref) curves are internal joins between subdomains of an
[`AbsCompositeDomain`](@ref), [`PeriodicX`](@ref) curves carry a periodic
condition in `x` and [`QuantumSolverIgnore`](@ref) curves are virtual edges
excluded from the solver boundary.
"""
abstract type AbsBoundaryCondition end

"""
    AbsCurve{BC<:AbsBoundaryCondition}

Supertype of parametrized planar boundary curves `r(t)`, `t ∈ [0,1]`.

## Description
Concrete subtypes carry the fields `length` (arc length), `orientation`
(`±1`, selecting which side of the curve is the inside), `bc::BC` (boundary
condition), `domain_id` and `segment_id` (identify the curve within a billiard,
see [`get_curve`](@ref)).

## API
Every subtype must implement `curve(crv, t)` returning an `SVector{2}` and
[`domain_fun`](@ref), which is negative inside the region bounded by the curve.
The following have generic fallbacks based on automatic differentiation and
quadrature and may be specialized:
- [`tangent`](@ref), [`tangent_2`](@ref), [`tangent_vec`](@ref), [`normal_vec`](@ref), [`curvature`](@ref)
- [`arc_length`](@ref), [`invert_curve`](@ref), [`area`](@ref), [`is_inside`](@ref)
"""
abstract type AbsCurve{BC} end 

"""
    AbsPolarCurve{BC} <: AbsCurve{BC}

Supertype of curves given in polar form `r(φ)` about a `center`, with
`φ = shift_angle + t*arc_angle`. Subtypes implement [`polar_radius`](@ref); the
curve, [`domain_fun`](@ref) and arc length then follow generically. See
[`FourierCoeffPolarSegment`](@ref) and [`PolarSegment`](@ref).
"""
abstract type AbsPolarCurve{BC} <: AbsCurve{BC} end 

"""
    AbsCompositeCurve

Supertype of curves assembled from consecutive [`AbsCurve`](@ref) pieces, stored
in the field `subcurves`, with the cumulative lengths `end_lengths` and the total
`length`. See [`CompositeCurve`](@ref).
"""
abstract type AbsCompositeCurve end

"""
    AbsSymmetry

Supertype of the discrete geometric symmetries of a billiard.

## Description
A concrete symmetry carries a `sym_id::Int` assigned by
[`register_symmetries`](@ref). It is purely geometric: no representation or
parity data is stored. Subtypes implement [`apply_symmetry`](@ref),
[`orientation_reversing`](@ref) and [`symmetry_node_multiple`](@ref).
"""
abstract type AbsSymmetry end

"""
    AbsCoords

Supertype of coordinates on the billiard boundary, e.g. [`PoincareBirkhoff`](@ref).
"""
abstract type AbsCoords end

"""
    AbsDomain

Supertype of planar regions bounded by [`AbsCurve`](@ref) segments.

## API
The following functions are defined for every domain:
- [`get_boundary_curves`](@ref), [`get_all_curves`](@ref), [`get_all_domains`](@ref)
- [`boundary_components`](@ref), [`geometric_boundary_components`](@ref)
- [`is_inside`](@ref), [`area`](@ref), [`genus`](@ref)
"""
abstract type AbsDomain end

"""
    AbsSimpleDomain <: AbsDomain

Supertype of domains bounded by a single list of curves, stored in the fields
`boundary` (curves), `corners` and `id`. See [`SimpleDomain`](@ref),
[`Polygon`](@ref), [`CircleWedge`](@ref).
"""
abstract type AbsSimpleDomain <: AbsDomain end

"""
    AbsPolarDomain <: AbsSimpleDomain

Supertype of simple domains bounded by polar curves. See [`PolarDomain`](@ref).
"""
abstract type AbsPolarDomain <: AbsSimpleDomain end

"""
    AbsCompositeDomain <: AbsDomain

Supertype of domains that are the union of several [`AbsSimpleDomain`](@ref)
stored in the field `subdomains`. Subdomains are joined along curves with a
[`Transparent`](@ref) boundary condition. See [`CompositeDomain`](@ref).
"""
abstract type AbsCompositeDomain <: AbsDomain end

"""
    AbsMultiplyConnectedDomain <: AbsDomain

Supertype of domains with holes: the fields `subdomains` describe the outer
region and `holes` the removed regions. See [`MultiplyConnectedDomain`](@ref).
"""
abstract type AbsMultiplyConnectedDomain <: AbsDomain end

"""
    AbsDomainHole

Supertype of regions removed from an [`AbsMultiplyConnectedDomain`](@ref). A hole
lives in the full physical frame and is never reconstructed by the symmetry
registry. See [`DomainHole`](@ref).
"""
abstract type AbsDomainHole end

"""
    AbsBilliard

Supertype of billiards.

## Description
Concrete subtypes carry the fields `fundamental_domain::AbsDomain` and
`symmetries::SymmetryRegistry`. The complete physical boundary is obtained with
[`full_boundary`](@ref); the physical area with [`area`](@ref).

## API
The following functions can be evaluated for every billiard:
- [`get_boundary_curves`](@ref), [`get_all_curves`](@ref), [`get_curve`](@ref), [`update_boundary_condition`](@ref)
- [`full_boundary`](@ref), [`full_boundary_components`](@ref), [`boundary_components`](@ref)
- [`area`](@ref), [`fundamental_area`](@ref), [`corner_angles`](@ref), [`is_inside`](@ref)
- [`pb_coords`](@ref), [`pb_sectors`](@ref), [`get_symmetries`](@ref)
"""
abstract type AbsBilliard end

"""
    AbsSampler

Supertype of boundary parameter samplers. Subtypes implement
[`sample_points`](@ref), returning nodes in `t ∈ [0,1]` and quadrature weights.
"""
abstract type AbsSampler end

"""
    AbsReflection <: AbsSymmetry

Supertype of the reflection-type symmetries: [`XAxisReflection`](@ref),
[`YAxisReflection`](@ref), [`XYAxisReflection`](@ref),
[`DiagonalReflection`](@ref), [`AntiDiagonalReflection`](@ref) and the generator
set [`CompositeReflection`](@ref).
"""
abstract type AbsReflection <: AbsSymmetry end

export AbsCurve, AbsPolarCurve, AbsCompositeCurve, AbsCoords, AbsDomain, AbsPolarDomain, AbsCompositeDomain, AbsSimpleDomain, AbsMultiplyConnectedDomain, AbsDomainHole, AbsBilliard, AbsBoundaryCondition, AbsSampler, AbsSymmetry

include("geometry/utils.jl")
export is_overlaping, is_connected, is_closed, connect_curves, find_unique_elements, point_curve_parameter

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
include("geometry/domains/holes.jl")
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
export XAxisReflection, YAxisReflection, XYAxisReflection, DiagonalReflection, AntiDiagonalReflection, CompositeReflection, NFoldRotation, apply_symmetry, apply_symmetry_pb, D2_symmetry, Cn_symmetry, orientation_reversing
export SymmetryRegistry, register_symmetries, symmetry_of
export SymmetryOrbitMap, fundamental_size, full_size, orbit_size, symmetry_orbit, symmetry_node_multiple, symmetry_index_orbits, get_symmetries
export SpecularReflection, QuantumSolverIgnore, Transparent, PeriodicX, SymmetryWall, get_boundary_curves, get_all_curves, get_curve, get_all_domains, get_domain, update_boundary_condition, genus, boundary_components, geometric_boundary_components
export LineSegment
export PolarSegment, FourierCoeffPolarSegment, polar_radius
export CircleSegment
export CompositeCurve
export Polygon
export CircleWedge
export SimpleDomain, CompositeDomain, reset_ids!
export MultiplyConnectedDomain, DomainHole, hole_domains, hole_loop_count, hole_leaves, hole_domain_ids, is_hole_curve, component_of_domain_id, symmetry_compatible, validate_holes
export boundary_component_is_hole, get_outer_boundary_curves, get_hole_boundary_curves, BoundaryComponentMap, full_boundary_components
export TriangleBilliard
export StadiumBilliard
export MushroomBilliard
export PolarBilliard, PolarDomain
export LimaconBilliard
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

include("geometry/billiards/square_in_square.jl")
export SquareWithinSquareBilliard

end
