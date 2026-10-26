# BilliardGeometry

[![Build Status](https://github.com/clozej/BilliardGeometry.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/clozej/BilliardGeometry.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Julia](https://img.shields.io/badge/julia-%E2%89%A5%201.9-9558B2.svg?logo=julia)](https://julialang.org)
[![License: MIT](https://img.shields.io/badge/license-MIT-yellow.svg)](LICENSE)

BilliardGeometry.jl is a Julia library for describing the geometry of planar billiard tables. It is the shared geometric foundation of the billiard ecosystem: [QuantumBilliards.jl](https://github.com/Quantum-Chaos-Julia/QuantumBilliards.jl) uses it to build boundary discretizations for quantum eigenvalue problems, and `ClassicalBilliards.jl` (under development) uses it for classical ray dynamics and Poincaré-Birkhoff maps. The package is purely geometric: it knows nothing about wave equations or trajectories, only about curves, domains, boundaries, symmetries and the quadratures needed to sample them.

## Introduction
A billiard is a bounded planar domain $\mathcal{B} \subset \mathbb{R}^2$ with boundary $\partial \mathcal{B}$. Both classical (a point particle reflected specularly from the boundary) and quantum (the Helmholtz equation $(\nabla^2 + k^2)\psi = 0$ with a boundary condition on $\partial \mathcal{B}$) problems are defined entirely by this boundary. The package provides one consistent description of $\partial \mathcal{B}$ so that the same billiard object can be passed to a classical or a quantum solver without modification.

A billiard is assembled from three layers:

1. **Curves** - parametrized boundary pieces $\gamma(t),\ t \in [0,1]$ with analytic tangents, normals, curvature and arc length.
2. **Domains** - ordered curves forming a closed region, optionally composed of several sub-domains glued along interfaces, or multiply connected with explicit holes.
3. **Billiards** - a fundamental domain plus a registry of discrete symmetries. The full table is reconstructed from the fundamental domain by applying the symmetries.

Each boundary curve carries a **boundary role** (specular reflection, transparent interface, symmetry wall, ...) so that downstream packages can decide how to treat it.

# Package features
- Parametrized boundary curves: `LineSegment`, `CircleSegment`, `PolarSegment`, `FourierCoeffPolarSegment`, and `CompositeCurve`, with tangent, normal, curvature, arc length, area and inversion routines.
- Ready-made billiards: circle, ellipse, rectangle, stadium, Robnik limaçon, mushroom, Prosen, Sinai, annular, square-in-square, triangles, general polygons, polar/Fourier-ripple, star and `C3` billiards.
- Simple, composite and multiply connected domains. A `MultiplyConnectedDomain` is an outer region minus explicit `DomainHole`s (semantics $\Omega = \bigcup \text{subdomains} \setminus \bigcup \text{holes}$), with `genus`, boundary components and separate access to outer and hole boundaries.
- Boundary roles as types: `SpecularReflection`, `Transparent`, `SymmetryWall`, `PeriodicX`, `QuantumSolverIgnore`, so that interfaces between sub-domains and symmetry cuts are explicit.
- Discrete symmetries: axis and diagonal reflections, `D2`, `Cn` rotations, and composite reflections, collected in a `SymmetryRegistry` whose generators carry a stable `sym_id` (lookup with `symmetry_of`). Exact, floating-point independent symmetry orbit maps (`SymmetryOrbitMap`) permute boundary points between fundamental and full boundary.
- `full_boundary` reconstructs the complete physical boundary from the fundamental one using the registered symmetries, with consistently oriented traversal. `full_boundary_components` returns the closed loops (outer boundary and holes) together with a `BoundaryComponentMap`.
- Boundary samplers and quadratures: `LinearNodes`, `GaussLegendreNodes`, `FourierNodes`, and Kress graded nodes for boundary corner singularities (`kress_graded_nodes_data`, `multi_kress_graded_nodes_data`).
- Geometric quantities: `area`, `fundamental_area`, `corner_angles`, `is_inside`, `component_lengths`.
- Poincaré-Birkhoff phase space coordinates (`PoincareBirkhoff`, `pb_coords`) for the boundary map, shared by classical dynamics and by quantum Husimi functions.
- Generic in the number type `T<:Real`, built on `StaticArrays` for allocation-free point evaluation.
- Type-hierarchy based, multiple-dispatch design: new curves, domains and billiards plug into the existing infrastructure by implementing only the specialized methods.

# Type hierarchy

| Abstract type | Role | Concrete examples |
|---|---|---|
| `AbsCurve{BC}` | Parametrized boundary piece with boundary role `BC` | `LineSegment`, `CircleSegment` |
| `AbsPolarCurve{BC}` | Curve given as $r(\varphi)$ | `PolarSegment`, `FourierCoeffPolarSegment` |
| `AbsCompositeCurve` | Concatenation of curves | `CompositeCurve` |
| `AbsDomain` | Region bounded by curves | `SimpleDomain`, `Polygon`, `CircleWedge`, `PolarDomain` (`AbsSimpleDomain`); `CompositeDomain` (`AbsCompositeDomain`); `MultiplyConnectedDomain` (`AbsMultiplyConnectedDomain`) |
| `AbsDomainHole` | Region removed from a multiply connected domain | `DomainHole` |
| `AbsBilliard` | Fundamental domain plus symmetries | `StadiumBilliard`, `ProsenBilliard`, `PolygonBilliard`, `AnnularBilliard`, ... |
| `AbsSymmetry` | Discrete symmetry operation | `XAxisReflection`, `YAxisReflection`, `XYAxisReflection`, `DiagonalReflection`, `AntiDiagonalReflection`, `CompositeReflection`, `NFoldRotation` |
| `AbsBoundaryCondition` | Role of a boundary curve | `SpecularReflection`, `Transparent`, `SymmetryWall` |
| `AbsSampler` | Boundary parameter nodes and weights | `LinearNodes`, `GaussLegendreNodes`, `FourierNodes` |
| `AbsCoords` | Boundary phase space coordinates | `PoincareBirkhoff` |

# Quick start

```julia
using BilliardGeometry, StaticArrays

billiard = StadiumBilliard(1.0)          # half-width of the straight part

area(billiard)                            # area of the full table
fundamental_area(billiard)                # area of the symmetry-reduced domain
corner_angles(billiard)                   # corner angles of the fundamental domain

is_inside(billiard, SVector(0.3, 0.2))    # point-in-billiard test

crvs = get_boundary_curves(billiard)      # physical boundary of the fundamental domain
full = full_boundary(billiard)            # reconstructed full boundary

# sample the first curve with Gauss-Legendre nodes
t, dt = sample_points(GaussLegendreNodes(), 64)
pts   = curve(crvs[1], t)
nvec  = normal_vec(crvs[1], t)
kappa = curvature(crvs[1], t)
```

Visualization is provided by [BilliardGeometryPlotting.jl](https://github.com/Quantum-Chaos-Julia/BilliardGeometryPlotting.jl) (Makie based).

# Billiards
Available out of the box:

| Billiard | Constructor |
|---|---|
| Circle | `CircleBilliard(R)` |
| Ellipse | `EllipseBilliard(...)` |
| Rectangle | `RectangleBilliard(...)` |
| Stadium | `StadiumBilliard(half_width)` |
| Robnik limaçon | `LimaconBilliard(...)` |
| Mushroom | `MushroomBilliard(half_width, stem_height; R)` |
| Prosen | `ProsenBilliard(...)` |
| Sinai | `SinaiBilliard(a, R)` |
| Annular | `AnnularBilliard(R_inner, shift; R_outer)` |
| Square within square | `SquareWithinSquareBilliard(...)` |
| Triangle | `TriangleBilliard(gamma, chi)` |
| Polygon | `PolygonBilliard(vertices)` |
| Star | `StarBilliard(R, a, n)` |
| `C3` | `C3Billiard(a)` |
| Polar / Fourier ripple | `PolarBilliard(coef)` |

Custom billiards are assembled from curves, grouped into a `SimpleDomain`, `CompositeDomain` or `MultiplyConnectedDomain`, and given a `SymmetryRegistry` (`register_symmetries`, `D2_symmetry`, `Cn_symmetry`). Interfaces between sub-domains are marked `Transparent`, and cuts along symmetry axes `SymmetryWall`. See the files in `src/geometry/billiards/` for templates.

# Holes
A billiard with obstacles uses a `MultiplyConnectedDomain(outer, holes, id)` as its fundamental domain, where `holes` is a vector of `DomainHole`s. `AnnularBilliard`, `SinaiBilliard` and `SquareWithinSquareBilliard` are built this way, with `AnnularBilliard` the reference example:

```julia
b = AnnularBilliard(0.3, 0.1)                 # inner radius, shift of the hole centre

genus(b)                                       # number of hole loops
get_outer_boundary_curves(b)                   # outer boundary
get_hole_boundary_curves(b)                    # hole boundaries
boundary_components(b)                         # closed loops of the boundary
boundary_component_is_hole(b)                  # which loops are holes
cmap = full_boundary_components(b)            # BoundaryComponentMap: curves, is_hole, group
area(b)                                        # outer area minus hole areas
```

Holes live in the full physical frame and are never replicated by the symmetry registry. Every hole loop is traversed with the hole region on its left. `validate_holes` and `symmetry_compatible` check whether holes are consistent with a given registry; a billiard with holes and a nonempty registry is not supported by symmetry-reduced consumers.

# Symmetries
A billiard stores its fundamental domain only. The registered symmetries are used to
- reconstruct the full boundary (`full_boundary`),
- map boundary points between symmetry-related images (`symmetry_orbit`, `SymmetryOrbitMap`),
- map phase space points between sectors (`apply_symmetry_pb`, `pb_sectors`),
- supply the symmetry data needed to separate eigenstates into symmetry sectors in QuantumBilliards.jl.

Generators are registered with `register_symmetries(gens...)`, which assigns each a stable `sym_id`. `SymmetryWall` curves refer to their generator by this id, and `symmetry_of(registry, sym_id)` retrieves it independently of storage order.

# Related packages
- [QuantumBilliards.jl](https://github.com/Quantum-Chaos-Julia/QuantumBilliards.jl) - spectra and eigenstates of quantum billiards.
- `ClassicalBilliards.jl` - classical billiard dynamics, trajectories and Poincaré sections (under development).
- [BilliardGeometryPlotting.jl](https://github.com/Quantum-Chaos-Julia/BilliardGeometryPlotting.jl) - plotting of curves, domains and billiards.
- [QBPlotting.jl](https://github.com/Quantum-Chaos-Julia/QBPlotting.jl) - plotting of quantum states, boundary functions and Husimi functions.
- [SpectralStatistics.jl](https://github.com/Quantum-Chaos-Julia/SpectralStatistics.jl) - spectral statistics and random matrix models.

# Roadmap
- Billiards on the hyperbolic plane and in a constant magnetic field.
- Additional billiard families and more general domain composition tools.
- Stable shared interface for classical dynamics in `ClassicalBilliards.jl`.

# Citation
If you use BilliardGeometry.jl through QuantumBilliards.jl in your research, please cite the QuantumBilliards.jl paper, see the [QuantumBilliards.jl README](https://github.com/Quantum-Chaos-Julia/QuantumBilliards.jl#citation).

# License
See [LICENSE](LICENSE).
