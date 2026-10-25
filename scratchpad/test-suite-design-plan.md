# BilliardGeometry.jl Test Suite Redesign — Implementation Plan

## Objective
Replace the monolithic `test/runtests.jl` (564 lines, 11 testsets, 611 assertions) with the two-tier layout used by QuantumBilliards.jl (unit tier mirroring `src/`, integration tier organized by computed quantity, Aqua quality check), port what is worth keeping, and close the coverage gaps found below. Source defects G1-G13 found during analysis are **fixed in the same effort, test first**: each gets a failing regression test before its source fix. The package moves to `julia = "1.12"` (decided), which also enables `[sources]` in `test/Project.toml`.

## Context / Current State
Facts below come from a real `Pkg.test(coverage=true)` run and two probe scripts (`/tmp/bgp/probe*.jl`), not assumptions.

* `test/runtests.jl`: one file, no groups, no `test/Project.toml` (`Test` lives in `[extras]/[targets]` of the main `Project.toml`). All 611 assertions pass in about 34 s wall time, most of it compilation. `symmetryorbits.jl` (379 assertions, 7.4 s) and `poincarebirkhoff.jl` (4.6 s) are the slow testsets.
* `testscripts.jl` (repo root, 17 lines) is a scratch script using Revise/CairoMakie, not a test.
* Empty testset: `limaconsegment.jl` (0 assertions). The name hints at `LimaconSegment`, which is **exported but not defined** anywhere in `src/`.
* Line coverage by file (covered/total). Files not listed are 100 %.

| File | Coverage | Gap |
|---|---|---|
| `geometry/area.jl` | 0/51 | `area` (all 4 methods), `fundamental_area`, `corner_angles`, `symmetry_reduction_factor`, `maximal_symmetry` |
| `geometry/billiards/{annular,ellipse,limacon,mushroom,polygon,prosen,sinai,star}.jl` | 0 % | constructors and argument validation |
| `geometry/domains/multiplyconnecteddomains.jl` | 0/4 | constructor |
| `geometry/inversions.jl` | 10/47 (21 %) | `invert_curve` all methods, grid search, Nelder-Mead |
| `geometry/boundarytypes.jl` | 28/99 (28 %) | `get_boundary_curves` (composite, multiply connected), `get_all_domains`, `get_domain`, `get_all_curves`, `get_curve`, `update_boundary_condition`, `genus`, `boundary_components`, `geometric_boundary_components` |
| `geometry/geometry.jl` | 8/25 (32 %) | `is_inside` (billiard, domain, curve, arrays), `domain_gradient_vector` array forms |
| `geometry/domains/compositedomains.jl` | 2/7 | `reset_ids!` |
| `geometry/boundarycomponents.jl` | 28/72 (39 %) | `print_component_junctions`, `_global_t_to_segment_u`, `_eval_composite_geom_global_t`, `_component_corner_angles` |
| `geometry/arclength.jl` | 9/21 (43 %) | `arc_length(crv, ts)`, `arc_length(crv, pt)` for generic curves, `construct_arc_length_interpolation` |
| `quadrature/samplers.jl` | 27/39 (69 %) | `LinearNodes`, `GaussLegendreNodes`, `FourierNodes` with `primes` |
| `geometry/utils.jl` | 64/84 (76 %) | `connect_curves` (code carries `#TODO: Let's test this function.`), `find_unique_elements`, `point_curve_parameter`, `signed_angle` |
| `geometry/segments/polarcurves.jl` | 33/42 | array methods, `polar_radius(::AbstractArray)` |
| `geometry/curvederivatives.jl` | 60/78 | `PolarSegment` `tangent_2`, `FourierCoeffPolarSegment` second derivative path |

Well covered by existing tests: `symmetryorbits.jl` (95 %), `symmetryregistry.jl` (100 %), `poincarebirkhoff.jl` (97 %), `kressgrading.jl` (92 %), `fullboundary.jl` (82 %).

* Existing testsets violate the layout rules in `julia-test-suite-structure`: the symmetry-orbit and full-boundary testsets loop over several real billiards and rebuild boundary samples (these are integration checks sitting in a unit file); no tolerance constants; very long provenance comments referencing migration "Step 10.9 / Step 15"; hand-rolled `_point_at_sigma` duplicated twice.
* QuantumBilliards.jl already works around one source defect: `test/integration/high_k/beyn_high_k_spectra.jl` line 22 states "`area(billiard)` ignores holes" and recomputes hole-aware areas by hand. Downstream validation of geometry through spectra already exists there, so this package's integration tier stays solver-free.

### Source defects found by probing (fixed in this plan, see "Source Fixes")
Ranked by downstream impact. Observations come from probe runs on Julia 1.12.5 unless marked "from code reading".

| ID | Where | Observation | Expected |
|---|---|---|---|
| G1 | `area.jl` `area(billiard)` | `AnnularBilliard(0.3, 0.2)` gives 3.4243 (= π(1+0.09)); `SinaiBilliard(1,0.3)` gives 4.2827 (= 4+0.2827); `SquareWithinSquareBilliard(1.0)` gives 4.125 (= 4+2·0.0625). Hole area is added: hole curves are stored CCW with `orientation=-1`, and `_area_integral` ignores `orientation`. | 2.8588, 3.7173, 3.875 |
| G2 | `area.jl` `corner_angles`, `boundarycomponents.jl` | `SinaiBilliard` fundamental gives `[0.0, 0.0, π/2, π/2, π/2]`: the inner arc is parametrized p0→p1 (CCW, `orientation=-1`) but the fundamental loop reaches it at p1, so chaining by `connect_curves` leaves it disconnected and the junction tangents are inconsistent. `SquareWithinSquareBilliard` hole corners give 0.7854 and 1.5708 instead of 3π/4 and 3π/2 (`_junction_angle` is unsigned, so reflex corners cannot be represented). `fundamental=false` treats the whole `full_boundary` as one component, mixing outer and hole loops. | Sinai: `[π/2, π/2, π/2, π/2, π/2]`; SWS hole corners reflex |
| G3 | `utils.jl` `is_overlaping` | `isapprox` without `atol` fails near zero: `(-1.8e-16, -1.0)` vs `(6e-17, -1.0)` is `false`, so `is_closed(full_boundary(CircleBilliard(1.0)))` is `false` (`is_connected` gives `[1,1,1,0]`). Also affects `connect_curves`. | `true` |
| G4 | `arclength.jl` / `circlesegment.jl` `arc_length(::CircleSegment, pt)` | negative for `arc_angle<0` (−0.785 vs 0.785), for lower-half points (−1.571 vs 4.712), and when `shift_angle` pushes the arc across ±π (−4.712 vs 1.571). | non-negative, equals `R*|arc_angle|*t` |
| G5 | `inversions.jl` `invert_curve` (`:roots`, `:auto`) | returns 0.370032 for true 0.37 (error 3e-5, root-finding on a double root of squared distance). `:optim` and `:grid_optim` are exact to 1e-16. | agree with `:optim` |
| G6 | `geometry.jl` `is_inside(::SquareWithinSquareBilliard)` | `is_inside(billiard)` ANDs the half-plane test of every curve; the exterior of a convex polygonal hole is a UNION of half-planes, so the result is almost empty: Monte-Carlo estimate of the fundamental area 0.1436 vs about 1.94. Circular holes (single curve; Annular, Sinai) are fine (2.858 vs 2.859). | about 1.94 |
| G7 | `boundarytypes.jl` `update_boundary_condition` on composite-domain billiards | `StadiumBilliard`: `MethodError: no constructor StadiumBilliard(::CompositeDomain, ::SymmetryRegistry)`. Works for `SimpleDomain` billiards. | returns updated billiard |
| G8 | `boundarytypes.jl` `get_all_curves(::AbsSimpleDomain)` | `subdomains = domain` then iterates a domain: `MethodError iterate(::SimpleDomain)`. | curve vector |
| G9 | `utils.jl` `find_unique_elements` | writes `vector[i] = false` instead of `u[i] = false`, mutating the input (`[0,0,1,1,2]` becomes `[0,0,0,1,2]`). Result is correct only by accident. | input unchanged |
| G10 | `BilliardGeometry.jl` | `export angle` shadows `Base.angle`; `using BilliardGeometry` then `angle(a,b)` in `Main` raises an ambiguity error. `LimaconSegment` exported but undefined (Aqua `undefined_exports`). | rename / drop |
| G11 | `curvederivatives.jl` | `tangent_vec`, `normal_vec` have only array methods; `tangent`, `tangent_2`, `curvature` have scalar methods. | consistent API |
| G12 | `Project.toml` | `julia = "1.9"`; several deps have no `[compat]` (Accessors, CircularArrays, CoordinateTransformations, ForwardDiff, Rotations, StaticArrays, LinearAlgebra). Test deps live in `[extras]/[targets]`. Aqua not run yet, so the list is predicted. | `julia = "1.12"`, compat for every dep, test deps in `test/Project.toml` |
| G13 | `fullboundary.jl` `_apply_symmetry_to_curve(::LineSegment)` (from code reading, to confirm by probe) | always sets `orientation=-c.orientation`, even for orientation-preserving images (`XYAxisReflection`, `NFoldRotation`) where `CircleSegment`/polar images keep it. Reflections are fixed by the later `_reverse_curve` (second flip), so only det>0 images of lines carry a wrong `orientation`, e.g. the XY image of a `RectangleBilliard` wall. | `orientation` preserved for det>0, flipped for det<0 |

Hole convention (must not change): hole curves are CCW-parametrized with `orientation=-1`. `QuantumBilliards.jl/src/solvers/sweepmethods/compositebim.jl` (`_group_is_hole`, `_flip_component_orientation`, `evaluate_points`) relies on it and flips hole normals itself. Reversing hole traversal in the constructors would double-flip and break CompositeBIM, so G1/G2/G6 are fixed in the measure/point-location code, not in the constructors.

## Proposed Approach
1. Dismantle `runtests.jl` into the QuantumBilliards layout. Port every existing assertion that is a legitimate unit check; move cross-billiard geometric checks to the integration tier; drop the empty testset and `testscripts.jl`-style content.
2. Unit tier: one explicit `@testset` per function or type, mirroring `src/`, no loops over billiards, small deterministic fixtures with analytic answers (`@inferred` on hot paths).
3. Integration tier for a solver-free geometry package: validate **pipelines that downstream solvers rely on** against **analytic or cross-method references**, organized by computed quantity. Reference priority: analytic > cross-method agreement (Green area vs point-in-domain quadrature, full boundary vs orbit images) > fixtures.
4. Known defects G1-G13 are fixed test-first: write the regression test (fails on current `src/`), apply the fix from "Source Fixes", confirm green. Each regression testset carries a searchable `# G<n>: <one line>` comment. `@test_broken` is used only if a fix is deferred, and then it is listed in the final report.
5. Fixes are gated by the downstream suites: after each batch of source fixes run QuantumBilliards.jl `QB_TEST_GROUP=unit` (and the `bim_sweep`/`high_k` integration files for batches touching holes, `connect_curves`, `arc_length` or `invert_curve`).

## Source Fixes
Order is chosen so earlier fixes unblock later tests. Every fix is minimal and local; no new exports except where stated.

| Order | ID | File | Fix |
|---|---|---|---|
| 1 | G3 | `geometry/utils.jl` | `is_overlaping(pt1, pt2)`: compare with `isapprox(a, b; atol=sqrt(eps(T)))` per coordinate, `T` from the point element type (relative-only tolerance fails at 0). Check `connect_curves`/`is_closed` results on all 15 billiards are unchanged except the previously false ones. `connect_curves` is used by QuantumBilliards `boundarypoints.jl`. |
| 2 | G9 | `geometry/utils.jl` | `find_unique_elements`: write `u[i] = false`, add the same `atol` as G3 so `0.0` and `1e-17` merge (matters for `pb_sectors`). Behavior: keep the last of each run of near-equal values; input unchanged. |
| 3 | G10 | `BilliardGeometry.jl` | Remove `angle` from `export` (internal use in `polygons.jl` stays); remove `LimaconSegment` from `export`. Grep over QuantumBilliards, QBPlotting, BilliardGeometryPlotting, QCPlotting, SpectralStatistics and QuantumBilliardsTests found no users of either. |
| 4 | G12 | `Project.toml`, `test/Project.toml` | `julia = "1.12"`; `[compat]` for every non-stdlib dep (use the major versions resolved in `Manifest.toml`) and for stdlibs only if Aqua `deps_compat` demands it; drop `[extras]`/`[targets]`; create `test/Project.toml` with `[sources] BilliardGeometry = {path = ".."}`, `Aqua`, `Elliptic`, `Test`, compat `julia = "1.12"`. Cross-package: QuantumBilliards.jl depends on this package, so its `Project.toml` and `test/Project.toml` compat must also be raised to `julia = "1.12"` (currently `1.11` in the test project). |
| 5 | G8 | `geometry/boundarytypes.jl` | `get_all_curves(domain::AbsDomain)`: for non-composite use `subdomains = (domain,)` instead of iterating the domain itself. |
| 6 | G7 | `geometry/boundarytypes.jl`, `geometry/billiards/{stadium,mushroom}.jl` | `@set` on a billiard rebuilds it through the bare type name; `CompositeDomain` has no `T`, so no constructor matches. Add `StadiumBilliard(fd::CompositeDomain, sym::SymmetryRegistry)` and the same for `MushroomBilliard`, forwarding to `StadiumBilliard{T}` with `T` taken from the first subdomain (one small `_domain_eltype(::CompositeDomain)` helper in `compositedomains.jl`). Alternative if preferred: add `ConstructionBase.constructorof` methods; rejected, new dependency. Also covers the nested `@set` into `subdomains[dom_idx].boundary[seg_idx].bc`. |
| 7 | G13 | `geometry/fullboundary.jl` | `_apply_symmetry_to_curve(sym, ::LineSegment)`: `orientation = det(_sym_matrix(sym, T)) > 0 ? c.orientation : -c.orientation`. Confirm with a probe first (`full_boundary(RectangleBilliard(1.0,0.6))` line orientations per sector) and with `domain_fun` of an image curve evaluated at an interior point of the unfolded domain. Does not change any geometry or `_reverse_curve`. |
| 8 | G4 | `geometry/segments/circlesegment.jl` | `arc_length(::CircleSegment, pt)`: `Δθ = sign(arc_angle)*(atan(y-cy, x-cx) - shift_angle)`, wrapped into `[-ε, 2π-ε)` with `ε=sqrt(eps)`, return `radius*Δθ`. Result is non-negative and equals `R*|arc_angle|*t` for points on the arc. Known limit: the end point of a full 2π arc maps to 0, not to the length; the test pins this, callers that need `t=1` use the `t` method. Used by QuantumBilliards `fundamental_pb_coords` through `pb_coords`; the existing PB tests must stay green. |
| 9 | G5 | `geometry/inversions.jl` | `invert_curve_roots`: solve `g(t) = <γ(t)-pt, γ'(t)> = 0` (stationary point of the distance, simple root) with Newton and a ForwardDiff derivative instead of root-finding on the squared distance (double root, linear convergence). After the solve require `|γ(t)-pt| < sqrt(atol)+rtol*scale`, else throw so `method=:auto` falls back to `:optim`. Keeps the existing signatures and the `error` for unknown methods. |
| 10 | G11 | `geometry/curvederivatives.jl` | Add scalar methods `tangent_vec(crv, t::Real)` and `normal_vec(crv, t::Real)`; array methods unchanged. |
| 11 | G1 | `geometry/area.jl` | `area(billiard)`: compute the outer loop and the holes separately. Hole ids = `domain_id`s of curves with `orientation == -1` in `get_boundary_curves(billiard)` (same rule as CompositeBIM `_group_is_hole`, read from the fundamental curves because image orientations are unreliable). `area = |Σ I_outer|/2 - Σ_h |Σ I_h|/2` over `full_boundary` curves partitioned by those ids. No hole: identical to today. `_area_integral`, the curve/composite/vector methods and `fundamental_area` unchanged. Expected: Annular 2.8588, Sinai 3.7173, SWS 3.875. |
| 12 | G6 | `geometry/geometry.jl` | `is_inside(billiard, pt)` / `is_inside(domain, pt)`: reduce per-curve booleans by group. Curves with `orientation == -1` are grouped by `domain_id` and reduced with `any` (outside a hole = on the billiard side of at least one hole edge); all other curves are reduced with `all`; result is `all` of the outer group AND `all` over hole groups. A single circular hole is the one-element case, so Annular/Sinai are unchanged. The vector-returning `is_inside(domain, pt)` API keeps its meaning; add the reducer as a new internal `_domain_contains(domain, pt)` used by the `AbsBilliard` method. |
| 13 | G2 | `geometry/boundarycomponents.jl`, `geometry/area.jl`, `geometry/boundarytypes.jl` | (a) add `_signed_junction_angle` (sign of the 2D cross product); `_junction_angle` and `_is_true_corner` stay unsigned, so Kress corner locations are unchanged; (b) `_component_corner_angles` returns `π - τ_signed` for a component traversed with the billiard on the left; (c) new internal `_connect_oriented_curves` for the fundamental boundary: like `connect_curves` but, when no forward match exists, may append `_reverse_curve` of a curve with `orientation == -1`; used only by `_geometric_boundary` for domains that contain such curves, so `get_boundary_curves`/`connect_curves` behavior seen by QuantumBilliards is untouched; (d) `corner_angles(...; fundamental=false)`: split `full_boundary` into the outer loop and one loop per hole id (as in G1) instead of one component; hole loops are reversed (CW) before the angle computation. Expected: Sinai fundamental all `π/2`; SWS hole corners `3π/2` (full domain) and `3π/4` where a hole edge meets the symmetry wall. |

Regression tests per fix (all live, none `@test_broken`):
* G3 `is_closed(full_boundary(CircleBilliard(1.0)))`; G9 input unchanged and result sorted-unique; G10 Aqua `undefined_exports` and `angle` ambiguity in a fresh `Main`; G12 Aqua `deps_compat`, `julia` compat string.
* G8, G7 `get_all_curves` on a `SimpleDomain`, `update_boundary_condition` on Stadium and Mushroom with return type and untouched original; G13 line orientations of all sectors give correct `is_inside` side at a known interior point.
* G4 table of (R, arc_angle sign, shift, t) including a negative arc and a wrap across ±π; G5 `:roots`, `:auto`, `:optim`, `:grid_optim` agree to `1e-10` on line, circle, polar; G11 scalar vs array agreement.
* G1 analytic areas of the three hole billiards; G6 Monte-Carlo-free grid area of SWS within `INSIDE_AREA_RTOL_GRID`; G2 corner-angle tables above, plus `sum` identities (outer loop angles sum to `(n-2)π` for polygons).

## Target Layout
```
test/
  runtests.jl                       # BG_TEST_GROUP = unit | integration | all
  Project.toml                      # Aqua, Elliptic (ellipse perimeter), Test; [sources] BilliardGeometry = {path = ".."}
  unit/
    unittests.jl                    # include list only
    geometry/
      utils_tests.jl                # is_overlaping, connect_curves, is_connected, is_closed, signed_angle, find_unique_elements, point_curve_parameter
      symmetry_tests.jl             # _sym_matrix, apply_symmetry, orientation_reversing, NFoldRotation, CompositeReflection ==/hash, D2_symmetry, Cn_symmetry, apply_symmetry_pb
      symmetryregistry_tests.jl     # register_symmetries ids, symmetry_of, AbstractVector behaviour
      symmetryorbits_tests.jl       # symmetry_node_multiple, SymmetryOrbitMap accessors/partition, get_symmetries, error paths (port of pure parts)
      boundarytypes_tests.jl        # BC types, get_boundary_curves, get_all_*, get_curve/get_domain, update_boundary_condition, genus, boundary_components, geometric_boundary_components
      boundarycomponents_tests.jl   # component_lengths, junction angles, corner locations, global-t mapping, print_component_junctions (capture output)
      fullboundary_tests.jl         # _apply_symmetry_to_curve and _reverse_curve per curve type, _full_boundary_data provenance
      arclength_tests.jl            # arc_length scalar/array/point, construct_arc_length_interpolation
      curvederivatives_tests.jl     # tangent/tangent_2 vs finite differences, tangent_vec/normal_vec/curvature
      inversions_tests.jl           # invert_curve all methods on line/circle/polar, error branch
      geometry_tests.jl             # is_inside, domain_fun, domain_gradient_vector
      area_tests.jl                 # area (curve, composite, vector), symmetry_reduction_factor, maximal_symmetry, fundamental_area, corner_angles on polygons
      poincarebirkhoff_tests.jl     # fundamental_pb_coords, get_pb_curve, pb_coords, pb_sectors
      segments/{linesegment,circlesegment,polarcurves,compositecurves}_tests.jl
      domains/{polygons,circular,compositedomains,multiplyconnecteddomains}_tests.jl
      billiards/{triangle,stadium,mushroom,polar,limacon,circle,ellipse,rectangle,polygon,star,c3,prosen,annular,sinai,square_in_square}_tests.jl
    quadrature/{samplers,kressgrading}_tests.jl
    quality/aqua_tests.jl           # one testset per Aqua check; live after fixes G10/G12
  integration/
    common.jl                       # tolerances + analytic generators, no tests
    integrationtests.jl             # include list only
    measures/                       # area_perimeter.jl, area_with_holes.jl
    boundary/                       # fullboundary_closure.jl, fundamental_vs_full.jl
    point_location/                 # inside_area_quadrature.jl
    symmetry/                       # orbit_maps_billiards.jl, image_invariance.jl
    corners/                        # corner_angles.jl
    discretization/                 # graded_boundary_quadrature.jl, fourier_nodes_boundary.jl
    poincare_birkhoff/              # pb_unfolding.jl
    precision/                      # eltype_consistency.jl
```
Billiard unit files are deliberately tiny (struct fields, registry contents, `SymmetryWall` ids, `@test_throws ArgumentError` for each guard in the constructor, `eltype` promotion). Annular, Sinai and SquareWithinSquare unit files carry the structural (not area) checks.

## Port / Rework / Drop Map for Existing Testsets
| Old testset (`runtests.jl`) | Decision | Destination |
|---|---|---|
| `linesegment.jl` | port, extend (array methods, `arc_length`, orientation sign, `Float32`) | `unit/geometry/segments/linesegment_tests.jl` |
| `circlesegment.jl` | port, extend; add G4 regression tests | `unit/geometry/segments/circlesegment_tests.jl` |
| `limaconsegment.jl` (empty) | drop | none; `LimaconBilliard` structural tests go to `billiards/limacon_tests.jl` |
| `polarsegment.jl` | port; keep the type-promotion and finite-difference checks | `segments/polarcurves_tests.jl` |
| `curvederivatives.jl` | port nearly as is; tighten FD checks with `h = cbrt(eps)` and per-order `atol` | `unit/geometry/curvederivatives_tests.jl` |
| `boundarycomponents.jl` | port; add global-t mapping and junction printing | `unit/geometry/boundarycomponents_tests.jl` |
| `kressgrading.jl` | port; add `kress_R_even!`/`kress_R_odd!` symmetry (`R[i,j]=R[j,i]`, parity) and `_corner_interval` | `unit/quadrature/kressgrading_tests.jl` |
| `samplers.jl - FourierNodes` | port; add `LinearNodes` and `GaussLegendreNodes` (weights sum to 1, exact for polynomials) | `unit/quadrature/samplers_tests.jl` |
| `symmetry.jl - CompositeReflection` | port | `unit/geometry/symmetry_tests.jl` |
| `symmetryorbits.jl` (379 asserts) | **split**: the partition/size/`get_symmetries`/`symmetry_of` error checks stay unit; the "sample the real full boundary and verify each orbit geometrically" loops over `CircleBilliard`, `C3Billiard`, `SquareWithinSquareBilliard`, a polar billiard with `AntiDiagonalReflection` move to integration | `unit/geometry/symmetryorbits_tests.jl`, `integration/symmetry/orbit_maps_billiards.jl` |
| `fullboundary.jl` | **split**: `_apply_symmetry_to_curve` / `_reverse_curve` per curve type and trivial-symmetry identity stay unit; stadium/C3/polar closure and index-permutation checks move to integration | `unit/geometry/fullboundary_tests.jl`, `integration/boundary/fullboundary_closure.jl` |
| `boundarytypes.jl - SymmetryWall` | port | `unit/geometry/boundarytypes_tests.jl` |
| `poincarebirkhoff.jl` | **split**: `pb_sectors` / `pb_coords` structural invariant over `RectangleBilliard` and `C3Billiard` loops is integration; unit gets one-curve `fundamental_pb_coords` and `apply_symmetry_pb` identities | `unit/geometry/poincarebirkhoff_tests.jl`, `integration/poincare_birkhoff/pb_unfolding.jl` |
| root `testscripts.jl` | delete (scratch, Revise/CairoMakie); nothing to port | none |

Shared helper: the duplicated `_point_at_sigma` / `_sample_full_boundary` becomes one function in `integration/common.jl` (arc-length-uniform midpoint sampler of a boundary vector). Unit files needing it define a locally prefixed copy (`fb_`).

## Proposed Integration Tests
All solver-free, no eigen-solves, deterministic grids instead of RNG, each file runs in seconds. One testset per billiard where analytic references differ; no solver x geometry product. Each geometry class is hit by at least two files.

1. **`measures/area_perimeter.jl`**: `area(b)` and total `full_boundary` length vs closed forms, one testset per exported billiard without holes.
   * Circle `πR²`, `2πR`; Rectangle(a,b) `4ab`, `4(a+b)`; Ellipse `πab`, perimeter `4a E(e²)` via `Elliptic.E`; Stadium(hw) `π+4hw`, `2π+4hw`; Mushroom `πR²/2 + 2·hw·h`; Limaçon `πR²(1+a²/2)`; Star/C3/Prosen/polar: `½∫r²dφ` and `∫√(r²+r'²)dφ` by independent `QuadGK` on `r(φ)`; Triangle(π/2,1) `1/2`, `2+√2`; Polygon unit square.
   * `fundamental_area(b) * symmetry_reduction_factor(b) ≈ area(b)`.
   * Verified now: circle, rectangle, ellipse, stadium, mushroom, limaçon, star, polygon, triangle all match analytic values to quadrature accuracy; C3 and Prosen need an independent reference in the test.
2. **`measures/area_with_holes.jl`**: Annular `π(R²−r²)`, Sinai `4a²−πR²`, SquareWithinSquare `4a²−Σs²`; `genus`, `length(boundary_components)` = genus+1. Live assertions (G1 fixed). Same file cross-checks `area(b)` against the QuantumBilliards-style recomputation (outer minus holes grouped by `domain_id`), so the two code paths cannot drift apart.
3. **`point_location/inside_area_quadrature.jl`**: cross-method area. Deterministic midpoint grid over the bounding box, count `is_inside(billiard, pt)`, compare with `fundamental_area` within a grid-resolution tolerance (`O(h)` boundary error, tolerance scaled by perimeter·h). Covers Circle, Stadium, Mushroom, Limaçon, Star, Annular, Sinai and SquareWithinSquare (G6 fixed; the last is the case that exposed the defect). Independent of Green's theorem, so it exposes sign/orientation errors like G1/G6 from the other side.
4. **`boundary/fullboundary_closure.jl`**: for every billiard with a single component: `full_boundary` is a connected loop (endpoint distance < `CLOSURE_ATOL`), has CCW orientation (signed Green area > 0), every curve has positive length, total length is `(nsym+1)·` fundamental physical length, and `is_closed(full_boundary(CircleBilliard(1.0)))` is `true` (G3 regression). Covers D2, Cn, single-reflection, diagonal-reflection and no-symmetry billiards. Ported from the stadium/C3/polar parts of the old testset.
5. **`boundary/fundamental_vs_full.jl`**: points sampled on the full boundary are on the boundary of the unfolded domain (`min` over curves of `|domain_fun|` < tol) and on exactly one fundamental-domain image; fundamental-domain interior samples map by every registered symmetry into the full-domain interior. Cross-checks `get_boundary_curves`, `full_boundary`, `apply_symmetry` and `domain_fun` against each other.
6. **`symmetry/orbit_maps_billiards.jl`**: ported geometric orbit verification (X/Y/XY reflections on `CircleBilliard`, `NFoldRotation` on `C3Billiard`, `DiagonalReflection` on SquareWithinSquare, `AntiDiagonalReflection` on the shifted polar circle, `CompositeReflection` closure equals `XYAxisReflection`). Single parametrized helper instead of five copies.
7. **`symmetry/image_invariance.jl`**: for each registered generator of each symmetric billiard, `apply_symmetry` of full-boundary samples lies on the full boundary (distance via `invert_curve(...; method=:optim)` below tolerance); `orientation_reversing` agrees with the sign of the Jacobian determinant of `_sym_matrix`; `|G| = length(registry)+1` equals `symmetry_reduction_factor` for pure reflection/rotation billiards.
8. **`corners/corner_angles.jl`**: Polygon/Rectangle/Triangle: angle sum `(n−2)π`, known angles (`TriangleBilliard(γ,χ)` returns `(α, β, γ)` with `α+β+γ=π`); smooth billiards (Circle `fundamental=false`, Ellipse, Stadium, Limaçon full) have no corners; fundamental-domain corners of the D2 billiards are `π/2`, C3 wedge `2π/3`. `SinaiBilliard` (all `π/2` in the fundamental domain) and SquareWithinSquare hole corners (`3π/2` full, `3π/4` at the symmetry wall) are live (G2 fixed). This is the geometric input of the Weyl law used by the solvers.
9. **`discretization/graded_boundary_quadrature.jl`**: the pipeline `component_lengths` → `multi_kress_graded_nodes_data` with `_component_corner_locations` → `_eval_composite_geom_global_t`. Using the Kress-graded global parameter on rectangle, stadium fundamental boundary and a polygon, integrate `∮ ds` (perimeter), `½∮(x dy − y dx)` (area) and a corner-singular smooth integrand; assert error decreases geometrically or reaches `GRADED_QUAD_ATOL` at fixed N, and that graded nodes beat uniform nodes at the corner integrand. This is exactly what BIM solvers consume.
10. **`discretization/fourier_nodes_boundary.jl`**: `FourierNodes(lengths)` with and without `primes` on the real fundamental boundary of several billiards: weights sum to total length, node spacing uniform across curve junctions (generalizes the existing three-length case), trapezoid rule of `∮ ds` and of a smooth periodic function converges spectrally with N.
11. **`poincare_birkhoff/pb_unfolding.jl`**: ported `pb_sectors`/`pb_coords` invariants plus (a) a velocity tangent to the boundary gives `|p|=1`, along the inward normal gives `p=0`; (b) reflection sectors reverse `s` and negate `p`, rotation sectors preserve them; (c) `s` is monotone over a full boundary sweep and continuous at sector joins; (d) `pb_sectors` breakpoints agree with cumulative `component_lengths` of `full_boundary`.
12. **`precision/eltype_consistency.jl`**: Float32, Float64 and BigFloat for Rectangle, Stadium, Star, C3 and Limaçon: result types are preserved (`@inferred` on `area`, `curve`, `tangent`), values agree across precisions within `eps(Float32)`-scaled tolerance, `NFoldRotation(...; T)` matches billiard `T`.

Not proposed: spectral/Weyl-law counting against analytic rectangle or right-triangle spectra. It would duplicate QuantumBilliards.jl's `integration/symmetry` and `high_k` tiers and needs `N(k)` remainder tolerances with no tight analytic bound; revisit only if a geometry-only Weyl check is wanted.

### Integration `common.jl` contents
* Tolerances (each commented `# <family> (worst observed <value>)`, at least 10x the worst observed, calibrated by a probe run and logged in `scratchpad/test-suite-calibration.md`): `AREA_RTOL_ANALYTIC`, `PERIMETER_RTOL_ANALYTIC`, `CLOSURE_ATOL`, `INSIDE_AREA_RTOL_GRID`, `ORBIT_ATOL`, `GRADED_QUAD_ATOL`, `PRECISION_RTOL_F32`.
* Reference generators with docstrings: `analytic_area(b)`, `analytic_perimeter(b)` dispatching on billiard type (hole billiards subtract), `polar_area(r)`, `polar_perimeter(r, dr)`, `uniform_boundary_samples(curves, N)`, `bounding_box(curves)`, `grid_inside_area(billiard, n)`.
* Fixture list `NO_HOLE_BILLIARDS` / `HOLE_BILLIARDS` as plain tuples of `(name, constructor)`; testsets stay explicit, the list is only used by helper functions inside one testset.

## Steps
Test-first throughout: scaffold, port (all green on unmodified `src/`), then for each source fix write the failing regression test, apply the fix, re-run.

1. Scaffold: apply fix 4 (G12) first so the environment is final: `julia = "1.12"`, full `[compat]`, `test/Project.toml` with `[sources]`, remove `[extras]/[targets]`. Rewrite `test/runtests.jl` as pure group dispatch with `BG_TEST_GROUP`. Create `test/unit/unittests.jl` and `test/integration/{common.jl,integrationtests.jl}` skeletons, files commented until they exist. Run the old testsets once under the new runner to confirm 611/611 before dismantling.
2. Port batch A (mechanical, no new logic): linesegment, circlesegment, polarcurves, curvederivatives, boundarycomponents, kressgrading, samplers, symmetry (CompositeReflection), boundarytypes (SymmetryWall). Assertion count must not drop.
3. Split batch B: symmetryorbits, fullboundary, poincarebirkhoff into unit and integration halves; delete the old testsets only after the new files reproduce every assertion.
4. Probe-confirm G13 (line orientations of `full_boundary` images) and the QuantumBilliards dependency map (grep already done: only `connect_curves` of the touched functions is used downstream; `angle`/`LimaconSegment` unused).
5. Source-fix batch 1, simple and local (fixes 1-3, 5, 6, 10 = G3, G9, G10, G8, G7, G11): write regression tests in `unit/geometry/{utils,boundarytypes,curvederivatives}_tests.jl`, `unit/quality/aqua_tests.jl`, `unit/geometry/billiards/{stadium,mushroom}_tests.jl`; apply fixes. Gate: BilliardGeometry unit tier + QuantumBilliards `QB_TEST_GROUP=unit`.
6. Source-fix batch 2, curve-level (fixes 7-9 = G13, G4, G5): tests in `unit/geometry/fullboundary_tests.jl`, `segments/circlesegment_tests.jl`, `inversions_tests.jl`; apply fixes. Gate: BilliardGeometry unit tier + `poincarebirkhoff` integration file + QuantumBilliards unit tier (uses `arc_length`).
7. New unit files for the remaining gap files (area, geometry, arclength, compositecurves, domains, symmetryregistry, the 15 billiard constructor files, the quadrature additions). `area_tests.jl` and `geometry_tests.jl` stay failing for hole cases until step 8.
8. Source-fix batch 3, hole handling (fixes 11-13 = G1, G6, G2): finish `area_tests.jl`, `geometry_tests.jl`, `boundarycomponents_tests.jl` regression tests; apply fixes in the order G1, G6, G2. Gate: QuantumBilliards `bim_sweep/compositebim_wavenumber.jl`, `high_k` hole testsets and `symmetry/reflection_spectra.jl` hole cases (they run Annular, Sinai, SWS through CompositeBIM and must not change).
9. Integration tier: write `common.jl` first, run each probe once to record worst observed values, set tolerances, then add files in the order measures, boundary, point_location, symmetry, corners, discretization, poincare_birkhoff, precision. With the fixes in, hole billiards are covered by live assertions everywhere.
10. Cleanup: remove old `runtests.jl` content and root `testscripts.jl`; add a "Running tests" section to `README.md` with the `BG_TEST_GROUP` commands; add a CI workflow (`julia-ci-setup`) with Julia 1.12 only; raise QuantumBilliards.jl julia compat to 1.12 (see fix 4).
11. Final: full run of `unit`, `integration`, `all` here and QuantumBilliards `QB_TEST_GROUP=all`; coverage run; record per-file coverage and tolerance calibration in `scratchpad/test-suite-calibration.md`; report any deferred `@test_broken`.

## Risks & Open Questions
* Julia 1.12 is decided. Consequence: QuantumBilliards.jl (depends on this package) must also move to 1.12; any other package in the workspace that depends on BilliardGeometry needs the same compat bump.
* Cross-package convention risk (highest): CompositeBIM treats `orientation == -1` curves as CCW-parametrized holes and flips their normals. Fixes 11-13 deliberately keep that convention; any later change of hole traversal must change both packages together. Gate in step 8 protects this.
* `connect_curves` tolerance change (fix 1) can alter which curves chain in edge cases; the check over all 15 billiards plus the QuantumBilliards `boundarypoints.jl` users is part of batch 1.
* Fix 13(c) adds a second chaining routine next to `connect_curves`. If review prefers one routine, extend `connect_curves` with an opt-in keyword instead, accepting a larger diff in a function QuantumBilliards calls.
* G13 is from code reading only; if the probe shows the field is never consumed for images, fix 7 may shrink to a test plus a one-line correctness fix.
* `is_inside` reducer (fix 12) assumes every `orientation == -1` group is a hole. True for all 15 shipped billiards; user-built domains violating it would change behavior, documented in the docstring.
* Grid-based `inside_area_quadrature.jl` tolerance is first order in the grid step; if too coarse to separate G1-sized errors (0.28 on 4.0 for Sinai) use Richardson or a finer grid on small billiards.
* Aqua has not been run; the G10/G12 expectations are predicted and finalized in step 5.
* Slow items today: `symmetryorbits` (7.4 s) and `poincarebirkhoff` (4.6 s). Target after split: unit tier well under 30 s, integration under 2 min on one thread.

## Testing & Validation
* Delegate test authoring to the Julia Test Writer agent per batch (A, B, gap files, integration files), passing the G-id table and the verified probe values above as evidence.
* Acceptance: `BG_TEST_GROUP=unit`, `integration`, and `all` each pass on Julia 1.12; QuantumBilliards `QB_TEST_GROUP=unit` and the compositebim/hole integration files still pass; no tolerance inlined in integration files; every regression testset has a `# G<n>` comment; no `@test_broken` remains (or each remaining one is listed with reason); line coverage of every `src/` file at least 90 % except explicitly listed dead code (`chebyshev_nodes`/`fourier_nodes` commented-out blocks, `invert_curve_nelder_mead`, `print_component_junctions` if left untested).
* Calibration log in `scratchpad/test-suite-calibration.md` with worst observed values per tolerance.

## Out of Scope
* Changing the hole convention (CCW traversal, `orientation=-1`) or CompositeBIM's normal flipping.
* New geometry features, new billiards, API renames beyond the G10 export removals.
* Solver, spectrum or Weyl-law validation (covered by QuantumBilliards.jl and SpectralStatistics.jl).
* Plotting packages (`BilliardGeometryPlotting.jl`, `QBPlotting.jl`) and their tests.
