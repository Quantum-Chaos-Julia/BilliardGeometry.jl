# Shared tolerances and reference helpers for the integration tier. No testsets in this file.
# Tolerance values are calibrated in scratchpad/test-suite-calibration.md; each constant is at least 10x the worst
# observed value on the real fixtures (floored at 1e-12 where the observed deviation is exactly zero).

using LinearAlgebra: norm
using StaticArrays: SVector

# --- tolerances -------------------------------------------------------------------------------------------------------

# Maximum gap between the end of a curve and the start of the next one in a reconstructed boundary loop (PROVISIONAL)
const CLOSURE_ATOL = 1e-6
# Absolute error of a total boundary length against its analytic value (PROVISIONAL)
const PERIMETER_ATOL = 1e-6
# Absolute error of the shoelace area of a densely sampled boundary against its analytic value (PROVISIONAL)
const AREA_ATOL_SHOELACE = 1e-3
# Distance between a boundary node mapped by a symmetry and the node the orbit map assigns to it (PROVISIONAL)
const ORBIT_ATOL = 1e-6
# Distance between a point of a polar-curve image and the radial function of its source curve (PROVISIONAL)
const POLAR_IMAGE_ATOL = 1e-6
# Slack of the Poincare-Birkhoff coordinate s around the window [(k - 1) L, k L] of its sector (PROVISIONAL)
const PB_WINDOW_ATOL = 1e-6
# Difference between s of a sector image and the cumulative length of the full boundary at that image (PROVISIONAL)
const PB_ARCLENGTH_ATOL = 1e-6
# Difference between pb_sectors and the cumulative curve lengths of the full boundary (PROVISIONAL)
const PB_SECTOR_ATOL = 1e-6

# --- reference helpers ------------------------------------------------------------------------------------------------

"""
    uniform_boundary_samples(curves, N) → pts::Vector{SVector{2,Float64}}

Sample `N` points on a connected boundary at equal arc-length spacing, at the midpoints of the `N` periodic cells,
that is at arc length `L (k - 1/2) / N` for `k = 1, ..., N` with `L` the total length of `curves`.

## Description
This is the node grid used by boundary-integral discretizations of a complete boundary component, built here only from
the curve geometry (`component_lengths` and `curve`), independently of the symmetry-orbit bookkeeping under test.

## Arguments
* `curves`: Connected boundary curves in traversal order, e.g. `full_boundary(billiard)`.
* `N::Int`: Number of sample points.

## Returns
* `pts`: The `N` boundary points in traversal order.
"""
function uniform_boundary_samples(curves, N::Int)
    _, cum, Ltot = component_lengths(curves)
    pts = Vector{SVector{2,Float64}}(undef, N)
    for k in 1:N
        target = Ltot * (k - 0.5) / N
        j = clamp(searchsortedlast(cum, target), 1, length(curves))
        u = clamp((target - cum[j]) / curves[j].length, 0.0, 1.0)
        pts[k] = curve(curves[j], u)
    end
    return pts
end

"""
    closure_gaps(curves) → gaps::Vector{Float64}

Distance between the end point of every curve and the start point of the next one, cyclically (the last curve is
followed by the first). All gaps vanish for a closed, continuously traversed boundary loop.
"""
function closure_gaps(curves)
    n = length(curves)
    return [norm(curve(curves[i], 1.0) - curve(curves[mod1(i + 1, n)], 0.0)) for i in 1:n]
end

"""
    boundary_signed_area(curves; n = 500) → A::Float64

Shoelace area of the polygon through `n` equally spaced parameter values of every curve. It is positive for a
counterclockwise traversal and converges to the enclosed area as `n` grows (error of order `1/n^2`).
"""
function boundary_signed_area(curves; n::Int = 500)
    pts = SVector{2,Float64}[]
    for c in curves, t in range(0.0, 1.0; length = n + 1)[1:end-1]
        push!(pts, curve(c, t))
    end
    A = 0.0
    for i in eachindex(pts)
        p, q = pts[i], pts[mod1(i + 1, length(pts))]
        A += p[1] * q[2] - q[1] * p[2]
    end
    return A / 2
end
