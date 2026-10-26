# Shared fixtures for the hole / multiply connected domain unit tests (Phase 1 of the multiply connected domain rework).
# Not a test file: it defines helpers and a minimal billiard type only, no @testset. It is included first in unittests.jl, so every
# file after it may use the `hf_` helpers. All fixtures are cheap (no eigen solves, no quadrature beyond Green-theorem areas).
#
# Reference billiards (all in the full physical frame unless stated):
#   hf_concentric_annulus()  AnnularBilliard(0.3, 0.0): unit disk without the disk of radius 0.3 about the origin, area pi (1 - 0.09)
#   hf_eccentric_annulus()   AnnularBilliard(0.2, 0.3): unit disk without the disk of radius 0.2 about (0.3, 0), area pi (1 - 0.04)
#   hf_rect_with_hole()      rectangle [0,3] x [0,2] without the rectangle [1,1.5] x [0.5,1.5] (off center), no symmetry, area 6 - 0.5
#   hf_bitten_square()       outer-only SimpleDomain: [0,1]^2 without the quarter disk of radius 0.5 about 0 (inverted circle, orientation -1), area 1 - pi/16
#   hf_d2_with_hole()        D2 billiard: fundamental frame is the quarter [0,1]^2 (full frame [-1,1]^2) with the asymmetric circular
#                            DomainHole of radius 0.2 about (0.3, 0.4) (domain_id 5), area 4 - pi 0.04

using StaticArrays: SVector

# Minimal billiard holding a MultiplyConnectedDomain and a symmetry registry (the production holes only exist in AnnularBilliard).
struct hf_HoleBilliard{T} <: AbsBilliard where {T<:Real}
    fundamental_domain::MultiplyConnectedDomain{T}
    symmetries::SymmetryRegistry
end

# Closed counterclockwise rectangle [x0, x1] x [y0, y1] of four specular LineSegments (segment ids 1..4) as a SimpleDomain with the
# given domain id; traversed with the rectangle on the left, as required for outer regions and holes alike.
function hf_rect_domain(x0, y0, x1, y1, id)
    p = (SVector(x0, y0), SVector(x1, y0), SVector(x1, y1), SVector(x0, y1))
    curves = AbsCurve[LineSegment(p[k], p[mod1(k + 1, 4)]; domain_id = id, segment_id = k) for k in 1:4]
    return SimpleDomain{typeof(float(x0))}(curves, collect(p), id)
end

# Full counterclockwise circle of radius R about (x0, y0) as a SimpleDomain with the given domain id (segment id 1).
function hf_circle_domain(R, x0, y0, id)
    T = typeof(float(R))
    crv = CircleSegment(R, T(2 * pi), zero(T), [x0, y0]; domain_id = id, segment_id = 1)
    return SimpleDomain{T}(AbsCurve[crv], SVector{2,T}[SVector(x0 + R, y0)], id)
end

# Composite hole of two disjoint rectangles [0,1]x[0,1] (domain id `id`) and [2,3]x[0,1] (domain id `id + 1`): one hole object, two loops.
hf_two_loop_hole(id = 3) = DomainHole(CompositeDomain(AbsSimpleDomain[hf_rect_domain(0.0, 0.0, 1.0, 1.0, id), hf_rect_domain(2.0, 0.0, 3.0, 1.0, id + 1)]))

hf_concentric_annulus() = AnnularBilliard(0.3, 0.0)
hf_eccentric_annulus() = AnnularBilliard(0.2, 0.3)

function hf_rect_with_hole()
    outer = hf_rect_domain(0.0, 0.0, 3.0, 2.0, 1)
    hole = DomainHole(hf_rect_domain(1.0, 0.5, 1.5, 1.5, 2))
    return hf_HoleBilliard{Float64}(MultiplyConnectedDomain(outer, [hole], 1), register_symmetries())
end

# Outer region: the quarter square [0,1]^2 closed by the two symmetry walls of RectangleBilliard(1.0, 1.0).
hf_quarter_square_domain() = RectangleBilliard(1.0, 1.0).fundamental_domain

function hf_d2_with_hole()
    outer = hf_quarter_square_domain()
    hole = DomainHole(hf_circle_domain(0.2, 0.3, 0.4, 5))
    return hf_HoleBilliard{Float64}(MultiplyConnectedDomain(outer, [hole], 1), D2_symmetry())
end

# Unit square [0,1]^2 without the quarter disk of radius 0.5 about the origin: four lines and one inverted (orientation -1) arc.
# The arc runs clockwise from (0, 0.5) to (0.5, 0) so that the region is on the left; the region lies outside the circle.
function hf_bitten_square()
    l1 = LineSegment([0.5, 0.0], [1.0, 0.0]; domain_id = 1, segment_id = 1)
    l2 = LineSegment([1.0, 0.0], [1.0, 1.0]; domain_id = 1, segment_id = 2)
    l3 = LineSegment([1.0, 1.0], [0.0, 1.0]; domain_id = 1, segment_id = 3)
    l4 = LineSegment([0.0, 1.0], [0.0, 0.5]; domain_id = 1, segment_id = 4)
    arc = CircleSegment(0.5, -pi / 2, pi / 2, [0.0, 0.0]; orientation = -1, domain_id = 1, segment_id = 5)
    return SimpleDomain{Float64}(AbsCurve[l1, l2, l3, l4, arc], SVector{2,Float64}[[0.5, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.5]], 1)
end
