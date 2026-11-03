# Unit tests for the annular billiard constructor and structure:
#   src/geometry/billiards/annular.jl   (AnnularBilliard{T}, AnnularBilliard(R_inner, shift; R_outer, center))
# Reference objects (the hole center is displaced by `shift` along the x axis, so the mirror y -> -y is always a symmetry):
#   * center = (0, 0) (default): symmetric construction. The outer region is the upper half disc, a SimpleDomain of the CircleSegment arc
#     (angle pi, domain id 1) and a SymmetryWall(1, 2) LineSegment on the x axis; the hole is a SymmetricHole whose piece is the upper
#     half disc of the hole (CircleWedge, domain id 2) and whose registry is {XAxisReflection}. The registry of the billiard is
#     {XAxisReflection}, so G' has order 2, the physical full boundary is two loops (outer circle, hole circle) of two blocks each.
#   * nonzero center: unreduced Phase 1 construction. The outer region is a SimpleDomain holding the full counterclockwise CircleSegment
#     (orientation +1, domain id 1), the hole is a DomainHole around a SimpleDomain holding the full counterclockwise CircleSegment
#     (domain id 2, segment id 1); the registry is empty.
# Closed forms: area pi (R_outer^2 - R_inner^2) (unchanged by the reduction), fundamental area = area / 2, loop lengths 2 pi R_outer and
# 2 pi R_inner. Guards: R_outer > 0, R_inner > 0, R_inner < R_outer, and the hole must lie inside the outer circle (validate_holes).
# get_boundary_curves rejects the symmetric billiard (holes plus a non-empty registry). No eigen solves.

using StaticArrays: SVector

# Total arc length of a loop of LineSegments and CircleSegments.
ann_loop_length(loop) = sum(c -> c isa LineSegment ? c.length : c.radius * c.arc_angle, loop)

@testset "AnnularBilliard fields" begin
    ann = AnnularBilliard(0.3, 0.2)
    @test ann isa AnnularBilliard{Float64}
    @test ann isa AbsBilliard
    @test ann.fundamental_domain isa MultiplyConnectedDomain{Float64}
    @test ann.symmetries isa SymmetryRegistry
    @test fieldnames(AnnularBilliard) == (:fundamental_domain, :symmetries)
    @test @inferred(AnnularBilliard(0.3, 0.2)) isa AnnularBilliard{Float64}
end

@testset "AnnularBilliard symmetries" begin
    # origin center: registry {XAxisReflection}, G' = G of order 2
    ann = AnnularBilliard(0.3, 0.2)
    @test length(ann.symmetries) == 1
    @test ann.symmetries[1] isa XAxisReflection
    @test typeof.(get_symmetries(ann)) == (XAxisReflection,)
    @test length(get_symmetry_group(ann)) == 1
    @test length(AnnularBilliard(0.3, 0.0).symmetries) == 1
    @test length(AnnularBilliard(0.3, 0.2; center = [0.0, 0.0]).symmetries) == 1
    # nonzero center: no symmetry (Phase 1)
    shifted = AnnularBilliard(0.3, 0.2; center = [0.5, 0.1])
    @test isempty(shifted.symmetries)
    @test get_symmetries(shifted) == ()
    @test length(get_symmetry_group(shifted)) == 0
end

@testset "AnnularBilliard curves" begin
    dom = AnnularBilliard(0.3, 0.2).fundamental_domain
    @test dom.id == 1
    # one outer subdomain (half disc: arc + symmetry wall), one SymmetricHole with two half-disc pieces
    @test length(dom.subdomains) == 1
    @test length(dom.holes) == 1
    @test dom.holes[1] isa SymmetricHole{Float64}
    outer_dom = dom.subdomains[1]
    @test outer_dom isa SimpleDomain{Float64}
    @test outer_dom.id == 1
    @test length(outer_dom.boundary) == 2
    arc, wall = outer_dom.boundary
    @test arc isa CircleSegment{Float64}
    @test wall isa LineSegment{Float64}
    @test arc.bc isa SpecularReflection
    @test wall.bc isa SymmetryWall
    @test [arc.orientation, wall.orientation] == [1, 1]
    @test [arc.domain_id, wall.domain_id] == [1, 1]
    @test [arc.segment_id, wall.segment_id] == [1, 2]
    # upper half of the outer circle of radius one about the origin, closed by the x axis from (-1, 0) to (1, 0)
    @test arc.radius == 1.0
    @test arc.center == SVector(0.0, 0.0)
    @test arc.arc_angle ≈ pi atol = 1e-14
    @test wall.pt0 == SVector(-1.0, 0.0)
    @test wall.pt1 == SVector(1.0, 0.0)

    # the hole: upper half disc piece of radius R_inner about (shift, 0), registry {XAxisReflection} acting about the hole center
    hole = dom.holes[1]
    @test hole.piece isa CircleWedge{Float64}
    @test hole.piece.id == 2
    @test length(hole.symmetries) == 1
    @test hole.symmetries[1] isa XAxisReflection
    @test hole.center == SVector(0.2, 0.0)
    @test length(hole.pieces) == 2
    @test hole.id_stride == 100
    hole_arc = hole.piece.boundary[1]
    @test hole_arc isa CircleSegment{Float64}
    @test hole_arc.radius == 0.3
    @test hole_arc.center == SVector(0.2, 0.0)
    @test hole_arc.arc_angle ≈ pi atol = 1e-14
    @test hole_arc.bc isa SpecularReflection
    @test count(c -> c.bc isa SymmetryWall, hole.piece.boundary) == 2
end

@testset "AnnularBilliard curves (shifted center)" begin
    # nonzero center: the Phase 1 construction with full circles
    dom = AnnularBilliard(0.3, 0.2; center = [0.5, 0.1]).fundamental_domain
    @test dom.id == 1
    @test length(dom.subdomains) == 1
    @test length(dom.holes) == 1
    @test dom.holes[1] isa DomainHole
    outer_dom = dom.subdomains[1]
    hole_dom = dom.holes[1].domain
    @test outer_dom isa SimpleDomain{Float64}
    @test hole_dom isa SimpleDomain{Float64}
    @test [outer_dom.id, hole_dom.id] == [1, 2]
    @test length(outer_dom.boundary) == 1
    @test length(hole_dom.boundary) == 1
    outer = outer_dom.boundary[1]
    inner = hole_dom.boundary[1]
    @test outer isa CircleSegment{Float64}
    @test inner isa CircleSegment{Float64}
    @test outer.bc isa SpecularReflection
    @test inner.bc isa SpecularReflection
    # both circles are counterclockwise (orientation +1, hole region on the left), with their own domain id
    @test [outer.orientation, inner.orientation] == [1, 1]
    @test [outer.domain_id, inner.domain_id] == [1, 2]
    @test [outer.segment_id, inner.segment_id] == [1, 1]
    @test outer.arc_angle == 2 * pi
    @test inner.arc_angle == 2 * pi
    @test outer.radius == 1.0
    @test outer.center == SVector(0.5, 0.1)
    @test inner.radius == 0.3
    @test inner.center ≈ SVector(0.7, 0.1) atol = 1e-14
end

@testset "AnnularBilliard topology" begin
    ann = AnnularBilliard(0.3, 0.2)
    @test genus(ann) == 1
    @test genus(ann.fundamental_domain) == 1
    @test !hasfield(MultiplyConnectedDomain, :genus)
    # boundary components of the fundamental domain: the physical outer arc and the physical arcs of the two hole pieces
    comps = boundary_components(ann)
    @test length(comps) == genus(ann) + 1
    @test [length(c) for c in comps] == [1, 2]
    @test [[c.domain_id for c in comp] for comp in comps] == [[1], [2, 102]]
    @test boundary_component_is_hole(ann) == [false, true]
    @test length(geometric_boundary_components(ann)) == 3
    # the two pieces of the hole are two hole domain ids of one component
    @test hole_domain_ids(ann) == [2, 102]
    @test component_of_domain_id(ann) == Dict(1 => 0, 2 => 1, 102 => 1)
    # the half-disc outer region has two curves, each half-disc hole piece three
    @test length(get_all_curves(ann)) == 2 + 2 * 3

    # shifted center: one component for the outer circle and one for the hole, one curve each
    shifted = AnnularBilliard(0.3, 0.2; center = [0.5, 0.1])
    comps_s = boundary_components(shifted)
    @test length(comps_s) == genus(shifted) + 1
    @test [length(c) for c in comps_s] == [1, 1]
    @test [c[1].domain_id for c in comps_s] == [1, 2]
    @test boundary_component_is_hole(shifted) == [false, true]
    @test length(geometric_boundary_components(shifted)) == 2
    @test hole_domain_ids(shifted) == [2]
    @test component_of_domain_id(shifted) == Dict(1 => 0, 2 => 1)
end

@testset "AnnularBilliard geometry" begin
    # area pi (R_outer^2 - R_inner^2) of the full physical billiard, concentric and eccentric; unchanged by the symmetry reduction
    ann = AnnularBilliard(0.3, 0.1)
    shifted = AnnularBilliard(0.3, 0.1; center = [0.5, 0.1])
    @test area(ann) ≈ pi * (1 - 0.3^2) rtol = 1e-12
    @test area(AnnularBilliard(0.3, 0.0)) ≈ pi * (1 - 0.3^2) rtol = 1e-12
    @test area(AnnularBilliard(0.2, 0.3)) ≈ pi * (1 - 0.2^2) rtol = 1e-12
    @test area(AnnularBilliard(0.3, 0.1; R_outer = 2.0)) ≈ pi * (4 - 0.3^2) rtol = 1e-12
    @test area(shifted) ≈ pi * (1 - 0.3^2) rtol = 1e-12
    # the fundamental area is the reduced-domain area: area / group order
    @test fundamental_area(ann) ≈ area(ann) / 2 rtol = 1e-12
    @test fundamental_area(shifted) ≈ area(shifted) rtol = 1e-12
    # membership refers to the reduced domain (upper half plane) of the symmetric billiard: ring, hole, outside
    @test is_inside(ann, SVector(0.8, 0.0)) === true
    @test is_inside(ann, SVector(0.1, 0.0)) === false      # hole center (shift 0.1)
    @test is_inside(ann, SVector(1.2, 0.0)) === false
    # a ring point close to the hole is inside, a hole point is not
    @test is_inside(ann, SVector(0.1, 0.35)) === true
    @test is_inside(ann, SVector(0.1, 0.25)) === false
    # corners of the reduced domain: two right angles where the wall meets the outer arc, two 3 pi / 2 where it meets the hole
    @test sort(corner_angles(ann)) ≈ [pi / 2, pi / 2, 3pi / 2, 3pi / 2] atol = 1e-12
    # no true corners in the full physical billiard
    @test isempty(corner_angles(ann; fundamental = false))
    @test isempty(corner_angles(shifted))
    # full physical boundary: outer loop (two half circles) then hole loop (two half circles), blocks [2, 2], loop lengths 2 pi R
    @test length(full_boundary(ann)) == 4
    fbc = full_boundary_components(ann)
    @test fbc.is_hole == [false, true]
    @test [length(b) for b in fbc.blocks] == [2, 2]
    @test [ann_loop_length(loop) for loop in fbc.curves] ≈ [2pi, 2pi * 0.3] rtol = 1e-12
    big = full_boundary_components(AnnularBilliard(0.3, 0.1; R_outer = 2.0))
    @test [ann_loop_length(loop) for loop in big.curves] ≈ [4pi, 2pi * 0.3] rtol = 1e-12
    # shifted center: two loops of one curve each, one block each
    fbc_s = full_boundary_components(shifted)
    @test length(full_boundary(shifted)) == 2
    @test fbc_s.is_hole == [false, true]
    @test [length(b) for b in fbc_s.blocks] == [1, 1]
    @test [ann_loop_length(loop) for loop in fbc_s.curves] ≈ [2pi, 2pi * 0.3] rtol = 1e-12
end

@testset "AnnularBilliard get_boundary_curves" begin
    # holes plus a non-empty registry: the outer curves live in the fundamental frame, the hole curves in the full frame
    @test_throws ArgumentError get_boundary_curves(AnnularBilliard(0.3, 0.2))
    # unreduced construction: outer circle then hole circle
    @test length(get_boundary_curves(AnnularBilliard(0.3, 0.2; center = [0.5, 0.1]))) == 2
end

@testset "AnnularBilliard keywords" begin
    ann = AnnularBilliard(0.3, 0.2; R_outer = 2.0, center = [0.1, 0.0])
    dom = ann.fundamental_domain
    outer = dom.subdomains[1].boundary[1]
    inner = dom.holes[1].domain.boundary[1]
    @test outer.radius == 2.0
    @test outer.center == SVector(0.1, 0.0)
    @test inner.radius == 0.3
    @test inner.center ≈ SVector(0.3, 0.0) atol = 1e-14
    @test area(ann) ≈ pi * (4 - 0.09) rtol = 1e-12
    # R_outer scales the symmetric construction too
    sym = AnnularBilliard(0.3, 0.2; R_outer = 2.0)
    @test sym.fundamental_domain.subdomains[1].boundary[1].radius == 2.0
    @test sym.fundamental_domain.subdomains[1].boundary[2].pt0 == SVector(-2.0, 0.0)
    @test length(sym.symmetries) == 1
end

@testset "AnnularBilliard argument guards" begin
    # R_outer > 0
    @test_throws ArgumentError AnnularBilliard(0.3, 0.2; R_outer = 0.0)
    @test_throws ArgumentError AnnularBilliard(0.3, 0.2; R_outer = -1.0)
    # R_inner > 0
    @test_throws ArgumentError AnnularBilliard(0.0, 0.2)
    @test_throws ArgumentError AnnularBilliard(-0.3, 0.2)
    # R_inner < R_outer
    @test_throws ArgumentError AnnularBilliard(1.0, 0.0)
    @test_throws ArgumentError AnnularBilliard(1.5, 0.0)
    @test_throws ArgumentError AnnularBilliard(0.5, 0.0; R_outer = 0.5)
    # the sample point of the hole loop must lie inside the outer circle (validate_holes)
    @test_throws ArgumentError AnnularBilliard(0.3, 3.0)
    @test_throws ArgumentError AnnularBilliard(0.3, 3.0; center = [0.5, 0.1])
end

@testset "AnnularBilliard eltype promotion" begin
    ann32 = AnnularBilliard(0.3f0, 0.2f0)
    @test ann32 isa AnnularBilliard{Float32}
    @test ann32.fundamental_domain isa MultiplyConnectedDomain{Float32}
    @test all(c -> c isa Union{CircleSegment{Float32},LineSegment{Float32}}, get_all_curves(ann32))
    @test genus(ann32) == 1
    @test area(ann32) isa Float32
    @test fundamental_area(ann32) isa Float32
    shifted32 = AnnularBilliard(0.3f0, 0.2f0; center = Float32[0.5, 0.1])
    @test all(c -> c isa CircleSegment{Float32}, get_all_curves(shifted32))
    @test area(shifted32) isa Float32
end
