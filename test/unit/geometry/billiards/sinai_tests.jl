# Unit tests for the Sinai billiard constructor and structure:
#   src/geometry/billiards/sinai.jl   (SinaiBilliard{T}, SinaiBilliard(a, R; center))
# Reference objects (square [-a, a]^2 with a concentric circular obstacle of radius R):
#   * center = (0, 0) (default): symmetric construction. The outer region is the fundamental triangle (0,0), (a,0), (a,a): a SimpleDomain
#     of one physical side (x = a, domain id 1), SymmetryWall(1, 2) on the diagonal and SymmetryWall(7, 8) on the x axis. The obstacle is
#     a SymmetricHole whose piece is the eighth disc (CircleWedge, angle pi/4, domain id 2) with the same two walls, registry
#     D4_symmetry() (7 non-identity elements, G' = D4 of order 8, id_stride 100, 8 pieces). Full physical boundary via
#     full_boundary_components: loop lengths [8 a, 2 pi R], blocks [8, 8].
#   * nonzero center: unreduced Phase 1 construction: SimpleDomain of the four counterclockwise LineSegments bottom, right, top, left
#     (orientation 1, domain id 1, segment ids 1..4) and a DomainHole around the full counterclockwise CircleSegment of radius R about
#     `center` (orientation +1, domain id 2, segment id 1); empty registry.
# Closed forms: area 4 a^2 - pi R^2 (unchanged by the reduction), fundamental area = area / 8 for the symmetric billiard, membership in
# the square / hole / outside. Guards: a > 0, R > 0, R < a. get_boundary_curves rejects the symmetric billiard. No eigen solves.

using StaticArrays: SVector

# Total arc length of a loop of LineSegments and CircleSegments.
sin_loop_length(loop) = sum(c -> c isa LineSegment ? c.length : c.radius * c.arc_angle, loop)

@testset "SinaiBilliard fields" begin
    sin_b = SinaiBilliard(1.0, 0.3)
    @test sin_b isa SinaiBilliard{Float64}
    @test sin_b isa AbsBilliard
    @test sin_b.fundamental_domain isa MultiplyConnectedDomain{Float64}
    @test sin_b.symmetries isa SymmetryRegistry
    @test fieldnames(SinaiBilliard) == (:fundamental_domain, :symmetries)
    @test @inferred(SinaiBilliard(1.0, 0.3)) isa SinaiBilliard{Float64}
end

@testset "SinaiBilliard symmetries" begin
    # origin center: registry D4_symmetry(), G' = G of order 8
    sin_b = SinaiBilliard(1.0, 0.3)
    @test length(sin_b.symmetries) == 7
    @test [g.sym_id for g in sin_b.symmetries] == [g.sym_id for g in D4_symmetry()]
    @test length(get_symmetry_group(sin_b)) == 7
    @test typeof.(get_symmetries(sin_b)) == (DiagonalReflection, YAxisReflection)
    @test length(SinaiBilliard(1.0, 0.3; center = [0.0, 0.0]).symmetries) == 7
    # nonzero center: no symmetry (Phase 1)
    shifted = SinaiBilliard(1.0, 0.3; center = [0.5, -0.2])
    @test isempty(shifted.symmetries)
    @test get_symmetries(shifted) == ()
    @test length(get_symmetry_group(shifted)) == 0
end

@testset "SinaiBilliard curves" begin
    dom = SinaiBilliard(1.0, 0.3).fundamental_domain
    @test dom.id == 1
    # one outer subdomain (fundamental triangle), one SymmetricHole of eight pieces
    @test length(dom.subdomains) == 1
    @test length(dom.holes) == 1
    @test dom.holes[1] isa SymmetricHole{Float64}
    outer_dom = dom.subdomains[1]
    @test outer_dom isa SimpleDomain{Float64}
    @test outer_dom.id == 1
    @test length(outer_dom.boundary) == 3
    @test all(c -> c isa LineSegment{Float64}, outer_dom.boundary)
    side, diagonal, axis = outer_dom.boundary
    # one physical side and two symmetry walls (diagonal and x axis), all counterclockwise with domain id 1
    @test side.bc isa SpecularReflection
    @test diagonal.bc isa SymmetryWall
    @test axis.bc isa SymmetryWall
    @test [c.orientation for c in outer_dom.boundary] == [1, 1, 1]
    @test [c.domain_id for c in outer_dom.boundary] == [1, 1, 1]
    @test [c.segment_id for c in outer_dom.boundary] == [1, 2, 3]
    # triangle (0,0), (a,0), (a,a) traversed (a,0) -> (a,a) -> (0,0) -> (a,0)
    @test [side.pt0, side.pt1] == [SVector(1.0, 0.0), SVector(1.0, 1.0)]
    @test [diagonal.pt0, diagonal.pt1] == [SVector(1.0, 1.0), SVector(0.0, 0.0)]
    @test [axis.pt0, axis.pt1] == [SVector(0.0, 0.0), SVector(1.0, 0.0)]
    @test outer_dom.corners == [SVector(1.0, 0.0), SVector(1.0, 1.0), SVector(0.0, 0.0)]

    # obstacle: eighth disc of radius R about the origin (domain id 2) and its 8 images under D4 (id_stride 100)
    hole = dom.holes[1]
    @test hole.piece isa CircleWedge{Float64}
    @test hole.piece.id == 2
    @test length(hole.symmetries) == 7
    @test length(hole.pieces) == 8
    @test hole.id_stride == 100
    @test hole.center == SVector(0.0, 0.0)
    arc = hole.piece.boundary[1]
    @test arc isa CircleSegment{Float64}
    @test arc.radius == 0.3
    @test arc.center == SVector(0.0, 0.0)
    @test arc.arc_angle ≈ pi / 4 atol = 1e-14
    @test arc.bc isa SpecularReflection
    @test count(c -> c.bc isa SymmetryWall, hole.piece.boundary) == 2
end

@testset "SinaiBilliard curves (shifted center)" begin
    # nonzero center: the Phase 1 construction
    dom = SinaiBilliard(1.0, 0.3; center = [0.5, -0.2]).fundamental_domain
    @test dom.id == 1
    @test length(dom.subdomains) == 1
    @test length(dom.holes) == 1
    @test dom.holes[1] isa DomainHole
    outer_dom = dom.subdomains[1]
    hole_dom = dom.holes[1].domain
    @test outer_dom isa SimpleDomain{Float64}
    @test hole_dom isa SimpleDomain{Float64}
    @test [outer_dom.id, hole_dom.id] == [1, 2]
    @test length(outer_dom.boundary) == 4
    @test length(hole_dom.boundary) == 1
    @test all(c -> c isa LineSegment{Float64,SpecularReflection}, outer_dom.boundary)
    inner = hole_dom.boundary[1]
    @test inner isa CircleSegment{Float64,SpecularReflection}
    # the square and the circle are counterclockwise (orientation +1, hole region on the left), each with its own domain id
    @test [c.orientation for c in outer_dom.boundary] == [1, 1, 1, 1]
    @test [c.domain_id for c in outer_dom.boundary] == [1, 1, 1, 1]
    @test [c.segment_id for c in outer_dom.boundary] == [1, 2, 3, 4]
    @test inner.orientation == 1
    @test inner.domain_id == 2
    @test inner.segment_id == 1
    @test inner.arc_angle == 2 * pi
    # sides in the order bottom, right, top, left, each starting where the previous one ends
    sides = outer_dom.boundary
    @test sides[1].pt0 ≈ SVector(-0.5, -1.2) atol = 1e-14
    @test sides[3].pt0 ≈ SVector(1.5, 0.8) atol = 1e-14
    @test all(i -> sides[i].pt1 == sides[mod1(i + 1, 4)].pt0, 1:4)
    @test all(s -> s.length == 2.0, sides)
    @test outer_dom.corners == [s.pt1 for s in sides]
    @test inner.center ≈ SVector(0.5, -0.2) atol = 1e-14
    @test inner.radius == 0.3
    @test hole_dom.corners ≈ [SVector(0.8, -0.2)] atol = 1e-14
end

@testset "SinaiBilliard topology" begin
    sin_b = SinaiBilliard(1.0, 0.3)
    @test genus(sin_b) == 1
    @test genus(sin_b.fundamental_domain) == 1
    # components of the fundamental domain: the physical side alone, and the physical arcs of the eight obstacle pieces
    comps = boundary_components(sin_b)
    @test length(comps) == genus(sin_b) + 1
    @test [length(c) for c in comps] == [1, 8]
    @test [[c.domain_id for c in comp] for comp in comps] == [[1], [2 + 100k for k in 0:7]]
    @test boundary_component_is_hole(sin_b) == [false, true]
    @test length(geometric_boundary_components(sin_b)) == 3
    @test hole_domain_ids(sin_b) == [2 + 100k for k in 0:7]
    @test component_of_domain_id(sin_b) == Dict(1 => 0, [2 + 100k => 1 for k in 0:7]...)
    # three outer curves (side + two walls) and three curves per obstacle piece; only the obstacle curves are hole curves
    @test [is_hole_curve(sin_b, c) for c in get_all_curves(sin_b)] == [fill(false, 3); fill(true, 24)]

    # shifted center: one component for the square (four curves) and one for the circle
    shifted = SinaiBilliard(1.0, 0.3; center = [0.5, -0.2])
    comps_s = boundary_components(shifted)
    @test genus(shifted) == 1
    @test [length(c) for c in comps_s] == [4, 1]
    @test [[c.domain_id for c in comp] for comp in comps_s] == [[1, 1, 1, 1], [2]]
    @test boundary_component_is_hole(shifted) == [false, true]
    @test length(geometric_boundary_components(shifted)) == 2
    @test hole_domain_ids(shifted) == [2]
    @test component_of_domain_id(shifted) == Dict(1 => 0, 2 => 1)
    @test [is_hole_curve(shifted, c) for c in get_all_curves(shifted)] == [false, false, false, false, true]
end

@testset "SinaiBilliard geometry" begin
    sin_b = SinaiBilliard(1.0, 0.3)
    shifted = SinaiBilliard(1.0, 0.3; center = [0.5, -0.2])
    # area 4 a^2 - pi R^2 of the full physical billiard (unchanged by the reduction); fundamental area is area / 8
    @test area(sin_b) ≈ 4 - pi * 0.3^2 rtol = 1e-12
    @test area(SinaiBilliard(2.0, 0.5)) ≈ 16 - pi * 0.25 rtol = 1e-12
    @test area(shifted) ≈ 4 - pi * 0.3^2 rtol = 1e-12
    @test fundamental_area(sin_b) ≈ area(sin_b) / 8 rtol = 1e-12
    @test fundamental_area(shifted) ≈ area(shifted) rtol = 1e-14
    # membership of the symmetric billiard refers to the fundamental triangle 0 <= y <= x <= a minus the obstacle piece
    @test is_inside(sin_b, SVector(0.7, 0.3)) === true
    @test is_inside(sin_b, SVector(0.5, 0.2)) === true
    @test is_inside(sin_b, SVector(0.9, 0.05)) === true
    @test is_inside(sin_b, SVector(0.2, 0.1)) === false     # inside the hole piece
    @test is_inside(sin_b, SVector(0.25, 0.1)) === false
    @test is_inside(sin_b, SVector(0.3, 0.7)) === false     # outside the fundamental triangle (y > x)
    @test is_inside(sin_b, SVector(-0.5, -0.5)) === false   # other quadrants lie outside the fundamental domain
    @test is_inside(sin_b, SVector(1.1, 0.5)) === false
    # membership of the unreduced billiard: the full square minus the disk, in all four quadrants
    @test is_inside(shifted, SVector(0.5, 0.5)) === true
    @test is_inside(shifted, SVector(-0.3, -0.5)) === true
    @test is_inside(shifted, SVector(-0.4, 0.1)) === true
    @test is_inside(shifted, SVector(0.5, -0.2)) === false  # hole center
    @test is_inside(shifted, SVector(0.6, -0.1)) === false
    @test is_inside(shifted, SVector(2.0, 0.5)) === false
    # corners of the full physical square: four right angles; unreduced billiard likewise
    @test corner_angles(sin_b; fundamental = false) ≈ fill(pi / 2, 4) atol = 1e-12
    @test corner_angles(shifted) ≈ fill(pi / 2, 4) atol = 1e-12
    @test corner_angles(shifted; fundamental = false) ≈ fill(pi / 2, 4) atol = 1e-12
    # corners of the reduced domain: includes the pi/4 angles of the fundamental triangle and the corners where walls meet the obstacle
    @test length(corner_angles(sin_b)) == 6
    @test count(a -> isapprox(a, pi / 4; atol = 1e-12), corner_angles(sin_b)) == 2
    # full physical boundary: square loop (8 pieces) and circle loop (8 arcs), blocks [8, 8], loop lengths 8 a and 2 pi R
    fbc = full_boundary_components(sin_b)
    @test length(full_boundary(sin_b)) == 16
    @test fbc.is_hole == [false, true]
    @test [length(b) for b in fbc.blocks] == [8, 8]
    @test [sin_loop_length(loop) for loop in fbc.curves] ≈ [8.0, 2pi * 0.3] rtol = 1e-12
    fbc2 = full_boundary_components(SinaiBilliard(2.0, 0.5))
    @test [sin_loop_length(loop) for loop in fbc2.curves] ≈ [16.0, 2pi * 0.5] rtol = 1e-12
    # shifted center: square loop of four sides and circle loop of one curve, one block each
    fbc_s = full_boundary_components(shifted)
    @test length(full_boundary(shifted)) == 5
    @test fbc_s.is_hole == [false, true]
    @test [length(b) for b in fbc_s.blocks] == [1, 1]
    @test [sin_loop_length(loop) for loop in fbc_s.curves] ≈ [8.0, 2pi * 0.3] rtol = 1e-12
end

@testset "SinaiBilliard get_boundary_curves" begin
    # holes plus a non-empty registry: the outer curves live in the fundamental frame, the hole curves in the full frame
    @test_throws ArgumentError get_boundary_curves(SinaiBilliard(1.0, 0.3))
    # unreduced construction: four sides then the circle
    @test length(get_boundary_curves(SinaiBilliard(1.0, 0.3; center = [0.5, -0.2]))) == 5
end

@testset "SinaiBilliard keywords" begin
    # the square and the obstacle share an arbitrary center
    sin_b = SinaiBilliard(1.0, 0.3; center = [0.5, -0.2])
    dom = sin_b.fundamental_domain
    inner = dom.holes[1].domain.boundary[1]
    @test inner.center ≈ SVector(0.5, -0.2) atol = 1e-14
    @test inner.radius == 0.3
    @test area(sin_b) ≈ 4 - pi * 0.09 rtol = 1e-12
    @test is_inside(sin_b, SVector(0.5, -0.2)) === false
    @test is_inside(sin_b, SVector(1.2, -0.2)) === true
    @test is_inside(sin_b, SVector(0.0, 0.0)) === true
    @test is_inside(sin_b, SVector(2.0, 0.0)) === false
    # the explicit origin is the symmetric default
    @test SinaiBilliard(1.0, 0.3; center = [0.0, 0.0]) isa SinaiBilliard{Float64}
    @test length(SinaiBilliard(1.0, 0.3; center = [0.0, 0.0]).symmetries) == 7
end

@testset "SinaiBilliard argument guards" begin
    # a > 0
    @test_throws ArgumentError SinaiBilliard(0.0, 0.3)
    @test_throws ArgumentError SinaiBilliard(-1.0, 0.3)
    # R > 0
    @test_throws ArgumentError SinaiBilliard(1.0, 0.0)
    @test_throws ArgumentError SinaiBilliard(1.0, -0.1)
    # R < a
    @test_throws ArgumentError SinaiBilliard(1.0, 1.0)
    @test_throws ArgumentError SinaiBilliard(1.0, 1.5)
end

@testset "SinaiBilliard eltype promotion" begin
    sin32 = SinaiBilliard(1.0f0, 0.3f0)
    @test sin32 isa SinaiBilliard{Float32}
    @test sin32.fundamental_domain isa MultiplyConnectedDomain{Float32}
    @test sin32.fundamental_domain.subdomains[1].boundary[1] isa LineSegment{Float32}
    @test sin32.fundamental_domain.holes[1] isa SymmetricHole{Float32}
    @test sin32.fundamental_domain.holes[1].piece.boundary[1] isa CircleSegment{Float32}
    @test eltype(sin32.fundamental_domain.subdomains[1].corners) == SVector{2,Float32}
    @test all(c -> c isa Union{LineSegment{Float32},CircleSegment{Float32}}, get_all_curves(sin32))
    @test genus(sin32) == 1
    @test area(sin32) isa Float32
    @test fundamental_area(sin32) isa Float32
    shifted32 = SinaiBilliard(1.0f0, 0.3f0; center = Float32[0.5, -0.2])
    @test shifted32.fundamental_domain.subdomains[1].boundary[1] isa LineSegment{Float32}
    @test shifted32.fundamental_domain.holes[1].domain.boundary[1] isa CircleSegment{Float32}
    @test area(shifted32) isa Float32
end
