# Unit tests for the square-within-square billiard constructor and structure:
#   src/geometry/billiards/square_in_square.jl   (SquareWithinSquareBilliard{T}, SquareWithinSquareBilliard(a; positions, sizes))
# Reference object: a MultiplyConnectedDomain with id 1, genus n and an empty symmetry registry, made of an outer SimpleDomain
# (the four counterclockwise LineSegments of [-a, a]^2, domain id 1) and n DomainHoles, one per obstacle square of side
# `sizes[i]` centred at (d_i, d_i) with the positions sorted ascending. Each hole is a SimpleDomain of four counterclockwise
# LineSegments (bottom, right, top, left; segment ids 1..4, orientation +1) with domain id i + 1, i counted from the lower left.
# Closed forms: area 4 a^2 - sum(sizes^2), membership, corner angles (pi/2 outer, 3 pi/2 at each hole corner).
# Guards: a > 0, equal lengths of positions and sizes, positive sizes, holes strictly inside, no overlap. No eigen solves.

using StaticArrays: SVector

@testset "SquareWithinSquareBilliard fields" begin
    sws = SquareWithinSquareBilliard(1.0)
    @test sws isa SquareWithinSquareBilliard{Float64}
    @test sws isa AbsBilliard
    @test sws.fundamental_domain isa MultiplyConnectedDomain{Float64}
    @test sws.symmetries isa SymmetryRegistry
    @test fieldnames(SquareWithinSquareBilliard) == (:fundamental_domain, :symmetries)
    @test @inferred(SquareWithinSquareBilliard(1.0)) isa SquareWithinSquareBilliard{Float64}
end

@testset "SquareWithinSquareBilliard symmetries" begin
    sws = SquareWithinSquareBilliard(1.0)
    @test isempty(sws.symmetries)
    @test get_symmetries(sws) == ()
end

@testset "SquareWithinSquareBilliard curves" begin
    dom = SquareWithinSquareBilliard(1.0).fundamental_domain
    @test dom.id == 1
    # one outer subdomain and two holes by default, every domain a closed loop of four lines
    @test length(dom.subdomains) == 1
    @test length(dom.holes) == 2
    @test all(h -> h isa DomainHole, dom.holes)
    doms = [dom.subdomains[1], dom.holes[1].domain, dom.holes[2].domain]
    @test all(d -> d isa SimpleDomain{Float64}, doms)
    @test [d.id for d in doms] == [1, 2, 3]
    for d in doms
        @test length(d.boundary) == 4
        @test all(c -> c isa LineSegment{Float64,SpecularReflection}, d.boundary)
        # counterclockwise, one domain id per loop, segment ids bottom, right, top, left
        @test all(c -> c.orientation == 1, d.boundary)
        @test all(c -> c.domain_id == d.id, d.boundary)
        @test [c.segment_id for c in d.boundary] == [1, 2, 3, 4]
        # every curve starts where the previous one ends, and the loop closes
        @test all(i -> d.boundary[i].pt1 == d.boundary[mod1(i + 1, 4)].pt0, 1:4)
        @test d.corners == [c.pt1 for c in d.boundary]
    end

    # outer square [-1, 1]^2, holes of side 1/4 centred at (-1/2, -1/2) and (1/2, 1/2), lower left one first
    @test doms[1].corners == [SVector(1.0, -1.0), SVector(1.0, 1.0), SVector(-1.0, 1.0), SVector(-1.0, -1.0)]
    @test doms[2].boundary[1].pt0 == SVector(-0.625, -0.625)
    @test doms[2].boundary[3].pt0 == SVector(-0.375, -0.375)
    @test doms[3].boundary[1].pt0 == SVector(0.375, 0.375)
    @test doms[3].boundary[3].pt0 == SVector(0.625, 0.625)
end

@testset "SquareWithinSquareBilliard hole count" begin
    one_hole = SquareWithinSquareBilliard(1.0; positions=[0.0], sizes=[0.4]).fundamental_domain
    @test length(one_hole.holes) == 1
    @test one_hole.holes[1].domain.id == 2
    @test one_hole.holes[1].domain.boundary[1].pt0 ≈ SVector(-0.2, -0.2) atol = 1e-14

    three_holes = SquareWithinSquareBilliard(1.0; positions=[0.5, -0.5, 0.0], sizes=[0.2, 0.2, 0.2]).fundamental_domain
    @test length(three_holes.holes) == 3
    @test [h.domain.id for h in three_holes.holes] == [2, 3, 4]
    # the ids follow the sorted positions: hole 2 sits at -0.5, hole 4 at 0.5
    @test [h.domain.boundary[1].pt0[1] for h in three_holes.holes] ≈ [-0.6, -0.1, 0.4] atol = 1e-14

    # unsorted positions are sorted with their sizes: the hole at 0.5 has side 0.1 and the one at -0.5 side 0.2
    unsorted = SquareWithinSquareBilliard(1.0; positions=[0.5, -0.5], sizes=[0.1, 0.2]).fundamental_domain
    @test unsorted.holes[1].domain.boundary[1].pt0 ≈ SVector(-0.6, -0.6) atol = 1e-14
    @test unsorted.holes[2].domain.boundary[1].pt0 ≈ SVector(0.45, 0.45) atol = 1e-14
end

@testset "SquareWithinSquareBilliard topology" begin
    sws = SquareWithinSquareBilliard(1.0)
    # genus equals the number of holes
    @test genus(sws) == 2
    @test genus(sws.fundamental_domain) == 2
    @test genus(SquareWithinSquareBilliard(1.0; positions=[0.0], sizes=[0.4])) == 1
    @test genus(SquareWithinSquareBilliard(1.0; positions=[0.5, -0.5, 0.0], sizes=[0.2, 0.2, 0.2])) == 3
    # one component of four curves for the outer square and for each hole
    comps = boundary_components(sws)
    @test length(comps) == genus(sws) + 1
    @test [length(c) for c in comps] == [4, 4, 4]
    @test [[c.domain_id for c in comp] for comp in comps] == [fill(1, 4), fill(2, 4), fill(3, 4)]
    @test boundary_component_is_hole(sws) == [false, true, true]
    @test length(geometric_boundary_components(sws)) == genus(sws) + 1
    @test hole_domain_ids(sws) == [2, 3]
    @test component_of_domain_id(sws) == Dict(1 => 0, 2 => 1, 3 => 2)
    @test [is_hole_curve(sws, c) for c in get_all_curves(sws)] == [fill(false, 4); fill(true, 8)]
end

@testset "SquareWithinSquareBilliard geometry" begin
    sws = SquareWithinSquareBilliard(1.0)
    # area 4 a^2 - sum of the squared side lengths
    @test area(sws) ≈ 4 - 2 * 0.25^2 rtol = 1e-12
    @test area(SquareWithinSquareBilliard(1.0; positions=[0.0], sizes=[0.4])) ≈ 4 - 0.16 rtol = 1e-12
    @test area(SquareWithinSquareBilliard(2.0)) ≈ 16 - 2 * 0.5^2 rtol = 1e-12
    # without symmetry the fundamental area is the area
    @test fundamental_area(sws) ≈ area(sws) rtol = 1e-14
    # membership: the whole square minus the holes
    @test is_inside(sws, SVector(0.0, 0.0)) === true
    @test is_inside(sws, SVector(0.9, 0.2)) === true
    @test is_inside(sws, SVector(-0.9, 0.1)) === true
    @test is_inside(sws, SVector(0.2, 0.9)) === true
    @test is_inside(sws, SVector(0.5, 0.5)) === false        # centre of the upper hole
    @test is_inside(sws, SVector(-0.5, -0.5)) === false      # centre of the lower hole
    @test is_inside(sws, SVector(1.2, 0.1)) === false
    # corners: four right angles of the outer square, then 3 pi / 2 at each of the 8 hole corners
    @test corner_angles(sws) ≈ [fill(pi / 2, 4); fill(3pi / 2, 8)] atol = 1e-12
    @test corner_angles(sws; fundamental = false) ≈ [fill(pi / 2, 4); fill(3pi / 2, 8)] atol = 1e-12
    # full boundary: the outer loop, then one loop per hole
    @test length(full_boundary(sws)) == 12
    @test full_boundary_components(sws).is_hole == [false, true, true]
end

@testset "SquareWithinSquareBilliard argument guards" begin
    # a > 0
    @test_throws ArgumentError SquareWithinSquareBilliard(0.0)
    @test_throws ArgumentError SquareWithinSquareBilliard(-1.0)
    # positions and sizes have the same length
    @test_throws ArgumentError SquareWithinSquareBilliard(1.0; positions=[0.0], sizes=[0.2, 0.3])
    # sizes are positive
    @test_throws ArgumentError SquareWithinSquareBilliard(1.0; positions=[0.0], sizes=[0.0])
    @test_throws ArgumentError SquareWithinSquareBilliard(1.0; positions=[0.0], sizes=[-0.2])
    # holes lie strictly inside the outer square, on either side
    @test_throws ArgumentError SquareWithinSquareBilliard(1.0; positions=[0.9], sizes=[0.3])
    @test_throws ArgumentError SquareWithinSquareBilliard(1.0; positions=[-0.9], sizes=[0.3])
    @test_throws ArgumentError SquareWithinSquareBilliard(1.0; positions=[0.9], sizes=[0.2])
    # holes may neither overlap nor touch
    @test_throws ArgumentError SquareWithinSquareBilliard(1.0; positions=[-0.5, -0.4], sizes=[0.4, 0.4])
    @test_throws ArgumentError SquareWithinSquareBilliard(1.0; positions=[-0.5, -0.3], sizes=[0.2, 0.2])
end

@testset "SquareWithinSquareBilliard eltype promotion" begin
    sws32 = SquareWithinSquareBilliard(1.0f0)
    @test sws32 isa SquareWithinSquareBilliard{Float32}
    @test sws32.fundamental_domain isa MultiplyConnectedDomain{Float32}
    @test all(c -> c isa LineSegment{Float32}, get_all_curves(sws32))
    @test eltype(sws32.fundamental_domain.subdomains[1].corners) == SVector{2,Float32}
    @test genus(sws32) == 2
    @test area(sws32) isa Float32
    # Float64 positions and sizes are converted to the element type of a
    sws32 = SquareWithinSquareBilliard(1.0f0; positions=[-0.5, 0.5], sizes=[0.25, 0.25])
    @test sws32 isa SquareWithinSquareBilliard{Float32}
    @test all(c -> c isa LineSegment{Float32}, get_all_curves(sws32))
end
