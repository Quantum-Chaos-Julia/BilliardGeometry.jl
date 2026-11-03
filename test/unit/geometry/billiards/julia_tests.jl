# Unit tests for the Julia-logo billiard constructor and structure:
#   src/geometry/billiards/julia.jl   (JuliaBilliard{T}, JuliaBilliard(r, d; R))
# Reference object: disc of radius R = 1 with three obstacles of radius r at distance d from the origin (polar angles pi/2, 7 pi/6,
#   11 pi/6), symmetry group D3 of order 6 = Dn_symmetry(3; shift = 1//6) (5 non-identity elements, sym_ids 1..5, mirrors at pi/2, pi/6,
#   5 pi/6). The outer region is the sector pi/6 <= phi <= pi/2: one physical arc (domain id 1) and SymmetryWall(1, 2), SymmetryWall(5, 2).
#   The obstacles are one SymmetricHole whose piece is the half disc about (0, d) with x >= 0 (domain id 2, arc + mirror SymmetryWall(1, 2)),
#   6 pieces with ids 2, 102, ..., 502 forming 3 closed loops of 2 arcs each. Closed forms: area pi R^2 - 3 pi r^2, loop lengths
#   [2 pi R, 3 x 2 pi r], genus 3, blocks [6, 6]. Guards: r, d, R > 0, d + r < R, sqrt(3) d > 2 r. No eigen solves.

using StaticArrays: SVector

const jb_ATOL = 1e-12

@testset "JuliaBilliard fields" begin
    jb = JuliaBilliard(0.2, 0.5)
    @test jb isa JuliaBilliard{Float64}
    @test jb isa AbsBilliard
    @test jb.fundamental_domain isa MultiplyConnectedDomain{Float64}
    @test jb.symmetries isa SymmetryRegistry
    @test fieldnames(JuliaBilliard) == (:fundamental_domain, :symmetries)
    @test @inferred(JuliaBilliard(0.2, 0.5)) isa JuliaBilliard{Float64}
end

@testset "JuliaBilliard symmetries" begin
    jb = JuliaBilliard(0.2, 0.5)
    reg = jb.symmetries
    @test length(reg) == 5
    @test [s.sym_id for s in reg] == collect(1:5)
    @test typeof.(collect(reg)) == typeof.(collect(Dn_symmetry(3; shift = 1 // 6)))
    # mirrors at pi/2 (sym_id 1) and pi/6 (sym_id 5) bound the wedge; rotations by 2 pi/3 and 4 pi/3 sit at even ids
    @test reg[1] isa YAxisReflection
    @test reg[5] isa AxisReflection{1 // 6}
    @test reg[2] isa NFoldRotation{Float64}
    @test reg[2].m == 1
    @test reg[4].m == 2
    @test length(get_symmetry_group(jb)) == 5
    @test typeof.(get_symmetries(jb)) == (YAxisReflection, AxisReflection{5 // 6})
    @test symmetry_compatible(jb.fundamental_domain.holes[1], reg)
end

@testset "JuliaBilliard outer sector" begin
    dom = JuliaBilliard(0.2, 0.5).fundamental_domain
    @test length(dom.subdomains) == 1
    outer = dom.subdomains[1]
    @test outer.id == 1
    @test length(outer.boundary) == 3
    arc, wall1, wall0 = outer.boundary
    @test arc isa CircleSegment{Float64}
    @test wall1 isa LineSegment{Float64}
    @test wall0 isa LineSegment{Float64}
    @test [c.bc for c in outer.boundary] == [SpecularReflection(), SymmetryWall(1, 2), SymmetryWall(5, 2)]
    @test arc.radius == 1.0
    @test arc.arc_angle ≈ pi / 3 atol = jb_ATOL
    @test outer.corners[1] ≈ SVector(cos(pi / 6), sin(pi / 6)) atol = jb_ATOL
    @test outer.corners[2] ≈ SVector(0.0, 1.0) atol = jb_ATOL
    @test outer.corners[3] == SVector(0.0, 0.0)
end

@testset "JuliaBilliard obstacles (SymmetricHole)" begin
    jb = JuliaBilliard(0.2, 0.5)
    dom = jb.fundamental_domain
    @test length(dom.holes) == 1
    h = dom.holes[1]
    @test h isa SymmetricHole{Float64}
    @test h.center == SVector(0.0, 0.0)
    @test length(h.symmetries) == 5
    # half disc piece plus its 5 images: three obstacles cut along the mirror lines
    @test length(h.pieces) == 6
    @test [p.id for p in h.pieces] == [2, 102, 202, 302, 402, 502]
    piece = h.pieces[1]
    @test length(piece.boundary) == 2
    @test piece.boundary[1] isa CircleSegment{Float64}
    @test piece.boundary[1].radius == 0.2
    @test piece.boundary[1].center == SVector(0.0, 0.5)
    @test piece.boundary[2].bc == SymmetryWall(1, 2)
    # three closed loops of two arcs each, hole on the left, one around each obstacle centre
    @test hole_loop_count(h) == 3
    comps = boundary_components(h)
    @test [length(l) for l in comps] == [2, 2, 2]
    for loop in comps
        @test all(c -> c isa CircleSegment{Float64} && c.radius == 0.2, loop)
        @test curve(loop[1], 1.0) ≈ curve(loop[2], 0.0) atol = jb_ATOL
        @test curve(loop[2], 1.0) ≈ curve(loop[1], 0.0) atol = jb_ATOL
        @test sum(area, loop) ≈ pi * 0.04 atol = jb_ATOL
    end
    centers = [(curve(l[1], 0.0) + curve(l[1], 1.0)) / 2 for l in comps]
    expected = [0.5 * SVector(cos(a), sin(a)) for a in (pi / 2, 7pi / 6, 11pi / 6)]
    @test all(c -> any(e -> isapprox(c, e; atol = jb_ATOL), expected), centers)
    @test area(h) ≈ 3 * pi * 0.04 atol = jb_ATOL
    @test is_inside(h, SVector(0.05, 0.5)) === true
    @test is_inside(h, SVector(0.5, 0.0)) === false
    # the obstacle group is the billiard group
    @test symmetry_compatible(h, jb.symmetries)
end

@testset "JuliaBilliard area, genus and full boundary" begin
    jb = JuliaBilliard(0.2, 0.5)
    @test genus(jb.fundamental_domain) == 3
    @test area(jb) ≈ pi - 3 * pi * 0.04 atol = jb_ATOL
    @test area(JuliaBilliard(0.15, 0.4)) ≈ pi - 3 * pi * 0.15^2 atol = jb_ATOL
    @test area(JuliaBilliard(0.2, 0.5; R = 2.0)) ≈ 4pi - 3 * pi * 0.04 atol = jb_ATOL
    fb = full_boundary_components(jb)
    @test fb.is_hole == [false, true, true, true]
    @test [length(c) for c in fb.curves] == [6, 2, 2, 2]
    @test [sum(c.length for c in loop) for loop in fb.curves] ≈ [2pi, 0.4pi, 0.4pi, 0.4pi] atol = jb_ATOL
    @test all(b -> length(b) == 6, fb.blocks)
    @test length(fb.blocks) == 2
end

@testset "JuliaBilliard membership" begin
    dom = JuliaBilliard(0.2, 0.5).fundamental_domain
    @test is_inside(dom, SVector(0.2, 0.8)) === true
    @test is_inside(dom, SVector(0.1, 0.5)) === false       # inside the obstacle
    @test is_inside(dom, SVector(0.8, 0.8)) === false       # outside the disc
    @test is_inside(dom, SVector(0.9, 0.0)) === false       # outside the sector
end

@testset "JuliaBilliard argument guards" begin
    @test_throws ArgumentError JuliaBilliard(-0.1, 0.5)
    @test_throws ArgumentError JuliaBilliard(0.0, 0.5)
    @test_throws ArgumentError JuliaBilliard(0.2, -0.5)
    @test_throws ArgumentError JuliaBilliard(0.2, 0.5; R = -1.0)
    # obstacle must lie strictly inside the disc
    @test_throws ArgumentError JuliaBilliard(0.2, 0.9)
    @test_throws ArgumentError JuliaBilliard(0.2, 0.8)
    # obstacles must not touch: sqrt(3) d > 2 r
    @test_throws ArgumentError JuliaBilliard(0.2, 0.2)
end

@testset "JuliaBilliard eltype promotion" begin
    jb32 = JuliaBilliard(0.2f0, 0.5f0)
    @test jb32 isa JuliaBilliard{Float32}
    @test jb32.fundamental_domain isa MultiplyConnectedDomain{Float32}
    @test length(jb32.symmetries) == 5
end
