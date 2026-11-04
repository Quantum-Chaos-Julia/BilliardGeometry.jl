# Unit tests for the C6C3HoleBilliard constructor and structure:
#   src/geometry/billiards/c6_c3_hole.jl   (C6C3HoleBilliard{T}, C6C3HoleBilliard(a, b; R, R_hole))
# Reference object: outer boundary r = R (1 + a cos 6 phi) (C6), hole boundary r = R_hole (1 + b cos 3 phi) (C3). The billiard symmetry
# group is C6 ∩ C3 = C3, a proper subgroup of the outer symmetry group (the constructor logs this with @info). The fundamental domain is
# a MultiplyConnectedDomain: the outer sector 0 <= phi <= 2 pi / 3 (one arc, two SymmetryWalls, domain id 1) and one SymmetricHole whose
# piece is the matching hole sector (domain id 2) with group Cn_symmetry(3) about the origin; registry Cn_symmetry(3).
# Closed forms: area = pi R^2 (1 + a^2 / 2) - pi R_hole^2 (1 + b^2 / 2) (polar area (1/2) ∫ r^2 dphi; cross terms vanish).
# Guards: R > 0, R_hole > 0, |a| < 1, |b| < 1, hole inside the outer boundary. No eigen solves.

c6c3_area(a, b, R, Rh) = pi * R^2 * (1 + a^2 / 2) - pi * Rh^2 * (1 + b^2 / 2)

@testset "C6C3HoleBilliard fields" begin
    b6 = C6C3HoleBilliard(0.1, 0.15)
    @test b6 isa C6C3HoleBilliard{Float64}
    @test b6 isa AbsBilliard
    @test b6.fundamental_domain isa MultiplyConnectedDomain{Float64}
    @test b6.symmetries isa SymmetryRegistry
    @test fieldnames(C6C3HoleBilliard) == (:fundamental_domain, :symmetries)
end

@testset "C6C3HoleBilliard symmetries" begin
    b6 = C6C3HoleBilliard(0.1, 0.15)
    # rotations by 2 pi / 3 and 4 pi / 3 (non-identity elements of C3)
    @test length(b6.symmetries) == 2
    @test [g.sym_id for g in b6.symmetries] == [g.sym_id for g in Cn_symmetry(3)]
    @test length(get_symmetry_group(b6)) == 2
    @test_logs C6C3HoleBilliard(0.1, 0.15)
end

@testset "C6C3HoleBilliard curves" begin
    dom = C6C3HoleBilliard(0.1, 0.15).fundamental_domain
    @test dom.id == 1
    @test length(dom.subdomains) == 1
    @test length(dom.holes) == 1
    @test dom.holes[1] isa SymmetricHole{Float64}
    outer_dom = dom.subdomains[1]
    @test outer_dom isa SimpleDomain{Float64}
    @test outer_dom.id == 1
    @test length(outer_dom.boundary) == 3
    arc, wall1, wall0 = outer_dom.boundary
    @test arc.bc isa SpecularReflection
    @test wall1.bc isa SymmetryWall
    @test wall0.bc isa SymmetryWall
    @test [c.domain_id for c in outer_dom.boundary] == [1, 1, 1]
    @test [c.segment_id for c in outer_dom.boundary] == [1, 2, 3]
    # walls are radii through the origin, closing the sector
    @test wall1.pt1 == SVector(0.0, 0.0)
    @test wall0.pt0 == SVector(0.0, 0.0)
    # outer sector corners: radius R (1 + a) at phi = 0 (cos 0 = 1), R (1 + a cos 4 pi) = R (1 + a) at phi = 2 pi / 3
    @test hypot(outer_dom.corners[1]...) ≈ 1.1 atol = 1e-12
    @test hypot(outer_dom.corners[2]...) ≈ 1.1 atol = 1e-12
    @test outer_dom.corners[1] ≈ SVector(1.1, 0.0) atol = 1e-12
    @test outer_dom.corners[2] ≈ 1.1 * SVector(cos(2pi / 3), sin(2pi / 3)) atol = 1e-12
    @test outer_dom.corners[3] == SVector(0.0, 0.0)

    hole = dom.holes[1]
    @test hole.center == SVector(0.0, 0.0)
    @test length(hole.symmetries) == 2
    @test length(hole.pieces) == 3
    @test hole.piece.id == 2
    @test count(c -> c.bc isa SymmetryWall, hole.piece.boundary) == 2
    # hole sector corners: radius R_hole (1 + b) = (1/3)(1.15) at phi = 0 and 2 pi / 3 (cos 2 pi = 1)
    @test hole.piece.corners[1] ≈ (1 / 3) * 1.15 * SVector(1.0, 0.0) atol = 1e-12
    @test hole.piece.corners[2] ≈ (1 / 3) * 1.15 * SVector(cos(2pi / 3), sin(2pi / 3)) atol = 1e-12
end

@testset "C6C3HoleBilliard boundary components" begin
    b6 = C6C3HoleBilliard(0.1, 0.15)
    @test length(boundary_components(b6)) == 2
    # three rotated copies of (arc) outer and hole, each plus ... : 6 curves in total
    @test length(full_boundary(b6)) == 6
end

@testset "C6C3HoleBilliard area" begin
    @test area(C6C3HoleBilliard(0.1, 0.15)) ≈ c6c3_area(0.1, 0.15, 1.0, 1 / 3) rtol = 1e-10
    @test area(C6C3HoleBilliard(0.2, -0.3; R = 2.0, R_hole = 0.5)) ≈ c6c3_area(0.2, -0.3, 2.0, 0.5) rtol = 1e-10
    # a = b = 0: annulus
    @test area(C6C3HoleBilliard(0.0, 0.0)) ≈ pi * (1 - 1 / 9) rtol = 1e-10
    # one third of the physical area lies in the fundamental domain
    b6 = C6C3HoleBilliard(0.1, 0.15)
    @test fundamental_area(b6) ≈ area(b6) / 3 rtol = 1e-10
end

@testset "C6C3HoleBilliard Float32" begin
    b32 = C6C3HoleBilliard(0.1f0, 0.15f0)
    @test b32 isa C6C3HoleBilliard{Float32}
    @test b32.fundamental_domain isa MultiplyConnectedDomain{Float32}
    @test area(b32) ≈ c6c3_area(0.1f0, 0.15f0, 1.0f0, 1.0f0 / 3) rtol = 1e-4
end

@testset "C6C3HoleBilliard argument guards" begin
    @test_throws ArgumentError C6C3HoleBilliard(0.1, 0.15; R = 0.0)
    @test_throws ArgumentError C6C3HoleBilliard(0.1, 0.15; R = -1.0)
    @test_throws ArgumentError C6C3HoleBilliard(0.1, 0.15; R_hole = 0.0)
    @test_throws ArgumentError C6C3HoleBilliard(0.1, 0.15; R_hole = -0.2)
    @test_throws ArgumentError C6C3HoleBilliard(1.0, 0.15)
    @test_throws ArgumentError C6C3HoleBilliard(-1.5, 0.15)
    @test_throws ArgumentError C6C3HoleBilliard(0.1, 1.0)
    @test_throws ArgumentError C6C3HoleBilliard(0.1, -1.2)
    # hole not inside the outer boundary
    @test_throws ArgumentError C6C3HoleBilliard(0.1, 0.15; R_hole = 0.9)
end
