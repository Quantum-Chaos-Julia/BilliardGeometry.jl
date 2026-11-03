# Unit tests for the generalized Sinai billiard:
#   src/geometry/billiards/generalized_sinai.jl   (GeneralizedSinaiBilliard{T}, GeneralizedSinaiBilliard(a, b))
# Reference object (a, b) = (0.4, 0.6): first-quadrant fundamental domain bounded by xwall (0,0)->(a,0), concave arc_x through (a,0) and
# (1,1) centered at (cx, 0), concave arc_y through (1,1) and (0,b) centered at (0, cy), and ywall (0,b)->(0,0).
# Evidence: analytic (arcs cross axes at right angles, so cx = (2 - a^2) / (2 (1 - a)); circle passes through (a,0) and (1,1); curve
# continuity), plus values observed from a run for a = 0.4, b = 0.6 (area, corner_angles). Guards: 0 < a, b < 1. No eigen solves.

using StaticArrays: SVector

const gsb_a = 0.4
const gsb_b = 0.6
const gsb_cx = (2 - gsb_a^2) / (2 * (1 - gsb_a))
const gsb_cy = (2 - gsb_b^2) / (2 * (1 - gsb_b))

@testset "GeneralizedSinaiBilliard fields" begin
    gsb = GeneralizedSinaiBilliard(gsb_a, gsb_b)
    @test gsb isa GeneralizedSinaiBilliard{Float64}
    @test gsb isa AbsBilliard
    @test gsb.fundamental_domain isa SimpleDomain{Float64}
    @test gsb.fundamental_domain.id == 1
    @test gsb.symmetries isa SymmetryRegistry
    @test fieldnames(GeneralizedSinaiBilliard) == (:fundamental_domain, :symmetries)
    @test GeneralizedSinaiBilliard(1 // 2, 1 // 4) isa GeneralizedSinaiBilliard{Float64}
end

@testset "GeneralizedSinaiBilliard curves" begin
    cs = GeneralizedSinaiBilliard(gsb_a, gsb_b).fundamental_domain.boundary
    @test length(cs) == 4
    arc_x, arc_y, ywall, xwall = cs
    @test xwall isa LineSegment
    @test arc_x isa CircleSegment
    @test arc_y isa CircleSegment
    @test ywall isa LineSegment
    @test xwall.bc isa SymmetryWall
    @test ywall.bc isa SymmetryWall
    @test arc_x.bc isa SpecularReflection
    @test arc_y.bc isa SpecularReflection
    @test [c.segment_id for c in cs] == [1, 2, 3, 4]
    @test all(c -> c.domain_id == 1, cs)
    # end-to-start continuity, closed loop
    for i in 1:4
        @test curve(cs[i], 1.0) ≈ curve(cs[mod1(i + 1, 4)], 0.0) atol = 1e-12
    end
    # end points
    @test curve(xwall, 0.0) ≈ SVector(0.0, 0.0) atol = 1e-12
    @test curve(xwall, 1.0) ≈ SVector(gsb_a, 0.0) atol = 1e-12
    @test curve(arc_x, 1.0) ≈ SVector(1.0, 1.0) atol = 1e-12
    @test curve(arc_y, 0.0) ≈ SVector(1.0, 1.0) atol = 1e-12
    @test curve(arc_y, 1.0) ≈ SVector(0.0, gsb_b) atol = 1e-12
    # arc centers on the axes, radii reach the axis crossings
    @test arc_x.center ≈ SVector(gsb_cx, 0.0)
    @test arc_y.center ≈ SVector(0.0, gsb_cy)
    @test arc_x.radius ≈ gsb_cx - gsb_a
    @test arc_y.radius ≈ gsb_cy - gsb_b
    # circle through (a,0) and (1,1)
    @test hypot(gsb_a - gsb_cx, 0.0) ≈ arc_x.radius
    @test hypot(1.0 - gsb_cx, 1.0) ≈ arc_x.radius
    @test hypot(1.0, 1.0 - gsb_cy) ≈ arc_y.radius
    @test hypot(0.0, gsb_b - gsb_cy) ≈ arc_y.radius
    @test curve(arc_x, 0.5) isa SVector{2,Float64}
    @test hypot((curve(arc_x, 0.5) - arc_x.center)...) ≈ arc_x.radius
    @test hypot((curve(arc_y, 0.5) - arc_y.center)...) ≈ arc_y.radius
    corners = GeneralizedSinaiBilliard(gsb_a, gsb_b).fundamental_domain.corners
    @test corners == [SVector(gsb_a, 0.0), SVector(1.0, 1.0), SVector(0.0, gsb_b), SVector(0.0, 0.0)]
end

@testset "GeneralizedSinaiBilliard is_inside" begin
    gsb = GeneralizedSinaiBilliard(gsb_a, gsb_b)
    @test is_inside(gsb, SVector(0.2, 0.2))
    # inside the concave circles, hence outside the billiard
    @test !is_inside(gsb, SVector(0.95, 0.95))
    @test !is_inside(gsb, SVector(0.9, 0.1))
end

@testset "GeneralizedSinaiBilliard area" begin
    # observed from a run
    @test area(GeneralizedSinaiBilliard(gsb_a, gsb_b)) ≈ 1.1900527359836985
end

@testset "GeneralizedSinaiBilliard corner_angles" begin
    # observed from a run; the three right angles are analytic, the concave corner at (1,1) is observed
    angles = corner_angles(GeneralizedSinaiBilliard(gsb_a, gsb_b))
    @test length(angles) == 4
    @test angles ≈ [π / 2, 0.27105542797100135, π / 2, π / 2]
end

@testset "GeneralizedSinaiBilliard symmetries" begin
    gsb = GeneralizedSinaiBilliard(gsb_a, gsb_b)
    @test [g.sym_id for g in gsb.symmetries] == [g.sym_id for g in D2_symmetry()]
    @test length(gsb.symmetries) == 3
end

@testset "GeneralizedSinaiBilliard argument guards" begin
    @test_throws ArgumentError GeneralizedSinaiBilliard(0.0, 0.5)
    @test_throws ArgumentError GeneralizedSinaiBilliard(1.0, 0.5)
    @test_throws ArgumentError GeneralizedSinaiBilliard(-0.1, 0.5)
    @test_throws ArgumentError GeneralizedSinaiBilliard(0.5, 0.0)
    @test_throws ArgumentError GeneralizedSinaiBilliard(0.5, 1.0)
    @test_throws ArgumentError GeneralizedSinaiBilliard(0.5, 1.5)
end

# a == b: D4 symmetry, fundamental domain is the eighth 0 <= y <= x (arc_x, diagonal, xwall).
const gsb_s = 0.5

@testset "GeneralizedSinaiBilliard a == b curves and corners" begin
    gsb = GeneralizedSinaiBilliard(gsb_s, gsb_s)
    @test gsb isa GeneralizedSinaiBilliard{Float64}
    cs = gsb.fundamental_domain.boundary
    @test length(cs) == 3
    arc_x, diagonal, xwall = cs
    @test arc_x isa CircleSegment
    @test diagonal isa LineSegment
    @test xwall isa LineSegment
    @test arc_x.bc isa SpecularReflection
    @test diagonal.bc isa SymmetryWall
    @test xwall.bc isa SymmetryWall
    @test (diagonal.bc.sym_id, diagonal.bc.sector_id) == (1, 2)
    @test (xwall.bc.sym_id, xwall.bc.sector_id) == (7, 8)
    @test [c.segment_id for c in cs] == [1, 2, 3]
    @test all(c -> c.domain_id == 1, cs)
    for i in 1:3
        @test curve(cs[i], 1.0) ≈ curve(cs[mod1(i + 1, 3)], 0.0) atol = 1e-12
    end
    @test curve(arc_x, 0.0) ≈ SVector(gsb_s, 0.0) atol = 1e-12
    @test curve(arc_x, 1.0) ≈ SVector(1.0, 1.0) atol = 1e-12
    @test curve(diagonal, 0.0) ≈ SVector(1.0, 1.0) atol = 1e-12
    @test curve(diagonal, 1.0) ≈ SVector(0.0, 0.0) atol = 1e-12
    @test curve(xwall, 0.0) ≈ SVector(0.0, 0.0) atol = 1e-12
    @test curve(xwall, 1.0) ≈ SVector(gsb_s, 0.0) atol = 1e-12
    @test gsb.fundamental_domain.corners == [SVector(gsb_s, 0.0), SVector(1.0, 1.0), SVector(0.0, 0.0)]
end

@testset "GeneralizedSinaiBilliard a == b symmetries" begin
    gsb = GeneralizedSinaiBilliard(gsb_s, gsb_s)
    @test length(gsb.symmetries) == 7
    @test [g.sym_id for g in gsb.symmetries] == [g.sym_id for g in D4_symmetry()]
end

@testset "GeneralizedSinaiBilliard a == b area and is_inside" begin
    eighth = GeneralizedSinaiBilliard(gsb_s, gsb_s)
    # area(billiard) is the full-domain area (fundamental domain times symmetry order), so it matches the nearby
    # quadrant billiard (a != b) up to O(eps); the eighth fundamental domain is half the quadrant one.
    quadrant = GeneralizedSinaiBilliard(gsb_s, gsb_s + 1e-9)
    @test area(eighth) ≈ area(quadrant) rtol = 1e-6
    @test area(eighth.fundamental_domain) ≈ area(quadrant.fundamental_domain) / 2 rtol = 1e-6
    # observed from a run
    @test area(eighth) ≈ 1.2044048874899227
    for p in (SVector(0.2, 0.1), SVector(0.95, 0.95), SVector(0.4, 0.3), SVector(0.9, 0.1))
        @test is_inside(eighth, p) == is_inside(quadrant, p)
    end
    @test is_inside(eighth, SVector(0.2, 0.1))
    @test !is_inside(eighth, SVector(0.95, 0.95))
end
