# Integration tests for the full boundary reconstructed from the fundamental domain and its symmetry group:
#   full_boundary, component_lengths, curve, apply_symmetry, get_boundary_curves of real billiards.
# The reconstructed boundary must be a connected counterclockwise loop whose length and enclosed area agree with closed
# forms: the stadium of half width hw (hw = 1/2 here) has perimeter 2 pi + 4 hw and area pi + 4 hw, and the C3 billiard with
# r(phi) = (1 + a (cos 3 phi - sin 6 phi)) / 2 has area (pi / 4)(1 + a^2) (the radial integral of r^2 / 2, with cos 3 phi and
# sin 6 phi orthogonal of squared norm pi each) and perimeter the integral of sqrt(r^2 + r'^2), evaluated here by an
# independent adaptive quadrature. The area is the shoelace area of 500 points per curve, an O(1/n^2) estimate, hence its own
# tolerance. The node permutation check samples the stadium boundary at the periodic midpoint grid (uniform_boundary_samples)
# and compares the geometric image of every node with the node assigned to it by the integer reflection maps. The polar
# D2 billiard has a full-circle fundamental curve, so only the image identity of its reconstructed curves is checked (its
# loop is not a proper sector reconstruction). Tolerances and reference helpers come from test/integration/common.jl.
# No eigen solves, no pinned wavefunction values. Helper functions and constants use the prefix `fc_`.

using QuadGK: quadgk
using BilliardGeometry: _apply_symmetry_to_curve

# Half width of the stadium and shape parameter of the C3 billiard.
const fc_HW = 0.5
const fc_A = 0.2

# Number of midpoint nodes of the index-permutation check (a multiple of 4).
const fc_N = 40

# Parameter values at which the polar images are compared.
const fc_TS = (0.0, 0.2, 0.5, 0.8, 1.0)

# Index maps of the reflections on the periodic midpoint grid of N nodes: X is q -> N + 1 - q, Y is q -> N/2 + 1 - q (cyclic).
fc_idx_reflect_x(q, N) = mod1(N - q + 1, N)
fc_idx_reflect_y(q, N) = mod1(N ÷ 2 - q + 1, N)

# Radial function of the C3 billiard and its derivative.
fc_c3_r(phi) = 0.5 * (1 + fc_A * (cos(3phi) - sin(6phi)))
fc_c3_dr(phi) = 0.5 * fc_A * (-3sin(3phi) - 6cos(6phi))

@testset "full_boundary closure (StadiumBilliard)" begin
    stad = StadiumBilliard(fc_HW)
    fund = get_boundary_curves(stad)
    fb = full_boundary(stad)

    # four copies of the fundamental curves form a connected loop of positive-length curves
    @test length(fb) == 4 * length(fund)
    @test all(g -> g < CLOSURE_ATOL, closure_gaps(fb))
    @test all(c -> c.length > 0, fb)

    # counterclockwise traversal enclosing the area pi + 4 hw
    @test boundary_signed_area(fb) > 0
    @test isapprox(boundary_signed_area(fb), pi + 4 * fc_HW; atol = AREA_ATOL_SHOELACE)

    # total length is the perimeter 2 pi + 4 hw and four times the fundamental physical length
    _, _, Ltot = component_lengths(fb)
    @test isapprox(Ltot, 2pi + 4 * fc_HW; atol = PERIMETER_ATOL)
    @test isapprox(Ltot, 4 * sum(c -> c.length, fund); atol = PERIMETER_ATOL)
end

@testset "full_boundary index permutation (StadiumBilliard)" begin
    fb = full_boundary(StadiumBilliard(fc_HW))
    xy = uniform_boundary_samples(fb, fc_N)
    # the geometric image of every midpoint node is the node given by the integer reflection map
    @test all(q -> norm(apply_symmetry(XAxisReflection(), xy[q]) - xy[fc_idx_reflect_x(q, fc_N)]) < ORBIT_ATOL, 1:fc_N)
    @test all(q -> norm(apply_symmetry(YAxisReflection(), xy[q]) - xy[fc_idx_reflect_y(q, fc_N)]) < ORBIT_ATOL, 1:fc_N)
end

@testset "full_boundary closure (C3Billiard)" begin
    c3 = C3Billiard(fc_A)
    fund = get_boundary_curves(c3)
    fb = full_boundary(c3)

    # three rotated copies of the single fundamental arc form a connected loop
    @test length(fb) == 3 * length(fund)
    @test all(g -> g < CLOSURE_ATOL, closure_gaps(fb))
    @test all(c -> c.length > 0, fb)

    # counterclockwise traversal enclosing the area (pi / 4)(1 + a^2)
    @test boundary_signed_area(fb) > 0
    @test isapprox(boundary_signed_area(fb), pi / 4 * (1 + fc_A^2); atol = AREA_ATOL_SHOELACE)

    # total length is the integral of sqrt(r^2 + r'^2) over the full turn and three times the fundamental length
    _, _, Ltot = component_lengths(fb)
    L_ref, _ = quadgk(phi -> sqrt(fc_c3_r(phi)^2 + fc_c3_dr(phi)^2), 0.0, 2pi; rtol = 1e-14)
    @test isapprox(Ltot, L_ref; atol = PERIMETER_ATOL)
    @test isapprox(Ltot, 3 * sum(c -> c.length, fund); atol = PERIMETER_ATOL)
end

@testset "full_boundary polar images (D2 PolarBilliard)" begin
    polar = PolarBilliard([0.0, 0.0, 0.0, 0.3])
    polar_d2 = PolarBilliard{Float64}(polar.fundamental_domain, register_symmetries(YAxisReflection(), XYAxisReflection(),
                                                                                    XAxisReflection()))
    fb = full_boundary(polar_d2)
    @test length(fb) == 4
    @test all(c -> c isa FourierCoeffPolarSegment, fb)

    # an image curve evaluates to the symmetry applied to the point of the source curve
    src = polar.fundamental_domain.boundary[1]
    for sym in (XAxisReflection(), YAxisReflection(), XYAxisReflection())
        g = _apply_symmetry_to_curve(sym, src)
        @test g isa FourierCoeffPolarSegment
        @test all(t -> norm(curve(g, t) - apply_symmetry(sym, curve(src, t))) < POLAR_IMAGE_ATOL, fc_TS)
    end
end
