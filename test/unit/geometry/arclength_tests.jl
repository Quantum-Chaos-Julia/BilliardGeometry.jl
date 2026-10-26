# Unit tests for src/geometry/arclength.jl: arc_length(crv, t), arc_length(crv, ts), arc_length(crv, pt) and construct_arc_length_interpolation.
# Reference values are closed forms (L t for a line, R |a| t for a circular arc and a polar circle) and an independent QuadGK integral of sqrt(r^2 + r'^2) for a polar curve; helpers use the prefix `al_`.

using StaticArrays: SVector
using QuadGK: quadgk

# Absolute tolerance for the closed-form lengths (worst observed 2.2e-16, margin 450x).
const al_ATOL = 1e-13
# Absolute tolerance against the independent polar integral and between the generic point and parameter forms (worst observed 7.3e-15, margin 13x).
const al_ATOL_POLAR = 1e-13
# Absolute tolerance for the cubic spline interpolation of the polar arc with the default 100 samples (worst observed 2.9e-9, margin 17x).
const al_ATOL_SPLINE = 5e-8
# Absolute tolerance for the spline interpolation of curves whose arc length is linear in t (worst observed 4.4e-16, margin 220x).
const al_ATOL_SPLINE_LINEAR = 1e-13

al_line() = LineSegment([0.2, -0.1], [1.4, 0.9])
al_arc() = CircleSegment(1.2, 1.3, 0.4, [0.1, -0.2])
# r(phi) = 1 + 0.2 cos(3 phi) on [0.3, 1.7]
al_polar() = PolarSegment(φ -> 1.0 + 0.2 * cos(3φ); arc_angle = 1.4, shift_angle = 0.3, center = [0.1, 0.2])
al_fourier() = FourierCoeffPolarSegment([0.0, 0.0, 0.0, 0.3]; arc_angle = 1.3, shift_angle = 0.4, center = [-0.2, 0.1])
# independent reference: integral of sqrt(r^2 + r'^2) dphi from the start of al_polar() to parameter t
al_polar_length(t) = quadgk(φ -> hypot(1.0 + 0.2 * cos(3φ), -0.6 * sin(3φ)), 0.3, 0.3 + 1.4 * t; rtol = 1e-13)[1]

@testset "arc_length (parameter)" begin
    # line: L t
    line = al_line()
    for t in (0.0, 0.25, 0.5, 1.0)
        @test arc_length(line, t) ≈ t * line.length atol = al_ATOL
    end
    @test line.length ≈ hypot(1.2, 1.0) atol = al_ATOL
    @test (@inferred arc_length(line, 0.5)) isa Float64

    # circular arc: R |a| t
    arc = al_arc()
    for t in (0.0, 0.2, 0.5, 1.0)
        @test arc_length(arc, t) ≈ 1.2 * 1.3 * t atol = al_ATOL
    end
    @test (@inferred arc_length(arc, 0.5)) isa Float64

    # polar circle r = 2 on an angle of 1.5: 2 * 1.5 * t
    pc = PolarSegment(φ -> 2.0; arc_angle = 1.5, shift_angle = 0.2, center = [0.0, 0.0])
    @test arc_length(pc, 0.4) ≈ 2.0 * 1.5 * 0.4 atol = al_ATOL
    @test arc_length(pc, 1.0) ≈ pc.length atol = al_ATOL

    # generic polar curve against the independent integral, and the full parameter range gives the stored length
    pol = al_polar()
    for t in (0.2, 0.5, 0.9)
        @test arc_length(pol, t) ≈ al_polar_length(t) atol = al_ATOL_POLAR
    end
    @test arc_length(pol, 1.0) ≈ pol.length atol = al_ATOL_POLAR
    @test (@inferred arc_length(pol, 0.5)) isa Float64

    # Fourier polar curve: the full range gives the stored length
    @test arc_length(al_fourier(), 1.0) ≈ al_fourier().length atol = al_ATOL_POLAR

    # additivity: s(0.9) - s(0.5) is the length of the middle part, positive and increasing
    @test arc_length(pol, 0.9) > arc_length(pol, 0.5) > arc_length(pol, 0.2) > 0

    # the tolerances are forwarded to the quadrature
    @test arc_length(line, 0.3; rtol = 1e-12) ≈ 0.3 * line.length atol = al_ATOL
    @test arc_length(line, 0.3; rtol = 1e-12, atol = 1e-14) ≈ 0.3 * line.length atol = al_ATOL

    # Float32 line of length 5: the result stays Float32
    l32 = LineSegment(Float32[0, 0], Float32[3, 4])
    @test arc_length(l32, 0.5f0) isa Float32
    @test arc_length(l32, 0.5f0) ≈ 2.5f0 atol = 1f-6
end

@testset "arc_length (parameter array)" begin
    line = al_line()
    ts = [0.25, 0.5, 1.0]
    s = arc_length(line, ts)
    @test s isa Vector{Float64}
    @test s ≈ ts .* line.length atol = al_ATOL
    # same values as the scalar method
    @test s == [arc_length(line, t) for t in ts]

    # the shape of the input is kept
    @test size(arc_length(al_arc(), [0.1 0.2; 0.3 0.4])) == (2, 2)
    @test arc_length(al_arc(), [0.1 0.2; 0.3 0.4]) ≈ 1.2 * 1.3 .* [0.1 0.2; 0.3 0.4] atol = al_ATOL

    # generic polar curve
    pol = al_polar()
    @test arc_length(pol, [0.2, 0.5, 0.9]) ≈ al_polar_length.([0.2, 0.5, 0.9]) atol = al_ATOL_POLAR
    @test arc_length(pol, Float64[]) == Float64[]
end

@testset "arc_length (point)" begin
    # line: the distance from the start point, so the point at parameter t has arc length t L
    line = al_line()
    @test arc_length(line, SVector(0.2, -0.1)) == 0.0
    @test arc_length(line, curve(line, 0.3)) ≈ 0.3 * line.length atol = al_ATOL
    @test arc_length(line, curve(line, 1.0)) ≈ line.length atol = al_ATOL
    @test (@inferred arc_length(line, SVector(0.8, 0.4))) isa Float64

    # circular arc: R |a| t
    arc = al_arc()
    @test arc_length(arc, curve(arc, 0.6)) ≈ 0.6 * arc.length atol = al_ATOL
    @test arc_length(arc, curve(arc, 0.1)) ≈ 0.1 * arc.length atol = al_ATOL
    @test (@inferred arc_length(arc, curve(arc, 0.6))) isa Float64

    # generic polar curve (through the inversion of the parametrization): the point at parameter t has the length of the parameter form
    pol = al_polar()
    for t in (0.2, 0.5, 0.9)
        @test arc_length(pol, curve(pol, t)) ≈ al_polar_length(t) atol = al_ATOL_POLAR
        @test arc_length(pol, curve(pol, t)) ≈ arc_length(pol, t) atol = al_ATOL_POLAR
    end
    pc = PolarSegment(φ -> 2.0; arc_angle = 1.5, shift_angle = 0.2, center = [0.0, 0.0])
    @test arc_length(pc, curve(pc, 0.4)) ≈ 2.0 * 1.5 * 0.4 atol = al_ATOL_POLAR
    fc = al_fourier()
    @test arc_length(fc, curve(fc, 0.7)) ≈ arc_length(fc, 0.7) atol = al_ATOL_POLAR
end

@testset "construct_arc_length_interpolation" begin
    # a line and a circular arc have arc length linear in t: the spline reproduces it, and its inverse
    for crv in (al_line(), al_arc())
        s_of_t, t_of_s = construct_arc_length_interpolation(crv)
        @test s_of_t(0.0) == 0.0
        @test s_of_t(1.0) ≈ crv.length atol = al_ATOL_SPLINE_LINEAR
        for t in (0.13, 0.55, 0.9)
            @test s_of_t(t) ≈ t * crv.length atol = al_ATOL_SPLINE_LINEAR
            @test t_of_s(t * crv.length) ≈ t atol = al_ATOL_SPLINE_LINEAR
        end
    end

    # generic polar curve: s(t) against the independent integral, and t(s(t)) = t
    pol = al_polar()
    s_of_t, t_of_s = construct_arc_length_interpolation(pol)
    @test s_of_t(0.0) == 0.0
    @test s_of_t(1.0) ≈ pol.length atol = al_ATOL_SPLINE
    for t in (0.13, 0.37, 0.5, 0.81)
        @test s_of_t(t) ≈ al_polar_length(t) atol = al_ATOL_SPLINE
        @test t_of_s(al_polar_length(t)) ≈ t atol = al_ATOL_SPLINE
    end

    # Fourier polar curve
    fc = al_fourier()
    s_of_t, t_of_s = construct_arc_length_interpolation(fc)
    @test s_of_t(1.0) ≈ fc.length atol = al_ATOL_SPLINE
    @test t_of_s(s_of_t(0.37)) ≈ 0.37 atol = al_ATOL_SPLINE

    # fewer samples and a tighter quadrature tolerance still reproduce a linear arc length
    s_of_t, t_of_s = construct_arc_length_interpolation(al_line(); n_samples = 11, rtol = 1e-12)
    @test s_of_t(0.55) ≈ 0.55 * al_line().length atol = al_ATOL_SPLINE_LINEAR
    @test t_of_s(0.55 * al_line().length) ≈ 0.55 atol = al_ATOL_SPLINE_LINEAR

    # another interpolation method: linear interpolation is exact for a linear arc length
    s_lin, t_lin = construct_arc_length_interpolation(al_arc(); interp_method = BilliardGeometry.LinearInterpolation)
    @test s_lin isa BilliardGeometry.LinearInterpolation
    @test s_lin(0.55) ≈ 0.55 * al_arc().length atol = al_ATOL_SPLINE_LINEAR
    @test t_lin(0.55 * al_arc().length) ≈ 0.55 atol = al_ATOL_SPLINE_LINEAR
end
