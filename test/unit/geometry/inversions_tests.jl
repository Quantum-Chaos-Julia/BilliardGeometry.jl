# Unit tests for the inversion of curve parametrizations:
#   src/geometry/inversions.jl   (invert_curve and its :roots, :auto, :optim and :grid_optim methods)
# Reference objects: a line, a circular arc and two polar arcs, placed off the origin. For a point p = curve(crv, t0) on the
# curve the inversion must return t0, so the reference value is the parameter used to build the point and no second
# implementation enters the comparison. All four methods are checked against the same parameters; they differ only in
# robustness, not in the result. The point of a curve is a regular point of the squared distance only up to a double root,
# which is why the roots method solves for the stationary point of the distance instead (G5). Points that are not on the
# curve have no zero distance; the roots method must then fail and the automatic method must fall back to the optimizer.
# No eigen solves, no pinned wavefunction values. Helper constants and functions use the prefix `iv_`.

using StaticArrays: SVector
using Logging: Warn

# Absolute tolerance on the recovered parameter (the optimizer is exact to about 1e-16, the Newton solve to a few ulps).
const iv_T_ATOL = 1e-10

# Parameters at which the curves are inverted; none is at the ends, where the wrap of the parameter to [0, 1) is ambiguous.
const iv_TS = (0.07, 0.2, 0.37, 0.5, 0.9)

const iv_METHODS = (:roots, :auto, :optim, :grid_optim)

iv_line() = LineSegment([0.2, -0.1], [1.4, 0.9])
iv_arc() = CircleSegment(1.2, 1.3, 0.4, [0.1, -0.2])
iv_polar() = PolarSegment(φ -> 1.0 + 0.2 * cos(3φ); arc_angle = 1.4, shift_angle = 0.3, center = [0.1, 0.2])
iv_fourier() = FourierCoeffPolarSegment([0.0, 0.0, 0.0, 0.3]; arc_angle = 1.3, shift_angle = 0.4, center = [-0.2, 0.1])

@testset "invert_curve (point on the curve)" begin
    # G5: every method recovers the parameter of a point on a line, a circle and polar curves
    for (name, crv) in (("line", iv_line()), ("circle", iv_arc()), ("polar", iv_polar()), ("fourier", iv_fourier()))
        @testset "$name" begin
            for method in iv_METHODS, t0 in iv_TS
                t = invert_curve(crv, curve(crv, t0); method = method)
                @test t ≈ t0 atol = iv_T_ATOL
            end
        end
    end

    # the default method is :auto and returns the parameter without falling back (no warning)
    crv = iv_polar()
    @test (@test_logs invert_curve(crv, curve(crv, 0.9))) ≈ 0.9 atol = iv_T_ATOL
    @test (@test_logs invert_curve(crv, curve(crv, 0.2))) ≈ 0.2 atol = iv_T_ATOL

    # the starting parameter does not change the result
    for t_init in (0.1, 0.5, 0.8)
        @test invert_curve(crv, curve(crv, 0.37); method = :roots, t_init = t_init) ≈ 0.37 atol = iv_T_ATOL
        @test invert_curve(crv, curve(crv, 0.37); method = :optim, t_init = t_init) ≈ 0.37 atol = iv_T_ATOL
    end

    # the returned parameter lies in [0, 1)
    for method in iv_METHODS
        t = invert_curve(iv_arc(), curve(iv_arc(), 0.6); method = method)
        @test 0.0 <= t < 1.0
    end
end

@testset "invert_curve (point off the curve)" begin
    # the foot point of a point off a line is its orthogonal projection: for pt0 = (0,0), pt1 = (1,0) the parameter is x
    line = LineSegment([0.0, 0.0], [1.0, 0.0])
    p = SVector(0.3, 0.8)
    @test invert_curve(line, p; method = :optim) ≈ 0.3 atol = iv_T_ATOL
    @test invert_curve(line, p; method = :grid_optim) ≈ 0.3 atol = iv_T_ATOL

    # the roots method accepts only points on the curve; the automatic method falls back to the optimizer and warns
    @test_throws Exception invert_curve(line, p; method = :roots)
    t = @test_logs (:warn,) min_level = Warn invert_curve(line, p; method = :auto)
    @test t ≈ 0.3 atol = iv_T_ATOL

    # a point off a circle projects radially onto it: the arc of one radian starts at angle 0, so angle 0.25 is t = 0.25
    arc = CircleSegment(1.0, 1.0, 0.0, [0.0, 0.0])
    q = 2.0 * SVector(cos(0.25), sin(0.25))
    @test invert_curve(arc, q; method = :optim) ≈ 0.25 atol = iv_T_ATOL
    @test invert_curve(arc, q; method = :grid_optim) ≈ 0.25 atol = iv_T_ATOL
end

@testset "invert_curve (grid search)" begin
    # a nearly closed circle (arc angle 1.9 pi) is inverted at parameters far from the default start t = 0.5
    arc = CircleSegment(1.0, 1.9 * pi, 0.0, [0.0, 0.0])
    for t0 in (0.05, 0.95)
        @test invert_curve(arc, curve(arc, t0); method = :grid_optim) ≈ t0 atol = iv_T_ATOL
    end
    # n_grid is honored: a grid of two points still reaches the answer through the optimizer
    @test invert_curve(arc, curve(arc, 0.3); method = :grid_optim, n_grid = 2) ≈ 0.3 atol = iv_T_ATOL
end

@testset "invert_curve (unknown method)" begin
    crv = iv_line()
    @test_throws ErrorException invert_curve(crv, curve(crv, 0.3); method = :bogus)
end
