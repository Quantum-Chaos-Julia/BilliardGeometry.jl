# Unit tests for parameter derivatives of curves:
#   src/geometry/curvederivatives.jl   (tangent and tangent_2 of LineSegment, CircleSegment, FourierCoeffPolarSegment and
#                                       PolarSegment, the generic ForwardDiff fallback of tangent, the array methods, tangent_vec,
#                                       normal_vec and curvature)
# Reference object: the circle r(t) = c + R (cos(a t + s), sin(a t + s)) with tangent R a (-sin, cos), tangent_2 -R a^2 (cos, sin),
# unit tangent (-sin, cos), unit outward normal (cos, sin) and curvature 1/R (sign of a), a line with constant tangent pt1 - pt0
# and zero second derivative, and the parabola (t, t^2). The first and second derivatives of every concrete curve are also
# compared with central finite differences of `curve`, and the analytic Fourier polar derivatives with the ForwardDiff based
# derivatives of the PolarSegment of the same radial function. The scalar methods tangent_vec(crv, t::Real) and
# normal_vec(crv, t::Real) are compared elementwise with the array methods and with the closed form of the line, circle and a
# polar circle: unit length, normal = (ty, -tx), orthogonal to the tangent. No eigen solves, no pinned wavefunction values.
# Helper functions and constants use the prefix `cd_`.

using StaticArrays: SVector
using LinearAlgebra: norm, dot

# Finite-difference step (central difference, error ~ h^2, roundoff ~ eps/h, balanced at cbrt(eps)).
const cd_FD_H = cbrt(eps(Float64))
# Absolute tolerance of the central difference of the first derivative (worst observed 1.2e-9).
const cd_FD1_ATOL = 2e-8
# Absolute tolerance of the second difference with the same step (roundoff ~ eps/h^2, worst observed 3.5e-5).
const cd_FD2_ATOL = 4e-4
# Absolute tolerance for closed-form and ForwardDiff comparisons of trigonometric expressions (a few ulps of Float64 times
# the size R a^2 ~ 10 of the derivatives).
const cd_ATOL = 1e-10

# Deviation of the analytic first and second derivative of `crv` at `t0` from finite differences of `curve`.
function cd_fd_deviation(crv, t0; h=cd_FD_H)
    xp = curve(crv, t0 + h)
    xm = curve(crv, t0 - h)
    x0 = curve(crv, t0)
    fd1 = (xp .- xm) ./ (2h)
    fd2 = (xp .- 2 .* x0 .+ xm) ./ (h^2)
    return maximum(abs.(fd1 .- tangent(crv, t0))), maximum(abs.(fd2 .- tangent_2(crv, t0)))
end

# A curve without dedicated derivative methods, to exercise the generic ForwardDiff fallback of `tangent`.
struct cd_Parabola <: AbsCurve{SpecularReflection} end
BilliardGeometry.curve(::cd_Parabola, t) = SVector(t, t^2)

cd_rfunc(phi) = 1.0 + 0.3 * cos(2 * phi)

@testset "tangent and tangent_2 (::CircleSegment)" begin
    # compare against central finite differences of curve(t)
    circ = CircleSegment(2.0, 0.5 * pi, 0.3, [0.1, -0.2])
    for t0 in (0.13, 0.5, 0.87)
        dev1, dev2 = cd_fd_deviation(circ, t0)
        @test dev1 < cd_FD1_ATOL
        @test dev2 < cd_FD2_ATOL
    end

    # closed form: R a (-sin(phi), cos(phi)) and -R a^2 (cos(phi), sin(phi)) with phi = a t + s
    R, a, s = 2.0, 0.5 * pi, 0.3
    for t0 in (0.0, 0.4, 1.0)
        phi = a * t0 + s
        @test tangent(circ, t0) ≈ R * a * SVector(-sin(phi), cos(phi)) atol = cd_ATOL
        @test tangent_2(circ, t0) ≈ -R * a^2 * SVector(cos(phi), sin(phi)) atol = cd_ATOL
    end

    # a negative arc reverses the sense of traversal: the tangent flips sign, tangent_2 does not
    cw = CircleSegment(2.0, -0.5 * pi, 0.3, [0.1, -0.2])
    @test tangent(cw, 0.4) ≈ -R * a * SVector(-sin(-a * 0.4 + s), cos(-a * 0.4 + s)) atol = cd_ATOL
    @test tangent_2(cw, 0.4) ≈ -R * a^2 * SVector(cos(-a * 0.4 + s), sin(-a * 0.4 + s)) atol = cd_ATOL
    dev1, dev2 = cd_fd_deviation(cw, 0.4)
    @test dev1 < cd_FD1_ATOL
    @test dev2 < cd_FD2_ATOL
end

@testset "tangent and tangent_2 (::LineSegment)" begin
    # tangent is the constant chord vector, tangent_2 is zero
    line = LineSegment([0.0, 0.0], [2.0, -1.0])
    @test tangent(line, 0.0) == line.pt1 .- line.pt0
    @test tangent(line, 0.7) == line.pt1 .- line.pt0
    @test tangent_2(line, 0.4) == zero(line.pt0)
    for t0 in (0.2, 0.8)
        dev1, dev2 = cd_fd_deviation(line, t0)
        @test dev1 < cd_FD1_ATOL
        @test dev2 < cd_FD2_ATOL
    end
end

@testset "tangent (generic AbsCurve fallback)" begin
    # ForwardDiff derivative of `curve`: d/dt (t, t^2) = (1, 2t)
    parabola = cd_Parabola()
    @test tangent(parabola, 0.0) == SVector(1.0, 0.0)
    @test tangent(parabola, 0.3) ≈ SVector(1.0, 0.6) atol = cd_ATOL
    @test tangent(parabola, [0.0, 0.5]) ≈ [SVector(1.0, 0.0), SVector(1.0, 1.0)] atol = cd_ATOL
end

@testset "tangent and tangent_2 (::PolarSegment)" begin
    # derivatives via nested ForwardDiff, checked against finite differences and against the closed form of a circle
    pseg = PolarSegment(cd_rfunc; arc_angle=pi, shift_angle=0.2, center=[0.3, -0.4])
    for t0 in (0.1, 0.5, 0.9)
        dev1, dev2 = cd_fd_deviation(pseg, t0)
        @test dev1 < cd_FD1_ATOL
        @test dev2 < cd_FD2_ATOL
    end

    # r(phi) = 2 is a circle of radius 2: tangent_2 = -R a^2 (cos(phi), sin(phi))
    circle = PolarSegment(phi -> 2.0)
    a = 2.0 * pi
    for t0 in (0.0, 0.3, 0.7)
        phi = a * t0
        @test tangent(circle, t0) ≈ 2.0 * a * SVector(-sin(phi), cos(phi)) atol = cd_ATOL
        @test tangent_2(circle, t0) ≈ -2.0 * a^2 * SVector(cos(phi), sin(phi)) atol = cd_ATOL
    end
end

@testset "tangent and tangent_2 (::FourierCoeffPolarSegment)" begin
    # vanishing coefficients: circle of radius R, tangent_2 = -R a^2 (cos(phi), sin(phi)) with a = 2 pi
    R = 1.5
    circle = FourierCoeffPolarSegment([0.0, 0.0]; R=R)
    a = 2.0 * pi
    for t0 in (0.0, 0.3, 0.7)
        phi = a * t0
        @test tangent(circle, t0) ≈ R * a * SVector(-sin(phi), cos(phi)) atol = cd_ATOL
        @test tangent_2(circle, t0) ≈ -R * a^2 * SVector(cos(phi), sin(phi)) atol = cd_ATOL
    end

    # analytic first and second derivative path against finite differences, with shift, partial arc and center
    fseg = FourierCoeffPolarSegment([0.1, 0.2, 0.3, 0.4]; R=1.5, arc_angle=pi, shift_angle=0.2, center=[0.3, -0.4])
    for t0 in (0.1, 0.5, 0.9)
        dev1, dev2 = cd_fd_deviation(fseg, t0)
        @test dev1 < cd_FD1_ATOL
        @test dev2 < cd_FD2_ATOL
    end

    # the analytic formulas equal the ForwardDiff derivatives of a PolarSegment with the same radial function
    # r(phi) = R + 0.2 cos(phi) + 0.4 cos(2 phi) + 0.1 sin(phi) + 0.3 sin(2 phi)
    rfun(phi) = 1.5 + 0.2 * cos(phi) + 0.4 * cos(2 * phi) + 0.1 * sin(phi) + 0.3 * sin(2 * phi)
    pseg = PolarSegment(rfun; arc_angle=pi, shift_angle=0.2, center=[0.3, -0.4])
    for t0 in (0.0, 0.25, 0.6, 1.0)
        @test tangent(fseg, t0) ≈ tangent(pseg, t0) atol = cd_ATOL
        @test tangent_2(fseg, t0) ≈ tangent_2(pseg, t0) atol = cd_ATOL
    end
end

@testset "tangent and tangent_2 array methods" begin
    # array-argument methods broadcast the scalar method
    circ = CircleSegment(2.0, 0.5 * pi, 0.3, [0.1, -0.2])
    ts = [0.1, 0.4, 0.9]
    @test tangent(circ, ts) == [tangent(circ, t) for t in ts]
    @test tangent_2(circ, ts) == [tangent_2(circ, t) for t in ts]
    @test tangent(circ, Float64[]) == SVector{2,Float64}[]

    line = LineSegment([0.0, 0.0], [2.0, -1.0])
    @test tangent(line, ts) == [tangent(line, t) for t in ts]
    @test tangent_2(line, ts) == [tangent_2(line, t) for t in ts]

    pseg = PolarSegment(cd_rfunc)
    fseg = FourierCoeffPolarSegment([0.0, 0.0, 0.0, 0.3])
    @test tangent(pseg, ts) == [tangent(pseg, t) for t in ts]
    @test tangent_2(pseg, ts) == [tangent_2(pseg, t) for t in ts]
    @test tangent(fseg, ts) == [tangent(fseg, t) for t in ts]
    @test tangent_2(fseg, ts) == [tangent_2(fseg, t) for t in ts]
end

@testset "curvature, tangent_vec and normal_vec" begin
    # `tangent_vec`/`normal_vec`/`curvature`. A circle of radius R has constant curvature 1/R, and its unit outward normal at
    # parameter t is the radial direction (cos(phi), sin(phi)).
    R = 2.0
    circ0 = CircleSegment(R, 2 * pi, 0.0, [0.0, 0.0])
    ts0 = collect(range(0, 1, length=5))
    @test all(isapprox.(curvature(circ0, ts0), 1 / R; atol=1e-12))
    @test isapprox(curvature(circ0, ts0[2]), 1 / R; atol=1e-12)
    tv = tangent_vec(circ0, ts0)
    @test all(isapprox.(norm.(tv), 1.0; atol=1e-12))
    nv = normal_vec(circ0, ts0)
    phis = 2 * pi .* ts0
    @test all(isapprox.(nv, [SVector(cos(phi), sin(phi)) for phi in phis]; atol=1e-10))

    # the unit tangent is (-sin, cos) and is orthogonal to the unit normal
    @test all(isapprox.(tv, [SVector(-sin(phi), cos(phi)) for phi in phis]; atol=cd_ATOL))
    @test all(abs.(dot.(tv, nv)) .< cd_ATOL)

    # a clockwise circle has negative curvature, a line has zero curvature
    cw = CircleSegment(R, -0.5 * pi, 0.0, [0.0, 0.0])
    @test all(isapprox.(curvature(cw, [0.1, 0.6]), -1 / R; atol=cd_ATOL))
    line = LineSegment([0.0, 0.0], [2.0, -1.0])
    @test all(iszero, curvature(line, [0.1, 0.6]))
    @test curvature(line, 0.3) == 0.0

    # the polar circle of radius 2 has curvature 1/2 for both polar types
    @test curvature(PolarSegment(phi -> 2.0), 0.3) ≈ 0.5 atol = cd_ATOL
    @test all(isapprox.(curvature(FourierCoeffPolarSegment([0.0, 0.0]; R=2.0), [0.2, 0.7]), 0.5; atol=cd_ATOL))

    # a line has a constant unit tangent along the chord and a normal rotated clockwise by 90 degrees
    chord = SVector(3.0, 4.0)
    l2 = LineSegment([0.0, 0.0], chord)
    @test tangent_vec(l2, [0.2, 0.8]) == [chord / 5.0, chord / 5.0]
    @test normal_vec(l2, [0.2, 0.8]) == [SVector(0.8, -0.6), SVector(0.8, -0.6)]
end

# Checks of the scalar unit tangent and unit normal of `crv` at the parameters `ts`: agreement with the array methods, unit
# length, normal = (ty, -tx) and orthogonality.
function cd_check_scalar_unit_vectors(crv, ts)
    for t0 in ts
        tv = tangent_vec(crv, t0)
        nv = normal_vec(crv, t0)
        @test tv isa SVector{2}
        @test nv isa SVector{2}
        @test norm(tv) ≈ 1.0 atol = cd_ATOL
        @test norm(nv) ≈ 1.0 atol = cd_ATOL
        @test nv ≈ SVector(tv[2], -tv[1]) atol = cd_ATOL
        @test abs(dot(tv, nv)) < cd_ATOL
        @test tv ≈ tangent_vec(crv, [t0])[1] atol = cd_ATOL
        @test nv ≈ normal_vec(crv, [t0])[1] atol = cd_ATOL
    end
end

@testset "tangent_vec and normal_vec (scalar, ::LineSegment)" begin
    # G11: scalar tangent_vec/normal_vec do not exist, only the array methods do
    chord = SVector(3.0, 4.0)
    line = LineSegment([0.0, 0.0], chord)
    for t0 in (0.0, 0.3, 1.0)
        @test tangent_vec(line, t0) ≈ SVector(0.6, 0.8) atol = cd_ATOL
        @test normal_vec(line, t0) ≈ SVector(0.8, -0.6) atol = cd_ATOL
    end
    cd_check_scalar_unit_vectors(line, (0.0, 0.3, 1.0))
    # elementwise agreement with the array method
    ts = [0.1, 0.5, 0.9]
    @test [tangent_vec(line, t) for t in ts] ≈ tangent_vec(line, ts)
    @test [normal_vec(line, t) for t in ts] ≈ normal_vec(line, ts)

    # type stability and eltype
    @test @inferred(tangent_vec(line, 0.3)) isa SVector{2,Float64}
    @test @inferred(normal_vec(line, 0.3)) isa SVector{2,Float64}
    line32 = LineSegment(Float32[0, 0], Float32[3, 4])
    @test tangent_vec(line32, 0.3) isa SVector{2,Float32}
    @test normal_vec(line32, 0.3) isa SVector{2,Float32}
end

@testset "tangent_vec and normal_vec (scalar, ::CircleSegment)" begin
    # G11: counter-clockwise circle, unit tangent (-sin(phi), cos(phi)) and outward unit normal (cos(phi), sin(phi))
    a, s = 2 * pi, 0.3
    circ = CircleSegment(2.0, a, s, [0.1, -0.2])
    for t0 in (0.0, 0.2, 0.55, 0.9)
        phi = a * t0 + s
        @test tangent_vec(circ, t0) ≈ SVector(-sin(phi), cos(phi)) atol = cd_ATOL
        @test normal_vec(circ, t0) ≈ SVector(cos(phi), sin(phi)) atol = cd_ATOL
    end
    cd_check_scalar_unit_vectors(circ, (0.0, 0.2, 0.55, 0.9))
    ts = [0.1, 0.5, 0.9]
    @test [tangent_vec(circ, t) for t in ts] ≈ tangent_vec(circ, ts)
    @test [normal_vec(circ, t) for t in ts] ≈ normal_vec(circ, ts)

    # a clockwise arc flips the tangent, and with it the normal
    cw = CircleSegment(2.0, -0.5 * pi, 0.3, [0.1, -0.2])
    phi = -0.5 * pi * 0.4 + 0.3
    @test tangent_vec(cw, 0.4) ≈ SVector(sin(phi), -cos(phi)) atol = cd_ATOL
    @test normal_vec(cw, 0.4) ≈ SVector(-cos(phi), -sin(phi)) atol = cd_ATOL

    @test @inferred(tangent_vec(circ, 0.3)) isa SVector{2,Float64}
    @test @inferred(normal_vec(circ, 0.3)) isa SVector{2,Float64}
end

@testset "tangent_vec and normal_vec (scalar, polar segments)" begin
    # G11: polar circle of radius 2, same closed form as the CircleSegment
    a = 2.0 * pi
    pcirc = PolarSegment(phi -> 2.0)
    for t0 in (0.1, 0.45, 0.8)
        phi = a * t0
        @test tangent_vec(pcirc, t0) ≈ SVector(-sin(phi), cos(phi)) atol = cd_ATOL
        @test normal_vec(pcirc, t0) ≈ SVector(cos(phi), sin(phi)) atol = cd_ATOL
    end

    # a non-circular radial function: unit length, orthogonality and agreement with the array method
    pseg = PolarSegment(cd_rfunc; arc_angle=pi, shift_angle=0.2, center=[0.3, -0.4])
    cd_check_scalar_unit_vectors(pseg, (0.1, 0.5, 0.9))
    fseg = FourierCoeffPolarSegment([0.1, 0.2, 0.3, 0.4]; R=1.5, arc_angle=pi, shift_angle=0.2, center=[0.3, -0.4])
    cd_check_scalar_unit_vectors(fseg, (0.1, 0.5, 0.9))
    ts = [0.1, 0.5, 0.9]
    @test [tangent_vec(pseg, t) for t in ts] ≈ tangent_vec(pseg, ts)
    @test [normal_vec(fseg, t) for t in ts] ≈ normal_vec(fseg, ts)
end
