# Unit tests for polar curves:
#   src/geometry/segments/polarcurves.jl      (FourierCoeffPolarSegment and PolarSegment constructors, polar_radius for a scalar and
#                                              an array of angles, curve for a scalar and an array of parameters, domain_fun)
#   src/geometry/curvederivatives.jl          (tangent and tangent_2 of both polar types, only as the finite-difference contract
#                                              of their curve; the closed-form derivative checks are in curvederivatives_tests.jl)
# Reference object: r(phi) = R + sum_n a_n cos(n phi) + sum_n b_n sin(n phi) with coef = [b1, a1, b2, a2, ...] (sine coefficient
# first), and the closure r(phi) = 1 + 0.3 cos(2 phi). A curve point is center + r(phi) (cos phi, sin phi) with
# phi = shift_angle + t * arc_angle, and the signed domain function is (|x - center| - r(theta)) * orientation. With vanishing
# coefficients the curve is a circle arc of radius R, so lengths are R |arc_angle|. No eigen solves, no pinned wavefunction
# values. Helper constants use the prefix `pc_`.

using StaticArrays: SVector
using LinearAlgebra: norm

# Finite-difference step for the first derivative (central difference, error ~ h^2, roundoff ~ eps/h, balanced at cbrt(eps)).
const pc_FD_H = cbrt(eps(Float64))
# Absolute tolerance of the central difference of the first derivative (worst observed 1.2e-8).
const pc_FD1_ATOL = 2e-7
# Absolute tolerance of the second difference with the same step (roundoff ~ eps/h^2, worst observed 1.7e-5).
const pc_FD2_ATOL = 2e-4
# Absolute tolerance for closed-form comparisons of trigonometric expressions (a few ulps of Float64).
const pc_ATOL = 1e-12

pc_rfunc(phi) = 1.0 + 0.3 * cos(2 * phi)

@testset "FourierCoeffPolarSegment constructor" begin
    # `PolarSegment`->`FourierCoeffPolarSegment` rename, with a `promote_type` fix so mismatched literal types
    # (Int `R`, Float32 `coef`, Float64 `center`) still produce a single consistent `T`.
    seg = FourierCoeffPolarSegment(Float32[0.0, 0.0, 0.0, 0.3]; R=1, center=[0.0, 0.0])
    @test typeof(seg).parameters[1] == Float64

    crv = FourierCoeffPolarSegment([0.1, 0.2, 0.3, 0.4]; R=1.5, arc_angle=pi, shift_angle=0.25, center=[1.0, -1.0],
        orientation=-1, bc=SymmetryWall(1, 2), domain_id=3, segment_id=4)
    @test crv isa FourierCoeffPolarSegment{Float64,SymmetryWall,4}
    @test crv.R == 1.5
    @test crv.coef == SVector(0.1, 0.2, 0.3, 0.4)
    @test crv.arc_angle == Float64(pi)
    @test crv.shift_angle == 0.25
    @test crv.center == SVector(1.0, -1.0)
    @test crv.orientation == -1
    @test crv.bc == SymmetryWall(1, 2)
    @test crv.domain_id == 3
    @test crv.segment_id == 4

    # defaults: a full 2 pi arc without shift about the origin
    default = FourierCoeffPolarSegment([0.0, 0.0])
    @test default.R == 1.0
    @test default.arc_angle == 2.0 * pi
    @test default.shift_angle == 0.0
    @test default.center == SVector(0.0, 0.0)
    @test default.orientation == 1
    @test default.bc isa SpecularReflection

    # vanishing coefficients give a circle arc, so the stored length is R |arc_angle|
    @test default.length ≈ 2.0 * pi atol = pc_ATOL
    @test FourierCoeffPolarSegment([0.0, 0.0]; R=2.0, arc_angle=pi).length ≈ 2.0 * pi atol = pc_ATOL
    @test FourierCoeffPolarSegment(Float32[0.0, 0.0]; R=1.0f0, arc_angle=1.0f0, shift_angle=0.0f0, center=Float32[0, 0]) isa
          FourierCoeffPolarSegment{Float32}
end

@testset "PolarSegment constructor" begin
    seg = PolarSegment(pc_rfunc)
    @test seg isa PolarSegment{Float64,SpecularReflection}
    @test seg.r_func === pc_rfunc
    @test seg.R == 1.0
    @test seg.arc_angle == 2.0 * pi
    @test seg.shift_angle == 0.0
    @test seg.center == SVector(0.0, 0.0)
    @test seg.orientation == 1
    @test seg.bc isa SpecularReflection

    tagged = PolarSegment(pc_rfunc; arc_angle=pi, shift_angle=0.5, center=[0.5, 0.5], orientation=-1, bc=SymmetryWall(2, 3),
        domain_id=7, segment_id=8)
    @test tagged.arc_angle == Float64(pi)
    @test tagged.shift_angle == 0.5
    @test tagged.center == SVector(0.5, 0.5)
    @test tagged.orientation == -1
    @test tagged.bc == SymmetryWall(2, 3)
    @test tagged.domain_id == 7
    @test tagged.segment_id == 8

    # an Int radius still promotes to a single Float64 type parameter
    @test PolarSegment(pc_rfunc; R=1) isa PolarSegment{Float64}

    # a constant radial function is a circle: the length is 2 pi R, the half arc has length pi R
    @test PolarSegment(phi -> 2.0).length ≈ 4.0 * pi atol = pc_ATOL
    @test PolarSegment(phi -> 2.0; arc_angle=pi).length ≈ 2.0 * pi atol = pc_ATOL
end

@testset "polar_radius" begin
    # coef = [b1, a1, b2, a2]: the sine coefficient comes first
    crv = FourierCoeffPolarSegment([0.1, 0.2, 0.3, 0.4]; R=1.5)
    expected(phi) = 1.5 + 0.2 * cos(phi) + 0.4 * cos(2 * phi) + 0.1 * sin(phi) + 0.3 * sin(2 * phi)
    @test polar_radius(crv, 0.0) ≈ 1.5 + 0.2 + 0.4 atol = pc_ATOL
    @test polar_radius(crv, pi / 2) ≈ 1.5 + 0.1 - 0.4 atol = pc_ATOL
    @test polar_radius(crv, pi) ≈ 1.5 - 0.2 + 0.4 atol = pc_ATOL
    @test polar_radius(crv, 0.7) ≈ expected(0.7) atol = pc_ATOL
    # r is 2 pi periodic
    @test polar_radius(crv, 0.7 + 2 * pi) ≈ polar_radius(crv, 0.7) atol = pc_ATOL

    # a sine-only and a cosine-only second harmonic isolate the coefficient order
    @test polar_radius(FourierCoeffPolarSegment([0.0, 0.0, 0.5, 0.0]; R=1.0), pi / 4) ≈ 1.5 atol = pc_ATOL
    @test polar_radius(FourierCoeffPolarSegment([0.0, 0.0, 0.0, 0.5]; R=1.0), pi / 4) ≈ 1.0 atol = pc_ATOL

    # array method evaluates the scalar method elementwise
    phis = [0.0, 0.3, 1.2, 2.5, -1.0]
    @test polar_radius(crv, phis) == [polar_radius(crv, phi) for phi in phis]
    @test polar_radius(crv, phis) ≈ expected.(phis) atol = pc_ATOL
    @test polar_radius(crv, Float64[]) == Float64[]

    # PolarSegment forwards to its radial function
    pseg = PolarSegment(pc_rfunc)
    @test polar_radius(pseg, 0.4) == pc_rfunc(0.4)
    @test polar_radius(pseg, phis) == pc_rfunc.(phis)
end

@testset "curve(::AbsPolarCurve)" begin
    crv = FourierCoeffPolarSegment([0.0, 0.0, 0.0, 0.3]; R=1.0, center=[0.5, -0.5])
    # t = 0 is phi = 0: r = 1 + 0.3, t = 1/4 is phi = pi/2: r = 1 - 0.3
    @test curve(crv, 0.0) ≈ SVector(0.5 + 1.3, -0.5) atol = pc_ATOL
    @test curve(crv, 0.25) ≈ SVector(0.5, -0.5 + 0.7) atol = pc_ATOL
    # closed curve: the end point equals the start point
    @test curve(crv, 1.0) ≈ curve(crv, 0.0) atol = 1e-10

    # shift and arc angle select the angular range, phi = shift + t * arc_angle
    part = FourierCoeffPolarSegment([0.0, 0.0]; R=2.0, arc_angle=pi / 2, shift_angle=pi / 2)
    @test curve(part, 0.0) ≈ SVector(0.0, 2.0) atol = pc_ATOL
    @test curve(part, 1.0) ≈ SVector(-2.0, 0.0) atol = pc_ATOL

    # array method evaluates the scalar method elementwise, for both polar types
    ts = [0.0, 0.1, 0.35, 0.8, 1.0]
    pts = curve(crv, ts)
    @test pts isa Vector
    @test all(isapprox.(pts, [curve(crv, t) for t in ts]; atol=pc_ATOL))
    pseg = PolarSegment(pc_rfunc; center=[0.5, -0.5])
    @test all(isapprox.(curve(pseg, ts), [curve(pseg, t) for t in ts]; atol=pc_ATOL))
    # PolarSegment and FourierCoeffPolarSegment describe the same curve r = 1 + 0.3 cos(2 phi)
    @test all(isapprox.(curve(pseg, ts), pts; atol=pc_ATOL))
end

@testset "domain_fun (::AbsPolarCurve)" begin
    crv = FourierCoeffPolarSegment([0.0, 0.0, 0.0, 0.3]; R=1.0, center=[0.5, -0.5])
    # points at angle theta and radius rho about the center have domain function rho - r(theta)
    pt(rho, theta) = SVector(0.5 + rho * cos(theta), -0.5 + rho * sin(theta))
    @test domain_fun(crv, pt(0.5, 0.0)) ≈ 0.5 - 1.3 atol = pc_ATOL
    @test domain_fun(crv, pt(2.0, pi / 2)) ≈ 2.0 - 0.7 atol = pc_ATOL
    @test domain_fun(crv, pt(1.3, 0.0)) ≈ 0.0 atol = pc_ATOL
    @test is_inside(crv, pt(0.5, 0.0))
    @test !is_inside(crv, pt(2.0, pi / 2))

    # the curve itself lies on the zero level set
    for t in (0.0, 0.2, 0.55, 0.9)
        @test domain_fun(crv, curve(crv, t)) ≈ 0.0 atol = pc_ATOL
    end

    # array method and orientation sign
    testpts = [pt(0.5, 0.0), pt(2.0, pi / 2)]
    @test domain_fun(crv, testpts) ≈ [0.5 - 1.3, 2.0 - 0.7] atol = pc_ATOL
    hole = FourierCoeffPolarSegment([0.0, 0.0, 0.0, 0.3]; R=1.0, center=[0.5, -0.5], orientation=-1)
    @test domain_fun(hole, testpts) ≈ -domain_fun(crv, testpts) atol = pc_ATOL
    @test is_inside(hole, testpts) == [false, true]

    # PolarSegment gives the same domain function as the equivalent Fourier curve
    pseg = PolarSegment(pc_rfunc; center=[0.5, -0.5])
    @test domain_fun(pseg, testpts) ≈ domain_fun(crv, testpts) atol = pc_ATOL
end

@testset "polar curve derivatives against finite differences" begin
    # New function-based `PolarSegment{T,BC,F}`: derivatives via ForwardDiff, checked against central finite differences of
    # `curve(t)`. The Fourier curve uses the analytic formulas of curvederivatives.jl and is checked the same way.
    h = pc_FD_H
    pseg = PolarSegment(pc_rfunc)
    fseg = FourierCoeffPolarSegment([0.1, 0.2, 0.3, 0.4]; R=1.5)
    for crv in (pseg, fseg), t0 in (0.1, 0.4, 0.7)
        xp = curve(crv, t0 + h)
        xm = curve(crv, t0 - h)
        x0 = curve(crv, t0)
        fd1 = (xp .- xm) ./ (2h)
        fd2 = (xp .- 2 .* x0 .+ xm) ./ (h^2)
        @test all(isapprox.(fd1, tangent(crv, t0); atol=pc_FD1_ATOL))
        @test all(isapprox.(fd2, tangent_2(crv, t0); atol=pc_FD2_ATOL))
    end
    # array method of tangent broadcasts the scalar method
    ts = [0.1, 0.4, 0.7]
    @test tangent(pseg, ts) == [tangent(pseg, t) for t in ts]
    @test tangent(fseg, ts) == [tangent(fseg, t) for t in ts]
end
