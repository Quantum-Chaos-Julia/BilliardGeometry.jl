# Unit tests for circular arcs:
#   src/geometry/segments/circlesegment.jl    (CircleSegment constructors, curve, domain_fun, arc_length(circle, pt))
#   src/geometry/arclength.jl                 (generic arc_length(crv, t) and arc_length(crv, ts) applied to a CircleSegment)
#   src/geometry/geometry.jl                  (is_inside and domain_gradient_vector of a single curve)
# Reference object: arcs of the circle of radius R about a center c, r(t) = c + R (cos(a t + s), sin(a t + s)) with arc
# angle a and shift angle s. Every checked quantity has a closed form: the arc length is R |a| t, the signed domain function is
# (|x - c| - R) * orientation and its gradient is the radial unit vector. The G4 regression for arc_length(circle, pt) on
# negative arcs is not part of this file. No eigen solves, no pinned wavefunction values. Helper constants use the prefix `cs_`.

using StaticArrays: SVector
using LinearAlgebra: norm

# Absolute tolerance for closed-form comparisons of trigonometric expressions (a few ulps of Float64).
const cs_ATOL = 1e-12
# Absolute tolerance for quadrature-based arc lengths (quadgk with its default rtol = sqrt(eps) of the element type).
const cs_ARCLEN_ATOL = 1e-7

@testset "CircleSegment constructors" begin
    crv = CircleSegment(2.0, 0.5 * pi, 0.3, [0.1, -0.2])
    @test crv isa CircleSegment{Float64,SpecularReflection}
    @test crv.radius == 2.0
    @test crv.arc_angle == 0.5 * pi
    @test crv.shift_angle == 0.3
    @test crv.center == SVector(0.1, -0.2)
    @test crv.length == 2.0 * 0.5 * pi
    @test crv.orientation == 1
    @test crv.domain_id == 1
    @test crv.segment_id == 1
    @test crv.bc isa SpecularReflection

    # the (x0, y0) method builds the same arc as the center-vector method
    @test CircleSegment(2.0, 0.5 * pi, 0.3, 0.1, -0.2) == crv
    tagged = CircleSegment(1.0, pi, 0.0, 0.5, 0.5; bc=SymmetryWall(1, 2), orientation=-1, domain_id=3, segment_id=4)
    @test tagged.bc == SymmetryWall(1, 2)
    @test tagged.orientation == -1
    @test tagged.domain_id == 3
    @test tagged.segment_id == 4

    # mixed argument types promote to a single T
    @test CircleSegment(1, 0.5 * pi, 0.0, [0.0, 0.0]) isa CircleSegment{Float64}
    @test CircleSegment(1.0f0, 0.5f0, 0.0f0, Float32[0, 0]) isa CircleSegment{Float32}
    @test CircleSegment(1.0f0, 0.5f0, 0.0f0, 0.0f0, 0.0f0) isa CircleSegment{Float32}
    @test CircleSegment(1.0f0, 0.5f0, 0.0, 0.0f0, 0.0f0) isa CircleSegment{Float64}

    # the length is R |arc_angle|, also for a negative arc
    @test CircleSegment(2.0, -0.5 * pi, 0.0, [0.0, 0.0]).length == 2.0 * 0.5 * pi
    @test CircleSegment(1.5, 2 * pi, 0.0, [0.0, 0.0]).length ≈ 3.0 * pi
end

@testset "curve(::CircleSegment)" begin
    R = 1.0
    arc_angle = 0.5 * pi
    shift_angle = 0.0
    center = [0.0, 0.0]
    crv = CircleSegment(R, arc_angle, shift_angle, center)
    @test all(isapprox.(curve(crv, 0.0), SVector{2,Float64}([1.0, 0.0]); atol=1e-8))
    @test all(isapprox.(curve(crv, 1.0), SVector{2,Float64}([0.0, 1.0]); atol=1e-8))
    @test curve(crv, 0.5) ≈ SVector(sqrt(2) / 2, sqrt(2) / 2) atol = cs_ATOL

    # points lie at distance R from the center for any t
    shifted = CircleSegment(2.0, 1.3, 0.4, [0.1, -0.2])
    for t in (0.0, 0.2, 0.5, 0.9, 1.0)
        @test norm(curve(shifted, t) - shifted.center) ≈ 2.0 atol = cs_ATOL
    end

    # array method evaluates the scalar method elementwise
    ts = [0.0, 0.25, 0.5, 1.0]
    pts = curve(shifted, ts)
    @test pts isa Vector{<:SVector{2,Float64}}
    @test pts == [curve(shifted, t) for t in ts]
end

@testset "CircleSegment shift_angle" begin
    # shift_angle rotates the start point: curve(0) = c + R (cos s, sin s)
    s = 0.3
    crv = CircleSegment(1.0, 0.5 * pi, s, [0.0, 0.0])
    @test curve(crv, 0.0) ≈ SVector(cos(s), sin(s)) atol = cs_ATOL
    @test curve(crv, 1.0) ≈ SVector(cos(s + 0.5 * pi), sin(s + 0.5 * pi)) atol = cs_ATOL

    # a shift of a quarter turn maps the start point (R, 0) to (0, R)
    quarter = CircleSegment(1.0, 0.5 * pi, 0.5 * pi, [0.0, 0.0])
    @test curve(quarter, 0.0) ≈ SVector(0.0, 1.0) atol = cs_ATOL
    @test curve(quarter, 1.0) ≈ SVector(-1.0, 0.0) atol = cs_ATOL

    # the shift does not change the length nor the domain function
    unshifted = CircleSegment(1.0, 0.5 * pi, 0.0, [0.0, 0.0])
    @test crv.length == unshifted.length
    @test domain_fun(crv, SVector(0.1, 0.2)) == domain_fun(unshifted, SVector(0.1, 0.2))

    # arc_length(circle, pt) for a positive arc that does not cross the branch cut of atan: R * a * t
    for t in (0.0, 0.25, 0.5, 1.0)
        @test arc_length(crv, curve(crv, t)) ≈ crv.length * t atol = cs_ARCLEN_ATOL
    end
end

@testset "CircleSegment negative arc_angle" begin
    R = 2.0
    crv = CircleSegment(R, -0.5 * pi, 0.0, [0.0, 0.0])
    # clockwise traversal from (R, 0) to (0, -R)
    @test curve(crv, 0.0) ≈ SVector(R, 0.0) atol = cs_ATOL
    @test curve(crv, 1.0) ≈ SVector(0.0, -R) atol = cs_ATOL
    @test curve(crv, 0.5) ≈ SVector(R * sqrt(2) / 2, -R * sqrt(2) / 2) atol = cs_ATOL
    @test crv.length == R * 0.5 * pi

    # the tangent points clockwise: at t = 0 it is (0, -R |a|)
    @test tangent(crv, 0.0) ≈ SVector(0.0, -R * 0.5 * pi) atol = cs_ATOL
    @test norm(tangent(crv, 0.37)) ≈ crv.length atol = cs_ATOL
    @test curvature(crv, 0.3) ≈ -1 / R atol = cs_ATOL

    # the arc length integrates |tangent| and is therefore non-negative: R |a| t
    @test arc_length(crv, 1.0) ≈ crv.length atol = cs_ARCLEN_ATOL
    @test arc_length(crv, 0.5) ≈ 0.5 * crv.length atol = cs_ARCLEN_ATOL

    # the domain function only depends on the distance from the center, not on the sense of traversal
    pos = CircleSegment(R, 0.5 * pi, 0.0, [0.0, 0.0])
    inside, outside = SVector(0.5, 0.5), SVector(3.0, 0.0)
    @test domain_fun(crv, inside) == domain_fun(pos, inside)
    @test domain_fun(crv, outside) == domain_fun(pos, outside)
    @test is_inside(crv, [inside, outside]) == [true, false]
end

@testset "arc_length(::CircleSegment, t) and arc_length(::CircleSegment, pt)" begin
    R = 1.5
    crv = CircleSegment(R, 0.5 * pi, 0.0, [0.0, 0.0])
    L = crv.length
    @test L == R * 0.5 * pi

    # point method with the end point of the arc (existing positive-arc behaviour)
    unit_arc = CircleSegment(1.0, 0.5 * pi, 0.0, [0.0, 0.0])
    @test arc_length(unit_arc, curve(unit_arc, 1.0)) ≈ 0.5 * pi
    @test arc_length(crv, curve(crv, 1.0)) ≈ L atol = cs_ARCLEN_ATOL
    @test arc_length(crv, curve(crv, 0.0)) == 0.0

    # parameter method integrates |tangent| = R |a|, so arc_length(crv, t) = t L
    @test arc_length(crv, 1.0) ≈ L atol = cs_ARCLEN_ATOL
    @test arc_length(crv, 0.0) == 0.0
    ts = [0.0, 0.2, 0.5, 1.0]
    @test arc_length(crv, ts) ≈ ts .* L atol = cs_ARCLEN_ATOL

    # a full circle has length 2 pi R
    full = CircleSegment(R, 2 * pi, 0.0, [0.3, -0.1])
    @test arc_length(full, 1.0) ≈ 2 * pi * R atol = cs_ARCLEN_ATOL
    @test full.length ≈ 2 * pi * R atol = cs_ATOL
end

@testset "domain_fun and is_inside (::CircleSegment)" begin
    crv = CircleSegment(1.0, 0.5 * pi, 0.0, [0.0, 0.0])
    @test domain_gradient_vector(crv, curve(crv, 0.5)) ≈ SVector{2,Float64}([sqrt(2) / 2, sqrt(2) / 2])

    testpt1 = SVector{2,Float64}([0.1, 0.2])
    testpt2 = SVector{2,Float64}([1.1, 1.0])
    @test domain_fun(crv, testpt1) < 0.0
    @test domain_fun(crv, testpt2) > 0.0
    @test is_inside(crv, [testpt1, testpt2]) == [true, false]

    # signed distance to the circle: |x - c| - R, exactly zero on it
    @test domain_fun(crv, SVector(3.0, 4.0)) ≈ 4.0 atol = cs_ATOL
    @test domain_fun(crv, SVector(0.6, 0.0)) ≈ -0.4 atol = cs_ATOL
    @test domain_fun(crv, SVector(1.0, 0.0)) == 0.0
    @test is_inside(crv, SVector(1.0, 0.0))
    @test is_inside(crv, testpt1)
    @test !is_inside(crv, testpt2)

    # array methods
    @test domain_fun(crv, [testpt1, testpt2]) ≈ [hypot(0.1, 0.2) - 1.0, hypot(1.1, 1.0) - 1.0] atol = cs_ATOL
    gs = domain_gradient_vector(crv, [SVector(3.0, 4.0), SVector(0.0, -2.0)])
    @test gs[1] ≈ SVector(0.6, 0.8) atol = cs_ATOL
    @test gs[2] ≈ SVector(0.0, -1.0) atol = cs_ATOL

    # a non-zero center shifts the circle
    off = CircleSegment(1.0, 0.5 * pi, 0.0, [2.0, 3.0])
    @test domain_fun(off, SVector(2.0, 3.0)) == -1.0
    @test domain_fun(off, SVector(3.0, 3.0)) == 0.0

    # orientation -1 flips the sign and therefore the inside region (outside of the circle, the hole convention)
    hole = CircleSegment(1.0, 0.5 * pi, 0.0, [0.0, 0.0]; orientation=-1)
    @test domain_fun(hole, testpt1) == -domain_fun(crv, testpt1)
    @test is_inside(hole, [testpt1, testpt2]) == [false, true]
end

@testset "CircleSegment Float32 element type" begin
    crv = CircleSegment(1.0f0, 0.5f0 * Float32(pi), 0.0f0, Float32[0, 0])
    t = 0.5f0
    pt = SVector(0.5f0, 0.5f0)
    @test @inferred(curve(crv, t)) isa SVector{2,Float32}
    @test @inferred(domain_fun(crv, pt)) isa Float32
    @test @inferred(tangent(crv, t)) isa SVector{2,Float32}
    @test @inferred(arc_length(crv, t)) isa Float32
    # Float32 quadgk uses rtol = sqrt(eps(Float32)) ~ 3.4e-4
    @test arc_length(crv, 1.0f0) ≈ crv.length rtol = 1.0f-3
end
