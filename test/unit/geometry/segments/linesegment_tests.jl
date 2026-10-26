# Unit tests for line segments:
#   src/geometry/segments/linesegment.jl      (LineSegment constructor, curve, domain_fun, arc_length(line, pt))
#   src/geometry/curvederivatives.jl          (tangent, tangent_2 of an AbsLine, scalar and array methods)
#   src/geometry/arclength.jl                 (generic arc_length(crv, t) and arc_length(crv, ts) applied to a line)
#   src/geometry/geometry.jl                  (is_inside and domain_gradient_vector of a single curve)
# Reference object: the chord (0,0) -> (3,4) of length 5 and the diagonal (0,0) -> (1,1). A line is affine in t, so every
# quantity has a closed form: curve = pt0 + t (pt1 - pt0), tangent = pt1 - pt0, tangent_2 = 0, arc length = t L, and the
# signed domain function of the diagonal is (x - y) * orientation. No eigen solves, no pinned wavefunction values.
# Helper constants use the prefix `ls_`.

using StaticArrays: SVector
using LinearAlgebra: norm

# Absolute tolerance for quadrature-based arc lengths (quadgk with its default rtol = sqrt(eps) of the element type).
const ls_ARCLEN_ATOL = 1e-7

@testset "LineSegment constructor" begin
    crv = LineSegment([0.0, 0.0], [3.0, 4.0])
    @test crv isa LineSegment{Float64,SpecularReflection}
    @test crv.pt0 == SVector(0.0, 0.0)
    @test crv.pt1 == SVector(3.0, 4.0)
    @test crv.length == 5.0
    @test crv.orientation == 1
    @test crv.domain_id == 1
    @test crv.segment_id == 1
    @test crv.bc isa SpecularReflection

    tagged = LineSegment([0.0, 0.0], [1.0, 0.0]; bc=SymmetryWall(2, 3), orientation=-1, domain_id=4, segment_id=5)
    @test tagged.bc == SymmetryWall(2, 3)
    @test tagged.orientation == -1
    @test tagged.domain_id == 4
    @test tagged.segment_id == 5

    # mixed element types promote to a single T
    @test LineSegment([0.0, 0.0], Float32[1, 1]) isa LineSegment{Float64}
    @test LineSegment(Float32[0, 0], Float32[3, 4]) isa LineSegment{Float32}
    @test LineSegment(Float32[0, 0], Float32[3, 4]).length === 5.0f0
end

@testset "curve(::LineSegment)" begin
    pt0, pt1 = [0.0, 0.0], [1.0, 1.0]
    crv = LineSegment(pt0, pt1)
    @test curve(crv, 0.0) == pt0
    @test curve(crv, 1.0) == pt1
    @test curve(crv, 0.25) ≈ SVector(0.25, 0.25)

    # array method evaluates the scalar method elementwise
    ts = [0.0, 0.25, 0.5, 1.0]
    pts = curve(crv, ts)
    @test pts isa Vector{<:SVector{2,Float64}}
    @test pts == [curve(crv, t) for t in ts]
    @test length(curve(crv, Float64[])) == 0

    # extrapolation outside [0,1] stays on the supporting line
    @test curve(crv, 2.0) ≈ SVector(2.0, 2.0)
    @test curve(crv, -1.0) ≈ SVector(-1.0, -1.0)
end

@testset "tangent and tangent_2 (::LineSegment)" begin
    crv = LineSegment([0.0, 0.0], [2.0, -1.0])
    chord = crv.pt1 .- crv.pt0
    @test tangent(crv, 0.0) == chord
    @test tangent(crv, 0.7) == chord
    @test tangent_2(crv, 0.4) == zero(crv.pt0)

    ts = [0.1, 0.4, 0.9]
    @test tangent(crv, ts) == [tangent(crv, t) for t in ts]
    @test tangent(crv, ts) == [chord, chord, chord]
    @test tangent_2(crv, ts) == [tangent_2(crv, t) for t in ts]
    @test all(iszero, tangent_2(crv, ts))
    @test norm(tangent(crv, 0.3)) ≈ crv.length
end

@testset "arc_length(::LineSegment, t) and arc_length(::LineSegment, pt)" begin
    crv = LineSegment([0.0, 0.0], [3.0, 4.0])
    L = crv.length

    # point method is the Euclidean distance from pt0
    @test arc_length(crv, crv.pt1) == L
    @test arc_length(crv, crv.pt0) == 0.0
    @test arc_length(crv, curve(crv, 0.3)) ≈ 0.3 * L

    # parameter method integrates |tangent| = L, so arc_length(crv, t) = t L
    @test arc_length(crv, 1.0) ≈ L atol = ls_ARCLEN_ATOL
    @test arc_length(crv, 0.25) ≈ 0.25 * L atol = ls_ARCLEN_ATOL
    @test arc_length(crv, 0.0) == 0.0
    ts = [0.0, 0.2, 0.5, 1.0]
    @test arc_length(crv, ts) ≈ ts .* L atol = ls_ARCLEN_ATOL
    # both methods agree on the same point of the curve
    @test arc_length(crv, 0.6) ≈ arc_length(crv, curve(crv, 0.6)) atol = ls_ARCLEN_ATOL
end

@testset "domain_fun and is_inside (::LineSegment)" begin
    crv = LineSegment([0.0, 0.0], [1.0, 1.0])
    # gradient of the signed domain function x - y is the constant normal (1, -1)
    @test domain_gradient_vector(crv, curve(crv, 0.5)) == SVector{2,Float64}([1.0, -1.0])

    testpt1 = SVector{2,Float64}([0.1, 0.2])
    testpt2 = SVector{2,Float64}([0.1, 0.0])
    @test domain_fun(crv, testpt1) < 0.0
    @test domain_fun(crv, testpt2) > 0.0
    @test domain_fun(crv, testpt1) ≈ -0.1
    @test domain_fun(crv, testpt2) ≈ 0.1
    @test is_inside(crv, [testpt1, testpt2]) == [true, false]
    @test is_inside(crv, testpt1)
    @test !is_inside(crv, testpt2)

    # array methods
    @test domain_fun(crv, [testpt1, testpt2]) ≈ [-0.1, 0.1]
    @test domain_gradient_vector(crv, [testpt1, testpt2]) == [SVector(1.0, -1.0), SVector(1.0, -1.0)]
    # points on the line are on the boundary: zero domain function, counted as inside
    @test domain_fun(crv, SVector(0.4, 0.4)) == 0.0
    @test is_inside(crv, SVector(0.4, 0.4))
end

@testset "LineSegment orientation sign" begin
    ccw = LineSegment([0.0, 0.0], [1.0, 1.0]; orientation=1)
    cw = LineSegment([0.0, 0.0], [1.0, 1.0]; orientation=-1)
    testpt1 = SVector(0.1, 0.2)
    testpt2 = SVector(0.1, 0.0)
    # orientation -1 flips the sign of the domain function and therefore the inside half-plane
    @test domain_fun(cw, testpt1) == -domain_fun(ccw, testpt1)
    @test domain_fun(cw, testpt2) == -domain_fun(ccw, testpt2)
    @test is_inside(cw, [testpt1, testpt2]) == [false, true]
    @test domain_gradient_vector(cw, testpt1) == -domain_gradient_vector(ccw, testpt1)
    # the geometry itself does not depend on orientation
    @test curve(cw, 0.3) == curve(ccw, 0.3)
    @test tangent(cw, 0.3) == tangent(ccw, 0.3)
    @test cw.length == ccw.length
end

@testset "LineSegment Float32 element type" begin
    crv = LineSegment(Float32[0, 0], Float32[3, 4])
    t = 0.5f0
    pt = SVector(1.0f0, 2.0f0)
    @test @inferred(curve(crv, t)) isa SVector{2,Float32}
    @test @inferred(curve(crv, t)) == SVector(1.5f0, 2.0f0)
    @test @inferred(tangent(crv, t)) == SVector(3.0f0, 4.0f0)
    @test @inferred(tangent_2(crv, t)) == SVector(0.0f0, 0.0f0)
    @test @inferred(domain_fun(crv, pt)) isa Float32
    @test @inferred(arc_length(crv, curve(crv, t))) ≈ 2.5f0
    @test @inferred(arc_length(crv, t)) isa Float32
    # Float32 quadgk uses rtol = sqrt(eps(Float32)) ~ 3.4e-4
    @test arc_length(crv, t) ≈ 2.5f0 rtol = 1.0f-3
    @test eltype(curve(crv, Float32[0, 1])) == SVector{2,Float32}
end
