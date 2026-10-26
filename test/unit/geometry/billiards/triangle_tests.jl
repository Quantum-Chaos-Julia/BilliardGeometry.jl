# Unit tests for the triangle billiard constructor (structural checks only):
#   src/geometry/billiards/triangle.jl   (TriangleBilliard{T}, TriangleBilliard(gamma, chi; bcs, x0, y0, h))
# Reference object: the triangle with interior angles alpha = (pi - gamma)/(1 + chi), beta = chi*alpha and gamma, with the gamma
# corner at (-x0, -y0) and height h; no symmetries. No area, spectra or boundary samples are computed.

using StaticArrays: SVector

@testset "TriangleBilliard fields" begin
    tri = TriangleBilliard(pi / 2, 1.0)
    @test tri isa TriangleBilliard{Float64}
    @test tri isa AbsBilliard
    @test tri.fundamental_domain isa Polygon{Float64}
    @test tri.symmetries isa SymmetryRegistry
    @test fieldnames(TriangleBilliard) == (:fundamental_domain, :symmetries)
    @test @inferred(TriangleBilliard(1.0, 2.0)) isa TriangleBilliard{Float64}

    # three physical sides, ids 1..3 in the polygon with id 1
    dom = tri.fundamental_domain
    @test dom.id == 1
    @test length(dom.boundary) == 3
    @test all(c -> c isa LineSegment{Float64}, dom.boundary)
    @test [c.segment_id for c in dom.boundary] == [1, 2, 3]
    @test all(c -> c.domain_id == 1, dom.boundary)
    @test [c.bc for c in dom.boundary] == fill(SpecularReflection(), 3)
end

@testset "TriangleBilliard symmetries" begin
    @test isempty(TriangleBilliard(1.0, 2.0).symmetries)
end

@testset "TriangleBilliard corners" begin
    # gamma = pi/2, chi = 1 is the right isosceles triangle with the right angle at the origin and legs of length h = 1
    dom = TriangleBilliard(pi / 2, 1.0).fundamental_domain
    @test length(dom.corners) == 3
    @test dom.corners[1] ≈ SVector(1.0, 0.0) atol = 1e-14
    @test dom.corners[2] ≈ SVector(0.0, 1.0) atol = 1e-14
    @test dom.corners[3] ≈ SVector(0.0, 0.0) atol = 1e-14
    @test dom.angles ≈ [pi / 4, pi / 4, pi / 2] atol = 1e-14

    # the gamma corner sits at (-x0, -y0) and the angle there equals gamma
    gamma = 1.0
    shifted = TriangleBilliard(gamma, 2.0; x0=0.3, y0=-0.2, h=1.5).fundamental_domain
    @test shifted.corners[3] ≈ SVector(-0.3, 0.2) atol = 1e-14
    @test shifted.angles[3] ≈ gamma atol = 1e-14
    # the angles of a triangle sum to pi, alpha = (pi - gamma)/(1 + chi) and beta = chi*alpha
    @test sum(shifted.angles) ≈ pi atol = 1e-13
    @test shifted.angles[1] ≈ (pi - gamma) / 3 atol = 1e-13
    @test shifted.angles[2] ≈ 2 * (pi - gamma) / 3 atol = 1e-13
end

@testset "TriangleBilliard boundary conditions keyword" begin
    bcs = [SpecularReflection(), Transparent(2), SymmetryWall(1, 2)]
    dom = TriangleBilliard(1.0, 2.0; bcs=bcs).fundamental_domain
    @test [c.bc for c in dom.boundary] == bcs
end

@testset "TriangleBilliard argument guards" begin
    # gamma must lie strictly between 0 and pi
    @test_throws ArgumentError TriangleBilliard(0.0, 1.0)
    @test_throws ArgumentError TriangleBilliard(-0.1, 1.0)
    @test_throws ArgumentError TriangleBilliard(pi, 1.0)
    @test_throws ArgumentError TriangleBilliard(4.0, 1.0)
    # chi must be positive so that alpha and beta stay positive
    @test_throws ArgumentError TriangleBilliard(1.0, 0.0)
    @test_throws ArgumentError TriangleBilliard(1.0, -1.0)
end

@testset "TriangleBilliard eltype promotion" begin
    tri32 = TriangleBilliard(1.0f0, 2.0f0)
    @test tri32 isa TriangleBilliard{Float32}
    @test tri32.fundamental_domain isa Polygon{Float32}
    @test eltype(tri32.fundamental_domain.corners) == SVector{2,Float32}
    @test tri32.fundamental_domain.boundary[1] isa LineSegment{Float32}
    @test tri32.fundamental_domain.angles isa Vector{Float32}

    # mixed argument types promote to the wider type
    @test TriangleBilliard(1.0f0, 2.0) isa TriangleBilliard{Float64}
    @test TriangleBilliard(1.0, 2) isa TriangleBilliard{Float64}
end
