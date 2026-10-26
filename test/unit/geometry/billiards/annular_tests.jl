# Unit tests for the annular billiard constructor (structural checks only):
#   src/geometry/billiards/annular.jl   (AnnularBilliard{T}, AnnularBilliard(R_inner, shift; R_outer, center))
# Reference object: a MultiplyConnectedDomain with id 1, genus 1 and no symmetries, made of an outer full CircleSegment (orientation 1,
# domain id 1) and an inner full CircleSegment (orientation -1, domain id 2) whose centre is displaced by `shift` along the x axis.
# Only structure is checked: curve counts, orientations, ids, genus and components. Guards: R_outer > 0, R_inner > 0, R_inner < R_outer.

using StaticArrays: SVector

@testset "AnnularBilliard fields" begin
    ann = AnnularBilliard(0.3, 0.2)
    @test ann isa AnnularBilliard{Float64}
    @test ann isa AbsBilliard
    @test ann.fundamental_domain isa MultiplyConnectedDomain{Float64}
    @test ann.symmetries isa SymmetryRegistry
    @test fieldnames(AnnularBilliard) == (:fundamental_domain, :symmetries)
    @test @inferred(AnnularBilliard(0.3, 0.2)) isa AnnularBilliard{Float64}
end

@testset "AnnularBilliard symmetries" begin
    @test isempty(AnnularBilliard(0.3, 0.2).symmetries)
end

@testset "AnnularBilliard curves" begin
    dom = AnnularBilliard(0.3, 0.2).fundamental_domain
    @test dom.id == 1
    @test length(dom.boundary) == 2
    outer, inner = dom.boundary
    @test outer isa CircleSegment{Float64}
    @test inner isa CircleSegment{Float64}
    @test all(c -> c.bc isa SpecularReflection, dom.boundary)

    # the hole curve is stored with orientation -1 and its own domain id
    @test [c.orientation for c in dom.boundary] == [1, -1]
    @test [c.domain_id for c in dom.boundary] == [1, 2]
    @test [c.segment_id for c in dom.boundary] == [1, 1]
    @test outer.arc_angle == 2 * pi
    @test inner.arc_angle == 2 * pi

    # outer circle of radius one about the origin, inner circle of radius R_inner about (shift, 0)
    @test outer.radius == 1.0
    @test outer.center == SVector(0.0, 0.0)
    @test inner.radius == 0.3
    @test inner.center == SVector(0.2, 0.0)
end

@testset "AnnularBilliard topology" begin
    ann = AnnularBilliard(0.3, 0.2)
    @test genus(ann) == 1
    @test genus(ann.fundamental_domain) == 1
    @test ann.fundamental_domain.genus == 1
    # one component for the outer circle and one for the hole, one curve each
    comps = boundary_components(ann)
    @test length(comps) == genus(ann) + 1
    @test [length(c) for c in comps] == [1, 1]
    @test [c[1].domain_id for c in comps] == [1, 2]
    @test length(geometric_boundary_components(ann)) == genus(ann) + 1
end

@testset "AnnularBilliard keywords" begin
    dom = AnnularBilliard(0.3, 0.2; R_outer=2.0, center=[0.1, 0.0]).fundamental_domain
    @test dom.boundary[1].radius == 2.0
    @test dom.boundary[1].center == SVector(0.1, 0.0)
    @test dom.boundary[2].radius == 0.3
    @test dom.boundary[2].center ≈ SVector(0.3, 0.0) atol = 1e-14
end

@testset "AnnularBilliard argument guards" begin
    # R_outer > 0
    @test_throws ArgumentError AnnularBilliard(0.3, 0.2; R_outer=0.0)
    @test_throws ArgumentError AnnularBilliard(0.3, 0.2; R_outer=-1.0)
    # R_inner > 0
    @test_throws ArgumentError AnnularBilliard(0.0, 0.2)
    @test_throws ArgumentError AnnularBilliard(-0.3, 0.2)
    # R_inner < R_outer
    @test_throws ArgumentError AnnularBilliard(1.0, 0.0)
    @test_throws ArgumentError AnnularBilliard(1.5, 0.0)
    @test_throws ArgumentError AnnularBilliard(0.5, 0.0; R_outer=0.5)
end

@testset "AnnularBilliard eltype promotion" begin
    ann32 = AnnularBilliard(0.3f0, 0.2f0)
    @test ann32 isa AnnularBilliard{Float32}
    @test ann32.fundamental_domain isa MultiplyConnectedDomain{Float32}
    @test all(c -> c isa CircleSegment{Float32}, ann32.fundamental_domain.boundary)
    @test genus(ann32) == 1
end
