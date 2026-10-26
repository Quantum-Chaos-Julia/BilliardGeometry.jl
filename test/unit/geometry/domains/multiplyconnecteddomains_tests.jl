# Unit tests for src/geometry/domains/multiplyconnecteddomains.jl: the MultiplyConnectedDomain type and its constructor.
# Fixtures are an outer boundary of two line pieces and hole curve groups of full circles stored with orientation -1; helpers use the prefix `mc_`.

using StaticArrays: SVector

mc_outer() = [LineSegment([0.0, 0.0], [2.0, 0.0]), LineSegment([2.0, 0.0], [2.0, 2.0])]
mc_hole(R, x0, y0, id) = CircleSegment(R, 2pi, 0.0, [x0, y0]; orientation = -1, domain_id = id)
mc_corners() = [SVector(0.0, 0.0), SVector(2.0, 0.0)]

@testset "MultiplyConnectedDomain fields" begin
    outer = mc_outer()
    hole_a = AbsCurve[mc_hole(0.3, 1.0, 1.0, 2)]
    hole_b = AbsCurve[mc_hole(0.1, 0.5, 0.5, 3), mc_hole(0.1, 1.5, 0.5, 3)]
    corners = mc_corners()
    dom = MultiplyConnectedDomain(outer, [hole_a, hole_b], corners, 4)

    @test dom isa MultiplyConnectedDomain{Float64}
    @test dom isa AbsMultiplyConnectedDomain
    @test !(dom isa AbsSimpleDomain)
    @test fieldnames(MultiplyConnectedDomain) == (:boundary, :corners, :id, :genus)
    @test dom.id == 4
    @test dom.corners == corners
    @test dom.genus == 2
    @test dom.boundary isa Vector{AbsCurve}

    # outer curves come first, then the curves of each hole group in order, all in one flat vector
    @test length(dom.boundary) == 5
    @test dom.boundary[1] === outer[1]
    @test dom.boundary[2] === outer[2]
    @test dom.boundary[3] === hole_a[1]
    @test dom.boundary[4] === hole_b[1]
    @test dom.boundary[5] === hole_b[2]
    @test [c.domain_id for c in dom.boundary] == [1, 1, 2, 3, 3]

    # hole curves keep their reversed orientation
    @test [c.orientation for c in dom.boundary] == [1, 1, -1, -1, -1]

    # the inputs are not modified
    @test length(outer) == 2
    @test length(hole_a) == 1
    @test length(hole_b) == 2
end

@testset "MultiplyConnectedDomain genus" begin
    outer = mc_outer()
    hole = AbsCurve[mc_hole(0.3, 1.0, 1.0, 2)]

    # genus is the number of hole groups, not of hole curves
    @test MultiplyConnectedDomain(outer, Vector{AbsCurve}[], mc_corners(), 1).genus == 0
    @test MultiplyConnectedDomain(outer, [hole], mc_corners(), 1).genus == 1
    @test MultiplyConnectedDomain(outer, [hole, hole], mc_corners(), 1).genus == 2
    two_curve_hole = AbsCurve[mc_hole(0.1, 0.5, 0.5, 2), mc_hole(0.1, 1.5, 0.5, 2)]
    @test MultiplyConnectedDomain(outer, [two_curve_hole], mc_corners(), 1).genus == 1

    # without holes the boundary is the outer curves only
    no_holes = MultiplyConnectedDomain(outer, Vector{AbsCurve}[], mc_corners(), 1)
    @test length(no_holes.boundary) == 2
    @test no_holes.boundary isa Vector{AbsCurve}
    @test MultiplyConnectedDomain(outer, Vector{Vector{AbsCurve}}(), mc_corners(), 1).boundary isa Vector{AbsCurve}
end

@testset "MultiplyConnectedDomain curve vector types" begin
    # outer curves and hole groups may be concretely or abstractly typed
    concrete_outer = mc_outer()
    @test concrete_outer isa Vector{<:LineSegment}
    abstract_outer = AbsCurve[concrete_outer...]
    concrete_hole = [mc_hole(0.3, 1.0, 1.0, 2)]
    @test concrete_hole isa Vector{<:CircleSegment}
    from_concrete = MultiplyConnectedDomain(concrete_outer, [concrete_hole], mc_corners(), 1)
    from_abstract = MultiplyConnectedDomain(abstract_outer, [AbsCurve[concrete_hole...]], mc_corners(), 1)
    @test from_concrete.genus == from_abstract.genus == 1
    @test length(from_concrete.boundary) == length(from_abstract.boundary) == 3
    @test from_concrete.boundary isa Vector{AbsCurve}
end

@testset "MultiplyConnectedDomain element type" begin
    # the numeric type is that of the corners
    outer32 = [LineSegment(Float32[0, 0], Float32[1, 0])]
    hole32 = AbsCurve[CircleSegment(0.1f0, 2.0f0 * Float32(pi), 0.0f0, Float32[0.5, 0.5]; orientation = -1, domain_id = 2)]
    corners32 = SVector{2,Float32}[SVector(0.0f0, 0.0f0)]
    dom32 = MultiplyConnectedDomain(outer32, [hole32], corners32, 1)
    @test dom32 isa MultiplyConnectedDomain{Float32}
    @test dom32.corners isa Vector{SVector{2,Float32}}
    @test dom32.genus == 1
    @test length(dom32.boundary) == 2
end

@testset "MultiplyConnectedDomain argument errors" begin
    outer = mc_outer()
    hole = AbsCurve[mc_hole(0.3, 1.0, 1.0, 2)]
    corners = mc_corners()

    # the id must be an Int64 and the corners a vector of static 2-vectors
    @test_throws MethodError MultiplyConnectedDomain(outer, [hole], corners, Int32(1))
    @test_throws MethodError MultiplyConnectedDomain(outer, [hole], corners, 1.0)
    @test_throws MethodError MultiplyConnectedDomain(outer, [hole], [[0.0, 0.0]], 1)

    # the outer curves are a vector and the holes a vector of curve vectors
    @test_throws MethodError MultiplyConnectedDomain(Tuple(outer), [hole], corners, 1)
    @test_throws MethodError MultiplyConnectedDomain(outer, hole, corners, 1)
end
