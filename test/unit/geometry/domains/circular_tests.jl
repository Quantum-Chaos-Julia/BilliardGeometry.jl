# Unit tests for src/geometry/domains/circular.jl: the CircleWedge domain type and its constructor.
# Reference wedges have corners at center + R (cos phi, sin phi) for phi the start and end angle of the arc; helpers use the prefix `cw_`.

using StaticArrays: SVector

# Absolute tolerance for corner coordinates (worst observed 1.1e-16, margin 9000x).
const cw_POINT_ATOL = 1e-12
# Tolerances for Float32 coordinates and lengths (worst observed about 1 ulp = 6e-8 relative, 1.2e-7 absolute, margin 8x).
const cw_RTOL_F32 = 1e-6
const cw_ATOL_F32 = 1e-6

const cw_BCS = AbsBoundaryCondition[SpecularReflection(), Transparent(1), SymmetryWall(3, 4)]

# Quarter wedge of the unit circle centred at (0.5, 0): the arc runs from (1.5, 0) to (0.5, 1).
cw_quarter(id = 2) = CircleWedge(1.0, pi / 2, 0.0, 0.5, 0.0, id; bcs = cw_BCS)

@testset "CircleWedge fields" begin
    w = cw_quarter()
    @test w isa CircleWedge{Float64}
    @test w isa AbsSimpleDomain
    @test fieldnames(CircleWedge) == (:boundary, :corners, :center, :id)
    @test w.id == 2
    @test w.center == SVector(0.5, 0.0)
    @test w.boundary isa Vector{AbsCurve}
    @test w.corners isa Vector{SVector{2,Float64}}
    @test length(w.boundary) == 3
    @test length(w.corners) == 3

    # corners: start of the arc, end of the arc, centre
    @test w.corners[1] ≈ SVector(1.5, 0.0) atol = cw_POINT_ATOL
    @test w.corners[2] ≈ SVector(0.5, 1.0) atol = cw_POINT_ATOL
    @test w.corners[3] == w.center
end

@testset "CircleWedge boundary curves" begin
    w = cw_quarter()
    arc, radial_out, radial_in = w.boundary
    @test arc isa CircleSegment{Float64}
    @test radial_out isa LineSegment{Float64}
    @test radial_in isa LineSegment{Float64}

    # the arc carries the constructor data
    @test arc.radius == 1.0
    @test arc.arc_angle ≈ pi / 2
    @test arc.shift_angle == 0.0
    @test arc.center == SVector(0.5, 0.0)

    # segments are numbered 1, 2, 3 and all carry the domain id
    @test [c.segment_id for c in w.boundary] == [1, 2, 3]
    @test all(c -> c.domain_id == 2, w.boundary)
    @test cw_quarter(7).boundary[1].domain_id == 7

    # the boundary is a loop: arc end -> centre -> arc start -> arc end
    @test radial_out.pt1 == w.center
    @test radial_in.pt0 == w.center
    @test radial_out.pt0 ≈ w.corners[2] atol = cw_POINT_ATOL
    @test radial_in.pt1 ≈ w.corners[1] atol = cw_POINT_ATOL
    @test all(is_connected(w.boundary))

    # lengths: arc R * arc_angle, two radii R; perimeter R (2 + arc_angle)
    @test arc.length ≈ pi / 2
    @test radial_out.length ≈ 1.0
    @test radial_in.length ≈ 1.0
    @test sum(c.length for c in w.boundary) ≈ 2.0 + pi / 2

    # a second wedge with other radius, opening, orientation and centre: corners follow from the closed forms
    R, a, s, x0, y0 = 2.0, pi / 3, pi / 6, -1.0, 0.5
    w2 = CircleWedge(R, a, s, x0, y0, 1; bcs = cw_BCS)
    @test w2.corners[1] ≈ SVector(x0 + R * cos(s), y0 + R * sin(s)) atol = cw_POINT_ATOL
    @test w2.corners[2] ≈ SVector(x0 + R * cos(s + a), y0 + R * sin(s + a)) atol = cw_POINT_ATOL
    @test w2.center == SVector(x0, y0)
    @test w2.boundary[1].length ≈ R * a
    @test all(is_connected(w2.boundary))
end

@testset "CircleWedge boundary conditions" begin
    # the three conditions go to the arc, the radius leaving the arc end and the radius returning to the arc start
    w = cw_quarter()
    @test [c.bc for c in w.boundary] == cw_BCS
    @test w.boundary[1] isa CircleSegment{Float64,SpecularReflection}
    @test w.boundary[2] isa LineSegment{Float64,Transparent}
    @test w.boundary[3] isa LineSegment{Float64,SymmetryWall}

    # the keyword is typed as Vector{AbsBoundaryCondition}, so a vector of a single concrete condition type is rejected
    @test_throws TypeError CircleWedge(1.0, pi / 2, 0.0, 0.5, 0.0, 2; bcs = [Transparent(1), Transparent(2), Transparent(3)])
    # fewer than three conditions cannot be matched to the three curves
    @test_throws BoundsError CircleWedge(1.0, pi / 2, 0.0, 0.5, 0.0, 2; bcs = AbsBoundaryCondition[SpecularReflection()])
end

@testset "CircleWedge element type" begin
    w32 = CircleWedge(1.0f0, Float32(pi / 2), 0.0f0, 0.5f0, 0.0f0, 2; bcs = cw_BCS)
    @test w32 isa CircleWedge{Float32}
    @test w32.corners isa Vector{SVector{2,Float32}}
    @test w32.center isa SVector{2,Float32}
    @test w32.boundary[1] isa CircleSegment{Float32}
    @test w32.corners[1] ≈ SVector(1.5f0, 0.0f0) atol = cw_ATOL_F32
    @test w32.corners[2] ≈ SVector(0.5f0, 1.0f0) atol = cw_ATOL_F32
    @test w32.boundary[1].length ≈ Float32(pi / 2) rtol = cw_RTOL_F32
end
