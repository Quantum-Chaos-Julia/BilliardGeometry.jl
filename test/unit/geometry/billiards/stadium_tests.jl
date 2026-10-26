# Unit tests for the stadium billiard constructor (structural checks only):
#   src/geometry/billiards/stadium.jl   (StadiumBilliard{T}, StadiumBilliard(half_width))
# Reference object: the quarter stadium of radius 1 and half-width hw, a composite domain of two subdomains, the circle wedge with
# id 1 centred at (hw, 0) and the rectangle with id 2 spanning [0, hw] x [0, 1], joined by a Transparent curve and closed by
# symmetry walls; D2 symmetry registry (YAxisReflection, XYAxisReflection, XAxisReflection with sym_ids 1, 2, 3). No area, spectra
# or boundary samples are computed, so no eigen solves and no pinned wavefunction values. Helper constants use the prefix `sb_`.

using StaticArrays: SVector

const sb_HW = 0.7

@testset "StadiumBilliard fields" begin
    stad = StadiumBilliard(sb_HW)
    @test stad isa StadiumBilliard{Float64}
    @test stad isa AbsBilliard
    @test stad.fundamental_domain isa CompositeDomain
    @test stad.symmetries isa SymmetryRegistry
    @test fieldnames(StadiumBilliard) == (:fundamental_domain, :symmetries)
end

@testset "StadiumBilliard symmetries" begin
    reg = StadiumBilliard(sb_HW).symmetries
    @test length(reg) == 3
    @test reg[1] isa YAxisReflection
    @test reg[2] isa XYAxisReflection
    @test reg[3] isa XAxisReflection
    @test [s.sym_id for s in reg] == [1, 2, 3]
    @test symmetry_of(reg, 2) isa XYAxisReflection
    # the registry equals the D2 registry used by the other D2-symmetric billiards
    @test reg == D2_symmetry()
end

@testset "StadiumBilliard subdomains" begin
    subs = StadiumBilliard(sb_HW).fundamental_domain.subdomains
    @test length(subs) == 2
    @test [d.id for d in subs] == [1, 2]
    @test subs[1] isa CircleWedge{Float64}
    @test subs[2] isa Polygon{Float64}
    @test [length(d.boundary) for d in subs] == [3, 4]
    # segment ids number the curves of each subdomain from 1, and every curve carries its subdomain id
    for d in subs
        @test [c.segment_id for c in d.boundary] == collect(1:length(d.boundary))
        @test all(c -> c.domain_id == d.id, d.boundary)
    end

    # boundary conditions: the semicircle cap and the straight wall are physical, the rest are cuts
    @test [c.bc for c in subs[1].boundary] == [SpecularReflection(), Transparent(2), SymmetryWall(3, 4)]
    @test [c.bc for c in subs[2].boundary] == [SpecularReflection(), SymmetryWall(1, 4), SymmetryWall(3, 4), Transparent(1)]

    # geometry of the corners: the wedge is centred at (hw, 0), the rectangle spans [0, hw] x [0, 1]
    @test subs[1].center == SVector(sb_HW, 0.0)
    @test subs[1].boundary[1] isa CircleSegment
    @test subs[1].boundary[1].radius == 1.0
    @test subs[2].corners == [SVector(sb_HW, 1.0), SVector(0.0, 1.0), SVector(0.0, 0.0), SVector(sb_HW, 0.0)]
    # the cap arc runs from (hw + 1, 0) to (hw, 1), the top of the Transparent joint with the rectangle
    @test curve(subs[1].boundary[1], 0.0) ≈ SVector(sb_HW + 1.0, 0.0) atol = 1e-14
    @test curve(subs[1].boundary[1], 1.0) ≈ SVector(sb_HW, 1.0) atol = 1e-14
end

@testset "StadiumBilliard eltype promotion" begin
    stad32 = StadiumBilliard(0.7f0)
    @test stad32 isa StadiumBilliard{Float32}
    subs = stad32.fundamental_domain.subdomains
    @test subs[1] isa CircleWedge{Float32}
    @test subs[2] isa Polygon{Float32}
    @test subs[1].boundary[1].radius isa Float32
    @test eltype(subs[2].corners) == SVector{2,Float32}
    @test subs[2].boundary[1] isa LineSegment{Float32}
end
