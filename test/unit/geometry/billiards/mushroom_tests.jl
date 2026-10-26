# Unit tests for the mushroom billiard constructor (structural checks only):
#   src/geometry/billiards/mushroom.jl   (MushroomBilliard{T}, MushroomBilliard(half_width, stem_heigth; R, origin))
# Reference object: the half mushroom with cap radius R, stem half-width hw < R and stem height h, a composite domain of three
# subdomains (cap sector with id 1, triangle with id 2, stem rectangle with id 3) joined by Transparent curves, one YAxisReflection
# symmetry (sym_id 1), and the constructor guard that rejects a cap centre away from the origin because the YAxisReflection
# requires it. No area, spectra or boundary samples are computed, so no eigen solves and no pinned wavefunction values. Helper
# constants use the prefix `mb_`.

using StaticArrays: SVector

const mb_HW = 0.5
const mb_STEM = 2.0

@testset "MushroomBilliard fields" begin
    mush = MushroomBilliard(mb_HW, mb_STEM)
    @test mush isa MushroomBilliard{Float64}
    @test mush isa AbsBilliard
    @test mush.fundamental_domain isa CompositeDomain
    @test mush.symmetries isa SymmetryRegistry
    @test fieldnames(MushroomBilliard) == (:fundamental_domain, :symmetries)
    # the stem height defaults to one
    @test MushroomBilliard(mb_HW).fundamental_domain.subdomains[3].corners[2] == SVector(0.0, -1.0)
end

@testset "MushroomBilliard symmetries" begin
    reg = MushroomBilliard(mb_HW, mb_STEM).symmetries
    @test length(reg) == 1
    @test reg[1] isa YAxisReflection
    @test reg[1].sym_id == 1
    @test symmetry_of(reg, 1) isa YAxisReflection
end

@testset "MushroomBilliard subdomains" begin
    subs = MushroomBilliard(mb_HW, mb_STEM).fundamental_domain.subdomains
    @test length(subs) == 3
    @test [d.id for d in subs] == [1, 2, 3]
    @test subs[1] isa SimpleDomain{Float64}
    @test subs[2] isa Polygon{Float64}
    @test subs[3] isa Polygon{Float64}
    @test [length(d.boundary) for d in subs] == [3, 3, 4]
    for d in subs
        @test [c.segment_id for c in d.boundary] == collect(1:length(d.boundary))
        @test all(c -> c.domain_id == d.id, d.boundary)
    end

    # boundary conditions: the cap and the stem walls are physical, the cuts between subdomains are Transparent
    @test [c.bc for c in subs[1].boundary] == [SpecularReflection(), SpecularReflection(), Transparent(2)]
    @test [c.bc for c in subs[2].boundary] == [SymmetryWall(1, 2), Transparent(3), Transparent(1)]
    @test [c.bc for c in subs[3].boundary] == [SymmetryWall(1, 2), SpecularReflection(), SpecularReflection(), Transparent(2)]

    # geometry: cap sector corners, triangle corners and stem corners
    @test subs[1].corners == [SVector(1.0, 0.0), SVector(0.0, 1.0), SVector(mb_HW, 0.0)]
    @test subs[1].boundary[2] isa CircleSegment
    @test subs[1].boundary[2].radius == 1.0
    @test subs[2].corners == [SVector(0.0, 1.0), SVector(0.0, 0.0), SVector(mb_HW, 0.0)]
    @test subs[3].corners == [SVector(0.0, 0.0), SVector(0.0, -mb_STEM), SVector(mb_HW, -mb_STEM), SVector(mb_HW, 0.0)]
end

@testset "MushroomBilliard cap radius keyword" begin
    subs = MushroomBilliard(mb_HW, mb_STEM; R=2.0).fundamental_domain.subdomains
    @test subs[1].boundary[2].radius == 2.0
    @test subs[1].corners[1] == SVector(2.0, 0.0)
    @test subs[1].corners[2] == SVector(0.0, 2.0)
    @test subs[2].corners[1] == SVector(0.0, 2.0)
end

@testset "MushroomBilliard origin guard" begin
    # the YAxisReflection symmetry requires the cap centre at the origin
    @test_throws ArgumentError MushroomBilliard(mb_HW, mb_STEM; origin=[0.1, 0.0])
    @test_throws ArgumentError MushroomBilliard(mb_HW, mb_STEM; origin=[0.0, 0.2])
    @test_throws ArgumentError MushroomBilliard(mb_HW; origin=[-0.3, 0.4])
    # an explicit zero origin is accepted and gives the default billiard
    @test MushroomBilliard(mb_HW, mb_STEM; origin=[0.0, 0.0]) isa MushroomBilliard{Float64}
end

@testset "MushroomBilliard eltype promotion" begin
    mush32 = MushroomBilliard(0.5f0, 2.0f0)
    @test mush32 isa MushroomBilliard{Float32}
    subs = mush32.fundamental_domain.subdomains
    @test subs[1] isa SimpleDomain{Float32}
    @test subs[2] isa Polygon{Float32}
    @test subs[3] isa Polygon{Float32}
    @test subs[1].boundary[2].radius isa Float32
    @test eltype(subs[3].corners) == SVector{2,Float32}
    # keyword arguments given in Float64 are converted to the element type of the half-width
    mush32 = MushroomBilliard(0.5f0, 2.0f0; R=1.5, origin=[0.0, 0.0])
    @test mush32 isa MushroomBilliard{Float32}
    @test mush32.fundamental_domain.subdomains[1].boundary[2].radius === 1.5f0
    # the default stem height of the one-argument method takes the element type of the half-width
    @test MushroomBilliard(0.5f0) isa MushroomBilliard{Float32}
end
