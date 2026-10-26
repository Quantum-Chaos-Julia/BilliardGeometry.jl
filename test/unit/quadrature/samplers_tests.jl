# Unit tests for the one-dimensional node samplers on the unit parameter interval:
#   src/quadrature/samplers.jl   (LinearNodes, GaussLegendreNodes, FourierNodes and their sample_points methods)
# Reference objects: a sampler returns nodes t_i in [0, 1] and weights dt_i, so that sum_i f(t_i) dt_i approximates the
# integral of f over [0, 1]. LinearNodes is the midpoint rule on N equal cells (exact for degree <= 1), GaussLegendreNodes
# the N-point Gauss-Legendre rule mapped to [0, 1] (exact for degree <= 2N-1). FourierNodes lays out equispaced nodes
# over consecutive curves of the given lengths, the node spacing along the concatenated boundary being L / M with
# L = sum(lengths) and M = N, or M = nextprod(primes, N) when `primes` is given; the weights dts[i] are in units of the
# parameter of curve i, so dts[i] * lengths[i] is arc length. No eigen solves, no pinned wavefunction values. Helper
# constants use the prefix `sm_`.

using BilliardGeometry: LinearNodes, GaussLegendreNodes, FourierNodes, sample_points

# Absolute tolerance for quadrature sums of polynomials (a few hundred ulps of Float64).
const sm_ATOL = 1e-13
# Absolute tolerance for Float32 comparisons.
const sm_ATOL_F32 = 1f-6

@testset "LinearNodes" begin
    @test LinearNodes() isa LinearNodes{Float64}
    @test LinearNodes(Float32) isa LinearNodes{Float32}
    @test LinearNodes() isa AbsSampler

    for N in (1, 4, 10)
        t, dt = sample_points(LinearNodes(), N)
        @test length(t) == N && length(dt) == N
        # midpoints of N equal cells
        @test collect(t) ≈ [(i - 0.5) / N for i in 1:N] atol = sm_ATOL
        @test all(dt .≈ 1 / N)
        @test sum(dt) ≈ 1 atol = sm_ATOL
        # the midpoint rule is exact for degree <= 1
        @test sum(dt .* t) ≈ 1 / 2 atol = sm_ATOL
    end

    # the midpoint rule error for t^2 is h^2 / 12 = 1 / (12 N^2) with the exact value 1/3
    t, dt = sample_points(LinearNodes(), 5)
    @test sum(dt .* t .^ 2) ≈ 1 / 3 - 1 / (12 * 5^2) atol = sm_ATOL

    # element type follows the sampler
    t32, dt32 = sample_points(LinearNodes(Float32), 4)
    @test eltype(t32) == Float32 && eltype(dt32) == Float32
    @test collect(t32) ≈ Float32[0.125, 0.375, 0.625, 0.875] atol = sm_ATOL_F32
end

@testset "GaussLegendreNodes" begin
    @test GaussLegendreNodes() isa GaussLegendreNodes{Float64}
    @test GaussLegendreNodes(Float32) isa GaussLegendreNodes{Float32}
    @test GaussLegendreNodes() isa AbsSampler

    for N in (1, 3, 5, 12)
        t, dt = sample_points(GaussLegendreNodes(), N)
        @test length(t) == N && length(dt) == N
        @test all(0 .< t .< 1)
        @test issorted(t)
        @test all(dt .> 0)
        # weights sum to the length of [0, 1]
        @test sum(dt) ≈ 1 atol = sm_ATOL
        # exact for every monomial up to degree 2N-1
        @test all(p -> isapprox(sum(dt .* t .^ p), 1 / (p + 1); atol = sm_ATOL), 0:(2*N-1))
    end

    # one point: the midpoint with unit weight
    t, dt = sample_points(GaussLegendreNodes(), 1)
    @test t ≈ [0.5] atol = sm_ATOL
    @test dt ≈ [1.0] atol = sm_ATOL

    # nodes are symmetric about 1/2 and the weights are palindromic
    t, dt = sample_points(GaussLegendreNodes(), 6)
    @test t ≈ reverse(1 .- t) atol = sm_ATOL
    @test dt ≈ reverse(dt) atol = sm_ATOL

    # three-point rule: nodes 1/2 -+ sqrt(3/5)/2, weights 5/18, 4/9, 5/18
    t, dt = sample_points(GaussLegendreNodes(), 3)
    @test t ≈ [0.5 - sqrt(3 / 5) / 2, 0.5, 0.5 + sqrt(3 / 5) / 2] atol = sm_ATOL
    @test dt ≈ [5 / 18, 4 / 9, 5 / 18] atol = sm_ATOL

    # element type follows the sampler
    t32, dt32 = sample_points(GaussLegendreNodes(Float32), 3)
    @test eltype(t32) == Float32 && eltype(dt32) == Float32
    @test sum(dt32) ≈ 1f0 atol = sm_ATOL_F32
end

@testset "FourierNodes" begin
    @test FourierNodes([1.0, 0.5]) isa FourierNodes{Float64}
    @test FourierNodes([1.0, 0.5]).primes === nothing
    @test FourierNodes([1.0, 0.5]) isa AbsSampler

    @testset "weights across curve junctions" begin
        for lengths in ([1.0, 0.5], [0.3, 1.1, 0.7], [2.0, 1.0, 1.0, 0.25]), N in (16, 64, 512)
            ts, dts = sample_points(FourierNodes(lengths), N)
            L = sum(lengths)
            weights = reduce(vcat, [dts[i] .* lengths[i] for i in eachindex(lengths)])
            @test sum(weights) ≈ L rtol = 1e-12
            @test all(w -> isapprox(w, L / N; rtol = 1e-9), weights[1:end-1])
            @test all(i -> length(ts[i]) == length(dts[i]), eachindex(lengths))
        end
    end

    @testset "node layout" begin
        lengths = [0.3, 1.1, 0.7]
        N = 64
        L = sum(lengths)
        ts, dts = sample_points(FourierNodes(lengths), N)
        @test length(ts) == 3 && length(dts) == 3
        # every curve starts at or after t = 0 and stays in [0, 1]
        @test ts[1][1] == 0.0
        @test all(t -> all(0 .<= t .<= 1 + 1e-12), ts)
        # equal arc-length spacing L / N along the whole concatenated boundary
        @test all(i -> all(isapprox.(diff(ts[i]) .* lengths[i], L / N; rtol = 1e-9)), eachindex(lengths))
        # across a junction the arc length from the last node of a curve to the first node of the next is L / N
        for i in 1:2
            gap = (1 - ts[i][end]) * lengths[i] + ts[i+1][1] * lengths[i+1]
            @test gap ≈ L / N rtol = 1e-9
        end
    end

    @testset "primes" begin
        # N already a product of the given primes: identical to the unrestricted sampler
        lengths = [1.0, 0.5]
        @test sample_points(FourierNodes([2, 3, 5], lengths), 16) == sample_points(FourierNodes(lengths), 16)

        # otherwise M = nextprod(primes, N) nodes per boundary length L, here N = 17 -> M = 18
        @test nextprod([2, 3, 5], 17) == 18
        for lengths in ([1.0], [1.0, 0.5], [0.3, 1.1, 0.7])
            L = sum(lengths)
            M = 18
            ts, dts = sample_points(FourierNodes([2, 3, 5], lengths), 17)
            weights = reduce(vcat, [dts[i] .* lengths[i] for i in eachindex(lengths)])
            @test sum(weights) ≈ L rtol = 1e-12
            @test all(w -> isapprox(w, L / M; rtol = 1e-9), weights[1:end-1])
            @test all(i -> all(isapprox.(diff(ts[i]) .* lengths[i], L / M; rtol = 1e-9)), eachindex(lengths))
        end

        # a prime list that does not contain 2 can force a larger M: nextprod([3], 17) = 27
        ts, dts = sample_points(FourierNodes([3], [1.0]), 17)
        @test nextprod([3], 17) == 27
        @test all(w -> isapprox(w, 1 / 27; rtol = 1e-9), dts[1][1:end-1])
    end
end
