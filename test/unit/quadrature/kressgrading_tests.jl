# Unit tests for Kress logarithmic quadrature and boundary grading:
#   src/quadrature/kressgrading.jl   (kress_R_even!, kress_R_odd!, kress_R!, s_mid, kress_graded_nodes_data,
#                                     multi_kress_graded_nodes_data, _corner_interval)
# Reference objects: the periodic grid t_j = 2 pi j / N. The Kress matrix R is circulant with first column
# 2 pi R_j, where for N = 2n  R_j = -(1/n) sum_{m=1}^{n-1} cos(m j pi / n) / m - (-1)^j / (2 n^2)  and for N = 2n-1
# R_j = -(2/N) sum_{m=1}^{n-1} cos(2 pi m j / N) / m. It integrates the kernel log(4 sin^2((t - tau) / 2)) exactly
# against trigonometric polynomials: R * cos(m t) = -(2 pi / m) cos(m t), and in particular annihilates constants (m = 0).
# The graded maps t(sigma) are monotone, satisfy t(2 pi) - t(0) = 2 pi (so the weights h * t'(sigma) sum to 2 pi) and
# cluster nodes at the corners. A single corner at x > 0 is only probed for nodes after the corner. No eigen
# solves, no pinned wavefunction values. Helper functions and constants use the prefix `kg_`.

using BilliardGeometry: kress_R_even!, kress_R_odd!, kress_R!, s_mid, kress_graded_nodes_data,
                        multi_kress_graded_nodes_data, _corner_interval

# Absolute tolerance for closed-form comparisons of the FFT-built Kress matrix (a few hundred ulps of Float64).
const kg_ATOL = 1e-12
# Absolute tolerance for Float32 comparisons.
const kg_ATOL_F32 = 1f-5
# Absolute tolerance for sum(wq) = 2 pi; the periodic midpoint rule on the graded Jacobian is only algebraically accurate
# (observed error about 4e-7 for N = 64, q = 4).
const kg_WSUM_ATOL = 1e-5

# First column of the Kress matrix for an even node count N = 2n, from the closed form (no FFT).
kg_even_column(N) = (n = N ÷ 2; [-2 * pi * (sum(cos(m * j * pi / n) / m for m in 1:(n-1)) / n + (-1)^j / (2 * n^2)) for j in 0:(N-1)])
# First column of the Kress matrix for an odd node count N = 2n+1, from the closed form (no FFT).
kg_odd_column(N) = (n = (N - 1) ÷ 2; [-2 * pi * (2 / N) * sum(cos(2 * pi * m * j / N) / m for m in 1:n) for j in 0:(N-1)])

@testset "kress_R_even!" begin
    for N in (4, 8, 16)
        R = zeros(N, N)
        @test kress_R_even!(R) === nothing
        @test R[:, 1] ≈ kg_even_column(N) atol = kg_ATOL
        # circulant and symmetric
        @test R ≈ R' atol = kg_ATOL
        @test all(isapprox.(R[:, 2], circshift(R[:, 1], 1); atol = kg_ATOL))
        # constants are annihilated
        @test all(abs.(vec(sum(R, dims = 2))) .< kg_ATOL)
    end

    # exact for trigonometric polynomials of degree m < n
    N = 16
    R = zeros(N, N)
    kress_R_even!(R)
    t = 2 * pi .* (0:(N-1)) ./ N
    @test all(m -> isapprox(R * cos.(m .* t), -(2 * pi / m) .* cos.(m .* t); atol = kg_ATOL), 1:(N÷2-1))
    @test all(m -> isapprox(R * sin.(m .* t), -(2 * pi / m) .* sin.(m .* t); atol = kg_ATOL), 1:(N÷2-1))

    # element type is preserved
    R32 = zeros(Float32, 16, 16)
    kress_R_even!(R32)
    @test eltype(R32) == Float32
    @test R32 ≈ R atol = kg_ATOL_F32
end

@testset "kress_R_odd!" begin
    for N in (5, 9, 17)
        R = zeros(N, N)
        @test kress_R_odd!(R) === nothing
        @test R[:, 1] ≈ kg_odd_column(N) atol = kg_ATOL
        @test R ≈ R' atol = kg_ATOL
        @test all(isapprox.(R[:, 2], circshift(R[:, 1], 1); atol = kg_ATOL))
        @test all(abs.(vec(sum(R, dims = 2))) .< kg_ATOL)
    end

    # exact for trigonometric polynomials of degree m <= n = (N-1)/2
    N = 17
    R = zeros(N, N)
    kress_R_odd!(R)
    t = 2 * pi .* (0:(N-1)) ./ N
    @test all(m -> isapprox(R * cos.(m .* t), -(2 * pi / m) .* cos.(m .* t); atol = kg_ATOL), 1:((N-1)÷2))
    @test all(m -> isapprox(R * sin.(m .* t), -(2 * pi / m) .* sin.(m .* t); atol = kg_ATOL), 1:((N-1)÷2))

    R32 = zeros(Float32, 17, 17)
    kress_R_odd!(R32)
    @test eltype(R32) == Float32
    @test R32 ≈ R atol = kg_ATOL_F32
end

@testset "kress_R!" begin
    # kress_R! circulant/symmetry structural invariants (Kress logarithmic kernel is even)
    for N in (16, 17, 32, 33)
        R0 = zeros(Float64, N, N)
        kress_R!(R0)
        @test all(isfinite, R0)
        @test isapprox(R0, R0'; atol = 1e-12)
        # circulant: each row is a cyclic shift of the previous
        @test all(isapprox.(R0[:, 2], circshift(R0[:, 1], 1); atol = 1e-10))
    end

    # dispatches on the parity of the matrix size
    Re = zeros(16, 16); kress_R_even!(Re)
    Ro = zeros(17, 17); kress_R_odd!(Ro)
    R16 = zeros(16, 16); kress_R!(R16)
    R17 = zeros(17, 17); kress_R!(R17)
    @test R16 == Re
    @test R17 == Ro
    @test kress_R!(zeros(8, 8)) === nothing
end

@testset "s_mid" begin
    @test s_mid(1, 4) ≈ 2 * pi * 0.5 / 4
    # midpoints of N equal arcs of [0, 2 pi]
    @test [s_mid(k, 4) for k in 1:4] ≈ [pi / 4, 3pi / 4, 5pi / 4, 7pi / 4] atol = kg_ATOL
    @test s_mid(1, 1) ≈ pi atol = kg_ATOL
end

@testset "kress_graded_nodes_data" begin
    # single-corner grading: nodes cluster near sigma=0≡2pi, jacobian vanishes there
    N = 64
    σ, s, jac, jac2, wq = kress_graded_nodes_data(Float64, N; q = 4)
    @test length(σ) == N && length(s) == N
    @test all(isfinite, s) && all(isfinite, jac)
    @test issorted(s) # graded map remains monotone increasing
    @test jac[1] < jac[N÷2] # more clustering (smaller jacobian) near the corner than mid-arc
    @test all(wq .≈ (2 * pi / N) .* jac)

    # the computational nodes are the periodic midpoints and t(sigma) maps [0, 2 pi] onto itself
    @test σ ≈ [s_mid(k, N) for k in 1:N] atol = kg_ATOL
    @test all(0 .< s .< 2 * pi)
    @test sum(wq) ≈ 2 * pi atol = kg_WSUM_ATOL

    # the grading is antisymmetric about sigma = pi: t(2 pi - sigma) = 2 pi - t(sigma)
    @test s ≈ reverse(2 * pi .- s) atol = kg_ATOL

    # the grading order must exceed 1
    @test_throws ErrorException kress_graded_nodes_data(Float64, 16; q = 1)

    # element type is preserved
    σ32, s32, jac32, jac232, wq32 = kress_graded_nodes_data(Float32, 32; q = 3)
    @test eltype(σ32) == Float32 && eltype(s32) == Float32 && eltype(wq32) == Float32
    @test sum(wq32) ≈ 2f0 * Float32(pi) atol = kg_ATOL_F32
end

@testset "multi_kress_graded_nodes_data" begin
    N = 64
    # multi-corner grading with no corners falls back to the uniform grid
    σu, tmap, jacu, jac2u, wqu = multi_kress_graded_nodes_data(Float64, N, Float64[])
    @test tmap == σu
    @test all(jacu .== 1.0)
    @test all(jac2u .== 0.0)
    @test all(wqu .≈ 2 * pi / N)

    # multi-corner grading clusters near each supplied corner location
    corners = [0.0, pi]
    σc, tmapc, jacc, jac2c, wqc = multi_kress_graded_nodes_data(Float64, N, corners; q = 4)
    @test all(isfinite, tmapc)
    icorner = argmin(abs.(σc .- 0.0))
    imid = argmin(abs.(σc .- pi / 2))
    @test jacc[icorner] < jacc[imid]

    # the jacobian is small at the second corner too, and the map is monotone and total-length preserving
    icorner2 = argmin(abs.(σc .- pi))
    @test jacc[icorner2] < jacc[imid]
    @test issorted(tmapc)
    @test sum(wqc) ≈ 2 * pi atol = kg_WSUM_ATOL
    @test all(wqc .≈ (2 * pi / N) .* jacc)

    # the grading order must exceed 1
    @test_throws ErrorException multi_kress_graded_nodes_data(Float64, 16, [0.0]; q = 1)
end

@testset "_corner_interval" begin
    corners = [1.0, 3.0, 5.0]
    # x between two consecutive corners: that pair
    @test _corner_interval(corners, 2.0) == (1.0, 3.0)
    @test _corner_interval(corners, 4.0) == (3.0, 5.0)
    # x on a corner starts the interval at that corner
    @test _corner_interval(corners, 3.0) == (3.0, 5.0)
    # x before the first corner / after the last corner wraps around through 2 pi
    left, right = _corner_interval(corners, 0.5)
    @test left ≈ 5.0 - 2 * pi atol = kg_ATOL
    @test right == 1.0
    left, right = _corner_interval(corners, 6.0)
    @test left == 5.0
    @test right ≈ 1.0 + 2 * pi atol = kg_ATOL

    # a single corner gives the full period starting at the corner
    left, right = _corner_interval([0.0], 2.0)
    @test left == 0.0
    @test right ≈ 2 * pi atol = kg_ATOL
    left, right = _corner_interval([1.0], 2.0)
    @test left == 1.0
    @test right ≈ 1.0 + 2 * pi atol = kg_ATOL
end
