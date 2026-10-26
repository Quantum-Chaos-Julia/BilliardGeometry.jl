# Integration tests for the node orbits that the registered symmetry group of a real billiard induces on its boundary grid:
#   symmetry_index_orbits, full_boundary, component_lengths, curve, fundamental_size, orbit_size, length of SymmetryOrbitMap.
# The N midpoint nodes of the complete boundary are sampled geometrically (uniform_boundary_samples, independent of the orbit
# bookkeeping) and the orbit map is checked against the point maps of the symmetry: every orbit has exactly the expected
# number of members, and the members are, up to ORBIT_ATOL, the base node and its images under the point maps. The cases are
# X, Y, XY and CompositeReflection(X, Y) on the unit circle (D2), the rotation by 2 pi / 3 on the C3 billiard, and the diagonal and
# anti-diagonal reflections on circles whose angular origin is shifted by +pi/4 and -pi/4 (circular PolarBilliards with that reflection
# as their only registered symmetry, since no shipped billiard registers them). A single order-2 reflection also gets the irrep character check: an explicit trailing character -1 gives the
# phase 1 on the fundamental nodes and -1 on their partners. Tolerances come from test/integration/common.jl.
# No eigen solves, no pinned wavefunction values. Helper functions and constants use the prefix `om_`.

using StaticArrays: SVector
using BilliardGeometry: PolarDomain

# Number of boundary nodes: a multiple of 4 and 8 (D2 and diagonal reductions) and the C3 grid, a multiple of 3.
const om_N = 40
const om_N3 = 39

# Point maps of the symmetries as functions of a boundary point.
om_reflect_x(p) = SVector(p[1], -p[2])
om_reflect_y(p) = SVector(-p[1], p[2])
om_reflect_xy(p) = SVector(-p[1], -p[2])
om_reflect_diagonal(p) = SVector(p[2], p[1])
om_reflect_antidiagonal(p) = SVector(-p[2], -p[1])
om_rotate(p, angle) = SVector(cos(angle) * p[1] - sin(angle) * p[2], sin(angle) * p[1] + cos(angle) * p[2])

# Circle of unit radius whose angular origin is shifted by `shift`, with `sym` as the only registered symmetry.
function om_shifted_circle(sym, shift)
    seg = FourierCoeffPolarSegment(Float64[]; shift_angle = shift)
    dom = PolarDomain{Float64}([seg], [curve(seg, 0.0)], 1)
    return PolarBilliard{Float64}(dom, register_symmetries(sym))
end

"""
    om_check_orbit_map(billiard, sym, N, maps; character = nothing)

Check `symmetry_index_orbits(Float64, billiard, N, sym)` against the geometric samples of `full_boundary(billiard)`.
`maps` are the point maps of the non-identity group elements, so every orbit must have `length(maps) + 1` members. When
`character` is given, the orbits built with that trailing irrep character must have phase 1 on the fundamental nodes and
`character` on all other nodes.
"""
function om_check_orbit_map(billiard, sym, N::Int, maps; character = nothing)
    xy = uniform_boundary_samples(full_boundary(billiard), N)
    orbits = symmetry_index_orbits(Float64, billiard, N, sym)
    n_orbit = length(maps) + 1
    @test length(orbits) == N
    @test orbit_size(orbits) == n_orbit
    @test fundamental_size(orbits) == N ÷ n_orbit
    for b in 1:fundamental_size(orbits)
        members = findall(==(b), orbits.orbit_of)
        @test length(members) == n_orbit
        # the orbit is the base node and its images, each matched by exactly one member
        base = xy[orbits.fundamental_indices[b]]
        expected = vcat([base], [f(base) for f in maps])
        @test all(e -> count(m -> norm(e - xy[m]) < ORBIT_ATOL, members) == 1, expected)
    end
    if !isnothing(character)
        orbits_chi = symmetry_index_orbits(Float64, billiard, N, sym, ComplexF64(character))
        @test all(==(one(ComplexF64)), orbits_chi.phase[orbits_chi.fundamental_indices])
        partners = setdiff(1:N, orbits_chi.fundamental_indices)
        @test all(==(ComplexF64(character)), orbits_chi.phase[partners])
    end
    return nothing
end

@testset "orbit maps (XAxisReflection, CircleBilliard)" begin
    om_check_orbit_map(CircleBilliard(1.0), XAxisReflection(), om_N, [om_reflect_x]; character = -1)
end

@testset "orbit maps (YAxisReflection, CircleBilliard)" begin
    om_check_orbit_map(CircleBilliard(1.0), YAxisReflection(), om_N, [om_reflect_y]; character = -1)
end

@testset "orbit maps (XYAxisReflection, CircleBilliard)" begin
    # orbits {p, X(p), Y(p), XY(p)}
    om_check_orbit_map(CircleBilliard(1.0), XYAxisReflection(), om_N, [om_reflect_x, om_reflect_y, om_reflect_xy])
end

@testset "orbit maps (CompositeReflection, CircleBilliard)" begin
    # the closure of X and Y is the same D2 group as XYAxisReflection
    circ = CircleBilliard(1.0)
    comp = CompositeReflection(XAxisReflection(), YAxisReflection())
    om_check_orbit_map(circ, comp, om_N, [om_reflect_x, om_reflect_y, om_reflect_xy])
    orbits_comp = symmetry_index_orbits(Float64, circ, om_N, comp)
    orbits_xy = symmetry_index_orbits(Float64, circ, om_N, XYAxisReflection())
    @test fundamental_size(orbits_comp) == fundamental_size(orbits_xy)
    @test orbit_size(orbits_comp) == orbit_size(orbits_xy)
end

@testset "orbit maps (NFoldRotation, C3Billiard)" begin
    # orbits {p, R(p), R^2(p)} with R the rotation by 2 pi / 3
    maps = [p -> om_rotate(p, 2pi / 3), p -> om_rotate(p, 4pi / 3)]
    om_check_orbit_map(C3Billiard(0.2), NFoldRotation(3, 1), om_N3, maps)
end

@testset "orbit maps (DiagonalReflection, shifted circle)" begin
    om_check_orbit_map(om_shifted_circle(DiagonalReflection(), pi / 4), DiagonalReflection(), om_N, [om_reflect_diagonal]; character = -1)
end

@testset "orbit maps (AntiDiagonalReflection, shifted circle)" begin
    om_check_orbit_map(om_shifted_circle(AntiDiagonalReflection(), -pi / 4), AntiDiagonalReflection(), om_N, [om_reflect_antidiagonal]; character = -1)
end
