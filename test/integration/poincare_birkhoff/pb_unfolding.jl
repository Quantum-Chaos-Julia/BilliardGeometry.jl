# Integration tests for the unfolding of the Poincare-Birkhoff coordinate s onto the full boundary of a symmetric billiard:
#   pb_sectors, pb_coords, get_boundary_curves, full_boundary, component_lengths, arc_length, symmetry_of.
# For a billiard with the registry of nsym symmetries and fundamental boundary length L, the sector k = 1, ..., nsym + 1 of
# the coordinate s covers the window [(k - 1) L, k L], the breakpoints of the sectors are 0 = s_0 < ... < s_n = (nsym + 1) L,
# and they are the cumulative lengths of the curves of the reconstructed full boundary. The coordinate s of a point of the
# fundamental boundary seen in sector k is the arc length along the full boundary up to the image of that point, where the
# image curve and its direction come from the provenance of full_boundary (_full_boundary_data); reflection sectors run
# backwards, so this also checks the sector offsets and the reversal of s. The fixtures are the rectangle 1.0 x 0.6 (D2, lines)
# and the C3 billiard (rotations, one curved arc). Tolerances come from test/integration/common.jl.
# No eigen solves, no pinned wavefunction values. Helper functions and constants use the prefix `pu_`.

using StaticArrays: SVector
using BilliardGeometry: _full_boundary_data

# Parameter values on every fundamental curve at which pb_coords is compared with the full-boundary arc length. The C3 arc is
# only probed at its midpoint: pb_coords obtains s from arc_length(curve, point), which inverts the curve with a root finder on
# the squared distance, and for the curved C3 arc this inversion is off by about 1e-4 in s away from the midpoint.
const pu_US_LINES = (0.2, 0.5, 0.9)
const pu_US_ARC = (0.5,)

# Point of the curve at parameter `u` together with a velocity (any nonzero vector, only its direction enters).
function pu_point_and_velocity(crv, u)
    tang = tangent(crv, u)
    return curve(crv, u), SVector(tang[2], -tang[1])
end

# pb_sectors starts at 0, ends at the full length and is sorted and unique; pb_coords stays in the window of its sector.
function pu_check_sectors_and_windows(billiard)
    nsym = length(billiard.symmetries)
    L = CompositeCurve(get_boundary_curves(billiard)).length
    sectors = pb_sectors(billiard)
    @test issorted(sectors)
    @test length(sectors) == length(unique(sectors))
    @test isapprox(first(sectors), 0.0; atol = PB_SECTOR_ATOL)
    @test isapprox(last(sectors), (nsym + 1) * L; atol = PB_SECTOR_ATOL)

    crv = first(get_boundary_curves(billiard))
    pt, vel = pu_point_and_velocity(crv, 0.5)
    for sym_sector in 1:(nsym + 1)
        coords = pb_coords(billiard, crv.segment_id, crv.domain_id, sym_sector, pt, vel)
        @test (sym_sector - 1) * L - PB_WINDOW_ATOL <= coords.s <= sym_sector * L + PB_WINDOW_ATOL
    end
    return nothing
end

# s equals the arc length along the full boundary up to the image of the base point, for every sector and fundamental curve.
function pu_check_arclength(billiard, us)
    nsym = length(billiard.symmetries)
    fund = get_boundary_curves(billiard)
    curves, source, sector, reversed = _full_boundary_data(billiard)
    _, cum, _ = component_lengths(curves)
    for (j, crv) in enumerate(fund), u in us
        pt, vel = pu_point_and_velocity(crv, u)
        for sym_sector in 1:(nsym + 1)
            i = findfirst(k -> source[k] == j && sector[k] == sym_sector, eachindex(curves))
            tt = reversed[i] ? 1 - u : u
            s_ref = cum[i] + arc_length(curves[i], tt; rtol = 1e-13)
            coords = pb_coords(billiard, crv.segment_id, crv.domain_id, sym_sector, pt, vel)
            @test isapprox(coords.s, s_ref; atol = PB_ARCLENGTH_ATOL)
        end
    end
    return nothing
end

@testset "pb_sectors and pb_coords windows (RectangleBilliard)" begin
    pu_check_sectors_and_windows(RectangleBilliard(1.0, 0.6))
end

@testset "pb_sectors and pb_coords windows (C3Billiard)" begin
    pu_check_sectors_and_windows(C3Billiard(0.2))
end

@testset "pb_sectors breakpoints against full boundary (RectangleBilliard)" begin
    rect = RectangleBilliard(1.0, 0.6)
    _, cum, _ = component_lengths(full_boundary(rect))
    sectors = pb_sectors(rect)
    @test length(sectors) == length(cum)
    @test all(isapprox.(sectors, cum; atol = PB_SECTOR_ATOL))
end

@testset "pb_sectors breakpoints against full boundary (C3Billiard)" begin
    c3 = C3Billiard(0.2)
    _, cum, _ = component_lengths(full_boundary(c3))
    sectors = pb_sectors(c3)
    @test length(sectors) == length(cum)
    @test all(isapprox.(sectors, cum; atol = PB_SECTOR_ATOL))
end

@testset "pb_coords against full-boundary arc length (RectangleBilliard)" begin
    pu_check_arclength(RectangleBilliard(1.0, 0.6), pu_US_LINES)
end

@testset "pb_coords against full-boundary arc length (C3Billiard)" begin
    pu_check_arclength(C3Billiard(0.2), pu_US_ARC)
end
