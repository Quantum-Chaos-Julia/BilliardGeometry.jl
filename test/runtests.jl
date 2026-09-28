using BilliardGeometry
using Test
using StaticArrays
using LinearAlgebra
using CoordinateTransformations

@testset "linesegment.jl" begin
    pt0, pt1 = [0.0,0.0],  [1.0,1.0]
    crv = LineSegment(pt0, pt1) 
    # curve functions test
    @test curve(crv, 0.0) == pt0
    @test curve(crv, 1.0) == pt1
    @test arc_length(crv, crv.pt1) == crv.length
    # gradient functions test
    @test domain_gradient_vector(crv, curve(crv,0.5)) == SVector{2,Float64}([1.0,-1.0])
    # domain functions test
    testpt1 = SVector{2,Float64}([0.1,0.2])
    testpt2 = SVector{2,Float64}([0.1,0.0])
    @test domain_fun(crv, testpt1) < 0.0
    @test domain_fun(crv, testpt2) > 0.0
    @test is_inside(crv, [testpt1,testpt2]) == [true, false]

end

@testset "circlesegment.jl" begin
    R = 1.0
    arc_angle = 0.5*pi
    shift_angle = 0.0
    center = [0.0, 0.0]
    crv = CircleSegment(R, arc_angle, shift_angle, center)
    # curve functions test
    @test all(isapprox.(curve(crv, 0.0), SVector{2,Float64}([1.0,0.0]); atol=1e-8))
    @test all(isapprox.(curve(crv, 1.0), SVector{2,Float64}([0.0,1.0]); atol=1e-8))
    @test arc_length(crv, curve(crv,1.0)) ≈ crv.length
    # gradient functions test
    @test domain_gradient_vector(crv, curve(crv,0.5)) ≈ SVector{2,Float64}([sqrt(2)/2,sqrt(2)/2])
    # domain functions test
    testpt1 = SVector{2,Float64}([0.1,0.2])
    testpt2 = SVector{2,Float64}([ 1.1,1.0])
    @test domain_fun(crv, testpt1) < 0.0
    @test domain_fun(crv, testpt2) > 0.0
    @test is_inside(crv, [testpt1,testpt2]) == [true, false]
end

@testset "limaconsegment.jl" begin
    
end

@testset "polarsegment.jl" begin
    # Step 10.9 item C: `PolarSegment`->`FourierCoeffPolarSegment` rename, with
    # a `promote_type` fix so mismatched literal types (Int `R`, Float32
    # `coef`, Float64 `center`) still produce a single consistent `T`.
    seg = FourierCoeffPolarSegment(Float32[0.0,0.0,0.0,0.3]; R=1, center=[0.0,0.0])
    @test typeof(seg).parameters[1] == Float64

    # New function-based `PolarSegment{T,BC,F}`: derivatives via ForwardDiff,
    # checked against central finite differences of `curve(t)`.
    rfunc(phi) = 1.0 + 0.3*cos(2*phi)
    pseg = PolarSegment(rfunc)
    h = 1e-6
    for t0 in (0.1, 0.4, 0.7)
        xp = curve(pseg, t0+h)
        xm = curve(pseg, t0-h)
        x0 = curve(pseg, t0)
        fd1 = (xp .- xm) ./ (2h)
        fd2 = (xp .- 2 .* x0 .+ xm) ./ (h^2)
        @test all(isapprox.(fd1, tangent(pseg, t0); atol=1e-5))
        @test all(isapprox.(fd2, tangent_2(pseg, t0); atol=1e-3))
    end
end

@testset "curvederivatives.jl" begin
    h = 1e-6
    # CircleSegment: compare tangent/tangent_2 against central finite differences of curve(t)
    circ = CircleSegment(2.0, 0.5*pi, 0.3, [0.1,-0.2])
    for t0 in (0.13, 0.5, 0.87)
        xp = curve(circ, t0+h)
        xm = curve(circ, t0-h)
        x0 = curve(circ, t0)
        fd1 = (xp .- xm) ./ (2h)
        fd2 = (xp .- 2 .* x0 .+ xm) ./ (h^2)
        @test all(isapprox.(fd1, tangent(circ, t0); atol=1e-5))
        @test all(isapprox.(fd2, tangent_2(circ, t0); atol=1e-3))
    end
    # LineSegment: tangent is the constant chord vector, tangent_2 is zero
    line = LineSegment([0.0,0.0], [2.0,-1.0])
    @test tangent(line, 0.0) == line.pt1 .- line.pt0
    @test tangent(line, 0.7) == line.pt1 .- line.pt0
    @test tangent_2(line, 0.4) == zero(line.pt0)
    # array-argument methods broadcast the scalar method
    ts = [0.1,0.4,0.9]
    @test tangent(circ, ts) == [tangent(circ,t) for t in ts]
    @test tangent_2(circ, ts) == [tangent_2(circ,t) for t in ts]

    # Step 10.9 item D: `tangent_vec`/`normal_vec`/`curvature`. A circle of
    # radius R has constant curvature 1/R, and its unit outward normal at
    # parameter t is the radial direction (cos(phi),sin(phi)).
    R = 2.0
    circ0 = CircleSegment(R, 2*pi, 0.0, [0.0,0.0])
    ts0 = collect(range(0,1,length=5))
    @test all(isapprox.(curvature(circ0, ts0), 1/R; atol=1e-12))
    @test isapprox(curvature(circ0, ts0[2]), 1/R; atol=1e-12)
    tv = tangent_vec(circ0, ts0)
    @test all(isapprox.(norm.(tv), 1.0; atol=1e-12))
    nv = normal_vec(circ0, ts0)
    phis = 2*pi .* ts0
    @test all(isapprox.(nv, [SVector(cos(phi),sin(phi)) for phi in phis]; atol=1e-10))
end

@testset "boundarycomponents.jl" begin
    # two perpendicular unit line segments forming a right-angle corner
    line1 = LineSegment([0.0,0.0], [1.0,0.0])
    line2 = LineSegment([1.0,0.0], [1.0,1.0])
    comp = [line1, line2]
    lens, cum, Ltot = component_lengths(comp)
    @test lens == [1.0,1.0]
    @test cum == [0.0,1.0,2.0]
    @test Ltot == 2.0
    @test BilliardGeometry._is_true_corner(line1, line2, Float64)
    corners = BilliardGeometry._component_corner_locations(Float64, comp)
    @test length(corners) == 2
    @test 0.0 in corners
    @test isapprox(maximum(corners), pi; atol=1e-10)

    # two colinear segments have a smooth join (no corner)
    line3 = LineSegment([1.0,0.0], [2.0,0.0])
    @test !BilliardGeometry._is_true_corner(line1, line3, Float64)
    @test isempty(BilliardGeometry._component_corner_locations(Float64, [line1,line3]))

    # `_boundary_components` was renamed/replaced by `_group_curves_by_domain_id`
    # (groups a flat physical-boundary curve list into connected components by
    # curve `domain_id`, preserving first-seen order): curves with distinct
    # `domain_id`s become distinct components, curves sharing one `domain_id`
    # stay in a single component.
    line1b = LineSegment([0.0,0.0], [1.0,0.0]; domain_id=1)
    line2b = LineSegment([1.0,0.0], [1.0,1.0]; domain_id=2)
    comps = BilliardGeometry._group_curves_by_domain_id(AbsCurve[line1b,line2b])
    @test length(comps) == 2
    @test comps[1] == [line1b]
    @test comps[2] == [line2b]
    comps_same = BilliardGeometry._group_curves_by_domain_id(AbsCurve[line1,line2])
    @test length(comps_same) == 1
    @test comps_same[1] == [line1,line2]
end

@testset "kressgrading.jl" begin
    # kress_R! circulant/symmetry structural invariants (Kress logarithmic kernel is even)
    for N in (16,17,32,33)
        R0 = zeros(Float64, N, N)
        kress_R!(R0)
        @test all(isfinite, R0)
        @test isapprox(R0, R0'; atol=1e-12)
        # circulant: each row is a cyclic shift of the previous
        @test all(isapprox.(R0[:,2], circshift(R0[:,1],1); atol=1e-10))
    end

    # single-corner grading: nodes cluster near sigma=0≡2pi, jacobian vanishes there
    N = 64
    σ,s,jac,jac2,wq = kress_graded_nodes_data(Float64, N; q=4)
    @test length(σ) == N && length(s) == N
    @test all(isfinite, s) && all(isfinite, jac)
    @test issorted(s) # graded map remains monotone increasing
    @test jac[1] < jac[N÷2] # more clustering (smaller jacobian) near the corner than mid-arc
    @test all(wq .≈ (2*pi/N).*jac)

    # multi-corner grading with no corners falls back to the uniform grid
    σu,tmap,jacu,jac2u,wqu = multi_kress_graded_nodes_data(Float64, N, Float64[])
    @test tmap == σu
    @test all(jacu .== 1.0)
    @test all(jac2u .== 0.0)

    # multi-corner grading clusters near each supplied corner location
    corners = [0.0, pi]
    σc,tmapc,jacc,jac2c,wqc = multi_kress_graded_nodes_data(Float64, N, corners; q=4)
    @test all(isfinite, tmapc)
    icorner = argmin(abs.(σc .- 0.0))
    imid = argmin(abs.(σc .- pi/2))
    @test jacc[icorner] < jacc[imid]

    @test s_mid(1,4) ≈ 2*pi*0.5/4
end

@testset "symmetryorbits.jl" begin
    # `NFoldRotation(N,m)`'s convenience constructor bug (it used to build
    # `LinearMap(RotZ(...))`, a 3x3 rotation, which could not convert to the
    # struct's declared `LinearMap{SMatrix{2,2,Float64,4}}` field type) was
    # fixed as part of the Step 4 symmetry-infrastructure migration (now uses
    # `rotation_matrix_z`), so the convenience constructor can be used
    # directly here.
    _make_nfold(n; m=1) = NFoldRotation(n, m)

    @test symmetry_node_multiple(XAxisReflection()) == 4
    @test symmetry_node_multiple(YAxisReflection()) == 4
    @test symmetry_node_multiple(XYAxisReflection()) == 4
    @test symmetry_node_multiple(_make_nfold(6)) == 6

    # Step 15 migration: `symmetry_index_orbits` no longer takes a raw point
    # vector; it derives the exact midpoint-node permutation from a real
    # billiard's registered `SymmetryRegistry` (`billiard.symmetries`), never
    # from floating-point coordinates. Every sub-test below therefore builds
    # a real billiard with the relevant symmetry actually registered, then
    # samples true boundary points from its `full_boundary` reconstruction to
    # independently check the returned orbit permutation against the actual
    # geometric symmetry (mirrors the "fullboundary.jl" testset's index-
    # permutation consistency check below).
    function _sample_full_boundary(billiard, N)
        fb = full_boundary(billiard)
        _, _, Lt = component_lengths(fb)
        function _point_at_sigma(sigma)
            target = Lt*sigma/(2*pi)
            offset = 0.0
            for j in eachindex(fb)
                Lj = fb[j].length
                if target < offset+Lj || j==lastindex(fb)
                    u = clamp((target-offset)/Lj, 0.0, 1.0)
                    return curve(fb[j], u)
                end
                offset += Lj
            end
        end
        return [_point_at_sigma(2*pi*(k-0.5)/N) for k in 1:N]
    end

    # CircleBilliard: real D2-symmetric billiard (sym_id 1=Y, 2=XY, 3=X).
    circ = CircleBilliard(1.0)
    N = 40
    xy = _sample_full_boundary(circ, N)

    for (sym, reflect) in ((XAxisReflection(), pt->SVector(pt[1],-pt[2])),
                           (YAxisReflection(), pt->SVector(-pt[1],pt[2])))
        orbits = symmetry_index_orbits(Float64, circ, N, sym)
        @test fundamental_size(orbits) == N÷2
        @test length(orbits) == N
        for b in 1:fundamental_size(orbits)
            members = findall(==(b), orbits.orbit_of)
            @test length(members) == 2
            q1,q2 = members
            @test isapprox(reflect(xy[q1]), xy[q2]; atol=1e-8) || isapprox(reflect(xy[q2]), xy[q1]; atol=1e-8)
        end
        # Step 15 removed `symmetry_irrep_character`: the irrep character is
        # now an explicit trailing argument to `symmetry_index_orbits`
        # (defaulting to the trivial representation). Passing `-1` explicitly
        # reproduces the old default-antisymmetric-partner behavior.
        orbits_χ = symmetry_index_orbits(Float64, circ, N, sym, ComplexF64(-1))
        @test all(==(one(ComplexF64)), orbits_χ.phase[orbits_χ.fundamental_indices])
        partners = setdiff(1:N, orbits_χ.fundamental_indices)
        @test all(==(ComplexF64(-1)), orbits_χ.phase[partners])
    end

    # XYAxisReflection: built from the composition of X and Y, giving
    # 4-element orbits {p, X(p), Y(p), XY(p)}.
    orbits_xy = symmetry_index_orbits(Float64, circ, N, XYAxisReflection())
    @test fundamental_size(orbits_xy) == N÷4
    for b in 1:fundamental_size(orbits_xy)
        members = findall(==(b), orbits_xy.orbit_of)
        @test length(members) == 4
        p = xy[orbits_xy.fundamental_indices[b]]
        expected = Set([p, SVector(p[1],-p[2]), SVector(-p[1],p[2]), SVector(-p[1],-p[2])])
        actual = [xy[m] for m in members]
        @test all(e -> any(a -> isapprox(e,a;atol=1e-8), actual), expected)
    end

    # fund_to_full/fund_to_scale/symmetry_orbit: consistent with orbit_of/phase
    for b in 1:fundamental_size(orbits_xy)
        qs, χs = symmetry_orbit(orbits_xy, b)
        @test orbit_size(orbits_xy) == length(qs) == 4
        for (q,χ) in zip(qs,χs)
            @test orbits_xy.orbit_of[q] == b
            @test orbits_xy.phase[q] == χ
        end
    end
    @test full_size(orbits_xy) == N

    # NFoldRotation: real C3-symmetric billiard (sym_id 1,2 = rotate by
    # +2π/3, +4π/3).
    c3 = C3Billiard(0.2)
    n = 3
    N3 = 39
    xy3 = _sample_full_boundary(c3, N3)
    orbits_rot = symmetry_index_orbits(Float64, c3, N3, _make_nfold(n))
    @test fundamental_size(orbits_rot) == N3÷n
    rotate(pt) = BilliardGeometry.rotation_matrix_z(2*pi/n)*pt
    for b in 1:fundamental_size(orbits_rot)
        members = findall(==(b), orbits_rot.orbit_of)
        @test length(members) == n
        p = xy3[orbits_rot.fundamental_indices[b]]
        actual = [xy3[m] for m in members]
        pk = p
        for _ in 1:n
            @test any(a -> isapprox(pk,a;atol=1e-8), actual)
            pk = rotate(pk)
        end
    end

    # DiagonalReflection: real fixture (`SquareWithinSquareBilliard` registers
    # exactly this reflection, sym_id 1).
    sq = SquareWithinSquareBilliard(1.0)
    @test symmetry_node_multiple(DiagonalReflection()) == 8
    Nd = 40
    xy_diag = _sample_full_boundary(sq, Nd)
    orbits_diag = symmetry_index_orbits(Float64, sq, Nd, DiagonalReflection())
    @test fundamental_size(orbits_diag) == Nd÷2
    @test length(orbits_diag) == Nd
    for b in 1:fundamental_size(orbits_diag)
        members = findall(==(b), orbits_diag.orbit_of)
        @test length(members) == 2
        q1,q2 = members
        reflect_diag(pt) = SVector(pt[2],pt[1])
        @test isapprox(reflect_diag(xy_diag[q1]), xy_diag[q2]; atol=1e-8) || isapprox(reflect_diag(xy_diag[q2]), xy_diag[q1]; atol=1e-8)
    end

    # AntiDiagonalReflection: no shipped billiard registers this reflection
    # alone, so a real (perfectly circular) `PolarBilliard` fixture is built
    # here with it as the sole registered generator (a single non-identity
    # reflection is always a closed order-2 group, so this is a valid
    # `SymmetryRegistry` regardless of the underlying shape). The circle's
    # angular origin is shifted by `-π/4` so the periodic midpoint grid's
    # sector seams align with the anti-diagonal mirror line.
    @test symmetry_node_multiple(AntiDiagonalReflection()) == 8
    ad_seg = FourierCoeffPolarSegment(Float64[]; shift_angle=-pi/4)
    ad_dom = BilliardGeometry.PolarDomain{Float64}([ad_seg], [curve(ad_seg,0.0)], 1)
    circ_ad = PolarBilliard{Float64}(ad_dom, register_symmetries(AntiDiagonalReflection()))
    xy_ad = _sample_full_boundary(circ_ad, N)
    orbits_ad = symmetry_index_orbits(Float64, circ_ad, N, AntiDiagonalReflection())
    @test fundamental_size(orbits_ad) == N÷2
    @test length(orbits_ad) == N
    for b in 1:fundamental_size(orbits_ad)
        members = findall(==(b), orbits_ad.orbit_of)
        @test length(members) == 2
        q1,q2 = members
        reflect_antidiag(pt) = SVector(-pt[2],-pt[1])
        @test isapprox(reflect_antidiag(xy_ad[q1]), xy_ad[q2]; atol=1e-8) || isapprox(reflect_antidiag(xy_ad[q2]), xy_ad[q1]; atol=1e-8)
    end

    # CompositeReflection(XAxisReflection(),YAxisReflection()) generates the
    # same D2 group as XYAxisReflection() by closure, on the same
    # `CircleBilliard` fixture used above (its registry already contains X,Y,XY).
    comp = CompositeReflection(XAxisReflection(), YAxisReflection())
    @test symmetry_node_multiple(comp) == 4
    orbits_comp = symmetry_index_orbits(Float64, circ, N, comp)
    @test fundamental_size(orbits_comp) == N÷4
    for b in 1:fundamental_size(orbits_comp)
        members = findall(==(b), orbits_comp.orbit_of)
        @test length(members) == 4
        p = xy[orbits_comp.fundamental_indices[b]]
        expected = Set([p, SVector(p[1],-p[2]), SVector(-p[1],p[2]), SVector(-p[1],-p[2])])
        actual = [xy[m] for m in members]
        @test all(e -> any(a -> isapprox(e,a;atol=1e-8), actual), expected)
    end
    # The generated closure and the native `XYAxisReflection` registration
    # reach the same group order/reduction on the same billiard: two paths
    # to the same D2 group element agree on a real fixture.
    @test fundamental_size(orbits_comp) == fundamental_size(orbits_xy)
    @test orbit_size(orbits_comp) == orbit_size(orbits_xy)

    # `get_symmetries`: minimal generating subset, not the full registry.
    # A D2 billiard registers 3 elements (Y,XY,X) but only 2 are needed to
    # generate the group by composition.
    rect = RectangleBilliard(1.0, 0.6)
    gens_d2 = get_symmetries(rect)
    @test length(gens_d2) == 2
    @test all(g -> g isa YAxisReflection || g isa XAxisReflection, gens_d2)
    @test !any(g -> g isa XYAxisReflection, gens_d2)
    # A cyclic Cn billiard registers n-1 rotation images but is generated by
    # the single m=1 element alone.
    gens_c3 = get_symmetries(c3)
    @test length(gens_c3) == 1
    @test only(gens_c3) isa NFoldRotation
    @test only(gens_c3).m == 1
    # No registered symmetry: empty generating set.
    tri_empty = TriangleBilliard(1.0, 1.0)
    @test get_symmetries(tri_empty) == ()

    # `symmetry_of` error path: an unknown sym_id must raise ArgumentError,
    # not silently return `nothing` or the wrong generator.
    @test_throws ArgumentError symmetry_of(rect.symmetries, 999)

    # `SymmetryOrbitMap` internal consistency (partition/counting sanity
    # checks computable from the orbit map alone, without any externally
    # known reference value): every fundamental orbit has exactly
    # `orbit_size` members, the orbits partition `1:full_size` exactly once,
    # and `fundamental_size * orbit_size == full_size`.
    @test fundamental_size(orbits_xy) * orbit_size(orbits_xy) == full_size(orbits_xy)
    @test sum(b -> count(==(b), orbits_xy.orbit_of), 1:fundamental_size(orbits_xy)) == full_size(orbits_xy)
    all_members = Int[]
    for b in 1:fundamental_size(orbits_xy)
        qs, _ = symmetry_orbit(orbits_xy, b)
        append!(all_members, qs)
    end
    @test sort(all_members) == collect(1:full_size(orbits_xy))
    @test length(all_members) == length(unique(all_members))
end

@testset "fullboundary.jl" begin
    # Trivial symmetry: full_boundary reproduces get_boundary_curves exactly.
    tri = TriangleBilliard(1.0, 1.0)
    @test full_boundary(tri) == get_boundary_curves(tri)

    # D2-symmetric stadium: full_boundary reconstructs the complete closed
    # physical boundary from the quarter fundamental domain.
    hw = 0.5
    stad = StadiumBilliard(hw)
    fb = full_boundary(stad)
    @test length(fb) == 4*length(get_boundary_curves(stad))

    # Total arc length equals the full stadium perimeter (two straight edges
    # of length 2*hw each, plus the full circle circumference 2*pi).
    Ltot = sum(c.length for c in fb)
    @test isapprox(Ltot, 4*hw + 2*pi; atol=1e-10)

    # Closed, continuous CCW loop: each curve's endpoint matches the next
    # curve's start point.
    for i in eachindex(fb)
        p_end = curve(fb[i], 1.0)
        p_start = curve(fb[mod1(i+1,length(fb))], 0.0)
        @test isapprox(p_end, p_start; atol=1e-10)
    end

    # Exact index-permutation consistency: sampling full_boundary at N
    # midpoint nodes and applying apply_symmetry must land on the node
    # predicted by the canonical periodic reflection index maps.
    N = 40
    lens, cum, Lt = component_lengths(fb)
    function _point_at_sigma(fb, sigma, Lt)
        target = Lt*sigma/(2*pi)
        offset = 0.0
        for j in eachindex(fb)
            Lj = fb[j].length
            if target < offset+Lj || j==lastindex(fb)
                u = clamp((target-offset)/Lj, 0.0, 1.0)
                return curve(fb[j], u)
            end
            offset += Lj
        end
    end
    xy = [_point_at_sigma(fb, 2*pi*(k-0.5)/N, Lt) for k in 1:N]
    _idx_reflect_x(q,N) = mod1(N-q+1,N)
    _idx_reflect_y(q,N) = mod1(N÷2-q+1,N)
    for q in 1:N
        @test isapprox(apply_symmetry(XAxisReflection(), xy[q]), xy[_idx_reflect_x(q,N)]; atol=1e-8)
        @test isapprox(apply_symmetry(YAxisReflection(), xy[q]), xy[_idx_reflect_y(q,N)]; atol=1e-8)
    end

    # `full_boundary` on a rotation (Cn) billiard: generalizes the D2 stadium
    # check above to `NFoldRotation`. `C3Billiard`'s fundamental domain has a
    # single physical curve (the arc; the two wedge-cut walls are
    # `SymmetryWall`-tagged and excluded), so the reconstructed boundary has
    # exactly 3 copies of it, forming a closed CCW loop.
    c3 = C3Billiard(0.2)
    fb_c3 = full_boundary(c3)
    @test length(fb_c3) == 3*length(get_boundary_curves(c3))
    for i in eachindex(fb_c3)
        p_end = curve(fb_c3[i], 1.0)
        p_start = curve(fb_c3[mod1(i+1,length(fb_c3))], 0.0)
        @test isapprox(p_end, p_start; atol=1e-8)
    end

    # Step 10.9 item B: `full_boundary` support for `FourierCoeffPolarSegment`
    # under a D2-symmetric polar billiard. The core correctness identity is
    # `curve(_apply_symmetry_to_curve(sym,c), t) == apply_symmetry(sym, curve(c,t))`
    # (checked directly, white-box, since `PolarBilliard`'s default full-circle
    # fundamental curve is not itself a proper fundamental sector, so a
    # closed-loop reconstruction check like the Stadium one above does not
    # apply here -- that requires an actual quarter-sector polar billiard,
    # deferred to Step 11).
    polar_bil = PolarBilliard([0.0,0.0,0.0,0.3])
    orig_polar_curve = polar_bil.fundamental_domain.boundary[1]
    for sym in (XAxisReflection(), YAxisReflection(), XYAxisReflection())
        g = BilliardGeometry._apply_symmetry_to_curve(sym, orig_polar_curve)
        @test g isa FourierCoeffPolarSegment
        for t in (0.0, 0.2, 0.5, 0.8, 1.0)
            @test isapprox(curve(g, t), apply_symmetry(sym, curve(orig_polar_curve, t)); atol=1e-8)
        end
    end

    # `full_boundary` itself runs without error on a symmetric polar billiard
    # and returns the expected number/type of curves (fundamental curve plus
    # one image per non-identity symmetry).
    D2sym = [YAxisReflection(), XYAxisReflection(), XAxisReflection()]
    polar_bil_d2 = PolarBilliard{Float64}(polar_bil.fundamental_domain, register_symmetries(D2sym...))
    fb_polar = full_boundary(polar_bil_d2)
    @test length(fb_polar) == 4
    @test all(c isa FourierCoeffPolarSegment for c in fb_polar)
end

@testset "boundarytypes.jl - SymmetryWall" begin
    # Round-trip: fields are read back unchanged, and a `SymmetryWall`-tagged
    # curve is excluded from `get_boundary_curves`'s physical-boundary filter
    # (`typeof(crv.bc) <: SpecularReflection`), since it is a fundamental-
    # domain cut, not part of the physical boundary. `RectangleBilliard`
    # tags its two symmetry walls with `SymmetryWall(1,2)`/`SymmetryWall(3,2)`
    # (see rectangle.jl), giving a real fixture rather than a synthetic curve.
    wall = SymmetryWall(2, 3)
    @test wall.sym_id == 2
    @test wall.sector_id == 3

    rect = RectangleBilliard(1.0, 0.6)
    bc = get_boundary_curves(rect)
    @test length(bc) == 2 # only the two SpecularReflection curves; both SymmetryWall curves excluded
    @test all(c -> typeof(c.bc) <: SpecularReflection, bc)

    walls = filter(c -> c.bc isa SymmetryWall, rect.fundamental_domain.boundary)
    @test length(walls) == 2
    @test Set((w.bc.sym_id, w.bc.sector_id) for w in walls) == Set([(1,2),(3,2)])
end

@testset "poincarebirkhoff.jl" begin
    # Smoke test on a real D2 billiard (`RectangleBilliard`) and a real Cn
    # billiard (`C3Billiard`): `pb_sectors` must return a sorted, deduplicated
    # set of breakpoints starting at 0 and ending at the full (symmetry-
    # unfolded) perimeter; `pb_coords` for a point known to lie in the
    # `sym_sector`-th copy must map back into that copy's own global-
    # coordinate window `[(sym_sector-1)*L, sym_sector*L]` (an exact,
    # computable structural invariant of `apply_symmetry_pb`'s reconstruction
    # formula -- both the orientation-reversing and orientation-preserving
    # branches produce this same range -- not a fabricated numeric reference).
    for billiard in (RectangleBilliard(1.0, 0.6), C3Billiard(0.2))
        nsym = length(billiard.symmetries)
        L = CompositeCurve(get_boundary_curves(billiard)).length

        sectors = pb_sectors(billiard)
        @test issorted(sectors)
        @test length(sectors) == length(unique(sectors))
        @test isapprox(first(sectors), 0.0; atol=1e-8)
        @test isapprox(last(sectors), (nsym+1)*L; atol=1e-8)

        crv = first(get_boundary_curves(billiard))
        pt = curve(crv, 0.5)
        tang = tangent(crv, 0.5)
        vel = SVector(tang[2], -tang[1]) # any nonzero velocity; only its direction matters
        for sym_sector in 1:(nsym+1)
            coords = pb_coords(billiard, crv.segment_id, crv.domain_id, sym_sector, pt, vel)
            lo, hi = (sym_sector-1)*L, sym_sector*L
            @test lo - 1e-8 <= coords.s <= hi + 1e-8
        end
    end
end

