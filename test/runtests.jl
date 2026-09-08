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

    # _boundary_components normalizes bare curve vectors into single-segment components
    comps = BilliardGeometry._boundary_components([line1,line2])
    @test length(comps) == 2
    @test comps[1] == [line1]
    @test comps[2] == [line2]
    # already-nested input is passed through unchanged
    @test BilliardGeometry._boundary_components([comp]) == [comp]
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

    N = 40
    xy = [SVector(cos(2*pi*(k-0.5)/N), sin(2*pi*(k-0.5)/N)) for k in 1:N]

    for (sym, reflect) in ((XAxisReflection(), pt->SVector(pt[1],-pt[2])),
                           (YAxisReflection(), pt->SVector(-pt[1],pt[2])))
        orbits = symmetry_index_orbits(Float64, xy, sym)
        @test fundamental_size(orbits) == N÷2
        @test length(orbits) == N
        for b in 1:fundamental_size(orbits)
            members = findall(==(b), orbits.orbit_of)
            @test length(members) == 2
            q1,q2 = members
            @test isapprox(reflect(xy[q1]), xy[q2]; atol=1e-10)
        end
        @test all(==(one(ComplexF64)), orbits.phase)
    end

    orbits_xy = symmetry_index_orbits(Float64, xy, XYAxisReflection())
    @test fundamental_size(orbits_xy) == N÷4
    for b in 1:fundamental_size(orbits_xy)
        members = findall(==(b), orbits_xy.orbit_of)
        @test length(members) == 4
    end

    n = 5
    orbits_rot = symmetry_index_orbits(Float64, xy[1:N-mod(N,n)], _make_nfold(n))
    Nr = N-mod(N,n)
    @test fundamental_size(orbits_rot) == Nr÷n
    for b in 1:fundamental_size(orbits_rot)
        @test length(findall(==(b), orbits_rot.orbit_of)) == n
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

    # DiagonalReflection / AntiDiagonalReflection: two-element orbits
    @test symmetry_node_multiple(DiagonalReflection()) == 8
    @test symmetry_node_multiple(AntiDiagonalReflection()) == 8
    for (sym, reflect) in ((DiagonalReflection(), pt->SVector(pt[2],pt[1])),
                           (AntiDiagonalReflection(), pt->SVector(-pt[2],-pt[1])))
        orbits = symmetry_index_orbits(Float64, xy, sym)
        @test fundamental_size(orbits) == N÷2
        @test length(orbits) == N
        for b in 1:fundamental_size(orbits)
            members = findall(==(b), orbits.orbit_of)
            @test length(members) == 2
            q1,q2 = members
            @test isapprox(reflect(xy[q1]), xy[q2]; atol=1e-10)
        end
    end

    # CompositeReflection(XAxisReflection(),YAxisReflection()) generates the
    # same D2 group as XYAxisReflection() by closure.
    comp = CompositeReflection(XAxisReflection(), YAxisReflection())
    @test symmetry_node_multiple(comp) == 4
    orbits_comp = symmetry_index_orbits(Float64, xy, comp)
    @test fundamental_size(orbits_comp) == N÷4
    for b in 1:fundamental_size(orbits_comp)
        @test length(findall(==(b), orbits_comp.orbit_of)) == 4
    end
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
end

