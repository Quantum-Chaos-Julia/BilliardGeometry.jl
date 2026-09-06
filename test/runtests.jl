using BilliardGeometry
using Test
using StaticArrays
using LinearAlgebra

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
        @test issymmetric(R0)
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
    @test symmetry_node_multiple(XAxisReflection()) == 4
    @test symmetry_node_multiple(YAxisReflection()) == 4
    @test symmetry_node_multiple(XYAxisReflection()) == 4
    @test symmetry_node_multiple(NFoldRotation(6,1)) == 6

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
    orbits_rot = symmetry_index_orbits(Float64, xy[1:N-mod(N,n)], NFoldRotation(n,1))
    Nr = N-mod(N,n)
    @test fundamental_size(orbits_rot) == Nr÷n
    for b in 1:fundamental_size(orbits_rot)
        @test length(findall(==(b), orbits_rot.orbit_of)) == n
    end
end
