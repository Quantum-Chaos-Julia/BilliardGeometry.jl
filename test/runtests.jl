using BilliardGeometry
using Test
using StaticArrays
using LinearAlgebra

@testset "linesegment.jl" begin
    pt0,pt1=[0.0,0.0],[1.0,1.0]
    crv=LineSegment(pt0,pt1) 
    # curve functions test
    @test curve(crv,0.0)==pt0
    @test curve(crv,1.0)==pt1
    @test arc_length(crv,crv.pt1)==crv.length
    # gradient functions test
    @test domain_gradient_vector(crv,curve(crv,0.5))==SVector{2,Float64}([1.0,-1.0])
    # domain functions test
    testpt1=SVector{2,Float64}([0.1,0.2])
    testpt2=SVector{2,Float64}([0.1,0.0])
    @test domain_fun(crv,testpt1)<0.0
    @test domain_fun(crv,testpt2)>0.0
    @test is_inside(crv,[testpt1,testpt2])==[true,false]
    # arclength test
    s_of_t,t_of_s=construct_arc_length_interpolation(crv)
    @test s_of_t(0.5)≈0.5*crv.length atol=1e-8
    @test t_of_s(0.5*crv.length)≈0.5 atol=1e-8
end

@testset "circlesegment.jl" begin
    R=1.0
    arc_angle=0.5*pi
    shift_angle=0.0
    center=[0.0,0.0]
    crv=CircleSegment(R,arc_angle,shift_angle,center)
    # curve functions test
    @test all(isapprox.(curve(crv,0.0),SVector{2,Float64}([1.0,0.0]);atol=1e-8))
    @test all(isapprox.(curve(crv,1.0),SVector{2,Float64}([0.0,1.0]);atol=1e-8))
    @test arc_length(crv,curve(crv,1.0))≈crv.length
    # gradient functions test
    @test domain_gradient_vector(crv,curve(crv,0.5))≈SVector{2,Float64}([sqrt(2)/2,sqrt(2)/2])
    # domain functions test
    testpt1=SVector{2,Float64}([0.1,0.2])
    testpt2=SVector{2,Float64}([1.1,1.0])
    @test domain_fun(crv,testpt1)<0.0
    @test domain_fun(crv,testpt2)>0.0
    @test is_inside(crv,[testpt1,testpt2])==[true,false]
    # arclength test
    s_of_t,t_of_s=construct_arc_length_interpolation(crv)
    @test s_of_t(0.5)≈0.5*crv.length atol=1e-8
    @test t_of_s(0.5*crv.length)≈0.5 atol=1e-8
end

@testset "limacon_polarsegment.jl" begin
    a=0.2
    rfun=φ->one(Float64)+a*cos(φ)
    crv=PolarSegment(Float64,rfun;arc_angle=2π,shift_angle=0.0,center=SVector(0.0,0.0),orientation=1)
    @test arc_length(crv,curve(crv,1-1e-4))≈crv.length atol=1e-1
    p=curve(crv,0.25);g=domain_gradient_vector(crv,p)
    ĝ=g/norm(g)
    @test norm(ĝ)≈1 atol=1e-2
    ϵ=1e-3
    @test domain_fun(crv,p-ϵ*ĝ)<0
    @test domain_fun(crv,p+ϵ*ĝ)>0
    @test is_inside(crv,[p-ϵ*ĝ,p+ϵ*ĝ])==[true,false]
    # chebyshev interpolation test
    s_of_t,t_of_s=construct_arc_length_interpolation(Float64,crv;q=5.0,p=12)
    @test s_of_t(0.5)≈0.5*crv.length atol=1e-8
    @test t_of_s(0.5*crv.length)≈0.5 atol=1e-8
    # cubic spline interpolation test
    s_of_t_cs,t_of_s_cs=construct_arc_length_interpolation(crv;method=:cubic_spline,n_samples=2000)
    @test s_of_t_cs(0.5)≈0.5*crv.length atol=1e-8
    @test t_of_s_cs(0.5*crv.length)≈0.5 atol=1e-8
end

@testset "ellipse_polarsegment.jl" begin
    a=1.0;b=0.6
    rfun=φ->(a*b)/sqrt((b*cos(φ))^2+(a*sin(φ))^2)
    crv=PolarSegment(Float64,rfun;arc_angle=2π,shift_angle=0.0,center=SVector(0.0,0.0),orientation=1)
    @test arc_length(crv,curve(crv,1-1e-4))≈crv.length atol=1e-1
    p=curve(crv,0.33);g=domain_gradient_vector(crv,p)
    ĝ=g/norm(g)
    @test norm(ĝ)≈1 atol=1e-2
    ϵ=1e-3
    @test domain_fun(crv,p-ϵ*ĝ)<0
    @test domain_fun(crv,p+ϵ*ĝ)>0
    @test is_inside(crv,[p-ϵ*ĝ,p+ϵ*ĝ])==[true,false]
    # chebyshev interpolation test
    s_of_t,t_of_s=construct_arc_length_interpolation(Float64,crv;q=5.0,p=12)
    @test s_of_t(0.5)≈0.5*crv.length atol=1e-8
    @test t_of_s(0.5*crv.length)≈0.5 atol=1e-8
    # cubic spline interpolation test
    s_of_t_cs,t_of_s_cs=construct_arc_length_interpolation(crv;method=:cubic_spline,n_samples=2000)
    @test s_of_t_cs(0.5)≈0.5*crv.length atol=1e-8
    @test t_of_s_cs(0.5*crv.length)≈0.5 atol=1e-8
end

@testset "prosen_polarsegment.jl" begin
    a=0.08
    rfun=φ->one(Float64)+a*cos(4*φ)
    crv=PolarSegment(Float64,rfun;arc_angle=2π,shift_angle=0.0,center=SVector(0.0,0.0),orientation=1)
    @test arc_length(crv,curve(crv,1-1e-4))≈crv.length atol=1e-1
    p=curve(crv,0.41);g=domain_gradient_vector(crv,p)
    ĝ=g/norm(g)
    @test norm(ĝ)≈1 atol=1e-2
    ϵ=1e-3
    @test domain_fun(crv,p-ϵ*ĝ)<0
    @test domain_fun(crv,p+ϵ*ĝ)>0
    @test is_inside(crv,[p-ϵ*ĝ,p+ϵ*ĝ])==[true,false]
    s_of_t,t_of_s=construct_arc_length_interpolation(Float64,crv;q=5.0,p=12)
    # chebyshev interpolation test
    @test s_of_t(0.5)≈0.5*crv.length atol=1e-8
    @test t_of_s(0.5*crv.length)≈0.5 atol=1e-8
    # cubic spline interpolation test
    s_of_t_cs,t_of_s_cs=construct_arc_length_interpolation(crv;method=:cubic_spline,n_samples=2000)
    @test s_of_t_cs(0.5)≈0.5*crv.length atol=1e-8
    @test t_of_s_cs(0.5*crv.length)≈0.5 atol=1e-8
end

@testset "polarsegment.jl" begin
    N=20;coeffs=0.1.*rand(Float64,N)
    seg=PolarSegment(coeffs)
    @test arc_length(seg,curve(seg,1-1e-4))≈seg.length atol=1e-1
    t=0.25;p=curve(seg,t);g=domain_gradient_vector(seg,p);ĝ=g/norm(g);ϵ=1e-3
    @test domain_fun(seg,p-ϵ*ĝ)<0
    @test domain_fun(seg,p+ϵ*ĝ)>0
    @test is_inside(seg,[p-ϵ*ĝ,p+ϵ*ĝ])==[true,false]
    s_of_t,t_of_s=construct_arc_length_interpolation(Float64,seg;q=5.0,p=12)
    @test s_of_t(0.0)≈0.0 atol=1e-8
    @test s_of_t(1.0)≈seg.length atol=1e-6
    @test t_of_s(0.0)≈0.0 atol=1e-8
    @test t_of_s(seg.length)≈1.0 atol=1e-6
    for tq in (0.1,0.25,0.5,0.77,0.9)
        sq=s_of_t(tq)
        @test t_of_s(sq)≈tq atol=1e-6
    end
    for sq in (0.1,0.33,0.5,0.8,0.95).*seg.length
        tq=t_of_s(sq)
        @test s_of_t(tq)≈sq atol=1e-6
    end
end

@testset "symmetry reflections" begin
    T=Float64
    x0=1.2
    y0=-0.7
    pt=SVector{2,T}(0.3,-2.0)
    n=SVector{2,T}(0.6,-0.8)
    sx=YAxisReflection(x0)          # reflect across x = x0 (flip x)
    sy=XAxisReflection(y0)          # reflect across y = y0 (flip y)
    sxy=XYAxisReflection(x0, y0)    # both flips
    # point reflections
    @test apply_symmetry(sx,pt)==SVector{2,T}(2*x0-pt[1],pt[2])
    @test apply_symmetry(sy,pt)==SVector{2,T}(pt[1],2*y0-pt[2])
    @test apply_symmetry(sxy,pt)==SVector{2,T}(2*x0-pt[1],2*y0-pt[2])
    # involution property: reflecting twice gives original (for points)
    @test apply_symmetry(sx,apply_symmetry(sx,pt))==pt
    @test apply_symmetry(sy,apply_symmetry(sy,pt))==pt
    @test apply_symmetry(sxy,apply_symmetry(sxy,pt))==pt
    # normals: shifts do NOT matter (only sign flips)
    @test apply_symmetry_normals(sx,n)==SVector{2,T}(-n[1],n[2])
    @test apply_symmetry_normals(sy,n)==SVector{2,T}(n[1],-n[2])
    @test apply_symmetry_normals(sxy,n)==-n
    # involution for normals too
    @test apply_symmetry_normals(sx,apply_symmetry_normals(sx,n))==n
    @test apply_symmetry_normals(sy,apply_symmetry_normals(sy,n))==n
    @test apply_symmetry_normals(sxy,apply_symmetry_normals(sxy,n))==n
end

@testset "symmetry rotations" begin
    T=Float64
    c0=SVector{2,T}(0.4,-0.9)
    pt=SVector{2,T}(1.3,2.1)
    n=SVector{2,T}(0.6,-0.8)
    # quarter-turn around c0: (dx,dy)->(-dy,dx)
    rot90=Rotation(4,1,c0)
    pt90=apply_symmetry(rot90,pt)
    dx=pt[1]-c0[1]
    dy=pt[2]-c0[2]
    @test pt90≈SVector{2,T}(c0[1]-dy,c0[2]+dx) atol=1e-14
    # normals rotate about origin (no c0 shift)
    n90=apply_symmetry_normals(rot90,n)
    @test n90≈SVector{2,T}(-n[2],n[1]) atol=1e-14
    # N-fold closure: applying rotation N times returns original
    N=7
    rot1=Rotation(N,1,c0)
    p=pt
    for _ in 1:N
        p=apply_symmetry(rot1,p)
    end
    @test p≈pt atol=1e-12
    v=n
    for _ in 1:N
        v=apply_symmetry_normals(rot1,v)
    end
    @test v≈n atol=1e-12
    # length preservation for normals (direction vectors)
    @test norm(apply_symmetry_normals(rot1,n))≈norm(n) atol=1e-14
end

@testset "vector versions" begin
    T=Float64
    pts=[SVector{2,T}(0.1,0.2),SVector{2,T}(-1.0,3.0),SVector{2,T}(2.0,-4.0)]
    ns=[SVector{2,T}(1.0,0.0),SVector{2,T}(0.0,1.0),SVector{2,T}(0.6,-0.8)]
    x0,y0=0.3,-0.4
    sx=YAxisReflection(x0)
    sy=XAxisReflection(y0)
    sxy=XYAxisReflection(x0,y0)
    pts2=apply_symmetry(sxy,pts)
    ns2=apply_symmetry_normals(sxy,ns)
    @test length(pts2)==length(pts)
    @test length(ns2)==length(ns)
    @test pts2[1]==apply_symmetry(sxy,pts[1])
    @test ns2[3]==apply_symmetry_normals(sxy,ns[3])
    # input vectors not mutated
    @test pts[1]==SVector{2,T}(0.1,0.2)
    @test ns[2]==SVector{2,T}(0.0,1.0)
    # D2 generator sanity (assuming D2_symmetry(x0,y0) returns [YAxisReflection(x0), XYAxisReflection(x0,y0), XAxisReflection(y0)])
    syms=D2_symmetry(x0,y0)
    @test length(syms)==3
    @test apply_symmetry(syms[1],pts[1])==apply_symmetry(YAxisReflection(x0),pts[1])
    @test apply_symmetry(syms[3],pts[1])==apply_symmetry(XAxisReflection(y0),pts[1])
end

@testset "Cn_symmetry generator" begin
    T=Float64
    c0=SVector{2,T}(0.0,0.0)
    n=6
    syms=Cn_symmetry(n,c0)
    @test length(syms)==n
    pt=SVector{2,T}(1.0, 0.0)
    # i=0 should be identity rotation
    @test apply_symmetry(syms[1],pt)≈pt atol=1e-14
    # last element should be rotation by 2π*(n-1)/n, and composing with rotation by 2π/n gives identity
    rot1=syms[2]
    rotm1=syms[end]
    p=apply_symmetry(rot1,apply_symmetry(rotm1,pt))
    @test p≈pt atol=1e-12
end
