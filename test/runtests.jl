using BilliardGeometry
using Test
using StaticArrays

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
    @test s_of_t(0.5)≈0.5*crv.length atol=1e-2
    @test t_of_s(0.5*crv.length)≈0.5 atol=1e-2
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
    @test s_of_t(0.5)≈0.5*crv.length atol=1e-2
    @test t_of_s(0.5*crv.length)≈0.5 atol=1e-2
end

@testset "limaconsegment.jl" begin
    a=0.2
    crv=LimaconSegment(a)
    @test arc_length(crv,curve(crv,1.0))≈crv.length
    # gradient functions test
    @test domain_gradient_vector(crv,curve(crv,0.25))≈SVector{2,Float64}([0.0,1.0])
    # domain functions test
    testpt1=SVector{2,Float64}([0.1,0.2])
    testpt2=SVector{2,Float64}([-0.3,0.0])
    @test domain_fun(crv,testpt1)<0.0
    @test domain_fun(crv,testpt2)>0.0
    @test is_inside(crv,[testpt1,testpt2])==[true,false]
    # arclength test
    s_of_t,t_of_s=construct_arc_length_interpolation(crv)
    @test s_of_t(0.5)≈0.5*crv.length atol=1e-2
    @test t_of_s(0.5*crv.length)≈0.5 atol=1e-2
end

@testset "polarsegment.jl" begin
    N=10
    center=SVector{2,Float64}([0.5,0.5])
    coeffs=[rand() for _ in 1:N]
    crv=PolarSegment(coeffs,center=center)
    # curve functions test
    @test arc_length(crv,curve(crv,1.0))≈crv.length
    # arclength test
    s_of_t,t_of_s=construct_arc_length_interpolation(crv)
    @test s_of_t(0.5)≈0.5*crv.length atol=1e-2
    @test t_of_s(0.5*crv.length)≈0.5 atol=1e-2
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
