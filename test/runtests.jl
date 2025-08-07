using BilliardGeometry
using Test
using StaticArrays

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
