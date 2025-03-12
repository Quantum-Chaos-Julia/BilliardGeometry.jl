using BilliardGeometry
using Test
using StaticArrays

@testset "linesegment.jl" begin
    line = LineSegment([0.0,0.0],[1.0,1.0])
    @test curve(line, [0.0,1.0]) == [SVector(0.0,0.0),SVector(1.0,1.0)]
    @test domain_fun(line, SVector{2,Float64}([0.1,0.2])) < 0.0
    @test domain_fun(line, SVector{2,Float64}([0.1,0.0])) > 0.0
    @test arc_length(line, line.pt1) == line.length
end

@testset "circlesegment.jl" begin
    circle = LineSegment([0.0,0.0],[1.0,1.0])
    @test curve(circle, [0.0,1.0]) == [SVector(0.0,0.0),SVector(1.0,1.0)]
    @test domain_fun(circle, SVector{2,Float64}([0.1,0.2])) < 0.0
    @test domain_fun(circle, SVector{2,Float64}([0.1,0.0])) > 0.0
    @test arc_length(circle, circle.pt1) == circle.length
end
