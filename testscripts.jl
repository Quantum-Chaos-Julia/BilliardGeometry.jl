using Revise
using BilliardGeometry
using StaticArrays
using CircularArrays
using BilliardGeometryPlotting
using CairoMakie

line = LineSegment([0.0,0.0],[1.0,1.0])
line.bc

curve(line, [0.0,1.0])
domain_fun(line, SVector{2,Float64}([0.1,0.2]))

triangle = Polygon([[0.0,0.0],[1.0,0.0],[0.0,1.0]],1)
triangle.boundary

boundary = Vector{AbsCurve}()
Polygon(boundary,corners,1)