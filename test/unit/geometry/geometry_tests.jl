# Unit tests for src/geometry/geometry.jl: is_inside (billiard, domain, curve, array forms), domain_fun, domain_gradient_vector.
# Fixtures: the quarter circle, the rectangle, a right triangle with closed-form half-planes, the stadium and the three hole billiards; helpers use the prefix `gm_`.

using StaticArrays: SVector

# Absolute tolerance for domain functions and their gradients (worst observed 5.6e-17, margin 180x).
const gm_ATOL = 1e-14
# Absolute tolerance for Float32 domain functions (worst observed 0, set to 8 eps(Float32)).
const gm_ATOL_F32 = 1f-6

# Quarter circle of radius 1 (CircleBilliard fundamental domain): arc, y-wall (x = 0), x-wall (y = 0).
gm_arc() = CircleBilliard(1.0).fundamental_domain.boundary[1]
# Right triangle (0,0), (2,0), (0,1), edges (0,0)->(2,0), (2,0)->(0,1), (0,1)->(0,0).
gm_triangle() = PolygonBilliard([SVector(0.0, 0.0), SVector(2.0, 0.0), SVector(0.0, 1.0)])

@testset "domain_fun (circle arc)" begin
    arc = gm_arc()
    # |p - c| - R: negative inside, zero on, positive outside
    @test domain_fun(arc, SVector(0.5, 0.5)) ≈ sqrt(0.5) - 1 atol = gm_ATOL
    @test domain_fun(arc, SVector(2.0, 0.0)) ≈ 1.0 atol = gm_ATOL
    @test domain_fun(arc, SVector(1.0, 0.0)) == 0.0
    @test domain_fun(arc, SVector(0.6, 0.8)) ≈ 0.0 atol = gm_ATOL
    # array form: one value per point, same shape
    @test domain_fun(arc, [SVector(0.5, 0.5), SVector(2.0, 0.0)]) ≈ [sqrt(0.5) - 1, 1.0] atol = gm_ATOL
    @test domain_fun(arc, reshape([SVector(2.0, 0.0), SVector(1.0, 0.0), SVector(0.5, 0.5), SVector(3.0, 0.0)], 2, 2)) ≈ [1.0 sqrt(0.5)-1; 0.0 2.0] atol = gm_ATOL
    @test (@inferred domain_fun(arc, SVector(0.6, 0.8))) isa Float64

    # off-center arc of radius 2 about (1, 1): the point (4, 1) is at distance 3
    off = CircleSegment(2.0, pi / 2, 0.0, [1.0, 1.0])
    @test domain_fun(off, SVector(4.0, 1.0)) ≈ 1.0 atol = gm_ATOL
    @test domain_fun(off, SVector(1.0, 1.0)) ≈ -2.0 atol = gm_ATOL

    # Float32
    arc32 = CircleBilliard(1f0).fundamental_domain.boundary[1]
    @test domain_fun(arc32, SVector(2f0, 0f0)) isa Float32
    @test domain_fun(arc32, SVector(2f0, 0f0)) ≈ 1f0 atol = gm_ATOL_F32
end

@testset "domain_fun (line)" begin
    # the value is the signed distance times the segment length, negative on the inner side;
    # RectangleBilliard(1, 0.6) fundamental domain [0, 1] x [0, 0.6]: right wall x = 1 (length 0.6), top y = 0.6 (length 1),
    # y-wall x = 0 (length 0.6), x-wall y = 0 (length 1)
    right, top, ywall, xwall = RectangleBilliard(1.0, 0.6).fundamental_domain.boundary
    p = SVector(0.5, 0.3)
    q = SVector(1.1, 0.3)
    @test domain_fun(right, p) ≈ (0.5 - 1.0) * 0.6 atol = gm_ATOL
    @test domain_fun(top, p) ≈ (0.3 - 0.6) * 1.0 atol = gm_ATOL
    @test domain_fun(ywall, p) ≈ -0.5 * 0.6 atol = gm_ATOL
    @test domain_fun(xwall, p) ≈ -0.3 * 1.0 atol = gm_ATOL
    # beyond the right wall the value changes sign
    @test domain_fun(right, q) ≈ 0.1 * 0.6 atol = gm_ATOL
    @test domain_fun(right, SVector(1.0, 0.2)) == 0.0

    # array form
    @test domain_fun(right, [p, q]) ≈ [-0.3, 0.06] atol = gm_ATOL
    @test domain_fun(xwall, [p, q]) ≈ [-0.3, -0.3] atol = gm_ATOL
    @test (@inferred domain_fun(right, p)) isa Float64

    # the orientation flips the sign
    flipped = LineSegment([0.0, 0.0], [1.0, 0.0]; orientation = -1)
    @test domain_fun(LineSegment([0.0, 0.0], [1.0, 0.0]), SVector(0.5, 0.5)) ≈ -0.5 atol = gm_ATOL
    @test domain_fun(flipped, SVector(0.5, 0.5)) ≈ 0.5 atol = gm_ATOL

    # right triangle: half-plane functions of the three edges at (0.5, 0.25)
    e1, e2, e3 = gm_triangle().fundamental_domain.boundary
    t = SVector(0.5, 0.25)
    @test domain_fun(e1, t) ≈ -0.5 atol = gm_ATOL
    @test domain_fun(e2, t) ≈ -1.0 atol = gm_ATOL
    @test domain_fun(e3, t) ≈ -0.5 atol = gm_ATOL
    @test domain_fun(e2, SVector(1.5, 0.5)) ≈ 0.5 atol = gm_ATOL

    # Float32
    line32 = LineSegment(Float32[0, 0], Float32[1, 0])
    @test domain_fun(line32, SVector(0.5f0, 0.5f0)) isa Float32
    @test domain_fun(line32, SVector(0.5f0, 0.5f0)) ≈ -0.5f0 atol = gm_ATOL_F32
end

@testset "is_inside (curve)" begin
    arc = gm_arc()
    @test is_inside(arc, SVector(0.5, 0.5)) === true
    @test is_inside(arc, SVector(0.8, 0.7)) === false
    @test is_inside(arc, SVector(2.0, 0.0)) === false
    # the boundary counts as inside (domain_fun <= 0)
    @test is_inside(arc, SVector(1.0, 0.0)) === true
    @test (@inferred is_inside(arc, SVector(0.6, 0.8))) isa Bool

    # array form: one Bool per point
    @test is_inside(arc, [SVector(0.5, 0.5), SVector(2.0, 0.0)]) == [true, false]
    @test is_inside(arc, reshape([SVector(0.5, 0.5), SVector(2.0, 0.0), SVector(0.8, 0.7), SVector(0.1, 0.1)], 2, 2)) == [true false; false true]

    # a line: with orientation 1 the left side of the directed edge is the inside, orientation -1 swaps the sides
    line = LineSegment([0.0, 0.0], [1.0, 0.0])
    @test is_inside(line, SVector(0.5, 0.5)) === true
    @test is_inside(line, SVector(0.5, -0.5)) === false
    @test is_inside(LineSegment([0.0, 0.0], [1.0, 0.0]; orientation = -1), SVector(0.5, 0.5)) === false
end

@testset "is_inside (domain)" begin
    # one entry per boundary curve of the quarter circle: arc, y-wall, x-wall
    dom = CircleBilliard(1.0).fundamental_domain
    @test is_inside(dom, SVector(0.5, 0.5)) == [true, true, true]
    @test is_inside(dom, SVector(0.8, 0.7)) == [false, true, true]
    @test is_inside(dom, SVector(-0.1, 0.2)) == [true, false, true]
    @test is_inside(dom, SVector(0.2, -0.1)) == [true, true, false]

    # array form: rows are the points, columns the curves
    pts = [SVector(0.5, 0.5), SVector(1.2, 0.1), SVector(-0.1, 0.1)]
    @test is_inside(dom, pts) == [true true true; false true true; true false true]

    # a rectangle domain has four curves; a point inside satisfies all, a point beyond the top wall only the top one fails
    rdom = RectangleBilliard(1.0, 0.6).fundamental_domain
    @test is_inside(rdom, SVector(0.5, 0.3)) == [true, true, true, true]
    @test is_inside(rdom, SVector(0.5, 0.7)) == [true, false, true, true]
end

@testset "is_inside (billiard)" begin
    # quarter circle: the intersection of the arc and the two walls
    c = CircleBilliard(1.0)
    @test is_inside(c, SVector(0.5, 0.5)) === true
    @test is_inside(c, SVector(0.3, 0.0)) === true
    @test is_inside(c, SVector(0.0, 0.0)) === true
    @test is_inside(c, SVector(0.8, 0.7)) === false      # beyond the arc
    @test is_inside(c, SVector(-0.1, 0.2)) === false     # beyond the y-wall
    @test is_inside(c, SVector(0.2, -0.1)) === false     # beyond the x-wall
    @test is_inside(c, SVector(1.2, 0.0)) === false

    # array form keeps the shape of the input
    @test is_inside(c, [SVector(0.5, 0.5), SVector(1.2, 0.0)]) == [true, false]
    @test is_inside(c, reshape([SVector(0.5, 0.5), SVector(-0.1, 0.2), SVector(0.8, 0.7), SVector(0.3, 0.1)], 2, 2)) == [true false; false true]

    # rectangle [0, 1] x [0, 0.6]
    r = RectangleBilliard(1.0, 0.6)
    @test is_inside(r, SVector(0.5, 0.3)) === true
    @test is_inside(r, SVector(1.1, 0.3)) === false
    @test is_inside(r, SVector(0.5, 0.7)) === false
    @test is_inside(r, SVector(-0.1, 0.3)) === false
    @test is_inside(r, SVector(0.5, -0.1)) === false

    # triangle (0,0), (2,0), (0,1)
    tri = gm_triangle()
    @test is_inside(tri, SVector(0.5, 0.25)) === true
    @test is_inside(tri, SVector(1.5, 0.5)) === false    # beyond the hypotenuse x + 2y = 2
    @test is_inside(tri, SVector(1.0, -0.5)) === false
    @test is_inside(tri, [SVector(0.5, 0.25), SVector(1.5, 0.5)]) == [true, false]

    # stadium of radius 1 and half-width 0.7: the quarter stadium is the union of the circle wedge and the rectangle
    st = StadiumBilliard(0.7)
    @test is_inside(st, SVector(0.3, 0.5)) === true      # rectangle part
    @test is_inside(st, SVector(1.2, 0.3)) === true      # wedge part
    @test is_inside(st, SVector(1.0, 0.9)) === true      # |p - (0.7, 0)| = 0.95 < 1
    @test is_inside(st, SVector(1.5, 0.9)) === false     # |p - (0.7, 0)| = 1.2 > 1
    @test is_inside(st, SVector(-0.1, 0.3)) === false
    @test is_inside(st, SVector(0.3, 1.1)) === false
end

@testset "is_inside (AnnularBilliard)" begin
    # AnnularBilliard(0.3, 0.2): unit disk without the disk of radius 0.3 centered at (0.2, 0)
    an = AnnularBilliard(0.3, 0.2)
    @test is_inside(an, SVector(0.9, 0.0)) === true
    @test is_inside(an, SVector(-0.7, 0.0)) === true
    @test is_inside(an, SVector(0.55, 0.0)) === true     # 0.35 from the hole center
    @test is_inside(an, SVector(0.2, 0.0)) === false     # hole center
    @test is_inside(an, SVector(0.0, 0.0)) === false     # 0.2 from the hole center
    @test is_inside(an, SVector(0.4, 0.0)) === false
    @test is_inside(an, SVector(1.1, 0.0)) === false     # outside the outer circle
    @test is_inside(an, [SVector(0.9, 0.0), SVector(0.2, 0.0)]) == [true, false]
end

@testset "is_inside (MultiplyConnectedDomain)" begin
    # rectangle [0,3] x [0,2] without [1,1.5] x [0.5,1.5]: ring, hole, outside
    rect = hf_rect_with_hole()
    dom = rect.fundamental_domain
    @test is_inside(dom, SVector(0.5, 0.5)) === true
    @test is_inside(dom, SVector(2.0, 1.0)) === true
    @test is_inside(dom, SVector(1.25, 1.0)) === false     # in the hole
    @test is_inside(dom, SVector(3.5, 1.0)) === false      # outside
    @test is_inside(dom, [SVector(0.5, 0.5), SVector(1.25, 1.0), SVector(3.5, 1.0)]) == [true, false, false]
    # the billiard method dispatches to the domain method: the hole is not a filled region
    @test is_inside(rect, SVector(1.25, 1.0)) === false
    @test is_inside(rect, SVector(0.5, 0.5)) === true
    @test is_inside(rect, [SVector(0.5, 0.5), SVector(1.25, 1.0)]) == [true, false]

    # a point in no outer subdomain is outside even if it is in no hole
    ann = hf_concentric_annulus()
    @test is_inside(ann.fundamental_domain, SVector(2.0, 0.0)) === false
    @test is_inside(ann.fundamental_domain, SVector(0.0, 0.0)) === false
    @test is_inside(ann.fundamental_domain, SVector(0.6, 0.0)) === true
    @test is_inside(hf_eccentric_annulus(), SVector(0.3, 0.0)) === false   # hole center
    @test is_inside(hf_eccentric_annulus(), SVector(-0.5, 0.0)) === true

    # D2 billiard: the outer region is the quarter square, the asymmetric hole of radius 0.2 about (0.3, 0.4) is in the full frame
    d2 = hf_d2_with_hole()
    @test is_inside(d2, SVector(0.3, 0.4)) === false       # hole center
    @test is_inside(d2, SVector(0.9, 0.1)) === true
    @test is_inside(d2, SVector(-0.9, -0.1)) === false     # outside the fundamental quarter
    # inside the fundamental quarter the billiard test is "not in the hole"
    for x in 0.05:0.1:0.95, y in 0.05:0.1:0.95
        p = SVector(x, y)
        @test is_inside(d2, p) == (hypot(x - 0.3, y - 0.4) > 0.2)
    end
    # the hole is tested at the point itself, not at its symmetry images: the mirror image of the hole center is no hole point
    hole = d2.fundamental_domain.holes[1]
    @test is_inside(hole, SVector(0.3, 0.4)) === true
    @test is_inside(hole, SVector(-0.3, 0.4)) === false
end

@testset "is_inside (SinaiBilliard)" begin
    # the square [-1, 1]^2 without the disk of radius 0.3 about the origin
    sn = SinaiBilliard(1.0, 0.3)
    @test is_inside(sn, SVector(0.5, 0.5)) === true
    @test is_inside(sn, SVector(-0.5, -0.5)) === true
    @test is_inside(sn, SVector(-0.9, 0.1)) === true
    @test is_inside(sn, SVector(0.9, 0.1)) === true
    @test is_inside(sn, SVector(0.25, 0.25)) === true    # distance 0.354 from the origin
    @test is_inside(sn, SVector(0.31, 0.0)) === true
    @test is_inside(sn, SVector(0.1, 0.1)) === false     # inside the hole
    @test is_inside(sn, SVector(-0.1, -0.2)) === false
    @test is_inside(sn, SVector(0.2, 0.2)) === false     # distance 0.283 from the origin
    @test is_inside(sn, SVector(0.29, 0.0)) === false
    @test is_inside(sn, SVector(1.1, 0.5)) === false
    @test is_inside(sn, [SVector(0.5, 0.5), SVector(0.1, 0.1)]) == [true, false]
end

@testset "is_inside (SquareWithinSquareBilliard)" begin
    # the square [-1, 1]^2 without the two square holes [0.375, 0.625]^2 and [-0.625, -0.375]^2
    sws = SquareWithinSquareBilliard(1.0)
    # points in a hole or outside the outer square are outside
    @test is_inside(sws, SVector(0.55, 0.45)) === false
    @test is_inside(sws, SVector(0.5, 0.4)) === false
    @test is_inside(sws, SVector(-0.5, -0.6)) === false
    @test is_inside(sws, SVector(1.2, 0.1)) === false
    @test is_inside(sws, SVector(-1.2, 0.1)) === false

    # points of the billiard, on both sides of the diagonal, are inside
    @testset "points of the billiard" begin
        @test is_inside(sws, SVector(0.9, 0.2)) === true
        @test is_inside(sws, SVector(0.7, 0.1)) === true
        @test is_inside(sws, SVector(0.3, 0.1)) === true
        @test is_inside(sws, SVector(0.9, 0.5)) === true
        @test is_inside(sws, SVector(0.5, 0.3)) === true
        @test is_inside(sws, SVector(-0.2, -0.5)) === true
        @test is_inside(sws, SVector(0.0, -0.9)) === true
        @test is_inside(sws, SVector(0.2, 0.9)) === true
        @test is_inside(sws, SVector(-0.9, 0.1)) === true
    end
    @testset "array form" begin
        @test is_inside(sws, [SVector(0.9, 0.2), SVector(0.55, 0.45), SVector(0.3, 0.1)]) == [true, false, true]
    end
end

@testset "domain_gradient_vector" begin
    arc = gm_arc()
    # gradient of |p| - R is the radial unit vector
    @test domain_gradient_vector(arc, SVector(0.6, 0.8)) ≈ SVector(0.6, 0.8) atol = gm_ATOL
    @test domain_gradient_vector(arc, SVector(0.3, 0.4)) ≈ SVector(0.6, 0.8) atol = gm_ATOL
    @test domain_gradient_vector(arc, SVector(0.0, 2.0)) ≈ SVector(0.0, 1.0) atol = gm_ATOL
    @test (@inferred domain_gradient_vector(arc, SVector(0.6, 0.8))) isa SVector{2,Float64}
    # off-center arc about (1, 1): unit vector from the center
    @test domain_gradient_vector(CircleSegment(2.0, pi / 2, 0.0, [1.0, 1.0]), SVector(4.0, 5.0)) ≈ SVector(0.6, 0.8) atol = gm_ATOL

    # array form: one gradient per point
    gs = domain_gradient_vector(arc, [SVector(0.6, 0.8), SVector(0.0, 2.0)])
    @test gs isa Vector{SVector{2,Float64}}
    @test gs ≈ [SVector(0.6, 0.8), SVector(0.0, 1.0)] atol = gm_ATOL

    # a line has the constant gradient (y1 - y0, -(x1 - x0)) * orientation, independent of the point;
    # right triangle (0,0), (2,0), (0,1)
    e1, e2, e3 = gm_triangle().fundamental_domain.boundary
    p = SVector(0.5, 0.25)
    @test domain_gradient_vector(e1, p) ≈ SVector(0.0, -2.0) atol = gm_ATOL
    @test domain_gradient_vector(e2, p) ≈ SVector(1.0, 2.0) atol = gm_ATOL
    @test domain_gradient_vector(e3, p) ≈ SVector(-1.0, 0.0) atol = gm_ATOL
    @test domain_gradient_vector(e2, SVector(9.0, -3.0)) ≈ SVector(1.0, 2.0) atol = gm_ATOL
    @test domain_gradient_vector(e2, [p, SVector(9.0, -3.0)]) ≈ [SVector(1.0, 2.0), SVector(1.0, 2.0)] atol = gm_ATOL
    # the gradient of the reversed orientation is the negative one
    @test domain_gradient_vector(LineSegment([0.0, 0.0], [2.0, 0.0]; orientation = -1), p) ≈ SVector(0.0, 2.0) atol = gm_ATOL

    # the gradient is normal to the line: orthogonal to the edge direction
    @test sum(domain_gradient_vector(e2, p) .* (e2.pt1 - e2.pt0)) ≈ 0.0 atol = gm_ATOL

    # rectangle right wall x = 1 (length 0.6)
    right = RectangleBilliard(1.0, 0.6).fundamental_domain.boundary[1]
    @test domain_gradient_vector(right, SVector(0.5, 0.3)) ≈ SVector(0.6, 0.0) atol = gm_ATOL

    # Float32
    arc32 = CircleBilliard(1f0).fundamental_domain.boundary[1]
    g32 = domain_gradient_vector(arc32, SVector(0.6f0, 0.8f0))
    @test g32 isa SVector{2,Float32}
    @test g32 ≈ SVector(0.6f0, 0.8f0) atol = gm_ATOL_F32
end
