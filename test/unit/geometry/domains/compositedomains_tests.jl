# Unit tests for src/geometry/domains/compositedomains.jl: SimpleDomain, CompositeDomain, _domain_eltype, reset_ids!.
# Fixtures are a square, a triangle and a circle wedge with chosen ids; helpers use the prefix `cd_`.

using StaticArrays: SVector
using BilliardGeometry: _domain_eltype

cd_square(id; T = Float64) = Polygon([T[0, 0], T[1, 0], T[1, 1], T[0, 1]], id)
cd_triangle(id; T = Float64) = Polygon([T[0, 0], T[1, 0], T[0, 1]], id)
cd_wedge(id; T = Float64) = CircleWedge(one(T), T(pi / 2), zero(T), T(0.5), zero(T), id;
    bcs = AbsBoundaryCondition[SpecularReflection(), Transparent(1), SymmetryWall(3, 4)])

@testset "SimpleDomain" begin
    sq = cd_square(3)
    dom = SimpleDomain{Float64}(sq.boundary, sq.corners, 5)
    @test dom isa SimpleDomain{Float64}
    @test dom isa AbsSimpleDomain
    @test fieldnames(SimpleDomain) == (:boundary, :corners, :id)
    @test dom.boundary === sq.boundary
    @test dom.corners === sq.corners
    @test dom.id == 5

    # the numeric type is inferred from the corners when it is not given
    @test SimpleDomain(sq.boundary, sq.corners, 5) isa SimpleDomain{Float64}
    sq32 = cd_square(1; T = Float32)
    @test SimpleDomain(sq32.boundary, sq32.corners, 1) isa SimpleDomain{Float32}

    # the id is an integer
    @test_throws InexactError SimpleDomain{Float64}(sq.boundary, sq.corners, 5.5)
end

@testset "CompositeDomain" begin
    sq, tri = cd_square(1), cd_triangle(2)
    dom = CompositeDomain(AbsSimpleDomain[sq, tri])
    @test dom isa CompositeDomain
    @test dom isa AbsCompositeDomain
    @test fieldnames(CompositeDomain) == (:subdomains,)
    @test dom.subdomains isa Vector{AbsSimpleDomain}
    @test length(dom.subdomains) == 2
    @test dom.subdomains[1] === sq
    @test dom.subdomains[2] === tri

    # subdomains of different types and of one type are both stored as Vector{AbsSimpleDomain}
    mixed = CompositeDomain([sq, cd_wedge(2)])
    @test mixed.subdomains isa Vector{AbsSimpleDomain}
    @test [typeof(d) for d in mixed.subdomains] == [Polygon{Float64}, CircleWedge{Float64}]
    @test CompositeDomain([sq, tri]).subdomains isa Vector{AbsSimpleDomain}
end

@testset "_domain_eltype" begin
    sq = cd_square(1)
    @test _domain_eltype(SimpleDomain{Float64}(sq.boundary, sq.corners, 1)) === Float64
    @test _domain_eltype(SimpleDomain{Float32}(sq.boundary, SVector{2,Float32}.(sq.corners), 1)) === Float32
    @test _domain_eltype(sq) === Float64
    @test _domain_eltype(cd_square(1; T = Float32)) === Float32
    @test _domain_eltype(cd_wedge(1)) === Float64
    @test _domain_eltype(cd_wedge(1; T = Float32)) === Float32
    @test (@inferred _domain_eltype(sq)) === Float64
    @test (@inferred _domain_eltype(cd_wedge(1; T = Float32))) === Float32

    # a composite domain takes the type of its first subdomain
    @test _domain_eltype(CompositeDomain([cd_wedge(1), cd_triangle(2)])) === Float64
    @test _domain_eltype(CompositeDomain(AbsSimpleDomain[cd_square(1; T = Float32), cd_triangle(2; T = Float32)])) === Float32
    @test _domain_eltype(CompositeDomain(AbsSimpleDomain[cd_triangle(1; T = Float32), cd_wedge(2)])) === Float32

    # the type is read from the first subdomain, so an empty composite domain has none
    @test_throws BoundsError _domain_eltype(CompositeDomain(AbsSimpleDomain[]))
end

@testset "reset_ids!" begin
    sq, tri, wedge = cd_square(9), cd_triangle(5), cd_wedge(2)
    dom = CompositeDomain(AbsSimpleDomain[sq, tri, wedge])
    subdomain_vector = dom.subdomains

    # ids become 1, 2, 3 in storage order, and the domain itself is returned
    @test reset_ids!(dom) === dom
    @test [d.id for d in dom.subdomains] == [1, 2, 3]
    # the subdomain vector is updated in place
    @test dom.subdomains === subdomain_vector

    # subdomain types and geometry are kept
    @test [typeof(d) for d in dom.subdomains] == [Polygon{Float64}, Polygon{Float64}, CircleWedge{Float64}]
    @test dom.subdomains[1].corners == sq.corners
    @test dom.subdomains[2].corners == tri.corners
    @test dom.subdomains[2].angles == tri.angles
    @test dom.subdomains[3].corners == wedge.corners
    @test dom.subdomains[3].center == wedge.center
    @test length(dom.subdomains[1].boundary) == 4

    # the original subdomain values are immutable and keep their own ids
    @test (sq.id, tri.id, wedge.id) == (9, 5, 2)

    # resetting twice changes nothing
    reset_ids!(dom)
    @test [d.id for d in dom.subdomains] == [1, 2, 3]

    # a composite domain with one subdomain gets id 1, an empty one stays empty
    single = CompositeDomain(AbsSimpleDomain[cd_wedge(4)])
    reset_ids!(single)
    @test single.subdomains[1].id == 1
    @test isempty(reset_ids!(CompositeDomain(AbsSimpleDomain[])).subdomains)

    # element types of the subdomains are kept
    dom32 = CompositeDomain(AbsSimpleDomain[cd_square(7; T = Float32), cd_wedge(8; T = Float32)])
    reset_ids!(dom32)
    @test [d.id for d in dom32.subdomains] == [1, 2]
    @test _domain_eltype(dom32) === Float32

    # only composite domains are renumbered
    @test_throws MethodError reset_ids!(cd_square(3))
end
