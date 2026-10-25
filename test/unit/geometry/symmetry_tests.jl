# Unit tests for the discrete symmetry generators (reflections and rotations) of a billiard:
#   src/geometry/symmetry.jl      (CompositeReflection ==/hash, _sym_matrix, apply_symmetry, apply_symmetry_pb, D2_symmetry,
#                                  NFoldRotation, Cn_symmetry, _with_sym_id)
#   src/geometry/fullboundary.jl  (orientation_reversing)
# Reference objects: the five coordinate reflections acting on the point (x, y), namely X: (x, -y), Y: (-x, y),
# XY: (-x, -y), Diagonal: (y, x) and AntiDiagonal: (-y, -x); and the counterclockwise rotation by 2 pi m / N, which maps
# (1, 0) to (cos(2 pi m / N), sin(2 pi m / N)). Pure reflections have determinant -1 (they reverse the boundary
# orientation), XY and rotations have determinant +1. The registry itself (register_symmetries, symmetry_of) is not
# tested here. No eigen solves, no pinned wavefunction values. Helper constants use the prefix `sy_`.

using StaticArrays: SVector, SMatrix
using LinearAlgebra: det, I
using BilliardGeometry: AbsReflection, _sym_matrix, _with_sym_id

# Absolute tolerance for rotation images (cos(pi/2) is 6e-17 in floating point).
const sy_ATOL = 1e-12
# Absolute tolerance for Float32 comparisons.
const sy_ATOL_F32 = 1f-6

# The five reflection generators with their action on the point (1, 2).
const sy_REFLECTIONS = [(XAxisReflection(), SVector(1.0, -2.0)),
                        (YAxisReflection(), SVector(-1.0, 2.0)),
                        (XYAxisReflection(), SVector(-1.0, -2.0)),
                        (DiagonalReflection(), SVector(2.0, 1.0)),
                        (AntiDiagonalReflection(), SVector(-2.0, -1.0))]

@testset "CompositeReflection equality and hash" begin
    a = CompositeReflection(XAxisReflection(), YAxisReflection())
    b = CompositeReflection(XAxisReflection(), YAxisReflection())
    c = CompositeReflection(YAxisReflection(), XAxisReflection())
    d = CompositeReflection(XAxisReflection())
    @test a !== b
    @test a == b
    @test isequal(a, b)
    @test hash(a) == hash(b)
    @test a != c                       # generator order is part of the identity
    @test a != d
    @test length(unique([a, b])) == 1
    @test Dict(a => 1)[b] == 1

    # constructors: vararg and vector forms build the same generator list
    @test CompositeReflection([XAxisReflection(), YAxisReflection()]) == a
    @test a.reflections isa Vector{AbsReflection}
    @test length(a.reflections) == 2
    # the sym_id is part of each generator's identity
    @test CompositeReflection(XAxisReflection(1)) != CompositeReflection(XAxisReflection(2))
    @test a isa AbsReflection
end

@testset "reflection constructors" begin
    # default sym_id is 0, explicit sym_id is stored
    for T in (XAxisReflection, YAxisReflection, XYAxisReflection, DiagonalReflection, AntiDiagonalReflection)
        @test T().sym_id == 0
        @test T(3).sym_id == 3
        @test T() isa AbsReflection
        @test T() isa AbsSymmetry
    end
end

@testset "_with_sym_id" begin
    for (sym, _) in sy_REFLECTIONS
        tagged = _with_sym_id(sym, 7)
        @test typeof(tagged) == typeof(sym)
        @test tagged.sym_id == 7
    end
    rot = NFoldRotation(5, 2)
    tagged = _with_sym_id(rot, 4)
    @test tagged isa NFoldRotation{Float64}
    @test tagged.sym_id == 4
    @test tagged.order == 5 && tagged.m == 2
    @test tagged.sym_map == rot.sym_map
    @test _with_sym_id(NFoldRotation(5, 2; T = Float32), 4) isa NFoldRotation{Float32}
end

@testset "_sym_matrix" begin
    @test _sym_matrix(XAxisReflection()) == [1.0 0.0; 0.0 -1.0]
    @test _sym_matrix(YAxisReflection()) == [-1.0 0.0; 0.0 1.0]
    @test _sym_matrix(XYAxisReflection()) == [-1.0 0.0; 0.0 -1.0]
    @test _sym_matrix(DiagonalReflection()) == [0.0 1.0; 1.0 0.0]
    @test _sym_matrix(AntiDiagonalReflection()) == [0.0 -1.0; -1.0 0.0]

    # the numeric type follows the requested type, default Float64
    @test _sym_matrix(XAxisReflection()) isa SMatrix{2,2,Float64,4}
    @test _sym_matrix(XAxisReflection(), Float32) isa SMatrix{2,2,Float32,4}
    @test _sym_matrix(DiagonalReflection(), BigFloat) isa SMatrix{2,2,BigFloat,4}

    # every generator is an involution with determinant -1 (pure reflections) or +1 (XY)
    for (sym, _) in sy_REFLECTIONS
        M = _sym_matrix(sym)
        @test M * M == Matrix(I, 2, 2)
        @test det(M) == (sym isa XYAxisReflection ? 1.0 : -1.0)
    end
end

@testset "apply_symmetry (reflections)" begin
    pt = SVector(1.0, 2.0)
    for (sym, image) in sy_REFLECTIONS
        @test apply_symmetry(sym, pt) == image
        # applying a reflection twice is the identity
        @test apply_symmetry(sym, apply_symmetry(sym, pt)) == pt
    end

    # the point type is preserved
    @test apply_symmetry(XAxisReflection(), SVector(1.0f0, 2.0f0)) isa SVector{2,Float32}
    @test apply_symmetry(XAxisReflection(), SVector(big(1.0), big(2.0))) == SVector(big(1.0), big(-2.0))

    # a collection of points is mapped pointwise
    pts = [SVector(1.0, 2.0), SVector(-3.0, 0.5), SVector(0.0, 0.0)]
    @test apply_symmetry(YAxisReflection(), pts) == [SVector(-1.0, 2.0), SVector(3.0, 0.5), SVector(0.0, 0.0)]

    # XY is the composition of the two axis reflections
    @test apply_symmetry(YAxisReflection(), apply_symmetry(XAxisReflection(), pt)) == apply_symmetry(XYAxisReflection(), pt)
end

@testset "orientation_reversing" begin
    @test orientation_reversing(XAxisReflection())
    @test orientation_reversing(YAxisReflection())
    @test orientation_reversing(DiagonalReflection())
    @test orientation_reversing(AntiDiagonalReflection())
    @test !orientation_reversing(XYAxisReflection())
    @test !orientation_reversing(NFoldRotation(4, 1))
    # agrees with the sign of the determinant of the linear map
    for (sym, _) in sy_REFLECTIONS
        @test orientation_reversing(sym) == (det(_sym_matrix(sym)) < 0)
    end
end

@testset "apply_symmetry_pb" begin
    L = 2.5
    s, p = 0.7, 0.3
    # orientation-reversing image k = 1 (sector 2) runs backwards: s -> 2L - s, p -> -p
    @test apply_symmetry_pb(XAxisReflection(), 2, s, p, L) == (2 * L - s, -p)
    # orientation-reversing image k = 3 (sector 4)
    @test apply_symmetry_pb(YAxisReflection(), 4, s, p, L) == (4 * L - s, -p)
    # orientation-preserving images continue the traversal: s -> k L + s, p unchanged
    @test apply_symmetry_pb(XYAxisReflection(), 3, s, p, L) == (2 * L + s, p)
    @test apply_symmetry_pb(NFoldRotation(3, 1), 2, s, p, L) == (L + s, p)
end

@testset "D2_symmetry" begin
    reg = D2_symmetry()
    @test reg isa SymmetryRegistry
    @test length(reg) == 3
    # order Y, XY, X with ids 1, 2, 3
    @test reg[1] isa YAxisReflection && reg[2] isa XYAxisReflection && reg[3] isa XAxisReflection
    @test [s.sym_id for s in reg] == [1, 2, 3]
    # the registered elements act like their untagged counterparts
    pt = SVector(1.0, 2.0)
    @test apply_symmetry(reg[1], pt) == SVector(-1.0, 2.0)
    @test apply_symmetry(reg[2], pt) == SVector(-1.0, -2.0)
    @test apply_symmetry(reg[3], pt) == SVector(1.0, -2.0)
end

@testset "NFoldRotation" begin
    rot = NFoldRotation(4, 1)
    @test rot isa NFoldRotation{Float64}
    @test rot isa AbsSymmetry
    @test rot.order == 4
    @test rot.m == 1
    @test rot.sym_id == 0
    @test rot.angle == 2 * pi / 4
    @test NFoldRotation(4, 1, 6).sym_id == 6

    # the power m is reduced modulo N
    @test NFoldRotation(4, 5).m == 1
    @test NFoldRotation(4, -1).m == 3
    @test NFoldRotation(4, 4).m == 0

    # rotation by 2 pi m / N maps (1, 0) to (cos, sin)
    e1 = SVector(1.0, 0.0)
    @test apply_symmetry(rot, e1) ≈ SVector(0.0, 1.0) atol = sy_ATOL
    @test apply_symmetry(NFoldRotation(3, 1), e1) ≈ SVector(-0.5, sqrt(3) / 2) atol = sy_ATOL
    @test apply_symmetry(NFoldRotation(3, 2), e1) ≈ SVector(-0.5, -sqrt(3) / 2) atol = sy_ATOL
    @test apply_symmetry(NFoldRotation(4, 2), SVector(1.0, 2.0)) ≈ SVector(-1.0, -2.0) atol = sy_ATOL
    @test apply_symmetry(NFoldRotation(4, 0), SVector(1.0, 2.0)) ≈ SVector(1.0, 2.0) atol = sy_ATOL

    # the rotation matrix is orthogonal with determinant +1
    M = rot.sym_map.linear
    @test M' * M ≈ Matrix(I, 2, 2) atol = sy_ATOL
    @test det(M) ≈ 1.0 atol = sy_ATOL

    # N applications of the generator give the identity; rotations compose by adding the powers
    pt = SVector(0.3, -0.8)
    img = pt
    for _ in 1:5
        img = apply_symmetry(NFoldRotation(5, 1), img)
    end
    @test img ≈ pt atol = sy_ATOL
    @test apply_symmetry(NFoldRotation(7, 2), apply_symmetry(NFoldRotation(7, 3), pt)) ≈ apply_symmetry(NFoldRotation(7, 5), pt) atol = sy_ATOL
    # rotations preserve the length
    @test sqrt(sum(abs2, apply_symmetry(NFoldRotation(7, 3), pt))) ≈ sqrt(sum(abs2, pt)) atol = sy_ATOL

    # a collection of points is mapped pointwise
    pts = [SVector(1.0, 0.0), SVector(0.0, 1.0)]
    @test apply_symmetry(rot, pts) ≈ [SVector(0.0, 1.0), SVector(-1.0, 0.0)] atol = sy_ATOL

    # the numeric type of the angle and the matrix follows the keyword T
    rot32 = NFoldRotation(4, 1; T = Float32)
    @test rot32 isa NFoldRotation{Float32}
    @test rot32.angle isa Float32
    @test eltype(rot32.sym_map.linear) == Float32
    @test apply_symmetry(rot32, SVector(1.0f0, 0.0f0)) ≈ SVector(0.0f0, 1.0f0) atol = sy_ATOL_F32
end

@testset "Cn_symmetry" begin
    reg = Cn_symmetry(3)
    @test reg isa SymmetryRegistry
    # n-1 nontrivial images with ids 1..n-1 and powers 1..n-1
    @test length(reg) == 2
    @test all(s -> s isa NFoldRotation{Float64}, reg)
    @test [s.sym_id for s in reg] == [1, 2]
    @test [s.m for s in reg] == [1, 2]
    @test all(s -> s.order == 3, reg)
    @test apply_symmetry(reg[1], SVector(1.0, 0.0)) ≈ SVector(-0.5, sqrt(3) / 2) atol = sy_ATOL

    @test length(Cn_symmetry(6)) == 5
    @test [s.m for s in Cn_symmetry(6)] == 1:5
    @test length(Cn_symmetry(1)) == 0

    # the keyword T selects the numeric type of every rotation
    reg32 = Cn_symmetry(4; T = Float32)
    @test length(reg32) == 3
    @test all(s -> s isa NFoldRotation{Float32}, reg32)
end
