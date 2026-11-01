# Unit tests for the discrete symmetry generators (reflections and rotations) of a billiard:
#   src/geometry/symmetry.jl      (AxisReflection and its aliases, show, symmetry_key, CompositeReflection ==/hash, _sym_matrix,
#                                  apply_symmetry, apply_symmetry_pb, D2_symmetry, D4_symmetry, NFoldRotation, Cn_symmetry,
#                                  _with_sym_id)
#   src/geometry/fullboundary.jl  (orientation_reversing)
# Reference objects: the five coordinate reflections acting on the point (x, y), namely X: (x, -y), Y: (-x, y),
# XY: (-x, -y), Diagonal: (y, x) and AntiDiagonal: (-y, -x); and the counterclockwise rotation by 2 pi m / N, which maps
# (1, 0) to (cos(2 pi m / N), sin(2 pi m / N)). Pure reflections have determinant -1 (they reverse the boundary
# orientation), XY and rotations have determinant +1. The registry itself (register_symmetries, symmetry_of) is not
# tested here. AxisReflection{A} is the mirror line through the origin at angle pi A (A in half-turns), with matrix
# Rot(2 pi A) * diag(1, -1); the four special axes A = 0, 1/4, 1/2, 3/4 are the const aliases X, Diagonal, Y, AntiDiagonal.
# The dihedral key (reflects, a) of an element is the exact label of its action (rotation by 2 pi a or mirror at angle pi a).
# No eigen solves, no pinned wavefunction values. Helper constants use the prefix `sy_`.

using StaticArrays: SVector, SMatrix
using LinearAlgebra: det, I, norm
using BilliardGeometry: AbsReflection, _sym_matrix, _with_sym_id, symmetry_key, _reconstruct_curves, _keys_closure, _keys_closure

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

@testset "AxisReflection aliases" begin
    @test XAxisReflection === AxisReflection{0 // 1}
    @test DiagonalReflection === AxisReflection{1 // 4}
    @test YAxisReflection === AxisReflection{1 // 2}
    @test AntiDiagonalReflection === AxisReflection{3 // 4}
    # XY is a separate struct: a rotation by pi, not a mirror
    @test XYAxisReflection !== AxisReflection{1 // 2}
    @test !(XYAxisReflection() isa AxisReflection)

    # constructors: sym_id defaults to 0 for every axis, also a generic one
    @test AxisReflection{1 // 6}().sym_id == 0
    @test AxisReflection{1 // 6}(4).sym_id == 4
    @test AxisReflection{1 // 6}() isa AbsReflection
    @test XAxisReflection(2) isa AxisReflection{0 // 1}
    @test typeof(YAxisReflection()) === AxisReflection{1 // 2}
    @test DiagonalReflection() isa AxisReflection
    # the mirror angle is part of the type
    @test AxisReflection{1 // 6}() isa AxisReflection{1 // 6}
    @test !(AxisReflection{1 // 6}() isa AxisReflection{1 // 3})

    # _with_sym_id keeps the generic axis
    tagged = _with_sym_id(AxisReflection{1 // 6}(), 9)
    @test typeof(tagged) === AxisReflection{1 // 6}
    @test tagged.sym_id == 9
end

@testset "AxisReflection show" begin
    # the four special axes print their alias names, with the sym_id
    @test sprint(show, XAxisReflection()) == "XAxisReflection(0)"
    @test sprint(show, YAxisReflection(3)) == "YAxisReflection(3)"
    @test sprint(show, DiagonalReflection(2)) == "DiagonalReflection(2)"
    @test sprint(show, AntiDiagonalReflection(5)) == "AntiDiagonalReflection(5)"
    # any other axis prints the parametric constructor
    @test sprint(show, AxisReflection{1 // 6}(7)) == "AxisReflection{1//6}(7)"
    @test sprint(show, AxisReflection{1 // 6}()) == "AxisReflection{1//6}(0)"
end

@testset "symmetry_key" begin
    @test symmetry_key(XAxisReflection()) == (true, 0 // 1)
    @test symmetry_key(DiagonalReflection()) == (true, 1 // 4)
    @test symmetry_key(YAxisReflection()) == (true, 1 // 2)
    @test symmetry_key(AntiDiagonalReflection()) == (true, 3 // 4)
    @test symmetry_key(XYAxisReflection()) == (false, 1 // 2)
    @test symmetry_key(NFoldRotation(4, 1)) == (false, 1 // 4)
    @test symmetry_key(NFoldRotation(6, 2)) == (false, 1 // 3)         # reduced fraction of a full turn
    @test symmetry_key(NFoldRotation(5, 0)) == (false, 0 // 1)
    # the angle is taken modulo 1: a mirror at angle pi (A = 1) is the x axis again
    @test symmetry_key(AxisReflection{1 // 1}()) == (true, 0 // 1)
    @test symmetry_key(AxisReflection{5 // 4}()) == (true, 1 // 4)
    @test symmetry_key(AxisReflection{-1 // 4}()) == (true, 3 // 4)
    # the sym_id does not enter the key
    @test symmetry_key(XAxisReflection(5)) == symmetry_key(XAxisReflection())
    # D4 registry: keys in registration (chain) order Diag, R1, Y, XY, Anti, R3, X
    @test [symmetry_key(s) for s in D4_symmetry()] ==
          [(true, 1 // 4), (false, 1 // 4), (true, 1 // 2), (false, 1 // 2), (true, 3 // 4), (false, 3 // 4), (true, 0 // 1)]
end

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

@testset "_sym_matrix (generic AxisReflection)" begin
    # Rot(2 pi A) * diag(1, -1) = [cos(2 pi A) sin(2 pi A); sin(2 pi A) -cos(2 pi A)]
    M = _sym_matrix(AxisReflection{1 // 6}())
    @test M isa SMatrix{2,2,Float64,4}
    @test M ≈ [0.5 sqrt(3) / 2; sqrt(3) / 2 -0.5] atol = sy_ATOL
    @test _sym_matrix(AxisReflection{1 // 6}(), Float32) isa SMatrix{2,2,Float32,4}
    @test _sym_matrix(AxisReflection{1 // 6}(), BigFloat) isa SMatrix{2,2,BigFloat,4}

    # the generic method agrees with the exact special matrices: A = 5/4 is the diagonal in a different representation
    @test _sym_matrix(AxisReflection{5 // 4}()) ≈ _sym_matrix(DiagonalReflection()) atol = sy_ATOL
    @test _sym_matrix(AxisReflection{-1 // 2}()) ≈ _sym_matrix(YAxisReflection()) atol = sy_ATOL

    # every axis is an involution with determinant -1 that fixes the mirror line (cos(pi A), sin(pi A))
    for A in (1 // 6, 1 // 3, 1 // 8, 5 // 12, 7 // 8)
        sym = AxisReflection{A}()
        Ms = _sym_matrix(sym)
        @test Ms * Ms ≈ Matrix(I, 2, 2) atol = sy_ATOL
        @test det(Ms) ≈ -1.0 atol = sy_ATOL
        axis = SVector(cospi(Float64(A)), sinpi(Float64(A)))
        @test Ms * axis ≈ axis atol = sy_ATOL
    end

    # special axes keep their exact entries
    @test _sym_matrix(DiagonalReflection()) == [0.0 1.0; 1.0 0.0]
    @test _sym_matrix(XAxisReflection()) == [1.0 0.0; 0.0 -1.0]
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

@testset "apply_symmetry (generic AxisReflection)" begin
    # mirror at angle pi/6 maps (1, 0) to (cos(pi/3), sin(pi/3))
    img = apply_symmetry(AxisReflection{1 // 6}(), SVector(1.0, 0.0))
    @test img ≈ SVector(0.5, sqrt(3) / 2) atol = sy_ATOL
    # a point on the mirror line is fixed, and the map is an involution
    on_axis = SVector(cospi(1 / 6), sinpi(1 / 6))
    @test apply_symmetry(AxisReflection{1 // 6}(), on_axis) ≈ on_axis atol = sy_ATOL
    pt = SVector(0.3, -0.8)
    @test apply_symmetry(AxisReflection{1 // 6}(), apply_symmetry(AxisReflection{1 // 6}(), pt)) ≈ pt atol = sy_ATOL
    # two mirrors compose to a rotation by twice the angle between them: X after the pi/6 mirror is rotation by -pi/3
    r = apply_symmetry(XAxisReflection(), apply_symmetry(AxisReflection{1 // 6}(), pt))
    @test r ≈ apply_symmetry(NFoldRotation(6, 5), pt) atol = sy_ATOL
    # a collection of points is mapped pointwise and the point type is preserved
    @test apply_symmetry(AxisReflection{1 // 6}(), [SVector(1.0, 0.0)]) ≈ [SVector(0.5, sqrt(3) / 2)] atol = sy_ATOL
    @test apply_symmetry(AxisReflection{1 // 6}(), SVector(1.0f0, 0.0f0)) isa SVector{2,Float32}
    @test_throws ArgumentError apply_symmetry(AxisReflection{1 // 6}(), [1.0, 0.0])

    # type stability of the special and the generic method
    @test (@inferred apply_symmetry(DiagonalReflection(), SVector(1.0, 2.0))) == SVector(2.0, 1.0)
    @test (@inferred apply_symmetry(XAxisReflection(), SVector(1.0, 2.0))) isa SVector{2,Float64}
    @test (@inferred apply_symmetry(AxisReflection{1 // 6}(), SVector(1.0, 2.0))) isa SVector{2,Float64}
end

@testset "apply_symmetry (plain vector point)" begin
    # G15: points are SVector{2}; a plain Vector of reals is rejected with an ArgumentError instead of recursing forever
    for sym in (XAxisReflection(), YAxisReflection(), XYAxisReflection(), DiagonalReflection(), AntiDiagonalReflection(),
                NFoldRotation(3, 1), NFoldRotation(4, 2))
        @test_throws ArgumentError apply_symmetry(sym, [1.0, 2.0])
        @test_throws ArgumentError apply_symmetry(sym, Float32[1, 2])
    end
    err = try
        apply_symmetry(XAxisReflection(), [1.0, 2.0])
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("SVector", err.msg)

    # collections of SVector points and single SVector points are unaffected
    @test apply_symmetry(XAxisReflection(), [SVector(1.0, 2.0)]) == [SVector(1.0, -2.0)]
    @test apply_symmetry(NFoldRotation(4, 1), [SVector(1.0, 0.0)]) ≈ [SVector(0.0, 1.0)] atol = 1e-12
end

@testset "orientation_reversing" begin
    @test orientation_reversing(XAxisReflection())
    @test orientation_reversing(YAxisReflection())
    @test orientation_reversing(DiagonalReflection())
    @test orientation_reversing(AntiDiagonalReflection())
    @test !orientation_reversing(XYAxisReflection())
    @test !orientation_reversing(NFoldRotation(4, 1))
    # every AxisReflection reverses the orientation, also a generic axis
    @test orientation_reversing(AxisReflection{1 // 6}())
    @test orientation_reversing(AxisReflection{5 // 12}(3))
    # agrees with the sign of the determinant of the linear map
    for (sym, _) in sy_REFLECTIONS
        @test orientation_reversing(sym) == (det(_sym_matrix(sym)) < 0)
    end
    @test orientation_reversing(AxisReflection{1 // 6}()) == (det(_sym_matrix(AxisReflection{1 // 6}())) < 0)
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

@testset "D4_symmetry" begin
    reg = D4_symmetry()
    @test reg isa SymmetryRegistry
    @test length(reg) == 7
    # chain order Diagonal, rotation by pi/2, Y, XY, AntiDiagonal, rotation by 3 pi/2, X with ids 1..7
    @test [typeof(s) for s in reg][[1, 3, 4, 5, 7]] ==
          [DiagonalReflection, YAxisReflection, XYAxisReflection, AntiDiagonalReflection, XAxisReflection]
    @test all(s -> s isa NFoldRotation, [reg[2], reg[6]])
    @test (reg[2].order, reg[2].m) == (4, 1)
    @test (reg[6].order, reg[6].m) == (4, 3)
    @test [s.sym_id for s in reg] == collect(1:7)
    # the seven elements and the identity are eight distinct geometric actions
    @test length(unique(vcat([(false, 0 // 1)], [symmetry_key(s) for s in reg]))) == 8
    # the registered elements act like their untagged counterparts
    pt = SVector(1.0, 2.0)
    @test apply_symmetry(reg[1], pt) == SVector(2.0, 1.0)
    @test apply_symmetry(reg[5], pt) == SVector(-2.0, -1.0)
    @test apply_symmetry(reg[2], pt) ≈ SVector(-2.0, 1.0) atol = sy_ATOL
    @test apply_symmetry(reg[6], pt) ≈ SVector(2.0, -1.0) atol = sy_ATOL
end

# Largest distance between the end of curve j and the start of curve j + 1 (cyclically) of a loop of curves.
sy_max_chain_gap(crvs) = maximum(norm(curve(crvs[j], 1.0) - curve(crvs[mod1(j + 1, length(crvs))], 0.0)) for j in eachindex(crvs))

@testset "D4_symmetry (chain order)" begin
    # one piece spanning the angles 0..pi/4 counter-clockwise; block 1 is the piece, block j + 1 the image under the j-th element
    piece = CircleSegment(1.0, pi / 4, 0.0, SVector(0.0, 0.0))
    keys = vcat([(false, 0 // 1)], [symmetry_key(s) for s in D4_symmetry()])
    @test length(unique(keys)) == 8
    curves, _, sector, _ = _reconstruct_curves([piece], D4_symmetry())
    @test length(curves) == 8
    @test sector == collect(1:8)
    # in registration order block j + 1 is entered where block j ends, so the eight images form a closed loop
    @test sy_max_chain_gap(curves) < 1e-12
    # the previous order (Y, XY, X, Diagonal, AntiDiagonal, R1, R3) does not chain
    old_order = register_symmetries(YAxisReflection(), XYAxisReflection(), XAxisReflection(), DiagonalReflection(),
                                    AntiDiagonalReflection(), NFoldRotation(4, 1), NFoldRotation(4, 3))
    old_curves, = _reconstruct_curves([piece], old_order)
    @test sy_max_chain_gap(old_curves) > 1
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

# Dn_symmetry: image k = 1..2n-1 of the wedge starting at angle pi*shift is the mirror at axis shift + (k + 1) / (2 n) for odd k
# and the rotation by 2 pi (k / 2) / n for even k (XYAxisReflection when that is the pi rotation). The registry is closed under
# composition, so the dihedral keys of its elements plus the identity form a group of order 2n generated by two mirrors. The
# D3 shift = 1/6 listing is an observed run.
@testset "Dn_symmetry" begin
    @testset "equals D2_symmetry and D4_symmetry" begin
        @test collect(Dn_symmetry(2)) == collect(D2_symmetry())
        @test string.(Dn_symmetry(4)) == string.(D4_symmetry())
        @test [s.sym_id for s in Dn_symmetry(4)] == [s.sym_id for s in D4_symmetry()]
    end

    @testset "element types and ids" begin
        reg = Dn_symmetry(3)
        @test reg isa SymmetryRegistry
        @test length(reg) == 5
        @test [s.sym_id for s in reg] == 1:5
        @test [symmetry_key(s) for s in reg] == [(true, 1 // 3), (false, 1 // 3), (true, 2 // 3), (false, 2 // 3), (true, 0 // 1)]
        @test all(s -> s isa NFoldRotation{Float64}, reg[[2, 4]])
        @test length(Dn_symmetry(6)) == 11
        @test collect(Dn_symmetry(1)) == [XAxisReflection(1)]
    end

    @testset "shift" begin
        reg = Dn_symmetry(3; shift = 1 // 6)
        @test reg[1] == YAxisReflection(1)
        @test reg[2] isa NFoldRotation && reg[2].order == 3 && reg[2].m == 1 && reg[2].sym_id == 2
        @test reg[3] == AxisReflection{5 // 6}(3)
        @test reg[4] isa NFoldRotation && reg[4].order == 3 && reg[4].m == 2 && reg[4].sym_id == 4
        @test reg[5] == AxisReflection{1 // 6}(5)
    end

    @testset "numeric type" begin
        reg = Dn_symmetry(3; T = Float32)
        @test all(s -> s isa NFoldRotation{Float32}, reg[[2, 4]])
    end

    @testset "invalid n" begin
        @test_throws ArgumentError Dn_symmetry(0)
        @test_throws ArgumentError Dn_symmetry(-2)
    end

    @testset "closure under composition and generators" begin
        for n in (3, 5, 6)
            reg = Dn_symmetry(n)
            keys = [symmetry_key(s) for s in reg]
            @test length(unique(keys)) == 2n - 1
            @test length(_keys_closure(keys)) == 2n
            gens = get_symmetries(reg)
            @test length(gens) == 2
            @test all(g -> g isa AbsReflection, gens)
            @test length(_keys_closure([symmetry_key(g) for g in gens])) == 2n
        end
    end
end

# Dn_symmetry (src/geometry/symmetry.jl). Reference: image k = 1..2n-1 of the wedge starting at angle pi*shift is the mirror at
# axis shift + (k + 1) / (2 n) for odd k and the rotation by 2 pi (k / 2) / n for even k (XYAxisReflection when that is the pi
# rotation); the registry is closed under composition, so the dihedral keys of its elements plus the identity form a group of
# order 2n, generated by two mirrors. The D3 shift = 1/6 listing is an observed run.
@testset "Dn_symmetry" begin
    @testset "equals D2_symmetry and D4_symmetry" begin
        @test collect(Dn_symmetry(2)) == collect(D2_symmetry())
        @test string.(Dn_symmetry(4)) == string.(D4_symmetry())
        @test [s.sym_id for s in Dn_symmetry(4)] == [s.sym_id for s in D4_symmetry()]
    end

    @testset "element types and ids" begin
        reg = Dn_symmetry(3)
        @test reg isa SymmetryRegistry
        @test length(reg) == 5
        @test [s.sym_id for s in reg] == 1:5
        @test [symmetry_key(s) for s in reg] == [(true, 1 // 3), (false, 1 // 3), (true, 2 // 3), (false, 2 // 3), (true, 0 // 1)]
        @test all(s -> s isa NFoldRotation{Float64}, reg[[2, 4]])
        @test length(Dn_symmetry(1)) == 1
        @test collect(Dn_symmetry(1)) == [XAxisReflection(1)]
        @test length(Dn_symmetry(6)) == 11
    end

    @testset "shift" begin
        reg = Dn_symmetry(3; shift = 1 // 6)
        @test reg[1] == YAxisReflection(1)
        @test reg[2] isa NFoldRotation && reg[2].order == 3 && reg[2].m == 1 && reg[2].sym_id == 2
        @test reg[3] == AxisReflection{5 // 6}(3)
        @test reg[4] isa NFoldRotation && reg[4].order == 3 && reg[4].m == 2 && reg[4].sym_id == 4
        @test reg[5] == AxisReflection{1 // 6}(5)
    end

    @testset "numeric type" begin
        reg = Dn_symmetry(3; T = Float32)
        @test all(s -> s isa NFoldRotation{Float32}, reg[[2, 4]])
    end

    @testset "invalid n" begin
        @test_throws ArgumentError Dn_symmetry(0)
        @test_throws ArgumentError Dn_symmetry(-2)
    end

    @testset "closure under composition and generators" begin
        for n in (3, 5, 6)
            reg = Dn_symmetry(n)
            keys = [symmetry_key(s) for s in reg]
            @test length(unique(keys)) == 2n - 1
            @test length(_keys_closure(keys)) == 2n
            gens = get_symmetries(reg)
            @test length(gens) == 2
            @test all(g -> g isa AbsReflection, gens)
            @test length(_keys_closure([symmetry_key(g) for g in gens])) == 2n
        end
    end
end
