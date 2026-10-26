# Unit tests for src/geometry/symmetryregistry.jl: SymmetryRegistry, register_symmetries, symmetry_of.
# Reference data are small registries of reflections and rotations whose ids follow from the argument order; helpers use the prefix `sr_`.

sr_ids(reg) = [s.sym_id for s in reg]

@testset "SymmetryRegistry" begin
    reg = register_symmetries(XAxisReflection(), YAxisReflection())
    @test reg isa SymmetryRegistry
    @test reg isa AbstractVector{AbsSymmetry}
    @test fieldnames(SymmetryRegistry) == (:generators,)
    @test reg.generators isa Vector{AbsSymmetry}
    @test eltype(reg) === AbsSymmetry

    # wrapping a vector by hand keeps the generators and their ids as given
    hand = SymmetryRegistry(AbsSymmetry[XAxisReflection(4), YAxisReflection(9)])
    @test sr_ids(hand) == [4, 9]
    @test hand.generators[2] == YAxisReflection(9)

    # equal contents compare equal as vectors, distinct registries are distinct objects
    @test register_symmetries(XAxisReflection()) == register_symmetries(XAxisReflection())
    @test register_symmetries(XAxisReflection()) != register_symmetries(YAxisReflection())
    @test register_symmetries(XAxisReflection()) !== register_symmetries(XAxisReflection())
end

@testset "SymmetryRegistry (vector interface)" begin
    reg = register_symmetries(XAxisReflection(), NFoldRotation(3, 1), DiagonalReflection())
    @test length(reg) == 3
    @test size(reg) == (3,)
    @test axes(reg) == (Base.OneTo(3),)
    @test firstindex(reg) == 1
    @test lastindex(reg) == 3
    @test IndexStyle(typeof(reg)) == IndexLinear()

    # indexing returns the stored generators, iteration follows the declaration order
    @test reg[1] isa XAxisReflection
    @test reg[2] isa NFoldRotation
    @test reg[3] isa DiagonalReflection
    @test [typeof(s) for s in reg] == [XAxisReflection, NFoldRotation{Float64}, DiagonalReflection]
    @test reg[2] === reg.generators[2]

    # collections derived from the registry are plain vectors
    @test collect(reg) isa Vector{AbsSymmetry}
    @test reg[2:3] isa Vector{AbsSymmetry}
    @test sr_ids(reg[2:3]) == [2, 3]
    @test map(s -> s.sym_id, reg) == [1, 2, 3]

    # out-of-range indices raise a BoundsError; the registry is read only
    @test_throws BoundsError reg[0]
    @test_throws BoundsError reg[4]
    @test_throws Base.CanonicalIndexError (reg[1] = XAxisReflection(5))

    # the empty registry is a valid zero-length vector
    empty_reg = SymmetryRegistry(AbsSymmetry[])
    @test isempty(empty_reg)
    @test size(empty_reg) == (0,)
end

@testset "register_symmetries" begin
    # ids are assigned 1, 2, ... in argument order and replace whatever id the argument carried
    reg = register_symmetries(YAxisReflection(), XYAxisReflection(), XAxisReflection())
    @test sr_ids(reg) == [1, 2, 3]
    @test reg[1] == YAxisReflection(1)
    @test reg[2] == XYAxisReflection(2)
    @test reg[3] == XAxisReflection(3)
    @test sr_ids(register_symmetries(DiagonalReflection(7), AntiDiagonalReflection(7))) == [1, 2]
    @test register_symmetries(AntiDiagonalReflection(), DiagonalReflection()) == AbsSymmetry[AntiDiagonalReflection(1), DiagonalReflection(2)]

    # the same generator passed twice gets two distinct ids
    @test sr_ids(register_symmetries(XAxisReflection(), XAxisReflection())) == [1, 2]

    # the arguments are not modified (the generators are immutable, so the tagged copies are new values)
    original = DiagonalReflection(7)
    register_symmetries(original)
    @test original.sym_id == 7

    # no generators give an empty registry
    @test register_symmetries() isa SymmetryRegistry
    @test isempty(register_symmetries())
    @test length(register_symmetries()) == 0

    # rotations keep their order, power, angle and numeric type, and only the id changes
    rot = register_symmetries(XAxisReflection(), NFoldRotation(5, 2))[2]
    @test rot isa NFoldRotation{Float64}
    @test (rot.order, rot.m, rot.sym_id) == (5, 2, 2)
    @test rot.angle == NFoldRotation(5, 2).angle
    @test rot.sym_map == NFoldRotation(5, 2).sym_map
    rot32 = register_symmetries(NFoldRotation(4, 1; T = Float32))[1]
    @test rot32 isa NFoldRotation{Float32}
    @test rot32.sym_id == 1
    @test rot32.angle isa Float32

    # the registry type is known from the argument types
    @test (@inferred register_symmetries(XAxisReflection(), YAxisReflection())) isa SymmetryRegistry

    # CompositeReflection is a generator set for orbit bookkeeping, not a taggable generator
    @test_throws MethodError register_symmetries(CompositeReflection(XAxisReflection(), YAxisReflection()))
end

@testset "symmetry_of" begin
    reg = register_symmetries(YAxisReflection(), XYAxisReflection(), XAxisReflection())
    @test symmetry_of(reg, 1) isa YAxisReflection
    @test symmetry_of(reg, 2) isa XYAxisReflection
    @test symmetry_of(reg, 3) isa XAxisReflection
    @test symmetry_of(reg, 2) === reg[2]
    # every generator is found under its own id
    @test all(s -> symmetry_of(reg, s.sym_id) === s, reg)

    # the lookup is by id and not by storage position
    swapped = SymmetryRegistry(AbsSymmetry[XAxisReflection(3), YAxisReflection(1)])
    @test symmetry_of(swapped, 1) isa YAxisReflection
    @test symmetry_of(swapped, 3) isa XAxisReflection
    @test symmetry_of(swapped, 1) === swapped[2]

    # the first generator wins when ids repeat
    dup = SymmetryRegistry(AbsSymmetry[XAxisReflection(4), YAxisReflection(4)])
    @test symmetry_of(dup, 4) isa XAxisReflection

    # an unregistered id, including the default id 0 of an untagged generator, raises an ArgumentError naming the id
    @test_throws ArgumentError symmetry_of(reg, 99)
    @test_throws "sym_id=99" symmetry_of(reg, 99)
    @test_throws ArgumentError symmetry_of(reg, 0)
    @test_throws ArgumentError symmetry_of(register_symmetries(), 1)

    # ids are integers
    @test_throws MethodError symmetry_of(reg, 1.0)
end
