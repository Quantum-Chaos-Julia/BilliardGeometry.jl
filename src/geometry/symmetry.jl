
ident = IdentityTransformation()
reflect_x = LinearMap(SMatrix{2,2}([-1.0 0.0;0.0 1.0]))
reflect_y = LinearMap(SMatrix{2,2}([1.0 0.0;0.0 -1.0]))
reflect_diag = LinearMap(SMatrix{2,2}([0.0 1.0;1.0 0.0]))
reflect_antidiag = LinearMap(SMatrix{2,2}([0.0 -1.0;-1.0 0.0]))
reflect_xy = reflect_x ∘ reflect_y

"""
    XAxisReflection(sym_id=0)

Reflection across the x-axis, `(x,y) -> (x,-y)`.

## Arguments
* `sym_id`: Stable identifier assigned by [`register_symmetries`](@ref) (purely geometric; carries no representation/parity data — see [`SymmetrySector`](@ref) in `QuantumBilliards.jl` for the per-solve representation choice).
"""
struct XAxisReflection <: AbsReflection
    sym_id::Int
end
XAxisReflection() = XAxisReflection(0)

"""
    YAxisReflection(sym_id=0)

Reflection across the y-axis, `(x,y) -> (-x,y)`.

## Arguments
* `sym_id`: Stable identifier assigned by [`register_symmetries`](@ref).
"""
struct YAxisReflection <: AbsReflection
    sym_id::Int
end
YAxisReflection() = YAxisReflection(0)

"""
    XYAxisReflection(sym_id=0)

Reflection across both coordinate axes, `(x,y) -> (-x,-y)`.

## Arguments
* `sym_id`: Stable identifier assigned by [`register_symmetries`](@ref).
"""
struct XYAxisReflection <: AbsReflection 
    sym_id::Int
end
XYAxisReflection() = XYAxisReflection(0)

"""
    DiagonalReflection(sym_id=0)

Reflection across the `y=x` diagonal, `(x,y) -> (y,x)`.

## Arguments
* `sym_id`: Stable identifier assigned by [`register_symmetries`](@ref).
"""
struct DiagonalReflection <: AbsReflection
    sym_id::Int
end
DiagonalReflection() = DiagonalReflection(0)

"""
    AntiDiagonalReflection(sym_id=0)

Reflection across the `y=-x` anti-diagonal, `(x,y) -> (-y,-x)`.

## Arguments
* `sym_id`: Stable identifier assigned by [`register_symmetries`](@ref).
"""
struct AntiDiagonalReflection <: AbsReflection
    sym_id::Int
end
AntiDiagonalReflection() = AntiDiagonalReflection(0)

# Layer-1 (geometric) reconstruction helpers used only by `register_symmetries`
# to tag a fresh copy of a generator with its assigned `sym_id`.
_with_sym_id(::XAxisReflection, id::Int) = XAxisReflection(id)
_with_sym_id(::YAxisReflection, id::Int) = YAxisReflection(id)
_with_sym_id(::XYAxisReflection, id::Int) = XYAxisReflection(id)
_with_sym_id(::DiagonalReflection, id::Int) = DiagonalReflection(id)
_with_sym_id(::AntiDiagonalReflection, id::Int) = AntiDiagonalReflection(id)

"""
    CompositeReflection(reflections)
    CompositeReflection(reflections...)

Combines several reflection constraints into one symmetry-group generator
set. Used only by [`symmetry_index_orbits`](@ref), which generates the
complete symmetry group by closure under composition; `CompositeReflection`
is not itself a single applyable geometric transformation (there is no
`apply_symmetry` method for it).

## Attributes
* `reflections`: Reflection constraints used as generators of the composite symmetry group.
"""
struct CompositeReflection <: AbsReflection
    reflections::Vector{AbsReflection}
end
CompositeReflection(reflections::AbstractVector{<:AbsReflection}) = CompositeReflection(AbsReflection[reflections...])
CompositeReflection(reflections::AbsReflection...) = CompositeReflection(AbsReflection[reflections...])

function apply_symmetry(sym::XAxisReflection, pt::SVector{2,T}) where T<:Real
    return reflect_y(pt)
end

function apply_symmetry(sym::XAxisReflection, pts)
    return [reflect_y(pt) for pt in pts]
end

function apply_symmetry(sym::YAxisReflection, pt::SVector{2,T})  where T<:Real
    return reflect_x(pt)
end
function apply_symmetry(sym::YAxisReflection, pts)
    return [reflect_x(pt) for pt in pts]
end

function apply_symmetry(sym::XYAxisReflection, pt::SVector{2,T})  where T<:Real
    return reflect_xy(pt)
end
function apply_symmetry(sym::XYAxisReflection, pts)
    return [reflect_xy(pt) for pt in pts]
end

function apply_symmetry(sym::DiagonalReflection, pt::SVector{2,T}) where T<:Real
    return reflect_diag(pt)
end
function apply_symmetry(sym::DiagonalReflection, pts)
    return [reflect_diag(pt) for pt in pts]
end

function apply_symmetry(sym::AntiDiagonalReflection, pt::SVector{2,T}) where T<:Real
    return reflect_antidiag(pt)
end
function apply_symmetry(sym::AntiDiagonalReflection, pts)
    return [reflect_antidiag(pt) for pt in pts]
end

function apply_symmetry_pb(sym::AbsReflection, sym_sector::Int64, s::T, p::T, L::T) where T<:Real
    if sym_sector == 1
        return s, p
    elseif sym_sector == 2
        return 2*L-s, -p
    elseif sym_sector == 3
        return 2*L+s, p
    elseif sym_sector == 4    
        return 4*L-s, -p
    end
end

"""
    D2_symmetry() → reg::SymmetryRegistry

Registers the three non-identity elements of the `D₂` reflection group
(`YAxisReflection`, `XYAxisReflection`, `XAxisReflection`, in that order) via
[`register_symmetries`](@ref), giving them fresh `sym_id`s `1,2,3`.
"""
D2_symmetry() = register_symmetries(YAxisReflection(), XYAxisReflection(), XAxisReflection())

abstract type AbsRotation <: AbsSymmetry end

# Counter-clockwise 2x2 rotation matrix, used instead of `LinearMap(RotZ(θ))`
# (Rotations.jl's `RotZ` is a 3x3 rotation and cannot convert to the
# `LinearMap{SMatrix{2,2,Float64,4}}` field type below).
@inline function rotation_matrix_z(θ::T) where T<:Real
    s, c = sincos(θ)
    return SMatrix{2,2,T}(c, s, -s, c)
end

"""
    NFoldRotation(N,m,sym_id=0)

One nontrivial image of an `N`-fold rotational symmetry: counter-clockwise
rotation by `2π*m/N`.

## Arguments
* `N`: Order of the rotational symmetry.
* `m`: Power of the fundamental rotation.
* `sym_id`: Stable identifier assigned by [`register_symmetries`](@ref) (purely geometric; carries no representation data \u2014 see `SymmetrySector` in `QuantumBilliards.jl` for the per-solve irrep-sector choice).
"""
struct NFoldRotation <: AbsRotation
    order::Int64
    m::Int64
    sym_id::Int64
    angle::Float64
    sym_map::LinearMap{SMatrix{2, 2, Float64, 4}}
end

function NFoldRotation(N, m, sym_id::Int=0)
    angle = 2*pi/N
    mm = mod(m, N)
    sym_map = LinearMap(rotation_matrix_z(angle*mm))
    return NFoldRotation(N, mm, sym_id, angle, sym_map)
end

_with_sym_id(sym::NFoldRotation, id::Int) = NFoldRotation(sym.order, sym.m, id)

function apply_symmetry(sym::NFoldRotation, pt::SVector{2,T}) where T<:Real
    return sym.sym_map(pt)
end
function apply_symmetry(sym::NFoldRotation, pts)
    return [sym.sym_map(pt) for pt in pts]
end

"""
    Cn_symmetry(n) → reg::SymmetryRegistry

Registers the `n-1` non-identity images of an `n`-fold rotational symmetry
(`NFoldRotation(n,1), ..., NFoldRotation(n,n-1)`) via
[`register_symmetries`](@ref), giving them fresh `sym_id`s `1,...,n-1`.
"""
Cn_symmetry(n) = register_symmetries((NFoldRotation(n,i) for i in 1:(n-1))...)