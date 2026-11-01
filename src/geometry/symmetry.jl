"""
    AxisReflection{A}(sym_id=0)

Reflection across the line through the origin at angle `π*A`, with `A::Rational{Int}`
the mirror angle in half-turns: matrix `Rot(2π*A) * diag(1,-1)`. The four special axes
have aliases [`XAxisReflection`](@ref) (`A=0`), [`DiagonalReflection`](@ref) (`A=1/4`),
[`YAxisReflection`](@ref) (`A=1/2`) and [`AntiDiagonalReflection`](@ref) (`A=3/4`).

## Arguments
* `sym_id`: Stable identifier assigned by [`register_symmetries`](@ref) (purely geometric; carries no representation/parity data — see [`SymmetrySector`](@ref) in `QuantumBilliards.jl` for the per-solve representation choice).
"""
struct AxisReflection{A} <: AbsReflection
    sym_id::Int
end
AxisReflection{A}() where {A} = AxisReflection{A}(0)

"""
    XAxisReflection(sym_id=0)

Reflection across the x-axis, `(x,y) -> (x,-y)`. Alias of `AxisReflection{0//1}`.
"""
const XAxisReflection = AxisReflection{0//1}

"""
    YAxisReflection(sym_id=0)

Reflection across the y-axis, `(x,y) -> (-x,y)`. Alias of `AxisReflection{1//2}`.
"""
const YAxisReflection = AxisReflection{1//2}

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

Reflection across the `y=x` diagonal, `(x,y) -> (y,x)`. Alias of `AxisReflection{1//4}`.
"""
const DiagonalReflection = AxisReflection{1//4}

"""
    AntiDiagonalReflection(sym_id=0)

Reflection across the `y=-x` anti-diagonal, `(x,y) -> (-y,-x)`. Alias of `AxisReflection{3//4}`.
"""
const AntiDiagonalReflection = AxisReflection{3//4}

function Base.show(io::IO, sym::AxisReflection{A}) where {A}
    name = A == 0 ? "XAxisReflection" : A == 1//2 ? "YAxisReflection" :
           A == 1//4 ? "DiagonalReflection" : A == 3//4 ? "AntiDiagonalReflection" :
           "AxisReflection{$A}"
    print(io, name, "(", sym.sym_id, ")")
end

# Reconstruction helpers used only by `register_symmetries`
# to tag a fresh copy of a generator with its assigned `sym_id`.
_with_sym_id(::AxisReflection{A}, id::Int) where {A} = AxisReflection{A}(id)
_with_sym_id(::XYAxisReflection, id::Int) = XYAxisReflection(id)
_with_sym_id(sym::AbsSymmetry, ::Int) = throw(ArgumentError("$(typeof(sym)) cannot be registered as a billiard symmetry element"))

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

Base.:(==)(a::CompositeReflection, b::CompositeReflection) = a.reflections == b.reflections
Base.hash(a::CompositeReflection, h::UInt) = hash(a.reflections, hash(:CompositeReflection, h))

"""
    _sym_matrix(sym::AbsReflection, ::Type{T}=Float64) where T<:Real → M::SMatrix{2,2,T}

2x2 linear-map matrix implementing [`apply_symmetry`](@ref) for a pure
reflection generator, built directly at the numeric type `T` of the point
being transformed (defaults to `Float64` for callers with no natural `T`,
e.g. symmetry-orbit bookkeeping that only inspects the matrix's sign
pattern). Entries are exact `±1`/`0` integer literals, so promotion to any
`T<:Real` (including `BigFloat`) is lossless. Reused by `apply_symmetry` and
by `fullboundary.jl`'s curve-image reconstruction.
"""
_sym_matrix(::XAxisReflection, ::Type{T}=Float64) where {T<:Real} = SMatrix{2,2,T}(1, 0, 0, -1)
_sym_matrix(::YAxisReflection, ::Type{T}=Float64) where {T<:Real} = SMatrix{2,2,T}(-1, 0, 0, 1)
_sym_matrix(::XYAxisReflection, ::Type{T}=Float64) where {T<:Real} = SMatrix{2,2,T}(-1, 0, 0, -1)
_sym_matrix(::DiagonalReflection, ::Type{T}=Float64) where {T<:Real} = SMatrix{2,2,T}(0, 1, 1, 0)
_sym_matrix(::AntiDiagonalReflection, ::Type{T}=Float64) where {T<:Real} = SMatrix{2,2,T}(0, -1, -1, 0)
function _sym_matrix(::AxisReflection{A}, ::Type{T}=Float64) where {A,T<:Real}
    s, c = sincospi(T(2 * A))
    return SMatrix{2,2,T}(c, s, s, -c)
end

"""
    apply_symmetry(sym::AbsSymmetry, pt::SVector{2,T}) → pt′::SVector{2,T}
    apply_symmetry(sym::AbsSymmetry, pts) → pts′::Vector

Apply the geometric transformation `sym` to a point, or to every point of a
collection.

## Description
Reflections act by the exact `±1`/`0` matrix of `_sym_matrix`,
[`NFoldRotation`](@ref) by its stored rotation. [`CompositeReflection`](@ref)
is a generator set and has no method. Plain vectors of reals are rejected with
an `ArgumentError`: convert them to `SVector{2}` first.

## Arguments
* `sym`: Symmetry element, for example [`XAxisReflection`](@ref) or [`NFoldRotation`](@ref).
* `pt`, `pts`: A planar point, or an iterable of points.

## Returns
* `pt′`, `pts′`: The image point(s), with the same numeric type as the input.
"""
function apply_symmetry(sym::AbsReflection, pt::SVector{2,T}) where T<:Real
    return SVector{2,T}(_sym_matrix(sym, T) * pt)
end
apply_symmetry(sym::AbsReflection, pts) = [apply_symmetry(sym, pt) for pt in pts]

_throw_not_svector_point(sym) = throw(ArgumentError("apply_symmetry(::$(nameof(typeof(sym))), pt) requires a point of type SVector{2}; got a plain vector of reals. Convert with SVector{2}(pt)."))
apply_symmetry(sym::AbsReflection, ::AbstractVector{<:Real}) = _throw_not_svector_point(sym)

"""
    apply_symmetry_pb(sym::AbsSymmetry, sym_sector::Int64, s::T, p::T, L::T) where T<:Real

Maps a Poincaré–Birkhoff (arc-length, sine-of-angle) coordinate `(s,p)`,
measured within the fundamental domain (of length `L`), to its global
coordinate on the `sym_sector`-th copy of the full physical boundary
(`sym_sector==1` is the fundamental domain itself, `sym_sector==k+1`
corresponds to the image under `billiard.symmetries[k]`).

## Description
Mirrors [`full_boundary`](@ref)'s reconstruction: the full physical boundary
is the concatenation of the fundamental domain (`[0,L]`) followed by, for
`k=1,...,length(billiard.symmetries)`, the image of the fundamental domain
under `billiard.symmetries[k]`. Orientation-preserving images
(`XYAxisReflection`, [`NFoldRotation`](@ref)) continue the same traversal
direction, so `s` maps to `k*L+s` with `p` unchanged; orientation-reversing
images (pure reflections) reverse it, so `s` maps to `(k+1)*L-s` with `p`
negated — see [`orientation_reversing`](@ref).

## Arguments
* `sym`: The symmetry generator whose image sector `(s,p)` is mapped into (only its orientation-reversal behavior is used; not called for `sym_sector==1`).
* `sym_sector`: 1-based index of the boundary copy, `k+1` where `k` is `sym`'s 1-based position in `billiard.symmetries`.
* `s`,`p`: Arc length / sine-of-angle coordinates within the fundamental domain.
* `L`: Length of the fundamental domain's physical boundary.
"""
function apply_symmetry_pb(sym::AbsSymmetry, sym_sector::Int64, s::T, p::T, L::T) where T<:Real
    k = sym_sector - 1
    if orientation_reversing(sym)
        return (k+1)*L - s, -p
    else
        return k*L + s, p
    end
    throw(ArgumentError("sym_sector must be in 1:4; received $sym_sector"))
end

"""
    D2_symmetry() → reg::SymmetryRegistry

Registers the three non-identity elements of the `D₂` reflection group
(`YAxisReflection`, `XYAxisReflection`, `XAxisReflection`, in that order) via
[`register_symmetries`](@ref), giving them fresh `sym_id`s `1,2,3`.
"""
D2_symmetry() = register_symmetries(YAxisReflection(), XYAxisReflection(), XAxisReflection())

"""
    D4_symmetry() → reg::SymmetryRegistry

Registers the seven non-identity elements of the `D₄` group (`DiagonalReflection`,
`NFoldRotation(4,1)`, `YAxisReflection`, `XYAxisReflection`,
`AntiDiagonalReflection`, `NFoldRotation(4,3)`, `XAxisReflection`, in that
order) via [`register_symmetries`](@ref), giving them `sym_id`s `1,...,7`.

The order is the reconstruction order of a fundamental piece that runs
counterclockwise from the angle `0` to the angle `π/4` (for example the
fundamental triangle `0 ≤ y ≤ x`): the images, in this order, form one
contiguous closed loop (see `_reconstruct_curves`).
"""
D4_symmetry() = register_symmetries(DiagonalReflection(), NFoldRotation(4, 1), YAxisReflection(), XYAxisReflection(),
    AntiDiagonalReflection(), NFoldRotation(4, 3), XAxisReflection())

"""
    AbsRotation <: AbsSymmetry

Supertype of rotational symmetries. See [`NFoldRotation`](@ref).
"""
abstract type AbsRotation <: AbsSymmetry end

# Counter-clockwise 2x2 rotation matrix, used instead of `LinearMap(RotZ(θ))`
# (Rotations.jl's `RotZ` is a 3x3 rotation and cannot convert to the
# `LinearMap{SMatrix{2,2,Float64,4}}` field type below).
@inline function rotation_matrix_z(θ::T) where T<:Real
    s, c = sincos(θ)
    return SMatrix{2,2,T}(c, s, -s, c)
end

"""
    NFoldRotation(N,m,sym_id=0; T::Type{<:Real}=Float64)

One nontrivial image of an `N`-fold rotational symmetry: counter-clockwise
rotation by `2π*m/N`.

## Arguments
* `N`: Order of the rotational symmetry.
* `m`: Power of the fundamental rotation.
* `sym_id`: Stable identifier assigned by [`register_symmetries`](@ref) (purely geometric; carries no representation data — see `SymmetrySector` in `QuantumBilliards.jl` for the per-solve irrep-sector choice).
* `T`: Numeric type of the rotation angle/matrix, matching the billiard's own `T<:Real` (defaults to `Float64`; pass explicitly for `BigFloat`/`Float32` billiards so the rotation isn't silently narrowed/widened to `Float64`).

## Attributes
* `order::Int64`: Order `N` of the rotational symmetry.
* `m::Int64`: Rotation power, reduced modulo `N`.
* `sym_id::Int64`: Identifier assigned by [`register_symmetries`](@ref).
* `angle::T`: Fundamental rotation angle `2π/N`.
* `sym_map::LinearMap`: Rotation by `m*angle`.
"""
struct NFoldRotation{T<:Real} <: AbsRotation
    order::Int64
    m::Int64
    sym_id::Int64
    angle::T
    sym_map::LinearMap{SMatrix{2, 2, T, 4}}
end

function NFoldRotation(N, m, sym_id::Int=0; T::Type{<:Real}=Float64)
    angle = T(2*pi/N)
    mm = mod(m, N)
    sym_map = LinearMap(rotation_matrix_z(angle*mm))
    return NFoldRotation(N, mm, sym_id, angle, sym_map)
end

_with_sym_id(sym::NFoldRotation, id::Int) = NFoldRotation(sym.order, sym.m, id; T=typeof(sym.angle))

function apply_symmetry(sym::NFoldRotation, pt::SVector{2,T}) where T<:Real
    return sym.sym_map(pt)
end
function apply_symmetry(sym::NFoldRotation, pts)
    return [sym.sym_map(pt) for pt in pts]
end
apply_symmetry(sym::NFoldRotation, ::AbstractVector{<:Real}) = _throw_not_svector_point(sym)

"""
    Cn_symmetry(n; T::Type{<:Real}=Float64) → reg::SymmetryRegistry

Registers the `n-1` non-identity images of an `n`-fold rotational symmetry
(`NFoldRotation(n,1), ..., NFoldRotation(n,n-1)`) via
[`register_symmetries`](@ref), giving them fresh `sym_id`s `1,...,n-1`.

`T` should match the billiard's own numeric type parameter (defaults to
`Float64`); see [`NFoldRotation`](@ref).
"""
Cn_symmetry(n; T::Type{<:Real}=Float64) = register_symmetries((NFoldRotation(n,i; T=T) for i in 1:(n-1))...)

"""
    Dn_symmetry(n; shift=0//1, T::Type{<:Real}=Float64) → reg::SymmetryRegistry

Registers the `2n-1` non-identity elements of the dihedral group `Dₙ` in the
reconstruction order of a fundamental wedge that runs counterclockwise from the
angle `π*shift` to `π*shift + π/n`, bounded by two mirror lines.

## Description
Image `k = 1,...,2n-1` of the wedge occupies the angles `π*shift + kπ/n` to
`π*shift + (k+1)π/n`, so the images form one contiguous closed loop. Odd `k` is the
reflection about the axis `A = shift + (k+1)/(2n)` and even `k` the rotation by `2π(k/2)/n`
(an [`XYAxisReflection`](@ref) when that is a rotation by `π`). The wedge wall at
`π*shift + π/n` borders image `sym_id = 1` and the wall at `π*shift` borders image
`sym_id = 2n-1`. For `n = 2, 4` and `shift = 0` the result equals [`D2_symmetry`](@ref) and
[`D4_symmetry`](@ref).

## Arguments
* `n`: Order of the rotation subgroup, positive.
* `shift`: Angle of the first mirror line in half-turns, a `Rational`.
* `T`: Numeric type of the rotation matrices, see [`NFoldRotation`](@ref).
"""
function Dn_symmetry(n::Integer; shift::Rational=0//1, T::Type{<:Real}=Float64)
    n >= 1 || throw(ArgumentError("n must be positive; received n=$n"))
    gens = AbsSymmetry[]
    for k in 1:2n-1
        if isodd(k)
            push!(gens, AxisReflection{mod(Rational{Int}(shift) + (k + 1) // (2n), 1 // 1)}())
        elseif 2*(k ÷ 2) == n
            push!(gens, XYAxisReflection())
        else
            push!(gens, NFoldRotation(n, k ÷ 2; T=T))
        end
    end
    return register_symmetries(gens...)
end
