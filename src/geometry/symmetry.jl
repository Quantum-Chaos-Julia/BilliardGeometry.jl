
ident = IdentityTransformation()
reflect_x = LinearMap(SMatrix{2,2}([-1.0 0.0;0.0 1.0]))
reflect_y = LinearMap(SMatrix{2,2}([1.0 0.0;0.0 -1.0]))
reflect_diag = LinearMap(SMatrix{2,2}([0.0 1.0;1.0 0.0]))
reflect_antidiag = LinearMap(SMatrix{2,2}([0.0 -1.0;-1.0 0.0]))
reflect_xy = reflect_x ∘ reflect_y

"""
    XAxisReflection(parity_y=-1)

Reflection across the x-axis, `(x,y) -> (x,-y)`.

## Arguments
* `parity_y`: Irreducible-representation parity associated with the reflected `y` coordinate.
"""
struct XAxisReflection <: AbsReflection
    parity_y::Int
end
XAxisReflection() = XAxisReflection(-1)

"""
    YAxisReflection(parity_x=-1)

Reflection across the y-axis, `(x,y) -> (-x,y)`.

## Arguments
* `parity_x`: Irreducible-representation parity associated with the reflected `x` coordinate.
"""
struct YAxisReflection <: AbsReflection
    parity_x::Int
end
YAxisReflection() = YAxisReflection(-1)

"""
    XYAxisReflection(parity_x=-1,parity_y=-1)

Reflection across both coordinate axes, `(x,y) -> (-x,-y)`.

## Arguments
* `parity_x`: Irreducible-representation parity associated with the reflected `x` coordinate.
* `parity_y`: Irreducible-representation parity associated with the reflected `y` coordinate.
"""
struct XYAxisReflection <: AbsReflection 
    parity_x::Int
    parity_y::Int
end
XYAxisReflection() = XYAxisReflection(-1,-1)

"""
    DiagonalReflection(parity=-1)

Reflection across the `y=x` diagonal, `(x,y) -> (y,x)`.
"""
struct DiagonalReflection <: AbsReflection
    parity::Int
end
DiagonalReflection() = DiagonalReflection(-1)

"""
    AntiDiagonalReflection(parity=-1)

Reflection across the `y=-x` anti-diagonal, `(x,y) -> (-y,-x)`.
"""
struct AntiDiagonalReflection <: AbsReflection
    parity::Int
end
AntiDiagonalReflection() = AntiDiagonalReflection(-1)

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

D2_symmetry = [YAxisReflection(), XYAxisReflection(), XAxisReflection()]

abstract type AbsRotation <: AbsSymmetry end

# Counter-clockwise 2x2 rotation matrix, used instead of `LinearMap(RotZ(θ))`
# (Rotations.jl's `RotZ` is a 3x3 rotation and cannot convert to the
# `LinearMap{SMatrix{2,2,Float64,4}}` field type below).
@inline function rotation_matrix_z(θ::T) where T<:Real
    s, c = sincos(θ)
    return SMatrix{2,2,T}(c, s, -s, c)
end

"""
    NFoldRotation(N,m,sector=0)

One nontrivial image of an `N`-fold rotational symmetry: counter-clockwise
rotation by `2π*m/N`.

## Arguments
* `N`: Order of the rotational symmetry.
* `m`: Power of the fundamental rotation.
* `sector`: Irreducible-representation sector.
"""
struct NFoldRotation <: AbsRotation
    order::Int64
    m::Int64
    sector::Int64
    angle::Float64
    sym_map::LinearMap{SMatrix{2, 2, Float64, 4}}
end

function NFoldRotation(N, m, sector::Int=0)
    angle = 2*pi/N
    mm = mod(m, N)
    sym_map = LinearMap(rotation_matrix_z(angle*mm))
    return NFoldRotation(N, mm, mod(sector,N), angle, sym_map)
end

function apply_symmetry(sym::NFoldRotation, pt::SVector{2,T}) where T<:Real
    return sym.sym_map(pt)
end
function apply_symmetry(sym::NFoldRotation, pts)
    return [sym.sym_map(pt) for pt in pts]
end

Cn_symmetry(n, sector::Int=0) = [NFoldRotation(n,i,sector) for i in 1:(n-1)]

"""
    symmetry_irrep_character(::Type{T}, sym) where T<:Real → χ::Complex{T}

Returns the one-dimensional irreducible-representation factor associated with
a symmetry image: the parity `±1` for a `<:AbsReflection`, or the character
factor `exp(2πi*sector*m/N)` for a `<:AbsRotation` image.
"""
symmetry_irrep_character(::Type{T}, sym::XAxisReflection) where T<:Real = Complex{T}(sym.parity_y)
symmetry_irrep_character(::Type{T}, sym::YAxisReflection) where T<:Real = Complex{T}(sym.parity_x)
symmetry_irrep_character(::Type{T}, sym::XYAxisReflection) where T<:Real = Complex{T}(sym.parity_x*sym.parity_y)
symmetry_irrep_character(::Type{T}, sym::DiagonalReflection) where T<:Real = Complex{T}(sym.parity)
symmetry_irrep_character(::Type{T}, sym::AntiDiagonalReflection) where T<:Real = Complex{T}(sym.parity)
symmetry_irrep_character(::Type{T}, sym::NFoldRotation) where T<:Real = cis(T(2*pi)*T(sym.sector*sym.m)/T(sym.order))