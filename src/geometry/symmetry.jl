const ident=IdentityTransformation()
const reflect_x=LinearMap(SMatrix{2,2}([-1.0 0.0;0.0 1.0]))
const reflect_y=LinearMap(SMatrix{2,2}([1.0 0.0;0.0 -1.0]))
const reflect_diag=LinearMap(SMatrix{2,2}([0.0 1.0;1.0 0.0]))
const reflect_antidiag=LinearMap(SMatrix{2,2}([0.0 -1.0;-1.0 0.0]))
const reflect_xy=reflect_x∘reflect_y
@inline function rotation_matrix_z(θ::T) where {T<:Real}
    s,c=sincos(θ)
    return SMatrix{2,2,T}(c,-s,s,c)
end

"""
    XAxisReflection(parity_y=-1)

Reflection across the x-axis.

## Arguments
* `parity_y`: Symmetry-sector parity associated with the reflected y-coordinate.

## Returns
An `XAxisReflection`.
"""
struct XAxisReflection<:AbsReflection
    parity_y::Int
end

"""
    YAxisReflection(parity_x=-1)

Reflection across the y-axis.

## Arguments
* `parity_x`: Symmetry-sector parity associated with the reflected x-coordinate.

## Returns
A `YAxisReflection`.
"""
struct YAxisReflection<:AbsReflection
    parity_x::Int
end

"""
    XYAxisReflection(parity_x=-1,parity_y=-1)

Reflection across both coordinate axes.

## Arguments
* `parity_x`: Symmetry-sector parity associated with x reflection.
* `parity_y`: Symmetry-sector parity associated with y reflection.

## Returns
An `XYAxisReflection`.
"""
struct XYAxisReflection<:AbsReflection
    parity_x::Int
    parity_y::Int
end

struct DiagonalReflection<:BilliardGeometry.AbsReflection
    parity::Int
end

struct AntiDiagonalReflection<:BilliardGeometry.AbsReflection
    parity::Int
end

"""
    CompositeReflection(reflections)
    CompositeReflection(reflections...)

Represent simultaneous reflection constraints whose generated symmetry group is
used for symmetry reduction.

Each constituent reflection carries its own one-dimensional irrep parity.
The complete symmetry group is obtained from the constituent actions by group
closure rather than by treating the collection as a single geometric
transformation.

This permits combinations such as an `XYAxisReflection` together with a
`DiagonalReflection`, which generate the full square `D₄` symmetry group when
supported by the geometry.

## Attributes
* `reflections`: Reflection constraints used as generators of the composite symmetry group.
"""
struct CompositeReflection<:AbsReflection
    reflections::Vector{AbsReflection}
end

CompositeReflection(reflections::AbstractVector{<:AbsReflection})=CompositeReflection(AbsReflection[reflections...])
CompositeReflection(reflections::AbsReflection...)=CompositeReflection(AbsReflection[reflections...])

# defaults 
XAxisReflection()=XAxisReflection(-1)
YAxisReflection()=YAxisReflection(-1)
XYAxisReflection()=XYAxisReflection(-1,-1)
DiagonalReflection()=DiagonalReflection(-1)
AntiDiagonalReflection()=AntiDiagonalReflection(-1)

"""
    apply_symmetry(sym,pt)
    apply_symmetry(sym,pts)

Apply a reflection symmetry to one point or a collection of points.

## Arguments
* `sym`: Reflection symmetry.
* `pt`: Point represented by an `SVector`.
* `pts`: Collection of points.

## Returns
The transformed point or collection of transformed points.
"""
@inline apply_symmetry(::XAxisReflection,pt::SVector{2,T}) where {T<:Real}=SVector{2,T}(reflect_y(pt))
@inline apply_symmetry(::YAxisReflection,pt::SVector{2,T}) where {T<:Real}=SVector{2,T}(reflect_x(pt))
@inline apply_symmetry(::XYAxisReflection,pt::SVector{2,T}) where {T<:Real}=SVector{2,T}(reflect_xy(pt))
@inline apply_symmetry(::DiagonalReflection,pt::SVector{2,T}) where {T<:Real}=SVector{2,T}(reflect_diag(pt))
@inline apply_symmetry(::AntiDiagonalReflection,pt::SVector{2,T}) where {T<:Real}=SVector{2,T}(reflect_antidiag(pt))
apply_symmetry(sym::AbsReflection,pts::AbstractVector)=[apply_symmetry(sym,pt) for pt in pts]

"""
    apply_symmetry_pb(sym,sym_sector,s,p,L)

Map Birkhoff coordinates between reflection-symmetry sectors.

## Arguments
* `sym`: Reflection symmetry.
* `sym_sector`: Symmetry sector in `1:4`.
* `s`: Boundary arclength coordinate.
* `p`: Momentum coordinate.
* `L`: Reference boundary length.

## Returns
The transformed pair `(s,p)`.
"""
function apply_symmetry_pb(::AbsReflection,sym_sector::Int,s::T,p::T,L::T) where {T<:Real}
    if sym_sector==1
        return s,p
    elseif sym_sector==2
        return 2*L-s,-p
    elseif sym_sector==3
        return 2*L+s,p
    elseif sym_sector==4
        return 4*L-s,-p
    end
    throw(ArgumentError("sym_sector must be in 1:4; received $sym_sector"))
end

const D2_symmetry=(YAxisReflection(),XYAxisReflection(),XAxisReflection())
abstract type AbsRotation<:AbsSymmetry end

"""
    NFoldRotation(N,m,sector=0)

One nontrivial image of an `N`-fold rotational symmetry.

## Arguments
* `N`: Order of the rotational symmetry.
* `m`: Power of the fundamental rotation.
* `sector`: Irreducible-representation sector.

## Returns
An `NFoldRotation`.
"""
struct NFoldRotation<:AbsRotation
    order::Int
    m::Int
    sector::Int
    angle::Float64
    sym_map::LinearMap{SMatrix{2,2,Float64,4}}
end

function NFoldRotation(N::Int,m::Int,sector::Int=0)
    angle=2pi/N
    mm=mod(m,N)
    return NFoldRotation(N,mm,mod(sector,N),angle,LinearMap(rotation_matrix_z(angle*mm)))
end

"""
    apply_symmetry(sym,pt)
    apply_symmetry(sym,pts)

Apply an `N`-fold rotational symmetry image.

## Arguments
* `sym`: Rotational symmetry image.
* `pt`: Point represented by an `SVector`.
* `pts`: Collection of points.

## Returns
The transformed point or collection of transformed points.
"""
@inline apply_symmetry(sym::NFoldRotation,pt::SVector{2,T}) where {T<:Real}=SVector{2,T}(sym.sym_map(pt))
apply_symmetry(sym::NFoldRotation,pts::AbstractVector)=[apply_symmetry(sym,pt) for pt in pts]

"""
    Cn_symmetry(n,sector=0)

Construct all nontrivial images of a cyclic `Cₙ` symmetry.

## Arguments
* `n`: Order of the cyclic symmetry.
* `sector`: Irreducible-representation sector.

## Returns
A vector containing rotations `1,…,n-1`.
"""
Cn_symmetry(n::Int,sector::Int=0)=[NFoldRotation(n,i,sector) for i in 1:n-1]

"""
    symmetry_irrep_character(::Type{T},sym)

Return the representation factor associated with a symmetry image.
For a `<:AbsReflection` this is the parity factor ±1, while for a `<:AbsRotation` this is the character factor `exp(2πim/N)`.

## Arguments
* `T`: Real scalar type.
* `sym`: Reflection or rotational symmetry image.

## Returns
The corresponding `Complex{T}` parity or character factor.
"""
@inline symmetry_irrep_character(::Type{T},sym::XAxisReflection) where {T<:Real}=Complex{T}(sym.parity_y)
@inline symmetry_irrep_character(::Type{T},sym::YAxisReflection) where {T<:Real}=Complex{T}(sym.parity_x)
@inline symmetry_irrep_character(::Type{T},sym::XYAxisReflection) where {T<:Real}=Complex{T}(sym.parity_x*sym.parity_y)
@inline symmetry_irrep_character(::Type{T},sym::NFoldRotation) where {T<:Real}=cis(T(2pi)*T(sym.sector*sym.m)/T(sym.order))
@inline symmetry_irrep_character(::Type{T},sym::DiagonalReflection) where {T<:Real}=Complex{T}(sym.parity)
@inline symmetry_irrep_character(::Type{T},sym::AntiDiagonalReflection) where {T<:Real}=Complex{T}(sym.parity)
