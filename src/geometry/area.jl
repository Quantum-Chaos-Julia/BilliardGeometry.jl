################################################################################
######################## AREA (GREEN'S THEOREM) ###############################
################################################################################

# Twice the oriented area contribution of the planar curve `crv`,
# ∫(x y' - y x')dt, used by `area` via Green's theorem.
@inline function _area_integral(crv::AbsCurve; rtol=1e-10)
    T = typeof(crv.length)
    f(t) = begin
        r = curve(crv, t)
        dr = tangent(crv, t)
        r[1]*dr[2] - r[2]*dr[1]
    end
    I, _ = quadgk(f, zero(T), one(T); rtol=rtol)
    return I
end

"""
    area(crv::AbsCurve; rtol=1e-10) → A::Real

Returns the geometric area enclosed by the planar boundary curve `crv`, via
Green's theorem `A = (1/2)|∫(x y' - y x')dt|`. The absolute value makes the
result independent of the boundary parametrization's orientation.
"""
@inline area(crv::AbsCurve; rtol=1e-10) = abs(_area_integral(crv; rtol=rtol))/2

"""
    area(crv::CompositeCurve; rtol=1e-10) → A::Real

Returns the geometric area enclosed by the composite boundary `crv`, summing
the Green-theorem contributions of every constituent curve before taking the
absolute value (preserving cancellation between oppositely oriented boundary
components).
"""
function area(crv::CompositeCurve; rtol=1e-10)
    T = typeof(crv.length)
    I = zero(T)
    @inbounds for subcrv in crv.subcurves
        I += _area_integral(subcrv; rtol=rtol)
    end
    return abs(I)/2
end

"""
    area(curves::AbstractVector{<:AbsCurve}; rtol=1e-10) → A::Real

Returns the geometric area enclosed by a collection of planar boundary
curves, summing oriented Green-theorem contributions before taking the
absolute value.
"""
function area(curves::AbstractVector{<:AbsCurve}; rtol=1e-10)
    isempty(curves) && return 0.0
    T = typeof(first(curves).length)
    I = zero(T)
    @inbounds for crv in curves
        I += _area_integral(crv; rtol=rtol)
    end
    return abs(I)/2
end

"""
    area(domain::AbsSimpleDomain; rtol=1e-10) → A::Real
    area(domain::AbsCompositeDomain; rtol=1e-10) → A::Real

Area of a closed simple domain from the Green-theorem integral over its full
`boundary` (symmetry walls included); for a composite domain the sum over its
subdomains (assumed non-overlapping).
"""
area(domain::AbsSimpleDomain; rtol=1e-10) = area(domain.boundary; rtol)
area(domain::AbsCompositeDomain; rtol=1e-10) = sum(area(sd; rtol) for sd in domain.subdomains)

"""
    area(hole::AbsDomainHole; rtol=1e-10) → A::Real

Positive area of a hole: the sum over its loops of the absolute Green-theorem area.
"""
area(hole::AbsDomainHole; rtol=1e-10) = sum(area(loop; rtol) for loop in boundary_components(hole))

# Floating-point type `T` of a `Billiard{T}`; read from the type so that it is
# inferable even though the boundary is stored as `Vector{AbsCurve}`.
@inline _billiard_eltype(::Bi) where {Bi<:AbsBilliard} = Bi.parameters[1]

"""
    area(billiard::AbsBilliard; rtol=1e-10) → A::Real

Returns the area of the complete physical billiard, using
[`full_boundary`](@ref) so the result is independent of any symmetry
reduction represented by the billiard's fundamental domain. For a billiard
with holes, the outer area minus the area of every hole.

## Arguments
* `billiard`: Billiard whose area is computed.

## Keyword arguments
* `rtol = 1e-10` : Relative tolerance of the Green-theorem quadrature.

## Returns
* `A`: Area of the physical region Ω.
"""
@inline area(billiard::AbsBilliard; rtol=1e-10) = _billiard_area(billiard.fundamental_domain, billiard; rtol)::_billiard_eltype(billiard)

_billiard_area(::AbsDomain, billiard::AbsBilliard; rtol=1e-10) = area(full_boundary(billiard); rtol)
function _billiard_area(domain::AbsMultiplyConnectedDomain, billiard::AbsBilliard; rtol=1e-10)
    outer = first(_full_boundary_data(billiard))
    return area(outer; rtol) - sum(area(h; rtol) for h in domain.holes; init=zero(_billiard_eltype(billiard)))
end

"""
    symmetry_reduction_factor(sym::AbsSymmetry) → n::Int
    symmetry_reduction_factor(billiard::AbsBilliard) → n::Int

Factor by which a discrete symmetry reduces the area of the billiard to that of
its fundamental domain: `2` for a single reflection, `4` for
[`XYAxisReflection`](@ref), the order for an [`NFoldRotation`](@ref). For a
billiard it is the largest factor among its registered symmetries, `1` without
symmetries. [`CompositeReflection`](@ref) is not reduced and returns `1`, which
is a conservative value for Weyl-law windows.
"""
@inline symmetry_reduction_factor(::AbsSymmetry) = 1
@inline symmetry_reduction_factor(::AxisReflection) = 2
@inline symmetry_reduction_factor(::XYAxisReflection) = 4
@inline symmetry_reduction_factor(sym::NFoldRotation) = sym.order
@inline symmetry_reduction_factor(::CompositeReflection) = 1

function maximal_symmetry(billiard::AbsBilliard)
    isempty(billiard.symmetries) && return nothing
    _, i = findmax(symmetry_reduction_factor, billiard.symmetries)
    return billiard.symmetries[i]
end

@inline symmetry_reduction_factor(billiard::AbsBilliard) = length(billiard.symmetries) + 1

"""
    fundamental_area(billiard::AbsBilliard; rtol=1e-10) → A::Real

Returns the area of the billiard's fundamental domain, obtained by dividing
the full physical [`area`](@ref) by the highest-order discrete symmetry's
[`symmetry_reduction_factor`](@ref). For a billiard with holes it is the
full area divided by the order of the symmetry group
[`get_symmetry_group`](@ref) (the outer group restricted to the elements that
leave every hole invariant), so it equals [`area`](@ref) when that group is
trivial, in particular for billiards whose only holes are `DomainHole`s.

## Arguments
* `billiard`: Billiard whose fundamental-domain area is computed.

## Keyword arguments
* `rtol = 1e-10` : Relative tolerance of the Green-theorem quadrature.

## Returns
* `A`: Area of the fundamental domain.
"""
@inline function fundamental_area(billiard::AbsBilliard; rtol=1e-10)
    return _fundamental_area(billiard.fundamental_domain, billiard; rtol)
end

@inline _fundamental_area(::AbsDomain, billiard::AbsBilliard; rtol=1e-10) = area(billiard; rtol)/symmetry_reduction_factor(billiard)
_fundamental_area(::AbsMultiplyConnectedDomain, billiard::AbsBilliard; rtol=1e-10) = area(billiard; rtol)/(length(get_symmetry_group(billiard)) + 1)

################################################################################
######################## CORNER ANGLES (WEYL LAW INPUT) #######################
################################################################################

"""
    corner_angles(billiard::AbsBilliard; fundamental::Bool=true, angle_tol=1e-8) → angles::Vector

Interior angles (in radians, matching e.g. [`Polygon`](@ref)'s `angles`
field convention) at every true corner of `billiard`'s boundary.

## Keyword arguments
* `fundamental::Bool = true`: When `true` (default), uses
  [`geometric_boundary_components`](@ref) (the fundamental-domain boundary
  actually diagonalized over by basis/BIM solvers, including any
  symmetry-wall cuts so corners where a physical curve meets a cut are
  correctly located); when `false`, uses [`full_boundary`](@ref) (the
  complete physical boundary), treated as a single component.
* `angle_tol = 1e-8`: Minimum turning angle (radians) classified as a true
  corner rather than a smooth join.
"""
function corner_angles(billiard::AbsBilliard; fundamental::Bool=true, angle_tol=1e-8)
    comps, is_hole = _corner_components(billiard.fundamental_domain, billiard, fundamental)
    T = _billiard_eltype(billiard)
    angles = T[]
    for (comp, hole) in zip(comps, is_hole)
        append!(angles, _component_corner_angles(T, comp; angle_tol=T(angle_tol), hole=hole))
    end
    return angles
end

function _corner_components(::AbsDomain, billiard::AbsBilliard, fundamental::Bool)
    comps = fundamental ? geometric_boundary_components(billiard) : [full_boundary(billiard)]
    return comps, falses(length(comps))
end
function _corner_components(::AbsMultiplyConnectedDomain, billiard::AbsBilliard, fundamental::Bool)
    fundamental && return geometric_boundary_components(billiard), boundary_component_is_hole(billiard)
    map = full_boundary_components(billiard)
    return map.curves, map.is_hole
end
