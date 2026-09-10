"""
    is_inside(billiard::B,pt::SVector{2,T}) where {B<:AbsBilliard,T<:Real}

Determine whether a point lies inside any simple domain of a billiard.

For each domain, the point must satisfy the inside condition for all of its
boundary curves. The point is considered inside the billiard if this holds for
at least one domain.

## Arguments
* `billiard::B`: Billiard.
* `pt::SVector{2,T}`: Point to test.

## Returns
* `Bool`: `true` if `pt` lies inside the billiard.
"""
function is_inside(billiard::B,pt::SVector{2,T}) where {B<:AbsBilliard,T<:Real}
    d=[all(is_inside(domain,pt)) for domain in get_all_domains(billiard)]
    return any(d)
end

"""
    is_inside(billiard::B,pts::AbstractArray) where {B<:AbsBilliard}

Determine whether each point in `pts` lies inside a billiard.

## Arguments
* `billiard::B`: Billiard.
* `pts::AbstractArray`: Points to test.

## Returns
* Array containing the result of `is_inside(billiard,pt)` for each point.
"""
function is_inside(billiard::B,pts::AbstractArray) where {B<:AbsBilliard}
    return [is_inside(billiard,pt) for pt in pts]
end

"""
    is_inside(domain::D,pt::SVector{2,T}) where {D<:AbsDomain,T<:Real}

Evaluate the inside condition of every boundary curve of a domain at a point.

## Arguments
* `domain::D`: Domain.
* `pt::SVector{2,T}`: Point to test.

## Returns
* Vector containing the inside condition for each boundary curve.
"""
function is_inside(domain::D,pt::SVector{2,T}) where {D<:AbsDomain,T<:Real}
    d=[is_inside(crv,pt) for crv in domain.boundary]
    return d
end

"""
    is_inside(domain::D,pts::AbstractArray) where {D<:AbsDomain}

Evaluate the inside condition of every boundary curve of a domain for multiple
points.

## Arguments
* `domain::D`: Domain.
* `pts::AbstractArray`: Points to test.

## Returns
* Matrix whose columns contain the inside conditions associated with the domain
  boundary curves.
"""
function is_inside(domain::D,pts::AbstractArray) where {D<:AbsDomain}
    d=[is_inside(crv,pts) for crv in domain.boundary]
    return reduce(hcat,d)
end

"""
    is_inside(curve::C,pt::SVector{2,T}) where {C<:AbsCurve,T<:Real}

Determine whether a point satisfies the domain condition associated with a
boundary curve.

A point is considered inside when

    domain_fun(curve,pt) <= 0.

## Arguments
* `curve::C`: Boundary curve.
* `pt::SVector{2,T}`: Point to test.

## Returns
* `Bool`: Result of the curve domain test.
"""
function is_inside(curve::C,pt::SVector{2,T}) where {C<:AbsCurve,T<:Real}
    return domain_fun(curve,pt).<=zero(eltype(pt))
end

"""
    is_inside(curve::C,pts::AbstractArray) where {C<:AbsCurve}

Determine whether multiple points satisfy the domain condition associated with a
boundary curve.

## Arguments
* `curve::C`: Boundary curve.
* `pts::AbstractArray`: Points to test.

## Returns
* Array of Boolean values obtained from the curve domain function.
"""
function is_inside(curve::C,pts::AbstractArray) where {C<:AbsCurve}
    let
        d=domain_fun(curve,pts)
        return d.<=zero(eltype(pts[1]))
    end
end

"""
    domain_gradient_vector(curve::C,pt::SVector{2,T}) where {C<:AbsCurve,T<:Real}

Compute the gradient of a curve's domain function at a point.

The gradient of the implicit domain function gives the local normal direction.

## Arguments
* `curve::C`: Boundary curve.
* `pt::SVector{2,T}`: Point at which to evaluate the gradient.

## Returns
* Gradient vector of `domain_fun(curve,pt)`.
"""
function domain_gradient_vector(curve::C,pt::SVector{2,T}) where {C<:AbsCurve,T<:Real}
    f(r)=domain_fun(curve,r)
    g=ForwardDiff.gradient(f,pt)
    return g
end

"""
    domain_gradient_vector(curve::C,pts::AbstractArray) where {C<:AbsCurve}

Compute the gradient of a curve's domain function at multiple points.

## Arguments
* `curve::C`: Boundary curve.
* `pts::AbstractArray`: Points at which to evaluate the gradient.

## Returns
* Vector containing the domain-function gradient at each point.
"""
function domain_gradient_vector(curve::C,pts::AbstractArray) where {C<:AbsCurve}
    f(r)=domain_fun(curve,r)
    gs=[ForwardDiff.gradient(f,pt) for pt in pts]
    return gs
end