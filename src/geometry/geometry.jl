"""
    is_inside(billiard::AbsBilliard, pt::SVector{2}) → flag::Bool
    is_inside(billiard::AbsBilliard, pts::AbstractArray) → flags::Vector{Bool}

`true` when `pt` lies in the fundamental domain of `billiard`, i.e. inside at
least one of its simple domains. Points are in the fundamental frame.
"""
function is_inside(billiard::B, pt::SVector{2,T}) where {B<:AbsBilliard, T<:Real}
    return _is_inside_domain(billiard.fundamental_domain, pt)
end

_is_inside_domain(domain::AbsDomain, pt) = any(all(is_inside(d, pt)) for d in get_all_domains(domain))
_is_inside_domain(domain::AbsMultiplyConnectedDomain, pt) = is_inside(domain, pt)

function is_inside(billiard::B, pts::AbstractArray) where {B<:AbsBilliard}
    return [is_inside(billiard, pt) for pt in pts]
end

"""
    is_inside(domain::AbsMultiplyConnectedDomain, pt) → flag::Bool

`true` when `pt` lies in some outer subdomain and in no hole. Holes are in the
full physical frame, so the test is valid for a fundamental-frame point and
for its symmetry images.
"""
function is_inside(domain::AbsMultiplyConnectedDomain, pt::SVector{2,T}) where {T<:Real}
    any(all(is_inside(sd, pt)) for sd in domain.subdomains) || return false
    return !any(is_inside(h, pt) for h in domain.holes)
end

function is_inside(domain::AbsMultiplyConnectedDomain, pts::AbstractArray)
    return [is_inside(domain, pt) for pt in pts]
end


"""
    is_inside(domain::AbsDomain, pt::SVector{2}) → flags::Vector{Bool}
    is_inside(domain::AbsDomain, pts::AbstractArray) → flags::Matrix{Bool}

For a simple `domain`, test `pt` against each curve of `domain.boundary`:
`flags[i]` is `true` when `pt` is on the inner side of the `i`-th curve
according to [`is_inside`](@ref) for curves. For an array of points the
result has one column per curve. A point is in the domain when all flags are
`true`; composite domains are handled by the billiard method.
"""
function is_inside(domain::D, pt::SVector{2,T}) where {D<:AbsDomain, T<:Real}
    d = [is_inside(crv, pt) for crv in domain.boundary]
    return d
end

function is_inside(domain::D, pts::AbstractArray) where {D<:AbsDomain}
    d = [is_inside(crv, pts) for crv in domain.boundary]
    return  reduce(hcat,d)
end


#check if points inside for general curves
"""
    is_inside(crv::AbsCurve, pt::SVector{2}) → flag::Bool
    is_inside(crv::AbsCurve, pts::AbstractArray) → flags::Vector{Bool}

`true` when [`domain_fun`](@ref)`(crv, pt) ≤ 0`, i.e. `pt` is on the inner side
of the curve (or on it).
"""
function is_inside(curve::C, pt::SVector{2,T}) where {C<:AbsCurve, T<:Real}
    return domain_fun(curve, pt) .<= zero(eltype(pt)) 
end

function is_inside(curve::C, pts::AbstractArray) where {C<:AbsCurve}
    let
    d = domain_fun(curve, pts)
    return d .<= zero(eltype(pts[1])) 
    end
end

#gradient of domain_function gives normal direcrion
"""
    domain_gradient_vector(crv::AbsCurve, pt::SVector{2}) → g::SVector{2}
    domain_gradient_vector(crv::AbsCurve, pts::AbstractArray) → gs::Vector

Gradient of [`domain_fun`](@ref) with respect to the point, computed by
automatic differentiation. Since `domain_fun` increases away from the interior,
it points along the outward normal of the curve.
"""
function domain_gradient_vector(curve::C, pt::SVector{2,T}) where {C<:AbsCurve, T<:Real}
    f(r) = domain_fun(curve, r)
    g = ForwardDiff.gradient(f, pt)
    return g
end

function domain_gradient_vector(curve::C, pts::AbstractArray) where {C<:AbsCurve}
    f(r) = domain_fun(curve, r)
    gs = [ForwardDiff.gradient(f, pt) for pt in pts]
    return gs
end

