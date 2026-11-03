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

"""
    verify_transparent_bcs(billiard::AbsBilliard; n = 9, δ = 1e-6, verbose = true) → ok::Bool

Check that every [`Transparent`](@ref) curve of `billiard` really joins the two
subdomains it claims to join.

## Description
For each transparent curve of a simple domain `A` with `bc = Transparent(next_id)`,
`n` interior parameter values `t` are sampled. With `p = curve(crv, t)` and the
outward unit normal `n̂ = normal_vec(crv, t)`, the test requires
* `p - δ n̂` lies in `A` (the normal points out of `A`),
* `p + δ n̂` lies in the domain `next_id` and not in `A`.

Failures (wrong `next_id`, missing domain, flipped curve orientation, overlapping
or gapped subdomains) are reported with `@warn` when `verbose`.

## Arguments
* `billiard`: Billiard to check; points are in the fundamental frame.

## Keyword arguments
* `n::Int = 9`: Number of samples per curve, at `t = i/(n+1)`.
* `δ = 1e-6`: Step along the normal.
* `verbose::Bool = true`: Emit a warning per failing sample.

## Returns
* `ok`: `true` if all samples pass (also when there are no transparent curves).
"""
function verify_transparent_bcs(billiard::B; n::Int=9, δ::Real=1e-6, verbose::Bool=true) where {B<:AbsBilliard}
    ok = true
    for dom in get_all_domains(billiard), crv in dom.boundary
        crv.bc isa Transparent || continue
        other = get_domain(billiard, crv.bc.next_id)
        if other === nothing
            verbose && @warn "Transparent curve (domain $(dom.id), segment $(crv.segment_id)) points to missing domain $(crv.bc.next_id)"
            ok = false
            continue
        end
        for i in 1:n
            t = i / (n + 1)
            p = curve(crv, t)
            nv = normal_vec(crv, t)
            pin = p - δ * nv
            pout = p + δ * nv
            in_self = all(is_inside(dom, pin))
            out_self = !all(is_inside(dom, pout))
            in_other = all(is_inside(other, pout))
            if !(in_self && out_self && in_other)
                verbose && @warn "Transparent curve (domain $(dom.id), segment $(crv.segment_id)) fails at t=$t" point=p inside_own=in_self outside_own=out_self inside_next=in_other next_id=crv.bc.next_id
                ok = false
            end
        end
    end
    return ok
end

