
"""
    LinearNodes{T} <: AbsSampler

Sampler with the midpoints of `N` equal intervals of `[0,1]`, i.e. the composite
midpoint rule. Construct with `LinearNodes(T = Float64)`.
"""
struct LinearNodes{T<:Real} <: AbsSampler 
end 
LinearNodes(::Type{T}=Float64) where {T<:Real} = LinearNodes{T}()

"""
    sample_points(sampler::AbsSampler, N::Int) → (t, dt)

Parameter nodes `t` in `[0,1]` and quadrature weights `dt` for a boundary curve.

## Description
* [`LinearNodes`](@ref): `N` midpoints with equal weights `1/N`.
* [`GaussLegendreNodes`](@ref): `N` Gauss–Legendre nodes and weights mapped from `[-1,1]` to `[0,1]`.
* [`FourierNodes`](@ref): nodes equally spaced in arc length along several consecutive curves, returned as one vector of nodes and weights per curve; `N` is the requested total number of nodes and may be increased, see [`FourierNodes`](@ref).

## Arguments
* `sampler`: Node distribution.
* `N`: Number of nodes (total over all curves for [`FourierNodes`](@ref)).

## Returns
* `t`: Nodes in `[0,1]`.
* `dt`: Weights, so that `sum(f.(t) .* dt)` approximates the integral of `f` over `[0,1]`. For [`FourierNodes`](@ref) they are spacings in the curve parameter.
"""
function sample_points(sampler::LinearNodes{T}, N::Int) where {T<:Real}
    t = midpoints(range(zero(T),one(T),length = (N+1)))
    dt = diff(range(zero(T),one(T),length =(N+1)))
    return t, dt
end

"""
    GaussLegendreNodes{T} <: AbsSampler

Sampler with the `N` Gauss–Legendre nodes and weights on `[0,1]`. Construct
with `GaussLegendreNodes(T = Float64)`; the nodes are computed in `Float64` and
converted to `T`.
"""
struct GaussLegendreNodes{T<:Real} <: AbsSampler 
end 
GaussLegendreNodes(::Type{T}=Float64) where {T<:Real} = GaussLegendreNodes{T}()

function sample_points(sampler::GaussLegendreNodes{T}, N::Int) where {T<:Real}
    # `gausslegendre` (FastGaussQuadrature) always returns `Float64` nodes/weights;
    # convert to `T` afterwards so callers get a type-stable result for the
    # sampler's declared precision.
    x, w = gausslegendre(N)
    t = T(0.5) .* T.(x)  .+ T(0.5)
    dt = T.(w) .* T(0.5) 
    return t, dt
end

"""
    FourierNodes{T} <: AbsSampler

Sampler placing nodes at equal arc-length spacing along a sequence of
consecutive curves, as used for periodic (trapezoidal) discretization of the
connected boundary.

## Description
For a requested count `N`, the total number of nodes is `M = N`, or
`M = nextprod(primes, N)` when `primes` is given, which makes `M` a product of
small primes suitable for FFTs. The spacing is `L/M` with `L` the total length,
and the first node of each curve continues the grid from the end of the
previous curve.

## Attributes
* `primes::Union{Vector{Int64},Nothing}`: Allowed prime factors of the node count, or `nothing` for no adjustment.
* `lengths::Vector{T}`: Arc length of each curve. The constructor `FourierNodes(lengths)` sets `primes = nothing`.
"""
struct FourierNodes{T<:Real} <: AbsSampler
    primes::Union{Vector{Int64},Nothing}
    lengths::Vector{T}
end 

FourierNodes(lengths::Vector{T}) where {T<:Real} = FourierNodes{T}(nothing, lengths)

function sample_points(sampler::FourierNodes, N::Int)
    T = typeof(sampler.lengths[1])
    if isnothing(sampler.primes) 
        M = N
    else
        M = nextprod(sampler.primes,N)
    end

    ts = Vector{Vector{T}}(undef,0)
    dts = Vector{Vector{T}}(undef,0)
    t::Vector{T} = Vector{T}(undef,0)
    dt::Vector{T} = Vector{T}(undef,0)

    crv_lengths::Vector{T} = sampler.lengths
    L::T = sum(crv_lengths)

    start::T = zero(T)
    dt_end::T = zero(T)
    ds::T = zero(T)
    ncurves = length(crv_lengths)
    for (i, l) in enumerate(crv_lengths)
        ds = L/(l*M) 
        t = collect(range(start*ds,one(T),step=ds))
        dt_end = one(T) - t[end]
        start = (ds - dt_end)/ds
        push!(ts,t)
        dt = diff(t)
        # The last node of a curve (except the final one) is followed by the first node of the next curve at
        # arclength distance L/M, so its weight is the full node spacing; the final curve wraps onto t = 0.
        push!(dt, i < ncurves ? ds : dt_end)
        push!(dts,dt)
    end
    return ts,dts
end
