

struct LinearNodes <: AbsSampler 
end 

function sample_points(sampler::LinearNodes, N::Int)
    t = midpoints(range(0,1.0,length = (N+1)))
    dt = diff(range(0,1.0,length =(N+1)))
    return t, dt
end

struct GaussLegendreNodes <: AbsSampler 
end 

function sample_points(sampler::GaussLegendreNodes, N::Int)
    x, w = gausslegendre(N)
    t = 0.5 .* x  .+ 0.5
    dt = w .* 0.5 
    return t, dt
end

#=a
#this one is not working yet
function chebyshev_nodes(N::Int)
    x = [cos((2*i-1)/(2*N)*pi) for i in 1:N]
    t = 0.5 .* x  .+ 0.5
    dt = ones(N)  #wrong
    return t, dt
end
=#

struct FourierNodes{T<:Real} <: AbsSampler
    primes::Union{Vector{Int64},Nothing}
    lengths::Vector{T}
end 

# Simple outer constructor
FourierNodes(lengths::Vector{T}) where {T<:Real} = FourierNodes{T}(nothing, lengths)


#TODO thread safe
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
    for l in crv_lengths
        ds = L/(l*M) 
        #println(start*ds)
        t = collect(range(start*ds,one(T),step=ds))
        #println(t)
        dt_end = one(T) - t[end]
        start = (ds - dt_end)/ds
        push!(ts,t)
        dt = diff(t)
        push!(dt,dt_end)
        push!(dts,dt)
    end
    return ts,dts
end


function random_interior_points(billiard::AbsBilliard, N::Int; grd::Int = 1000)
    xlim,ylim = boundary_limits(billiard.fundamental_boundary; grd=grd)
    dx =  xlim[2] - xlim[1]
    dy =  ylim[2] - ylim[1]
    pts = []
 
    #println(length(pts))
    while length(pts)<N
        x = (dx .* rand() .+ xlim[1]) 
        y = (dy .* rand() .+ ylim[1])
        pt = SVector(x,y)
        if is_inside(billiard, [pt])[1] #rework this
            push!(pts,pt)
        end
    end
    return pts
end

#=
#needs some work
function fourier_nodes(N::Int; primes=(2,3,5)) #starts at 0 ends at 
    if primes == false
        M = N
    else
        M = nextprod(primes,N)
    end
    t = collect(i/M for i in 0:(M-1))
    dt = diff(t)
    dt = push!(dt,dt[1])
    return t, dt
end

function fourier_nodes(N::Int, crv_lengths; primes=(2,3,5)) #starts at 0 ends at 
    if primes == false
        M = N
    else
        M = nextprod(primes,N)
    end
    L = sum(crv_lengths)
    ts =Vector{Vector{typeof(L)}}(undef,0)
    dts =Vector{Vector{typeof(L)}}(undef,0)
    start = 0.0
    for l in crv_lengths
        ds = L/(l*M) 
        println(start*ds)
        t = collect(range(start*ds,1.0,step=ds))
        #println(t)
        dt_end = 1.0 - t[end]
        start = (ds - dt_end)/ds
        push!(ts,t)
        dt = diff(t)
        push!(dt,dt_end)
        push!(dts,dt)
    end
    return ts,dts
end
=#