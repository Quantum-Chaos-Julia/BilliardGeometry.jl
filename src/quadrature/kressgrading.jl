#  Reference:
#  R. Kress, "Boundary integral equations in time-harmonic acoustic scattering",
#  Math. Comput. Modelling 15 (1991), 229-243.
#  FFT-based construction via Alex Barnett's mpspack (ifft -> circulant kernel).

"""
    kress_R_even!(R0::AbstractMatrix{T}) where {T<:Real}

Constructs the periodic Kress logarithmic correction matrix in place for an
even number of nodes `N = 2n`, corresponding to the periodic logarithmic
kernel `log(4sin²((t-τ)/2))`. The first column is computed spectrally via an
inverse FFT of the Fourier coefficients of the kernel (with the Nyquist
correction term added explicitly); the remaining columns follow by circulant
shifts.
"""
function kress_R_even!(R0::AbstractMatrix{T}) where {T<:Real}
    N = size(R0,1)
    n = N÷2
    twopi = T(2*pi)
    a = zeros(Complex{T}, N)
    @inbounds for m in 1:(n-1)
        a[m+1] = 1/m
        a[N-m+1] = 1/m
    end # a[n+1] stays 0 (no 1/n Nyquist term)
    rjn = real(FFTW.ifft(a))
    ks = 0:(N-1)
    alt = (-1).^ks
    @. R0[:,1] = -twopi*rjn - (2*twopi/(N^2))*alt
    @inbounds for j in 2:N
        @views R0[:,j] .= circshift(R0[:,j-1], 1)
    end
    return nothing
end

"""
    kress_R_odd!(R0::AbstractMatrix{T}) where {T<:Real}

Constructs the periodic Kress logarithmic correction matrix in place for an
odd number of nodes `N = 2n-1`. Since an odd periodic grid has no Nyquist
frequency, the first column is given by the symmetric Fourier sum alone; the
remaining columns follow by circulant shifts.
"""
function kress_R_odd!(R0::AbstractMatrix{T}) where {T<:Real}
    N = size(R0,1)
    n = (N-1)÷2
    twopi = T(2*pi)
    a = zeros(Complex{T}, N)
    @inbounds for m in 1:n
        a[m+1] = 1/m
        a[N-m+1] = 1/m
    end
    rjn = real(FFTW.ifft(a))
    @. R0[:,1] = -twopi*rjn
    @inbounds for j in 2:N
        @views R0[:,j] .= circshift(R0[:,j-1], 1)
    end
    return nothing
end

"""
    kress_R!(R0::AbstractMatrix{T}) where {T<:Real}

Constructs the periodic Kress logarithmic correction matrix in place,
dispatching to [`kress_R_even!`](@ref) or [`kress_R_odd!`](@ref) according to
the parity of `size(R0,1)`.
"""
function kress_R!(R0::AbstractMatrix{T}) where {T<:Real}
    iseven(size(R0,1)) ? kress_R_even!(R0) : kress_R_odd!(R0)
    return nothing
end

"""
    s_mid(k::Int, N::Int) → σ

Periodic midpoint node `σ_k = 2π(k-1/2)/N` used for uniform (ungraded)
periodic trapezoidal boundary discretizations.
"""
@inline s_mid(k::Int, N::Int) = 2*pi*(k-0.5)/N

################################################################################
################## SINGLE-CORNER KRESS GRADING TRANSFORM ######################
################################################################################
# Kress-type grading transform for a single closed curve with one corner at
# σ=0≡2π: t=w(σ) clusters nodes at the corner while remaining smooth and
# strictly increasing away from it.

@inline function _kress_v(s::T, q::T) where {T<:Real}
    x = (pi-s)/pi
    return (inv(q)-T(0.5))*x^3 + inv(q)*((s-pi)/pi) + T(0.5)
end
@inline function _kress_vprime(s::T, q::T) where {T<:Real}
    x = (pi-s)/pi
    return -(T(3)/pi)*(inv(q)-T(0.5))*x^2 + inv(q)/pi
end
@inline function _kress_vdoubleprime(s::T, q::T) where {T<:Real}
    x = (pi-s)/pi
    return T((6/(pi^2))*(inv(q)-0.5)*x)
end
@inline function _kress_w(s::T, q::T) where {T<:Real}
    twopi = T(2*pi)
    a = _kress_v(s,q)^q
    b = _kress_v(twopi-s,q)^q
    return twopi*a/(a+b)
end
@inline function _kress_wprime(s::T, q::T) where {T<:Real}
    twopi = T(2*pi)
    vs = _kress_v(s,q); vsp = _kress_vprime(s,q)
    vt = _kress_v(twopi-s,q); vtp = -_kress_vprime(twopi-s,q)
    a = vs^q; b = vt^q
    ap = q*vs^(q-1)*vsp; bp = q*vt^(q-1)*vtp
    den = a+b
    return twopi*(ap*den-a*(ap+bp))/den^2
end
@inline function _kress_wdoubleprime(s::T, q::T) where {T<:Real}
    twopi = T(2*pi)
    vs = _kress_v(s,q); vsp = _kress_vprime(s,q); vspp = _kress_vdoubleprime(s,q)
    vt = _kress_v(twopi-s,q); vtp = -_kress_vprime(twopi-s,q); vtpp = _kress_vdoubleprime(twopi-s,q)
    a = vs^q; b = vt^q
    ap = q*vs^(q-1)*vsp; bp = q*vt^(q-1)*vtp
    app = q*((q-one(T))*vs^(q-2)*vsp^2+vs^(q-1)*vspp)
    bpp = q*((q-one(T))*vt^(q-2)*vtp^2+vt^(q-1)*vtpp)
    den = a+b; denp = ap+bp
    num = ap*b-a*bp; nump = app*b-a*bpp
    return T(twopi*(nump*den-2*num*denp)/den^3)
end

@inline function _min_periodic_spacing_sorted(xs::Vector{T}) where {T<:Real}
    N = length(xs)
    N<=1 && return typemax(T)
    dmin = typemax(T)
    @inbounds for i in 1:N-1
        dmin = min(dmin, xs[i+1]-xs[i])
    end
    return min(dmin, T(2*pi)+xs[1]-xs[end])
end

"""
    kress_graded_nodes_data(::Type{T}, N::Int; q=3, minsep_tol=1e-12) where {T<:Real} → (σ, s, jac, jac2, wq)

Computes Kress-graded periodic nodes and quadrature weights for `N` points
and grading order `q`, clustering nodes at the single corner `σ=0≡2π`. If the
requested grading order would push neighboring mapped nodes below
`minsep_tol`, `q` is reduced (with a warning) until the minimum spacing
requirement is met.

## Returns
* `σ`: uniform computational nodes.
* `s`: graded physical parameter values `t(σ)`.
* `jac`: first derivative `dt/dσ`.
* `jac2`: second derivative `d²t/dσ²`.
* `wq`: quadrature weights `h*jac`.
"""
function kress_graded_nodes_data(::Type{T}, N::Int; q=3, minsep_tol=1e-12) where {T<:Real}
    qT = T(q)
    qT>one(T) || error("Require q>1 for Kress grading.")
    twopi = T(2*pi)
    h = twopi/T(N)
    δ = h/T(2)
    while qT>one(T)
        σ = Vector{T}(undef,N); s = Vector{T}(undef,N)
        a = Vector{T}(undef,N); a2 = Vector{T}(undef,N); wq = Vector{T}(undef,N)
        @inbounds for k in 1:N
            σ[k] = δ+T(k-1)*h
            s[k] = _kress_w(σ[k],qT)
            a[k] = _kress_wprime(σ[k],qT)
            a2[k] = _kress_wdoubleprime(σ[k],qT)
            wq[k] = h*a[k]
        end
        minsep = _min_periodic_spacing_sorted(s)
        minsep>=minsep_tol && return σ,s,a,a2,wq
        qnew = max(one(T),qT*0.9)
        @warn "Kress grading nodes too close; reducing q." q_old=qT q_new=qnew minsep=minsep minsep_tol=minsep_tol N=N
        qT = qnew
    end
    error("Kress grading is impossible: q reached 1 while min periodic spacing stayed below minsep_tol=$(minsep_tol).")
end

################################################################################
################### MULTI-CORNER KRESS GRADING TRANSFORM ######################
################################################################################
# Generalizes the single-corner map to a finite set of true geometric corners
# {c_j}: on each interval [c_j,c_{j+1}] a local Kress-type smoothstep clusters
# nodes at both endpoints, so that t(c_j)=c_j and dt/dσ→0 at every corner.

@inline _wrap_to_2pi(x::T) where {T<:Real} = (y=mod(x,T(2*pi)); y<zero(T) ? y+T(2*pi) : y)

@inline function _sort_unique_corners(::Type{T}, corners_in) where {T<:Real}
    isempty(corners_in) && return T[]
    cs = T[_wrap_to_2pi(T(c)) for c in corners_in]; sort!(cs)
    out = T[cs[1]]
    @inbounds for j in 2:length(cs)
        abs(cs[j]-out[end])>sqrt(eps(T)) && push!(out,cs[j])
    end
    return out
end

@inline function _kress_smoothstep(u::T, q::T) where {T<:Real}
    u = clamp(u,zero(T),one(T)); a=u^q; b=(one(T)-u)^q; return a/(a+b)
end
@inline function _kress_smoothstep_prime(u::T, q::T) where {T<:Real}
    u = clamp(u,eps(T),one(T)-eps(T))
    a = u^q; b = (one(T)-u)^q; d = a+b
    return q*u^(q-one(T))*(one(T)-u)^(q-one(T))/d^2
end
@inline function _kress_smoothstep_doubleprime(u::T, q::T) where {T<:Real}
    u = clamp(u,eps(T),one(T)-eps(T))
    vp = _kress_smoothstep_prime(u,q)
    d = u^q+(one(T)-u)^q
    dp = q*(u^(q-one(T))-(one(T)-u)^(q-one(T)))
    logder = (q-one(T))/u-(q-one(T))/(one(T)-u)-T(2)*dp/d
    return vp*logder
end

function _corner_interval(corners::Vector{T}, x::T) where {T<:Real}
    m = length(corners)
    if m==1
        return corners[1], corners[1]+T(2*pi)
    end
    j = searchsortedlast(corners,x)
    if j==0
        left = corners[end]-T(2*pi); right = corners[1]
    elseif j==m
        left = corners[end]; right = corners[1]+T(2*pi)
    else
        left = corners[j]; right = corners[j+1]
    end
    return left, right
end

"""
    multi_kress_graded_nodes_data(::Type{T}, N::Int, corners_in; q=3, minsep_tol=1e-12) where {T<:Real} → (σ, tmap, jac, jac2, wq)

Constructs a globally periodic Kress-type graded parametrization for a closed
boundary with the true geometric corners `corners_in` (given as global
periodic parameter locations in `[0,2π)`; smooth joins must not be included).
If `corners_in` is empty, returns the ungraded uniform periodic grid. As in
[`kress_graded_nodes_data`](@ref), `q` is reduced with a warning if the
requested grading pushes neighboring mapped nodes below `minsep_tol`.

## Returns
* `σ`: uniform computational nodes.
* `tmap`: graded physical parameter values `t(σ)`.
* `jac`: first derivative `dt/dσ`.
* `jac2`: second derivative `d²t/dσ²`.
* `wq`: quadrature weights `h*jac`.
"""
function multi_kress_graded_nodes_data(::Type{T}, N::Int, corners_in; q=3, minsep_tol=1e-12) where {T<:Real}
    qT = T(q); qT>one(T) || error("Require q>1.")
    corners = _sort_unique_corners(T,corners_in)
    twopi = T(2*pi)
    h = twopi/T(N)
    σ = Vector{T}(undef,N)
    δ = h/2
    @inbounds for k in 1:N
        σ[k] = _wrap_to_2pi(δ+T(k-1)*h)
    end
    sort!(σ)
    if isempty(corners)
        return σ, copy(σ), ones(T,N), zeros(T,N), fill(h,N)
    end
    while qT>one(T)
        tmap = Vector{T}(undef,N); jac = Vector{T}(undef,N)
        jac2 = Vector{T}(undef,N); wq = Vector{T}(undef,N)
        @inbounds for i in 1:N
            x = σ[i]
            left,right = _corner_interval(corners,x)
            L = right-left
            u = (x-left)/L
            v = _kress_smoothstep(u,qT)
            vp = _kress_smoothstep_prime(u,qT)
            vpp = _kress_smoothstep_doubleprime(u,qT)
            tmap[i] = _wrap_to_2pi(left+L*v)
            jac[i] = vp
            jac2[i] = vpp/L
            wq[i] = h*jac[i]
        end
        minsep = _min_periodic_spacing_sorted(tmap)
        minsep>=minsep_tol && return σ,tmap,jac,jac2,wq
        qnew = max(one(T),qT*0.9)
        @warn "Kress grading nodes too close; reducing q." q_old=qT q_new=qnew minsep=minsep minsep_tol=minsep_tol N=N
        qT = qnew
    end
    error("Kress grading is impossible: q reached 1 while min periodic spacing stayed below minsep_tol=$(minsep_tol).")
end
