# -------------------------
# Fixed Gauss–Legendre rules
# -------------------------
const GL16_X=(
    -0.9894009349916499,-0.9445750230732326,-0.8656312023878318,-0.7554044083550030,
    -0.6178762444026438,-0.4580167776572274,-0.2816035507792589,-0.09501250983763744,
     0.09501250983763744,0.2816035507792589,0.4580167776572274,0.6178762444026438,
     0.7554044083550030,0.8656312023878318,0.9445750230732326,0.9894009349916499)
const GL16_W=(
    0.027152459411754095,0.06225352393864789,0.09515851168249279,0.12462897125553387,
    0.14959598881657673,0.16915651939500254,0.1826034150449236,0.1894506104550685,
    0.1894506104550685,0.1826034150449236,0.16915651939500254,0.14959598881657673,
    0.12462897125553387,0.09515851168249279,0.06225352393864789,0.027152459411754095)
const GL32_X=(
    -0.9972638618494816,-0.9856115115452684,-0.9647622555875064,-0.9349060759377397,
    -0.8963211557660521,-0.8493676137325699,-0.7944837959679424,-0.7321821187402897,
    -0.6630442669302152,-0.5877157572407623,-0.5068999089322294,-0.4213512761306353,
    -0.3318686022821277,-0.2392873622521371,-0.1444719615827965,-0.048307665687738316,
     0.048307665687738316,0.1444719615827965,0.2392873622521371,0.3318686022821277,
     0.4213512761306353,0.5068999089322294,0.5877157572407623,0.6630442669302152,
     0.7321821187402897,0.7944837959679424,0.8493676137325699,0.8963211557660521,
     0.9349060759377397,0.9647622555875064,0.9856115115452684,0.9972638618494816)
const GL32_W=(
    0.007018610009470096,0.01627439473090567,0.02539206530926206,0.03427386291302143,
    0.04283589802222668,0.050998059262376176,0.05868409347853555,0.06582222277636185,
    0.07234579410884852,0.07819389578707031,0.08331192422694675,0.0876520930044038,
    0.09117387869576389,0.09384439908080457,0.09563872007927486,0.0965400885147278,
    0.0965400885147278,0.09563872007927486,0.09384439908080457,0.09117387869576389,
    0.0876520930044038,0.08331192422694675,0.07819389578707031,0.07234579410884852,
    0.06582222277636185,0.05868409347853555,0.050998059262376176,0.04283589802222668,
    0.03427386291302143,0.02539206530926206,0.01627439473090567,0.007018610009470096)

const GL64_X=(
    -0.9993050417357722,-0.9963401167719553,-0.9910133714767443,-0.983336253884626,
    -0.973326827789911,-0.9610087996520538,-0.9464113748584028,-0.9295691721319396,
    -0.9105221370785028,-0.8893154459951141,-0.8659993981540929,-0.8406292962525803,
    -0.8132653151227975,-0.7839723589433414,-0.7528199072605318,-0.7198818501716109,
    -0.6852363130542332,-0.6489654712546573,-0.6111553551723933,-0.571895646202634,
    -0.5312794640198947,-0.4894031457070531,-0.4463660172534642,-0.4022701579639918,
    -0.35722015833766824,-0.31132287199021125,-0.26468716220876765,-0.21742364374000725,
    -0.1696444204239929,-0.12146281929612088,-0.07299312178779943,-0.024350292663424446,
    0.024350292663424446,0.07299312178779943,0.12146281929612088,0.1696444204239929,
    0.21742364374000725,0.26468716220876765,0.31132287199021125,0.35722015833766824,
    0.4022701579639918,0.4463660172534642,0.4894031457070531,0.5312794640198947,
    0.571895646202634,0.6111553551723933,0.6489654712546573,0.6852363130542332,
    0.7198818501716109,0.7528199072605318,0.7839723589433414,0.8132653151227975,
    0.8406292962525803,0.8659993981540929,0.8893154459951141,0.9105221370785028,
    0.9295691721319396,0.9464113748584028,0.9610087996520538,0.973326827789911,
    0.983336253884626,0.9910133714767443,0.9963401167719553,0.9993050417357722)

const GL64_W=(
    0.0017832807216964326,0.004147033260562467,0.006504457968978363,0.008846759826363949,
    0.011168139460131126,0.013463047896718644,0.01572603047602472,0.01795171577569734,
    0.020134823153530212,0.02227017380838325,0.024352702568710864,0.026377469715054655,
    0.028339672614259476,0.030234657072402478,0.03205792835485155,0.03380516183714161,
    0.035472213256882386,0.03705512854024003,0.03855015317861562,0.039953741132720336,
    0.04126256324262353,0.04247351512365358,0.043583724529323443,0.044590558163756545,
    0.04549162792741815,0.0462847965813144,0.04696818281620999,0.0475401657148303,
    0.047999388596458276,0.04834476223480293,0.04857546744150339,0.04869095700913967,
    0.04869095700913967,0.04857546744150339,0.04834476223480293,0.047999388596458276,
    0.0475401657148303,0.04696818281620999,0.0462847965813144,0.04549162792741815,
    0.044590558163756545,0.043583724529323443,0.04247351512365358,0.04126256324262353,
    0.039953741132720336,0.03855015317861562,0.03705512854024003,0.035472213256882386,
    0.03380516183714161,0.03205792835485155,0.030234657072402478,0.028339672614259476,
    0.026377469715054655,0.024352702568710864,0.02227017380838325,0.020134823153530212,
    0.01795171577569734,0.01572603047602472,0.013463047896718644,0.011168139460131126,
    0.008846759826363949,0.006504457968978363,0.004147033260562467,0.0017832807216964326)

###############################################################
# _toT(::Type{T}, X::NTuple{N,Float64}) where {T<:Real,N}
#
# PURPOSE:
# - Convert a compile-time NTuple of Float64 constants into the same NTuple
#   but with element type T (e.g. Float32, BigFloat, etc.).
#
# INPUTS:
# - T: target scalar type (Real)
# - X: NTuple{N,Float64} of constants (e.g. Gauss-Legendre nodes/weights)
#
# OUTPUT:
# - NTuple{N,T} with entries T(X[i])
###############################################################
@inline _toT(::Type{T},X::NTuple{N,Float64}) where {T<:Real,N}=ntuple(i->T(X[i]),N)

###############################################################
# _gl_rule(f, a, b, X, W)
#
# FORMULA:
# - Map reference nodes x_j ∈ [-1,1] to t_j = c + h*x_j with:
#       c = (a+b)/2,  h = (b-a)/2
# - Then:
#       ∫_a^b f(t) dt  ≈  h * Σ_{j=1..N} W[j] * f(c + h*X[j])
#
# INPUTS:
# - f: integrand (callable), assumed reasonably smooth on [a,b]
# - a,b: integration limits (Real), a <= b
# - X: NTuple{N,TT} of Gauss–Legendre nodes on [-1,1]
# - W: NTuple{N,TT} of corresponding weights
#
# OUTPUT:
# - Approximation to ∫_a^b f(t) dt
###############################################################
@inline function _gl_rule(f::F,a::T,b::T,X::NTuple{N,TT},W::NTuple{N,TT}) where {F<:Function,T<:Real,N,TT<:Real}
    c=(a+b)/2
    h=(b-a)/2
    s=zero(promote_type(T,TT))
    @inbounds for j in 1:N
        s+=W[j]*f(c+h*X[j])
    end
    return h*s
end

###############################################################
# _gl16 / _gl32 / _gl64
#
# PURPOSE:
# - Convenience wrappers for fixed-order Gauss–Legendre quadrature on [a,b]
#   using pre-tabulated nodes/weights (GL16_*, GL32_*, GL64_*).
#
# INPUTS:
# - f: integrand
# - a,b: integration limits
#
# OUTPUT:
# - Quadrature estimate with the requested order.
###############################################################
@inline _gl16(f,a::T,b::T) where {T<:Real}=_gl_rule(f,a,b,_toT(T,GL16_X),_toT(T,GL16_W))
@inline _gl32(f,a::T,b::T) where {T<:Real}=_gl_rule(f,a,b,_toT(T,GL32_X),_toT(T,GL32_W))
@inline _gl64(f,a::T,b::T) where {T<:Real}=_gl_rule(f,a,b,_toT(T,GL64_X),_toT(T,GL64_W))

###############################################################
# _cheb_nodes_weights(::Type{T}, p::Int)
#
# PURPOSE:
# - Construct Chebyshev–Lobatto nodes and their barycentric weights for
#   interpolation on the reference interval ξ ∈ [-1,1].
#
# INPUTS:
# - T: scalar type of the returned arrays (e.g. Float64)
# - p: polynomial degree (>= 1). Number of nodes is (p+1).
#
# OUTPUT:
# - ξ::Vector{T} of length (p+1):
#       ξ[j+1] = cos(pi*j/p),   j=0..p
# - w::Vector{T} of length (p+1): barycentric weights for these nodes:
#       w0 =  1/2
#       wp =  (1/2)*(-1)^p
#       wj =  (-1)^j              for j=1..p-1
#
# WHY:
# - Chebyshev–Lobatto nodes cluster near endpoints, mitigating Runge effects.
###############################################################
@inline function _cheb_nodes_weights(::Type{T},p::Int) where {T<:Real}
    ξ=Vector{T}(undef,p+1)
    w=Vector{T}(undef,p+1)
    invp=one(T)/T(p)
    @inbounds for j in 0:p
        x=T(j)*invp
        ξ[j+1]=cospi(x)
        wj=(j==0 || j==p) ? T(0.5) : one(T)
        wj*=isodd(j) ? -one(T) : one(T)
        w[j+1]=wj
    end
    return ξ,w
end

###############################################################
# _cheb_nodes_weights(p::Int)
#
# PURPOSE:
# - Convenience overload defaulting to Float64 nodes/weights.
###############################################################
@inline _cheb_nodes_weights(p::Int)=_cheb_nodes_weights(Float64,p)

#################################################################
# _map_to_panel(a::T,b::T,ξ::T) where {T<:Real}
#
# INPUTS:
# - a,b: panel endpoints in parameter t
# - ξ: Chebyshev node in [-1,1]
# OUTPUT:
# - t in [a,b] corresponding to ξ
# LOGIC:
# 1) Affine map from [-1,1] to [a,b]
#################################################################
@inline _map_to_panel(a::T,b::T,ξ::T) where {T<:Real}=(a+b)/T(2)+(b-a)/T(2)*ξ

###############################################################
# _barycentric_eval(x::T, xn::AbstractVector{T}, yn::AbstractVector{T}, w::AbstractVector{T})
#
# FORMULA (barycentric interpolation):
#     p(x) =  ( Σ_j (w[j]/(x-xn[j])) * yn[j] ) / ( Σ_j (w[j]/(x-xn[j])) )
#
# PURPOSE:
# - Evaluate the barycentric-form Lagrange interpolant at a point x,
#   given interpolation nodes xn, node values yn, and barycentric weights w.
#
# INPUTS:
# - x:  query point
# - xn: interpolation nodes (length n)
# - yn: function values at nodes (length n)
# - w:  barycentric weights (length n)
#
# OUTPUT:
# - Interpolated value p(x) of the unique degree-(n-1) polynomial satisfying
#   p(xn[j]) = yn[j] for all j.
#
###############################################################
@inline function _barycentric_eval(x::T,xn::AbstractVector{T},yn::AbstractVector{T},w::AbstractVector{T}) where {T<:Real}
    num=zero(T);den=zero(T)
    @inbounds for j in eachindex(xn)
        dx=x-xn[j]
        if dx==zero(T)
            return yn[j]
        end
        tj=w[j]/dx
        num+=tj*yn[j]
        den+=tj
    end
    return num/den
end

###############################################################
# _pchip_slopes!(d::AbstractVector{T}, x::AbstractVector{T}, y::AbstractVector{T})
#
# PURPOSE:
# - Compute PCHIP (shape-preserving) node slopes d[i] for the data (x[i], y[i]).
#
# INPUTS:
# - d: preallocated vector (length n), overwritten in-place with slopes dy/dx
# - x: strictly increasing knot positions (length n)
# - y: values at knots (length n)
#
# OUTPUT:
# - d (also modified in-place)
#
# DEFINITIONS:
# - h[i] = x[i+1] - x[i]  interval widths
# - δ[i] = (y[i+1] - y[i]) / h[i] # secant slopes on intervals
#
# LOGIC:
# 1) Compute h and δ.
# 2) Endpoints:
#    - Use a one-sided weighted secant formula, then apply sign/magnitude limiting
#      to avoid overshoot and preserve local monotonicity.
# 3) Interior nodes i=2..n-1:
#    - If δ[i-1] and δ[i] differ in sign (or either is zero), set d[i]=0.
#    - Otherwise, set d[i] to a spacing-weighted harmonic mean of δ[i-1], δ[i],
#      which is stable and shape-preserving on nonuniform grids.
#
# NOTES:
# - Designed for monotone-safe interpolation (prevents cubic overshoot).
###############################################################
@inline function _pchip_slopes!(d::AbstractVector{T},x::AbstractVector{T},y::AbstractVector{T}) where {T<:Real}
    n=length(x)
    n<2 && return d
    h=Vector{T}(undef,n-1);δ=Vector{T}(undef,n-1)
    @inbounds for i in 1:n-1
        hi=x[i+1]-x[i]
        h[i]=hi
        δ[i]=(y[i+1]-y[i])/hi
    end
    if n==2
        d[1]=δ[1];d[2]=δ[1]
        return d
    end
    @inline function _edge_slope(h1,h2,δ1,δ2)
        dd=((T(2)*h1+h2)*δ1-h1*δ2)/(h1+h2)
        if sign(dd)!=sign(δ1);return zero(T);end
        if sign(δ1)!=sign(δ2) && abs(dd)>T(3)*abs(δ1);return T(3)*δ1;end
        return dd
    end
    d[1]=_edge_slope(h[1],h[2],δ[1],δ[2])
    d[n]=_edge_slope(h[n-1],h[n-2],δ[n-1],δ[n-2])
    @inbounds for i in 2:n-1
        δm=δ[i-1];δp=δ[i]
        if δm==zero(T) || δp==zero(T) || sign(δm)!=sign(δp)
            d[i]=zero(T)
        else
            hm=h[i-1];hp=h[i]
            w1=T(2)*hp+hm
            w2=hp+T(2)*hm
            d[i]=(w1+w2)/(w1/δm+w2/δp)
        end
    end
    return d
end

###############################################################
# _pchip_eval(xq::T, x::AbstractVector{T}, y::AbstractVector{T}, d::AbstractVector{T})
#
# PURPOSE:
# - Evaluate a *PCHIP* (Piecewise Cubic Hermite Interpolating Polynomial)
#   interpolant at a query point xq.
#
# CONTEXT (why this exists in your arc-length inversion):
# - On each panel you have monotone data (s_sorted, t_sorted).
# - You want t(s) quickly and safely:
#     - monotone-safe (no overshoot that breaks bracket / Newton)
#     - O(1) per query once you know the interval
# - PCHIP gives you a C^1 piecewise cubic, but with slopes chosen to preserve
#   monotonicity when the data are monotone (via your _pchip_slopes! routine).
#
# INPUTS:
# - xq::T:
#   Coord where we want the interpolated value.
#
# - x::AbstractVector{T}:
#   Knots x[1],...,x[n] (MUST be sorted strictly increasing): x == s_sorted
#
# - y::AbstractVector{T}:
#   Function values y[i] at knots x[i]: y == t_sorted
#
# - d::AbstractVector{T}:
#   Slopes at knots, i.e. d[i] ≈ dy/dx at x[i].
#   These slopes are not arbitrary; they are computed by
#   _pchip_slopes! to enforce shape/monotonicity (avoid oscillations).
#
# OUTPUT:
# - yq::T:
#   Interpolated value y(xq) given the piecewise cubic Hermite polynomial.
#
# LOGIC:
# 1) Clamp at boundaries:
#    - If xq <= x[1], return y[1]
#    - If xq >= x[n], return y[n]
#    This avoids extrapolation and guarantees safety.
#
# 2) Locate the interval:
#    - Find i such that x[i] <= xq < x[i+1]
#      using searchsortedlast (O(log n)).
#
# 3) Map xq to the local normalized coordinate t ∈ [0,1]:
#      h  = x1 - x0
#      t  = (xq - x0)/h
#
# 4) Evaluate the cubic Hermite interpolant on [x0,x1]:
#
#    The cubic is constructed to satisfy 4 constraints:
#      p(x0)  = y0
#      p(x1)  = y1
#      p'(x0) = d0
#      p'(x1) = d1
#
#    In normalized coordinate t ∈ [0,1], the standard Hermite form is:
#      p(t) = h00(t)*y0 + h10(t)*h*d0 + h01(t)*y1 + h11(t)*h*d1
#
#    where h00,h10,h01,h11 are the *Hermite basis functions*:
#
#      h00(t) =  2t^3 - 3t^2 + 1
#      h10(t) =    t^3 - 2t^2 + t
#      h01(t) = -2t^3 + 3t^2
#      h11(t) =    t^3 -   t^2
#
#    Each basis has Kronecker boundary behavior:
#
#    Values at endpoints:
#      h00(0)=1, h00(1)=0   (selects y0)
#      h01(0)=0, h01(1)=1   (selects y1)
#      h10(0)=0, h10(1)=0   (no contribution to endpoint values)
#      h11(0)=0, h11(1)=0
#
#    Derivatives at endpoints (with respect to t):
#      h00'(0)=0, h00'(1)=0
#      h01'(0)=0, h01'(1)=0
#      h10'(0)=1, h10'(1)=0   (selects slope at left endpoint)
#      h11'(0)=0, h11'(1)=1   (selects slope at right endpoint)
#
#    Because x = x0 + h*t, we have dp/dx = (1/h)*dp/dt.
#    That is why the slope terms appear as (h*d0) and (h*d1):
#    the Hermite basis expects derivatives w.r.t. x, but is written in t.
###############################################################
@inline function _pchip_eval(xq::T,x::AbstractVector{T},y::AbstractVector{T},d::AbstractVector{T}) where {T<:Real}
    n=length(x)
    xq<=x[1] && return y[1]
    xq>=x[n] && return y[n]
    i=searchsortedlast(x,xq)
    i>=n && (i=n-1)
    @inbounds begin
        x0=x[i];x1=x[i+1]
        y0=y[i];y1=y[i+1]
        d0=d[i];d1=d[i+1]
        h=x1-x0
        t=(xq-x0)/h
        t2=t*t;t3=t2*t
        h00=T(2)*t3-T(3)*t2+one(T)
        h10=t3-T(2)*t2+t
        h01=-T(2)*t3+T(3)*t2
        h11=t3-t2
        return h00*y0+h10*h*d0+h01*y1+h11*h*d1
    end
end

struct Panel{T<:Real}
    a::T
    b::T
    L::T
    qord::Int
end

struct PanelData{T<:Real}
    panel::Panel{T}
    tnodes::Vector{T}
    wbar::Vector{T}
    xs::Vector{T}
    ys::Vector{T}
    ss::Vector{T}
    vs::Vector{T}
    L16::T
    L32::T
    s_sorted::Vector{T}
    t_sorted::Vector{T}
    dt_ds::Vector{T}
end

###############################################################
# _panel_length_and_order(::Type{T},f::F,a::T,b::T,quad_rtol::T) where {T<:Real,F<:Function}
#
# INPUTS:
# - ::Type{T}: Floating-point type for all computations (e.g. Float64)
# - f::F: Function to integrate
# - a,b: panel endpoints
# - quad_rtol: relative tolerance for GL rule agreement
#
# OUTPUT:
# - L: estimated panel length
# - qord: chosen quadrature order (16/32/64)
# - ok_quad: Bool indicating whether panel passed quadrature criterion
#
# LOGIC: 
# - Compute GL16, GL32 estimates
# - If |L32-L16| <= quad_rtol*max(1,|L32|), accept L16, qord=16
# - Else compute GL64
# - If |L64-L32| <= quad_rtol*max(1,|L64|), accept L32, qord=32
# - Else reject panel (ok_quad=false), return L64, qord=64
###############################################################
@inline function _panel_length_and_order(::Type{T},f::F,a::T,b::T,quad_rtol::T) where {T<:Real,F<:Function}
    L16=_gl16(f,a,b);L32=_gl32(f,a,b)
    if abs(L32-L16)<=quad_rtol*max(one(T),abs(L32));return L16,16,true;end
    L64=_gl64(f,a,b)
    if abs(L64-L32)<=quad_rtol*max(one(T),abs(L64));return L32,32,true;end
    return L64,64,false
end

###############################################################
# _speed(crv::AbsCurve,t::T) where {T<:Real}
#
# INPUT:
# - crv<:AbsCurve
# - t::T
# OUTPUT:
# - speed |tangent(crv,t)|
###############################################################
@inline _speed(crv::AbsCurve,t)=norm(tangent(crv,t))

###############################################################
# _build_panels(::Type{T},crv;q=8,init_panels=8,max_panels=200000,nprobe=9,quad_rtol=1e-12,speed_ratio_max=3)
#
# PURPOSE:
# - Construct an adaptive panelization of the parameter interval t∈[0,1].
# - Panels are refined (bisected) until *all* acceptance criteria are met:
#     curvature Nyquist criterion (geometric resolution)
#     fixed Gauss–Legendre length convergence (quadrature reliability)
#     speed variation bound (parameterization regularity)
#
# INPUTS:
# - ::Type{T}: Floating-point type for all computations (e.g. Float64)
#
# - crv<:AbsCurve
#
# - q::T (default 8):
#   Curvature oversampling / “Nyquist” factor.
#   Enforces roughly:
#     L_panel <= 1/(q*κmax)   when κmax>0
#   Interpretation:
#   - κmax sets the smallest radius of curvature ~ 1/κmax.
#   - We require several samples across that radius along arc-length.
#   Larger q => more panels in high-curvature regions (safer).
#
# - init_panels::Int (default 8):
#   Start with a uniform partition of [0,1] into init_panels panels.
#   These are then refined adaptively.
#
# - max_panels::Int (default 200000):
#   Hard cap to avoid infinite refinement (e.g. cusps, bad parameterization).
#
# - nprobe::Int (default 9):
#   Number of probe points used to estimate:
#     κmax  (max curvature on panel)
#     vmin,vmax (speed range on panel)
#   Probes are uniform in t on [a,b].
#
# - quad_rtol::T (default 1e-12):
#   Relative tolerance for panel-length estimation by fixed GL rules.
#   Implemented in _panel_length_and_order via GL16/GL32/(GL64 fallback):
#     accept if |L32-L16| <= quad_rtol*max(1,|L32|)
#     else if |L64-L32| <= quad_rtol*max(1,|L64|)
#     else reject (force refinement)
#
# - speed_ratio_max::T (default 3):
#   Parameterization quality criterion:
#     vmax/vmin <= speed_ratio_max
#   Large speed variation means the mapping t→arc-length is “stiff”,
#   which can hurt interpolation and inversion accuracy.
#   Tightening this produces more panels where speed changes rapidly.
#
# OUTPUT:
# - panels::Vector{Panel{T}}:
#   Final accepted panels covering [0,1], in increasing order, with:
#     Panel(a,b,L,qord)
#   where:
#     a,b   panel parameter endpoints
#     L     estimated panel arc-length (from GL rule)
#     qord  chosen quadrature order (16/32/64) for that panel length estimate
#
# ACCEPTANCE CRITERIA (per panel [a,b]):
# - Curvature criterion:
#     ok_curv = (κmax==0) OR (L <= 1/(q*κmax))
# - Speed criterion:
#     ok_speed = (vmin != 0) AND (vmax/vmin <= speed_ratio_max)
# - Quadrature criterion:
#     ok_quad = true if GL16/GL32 (or GL32/GL64) agree to quad_rtol
#
# REFINEMENT STRATEGY:
# - If any criterion fails, bisect panel at midpoint m=(a+b)/2:
#     replace [a,b] by [a,m] and [m,b]
# - Continue until all panels pass or max_panels reached.
###############################################################
function _build_panels(::Type{T},crv;q::T=T(8),init_panels::Int=8,max_panels::Int=200000,nprobe::Int=9,quad_rtol::T=T(1e-12),speed_ratio_max::T=T(3)) where {T<:Real}
    panels=Vector{Panel{T}}(undef,init_panels)
    invN=one(T)/T(init_panels)
    @inbounds for i in 1:init_panels
        a=T(i-1)*invN;b=T(i)*invN
        panels[i]=Panel{T}(a,b,zero(T),0)
    end
    f_speed(t::T)=_speed(crv,t)
    i=1
    while i<=length(panels)
        p=panels[i];a=p.a;b=p.b
        κmax=zero(T);vmin=T(Inf);vmax=zero(T)
        @inbounds for k in 0:(nprobe-1)
            t=a+(b-a)*(T(k)/T(nprobe-1))
            κ=abs(curvature(crv,t));v=f_speed(t)
            κmax=max(κmax,κ);vmin=min(vmin,v);vmax=max(vmax,v)
        end
        ok_speed=(vmin!=zero(T))&&(vmax/vmin<=speed_ratio_max)
        L,qord,ok_quad=_panel_length_and_order(T,f_speed,a,b,quad_rtol)
        ok_curv=(κmax==zero(T))||L<=one(T)/(q*κmax)
        if ok_curv && ok_speed && ok_quad
            panels[i]=Panel{T}(a,b,L,qord)
            i+=1;continue
        end
        if length(panels)>=max_panels;error("Too many panels");end
        m=(a+b)/T(2)
        panels[i]=Panel{T}(a,m,zero(T),0)
        insert!(panels,i+1,Panel{T}(m,b,zero(T),0))
    end
    return panels
end

###############################################################
# _build_panel_data(crv,panels::Vector{Panel{T}};p_cheb::Int=16)
#
# PURPOSE:
# - Build all per-panel interpolation data needed for:
#     xy(t)   via Chebyshev barycentric interpolation
#     s(t)    via Chebyshev barycentric interpolation
#     v(t)    (speed) via Chebyshev barycentric interpolation
#     t(s)    via per-panel monotone PCHIP initial guess
#
# INPUTS:
# - crv:
#   Curve object providing (at least):
#     curve(crv,t) -> SVector{2,T}  (position)
#     _speed(crv,t) -> T            (speed = ||tangent||)
#
# - panels::Vector{Panel{T}}:
#   Panelization of [0,1] into sub-intervals [a,b].
#   Each Panel stores:
#     a::T, b::T      panel endpoints in parameter t
#     L::T            panel length 
#     qord::Int       quadrature order used 
#
# - p_cheb::Int (default 16):
#   Chebyshev polynomial degree per panel.
#   Uses p_cheb+1 Chebyshev–Lobatto nodes on each panel.
#
# OUTPUT:
# - data::Vector{PanelData{T}} of length = length(panels)
#   Each PanelData stores (typical fields):
#     panel     : Panel{T} (a,b,L,qord)
#     tnodes    : Chebyshev nodes mapped to [a,b]
#     wbar      : barycentric weights for those nodes
#     xs, ys    : x(tnodes), y(tnodes)
#     ss        : s(tnodes)  (GLOBAL arc length, monotone over whole curve)
#     vs        : v(tnodes)  (speed at nodes)
#     L16, L32  : GL16/GL32 panel length diagnostics
#     s_sorted  : node arc lengths in increasing t-order (for PCHIP)
#     t_sorted  : same nodes sorted increasingly in t
#     dt_ds     : PCHIP slopes for t(s) in that panel
#
# IDEA:
# - On each panel we tabulate geometry and arc length at Chebyshev nodes.
# - Accumulate arc length incrementally between consecutive nodes
#   using fixed Gauss–Legendre GL32 on each small subinterval:
#
#       s(t_j) = s(t_{j-1}) + ∫_{t_{j-1}}^{t_j} ||r'(u)|| du  (GL32)
#
# WHY SORTING IS NECESSARY:
# - Chebyshev–Lobatto nodes are not ordered in t.
# - The incremental accumulation needs t_{j-1} < t_j.
# - Therefore:
#     perm     = sortperm(tnodes)
#     t_sorted = tnodes[perm]
#   accumulate s_sorted along t_sorted, then permute back to original order
#   so barycentric interpolation (which uses the original tnodes ordering)
#   remains consistent.
#
# PER-PANEL LOGIC:
# 1) ξ,wbar = _cheb_nodes_weights(T,p_cheb) on [-1,1]
# 2) Map ξ -> tnodes ∈ [a,b] via _map_to_panel
# 3) Evaluate curve positions:
#       xs[j], ys[j] = curve(crv, tnodes[j])
# 4) Sort nodes by increasing t:
#       perm = sortperm(tnodes)
#       t_sorted = tnodes[perm]
# 5) Incrementally accumulate arc length:
#       s_sorted[1] = s0  (global offset at panel start)
#       for j=2..end:
#           s_sorted[j] = s_sorted[j-1] + _gl32(speed, t_sorted[j-1], t_sorted[j])
# 6) Permute back to Chebyshev order:
#       ss[perm[j]] = s_sorted[j]
# 7) Evaluate speed at nodes:
#       vs[j] = _speed(crv, tnodes[j])
# 8) Diagnostics:
#       L16 = _gl16(speed,a,b)
#       L32 = _gl32(speed,a,b)
# 9) Build PCHIP slopes for inverse initial guess t(s) on this panel:
#       dt_ds = slopes for monotone cubic Hermite on (s_sorted -> t_sorted)
#       _pchip_slopes!(dt_ds, s_sorted, t_sorted)
# 10) Store PanelData and advance global offset: s0 += L32
###############################################################
function _build_panel_data(crv,panels::Vector{Panel{T}};p_cheb::Int=16) where {T<:Real}
    ξ,wbar=_cheb_nodes_weights(T,p_cheb)
    f_speed(t::T)=_speed(crv,t)
    data=Vector{PanelData{T}}(undef,length(panels))
    s0=zero(T)
    for (ip,pan) in pairs(panels)
        a=pan.a;b=pan.b
        tnodes=Vector{T}(undef,length(ξ))
        @inbounds for j in eachindex(ξ)
            tnodes[j]=_map_to_panel(a,b,ξ[j])
        end
        xs=Vector{T}(undef,length(tnodes))
        ys=Vector{T}(undef,length(tnodes))
        @inbounds for j in eachindex(tnodes)
            r=curve(crv,tnodes[j])
            xs[j]=r[1];ys[j]=r[2]
        end
        perm=sortperm(tnodes);t_sorted=tnodes[perm]
        s_sorted=similar(t_sorted);s_sorted[1]=s0;acc=s0
        @inbounds for j in eachindex(t_sorted)[2:end]
            acc+=_gl32(f_speed,t_sorted[j-1],t_sorted[j])
            s_sorted[j]=acc
        end
        ss=Vector{T}(undef,length(tnodes))
        @inbounds for j in eachindex(perm)
            ss[perm[j]]=s_sorted[j]
        end
        vs=Vector{T}(undef,length(tnodes))
        @inbounds for j in eachindex(tnodes)
            vs[j]=_speed(crv,tnodes[j])
        end
        L16=_gl16(f_speed,a,b);L32=_gl32(f_speed,a,b)
        dt_ds=Vector{T}(undef,length(t_sorted))
        _pchip_slopes!(dt_ds,s_sorted,t_sorted)
        data[ip]=PanelData{T}(pan,tnodes,wbar,xs,ys,ss,vs,L16,L32,s_sorted,t_sorted,dt_ds)
        s0+=L32
    end
    return data
end

struct PanelChebArc{T<:Real,C}
    crv::C
    panels::Vector{Panel{T}}
    pdata::Vector{PanelData{T}}
    total_length::T
    send::Vector{T}
end

#################################################################
# build_panel_cheb_arc(::Type{T},crv;q::T=T(5),init_panels::Int=8,nprobe::Int=9,
#                      quad_rtol::T=T(1e-12),speed_ratio_max::T=T(3),p_cheb::Int=16)
# INPUTS:
# - T: real type
# - crv: curve object
# - q: curvature Nyquist factor (over-sampling)
# - init_panels: initial number of panels
# - nprobe: number of probe points per panel for curvature/speed estimation
# - quad_rtol: quadrature relative tolerance
# - speed_ratio_max: maximum allowed speed ratio per panel
# - p_cheb: number of Chebyshev nodes per panel
#
# OUTPUT:
# - PanelChebArc object
#
# LOGIC:
# 1) Build panels via _build_panels
# 2) Build panel data via _build_panel_data
# 3) Compute cumulative arc lengths at panel ends
#################################################################
function build_panel_cheb_arc(::Type{T},crv;q::T=T(5),init_panels::Int=8,nprobe::Int=9,quad_rtol::T=T(1e-12),speed_ratio_max::T=T(3),p_cheb::Int=16) where {T<:Real}
    panels=_build_panels(T,crv;q=q,init_panels=init_panels,nprobe=nprobe,quad_rtol=quad_rtol,speed_ratio_max=speed_ratio_max)
    pdata=_build_panel_data(crv,panels;p_cheb=p_cheb)
    send=Vector{T}(undef,length(pdata))
    acc=zero(T)
    @inbounds for i in eachindex(pdata)
        acc+=pdata[i].L32
        send[i]=acc
    end
    totalL=send[end]
    return PanelChebArc{T,typeof(crv)}(crv,panels,pdata,totalL,send)
end

###############################################################
# _find_panel_t(panels::Vector{Panel{T}},t::T)
# INPUTS:
# - panels: vector of Panel{T} objects
# - t: parameter t in [0,1]
# OUTPUT:
# - index of panel containing t
# LOGIC: 
# 1) Binary search on panel.b values in the following way:
#    - initialize lo=1, hi=length(panels)
#    - while lo<hi: # while panels remain to search
#        - mid=(lo+hi)>>>1 # integer divide by 2 for midpoint
#        - if t <= panels[mid].b: hi=mid # search lower half 
#          else: lo=mid+1 # search upper half
#    - return lo as the index of the panel containing t
###############################################################
@inline function _find_panel_t(panels::Vector{Panel{T}},t::T) where {T<:Real}
    lo=1
    hi=length(panels)
    while lo<hi
        mid=(lo+hi)>>>1
        if t<=panels[mid].b
            hi=mid
        else
            lo=mid+1
        end
    end
    return lo
end

###############################################################
# _find_panel_s(send::Vector{T},s::T)
#
# INPUTS:
# - send: vector of panel end arc lengths
# - s: arc length parameter s in [0, total_length]
#
# OUTPUT:
# - index of panel containing s
#
# LOGIC: 
# 1) Use searchsortedfirst to locate panel
###############################################################
@inline function _find_panel_s(send::Vector{T},s::T) where {T<:Real}
    return searchsortedfirst(send,s) 
end

###############################################################
# xy(obj::PanelChebArc,tq::Real)
#
# INPUTS:
# - obj: PanelChebArc object
# - tq: parameter t in [0,1]
#
# OUTPUT:
# - (x,y) coordinates at t
#
# LOGIC: 
# 1) Clamp t to [0,1]
# 2) Locate panel containing t via binary search on panel.b values
# 3) Evaluate x(t), y(t) via barycentric interpolation on (tnodes, xs) and (tnodes, ys) in that panel
###############################################################
function xy(obj::PanelChebArc{T},tq::Real) where {T<:Real}
    t=T(tq);t=t<zero(T) ? zero(T) : (t>one(T) ? one(T) : t)
    i=_find_panel_t(obj.panels,t)          # FIX
    pd=obj.pdata[i]
    x=_barycentric_eval(t,pd.tnodes,pd.xs,pd.wbar)
    y=_barycentric_eval(t,pd.tnodes,pd.ys,pd.wbar)
    return SVector{2,T}(x,y)
end

###############################################################
# _s_of_t(obj::PanelChebArc,tq::Real)
#
# INPUTS:
# - obj: PanelChebArc object
# - tq: parameter t in [0,1]
#
# OUTPUT:
# - s in [0, total_length]
#
# LOGIC: 
# 1) Clamp t to [0,1]
# 2) Locate panel containing t via binary search on panel.b values
# 3) Evaluate s(t) via barycentric interpolation on (tnodes, ss) in that panel
###############################################################
function _s_of_t(obj::PanelChebArc{T},tq::Real) where {T<:Real}
    t=T(tq);t=t<zero(T) ? zero(T) : (t>one(T) ? one(T) : t)
    i=_find_panel_t(obj.panels,t) 
    pd=obj.pdata[i]
    return _barycentric_eval(t,pd.tnodes,pd.ss,pd.wbar)
end

###############################################################
# _t_of_s(obj::PanelChebArc,sq::T;tol::T=1e-12)
#
# INPUTS:
# - obj: PanelChebArc object
# - sq: arc length parameter s in [0, total_length]
# - tol: tolerance for Newton iteration convergence (1 step Newton)
#
# OUTPUT:
# - t in [0,1] such that s(t) ≈ sq
#
# LOGIC: 
# 1) Handle boundary cases s<=0 and s>=total_length
# 2) Locate panel containing s via searchsortedfirst on send array !!! Binary search is slow here !!!
# 3) Use PCHIP (Piecewise Cubic Hermite Interpolating Polynomial) 
# interpolation on (s_sorted, t_sorted) in that panel to get initial t
# 4) Evaluate s(t) and v(t) via barycentric interpolation at that t
#   - s(t) is evaluated on (tnodes, ss), where ss are the arc lengths
#   at Chebyshev nodes that are precomputed
#   - v(t) is evaluated on (tnodes, vs), where vs are the speeds 
#   at Chebyshev nodes that are precomputed
# 5) One step Newton: t_new = t - (s(t)-sq)/v(t)
# 6) Clamp t_new to panel [a,b] and return
###############################################################
function _t_of_s(obj::PanelChebArc{T},sq::T;tol::T=T(1e-12)) where {T<:Real}
    s=T(sq)
    s<=zero(T) && return zero(T)
    s>=obj.total_length && return one(T)
    ip=searchsortedfirst(obj.send,s)
    pd=obj.pdata[ip]
    a=pd.panel.a;b=pd.panel.b
    t=_pchip_eval(s,pd.s_sorted,pd.t_sorted,pd.dt_ds)
    t=t<a ? a : (t>b ? b : t)
    st=_barycentric_eval(t,pd.tnodes,pd.ss,pd.wbar)
    f=st-s
    abs(f)<=tol && return t
    v=_barycentric_eval(t,pd.tnodes,pd.vs,pd.wbar)
    v<=eps(T) && return t
    tnew=t-f/v
    t=tnew<a ? a : (tnew>b ? b : tnew)
    return t
end