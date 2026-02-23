function _arc_length_integrand(crv::C,t::T) where {C<:AbsCurve,T<:Real}
    f(t)=curve(crv,t)
    return norm(ForwardDiff.derivative(f,t))
end

function arc_length(crv::C,t1::T;rtol=sqrt(eps(T)),atol=zero(T)) where {C<:AbsCurve,T<:Real}
    f(t)=_arc_length_integrand(crv,t)
    return quadgk(f,zero(T),t1;rtol,atol)[1]
end

function arc_length(crv::C,ts::AbstractArray;rtol=sqrt(eps(eltype(ts))),atol=zero(eltype(ts))) where {C<:AbsCurve}
    return [arc_length(crv,t;rtol,atol) for t in ts]
end

function arc_length(crv::C,pt::SVector{2,T};rtol=sqrt(eps(T)),atol=zero(T)) where {C<:AbsCurve,T<:Real}
    t1 = invert_curve(crv, pt)
    f(t)=_arc_length_integrand(crv,t)
    return quadgk(f,zero(T),t1;rtol,atol)[1]
end

##################################
#### ARCLENGTH INTERPOLATIONS ####
##################################

@inline _clamp_t(t::T) where {T<:Real}=t<zero(T) ? zero(T) : (t>one(T) ? one(T) : t)

function construct_arc_length_interpolation(crv::LineSegment{T,BC}) where {T<:Real,BC}
    L=crv.length
    s_of_t=(tq)->begin
        return L*_clamp_t(T(tq))
    end
    t_of_s=(sq)->begin
        s=T(sq)
        s<=zero(T) && return zero(T)
        s>=L && return one(T)
        return s/L
    end
    return s_of_t,t_of_s
end

function construct_arc_length_interpolation(crv::CircleSegment{T,BC}) where {T<:Real,BC}
    L=crv.length
    s_of_t=(tq)->begin
        return L*_clamp_t(T(tq))
    end

    t_of_s=(sq)->begin
        s=T(sq)
        s<=zero(T) && return zero(T)
        s>=L && return one(T)
        return s/L
    end
    return s_of_t,t_of_s
end

# q is defined as the ration of the maximum speed to the minimum speed along the curve, used for adaptive panel refinement.
# speed_ratio_max is the maximum allowed speed ratio for the adaptive panel refinement. If the actual speed ratio exceeds this value, the algorithm will refine the panels further to ensure better accuracy in the arc length interpolation.
# init_panels is the initial number of panels used for the adaptive refinement process. The algorithm will start with this number of panels and then refine them based on the speed ratio until the desired accuracy is achieved.
# nprobe is the number of probe points used to estimate the speed ratio along each panel. The algorithm will evaluate the speed at these probe points to determine if further refinement is needed.
# tol_newton is the tolerance for the Newton's method used to find the inverse mapping from arc length to parameter t. The algorithm will iterate until the difference between successive approximations is less than this tolerance, ensuring that the inverse mapping is accurate.
function _construct_arc_length_interpolation_CHEBYSHEV(::Type{T},crv::C;q=3.0,p::Int=8,quad_rtol=1e-8,speed_ratio_max=3.0,init_panels::Int=8,nprobe::Int=9,tol_newton=1e-11) where {BC,T<:Real,C<:AbsPolarCurve{BC}}
    obj=build_panel_cheb_arc(T,crv;q=q,init_panels=init_panels,nprobe=nprobe,quad_rtol=quad_rtol,speed_ratio_max=speed_ratio_max,p_cheb=p)
    s_of_t=(tq)->_s_of_t(obj,T(tq))
    t_of_s=(sq)->_t_of_s(obj,T(sq);tol=tol_newton)
    return s_of_t,t_of_s
end

# CubicSpline interpolation, reasonably fast butr not as accurate as the chebyshev adaptive panel refinement for the same number of panels. Not recommended for curves with large curvature variations.
function _construct_arc_length_interpolation_CUBIC_SPLINE(crv::C;target_prec=1e-8,test_points::Int=100,verbose=false) where {C<:AbsCurve}
    max_iters=50
    errors_s=Float64[]
    errors_t=Float64[]
    n_samples=100;best=Inf
    ts=collect(range(0.0,1.0;length=test_points))
    s_true=arc_length(crv,ts)
    local s_of_t,t_of_s,s_samples,t_samples
    for i in 0:max_iters-1
        t_samples=collect(range(0.0,1.0;length=n_samples))
        s_samples=zeros(Float64,n_samples)
        f(t)=_arc_length_integrand(crv,t)
        for k in 2:n_samples
            s_samples[k],_=quadgk(f,0.0,t_samples[k];rtol=1e-2*target_prec)
        end
        s_of_t=CubicSpline(s_samples,t_samples) # s(t)
        t_of_s=CubicSpline(t_samples,s_samples) # t(s)
        s_interp=t_of_s.(ts)
        s_clamp=clamp.(s_true,first(s_samples),last(s_samples))
        t_interp=s_of_t.(s_clamp)
        err_s=maximum(abs.(s_interp.-s_true))
        err_t=maximum(abs.(t_interp.-ts))
        running=max(err_s,err_t)
        push!(errors_s,err_s);push!(errors_t,err_t)
        verbose && println("iter $i: n=$n_samples  Δs=$err_s  Δt=$err_t")
        running<=target_prec && break
        running>=best && (verbose && @warn "No further improvement (best=$best, now=$running).";break)
        best=running;n_samples*=2
    end

    return verbose ? (t_of_s,s_of_t,errors_s,errors_t) : (t_of_s,s_of_t)
end

# Possible are :chebyshev, :cubic_spline
function construct_arc_length_interpolation(::Type{T},crv::C;method::Symbol=:chebyshev,n_samples=100,q=3.0,p::Int=8,quad_rtol=1e-10,speed_ratio_max=3.0,init_panels::Int=8,nprobe::Int=9,tol_newton=1e-11,verbose=false) where {BC,T<:Real,C<:AbsPolarCurve{BC}}
    if method==:chebyshev
        return _construct_arc_length_interpolation_CHEBYSHEV(T,crv;q=q,p=p,quad_rtol=quad_rtol,speed_ratio_max=speed_ratio_max,init_panels=init_panels,nprobe=nprobe,tol_newton=tol_newton)
    elseif method==:cubic_spline
        return _construct_arc_length_interpolation_CUBIC_SPLINE(crv;target_prec=100*quad_rtol,verbose=verbose) # this one is usually always less precise, but might be faster in some situtations
    else
        error("Unsupported interpolation method: $method")
    end
end