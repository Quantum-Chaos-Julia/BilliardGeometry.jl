"""
    invert_curve(crv, pt; method=:auto, t_init=0.5, rtol=sqrt(eps(T)), atol=1e-15, maxevals=50, n_grid=20)

Find parameter t such that curve(crv, t) ≈ pt.

# Methods
- `:auto` - Tries Roots.jl first, falls back to Optim.jl if needed
- `:roots` - Uses Roots.jl Newton method
- `:optim` - Uses Optim.jl BFGS (more robust)
- `:grid_optim` - Grid search + Optim refinement (most robust)
"""
function invert_curve(crv::C, pt::SVector{2,T}; 
                      method=:auto,
                      t_init=0.5, 
                      rtol=sqrt(eps(T)), 
                      atol=1e-15, 
                      maxevals=50, 
                      n_grid=20) where {C<:AbsCurve, T<:Real}
    
    if method == :auto
        # Try fast method first, fallback to robust
        try
            return invert_curve_roots(crv, pt, t_init, rtol, atol, maxevals)
        catch e
            @warn "Roots.jl failed, using Optim.jl" exception=e
            return invert_curve_optim(crv, pt, t_init, rtol, atol, maxevals)
        end
    elseif method == :roots
        return invert_curve_roots(crv, pt, t_init, rtol, atol, maxevals)
    elseif method == :optim
        return invert_curve_optim(crv, pt, t_init, rtol, atol, maxevals)
    elseif method == :grid_optim
        return invert_curve_grid_search(crv, pt; rtol, atol, maxevals, n_grid)
    else
        error("Unknown method: $method. Use :auto, :roots, :optim, or :grid_optim")
    end
end

"""
Roots.jl based inversion (fastest for well-behaved curves).
"""
function invert_curve_roots(crv::C, pt::SVector{2,T}, t_init, rtol, atol, maxevals) where {C<:AbsCurve, T<:Real}
    # Define objective function: distance squared
    function objective(t)
        pos = curve(crv, t)
        return sum((pos .- pt).^2)
    end
    
    result = find_zero(objective, t_init, atol=atol, rtol=rtol, maxevals=maxevals)
    # Wrap to [0, 1]
    t = mod(result, 1.0)
    return t
end

"""
Optim.jl based inversion (more robust, handles difficult cases better).
Uses BFGS with automatic differentiation.
"""
function invert_curve_optim(crv::C, pt::SVector{2,T}, t_init, rtol, atol, maxevals) where {C<:AbsCurve, T<:Real}
    # Objective: minimize distance squared
    function objective(t_vec)
        t = mod(t_vec[1], 1.0)  # Keep in [0, 1]
        pos = curve(crv, t)
        return sum((pos .- pt).^2)
    end
    
    # Gradient using ForwardDiff
    function gradient!(G, t_vec)
        t = mod(t_vec[1], 1.0)
        pos = curve(crv, t)
        vel = ForwardDiff.derivative(s -> curve(crv, s), t)
        residual = pos .- pt
        G[1] = 2 * dot(vel, residual)
    end
    
    # Use BFGS with gradient
    result = optimize(objective, gradient!, [t_init], BFGS(),
                     Optim.Options(
                         g_tol=atol,
                         f_abstol=atol,
                         f_reltol=atol,
                         iterations=maxevals,
                         show_trace=false
                     ))
    
    t = mod(result.minimizer[1], 1.0)
    return t
end

"""
Grid search followed by Optim refinement (most robust).
"""
function invert_curve_grid_search(crv::C, pt::SVector{2,T}; 
                                  rtol=sqrt(eps(T)), 
                                  atol=1e-15, 
                                  maxevals=50, 
                                  n_grid=20) where {C<:AbsCurve, T<:Real}
    # Coarse grid search
    t_samples = range(0.0, 1.0, length=n_grid+1)[1:end-1]
    min_dist = Inf
    best_t = 0.5
    
    for t in t_samples
        pos = curve(crv, t)
        dist_sq = sum((pos .- pt).^2)
        
        if dist_sq < min_dist
            min_dist = dist_sq
            best_t = t
        end
    end
    
    # Check if we're already very close
    if min_dist < atol^2
        return best_t
    end
    
    # Refine with Optim (more robust than Roots for difficult cases)
    return invert_curve_optim(crv, pt, best_t, rtol, atol, maxevals)
end

"""
Advanced: Use Optim's NelderMead (derivative-free) for very difficult curves.
Slower but works when derivatives are problematic.
"""
function invert_curve_nelder_mead(crv::C, pt::SVector{2,T}, t_init, atol, maxevals) where {C<:AbsCurve, T<:Real}
    function objective(t_vec)
        t = mod(t_vec[1], 1.0)
        pos = curve(crv, t)
        return sum((pos .- pt).^2)
    end
    
    result = optimize(objective, [t_init], NelderMead(),
                     Optim.Options(
                         f_tol=atol,
                         iterations=maxevals
                     ))
    
    t = mod(result.minimizer[1], 1.0)
    return t
end