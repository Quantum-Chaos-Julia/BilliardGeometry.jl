
"""
    is_overlaping(pt1, pt2) → flag::Bool

Return `true` when the points `pt1` and `pt2` coincide within the absolute
tolerance `sqrt(eps(T))` in both coordinates, where `T` is the floating point
type promoted from the point entries.

## Arguments
* `pt1`, `pt2`: Planar points, indexable by `1` and `2`.

## Returns
* `flag`: Whether the two points are the same up to tolerance.
"""
function is_overlaping(pt1, pt2)
    T = float(promote_type(typeof(pt1[1]), typeof(pt2[1]), typeof(pt1[2]), typeof(pt2[2])))
    atol = sqrt(eps(T))
    x = isapprox(pt1[1], pt2[1]; atol=atol)
    y = isapprox(pt1[2], pt2[2]; atol=atol)
    return x && y
end

#TODO: Let's test this function.
"""
    connect_curves(curves) → connected::Vector{AbsCurve}

Reorder `curves` so that the end point of every curve coincides with the start
point of the next one.

## Description
If `curves` is already closed in the sense of [`is_closed`](@ref) it is returned
unchanged. Otherwise curves are first chained greedily, each chain extended by a
remaining curve whose start point [`is_overlaping`](@ref) with the chain end,
and chains are then merged by the same rule. Curves are never reversed. Chains
that cannot be joined are concatenated in the order they were found, so the
result is always a permutation of the input.

## Arguments
* `curves`: Vector of [`AbsCurve`](@ref).

## Returns
* `connected`: The curves in chain order.
"""
function connect_curves(curves)
    if is_closed(curves)
        return curves
    end
    remaining_curves = copy(curves)
    chains = Vector{Vector{AbsCurve}}()
    while !isempty(remaining_curves)
        chain = Vector{AbsCurve}()
        push!(chain, popfirst!(remaining_curves))
        extended = true
        while extended
            extended = false
            end_pt = curve(chain[end],1.0)
            for j in eachindex(remaining_curves)
                start_pt = curve(remaining_curves[j],0.0)
                if is_overlaping(end_pt, start_pt)
                    push!(chain, splice!(remaining_curves,j))
                    extended = true
                    break
                end
            end
        end
        push!(chains, chain)
    end
    merged = true
    while merged && length(chains) > 1
        merged = false
        for i in eachindex(chains)
            end_pt = curve(chains[i][end], 1.0)
            for j in eachindex(chains)
                i == j && continue
                start_pt = curve(chains[j][1], 0.0)
                if is_overlaping(end_pt, start_pt)
                    chains[i] = vcat(chains[i], chains[j])
                    deleteat!(chains, j)
                    merged = true
                    break
                end
            end
            merged && break
        end
    end
    connected_curves = Vector{AbsCurve}()
    for chain in chains
        append!(connected_curves, chain)
    end
    return connected_curves
end

# Splits a chain-ordered curve list (as returned by `connect_curves`) into
# closed loops: a loop ends at the first curve whose end point coincides with
# the start point of the loop's first curve. Throws if a chain does not close.
function _split_closed_loops(curves)
    loops = Vector{Vector{AbsCurve}}()
    n = length(curves)
    i = 1
    while i <= n
        start_pt = curve(curves[i], zero(curves[i].length))
        j = i
        while j <= n && !is_overlaping(curve(curves[j], one(curves[j].length)), start_pt)
            j += 1
        end
        j > n && throw(ArgumentError("boundary curves do not form closed loops: the loop starting at curve $i does not close"))
        push!(loops, AbsCurve[curves[i:j]...])
        i = j + 1
    end
    return loops
end

"""
    is_connected(boundary; start_initial = true) → flags::Vector{Bool}

Test whether consecutive curves of `boundary` join end to start, treating the
list as periodic.

## Arguments
* `boundary`: Vector of [`AbsCurve`](@ref) in traversal order.

## Keyword arguments
* `start_initial::Bool = true` : If `true`, `flags[i]` tests whether the end of curve `i-1` meets the start of curve `i` (so `flags[1]` compares the last curve with the first). If `false`, `flags[i]` tests the end of curve `i` against the start of curve `i+1`.

## Returns
* `flags`: One entry per curve, see [`is_overlaping`](@ref).
"""
function is_connected(boundary; start_initial = true)
    start_pts = CircularArray([curve(crv,0.0) for crv in boundary])
    end_pts = CircularArray([curve(crv,1.0) for crv in boundary])    
    if start_initial
        test = [is_overlaping(end_pts[i-1], start_pts[i]) for i in 1:length(start_pts)]
    else
        test = [is_overlaping(end_pts[i], start_pts[i+1]) for i in 1:length(start_pts)]
    end
    return test
end

"""
    is_closed(boundary; check_periodic = false, start_initial = true) → flag::Bool

Return `true` when all consecutive curves of `boundary` are joined, see
[`is_connected`](@ref).

## Arguments
* `boundary`: Vector of [`AbsCurve`](@ref) in traversal order.

## Keyword arguments
* `check_periodic::Bool = false` : If `true`, the junction between the last and the first curve is checked as well. If `false`, only the junctions inside the list are checked, so an open chain is accepted.
* `start_initial::Bool = true` : Passed to [`is_connected`](@ref).

## Returns
* `flag`: Whether the curves form a connected chain (or closed loop if `check_periodic`).
"""
function is_closed(boundary; check_periodic=false,  start_initial = true)
    if check_periodic
        return all(is_connected(boundary; start_initial=start_initial))
    end

    if start_initial
        return all(is_connected(boundary; start_initial=start_initial)[2:end])
    else
        return all(is_connected(boundary; start_initial=start_initial)[1:end-1])
    end
end 


angle(a, b) = atan(norm(cross(a,b)),dot(a,b))

function signed_angle(a,b)
    return atan(a[1]*b[2] - a[2]*b[1], dot(a, b))
end

"""
    find_unique_elements(vector) → unique_values::Vector

Remove entries of a sorted `vector` that coincide, within `sqrt(eps)` of the
element type, with the following entry, and return the remaining distinct
values. Used to merge numerically duplicated boundary parameters such as
sector end points.

## Arguments
* `vector`: Real-valued vector, sorted so that near-duplicates are adjacent.

## Returns
* `unique_values`: The vector without near-duplicates.
"""
function find_unique_elements(vector)    
    atol = sqrt(eps(float(eltype(vector))))
    u = trues(length(vector))
    for i in 1:(length(vector)-1)
        if isapprox(vector[i], vector[i+1]; atol=atol)
            u[i] = false
        end
    end
    return unique(vector[u])
end

"""
    point_curve_parameter(crv::AbsCurve, pt) → t

Find the parameter `t ∈ [0,1]` for which `curve(crv, t) ≈ pt`.

## Description
The zeros of `y(t) - pt[2]` are located with `Roots.find_zeros`. If there is
more than one, the one that also solves the `x` equation is returned. Unlike
[`invert_curve`](@ref) it is a root-finding search on the whole interval and
does not require an initial guess. If `pt` is not on the curve the result is
undefined: `nothing` is returned when several `y` roots have no matching `x`
root, and an error is thrown when there is no `y` root at all.

## Arguments
* `crv`: Curve to be inverted.
* `pt`: Point on the curve.

## Returns
* `t`: Curve parameter of `pt`.
"""
function point_curve_parameter(crv::C, pt) where {C<:AbsCurve}
    inv_x(theta) = curve(crv,theta)[1] - pt[1] 
    inv_y(theta) = curve(crv,theta)[2] - pt[2]
    roots_y = find_zeros(inv_y, (0.0, 1.0))

    if length(roots_y)>1
        roots_x = find_zeros(inv_x, (0.0, 1.0))
        for t_y in roots_y
            for t_x in roots_x
                if isapprox(t_y, t_x)
                    return t_y
                end
            end
        end
    else
        return roots_y[1]
    end
end

