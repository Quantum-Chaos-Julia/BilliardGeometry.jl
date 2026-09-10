
export is_overlaping, is_connected, is_closed, angle, connect_curves, find_unique_elements, point_curve_parameter

function is_overlaping(pt1, pt2)
    x = isapprox(pt1[1],pt2[1])
    y = isapprox(pt1[2],pt2[2])
    return x && y
end

function connect_curves(curves)
    if is_closed(curves)
        return curves
    end
    # Greedily grow one connected chain at a time from whatever curves remain,
    # instead of assuming every input curve belongs to a single ring. This
    # correctly preserves multiply connected geometries (e.g. an outer wall
    # plus a disjoint circular obstacle) as separate closed components,
    # rather than silently dropping every component after the first.
    remaining_curves = copy(curves)
    connected_curves = Vector{AbsCurve}()
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
        append!(connected_curves, chain)
    end
    return connected_curves
end

function is_connected(boundary; start_intial = true)
    start_pts = CircularArray([curve(crv,0.0) for crv in boundary])
    end_pts = CircularArray([curve(crv,1.0) for crv in boundary])    
    if start_intial
        test = [is_overlaping(end_pts[i-1], start_pts[i]) for i in 1:length(start_pts)]
    else
        test = [is_overlaping(end_pts[i], start_pts[i+1]) for i in 1:length(start_pts)]
    end
    return test
end

function is_closed(boundary; check_periodic=false,  start_intial = true)
    if check_periodic
        return all(is_connected(boundary; start_intial))
    end

    if start_intial
        return all(is_connected(boundary; start_intial)[2:end])
    else
        return all(is_connected(boundary; start_intial)[1:end-1])
    end
end 


angle(a, b) = atan(norm(cross(a,b)),dot(a,b))

function signed_angle(a,b)
    return atan(a[1]*b[2] - a[2]*b[1], dot(a, b))
end

function find_unique_elements(vector)    
    u = [true for p in vector]
    for i in 1:(length(vector)-1)
        if isapprox(vector[i],vector[i+1])
            vector[i] = false
        end
    end
    return unique(vector[u])
end


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

