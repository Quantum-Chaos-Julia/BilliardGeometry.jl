
export is_overlaping, is_connected, is_closed, angle, connect_curves

function is_overlaping(pt1, pt2)
    x = isapprox(pt1[1],pt2[1])
    y = isapprox(pt1[2],pt2[2])
    return x && y
end

function connect_curves(curves)
    if is_closed(curves)
        return curves
    end
    connected_curves = Vector{AbsCurve}()
    remaining_curves = copy(curves)
    push!(connected_curves, popfirst!(remaining_curves))
    for i in 2:length(curves)
        end_pt = curve(connected_curves[end],1.0)
        N_remaining = length(remaining_curves)
        for j in 1:N_remaining
            start_pt = curve(remaining_curves[j],0.0)
            if is_overlaping(end_pt, start_pt)
                push!(connected_curves,splice!(remaining_curves,j))
                break
            end
        end
    end
    if length(remaining_curves) > 0
        second_section = connect_curves(remaining_curves)
        start_first, end_first = curve(connected_curves[1],0.0), curve(connected_curves[end],1.0)
        start_second, end_second = curve(second_section[1],0.0), curve(second_section[end],1.0)
        if is_overlaping(end_first, start_second)
            append!(connected_curves,second_section)
        elseif is_overlaping(end_second, start_first)
            prepend!(connected_curves,second_section)
        end
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

