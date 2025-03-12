
export is_overlaping, is_connected, is_closed, angle, circle_center, angle

function is_overlaping(pt1, pt2)
    x = isapprox(pt1[1],pt2[1])
    y = isapprox(pt1[2],pt2[2])
    return x && y
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

function is_closed(boundary; start_intial = true)
    return all(is_connected(boundary; start_intial))
end 

angle(a, b) = atan(norm(cross(a,b)),dot(a,b))

function circle_center(pt0, pt1, r)
    type = eltype(pt0)
    pt0 = SVector{2, type}(pt0)
    pt1 = SVector{2, type}(pt1)

    mid = (pt0 + pt1) ./ 2
    l = norm(pt1 - pt0)

    # Exit condition if circle is not possible
    if r < l/2
        println("Error: Radius is smaller than half the distance between points. Circle not possible.")
        return nothing
    end
    # The circle center will be along the perpendicular bisector of pt0 and pt1
    # Compute the slope of this bisector
    if pt0[1] == pt1[1]
        p = SVector(1.0,0.0)
    elseif pt0[2] == pt1[2]
        p = SVector(0.0,1.0) 
    else
        slope = (pt1[2] - pt0[2]) / (pt1[1] - pt0[1])
        p = normalize(SVector{2,type}(1.0, -1.0 / slope))
    end
    # Find the distance between M and the circle center C, which are located on the line with slope m_PM
    d = hypot(r,l/2)
    # The sign can be + or -, leading to two potential circle centers
    C1 = mid .+ d .* p
    C2 = mid .- d .* p
    return C1, C2
end