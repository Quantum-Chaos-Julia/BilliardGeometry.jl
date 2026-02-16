
function is_inside(billiard::B, pt::SVector{2,T}) where {B<:AbsBilliard, T<:Real}
    d = [all(is_inside(domain, pt)) for domain in get_all_domains(billiard)]
    return any(d)
end

function is_inside(billiard::B, pts::AbstractArray) where {B<:AbsBilliard}
    return [is_inside(billiard, pt) for pt in pts]
end


function is_inside(domain::D, pt::SVector{2,T}) where {D<:AbsDomain, T<:Real}
    d = [is_inside(crv, pt) for crv in domain.boundary]
    return d
end

function is_inside(domain::D, pts::AbstractArray) where {D<:AbsDomain}
    d = [is_inside(crv, pts) for crv in domain.boundary]
    return  reduce(hcat,d)
end


#check if points inside for general curves
function is_inside(curve::C, pt::SVector{2,T}) where {C<:AbsCurve, T<:Real}
    return domain_fun(curve, pt) .< zero(eltype(pt)) 
end

function is_inside(curve::C, pts::AbstractArray) where {C<:AbsCurve}
    let
    d = domain_fun(curve, pts)
    return d .< zero(eltype(pts[1])) 
    end
end

#gradient of domain_function gives normal direcrion
function domain_gradient_vector(curve::C, pt::SVector{2,T}) where {C<:AbsCurve, T<:Real}
    f(r) = domain_fun(curve, r)
    g = ForwardDiff.gradient(f, pt)
    return g
end

function domain_gradient_vector(curve::C, pts::AbstractArray) where {C<:AbsCurve}
    f(r) = domain_fun(curve, r)
    gs = [ForwardDiff.gradient(f, pt) for pt in pts]
    return gs
end

function curvature(crv::L,t) where {L<:AbsCurve}
    let 
        r(t)=curve(crv,t)
        dr(t)=ForwardDiff.derivative(r,t)
        ddr(t)=ForwardDiff.derivative(dr,t)
        der=dr(t)
        der2=ddr(t)
        norm=hypot(der[1],der[2])^3
        kap=der[1]*der2[2]-der[2]*der2[1]
        return kap/norm
    end
end
