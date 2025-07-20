struct SpecularReflection <:AbsBoundaryCondition end

struct Transparent <:AbsBoundaryCondition
    next_id::Int64
end

struct PeriodicX <:AbsBoundaryCondition
    next_id::Int64
end

struct ReflectionSymmetry{S} <:AbsBoundaryCondition where S<:AbsReflection
    symmetry::S
    N_sectors::Int64
end

function get_boundary_curves(domain::D) where D<:AbsSimpleDomain
    is_outer(crv) = typeof(crv.bc) <: SpecularReflection
    boundary = filter(is_outer, domain.boundary)
    return connect_curves(boundary)
end

function get_boundary_curves(composite_domain::D) where D<:AbsCompositeDomain
    boundary = Vector{AbsCurve}()
    for domain in composite_domain.subdomains
        subboundary = get_boundary_curves(domain)
        append!(boundary,subboundary)
    end
    return connect_curves(boundary)
end

function get_boundary_curves(billiard::B) where B<:AbsBilliard
    return get_boundary_curves(billiard.fundamental_domain)
end