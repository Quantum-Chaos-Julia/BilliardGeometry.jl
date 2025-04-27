struct SpecularReflection <:AbsBoundaryCondition end

struct Transparent <:AbsBoundaryCondition
    next_id::Int64
end

struct PeriodicX <:AbsBoundaryCondition
    next_id::Int64
end

struct ReflectionSymmetry{S} <:AbsBoundaryCondition where S<:AbsReflection
    symmetry::S
end

function get_boundary_curves(domain::D) where D<:BilliardGeometry.AbsSimpleDomain
    is_outer(crv) = typeof(crv.bc) <: BilliardGeometry.SpecularReflection
    return filter(is_outer, domain.boundary)
end
