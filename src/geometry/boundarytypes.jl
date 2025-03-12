struct SpecularReflection <:AbsBoundaryCondition end

struct Transparent <:AbsBoundaryCondition
    domain_exit_idx::Int64
end

struct ReflectionSymmetry{S} <:AbsBoundaryCondition where S<:AbsReflection
    symmetry::S
end