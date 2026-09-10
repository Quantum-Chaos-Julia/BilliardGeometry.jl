"""
    boundary_coords(crv::C,sampler::S,N::Int) where {C<:AbsCurve,S<:AbsSampler}

Sample a boundary curve using `sampler`.

## Arguments
* `crv::C`: Boundary curve.
* `sampler::S`: Parameter-space sampler.
* `N::Int`: Number of boundary nodes.

## Returns
* `xy`: Boundary coordinates.
* `normal`: Outward unit normals.
* `s`: Arc-length coordinates.
* `ds`: Physical quadrature weights.
"""
function boundary_coords(crv::C,sampler::S,N::Int) where {C<:AbsCurve,S<:AbsSampler}
    t,dt=sample_points(sampler,N)
    return boundary_coords(crv,t,dt)
end

"""
    boundary_coords(crv::C,t,dt) where {C<:AbsCurve}

Evaluate the boundary geometry at parameter nodes `t`.

## Arguments
* `crv::C`: Boundary curve.
* `t`: Parameter nodes.
* `dt`: Parameter-space quadrature weights.

## Returns
* `xy`: Boundary coordinates.
* `normal`: Outward unit normals.
* `s`: Arc-length coordinates.
* `ds`: Physical quadrature weights.
"""
function boundary_coords(crv::C,t,dt) where {C<:AbsCurve}
    xy=curve(crv,t)
    tang=tangent(crv,t)
    normal=domain_gradient_vector(crv,xy)
    normal./=norm.(normal)
    s=arc_length(crv,t)
    ds=norm.(tang).*dt
    return xy,normal,s,ds
end