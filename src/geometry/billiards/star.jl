"""
    StarBilliard{T} <: AbsBilliard

Smooth star-shaped billiard `r(φ) = R + a*cos(nφ)` with `n`-fold rotational
symmetry. The fundamental domain is the sector `0 ≤ φ ≤ 2π/n` bounded by the
boundary arc and two [`SymmetryWall`](@ref) radii. Only the `n-1` nontrivial
rotations of [`Cn_symmetry`](@ref) are registered, no reflections, to avoid
duplicated curves in [`full_boundary`](@ref).

## Attributes
* `fundamental_domain::SimpleDomain{T}`: The sector.
* `symmetries::SymmetryRegistry`: Rotations by `2πm/n`, `m = 1, ..., n-1`, with `sym_id = m`.
"""
struct StarBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::SimpleDomain{T}
    symmetries::SymmetryRegistry
end

"""
    StarBilliard(R, a, n; center = (0,0)) → billiard::StarBilliard

Construct a [`StarBilliard`](@ref).

## Arguments
* `R`: Mean radius, must satisfy `R > |a|` so that `r(φ) > 0`.
* `a`: Amplitude of the modulation.
* `n`: Order of the rotational symmetry, `n ≥ 2`.

## Keyword arguments
* `center` : Must be the origin, since [`NFoldRotation`](@ref) rotates about the origin.
"""
function StarBilliard(R::T, a::T, n::Int; center=SVector{2,T}(zero(T),zero(T))) where T<:Real
    n >= 2 || throw(ArgumentError("n must satisfy n>=2; received n=$n"))
    R > abs(a) || throw(ArgumentError("Require R>|a| so r(phi)>0; received R=$R, a=$a"))
    c = SVector{2,T}(center)
    (iszero(c[1]) && iszero(c[2])) || throw(ArgumentError("Cn symmetry requires center == (0,0) (NFoldRotation always rotates about the coordinate origin); received center=$c"))
    coef = SVector{2*n,T}(ntuple(i -> i==2*n ? a : zero(T), 2*n))
    bc = SpecularReflection()
    #only the n-1 nontrivial rotations are needed to tile the fundamental
    #2*pi/n sector back into the full 2*pi boundary; adding reflections on
    #top (even though r(phi) may also be reflection-symmetric for even n)
    #would duplicate already-covered curves in `full_boundary`.
    symmetries = Cn_symmetry(n; T=T) # sym_id 1 = rotate by +2π/n (m=1), ..., sym_id n-1 = rotate by +2π(n-1)/n (m=n-1)
    arc = FourierCoeffPolarSegment(coef; R=R, arc_angle=T(2*pi/n), shift_angle=zero(T), center=c, bc=bc, domain_id=1, segment_id=1)
    p0 = curve(arc, zero(T))
    p1 = curve(arc, one(T))
    # wall1 (t=1 end) borders the copy reached by rotating forward by +2π/n (sym_id 1);
    # wall0 (t=0 end) borders the copy reached by rotating backward, i.e. by +2π(n-1)/n (sym_id n-1).
    wall1 = LineSegment(p1, c; bc=SymmetryWall(1,2), domain_id=1, segment_id=2)
    wall0 = LineSegment(c, p0; bc=SymmetryWall(n-1,2), domain_id=1, segment_id=3)
    fundamental_boundary = AbsCurve[arc, wall1, wall0]
    vertices = SVector{2,T}[p0, p1, c]
    fundamental_domain = SimpleDomain{T}(fundamental_boundary, vertices, 1)
    return StarBilliard{T}(fundamental_domain, symmetries)
end
