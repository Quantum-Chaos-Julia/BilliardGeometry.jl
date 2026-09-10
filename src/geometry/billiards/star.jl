#smooth star-shaped billiard r(phi) = R + a*cos(n*phi), n-fold rotational symmetry
struct StarBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::SimpleDomain{T}
    symmetries::Vector{AbsSymmetry}
end

function StarBilliard(R::T, a::T, n::Int; center=SVector{2,T}(zero(T),zero(T))) where T<:Real
    n >= 2 || throw(ArgumentError("n must satisfy n>=2; received n=$n"))
    R > abs(a) || throw(ArgumentError("Require R>|a| so r(phi)>0; received R=$R, a=$a"))
    c = SVector{2,T}(center)
    coef = SVector{2*n,T}(ntuple(i -> i==2*n ? a : zero(T), 2*n))
    bc = SpecularReflection()
    arc = FourierCoeffPolarSegment(coef; R=R, arc_angle=T(2*pi/n), shift_angle=zero(T), center=c, bc=bc, domain_id=1, segment_id=1)
    p0 = curve(arc, zero(T))
    p1 = curve(arc, one(T))
    wall1 = LineSegment(p1, c; bc=QuantumSolverIgnore(), domain_id=1, segment_id=2)
    wall0 = LineSegment(c, p0; bc=QuantumSolverIgnore(), domain_id=1, segment_id=3)
    fundamental_boundary = AbsCurve[arc, wall1, wall0]
    vertices = SVector{2,T}[p0, p1, c]
    fundamental_domain = SimpleDomain{T}(fundamental_boundary, vertices, 1)
    #only the n-1 nontrivial rotations are needed to tile the fundamental
    #2*pi/n sector back into the full 2*pi boundary; adding reflections on
    #top (even though r(phi) may also be reflection-symmetric for even n)
    #would duplicate already-covered curves in `full_boundary`.
    symmetries = AbsSymmetry[Cn_symmetry(n)...]
    return StarBilliard{T}(fundamental_domain, symmetries)
end
