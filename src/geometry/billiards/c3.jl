#smooth threefold-rotationally symmetric billiard r(phi) = scale/2*(1+a*(cos(3*phi)-sin(6*phi)))
struct C3Billiard{T} <: AbsBilliard where T<:Real
    fundamental_domain::SimpleDomain{T}
    symmetries::SymmetryRegistry
end

function C3Billiard(a::T; scale::T=one(T)) where T<:Real
    scale > zero(T) || throw(ArgumentError("scale must be positive; received scale=$scale"))
    c = SVector{2,T}(zero(T),zero(T))
    R = scale/T(2)
    amp = a*scale/T(2)
    coef = SVector{12,T}(ntuple(i -> i==6 ? amp : i==11 ? -amp : zero(T), 12))
    bc = SpecularReflection()
    symmetries = Cn_symmetry(3) # sym_id 1 = rotate by +2π/3 (m=1), sym_id 2 = rotate by +4π/3 (m=2)
    arc = FourierCoeffPolarSegment(coef; R=R, arc_angle=T(2*pi/3), shift_angle=zero(T), center=c, bc=bc, domain_id=1, segment_id=1)
    p0 = curve(arc, zero(T))
    p1 = curve(arc, one(T))
    # wall1 (t=1 end) borders the copy reached by rotating forward by +2π/3 (sym_id 1);
    # wall0 (t=0 end) borders the copy reached by rotating backward, i.e. by +4π/3 (sym_id 2).
    wall1 = LineSegment(p1, c; bc=SymmetryWall(1,2), domain_id=1, segment_id=2)
    wall0 = LineSegment(c, p0; bc=SymmetryWall(2,2), domain_id=1, segment_id=3)
    fundamental_boundary = AbsCurve[arc, wall1, wall0]
    vertices = SVector{2,T}[p0, p1, c]
    fundamental_domain = SimpleDomain{T}(fundamental_boundary, vertices, 1)
    return C3Billiard{T}(fundamental_domain, symmetries)
end
