#Prosen billiard r(phi) = 1+a*cos(4*phi); fourfold symmetric shape reduced via
#reflection walls (D2 quadrant fundamental domain, matching the SymmetryWall
#wall boundary conditions below) rather than a rotational sector.
struct ProsenBilliard{T} <: AbsBilliard where T<:Real
    fundamental_domain::SimpleDomain{T}
    symmetries::SymmetryRegistry
end

function ProsenBilliard(a::T; center=SVector{2,T}(zero(T),zero(T))) where T<:Real
    c = SVector{2,T}(center)
    (iszero(c[1]) && iszero(c[2])) || throw(ArgumentError("D2 symmetry requires center == (0,0); received center=$c"))
    coef = SVector{8,T}(zero(T),zero(T),zero(T),zero(T),zero(T),zero(T),zero(T),a)
    bc = SpecularReflection()
    symmetries = D2_symmetry() # sym_id 1=YAxisReflection, 2=XYAxisReflection, 3=XAxisReflection
    arc = FourierCoeffPolarSegment(coef; R=one(T), arc_angle=T(pi/2), shift_angle=zero(T), center=c, bc=bc, domain_id=1, segment_id=1)
    p0 = curve(arc, zero(T))
    p1 = curve(arc, one(T))
    wall1 = LineSegment(c, p0; bc=SymmetryWall(3,2), domain_id=1, segment_id=2)
    wall2 = LineSegment(p1, c; bc=SymmetryWall(1,2), domain_id=1, segment_id=3)
    fundamental_boundary = AbsCurve[arc, wall2, wall1]
    vertices = SVector{2,T}[p0, p1, c]
    fundamental_domain = SimpleDomain{T}(fundamental_boundary, vertices, 1)
    return ProsenBilliard{T}(fundamental_domain, symmetries)
end
