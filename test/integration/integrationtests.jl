# Integration tier: solver-free geometry pipelines validated against analytic
# or cross-method references, grouped by computed quantity.

@group "full boundary" begin
    # Reconstructed full boundary: closure, orientation, length, area and node permutations.
    include("boundary/fullboundary_closure.jl")
end

@group "symmetry orbits" begin
    # Node orbits of the symmetry group against geometric point maps.
    include("symmetry/orbit_maps_billiards.jl")
end

@group "Poincare-Birkhoff" begin
    # Poincare-Birkhoff coordinate s against sector windows and full-boundary arc length.
    include("poincare_birkhoff/pb_unfolding.jl")
end
