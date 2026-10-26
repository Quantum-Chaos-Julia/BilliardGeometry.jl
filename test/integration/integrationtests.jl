# Integration tier: solver-free geometry pipelines validated against analytic
# or cross-method references, grouped by computed quantity.

# Reconstructed full boundary: closure, orientation, length, area and node permutations.
include("boundary/fullboundary_closure.jl")

# Node orbits of the symmetry group against geometric point maps.
include("symmetry/orbit_maps_billiards.jl")

# Poincare-Birkhoff coordinate s against sector windows and full-boundary arc length.
include("poincare_birkhoff/pb_unfolding.jl")

