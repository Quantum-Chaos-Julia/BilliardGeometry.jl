# Unit tier: one explicit @testset per function or type, mirroring src/.
# Geometric checks on real billiards (closure, perimeter, orbit maps, Poincare-Birkhoff unfolding) live in the integration tier.

include("holes_fixtures.jl")   # shared hole/multiply connected fixtures (no tests), must stay before the files that use them

@group "geometry: core" begin
    include("geometry/utils_tests.jl")
    include("geometry/arclength_tests.jl")
end

@group "geometry: segments" begin
    include("geometry/segments/linesegment_tests.jl")
    include("geometry/segments/circlesegment_tests.jl")
    include("geometry/segments/polarcurves_tests.jl")
    include("geometry/segments/compositecurves_tests.jl")
end

@group "geometry: curves and symmetry" begin
    include("geometry/curvederivatives_tests.jl")
    include("geometry/inversions_tests.jl")
    include("geometry/boundarycomponents_tests.jl")
    include("geometry/symmetry_tests.jl")
    include("geometry/symmetryregistry_tests.jl")
end

@group "geometry: domains" begin
    include("geometry/domains/polygons_tests.jl")
    include("geometry/domains/circular_tests.jl")
    include("geometry/domains/compositedomains_tests.jl")
    include("geometry/domains/holes_tests.jl")
    include("geometry/domains/multiplyconnecteddomains_tests.jl")
end

@group "geometry: billiards" begin
    include("geometry/billiards/stadium_tests.jl")
    include("geometry/billiards/mushroom_tests.jl")
    include("geometry/billiards/triangle_tests.jl")
    include("geometry/billiards/polar_tests.jl")
    include("geometry/billiards/limacon_tests.jl")
    include("geometry/billiards/circle_tests.jl")
    include("geometry/billiards/ellipse_tests.jl")
    include("geometry/billiards/rectangle_tests.jl")
    include("geometry/billiards/polygon_tests.jl")
    include("geometry/billiards/star_tests.jl")
    include("geometry/billiards/c3_tests.jl")
    include("geometry/billiards/prosen_tests.jl")
    include("geometry/billiards/annular_tests.jl")
    include("geometry/billiards/sinai_tests.jl")
    include("geometry/billiards/square_in_square_tests.jl")
end

@group "geometry: boundaries and orbits" begin
    include("geometry/boundarytypes_tests.jl")
    include("geometry/geometry_tests.jl")
    include("geometry/area_tests.jl")
    include("geometry/symmetryorbits_tests.jl")
    include("geometry/fullboundary_tests.jl")
    include("geometry/poincarebirkhoff_tests.jl")
end

@group "quadrature" begin
    include("quadrature/kressgrading_tests.jl")
    include("quadrature/samplers_tests.jl")
end

@group "quality" begin
    include("quality/aqua_tests.jl")
end
