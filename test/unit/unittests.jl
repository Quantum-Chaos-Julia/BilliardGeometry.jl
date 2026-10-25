# Unit tier: one explicit @testset per function or type, mirroring src/.
# legacy_tests.jl holds the testsets that are still to be split into unit and integration parts.

include("geometry/segments/linesegment_tests.jl")
include("geometry/segments/circlesegment_tests.jl")
include("geometry/segments/polarcurves_tests.jl")
include("geometry/curvederivatives_tests.jl")
include("geometry/boundarycomponents_tests.jl")
include("geometry/symmetry_tests.jl")
include("geometry/boundarytypes_tests.jl")
include("quadrature/kressgrading_tests.jl")
include("quadrature/samplers_tests.jl")

include("legacy_tests.jl")
