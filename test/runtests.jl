using Test
using BilliardGeometry

const BG_TEST_GROUPS = ("unit", "integration", "all")
const BG_TEST_GROUP = get(ENV, "BG_TEST_GROUP", "all")
BG_TEST_GROUP in BG_TEST_GROUPS ||
    error("Unknown BG_TEST_GROUP=\"$BG_TEST_GROUP\"; expected one of $(join(BG_TEST_GROUPS, ", "))")

if BG_TEST_GROUP in ("unit", "all")
    include("unit/unittests.jl")
end

if BG_TEST_GROUP in ("integration", "all")
    include("integration/common.jl")
    include("integration/integrationtests.jl")
end
