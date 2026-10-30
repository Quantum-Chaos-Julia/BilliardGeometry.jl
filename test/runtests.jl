using Test
using BilliardGeometry

const BG_TEST_GROUPS = ("unit", "integration", "all")
const BG_TEST_GROUP = get(ENV, "BG_TEST_GROUP", "all")
BG_TEST_GROUP in BG_TEST_GROUPS ||
    error("Unknown BG_TEST_GROUP=\"$BG_TEST_GROUP\"; expected one of $(join(BG_TEST_GROUPS, ", "))")

# BG_TEST_FILTER: comma-separated, case-insensitive substrings; only @group blocks whose name
# matches one of them run (e.g. BG_TEST_FILTER="segments,domains"). Unset or empty runs all groups.
const BG_TEST_FILTER = filter(!isempty, strip.(lowercase.(split(get(ENV, "BG_TEST_FILTER", ""), ','))))
group_enabled(name::AbstractString) =
    isempty(BG_TEST_FILTER) || any(f -> occursin(f, lowercase(name)), BG_TEST_FILTER)

# Each @group is a top-level verbose testset, so its summary prints as soon as it finishes.
# A failing group is recorded and the next group still runs; the run errors at the end.
const FAILED_GROUPS = String[]
macro group(name, body)
    return esc(quote
        if group_enabled($name)
            try
                @testset verbose=true $name $body
            catch e
                e isa Test.TestSetException || rethrow()
                push!(FAILED_GROUPS, $name)
            end
        end
    end)
end

if BG_TEST_GROUP in ("unit", "all")
    include("unit/unittests.jl")
end

if BG_TEST_GROUP in ("integration", "all")
    include("integration/common.jl")
    include("integration/integrationtests.jl")
end

isempty(FAILED_GROUPS) || error("Failed test groups: ", join(FAILED_GROUPS, ", "))
