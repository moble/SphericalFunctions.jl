# This file exists only for the things that insist on calling `Pkg.test`: the registry, CI
# actions like `julia-actions/julia-runtest`, downstream integration testing, and `] test`
# habits.  It is not the day-to-day runner — see docs/src/60-development/01-index.md.
#
# Day to day, use `juliati` (the TestItemApp command-line runner), the `julia` MCP server, or
# an editor's test-item support.  Those can filter by name, file or tag, run items in
# parallel, keep worker processes warm between runs, and report per-item results; none of
# that is worth reimplementing here.

using TestItemRunner

# `Pkg.test(...; test_args)` arrives here as `ARGS`.  The only filtering this shim supports is
# by tag — written as `:sometag` — because that is all the CI workflows use.
const CI = get(ENV, "CI", "false") == "true"
const requested_tags = Symbol[Symbol(a[2:end]) for a ∈ ARGS if startswith(a, ":")]

function testfilter(testitem)
    (; tags) = testitem
    if !isempty(requested_tags)
        # An explicit request wins, including over the `:skipci` rule below.
        return any(∈(tags), requested_tags)
    end
    # Items tagged `:skipci` need something CI does not have (a Python environment, say).
    !(CI && :skipci ∈ tags)
end

@run_package_tests verbose = true filter = testfilter

# Including the test files is not needed for discovery — `@run_package_tests` finds them on
# its own — but it makes `Pkg.test` parse each one, so a syntax error shows up here rather
# than as a silently missing test item.
include("aqua.jl")
include("complex_powers.jl")
include("haxis.jl")
include("hwedge.jl")
include("operators.jl")
include("strided.jl")
include("weights.jl")
include("mode_weights/indexing.jl")
include("mode_weights/containers.jl")
include("mode_weights/operations.jl")
include("mode_weights/mode_weights.jl")
include("sYlm/real_harmonics.jl")
include("sYlm/sYlm.jl")
include("ssht/map2salm.jl")
include("ssht/ssht.jl")
include("utilities/combinatorics.jl")
include("utilities/encoder.jl")
include("utilities/explicit_operators.jl")
include("utilities/explicit_wigner_matrices.jl")
include("utilities/naive_factorial.jl")
include("utilities/utilities.jl")
include("wigner/H_calculator.jl")
include("wigner/calculators.jl")
include("wigner/half_integer.jl")
include("wigner/half_integer_oracle.jl")
include("wigner/iteration.jl")
include("wigner/properties.jl")
include("wigner/recurrence.jl")
include("wigner/robustness.jl")
