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
        # An explicit request wins, including over the rules below.
        return any(∈(tags), requested_tags)
    end
    # Items tagged `:python` build a Conda environment on first use, which nobody running
    # `Pkg.test` — a user, a downstream package's CI — should get without asking for it.
    :python ∈ tags && return false
    # Items tagged `:skipci` need something CI does not have.
    !(CI && :skipci ∈ tags)
end

@run_package_tests verbose = true filter = testfilter

# Including the test files is not needed for discovery — `@run_package_tests` finds them on
# its own — but it makes `Pkg.test` parse each one, so a syntax error shows up here rather
# than as a silently missing test item.  Every file under `test/` is included, so that a new
# one cannot be forgotten.
for (root, _, files) ∈ walkdir(@__DIR__), file ∈ sort(files)
    path = joinpath(root, file)
    endswith(file, ".jl") && path != @__FILE__ && include(path)
end
