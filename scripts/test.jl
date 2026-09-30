# Call this from the top-level directory as
#   julia -t auto scripts/test.jl
# or to run with coverage as
#   julia -t auto scripts/test.jl --coverage
# See docs/src/60-development/01-index.md for more information.

import Dates
println("Running tests starting at ", Dates.format(Dates.now(), "HH:MM:SS"), ".")

using Pkg
cd((@__DIR__) * "/..")
Pkg.activate(".")

# Check for the `--coverage` flag
coverage = any(ARGS .== "--coverage")

# A failure is remembered rather than rethrown, so that coverage is still written; the exit
# status at the end reports it, so that `./runtests` and anything else calling this script
# can tell that the tests failed.
tests_passed = try
    Δt = @elapsed Pkg.test("SphericalFunctions"; coverage, test_args=ARGS)
    println("Running tests took $Δt seconds.")
    true
catch e
    if coverage
        println("Tests failed; proceeding to coverage")
    else
        println("Tests failed.")
    end
    false
end

if coverage
    cd((@__DIR__) * "/..")
    Pkg.activate("test")  # The test environment, which lists `Coverage`
    using Coverage
    coverage = Coverage.process_folder("src")
    Coverage.writefile("lcov.info", coverage)
    Coverage.clean_folder(".")
end

tests_passed || exit(1)
