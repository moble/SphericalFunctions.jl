# Call this from the top-level directory as
#   julia -t auto scripts/test.jl
# or to run with coverage as
#   julia -t auto scripts/test.jl --coverage
# See docs/src/development/index.md for more information.

import Dates
println("Running tests starting at ", Dates.format(Dates.now(), "HH:MM:SS"), ".")

using Pkg
cd((@__DIR__) * "/..")
Pkg.activate(".")

using LinearAlgebra
using Base.Threads
LinearAlgebra.BLAS.set_num_threads(nthreads())

# Check for the `--coverage` flag
coverage = any(ARGS .== "--coverage")

try
    Δt = @elapsed Pkg.test("SphericalFunctions"; coverage, test_args=ARGS)
    println("Running tests took $Δt seconds.")
catch e
    if coverage
        println("Tests failed; proceeding to coverage")
    else
        println("Tests failed.")
    end
end

if coverage
    Pkg.activate()  # Activate Julia's base (home) directory
    using Coverage
    cd((@__DIR__) * "/..")
    coverage = Coverage.process_folder("src")
    Coverage.writefile("lcov.info", coverage)
    Coverage.clean_folder(".")
end
