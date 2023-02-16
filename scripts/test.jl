# Call this from the top-level directory as
#   julia -t auto --project=. scripts/test.jl
# Optionally, specify the name of a top-level test group — e.g., ssht —
# at the end of this command to only run the tests in that group.  Or add
# --help at the end to see the possibilities.

import Dates
println("Running tests starting at ", Dates.format(Dates.now(), "HH:MM:SS"), ".")

using LinearAlgebra
using Base.Threads
LinearAlgebra.BLAS.set_num_threads(nthreads())

using Pkg
try
    Δt = @elapsed Pkg.test("SphericalFunctions"; coverage=true, test_args=ARGS)
    println("Running tests took $Δt seconds.")
catch e
    println("Tests failed; proceeding to coverage")
end

Pkg.activate()  # Activate Julia's base (home) directory
using Coverage
cd((@__DIR__) * "/..")
coverage = Coverage.process_folder()
Coverage.writefile("lcov.info", coverage)
Coverage.clean_folder(".")
