using SphericalFunctions
using Test
using Random, FFTW, FastTransforms
using OffsetArrays, StaticArrays, Quaternionic
using DoubleFloats, ProgressMeter
using LinearAlgebra, Base.Threads
import Hwloc: num_physical_cores

# Match the number of BLAS threads to the number of threads Julia was started with.
# This may not always be a great idea, but will be useful for these tests.
LinearAlgebra.BLAS.set_num_threads(nthreads())

include("test_utilities.jl")

enabled_tests = lowercase.(ARGS)

help = ("help" ∈ enabled_tests || "--help" ∈ enabled_tests)
helptests = []

# This block is cribbed from StaticArrays.jl/test/runtests.jl
function addtests(fname)
    key = lowercase(splitext(fname)[1])
    if help
        push!(helptests, key)
    else
        if isempty(enabled_tests) || key in enabled_tests
            println("Running $key.jl")
            Random.seed!(42)
            include(fname)
        end
    end
end


@testset verbose=true "All tests" begin
    addtests("complex_powers.jl")
    addtests("indexing.jl")
    addtests("iterators.jl")
    addtests("associated_legendre.jl")
    addtests("wigner_matrices.jl")
    addtests("weights.jl")
    addtests("map2salm.jl")
    addtests("ssht.jl")
end


if help
    println()
    println("Pass no args to run all tests, or select one or more of the following:")
    for helptest in helptests
        println("    ", helptest)
    end
end
