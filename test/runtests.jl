using Spherical
using Test, Random, FastTransforms


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
    addtests("alf.jl")
    addtests("weights.jl")
end


if help
    println()
    println("Pass no args to run all tests, or select one or more of the following:")
    for helptest in helptests
        println("    ", helptest)
    end
end
