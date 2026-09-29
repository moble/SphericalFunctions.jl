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
# by tag — written as `:sometag` — because that is all the CI workflows use.  Any other
# argument is ignored, with a warning, so that a request for part of the suite does not
# silently run all of it; filtering by name or file is left to `juliati` and the MCP runner.
const CI = get(ENV, "CI", "false") == "true"
const requested_tags = Symbol[Symbol(a[2:end]) for a ∈ ARGS if startswith(a, ":")]
for a ∈ ARGS
    startswith(a, ":") || @warn(
        "`runtests.jl` filters only by tag, written as `:sometag`, so this argument is "
        * "ignored.  To filter by name or file, use `juliati` or the `julia` MCP server.",
        argument=a
    )
end

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

# Facts about the machine and the run, printed before any item runs, so that they appear at
# the top of a CI log.  Whether `muladd` fuses is measured on compiled code, with the inputs
# hidden from constant folding; the case is the modulus of `cis(0.3)`, which is 1 - eps/2
# when rounded once and 1 when rounded twice.  The dynamic calls in the half-integer
# recurrence are counted in a fresh process with the same options, because asking for them
# here would fill this process's inference cache and could change what the items see.
let
    fused(x, y) = muladd(x, x, y * y)
    z = Base.inferencebarrier(cis(0.3))::ComplexF64
    have_fma = isdefined(Core.Intrinsics, :have_fma) ? Core.Intrinsics.have_fma(Float64) : missing
    @info("Test-run diagnostics", VERSION, Sys.CPU_NAME, Sys.ARCH, Threads.nthreads(),
        JULIA_CPU_TARGET=get(ENV, "JULIA_CPU_TARGET", "(unset)"),
        code_coverage=Base.JLOptions().code_coverage, check_bounds=Base.JLOptions().check_bounds,
        have_fma, muladd_fuses=(fused(z.re, z.im) == fma(z.re, z.re, z.im * z.im)),
        muladd=fused(z.re, z.im), fma=fma(z.re, z.re, z.im * z.im))
    script = """
        using SphericalFunctions: SphericalFunctions, HCalculator
        builtin(f::GlobalRef) = isdefined(f.mod, f.name) && getglobal(f.mod, f.name) isa Core.Builtin
        builtin(f) = f isa Core.Builtin
        for H ∈ (HCalculator(0.3, 7//2), HCalculator([0.3, 0.4], 7//2)),
                step! ∈ (SphericalFunctions.recurrence_step4!, SphericalFunctions.recurrence_step5!,
                         SphericalFunctions.recurrence_seed!)
            code = only(Base.code_typed(step!, (typeof(H),); optimize=true)).first.code
            calls = filter(ex -> Meta.isexpr(ex, :call) && !builtin(ex.args[1]), code)
            println("  ", nameof(step!), " for ", typeof(H), ": ", length(calls), " dynamic calls")
            foreach(c -> println("      ", c), calls)
        end
        """
    println("Dynamic calls in the half-integer recurrence, in a fresh process:")
    try
        run(`$(Base.julia_cmd()) --project=$(Base.active_project()) -e $script`)
    catch e
        @warn "The fresh process failed" exception=e
    end
    flush(stdout); flush(stderr)
end

@run_package_tests verbose = true filter = testfilter

# Including the test files is not needed for discovery — `@run_package_tests` finds them on
# its own — but it makes `Pkg.test` parse each one, so a syntax error shows up here rather
# than as a silently missing test item.  Every file under `test/` is included, so that a new
# one cannot be forgotten, except those in hidden directories: `test/.CondaPkg`, the Conda
# environment that the `:python` items build, holds thousands of files that are not part of
# the suite, and a `.jl` file among them must not be run here.  The comparisons with
# `@__DIR__` and `@__FILE__` are parenthesized because a macro called without parentheses
# takes the rest of the expression, `&&` and all, as its arguments.
for (root, _, files) ∈ walkdir(@__DIR__)
    (root != @__DIR__) && any(startswith("."), splitpath(relpath(root, @__DIR__))) && continue
    for file ∈ sort(files)
        path = joinpath(root, file)
        endswith(file, ".jl") && (path != @__FILE__) && include(path)
    end
end
