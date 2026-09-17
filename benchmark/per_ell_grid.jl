### The benchmark grid of the v3 design memo, section 9.
###
### The version-2 code filled one large array of H values for all ℓ at once, for a single
### rotor at a time.  The version-3 engine computes one ℓ at a time for a batch of Nᵣ rotors
### at once, which is what the transforms need and what keeps memory bounded at large ℓₘₐₓ
### (the whole-array form needed about 2.7 GB at ℓₘₐₓ = 1000).  The memo asked how much that
### costs for a single rotor and how much it gains for a batch.  The version-2 implementation
### was deleted in 3.0, so there is no longer a second implementation to time against; what
### this script measures instead is the *shape* of the cost, which answers the same design
### question:
###
###   * Absolute nanoseconds per H element per rotor, over the grid.  This is directly
###     comparable with the version-2 numbers recorded in memo section 9 (0.5 to 1.1 ns per
###     element per rotor at ℓₘₐₓ ≥ 64, 3.6 ns at ℓₘₐₓ = 8).
###   * The ratio of the single-rotor cost to the large-batch cost at the same ℓₘₐₓ.  Whatever
###     is left over at Nᵣ = 1 is per-ℓ fixed cost — recomputing recursion coefficients,
###     setting up the loop — and that is exactly what precomputing the coefficients would
###     remove.  A ratio near 1 means there is nothing to win.
###
### Run it on an otherwise idle machine, with one thread, from the package root:
###
###     julia --project=benchmark -t 1 benchmark/per_ell_grid.jl
###
### Optional arguments narrow the grid, e.g. `... per_ell_grid.jl 8,64 1,8`.

using SphericalFunctions
using SphericalFunctions: WignerHCalculator, recurrence!
using Printf

const T = Float64
const β = 1.1  # a generic angle: no pole, no symmetry

parse_list(s, default) = isnothing(s) ? default : parse.(Int, split(s, ","))
const ℓs = parse_list(get(ARGS, 1, nothing), [8, 64, 256, 1024])
const Nᵣs = parse_list(get(ARGS, 2, nothing), [1, 8, 64, 512])

"Best of `n` runs of `f`, in seconds.  The best is used rather than the mean because we are
after the cost of the computation, not of whatever else the machine was doing."
function best_time(f; n=5)
    t = Inf
    for _ ∈ 1:n
        t = min(t, @elapsed f())
    end
    t
end

"Number of H elements the wedge stores for this ℓₘₐₓ and m′ₘₐₓ, summed over ℓ."
function wedge_elements(ℓₘₐₓ, m′ₘₐₓ)
    total = 0
    for ℓ ∈ 0:ℓₘₐₓ
        mp = min(m′ₘₐₓ, ℓ)
        total += sum(ℓ - abs(m′) + 1 for m′ ∈ -mp:mp)
    end
    total
end

function bench(ℓₘₐₓ, m′ₘₐₓ, Nᵣ)
    angles = fill(T(β), Nᵣ)
    w = WignerHCalculator(angles, ℓₘₐₓ; m′ₘₐₓ)  # `angles::Vector{T}` fixes the element type
    function run()
        for ℓ ∈ 0:ℓₘₐₓ
            recurrence!(w, ℓ)
        end
    end
    run()  # warm up
    (best_time(run), @allocated(run()))
end

function main()
    println("threads = ", Threads.nthreads(), ", T = ", T, ", β = ", β)
    println("ns per H element per rotor, for a full sweep ℓ = 0 … ℓₘₐₓ\n")
    @printf("%6s %6s | %s\n", "ℓₘₐₓ", "m′ₘₐₓ",
            join((@sprintf("%12s", "Nᵣ=$Nᵣ") for Nᵣ ∈ Nᵣs), " "))
    worst_overhead = 0.0
    for ℓₘₐₓ ∈ ℓs, m′ₘₐₓ ∈ unique((ℓₘₐₓ, 2))
        elements = wedge_elements(ℓₘₐₓ, m′ₘₐₓ)
        cells = String[]
        percost = Float64[]
        for Nᵣ ∈ Nᵣs
            (t, alloc) = bench(ℓₘₐₓ, m′ₘₐₓ, Nᵣ)
            ns = 1e9t / (elements * Nᵣ)
            push!(percost, ns)
            push!(cells, @sprintf("%12.3f", ns))
            alloc == 0 || @warn "allocated $alloc bytes" ℓₘₐₓ m′ₘₐₓ Nᵣ
        end
        @printf("%6d %6d | %s\n", ℓₘₐₓ, m′ₘₐₓ, join(cells, " "))
        length(percost) > 1 && (worst_overhead = max(worst_overhead, percost[1] / percost[end]))
    end
    println()
    @printf("worst single-rotor overhead (Nᵣ=%d cost / Nᵣ=%d cost): %.2f×\n",
            first(Nᵣs), last(Nᵣs), worst_overhead)
    println("""
        Memo section 9 asked whether the single-rotor path needs precomputed recursion
        coefficients.  The per-rotor cost at Nᵣ = $(first(Nᵣs)) divided by the cost at
        Nᵣ = $(last(Nᵣs)) is the per-ℓ fixed overhead that precomputation would remove; the
        memo's threshold for acting on it was a factor of 2.""")
end

main()
