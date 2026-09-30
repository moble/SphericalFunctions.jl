### How the cost of the per-ℓ H recursion depends on ℓₘₐₓ and on the number of rotors.
###
### The engine computes one ℓ at a time for a batch of Nᵣ rotors at once, which is what the
### transforms need and what keeps memory bounded at large ℓₘₐₓ; the alternative, one large
### array of H values for all ℓ at once, for a single rotor at a time, needs about 2.7 GB at
### ℓₘₐₓ = 1000.  The design question is how much the per-ℓ form costs for a single rotor and
### how much it gains for a batch.  The package has no second implementation to time against,
### so what this script measures is the *shape* of the cost, which answers the same question:
###
###   * Absolute nanoseconds per H element per rotor, over the grid.  For comparison, the
###     whole-array code of the 2.x releases took 0.5 to 1.1 ns per element per rotor at
###     ℓₘₐₓ ≥ 64, and 3.6 ns at ℓₘₐₓ = 8.
###   * The ratio of the single-rotor cost to the large-batch cost at the same ℓₘₐₓ.  Whatever
###     is left over at Nᵣ = 1 is per-ℓ fixed cost: mostly the setup of the loop over rotors,
###     and then the square roots of the recursion coefficients.  The engine tabulates the
###     coefficients of steps 4 and 5 once per ℓ and writes their single-rotor case as one
###     statement, which together remove most of that cost for the full wedge.  For the axis
###     alone, which is all that spin weight 0 needs, the fixed cost is instead that of the
###     coefficients of step 2, which only a table kept for the lifetime of the calculator, of
###     O(ℓₘₐₓ²) entries, would remove.  A ratio near 1 means there is nothing to win.
###
### Run it on an otherwise idle machine, with one thread, from the package root:
###
###     julia --project=benchmark -t 1 benchmark/per_ell_grid.jl
###
### Optional arguments narrow the grid, e.g. `... per_ell_grid.jl 8,64 1,8`.  The default
### grid is large: its cell with ℓₘₐₓ = m′ₘₐₓ = 1024 and Nᵣ = 512 needs about 4.5 GB of memory
### for the calculator alone, and takes several minutes, because each cell runs seven full
### sweeps (one to warm up, the best of five, and one to count allocations).

using SphericalFunctions
using SphericalFunctions: HCalculator, recurrence!
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
    w = HCalculator(angles, ℓₘₐₓ; m′ₘₐₓ)  # `angles::Vector{T}` fixes the element type
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
        The per-rotor cost at Nᵣ = $(first(Nᵣs)) divided by the cost at Nᵣ = $(last(Nᵣs)) is the
        per-ℓ fixed overhead of the single-rotor path.  For the axis alone that overhead is the
        cost of the coefficients of step 2, which a table of O(ℓₘₐₓ²) coefficients kept for the
        lifetime of the calculator would remove; such a table is worth adding only if that
        ratio exceeds about 2.""")
end

main()
