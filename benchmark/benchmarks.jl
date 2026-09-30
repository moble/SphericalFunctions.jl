### These benchmarks can be run from the top-level directory of this repo with
###
###     julia --project=benchmark -e 'using PkgBenchmark; results=benchmarkpkg("SphericalFunctions"); export_markdown("benchmark/results.md", results)'
###
### This runs the benchmarks (possibly tuning them automatically first), and writes the
### results to a nice markdown file.  The `--project=benchmark` matters: PkgBenchmark runs
### the benchmarks in a new process that uses whichever project is active, and the benchmark
### project is the one that lists what they need.  The tuning is saved in `benchmark/tune.json`
### and reused on later runs; pass `retune=true` to `benchmarkpkg` after adding or changing
### benchmarks, since a benchmark missing from that file is never tuned.  The same suite can
### be run on GitHub with the "benchmarks" workflow.
###
### These are regression benchmarks: they are meant to be small enough to run often, and to
### cover each layer of the package once.  The separate script `per_ell_grid.jl` measures
### how the cost of the per-ℓ recursion depends on ℓₘₐₓ and on the number of rotors in a
### batch, over a much larger grid, and is not part of this suite.

using BenchmarkTools
using Random
using Quaternionic: Rotor
using SphericalFunctions

const SUITE = BenchmarkGroup()
const rng = Random.Xoshiro(1234)

### Complex powers: the innermost helper of every phase calculation.
SUITE["complex_powers"] = BenchmarkGroup(["recursions", "complex"])
for T in [big, Float64, Float32, Float16]
    z = exp(T(6)im/5)
    m = 10_000
    zpowers = zeros(typeof(z), m+1)
    SUITE["complex_powers"][T] = @benchmarkable complex_powers!($zpowers, $z)
end

### The Wigner engine, one ℓ at a time, for a single rotor and for a batch.  `m′ₘₐₓ = 2` is
### the spin-weighted case that the harmonics need; `m′ₘₐₓ = ℓₘₐₓ` is the full matrix.  A
### calculator built from a vector is batched even when the vector holds only one rotor, and
### one built from a single `Rotor` is not, so the case of one rotor is timed both ways.
SUITE["wigner"] = BenchmarkGroup(["recursions"])
for ℓₘₐₓ in (8, 64), m′ₘₐₓ in unique((ℓₘₐₓ, 2)), Nᵣ in (1, 64)
    R⃗ = randn(rng, Rotor{Float64}, Nᵣ)
    calc = DCalculator(R⃗, ℓₘₐₓ; m′ₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ)
    SUITE["wigner"]["D sweep", ℓₘₐₓ, m′ₘₐₓ, Nᵣ] = @benchmarkable begin
        for ℓ in 0:$ℓₘₐₓ
            recurrence!($calc, ℓ)
        end
    end
end
for ℓₘₐₓ in (8, 64), m′ₘₐₓ in unique((ℓₘₐₓ, 2))
    R = randn(rng, Rotor{Float64})
    calc = DCalculator(R, ℓₘₐₓ; m′ₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ)
    SUITE["wigner"]["D sweep (single Rotor)", ℓₘₐₓ, m′ₘₐₓ] = @benchmarkable begin
        for ℓ in 0:$ℓₘₐₓ
            recurrence!($calc, ℓ)
        end
    end
end
let R = randn(rng, Rotor{Float64})
    for ℓₘₐₓ in (8, 64)
        SUITE["wigner"]["D convenience", ℓₘₐₓ] = @benchmarkable D($R, $ℓₘₐₓ)
        SUITE["wigner"]["d convenience", ℓₘₐₓ] = @benchmarkable d(1.1, $ℓₘₐₓ)
    end
end

### The half-integer path, which shares the `m′` ladder with the integer one but reaches it
### through a Clebsch–Gordan seed rather than steps 1–3.  This is here as a guard rather than
### for its own sake: the indices are `HalfOddInteger`s, whose arithmetic lands back in `Int`,
### and if some unanticipated operation ever falls back to `Rational` instead, every answer
### stays correct while this gets roughly thirty times slower.  No correctness test can see
### that; this can.
for ℓₘₐₓ in (15//2, 127//2), Nᵣ in (1, 64)
    R⃗ = randn(rng, Rotor{Float64}, Nᵣ)
    calc = DCalculator(R⃗, ℓₘₐₓ)
    SUITE["wigner"]["D sweep (half-integer)", ℓₘₐₓ, Nᵣ] = @benchmarkable begin
        for ℓ in (1//2):($ℓₘₐₓ)
            recurrence!($calc, ℓ)
        end
    end
end
for ℓₘₐₓ in (15//2, 127//2)
    R = randn(rng, Rotor{Float64})
    calc = DCalculator(R, ℓₘₐₓ)
    SUITE["wigner"]["D sweep (half-integer, single Rotor)", ℓₘₐₓ] = @benchmarkable begin
        for ℓ in (1//2):($ℓₘₐₓ)
            recurrence!($calc, ℓ)
        end
    end
end

### Spin-weighted spherical harmonics: one rotor, a reused calculator, and the dense matrix.
SUITE["sYlm"] = BenchmarkGroup(["harmonics"])
let R = randn(rng, Rotor{Float64}), s = -2
    for ℓₘₐₓ in (8, 64)
        SUITE["sYlm"]["sYlm", ℓₘₐₓ] = @benchmarkable sYlm($R, $ℓₘₐₓ, $s)
        calc = sYlmCalculator(R, ℓₘₐₓ, s)
        SUITE["sYlm"]["sYlmCalculator reused with set_R!", ℓₘₐₓ] = @benchmarkable begin
            set_R!($calc, $R)
            for (ℓ, Yˡ) in $calc
            end
        end
    end
end
for ℓₘₐₓ in (8, 32)
    R⃗ = golden_ratio_spiral_rotors(0, ℓₘₐₓ)
    SUITE["sYlm"]["sYlm_matrix", ℓₘₐₓ] = @benchmarkable sYlm_matrix($R⃗, $ℓₘₐₓ, 0)
end

### Transforms: synthesis and analysis, for each method, at a size where all three are usable.
SUITE["ssht"] = BenchmarkGroup(["transforms"])
for method in ("RS", "Minimal", "Matrix"), ℓₘₐₓ in (8, 24)
    s = -2
    # Out of place throughout, so that repeated samples all see the same input.
    𝒯 = method == "RS" ? SSHT(s, ℓₘₐₓ; method) : SSHT(s, ℓₘₐₓ; method, inplace=false)
    f̃ = randn(rng, ComplexF64, Ysize(abs(s), ℓₘₐₓ))
    f = 𝒯 * f̃
    SUITE["ssht"]["synthesis", method, ℓₘₐₓ] = @benchmarkable $𝒯 * $f̃
    SUITE["ssht"]["analysis", method, ℓₘₐₓ] = @benchmarkable $𝒯 \ $f
end
