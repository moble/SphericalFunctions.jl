### These benchmarks can be run from the top-level directory of this repo with
###
###     julia -e 'using PkgBenchmark; results=benchmarkpkg("SphericalFunctions"); export_markdown("benchmark/results.md", results)'
###
### This runs the benchmarks (possibly tuning them automatically first), and writes the results to a nice markdown file.

using BenchmarkTools
using SphericalFunctions

const SUITE = BenchmarkGroup()

SUITE["complex_powers"] = BenchmarkGroup(["recursions", "complex"])
for T in [big, Float64, Float32, Float16]
    z = exp(T(6)im/5)
    m = 10_000
    zpowers = zeros(typeof(z), m+1)
    SUITE["complex_powers"][T] = @benchmarkable complex_powers!($zpowers, $z)
end


# SUITE["trigonometry"]["hyperbolic"] = BenchmarkGroup()
# for f in (sin, cos, tan)
#     for x in (0.0, pi)
#         SUITE["trigonometry"]["hyperbolic"][string(f), x] = @benchmarkable ($f)($x)
#     end
# end
