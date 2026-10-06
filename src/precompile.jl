# The precompile workload.  PrecompileTools runs the calls below while the package is being
# precompiled, so that the code they reach is compiled then and saved with the package,
# rather than compiled on the first call of each in a session.  The calls are those that
# most sessions begin with, in `Float64` and at small ℓ, since the code compiled does not
# depend on the size of ℓ:
#
#   - `D` and `d`, and the iteration of a `DCalculator` and an `sYlmCalculator`, each for
#     integer and for half-integer indices;
#   - the iteration of a `DCalculator` of a vector of rotors, for both kinds of index;
#   - `sYlm` for one rotor and for a vector of rotors, for both kinds of index, and `Ylm`;
#   - a round of the operations on `ModeWeights`: `setindex!`, evaluation `w(R)`, the
#     products `sYlm(R, ℓₘₐₓ, s) * w` and `D(R, ℓₘₐₓ) * w`, the operators `ð`, `ð̄` and
#     `L²`, `+`, and multiplication and division by a scalar;
#   - `SSHT(s, ℓₘₐₓ; method)` with `*` and `\` for the "RS", "Minimal" and "Matrix" methods,
#     and for the "RS" and "Matrix" methods at a half-integer spin weight.
#
# Setting the preference `precompile_workload = false` for this package, in a
# `LocalPreferences.toml`, skips the workload, so that the package precompiles quickly after
# each edit during development, at the cost of slower first calls.
PrecompileTools.@setup_workload begin
    R = from_spherical_coordinates(0.3, 0.7)
    R⃗ = [R, from_spherical_coordinates(1.1, 2.0)]
    PrecompileTools.@compile_workload begin
        D(R, 4)
        D(R, 7//2)
        d(0.3, 4)
        d(0.3, 7//2)
        for ℓₘₐₓ ∈ (4, 7//2)
            for (ℓ, 𝔇ˡ) ∈ DCalculator(R, ℓₘₐₓ)
            end
        end
        for ℓₘₐₓ ∈ (4, 7//2)
            for (ℓ, 𝔇ˡ) ∈ DCalculator(R⃗, ℓₘₐₓ)
            end
        end
        for (ℓₘₐₓ, s) ∈ ((4, -2), (7//2, 1//2))
            for (ℓ, Yˡ) ∈ sYlmCalculator(R, ℓₘₐₓ, s)
            end
        end
        sYlm(R, 4, -2)
        sYlm(R⃗, 4, -2)
        sYlm(R, 7//2, 1//2)
        sYlm(R⃗, 7//2, 1//2)
        Ylm(R, 4)

        w = ModeWeights(zeros(ComplexF64, Ysize(2, 4)), -2)
        w[2, 1] = 1
        w(R)
        sYlm(R, 4, -2) * w
        D(R, 4) * w
        ð * w
        ð̄ * w
        L² * w
        w + w
        2w
        w / 2

        wₕ = ModeWeights(zeros(ComplexF64, Ysize(1//2, 7//2)), 1//2)
        for method ∈ ("RS", "Minimal", "Matrix")
            𝒯 = SSHT(-2, 4; method)
            𝒯 \ (𝒯 * w)
        end
        for method ∈ ("RS", "Matrix")
            𝒯 = SSHT(1//2, 7//2; method)
            𝒯 \ (𝒯 * wₕ)
        end
    end
end
