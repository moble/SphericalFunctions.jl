@testset verbose=true "SSHT" begin

    function sYlm(s, ℓ, m, θϕ)
        NINJA.sYlm(s, ℓ, m, θϕ[1], θϕ[2])
    end

    # These test the ability of ssht to precisely reconstruct a pure `sYlm`.
    @testset "Synthesis: $T $method" for (method, T) in Iterators.product(
        ["Direct", "RS"], #["Direct", "EKKM", "RS"],
        [Double64, Float64, Float32]
    )

        # We can't go to very high ℓ, because NINJA.sYlm fails for low-precision numbers
        for ℓmax ∈ 3:7

            # We need ϵ to be huge, seemingly mostly due to the low-precision method
            # used for NINJA.sYlm; it is used because it is a simple reference method.
            ϵ = 500ℓmax^3 * eps(T)

            for s in -2:2
                𝒯 = SSHT(s, ℓmax; T=T, method=method)

                #for ℓmin in 0:abs(s)
                let ℓmin = abs(s)
                    for ℓ in abs(s):ℓmax
                        for m in -ℓ:ℓ
                            f = zeros(Complex{T}, SphericalFunctions.Ysize(ℓmin, ℓmax))
                            f[SphericalFunctions.Yindex(ℓ, m, ℓmin)] = one(T)
                            computed = 𝒯 * f
                            expected = sYlm.(s, ℓ, m, pixels(𝒯))
                            if ≉(computed, expected, atol=ϵ, rtol=ϵ)
                                @show method
                                @show T
                                @show ℓmax
                                @show s
                                @show ℓ
                                @show m
                                @show ϵ
                                comp = copy(computed)
                                @. comp[abs(comp)<ϵ]=0
                                @show comp
                                @show expected
                                println("max_diff = ", maximum(abs, computed .- expected), ";")
                                println()
                                error("")
                            end
                            @test computed ≈ expected atol=ϵ rtol=ϵ
                        end
                    end
                end
            end
        end
    end  # Synthesis

    # These test the ability of ssht to precisely decompose the results of `sYlm`.
    @testset "Analysis: $T $method" for (method, T) in Iterators.product(
        ["Direct", "RS"], #["Direct", "EKKM", "RS"],
        [Double64, Float64, Float32]
    )

        # We can't go to very high ℓ, because NINJA.sYlm fails for low-precision numbers
        for ℓmax ∈ 3:7

            # We need ϵ to be huge, seemingly mostly due to the low-precision method
            # used for NINJA.sYlm; it is used because it is a simple reference method.
            ϵ = 500ℓmax^3 * eps(T)

            for s in -2:2
                𝒯 = SSHT(s, ℓmax; T=T, method=method)

                #for ℓmin in 0:abs(s)
                let ℓmin = abs(s)
                    for ℓ in abs(s):ℓmax
                        for m in -ℓ:ℓ
                            f = sYlm.(s, ℓ, m, pixels(𝒯))
                            computed = 𝒯 \ f
                            expected = zeros(Complex{T}, size(computed))
                            expected[SphericalFunctions.Yindex(ℓ, m, ℓmin)] = one(T)
                            if ≉(computed, expected, atol=ϵ, rtol=ϵ)
                                @show method
                                @show T
                                @show ℓmax
                                @show s
                                @show ℓ
                                @show m
                                @show ϵ
                                comp = copy(computed)
                                @. comp[abs(comp)<ϵ]=0
                                @show comp
                                @show expected
                                println("max_diff = ", maximum(abs, computed .- expected), ";")
                                println()
                            end
                            @test computed ≈ expected atol=ϵ rtol=ϵ
                        end
                    end
                end
            end
        end
    end  # Analysis
end
