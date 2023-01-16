@testset verbose=true "ssht" begin

    # These test the ability of ssht to precisely decompose the results of `sYlm`.
    @testset "$method $T" for (method, T) in [("RS", Float64)]#Iterators.product(
    #     ["Direct", "RS"], #["Direct", "EKKM", "RS"],
    #     [Double64, Float64, Float32]
    # )
        if method == "RS" && T === BigFloat
            continue
        end

        function sYlm(s, ℓ, m, Rθϕ)
            (θ, ϕ) = to_spherical_coordinates(Rθϕ)
            NINJA.sYlm(s, ℓ, m, θ, ϕ)
        end

        for ℓmax ∈ [7, 15]
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
    end
end
