@testset verbose=true "ssht" begin

    @testset "$method $T" for (method, T) in Iterators.product(
        ["Direct", "RS"], #["Direct", "EKKM", "RS"],
        [BigFloat, Float64, Float32]
    )
        if method == "RS" && T === BigFloat
            continue
        end

        ϵ = 50eps(T)

        # These test the ability of ssht to precisely decompose the results of `sYlm`.
        ℓmax = 7

        function sYlm(s, ℓ, m, Rθϕ)
            θϕ = to_spherical_coordinates(Rθϕ)
            NINJA.sYlm(s, ℓ, m, θϕ[1], θϕ[2])
        end

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
                            @show T
                            @show ℓmax
                            @show s
                            @show ℓ
                            @show m
                            println("computed = $computed")
                            println("expected = $expected")
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
