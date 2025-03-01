# NOTE: Float16 irfft returns Float32, which leads to type conflicts, so we just don't
# test map2salm on Float16

@testitem "map2salm" setup=[Utilities] begin
    for T in [BigFloat, Float64, Float32]
        # These test the ability of map2salm to precisely decompose the results of `sYlm`.
        ℓmax = 7
        Nϑ = 2ℓmax + 1
        Nφ = 2ℓmax + 2
        for s in -2:2
            #for ℓmin in 0:abs(s)
            let ℓmin = 0
                for ℓ in abs(s):ℓmax
                    for m in -ℓ:ℓ
                        f = mapslices(
                            ϕθ -> sYlm(s, ℓ, m, ϕθ[2], ϕθ[1]),
                            phi_theta(Nφ, Nϑ, T),
                            dims=[3]
                        )
                        computed = map2salm(f, s, ℓmax)
                        expected = zeros(Complex{T}, size(computed))
                        expected[SphericalFunctions.Yindex(ℓ, m, ℓmin)] = one(T)
                        if ≉(computed, expected, atol=30eps(T), rtol=30eps(T))
                            @show T
                            @show ℓmax
                            @show Nϑ
                            @show Nφ
                            @show s
                            @show ℓmin
                            @show ℓ
                            @show m
                            println("computed = $computed")
                            println("expected = $expected")
                            println("max_diff = ", maximum(abs, computed .- expected), ";")
                            println()
                        end
                        @test computed ≈ expected atol=30eps(T) rtol=30eps(T)

                        plan = plan_map2salm(f, s, ℓmax)
                        computed2 = map2salm(f, plan)
                        @test array_equal(computed, computed2)
                    end
                end
            end
        end
    end
end
