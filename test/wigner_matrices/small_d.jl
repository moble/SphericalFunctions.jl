@testitem "Compare d expressions" setup=[ExplicitWignerMatrices,Utilities] begin
    @testset "$T" for T in [BigFloat, Float64, Float32]
        # This just compares the two versions of the `d` function from test_utilities
        # to ensure that later tests that use those functions are reliable
        tol = ifelse(T === BigFloat, 100, 1) * 30eps(T)
        for β in βrange(T)
            expiβ = cis(β)
            for ℓₘₐₓ in 0:2  # 2 is the max explicitly coded ℓ
                for m′ₘₐₓ in 0:ℓₘₐₓ
                    for n in 0:ℓₘₐₓ
                        for m′ in -min(n, m′ₘₐₓ):min(n, m′ₘₐₓ)
                            for m in -n:n
                                d_expl = ExplicitWignerMatrices.d_explicit(n, m′, m, expiβ)
                                d_form = ExplicitWignerMatrices.d_formula(n, m′, m, expiβ)
                                @test d_expl ≈ d_form atol=tol rtol=tol
                            end
                        end
                    end
                end
            end
        end
    end
end

@testitem "Compare d to formulaic d" setup=[ExplicitWignerMatrices,Utilities] begin
    @testset "$T" for T in [BigFloat, Float64, Float32]
        # Now, we're ready to check that d_{n}^{m′,m}(β) matches the expected values
        # for a range of β values
        tol = ifelse(T === BigFloat, 100, 1) * 30eps(T)
        for β in βrange(T)
            expiβ = cis(β)
            for ℓₘₐₓ in 0:4
                d_storage = d_prep(ℓₘₐₓ, T)
                d = d_matrices!(d_storage, expiβ)
                for n in 0:ℓₘₐₓ
                    for m′ in -n:n
                        for m in -n:n
                            d_formula = ExplicitWignerMatrices.d_formula(n, m′, m, expiβ)
                            d_recurrence = d[WignerDindex(n, m′, m)]
                            @test d_formula ≈ d_recurrence atol=tol rtol=tol
                        end
                    end
                end
            end
        end
    end
end

@testitem "Test d signatures" setup=[Utilities] begin
    @testset "$T" for T in [BigFloat, Float64, Float32]
        # 1 d_matrices(β, ℓₘₐₓ)
        # 2 d_matrices(expiβ, ℓₘₐₓ)
        # 3 d_matrices!(d_storage, β)
        # 4 d_matrices!(d_storage, expiβ)
        # 5 d_matrices!(d, β, ℓₘₐₓ)
        # 6 d_matrices!(d, expiβ, ℓₘₐₓ)
        ℓₘₐₓ = 8
        for β in βrange(T)
            expiβ = cis(β)
            dA = d_matrices(β, ℓₘₐₓ)  # 1
            dB = d_matrices(expiβ, ℓₘₐₓ)  # 2
            @test array_equal(dA, dB)
            dB .= 0
            d_matrices!(dB, β, ℓₘₐₓ)  # 5
            @test array_equal(dA, dB)
            dB .= 0
            d_matrices!(dB, expiβ, ℓₘₐₓ)  # 6
            @test array_equal(dA, dB)
            d_storage = d_prep(ℓₘₐₓ, T)
            dB = d_matrices!(d_storage, β)  # 3
            @test array_equal(dA, dB)
            dB .= 0
            dB = d_matrices!(d_storage, expiβ)  # 4
            @test array_equal(dA, dB)
        end
    end
end
