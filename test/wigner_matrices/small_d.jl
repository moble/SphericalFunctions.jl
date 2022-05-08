@testset verbose=true "d" begin

    @testset "Compare d expressions ($T)" for T in [BigFloat, Float64, Float32]
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

    @testset "Compare d to formulaic d ($T)" for T in [BigFloat, Float64, Float32]
        # Now, we're ready to check that d_{n}^{m′,m}(β) matches the expected values
        # for a range of β values
        tol = ifelse(T === BigFloat, 100, 1) * 30eps(T)
        for β in βrange(T)
            expiβ = cis(β)
            for ℓₘₐₓ in 0:4
                d, H_rec_coeffs = dprep(ℓₘₐₓ, T)
                d!(d, expiβ, ℓₘₐₓ, H_rec_coeffs)
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

    @testset "Test d signatures ($T)" for T in [BigFloat, Float64, Float32]
        # 1 d!(𝔡, expiβ::Complex{T}, ℓₘₐₓ) where {T<:Real}
        # 2 d!(𝔡, β::T, ℓₘₐₓ, (a,b,d)) where {T<:Real}
        # 3 d!(𝔡, β::T, ℓₘₐₓ) where {T<:Real}
        # 4 d(expiβ::Complex{T}, ℓₘₐₓ) where {T<:Real}
        # 5 d(β::Real, ℓₘₐₓ)
        ℓₘₐₓ = 8
        for β in βrange(T)
            expiβ = cis(β)
            dA = d(β, ℓₘₐₓ)
            dB = d(expiβ, ℓₘₐₓ)
            @test array_equal(dA, dB)
            dB .= 0
            d!(dB, β, ℓₘₐₓ)
            @test array_equal(dA, dB)
            dB .= 0
            d!(dB, expiβ, ℓₘₐₓ)
            @test array_equal(dA, dB)
            dB .= 0
            H_rec_coeffs = H_recursion_coefficients(ℓₘₐₓ, T)
            d!(dB, β, ℓₘₐₓ, H_rec_coeffs)
            @test array_equal(dA, dB)
        end
    end

end
