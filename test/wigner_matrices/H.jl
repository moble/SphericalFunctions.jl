@testset verbose=true "H" begin

    @testset "Test H(0) ($T)" for T in [BigFloat, Float64, Float32]
        # We have H_{n}^{m′,m}(0) = (-1)^m′ δ_{m′,m}
        let β = zero(T)
            expiβ = cis(β)
            for ℓₘₐₓ in 0:6  # Expect overflows for higher ℓ with Float32
                for m′ₘₐₓ in 0:ℓₘₐₓ
                    Hw = fill(T(NaN), WignerHsize(ℓₘₐₓ, m′ₘₐₓ))
                    H!(Hw, expiβ, ℓₘₐₓ, m′ₘₐₓ, H_recursion_coefficients(ℓₘₐₓ, T))
                    for n in 0:ℓₘₐₓ
                        for m′ in -min(n, m′ₘₐₓ):min(n, m′ₘₐₓ)
                            for m in abs(m′):n
                                H_rec = Hw[WignerHindex(n, m′, m, m′ₘₐₓ)]
                                if m′ == m
                                    @test H_rec ≈ T(-1)^m′ atol=eps(T)
                                else
                                    @test H_rec ≈ zero(T) atol=eps(T)
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    @testset "Test H(π) ($T)" for T in [BigFloat, Float64, Float32]
        # We have H_{n}^{m′,m}(0) = (-1)^m′ δ_{m′,m}
        let β = T(π)
            expiβ = cis(β)
            for ℓₘₐₓ in 0:6  # Expect overflows for higher ℓ with Float32
                for m′ₘₐₓ in 0:ℓₘₐₓ
                    Hw = fill(T(NaN), WignerHsize(ℓₘₐₓ, m′ₘₐₓ))
                    H!(Hw, expiβ, ℓₘₐₓ, m′ₘₐₓ, H_recursion_coefficients(ℓₘₐₓ, T))
                    for n in 0:ℓₘₐₓ
                        for m′ in -min(n, m′ₘₐₓ):min(n, m′ₘₐₓ)
                            for m in abs(m′):n
                                H_rec = Hw[WignerHindex(n, m′, m, m′ₘₐₓ)]
                                if m′ == -m
                                    @test H_rec ≈ T(-1)^(n+m′) atol=eps(T(π))
                                else
                                    @test H_rec ≈ zero(T) atol=1.5eps(T(π))
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    @testset "Compare H to explicit d ($T)" for T in [BigFloat, Float64, Float32]
        # This compares the H obtained via recurrence with the explicit Wigner d
        # d_{\ell}^{n,m} = \epsilon_n \epsilon_{-m} H_{\ell}^{n,m},
        for β in βrange(T)
            expiβ = cis(β)
            for ℓₘₐₓ in 0:2  # 2 is the max explicitly coded ℓ
                for m′ₘₐₓ in 0:ℓₘₐₓ
                    Hw = fill(T(NaN), WignerHsize(ℓₘₐₓ, m′ₘₐₓ))
                    H!(Hw, expiβ, ℓₘₐₓ, m′ₘₐₓ, H_recursion_coefficients(ℓₘₐₓ, T))
                    for n in 0:ℓₘₐₓ
                        for m′ in -min(n, m′ₘₐₓ):min(n, m′ₘₐₓ)
                            for m in -n:n
                                d_expl = ExplicitWignerMatrices.d_explicit(n, m′, m, expiβ)
                                d_rec = epsilon(m′) * epsilon(-m) * Hw[WignerHindex(n, m′, m, m′ₘₐₓ)]
                                @test d_rec ≈ d_expl atol=30eps(T) rtol=30eps(T)
                            end
                        end
                    end
                end

                # @show ℓₘₐₓ m′ₘₐₓ
                # for n in 0:ℓₘₐₓ
                #     @show n
                #     mat = fill(zero(T), 2n+1, 2n+1)
                #     for m′ in -min(n, m′ₘₐₓ):min(n, m′ₘₐₓ)
                #         for m in abs(m′):n
                #             mat[m′+n+1, m+n+1] =
                #                 Hw[SphericalFunctions.WignerHindex(n, m′, m, m′ₘₐₓ)].val
                #         end
                #     end
                #     println("H:")
                #     display(mat)
                #     println()
                #     mat = fill(zero(T), 2n+1, 2n+1)
                #     for m′ in -min(n, m′ₘₐₓ):min(n, m′ₘₐₓ)
                #         for m in abs(m′):n
                #             mat[m′+n+1, m+n+1] =
                #                 𝔇[SphericalFunctions.WignerDindex(n, m′, m, m′ₘₐₓ)].val
                #         end
                #     end
                #     println("𝔇 :")
                #     display(mat)
                #     println()
                #     println()
                # end
                # println()
            end
        end
    end

    @testset "Compare H to formulaic d ($T)" for T in [BigFloat, Float64, Float32]
        # This compares the H obtained via recurrence with the formulaic Wigner d
        # d_{\ell}^{n,m} = \epsilon_n \epsilon_{-m} H_{\ell}^{n,m},
        tol = ifelse(T === BigFloat, 100, 1) * 30eps(T)
        for β in βrange(T)
            expiβ = cis(β)
            for ℓₘₐₓ in 0:6  # Expect overflows for higher ℓ with Float32
                for m′ₘₐₓ in 0:ℓₘₐₓ
                    Hw = fill(T(NaN), WignerHsize(ℓₘₐₓ, m′ₘₐₓ))
                    H!(Hw, expiβ, ℓₘₐₓ, m′ₘₐₓ, H_recursion_coefficients(ℓₘₐₓ, T))
                    for n in 0:ℓₘₐₓ
                        for m′ in -min(n, m′ₘₐₓ):min(n, m′ₘₐₓ)
                            for m in -n:n
                                d_form = ExplicitWignerMatrices.d_formula(n, m′, m, expiβ)
                                d_rec = epsilon(m′) * epsilon(-m) * Hw[WignerHindex(n, m′, m, m′ₘₐₓ)]
                                @test d_rec ≈ d_form atol=tol rtol=tol
                            end
                        end
                    end
                end
            end
        end
    end

end
