@testset verbose=true "wigner_matrices" begin
    include("test_utilities.jl")
    using SphericalFunctions: WignerHsize, WignerDsize, WignerHindex, WignerDindex, abd
    αrange(T, n=15) = T[
        0; nextfloat(T(0)); rand(T(0):eps(T(π)):T(π), n÷2); prevfloat(T(π)); T(π);
        nextfloat(T(π)); rand(T(π):eps(2T(π)):2T(π), n÷2); prevfloat(T(π)); 2T(π)
    ]
    βrange(T, n=15) = T[
        0; nextfloat(T(0)); rand(T(0):eps(T(π)):T(π), n); prevfloat(T(π)); T(π)
    ]
    γrange(T, n=15) = αrange(T, n)
    epsilon(k) = ifelse(k>0 && isodd(k), -1, 1)

    @testset "Compare d expressions ($T)" for T in [BigFloat, Float64, Float32]
        # This just compares the two versions of the `d` function from test_utilities
        # to ensure that later tests that use those functions are reliable
        for β in βrange(T)
            expiβ = exp(im*β)
            for ℓₘₐₓ in 0:2  # 2 is the max explicitly coded ℓ
                for m′ₘₐₓ in 0:ℓₘₐₓ
                    for n in 0:ℓₘₐₓ
                        for m′ in -min(n, m′ₘₐₓ):min(n, m′ₘₐₓ)
                            for m in -n:n
                                d_expl = ExplicitWignerMatrices.d_explicit(n, m′, m, expiβ)
                                d_form = ExplicitWignerMatrices.d_formula(n, m′, m, expiβ)
                                @test d_expl ≈ d_form atol=30eps(T) rtol=30eps(T)
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
            expiβ = exp(im*β)
            for ℓₘₐₓ in 0:2  # 2 is the max explicitly coded ℓ
                for m′ₘₐₓ in 0:ℓₘₐₓ
                    Hw = fill(T(NaN), WignerHsize(m′ₘₐₓ, ℓₘₐₓ))
                    H!(Hw, expiβ, ℓₘₐₓ, m′ₘₐₓ, abd(ℓₘₐₓ, T))
                    for n in 0:ℓₘₐₓ
                        for m′ in -min(n, m′ₘₐₓ):min(n, m′ₘₐₓ)
                            for m in -n:n
                                d_expl = ExplicitWignerMatrices.d_explicit(n, m′, m, expiβ)
                                d_rec = epsilon(m′) * epsilon(-m) * Hw[WignerHindex(n, m′, m; m′ₘₐₓ=m′ₘₐₓ)]
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
                #                 Hw[SphericalFunctions.WignerHindex(n, m′, m; m′ₘₐₓ=m′ₘₐₓ)].val
                #         end
                #     end
                #     println("H:")
                #     display(mat)
                #     println()
                #     mat = fill(zero(T), 2n+1, 2n+1)
                #     for m′ in -min(n, m′ₘₐₓ):min(n, m′ₘₐₓ)
                #         for m in abs(m′):n
                #             mat[m′+n+1, m+n+1] =
                #                 𝔇[SphericalFunctions.WignerDindex(n, m′, m; m′ₘₐₓ=m′ₘₐₓ)].val
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
        for β in βrange(T)
            expiβ = exp(im*β)
            for ℓₘₐₓ in 0:6  # Expect overflows for higher ℓ with Float32
                for m′ₘₐₓ in 0:ℓₘₐₓ
                    Hw = fill(T(NaN), WignerHsize(m′ₘₐₓ, ℓₘₐₓ))
                    H!(Hw, expiβ, ℓₘₐₓ, m′ₘₐₓ, abd(ℓₘₐₓ, T))
                    for n in 0:ℓₘₐₓ
                        for m′ in -min(n, m′ₘₐₓ):min(n, m′ₘₐₓ)
                            for m in -n:n
                                d_form = ExplicitWignerMatrices.d_formula(n, m′, m, expiβ)
                                d_rec = epsilon(m′) * epsilon(-m) * Hw[WignerHindex(n, m′, m; m′ₘₐₓ=m′ₘₐₓ)]
                                @test d_rec ≈ d_form atol=30eps(T) rtol=30eps(T)
                            end
                        end
                    end
                end
            end
        end
    end

    @testset "Test H(0) ($T)" for T in [BigFloat, Float64, Float32]
        # We have H_{n}^{m′,m}(0) = (-1)^m′ δ_{m′,m}
        let β = zero(T)
            expiβ = exp(im*β)
            for ℓₘₐₓ in 0:6  # Expect overflows for higher ℓ with Float32
                for m′ₘₐₓ in 0:ℓₘₐₓ
                    Hw = fill(T(NaN), WignerHsize(m′ₘₐₓ, ℓₘₐₓ))
                    H!(Hw, expiβ, ℓₘₐₓ, m′ₘₐₓ, abd(ℓₘₐₓ, T))
                    for n in 0:ℓₘₐₓ
                        for m′ in -min(n, m′ₘₐₓ):min(n, m′ₘₐₓ)
                            for m in abs(m′):n
                                H_rec = Hw[WignerHindex(n, m′, m; m′ₘₐₓ=m′ₘₐₓ)]
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
            expiβ = exp(im*β)
            for ℓₘₐₓ in 0:6  # Expect overflows for higher ℓ with Float32
                for m′ₘₐₓ in 0:ℓₘₐₓ
                    Hw = fill(T(NaN), WignerHsize(m′ₘₐₓ, ℓₘₐₓ))
                    H!(Hw, expiβ, ℓₘₐₓ, m′ₘₐₓ, abd(ℓₘₐₓ, T))
                    for n in 0:ℓₘₐₓ
                        for m′ in -min(n, m′ₘₐₓ):min(n, m′ₘₐₓ)
                            for m in abs(m′):n
                                H_rec = Hw[WignerHindex(n, m′, m; m′ₘₐₓ=m′ₘₐₓ)]
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

    @testset "Compare H/D indexing ($T)" for T in [Float64, Float32]
        # Here, we check that we can pass in either an "H wedge" array to be used with
        # WignerHindex, or a full 𝔇 array used with WignerDindex, and obtain the same
        # H recurrence results
        ℓₘₐₓ = 8
        for m′ₘₐₓ in 0:ℓₘₐₓ
            expiβ = exp(im*rand(0:eps(T):π))
            expiβNaNCheck = complex(NaNCheck{T}(expiβ.re), NaNCheck{T}(expiβ.im))
            NCTN = NaNCheck{T}(NaN)
            Hw = fill(NCTN, WignerHsize(m′ₘₐₓ, ℓₘₐₓ))
            H!(Hw, expiβNaNCheck, ℓₘₐₓ, m′ₘₐₓ, abd(ℓₘₐₓ, T))
            𝔇 = fill(NCTN, WignerDsize(0, m′ₘₐₓ, ℓₘₐₓ))
            H!(𝔇, expiβNaNCheck, ℓₘₐₓ, m′ₘₐₓ, abd(ℓₘₐₓ, T), WignerDindex)
            for n in 0:ℓₘₐₓ
                for m′ in -min(n, m′ₘₐₓ):min(n, m′ₘₐₓ)
                    for m in abs(m′):n
                        Hnm′m = Hw[WignerHindex(n, m′, m; m′ₘₐₓ=m′ₘₐₓ)]
                        𝔇nm′m = 𝔇[WignerDindex(n, m′, m; m′ₘₐₓ=m′ₘₐₓ)]
                        @test Hnm′m == 𝔇nm′m
                    end
                end
            end
        end
    end

    @testset "Compare d to formulaic d ($T)" for T in [Float64]
        epsilon(k) = k>0 ? (-1)^k : 1
        for β in βrange(T, 2)#[2:end-1]
            expiβ = exp(im*β)
            @info "Fix ℓₘₐₓ and remove show_me"
            for ℓₘₐₓ in 0:4
                abd_vals = abd(ℓₘₐₓ, T)
                d = Array{T}(undef, WignerDsize(0, ℓₘₐₓ, ℓₘₐₓ))
                d!(d, expiβ, ℓₘₐₓ, abd_vals)
                for n in 0:ℓₘₐₓ

                    show_me = false

                    for m′ in -n:n
                        for m in -n:n
                            d_formula = ExplicitWignerMatrices.d_formula(n, m′, m, expiβ)
                            d_recurrence = d[WignerDindex(n, m′, m)]
                            if ≉(d_formula, d_recurrence, atol=30eps(T), rtol=30eps(T))
                                show_me = true
                                @show (n, m′, m) ℓₘₐₓ WignerDindex(n, m′, m) d_formula d_recurrence
                            end
                            @test d_formula ≈ d_recurrence atol=30eps(T) rtol=30eps(T)
                        end
                    end

                    if show_me
                        mat = fill(zero(T), 2n+1, 2n+1)
                        for m′ in -n:n
                            for m in -n:n
                                mat[m′+n+1, m+n+1] = d[WignerDindex(n, m′, m)]
                            end
                        end
                        println("d_recurrence{$n}:")
                        display(mat)
                        println()
                        mat = fill(zero(BigFloat), 2n+1, 2n+1)
                        for m′ in -n:n
                            for m in -n:n
                                mat[m′+n+1, m+n+1] = ExplicitWignerMatrices.d_formula(n, m′, m, big(expiβ))
                            end
                        end
                        println("d_formula{$n}:")
                        display(mat)
                        println()
                    end
                end
            end
        end
    end

    # @testset "Group characers $T" for T in [Float64, Float32]
    #     ℓₘₐₓ = 100
    #     m′ₘₐₓ = ℓₘₐₓ
    #     abd_vals = abd(ℓₘₐₓ, T)
    #     d = Array{T}(undef, WignerDsize(0, m′ₘₐₓ, ℓₘₐₓ))
    #     𝔇 = Array{Complex{T}}(undef, WignerDsize(0, m′ₘₐₓ, ℓₘₐₓ))
    #     for β in rand(T(0):eps(T(π)):T(π))
    #         expiβ = exp(im*β)
    #         @show typeof(expiβ)
    #         d!(d, expiβ, ℓₘₐₓ, abd_vals)
    #         #D!(𝔇, expiβ, ℓₘₐₓ, abd_vals)
    #         for j in 0:ℓₘₐₓ
    #             sin_ratio = sin((2j+1)*β/2) / sin(β/2)
    #             i1 = WignerDindex(j, -j, -j)
    #             i2 = WignerDindex(j, j, j)
    #             χʲ = sum(d[WignerDindex(j, m, m)] for m in -j:j)
    #             @test χʲ ≈ sin_ratio atol=30eps(T) rtol=30eps(T)
    #             #χʲ = sum(𝔇[i1:i2])
    #             #@test χʲ ≈ sin_ratio atol=30eps(T) rtol=30eps(T)
    #         end
    #     end
    # end

end
