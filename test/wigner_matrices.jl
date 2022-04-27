@testset verbose=true "wigner_matrices" begin
    include("test_utilities.jl")

    @testset "Compare explicit H{$T}" for T in [Float64, Float32]
        # d_{\ell}^{n,m} = \epsilon_n \epsilon_{-m} H_{\ell}^{n,m},
        epsilon(k) = k>0 ? (-1)^k : 1
        ℓₘₐₓ = 2  # This is the max explicitly coded ℓ
        for m′ₘₐₓ in 0:ℓₘₐₓ
            expiβ = exp(im*rand(0:eps(T):π))
            expiβNaNCheck = complex(NaNCheck{T}(expiβ.re), NaNCheck{T}(expiβ.im))
            NCTN = NaNCheck{T}(NaN)
            Hw = fill(NCTN, SphericalFunctions.WignerHsize(m′ₘₐₓ, ℓₘₐₓ))
            SphericalFunctions.H2!(Hw, expiβNaNCheck, ℓₘₐₓ, m′ₘₐₓ, SphericalFunctions.abd(ℓₘₐₓ, T))
            for n in 0:ℓₘₐₓ
                for m′ in -min(n, m′ₘₐₓ):min(n, m′ₘₐₓ)
                    for m in -n:n
                        d = ExplicitWignerMatrices.d(n, m′, m, expiβ)
                        expected = epsilon(m′) * epsilon(-m) * d
                        actual = Hw[SphericalFunctions.WignerHindex(n, m′, m; m′ₘₐₓ=m′ₘₐₓ)]
                        @test expected ≈ actual atol=30eps(T) rtol=30eps(T)
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

    @testset "Compare H/D indexing $T" for T in [Float64, Float32]
        # d_{\ell}^{n,m} = \epsilon_n \epsilon_{-m} H_{\ell}^{n,m},
        epsilon(k) = k>0 ? (-1)^k : 1
        ℓₘₐₓ = 8  # This is the max explicitly coded ℓ
        for m′ₘₐₓ in 0:ℓₘₐₓ
            expiβ = exp(im*rand(0:eps(T):π))
            expiβNaNCheck = complex(NaNCheck{T}(expiβ.re), NaNCheck{T}(expiβ.im))
            NCTN = NaNCheck{T}(NaN)
            Hw = fill(NCTN, SphericalFunctions.WignerHsize(m′ₘₐₓ, ℓₘₐₓ))
            SphericalFunctions.H2!(Hw, expiβNaNCheck, ℓₘₐₓ, m′ₘₐₓ, SphericalFunctions.abd(ℓₘₐₓ, T))
            𝔇 = fill(NCTN, SphericalFunctions.WignerDsize(0, m′ₘₐₓ, ℓₘₐₓ))
            SphericalFunctions.H2!(𝔇, expiβNaNCheck, ℓₘₐₓ, m′ₘₐₓ, SphericalFunctions.abd(ℓₘₐₓ, T), SphericalFunctions.WignerDindex)
            for n in 0:ℓₘₐₓ
                for m′ in -min(n, m′ₘₐₓ):min(n, m′ₘₐₓ)
                    for m in abs(m′):n
                        Hnm′m = Hw[SphericalFunctions.WignerHindex(n, m′, m; m′ₘₐₓ=m′ₘₐₓ)]
                        𝔇nm′m = 𝔇[SphericalFunctions.WignerDindex(n, m′, m; m′ₘₐₓ=m′ₘₐₓ)]
                        @test Hnm′m == 𝔇nm′m
                    end
                end
            end
        end
    end

    # @testset "Group characers $T" for T in [BigFloat, Float64, Float32]
    #     ℓₘₐₓ = 100
    #     m′ₘₐₓ = ℓₘₐₓ
    #     d = Array{T}(undef, WignerDsize(m′ₘₐₓ, ℓₘₐₓ))
    #     𝔇 = Array{Complex{T}}(undef, WignerDsize(m′ₘₐₓ, ℓₘₐₓ))
    #     for β in rand(0:eps(T):π)
    #         d!(d, expiβ, ℓₘₐₓ, m′ₘₐₓ, abd(ℓₘₐₓ, T), WignerDindex)
    #         D!(𝔇, expiβ, ℓₘₐₓ, m′ₘₐₓ, abd(ℓₘₐₓ, T), WignerDindex)
    #         for j in 0:ℓₘₐₓ
    #             sin_ratio = sin((2j+1)*β/2) / sin(β/2)
    #             i1 = WignerDindex(j, -j, -j)
    #             i2 = WignerDindex(j, j, j)
    #             χʲ = sum(d[i1:i2])
    #             @test χʲ ≈ sin_ratio atol=30eps(T) rtol=30eps(T)
    #             χʲ = sum(𝔇[i1:i2])
    #             @test χʲ ≈ sin_ratio atol=30eps(T) rtol=30eps(T)
    #         end
    #     end
    # end

end
