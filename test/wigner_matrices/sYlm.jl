@testitem "Issue #40" begin
    @test maximum(abs, sYlm_values(0.0, 0.0, 3, -2)) > 0
end

@testitem "Internal consistency" setup=[Utilities] begin
    using ProgressMeter
    using Quaternionic
    @testset "$T" for T in [Float64]
        ℓₘₐₓ = 8
        sₘₐₓ = 2
        ℓₘᵢₙ = 0
        tol = ℓₘₐₓ^2 * 2eps(T)
        sYlm_storage = sYlm_prep(ℓₘₐₓ, sₘₐₓ, T)

        let R = randn(Rotor{T})
            @test_throws ErrorException sYlm_values!(sYlm_storage, R, sₘₐₓ+1)
            @test_throws ErrorException sYlm_values!(sYlm_storage, R, -sₘₐₓ-1)
        end

        for spin in [0, -1, -2]
            for ι in βrange(T)
                for ϕ in αrange(T)
                    R = from_spherical_coordinates(ι, ϕ)
                    Y = sYlm_values!(sYlm_storage, R, spin)
                    Y′ = sYlm_values(ι, ϕ, ℓₘₐₓ, spin)  # Also tests issue #40
                    @test Y[(spin^2+1):end] ≈ Y′ atol=tol rtol=tol
                    i = 1
                    for ℓ in 0:abs(spin)-1
                        for m in -ℓ:ℓ
                            @test Y[i] == 0
                            i += 1
                        end
                    end
                end
            end
        end
    end
end

@testitem "Spin property" setup=[Utilities] begin
    using ProgressMeter
    using Quaternionic
    @testset "$T" for T in [Float64, Float32, BigFloat]
        # Test that ₛYₗₘ(R exp(γ*z/2)) = ₛYₗₘ(R) * exp(-im*s*γ)
        # See https://spherical.readthedocs.io/en/main/SWSHs/
        # for a more detailed explanation
        ℓₘₐₓ = 8
        sₘₐₓ = 2
        ℓₘᵢₙ = 0
        tol = 4ℓₘₐₓ * eps(T)
        sYlm_storage = sYlm_prep(ℓₘₐₓ, sₘₐₓ, T)
        @showprogress desc="Spin property ($T)" for spin in -sₘₐₓ:sₘₐₓ
            for ι in βrange(T)
                for ϕ in αrange(T)
                    for γ in γrange(T)
                        R = from_spherical_coordinates(ι, ϕ)
                        Y1 = copy(sYlm_values!(sYlm_storage, R * exp(γ*𝐤/2), spin))
                        Y2 = sYlm_values!(sYlm_storage, R, spin)
                        @test Y1 ≈ Y2 * cis(-spin*γ) atol=tol rtol=tol
                    end
                end
            end
        end
    end
end

@testitem "sYlm vs WignerD" setup=[Utilities] begin
    using ProgressMeter
    using Quaternionic
    @testset "$T" for T in [Float64, Float32, BigFloat]
        # ₛYₗₘ(R) = (-1)ˢ √((2ℓ+1)/(4π)) 𝔇ˡₘ₋ₛ(R)
        #        = (-1)ˢ √((2ℓ+1)/(4π)) 𝔇̄ˡ₋ₛₘ(R̄)
        ℓₘₐₓ = 8
        sₘₐₓ = 2
        ℓₘᵢₙ = 0
        tol = 4ℓₘₐₓ * eps(T)
        D_storage = D_prep(ℓₘₐₓ, T)
        sYlm_storage = sYlm_prep(ℓₘₐₓ, sₘₐₓ, T)
        @showprogress desc="sYlm vs WignerD ($T)" for s in -sₘₐₓ:sₘₐₓ
            for ι in βrange(T)
                for ϕ in αrange(T)
                    R = from_spherical_coordinates(ι, ϕ)
                    𝔇 = D_matrices!(D_storage, R)
                    Y = sYlm_values!(sYlm_storage, R, s)
                    i = 1
                    for ℓ in 0:abs(s)-1
                        for m in -ℓ:ℓ
                            @test Y[i] == 0
                            i += 1
                        end
                    end
                    for ℓ in abs(s):ℓₘₐₓ
                        for m in -ℓ:ℓ
                            sYlm1 = Y[i]
                            sYlm2 = (-1)^s * √((2ℓ+1)/(4T(π))) * 𝔇[WignerDindex(ℓ, m, -s)]
                            if ≉(sYlm1, sYlm2, atol=tol, rtol=tol)
                                println("Unequal at i=$i (s,ℓ,m)=$((s,ℓ,m)): ")
                                @show ι ϕ
                                @show sYlm1 sYlm2
                            end
                            @test sYlm1 ≈ sYlm2 atol=tol rtol=tol
                            i += 1
                        end
                    end
                end
            end
        end
    end
end

@testitem "sYlm conjugation" setup=[Utilities] begin
    using ProgressMeter
    using Quaternionic
    @testset "$T" for T in [Float64, Float32, BigFloat]
        # ₛȲₗₘ = (-1)ˢ⁺ᵐ ₋ₛYₗ₋ₘ
        ℓₘₐₓ = 8
        sₘₐₓ = 2
        ℓₘᵢₙ = 0
        tol = 4ℓₘₐₓ * eps(T)
        sYlm_storage = sYlm_prep(ℓₘₐₓ, sₘₐₓ, T)
        @showprogress desc="sYlm conjugation ($T)" for ι in βrange(T)
            for ϕ in αrange(T)
                for γ in γrange(T)
                    for s in -sₘₐₓ:sₘₐₓ
                        R = from_spherical_coordinates(ι, ϕ)
                        Y1 = copy(sYlm_values!(sYlm_storage, R, s))
                        Y2 = sYlm_values!(sYlm_storage, R, -s)
                        for ℓ in abs(s):ℓₘₐₓ
                            for m in -ℓ:ℓ
                                sYlm1 = conj(Y1[Yindex(ℓ, m)])
                                sYlm2 = (-1)^(s+m) * Y2[Yindex(ℓ, -m)]
                                @test sYlm1 ≈ sYlm2 atol=tol rtol=tol
                            end
                        end
                    end
                end
            end
        end
    end
end
