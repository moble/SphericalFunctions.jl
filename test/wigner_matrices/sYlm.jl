@testitem "Issue #40" begin
    @test maximum(abs, sYlm_values(0.0, 0.0, 3, -2)) > 0
end

@testitem "Test NINJA expressions" setup=[NINJA,Utilities] begin
    using ProgressMeter
    @testset "$T" for T in [Float64, Float32, BigFloat]
        ## This is just to test my implementation of the equations give in the paper.
        ## Note that this is a test of the testing code itself, not of the main code.
        tol = 2eps(T)
        @showprogress desc="Test NINJA expressions ($T)" for ι in βrange(T)
            for ϕ in αrange(T)
                @test NINJA.sYlm(-2, 2, 2, ι, ϕ) ≈ NINJA.m2Y22(ι, ϕ) atol=tol rtol=tol
                @test NINJA.sYlm(-2, 2, 1, ι, ϕ) ≈ NINJA.m2Y21(ι, ϕ) atol=tol rtol=tol
                @test NINJA.sYlm(-2, 2, 0, ι, ϕ) ≈ NINJA.m2Y20(ι, ϕ) atol=tol rtol=tol
                @test NINJA.sYlm(-2, 2, -1, ι, ϕ) ≈ NINJA.m2Y2m1(ι, ϕ) atol=tol rtol=tol
                @test NINJA.sYlm(-2, 2, -2, ι, ϕ) ≈ NINJA.m2Y2m2(ι, ϕ) atol=tol rtol=tol
            end
        end
    end
end

@testitem "Compare to NINJA expressions" setup=[NINJA,LAL,Utilities] begin
    using ProgressMeter
    using Quaternionic
    @testset "$T" for T in [Float64, Float32, BigFloat]
        ℓₘₐₓ = 8
        sₘₐₓ = 2
        ℓₘᵢₙ = 0
        tol = ℓₘₐₓ^2 * 2eps(T)  # Mostly because the NINJA.sYlm expressions are inaccurate
        sYlm_storage = sYlm_prep(ℓₘₐₓ, sₘₐₓ, T)

        let R = randn(Rotor{T})
            @test_throws ErrorException sYlm_values!(sYlm_storage, R, sₘₐₓ+1)
            @test_throws ErrorException sYlm_values!(sYlm_storage, R, -sₘₐₓ-1)
        end

        @showprogress desc="Compare to NINJA expressions ($T)" for spin in -sₘₐₓ:sₘₐₓ
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
                    for ℓ in abs(spin):ℓₘₐₓ
                        for m in -ℓ:ℓ
                            sYlm1 = Y[i]
                            sYlm2 = NINJA.sYlm(spin, ℓ, m, ι, ϕ)
                            @test sYlm1 ≈ sYlm2 atol=tol rtol=tol
                            if spin==-2 && T===Float64
                                sYlm3 = LAL.XLALSpinWeightedSphericalHarmonic(ι, ϕ, spin, ℓ, m)
                                @test sYlm1 ≈ sYlm3 atol=tol rtol=tol
                            end
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
