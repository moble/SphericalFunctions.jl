@testset verbose=true "sYlm" begin

    @testset "Test NINJA expressions ($T)" for T in [Float64, Float32, BigFloat]
        ## This is just to test my implementation of the equations give in the paper.
        ## Note that this is a test of the testing code itself, not of the main code.
        tol = 2eps(T)
        @showprogress "Test NINJA expressions ($T)" for ι in βrange(T)
            for ϕ in αrange(T)
                @test NINJA.sYlm(-2, 2, 2, ι, ϕ) ≈ NINJA.m2Y22(ι, ϕ) atol=tol rtol=tol
                @test NINJA.sYlm(-2, 2, 1, ι, ϕ) ≈ NINJA.m2Y21(ι, ϕ) atol=tol rtol=tol
                @test NINJA.sYlm(-2, 2, 0, ι, ϕ) ≈ NINJA.m2Y20(ι, ϕ) atol=tol rtol=tol
                @test NINJA.sYlm(-2, 2, -1, ι, ϕ) ≈ NINJA.m2Y2m1(ι, ϕ) atol=tol rtol=tol
                @test NINJA.sYlm(-2, 2, -2, ι, ϕ) ≈ NINJA.m2Y2m2(ι, ϕ) atol=tol rtol=tol
            end
        end
    end

    @testset "Compare to NINJA expressions ($T)" for T in [Float64, Float32, BigFloat]
        ℓₘₐₓ = 8
        sₘₐₓ = 2
        ℓₘᵢₙ = 0
        tol = ℓₘₐₓ^2 * 2eps(T)  # Mostly because the NINJA.sYlm expressions are inaccurate
        Y, H_rec_coeffs, Hwedge, expimϕ = Yprep(ℓₘₐₓ, sₘₐₓ, T)

        let R = randn(Rotor{T})
            @test_throws ErrorException Y!(Y, R, ℓₘₐₓ+1, 0, H_rec_coeffs, Hwedge, expimϕ)
            ℓ = ℓₘₐₓ - 1
            i = WignerHsize(ℓ, abs(0)) - 1
            @test_throws ErrorException Y!(Y, R, ℓ, 0, H_rec_coeffs, Hwedge[1:i], expimϕ)
        end

        @showprogress "Compare to NINJA expressions ($T)" for spin in -sₘₐₓ:sₘₐₓ
            for ι in βrange(T)
                for ϕ in αrange(T)
                    R = from_spherical_coordinates(ι, ϕ)
                    Y!(Y, R, ℓₘₐₓ, spin, H_rec_coeffs, Hwedge, expimϕ)
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
                            i += 1
                        end
                    end
                end
            end
        end
    end

    @testset "Spin property ($T)" for T in [Float64, Float32, BigFloat]
        # Test that ₛYₗₘ(R exp(γ*z/2)) = ₛYₗₘ(R) * exp(-im*s*γ)
        # See https://spherical.readthedocs.io/en/main/SWSHs/
        # for a more detailed explanation
        ℓₘₐₓ = 8
        sₘₐₓ = 2
        ℓₘᵢₙ = 0
        tol = 4ℓₘₐₓ * eps(T)
        Y1, H_rec_coeffs, Hwedge, expimϕ = Yprep(ℓₘₐₓ, sₘₐₓ, T)
        Y2 = Ystorage(ℓₘₐₓ, T)
        @showprogress "Spin property ($T)" for spin in -sₘₐₓ:sₘₐₓ
            for ι in βrange(T)
                for ϕ in αrange(T)
                    for γ in γrange(T)
                        R = from_spherical_coordinates(ι, ϕ)
                        Y!(Y1, R * exp(γ*𝐤/2), ℓₘₐₓ, spin, H_rec_coeffs, Hwedge, expimϕ)
                        Y!(Y2, R, ℓₘₐₓ, spin, H_rec_coeffs, Hwedge, expimϕ)
                        @test Y1 ≈ Y2 * cis(-spin*γ) atol=tol rtol=tol
                    end
                end
            end
        end
    end

    @testset "sYlm vs WignerD ($T)" for T in [Float64, Float32, BigFloat]
        # ₛYₗₘ(R) = (-1)ˢ √((2ℓ+1)/(4π)) 𝔇ˡₘ₋ₛ(R)
        #        = (-1)ˢ √((2ℓ+1)/(4π)) 𝔇̄ˡ₋ₛₘ(R̄)
        ℓₘₐₓ = 8
        sₘₐₓ = 2
        ℓₘᵢₙ = 0
        tol = 4ℓₘₐₓ * eps(T)
        𝔇, _, eⁱᵐᵅ, eⁱᵐᵞ = Dprep(ℓₘₐₓ, T)
        Y, H_rec_coeffs, Hwedge, expimϕ = Yprep(ℓₘₐₓ, sₘₐₓ, T)
        @showprogress "sYlm vs WignerD ($T)" for s in -sₘₐₓ:sₘₐₓ
            for ι in βrange(T)
                for ϕ in αrange(T)
                    R = from_spherical_coordinates(ι, ϕ)
                    D!(𝔇, R, ℓₘₐₓ, H_rec_coeffs, eⁱᵐᵅ, eⁱᵐᵞ)
                    Y!(Y, R, ℓₘₐₓ, s, H_rec_coeffs, Hwedge, expimϕ)
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

    @testset "sYlm conjugation ($T)" for T in [Float64, Float32, BigFloat]
        # ₛȲₗₘ = (-1)ˢ⁺ᵐ ₋ₛYₗ₋ₘ
        ℓₘₐₓ = 8
        sₘₐₓ = 2
        ℓₘᵢₙ = 0
        tol = 4ℓₘₐₓ * eps(T)
        Y1, H_rec_coeffs, Hwedge, expimϕ = Yprep(ℓₘₐₓ, sₘₐₓ, T)
        Y2 = Ystorage(ℓₘₐₓ, T)
        @showprogress "sYlm conjugation ($T)" for ι in βrange(T)
            for ϕ in αrange(T)
                for γ in γrange(T)
                    for s in -sₘₐₓ:sₘₐₓ
                        R = from_spherical_coordinates(ι, ϕ)
                        Y!(Y1, R, ℓₘₐₓ, s, H_rec_coeffs, Hwedge, expimϕ)
                        Y!(Y2, R, ℓₘₐₓ, -s, H_rec_coeffs, Hwedge, expimϕ)
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
