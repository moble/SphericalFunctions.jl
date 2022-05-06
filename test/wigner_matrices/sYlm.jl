@testset verbose=true "sYlm" begin

    @testset "Test NINJA expressions ($T)" for T in [Float64, Float32, BigFloat]
        ## This is just to test my implementation of the equations give in the paper.
        ## Note that this is a test of the testing code itself, not of the main code.
        tol = 2eps(T)
        for ι in βrange(T)
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
        tol = 100eps(T)
        ℓₘₐₓ = 8
        sₘₐₓ = 2
        ℓₘᵢₙ = 0
        Y, H_rec_coeffs, Hwedge, expimϕ = Yprep(ℓₘₐₓ, sₘₐₓ, T)
        for spin in -sₘₐₓ:sₘₐₓ
            for ι in βrange(T)
                for ϕ in αrange(T)
                    R = from_spherical_coordinates(ι, ϕ)
                    Y!(Y, R, ℓₘₐₓ, spin, H_rec_coeffs, Hwedge, expimϕ)
                    i = 1
                    for ℓ in 0:abs(spin)-1
                        for m in -ℓ:ℓ
                            if Y[i] != 0
                                println("Nonzero at i=$i (s,ℓ,m)=$((spin, ℓ,m))")
                            end
                            @test Y[i] == 0
                            i += 1
                        end
                    end
                    for ℓ in abs(spin):ℓₘₐₓ
                        for m in -ℓ:ℓ
                            sYlm1 = Y[i]
                            sYlm2 = NINJA.sYlm(spin, ℓ, m, ι, ϕ)
                            if ≉(sYlm1, sYlm2, atol=tol, rtol=tol)
                                println("Unequal at i=$i (s,ℓ,m)=$((spin,ℓ,m)): ")
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
