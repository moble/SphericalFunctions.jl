@testset verbose=true "Operators" begin

    @testset "Casimir $T" for T ∈ [Float32, Float64, Double64, BigFloat]
        ϵ = 100 * eps(T)
        for s ∈ -3:3
            for ℓₘₐₓ ∈ 4:7
                for ℓₘᵢₙ ∈ 0:min(abs(s)+1, ℓₘₐₓ)
                    let L²=L²(s, ℓₘᵢₙ, ℓₘₐₓ, T),
                        Lz=Lz(s, ℓₘᵢₙ, ℓₘₐₓ, T),
                        L₊=L₊(s, ℓₘᵢₙ, ℓₘₐₓ, T),
                        L₋=L₋(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                        L1 = L²
                        L2 = (L₊*L₋ .+ L₋*L₊ .+ 2Lz*Lz)/2
                        @test L1 ≈ L2 atol=ϵ rtol=ϵ
                    end
                    let L²=L²(s, ℓₘᵢₙ, ℓₘₐₓ, T),
                        R²=R²(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                        @test L² ≈ R² atol=ϵ rtol=ϵ
                    end
                    let
                        # R² = (2Rz² + R₊R₋ + R₋R₊)/2
                        R1 = R²(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                        R2 = T.(Array(
                            R₊(s+1, ℓₘᵢₙ, ℓₘₐₓ, T) * R₋(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            .+ R₋(s-1, ℓₘᵢₙ, ℓₘₐₓ, T) * R₊(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            .+ 2Rz(s, ℓₘᵢₙ, ℓₘₐₓ, T) * Rz(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                        ) / 2)
                        @test R1 ≈ R2 atol=ϵ rtol=ϵ
                    end
                end
            end
        end
    end

    @testset verbose=false "Formula $T" for T ∈ [Float32, Float64, Double64, BigFloat]
        ϵ = 100 * eps(T)
        @testset "$ℓₘₐₓ" for ℓₘₐₓ ∈ 4:7
            for s in -3:3
                let ℓₘᵢₙ = 0
                    𝒯₊ = SSHT(s+1, ℓₘₐₓ; T=T, method="Direct", inplace=false)
                    𝒯₋ = SSHT(s-1, ℓₘₐₓ; T=T, method="Direct", inplace=false)
                    i₊ = Yindex(abs(s+1), -abs(s+1), ℓₘᵢₙ)
                    i₋ = Yindex(abs(s-1), -abs(s-1), ℓₘᵢₙ)
                    Y = zeros(Complex{T}, Ysize(ℓₘᵢₙ, ℓₘₐₓ))
                    for ℓ in abs(s):ℓₘₐₓ
                        for m in -ℓ:ℓ
                            Y[:] .= zero(T)
                            Y[Yindex(ℓ, m, ℓₘᵢₙ)] = one(T)
                            ðY = 𝒯₊ * (ð(s, ℓₘᵢₙ, ℓₘₐₓ, T) * Y)[i₊:end]
                            Y₊ = 𝒯₊ * Y[i₊:end]
                            c₊ = ℓ < abs(s+1) ? zero(T) : √T((ℓ-s)*(ℓ+s+1))
                            @test ðY ≈ c₊ * Y₊ atol=ϵ rtol=ϵ
                            ð̄Y = 𝒯₋ * (ð̄(s, ℓₘᵢₙ, ℓₘₐₓ, T) * Y)[i₋:end]
                            Y₋ = 𝒯₋ * Y[i₋:end]
                            c₋ = ℓ < abs(s-1) ? zero(T) : -√T((ℓ+s)*(ℓ-s+1))
                            @test ð̄Y ≈ c₋ * Y₋ atol=ϵ rtol=ϵ
                        end
                    end
                end
            end
        end
    end

    @testset verbose=false "Commutators $T" for T ∈ [Float32, Float64, Double64, BigFloat]
        # [L², Lz] = 0     [L², L₊] = 0     [L², L₋] = 0
        # [R², Rz] = 0     [R², R₊] = 0     [R², R₋] = 0
        # [Lz, L₊] = L₊    [Lz, L₋] = -L₋   [L₊, L₋] = 2Lz
        # [Rz, R₊] = R₊    [Rz, R₋] = -R₋   [R₊, R₋] = 2Rz
        # [Rz, ð] = -ð     [Rz, ð̄] = ð̄      [ð, ð̄] = 2Rz
        ϵ = 100 * eps(T)
        @testset "$ℓₘₐₓ" for ℓₘₐₓ ∈ 4:7
            for s in -3:3
                let ℓₘᵢₙ = 0
                    Y = zeros(Complex{T}, Ysize(ℓₘᵢₙ, ℓₘₐₓ))
                    for Oᵢ ∈ [Lz, L₊, L₋, Rz, R₊, R₋]
                        for O² ∈ [L², R²]
                            let O²=O²(s, ℓₘᵢₙ, ℓₘₐₓ, T),
                                Oᵢ=Oᵢ(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                                # [O², Oᵢ] = 0
                                @test O²*Oᵢ-Oᵢ*O² ≈ 0*O² atol=ϵ rtol=ϵ
                            end
                        end
                    end
                    let Lz=Array(Lz(s, ℓₘᵢₙ, ℓₘₐₓ, T)),
                        L₊=Array(L₊(s, ℓₘᵢₙ, ℓₘₐₓ, T)),
                        L₋=Array(L₋(s, ℓₘᵢₙ, ℓₘₐₓ, T))
                        # [Lz, L₊] = L₊
                        @test Lz*L₊ - L₊*Lz ≈ L₊ atol=ϵ rtol=ϵ
                        # [Lz, L₋] = -L₋
                        @test Lz*L₋ - L₋*Lz ≈ -L₋ atol=ϵ rtol=ϵ
                        # [L₊, L₋] = 2Lz
                        @test L₊*L₋ - L₋*L₊ ≈ 2Lz atol=ϵ rtol=ϵ
                    end
                    let
                        # [Rz, R₊] = R₊
                        @test (
                            Rz(s-1, ℓₘᵢₙ, ℓₘₐₓ, T)*R₊(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            - R₊(s, ℓₘᵢₙ, ℓₘₐₓ, T)*Rz(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            ≈ R₊(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                        ) atol=ϵ rtol=ϵ
                        # [Rz, R₋] = -R₋
                        @test (
                            Rz(s+1, ℓₘᵢₙ, ℓₘₐₓ, T)*R₋(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            - R₋(s, ℓₘᵢₙ, ℓₘₐₓ, T)*Rz(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            ≈ -R₋(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                        ) atol=ϵ rtol=ϵ
                        # [R₊, R₋] = 2Rz
                        @test (
                            R₊(s+1, ℓₘᵢₙ, ℓₘₐₓ, T)*R₋(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            - R₋(s-1, ℓₘᵢₙ, ℓₘₐₓ, T)*R₊(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            ≈ 2Rz(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                        ) atol=ϵ rtol=ϵ
                        # [Rz, ð] = -ð
                        @test (
                            Rz(s+1, ℓₘᵢₙ, ℓₘₐₓ, T)*ð(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            - ð(s, ℓₘᵢₙ, ℓₘₐₓ, T)*Rz(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            ≈ -ð(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                        ) atol=ϵ rtol=ϵ
                        # [Rz, ð̄] = ð̄
                        @test (
                            Rz(s-1, ℓₘᵢₙ, ℓₘₐₓ, T)*ð̄(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            - ð̄(s, ℓₘᵢₙ, ℓₘₐₓ, T)*Rz(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            ≈ ð̄(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                        ) atol=ϵ rtol=ϵ
                        # [ð, ð̄] = 2Rz
                        @test (
                            ð(s-1, ℓₘᵢₙ, ℓₘₐₓ, T)*ð̄(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            -ð̄(s+1, ℓₘᵢₙ, ℓₘₐₓ, T)*ð(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            ≈ 2Rz(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                        ) atol=ϵ rtol=ϵ
                    end
                end
            end
        end
    end

end
