@testitem "Pretest ε and basis commutators" setup=[Utilities] begin
    using Quaternionic
    # Test that [eⱼ, eₖ] = 2∑ₗ ε(j,k,l) eₗ
    let e = [imx, imy, imz]
        for (j,eⱼ) ∈ enumerate(e)
            for (k,eₖ) ∈ enumerate(e)
                @test eⱼ*eₖ - eₖ*eⱼ == 2sum(ε(j,k,l)*e[l] for l ∈ 1:3)
            end
        end
    end
end

@testitem "Explicit definition" setup=[ExplicitOperators] begin
    using Quaternionic
    using DoubleFloats
    using Random
    Random.seed!(123)
    const L = ExplicitOperators.L
    const R = ExplicitOperators.R
    for T ∈ [Float32, Float64, Double64, BigFloat]
        # Test the `L` and `R` operators as defined above compared to eigenvalues on 𝔇
        ϵ = 100 * eps(T)
        for Q ∈ randn(Rotor{T}, 10)
            for ℓ ∈ 0:4
                for m ∈ -ℓ:ℓ
                    for m′ ∈ -ℓ:ℓ
                        f(Q) = D_matrices(Q, ℓ)[WignerDindex(ℓ, m, m′)]

                        @test R(imz, f)(Q) ≈ m′ * f(Q) atol=ϵ rtol=ϵ
                        @test L(imz, f)(Q) ≈ m * f(Q) atol=ϵ rtol=ϵ

                        if ℓ ≥ abs(m+1)
                            L₊1 = L(imx, f)(Q) + im * L(imy, f)(Q)
                            L₊2 = √T((ℓ-m)*(ℓ+m+1)) * D_matrices(Q, ℓ)[WignerDindex(ℓ, m+1, m′)]
                            @test L₊1 ≈ L₊2 atol=ϵ rtol=ϵ
                        end

                        if ℓ ≥ abs(m-1)
                            L₋1 = L(imx, f)(Q) - im * L(imy, f)(Q)
                            L₋2 = √T((ℓ+m)*(ℓ-m+1)) * D_matrices(Q, ℓ)[WignerDindex(ℓ, m-1, m′)]
                            @test L₋1 ≈ L₋2 atol=ϵ rtol=ϵ
                        end

                        if ℓ ≥ abs(m′+1)
                            K₊1 = R(imx, f)(Q) - im * R(imy, f)(Q)
                            K₊2 = √T((ℓ-m′)*(ℓ+m′+1)) * D_matrices(Q, ℓ)[WignerDindex(ℓ, m, m′+1)]
                            @test K₊1 ≈ K₊2 atol=ϵ rtol=ϵ
                        end

                        if ℓ ≥ abs(m′-1)
                            K₋1 = R(imx, f)(Q) + im * R(imy, f)(Q)
                            K₋2 = √T((ℓ+m′)*(ℓ-m′+1)) * D_matrices(Q, ℓ)[WignerDindex(ℓ, m, m′-1)]
                            @test K₋1 ≈ K₋2 atol=ϵ rtol=ϵ
                        end
                    end
                end
            end
        end
    end
end

@testitem "Composition" setup=[ExplicitOperators] begin
    # Test the order of operations:
    #   LₘLₙf(Q) = λ²∂ᵧ∂ᵨf(exp(ρn) exp(γm) Q)
    #   RₘRₙf(Q) = λ²∂ᵧ∂ᵨf(Q exp(γm) exp(ρn))
    using Quaternionic
    using DoubleFloats
    import ForwardDiff
    using Random
    Random.seed!(123)

    const L = ExplicitOperators.L
    const R = ExplicitOperators.R

    for T ∈ [Float32, Float64, Double64]
        z = zero(T)
        function LL(m, n, f, Q)
            - ForwardDiff.derivative(
                γ -> ForwardDiff.derivative(
                    ρ -> f((cos(ρ) + sin(ρ)*n) * (cos(γ) + sin(γ)*m) * Q),
                    z
                ),
                z
            ) / 4
        end
        function RR(m, n, f, Q)
            - ForwardDiff.derivative(
                γ -> ForwardDiff.derivative(
                    ρ -> f(Q * (cos(γ) + sin(γ)*m) * (cos(ρ) + sin(ρ)*n)),
                    z
                ),
                z
            ) / 4
        end

        ϵ = 100 * eps(T)
        M = randn(QuatVec{T}, 5)
        N = randn(QuatVec{T}, 5)
        for Q ∈ randn(Rotor{T}, 10)
            for ℓ ∈ 0:4
                for m ∈ -ℓ:ℓ
                    for m′ ∈ -ℓ:ℓ
                        f(Q) = D_matrices(Q, ℓ)[WignerDindex(ℓ, m, m′)]
                        for n ∈ N
                            for m ∈ M
                                @test L(m, L(n, f))(Q) ≈ LL(m, n, f, Q) atol=ϵ rtol=ϵ
                                @test R(m, R(n, f))(Q) ≈ RR(m, n, f, Q) atol=ϵ rtol=ϵ
                            end
                        end
                    end
                end
            end
        end
    end
end

@testitem "Scalar multiplication" setup=[ExplicitOperators] begin
    using Quaternionic
    using DoubleFloats
    const L = ExplicitOperators.L
    const R = ExplicitOperators.R
    for T ∈ [Float32, Float64, Double64]
        # Test L_{sg} = sL_{g} and R_{sg} = sR_{g}
        ϵ = 100 * eps(T)
        Ss = randn(T, 5)
        Gs = randn(QuatVec{T}, 5)
        for Q ∈ randn(Rotor{T}, 5)
            for ℓ ∈ 0:4
                for m ∈ -ℓ:ℓ
                    for m′ ∈ -ℓ:ℓ
                        f(Q) = D_matrices(Q, ℓ)[WignerDindex(ℓ, m, m′)]
                        for s ∈ Ss
                            for g ∈ Gs
                                @test L(s*g, f)(Q) ≈ s*L(g, f)(Q) atol=ϵ rtol=ϵ
                                @test R(s*g, f)(Q) ≈ s*R(g, f)(Q) atol=ϵ rtol=ϵ
                            end
                        end
                    end
                end
            end
        end
    end
end

@testitem "Additivity" setup=[ExplicitOperators] begin
    using Quaternionic
    using DoubleFloats
    const L = ExplicitOperators.L
    const R = ExplicitOperators.R
    for T ∈ [Float32, Float64, Double64]
        # Test L_{a+b} = L_{a}+L_{b} and R_{a+b} = R_{a}+R_{b}
        ϵ = 100 * eps(T)
        Gs = randn(QuatVec{T}, 5)
        for Q ∈ randn(Rotor{T}, 5)
            for ℓ ∈ 0:4
                for m ∈ -ℓ:ℓ
                    for m′ ∈ -ℓ:ℓ
                        f(Q) = D_matrices(Q, ℓ)[WignerDindex(ℓ, m, m′)]
                        for g₁ ∈ Gs
                            for g₂ ∈ Gs
                                @test L(g₁+g₂, f)(Q) ≈ L(g₁, f)(Q) + L(g₂, f)(Q) atol=ϵ rtol=ϵ
                                @test R(g₁+g₂, f)(Q) ≈ R(g₁, f)(Q) + R(g₂, f)(Q) atol=ϵ rtol=ϵ
                            end
                        end
                    end
                end
            end
        end
    end
end

@testitem "Basis commutators" setup=[ExplicitOperators] begin
    # [Lⱼ, Lₖ] =  im L_{[eⱼ,eₖ]/2} =  im ∑ₗ ε(j,k,l) Lₗ
    # [Rⱼ, Rₖ] = -im R_{[eⱼ,eₖ]/2} = -im ∑ₗ ε(j,k,l) Rₗ
    # [Lⱼ, Rₖ] = 0
    using Quaternionic
    using DoubleFloats
    import ForwardDiff
    using Random
    Random.seed!(1234)

    const L = ExplicitOperators.L
    const R = ExplicitOperators.R

    for T ∈ [Float32, Float64, Double64]
        ϵ = 400 * eps(T)
        E = QuatVec{T}[imx, imy, imz]
        for Q ∈ randn(Rotor{T}, 10)
            for ℓ ∈ 0:4
                for m ∈ -ℓ:ℓ
                    for m′ ∈ -ℓ:ℓ
                        f(Q) = D_matrices(Q, ℓ)[WignerDindex(ℓ, m, m′)]
                        for eⱼ ∈ E
                            for eₖ ∈ E
                                eⱼeₖ = QuatVec{T}(eⱼ * eₖ - eₖ * eⱼ) / 2
                                @test L(eⱼ, L(eₖ, f))(Q) - L(eₖ, L(eⱼ, f))(Q) ≈ im * L(eⱼeₖ, f)(Q) atol=ϵ rtol=ϵ
                                @test R(eⱼ, R(eₖ, f))(Q) - R(eₖ, R(eⱼ, f))(Q) ≈ -im * R(eⱼeₖ, f)(Q) atol=ϵ rtol=ϵ
                                @test L(eⱼ, R(eₖ, f))(Q) - R(eₖ, L(eⱼ, f))(Q) ≈ zero(T) atol=4ϵ
                            end
                        end
                    end
                end
            end
        end
    end
end

@testitem "Commutators" begin
    using DoubleFloats
    for T ∈ [Float32, Float64, Double64, BigFloat]
        # Test the following relations:
        # [L², Lz] = 0     [L², L₊] = 0     [L², L₋] = 0
        # [R², Rz] = 0     [R², R₊] = 0     [R², R₋] = 0
        # [Lz, L₊] = L₊    [Lz, L₋] = -L₋   [L₊, L₋] = 2Lz
        # [Rz, R₊] = R₊    [Rz, R₋] = -R₋   [R₊, R₋] = 2Rz
        # [Rz, ð] = -ð     [Rz, ð̄] = ð̄      [ð, ð̄] = 2Rz
        ϵ = 100 * eps(T)
        @testset "$ℓₘₐₓ" for ℓₘₐₓ ∈ 4:7
            for s in -3:3
                let ℓₘᵢₙ = 0
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

@testitem "Casimir" begin
    using DoubleFloats
    for T ∈ [Float32, Float64, Double64, BigFloat]
        # Test that L² = (L₊L₋ + L₋L₊ + 2Lz²)/2 = R² = (R₊R₋ + R₋R₊ + 2Rz²)/2
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
end

@testitem "Applied to ₛYₗₘ" begin
    using DoubleFloats
    for T ∈ [Float32, Float64, Double64, BigFloat]
        # Evaluate (on points) ðY = √((ℓ-s)(ℓ+s+1)) Y, and similarly for ð̄Y
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
end

## TODO: Add L_x, L_y, R_x, and R_y, then test these commutators.
## Note that R is harder because the basis in which all the matrices are returned
## assumes that you are dealing with a particular `s` eigenvalue.
# [Lⱼ, Lₖ] =  im L_{[eⱼ,eₖ]/2} =  im ∑ₗ ε(j,k,l) Lₗ
# [Rⱼ, Rₖ] = -im R_{[eⱼ,eₖ]/2} = -im ∑ₗ ε(j,k,l) Rₗ
# [Lⱼ, Rₖ] = 0
