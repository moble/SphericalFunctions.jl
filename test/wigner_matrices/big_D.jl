@testitem "Compare H/D indexing" setup=[NaNChecker] begin
    const NaNCheck = NaNChecker.NaNCheck
    @testset "$T" for T in [Float64, Float32]
        # Here, we check that we can pass in either an "H wedge" array to be used with
        # WignerHindex, or a full 𝔇 array used with WignerDindex, and obtain the same
        # H recurrence results
        ℓₘₐₓ = 8
        for m′ₘₐₓ in 0:ℓₘₐₓ
            expiβ = cis(rand(0:eps(T):π))
            expiβNaNCheck = complex(NaNCheck{T}(expiβ.re), NaNCheck{T}(expiβ.im))
            NCTN = NaNCheck{T}(NaN)
            Hw = fill(NCTN, WignerHsize(ℓₘₐₓ, m′ₘₐₓ))
            H!(Hw, expiβNaNCheck, ℓₘₐₓ, m′ₘₐₓ, H_recursion_coefficients(ℓₘₐₓ, T))
            𝔇 = fill(NCTN, WignerDsize(ℓₘₐₓ, m′ₘₐₓ))
            H!(
                𝔇, expiβNaNCheck, ℓₘₐₓ, m′ₘₐₓ,
                H_recursion_coefficients(ℓₘₐₓ, T), WignerDindex
            )
            for n in 0:ℓₘₐₓ
                for m′ in -min(n, m′ₘₐₓ):min(n, m′ₘₐₓ)
                    for m in abs(m′):n
                        Hnm′m = Hw[WignerHindex(n, m′, m, m′ₘₐₓ)]
                        𝔇nm′m = 𝔇[WignerDindex(n, m′, m, m′ₘₐₓ)]
                        @test Hnm′m == 𝔇nm′m
                    end
                end
            end
        end
    end
end

@testitem "Compare 𝔇 to formulaic d" setup=[ExplicitWignerMatrices,Utilities] begin
    using Quaternionic
    @testset "$T" for T in [BigFloat, Float64, Float32]
        # Now, we're ready to check that d_{n}^{m′,m}(β) matches the expected values
        # for a range of β values
        for ℓₘₐₓ in 0:4
            D_storage = D_prep(ℓₘₐₓ, T)
            expiα = complex(one(T))
            expiγ = complex(one(T))
            for β in βrange(T)
                expiβ = cis(β)
                R = from_euler_angles(zero(T), β, zero(T))
                𝔇 = D_matrices!(D_storage, R)
                for n in 0:ℓₘₐₓ
                    for m′ in -n:n
                        for m in -n:n
                            𝔇_formula = ExplicitWignerMatrices.D_formula(
                                n, m′, m, expiα, expiβ, expiγ
                            )
                            𝔇_recurrence = 𝔇[WignerDindex(n, m′, m)]
                            @test conj(𝔇_formula) ≈ 𝔇_recurrence atol=200eps(T) rtol=200eps(T)
                        end
                    end
                end
            end
        end
    end
end

@testitem "Compare 𝔇 to formulaic 𝔇" setup=[ExplicitWignerMatrices,Utilities] begin
    using Quaternionic
    using ProgressMeter
    using Random
    @testset "$T" for T in [BigFloat, Float64, Float32]
        # Now, we're ready to check that 𝔇_{n}^{m′,m}(β) matches the expected values
        # for a range of α, β, γ values
        Random.seed!(123)
        ℓₘₐₓ = T===BigFloat ? 4 : 8
        D_storage = D_prep(ℓₘₐₓ, T)
        @showprogress desc="Compare 𝔇 to formulaic 𝔇 ($T)" for α in αrange(T, 5)
            for β in βrange(T, 5)
                for γ in γrange(T, 5)
                    R = from_euler_angles(α, β, γ)
                    expiα, expiβ, expiγ = to_euler_phases(R)
                    𝔇 = D_matrices!(D_storage, R)
                    for n in 0:ℓₘₐₓ
                        for m′ in -n:n
                            for m in -n:n
                                𝔇_formula = ExplicitWignerMatrices.D_formula(
                                    n, m′, m, expiα, expiβ, expiγ
                                )
                                𝔇_recurrence = 𝔇[WignerDindex(n, m′, m)]
                                @test conj(𝔇_formula) ≈ 𝔇_recurrence atol=400eps(T) rtol=400eps(T)
                            end
                        end
                    end
                end
            end
        end
    end
end

@testitem "Group characters" setup=[Utilities] begin
    using Quaternionic
    using ProgressMeter
    using Random
    @testset "$T" for T in [BigFloat, Float64, Float32]
        # χʲ(β) ≔ Σₘ dʲₘₘ(β) ≡ Σₘ 𝔇ʲₘₘ(exp(v̂ β/2)) = sin((2j+1)β/2) / sin(β/2)
        # Here, v̂ is any unit vector; group characters are constant on conjugacy classes and
        # conjugacy classes of SO(3) are rotations through the same angle about any axis.
        Random.seed!(123)
        ℓₘₐₓ = T===BigFloat ? 10 : 20
        D_storage = D_prep(ℓₘₐₓ, T)
        d_storage = d_prep(ℓₘₐₓ, T)
        @showprogress desc="Group characters $T" for β in βrange(T)
            expiβ = cis(β)
            d = d_matrices!(d_storage, expiβ)
            for j in 0:ℓₘₐₓ
                sin_ratio = if abs(β) < 10eps(T)
                    T(2j+1)
                elseif abs(β-π) < 10eps(T)
                    T(-1)^j
                else
                    sin((2j+1)*β/2) / sin(β/2)
                end
                χʲ = sum(d[WignerDindex(j, m, m)] for m in -j:j)
                @test χʲ ≈ sin_ratio atol=500eps(T) rtol=500eps(T)
                for v̂ in v̂range(T)
                    R = exp(β/2 * v̂)
                    𝔇 = D_matrices!(D_storage, R)
                    χʲ = sum(𝔇[WignerDindex(j, m, m)] for m in -j:j)
                    @test χʲ ≈ sin_ratio atol=500eps(T) rtol=500eps(T)
                end
            end
        end
    end
end

@testitem "Representation property" setup=[Utilities] begin
    using Quaternionic
    using ProgressMeter
    using Random
    @testset "$T" for T in [Float64, Float32, BigFloat]
        # For each l, 𝔇ˡₙ,ₘ(R₁ R₂) = Σₚ 𝔇ˡₙ,ₚ(R₁) 𝔇ˡₚ,ₘ(R₂)
        Random.seed!(123)
        tol = 3eps(T)
        ℓₘₐₓ = 10
        D₁_storage = D_prep(ℓₘₐₓ, T)
        𝔇₁ = D₁_storage[1]
        𝔇₂ = similar(𝔇₁)
        D₂_storage = (𝔇₂, D₁_storage[2:end]...)
        𝔇₁₂ = similar(𝔇₁)
        D₁₂_storage = (𝔇₁₂, D₁_storage[2:end]...)
        @showprogress desc="Representation property ($T)" for R₁ in Rrange(T)
            for R₂ in Rrange(T)
                D_matrices!(D₁_storage, R₁)
                D_matrices!(D₂_storage, R₂)
                D_matrices!(D₁₂_storage, R₁*R₂)
                for ℓ in 0:ℓₘₐₓ
                    i = WignerDindex(ℓ, -ℓ, -ℓ)
                    j = WignerDindex(ℓ, ℓ, ℓ)
                    𝔇₁ˡ = transpose(reshape(𝔇₁[i:j], 2ℓ+1, 2ℓ+1))
                    𝔇₂ˡ = transpose(reshape(𝔇₂[i:j], 2ℓ+1, 2ℓ+1))
                    𝔇₁₂ˡ = transpose(reshape(𝔇₁₂[i:j], 2ℓ+1, 2ℓ+1))
                    @test 𝔇₁ˡ * 𝔇₂ˡ ≈ 𝔇₁₂ˡ atol=(2ℓ+1)^2*tol rtol=(2ℓ+1)^2*tol
                end
            end
        end
    end
end
