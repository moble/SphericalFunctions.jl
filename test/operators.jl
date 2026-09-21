# Tests of the angular-momentum operators in `src/utilities/operators.jl`, against the
# settled conventions documented in `docs/src/30-conventions/01-summary.md`:
#
#     L_𝐮 f(𝐑) =  i d/dϵ f(e^{-ϵ𝐮/2} 𝐑),      R_𝐮 f(𝐑) = -i d/dϵ f(𝐑 e^{-ϵ𝐮/2}),
#     L_± = L_x ± i L_y,   R_± = R_x ± i R_y,   [L_z, L_±] = ±L_±,   [R_z, R_±] = ±R_±,
#     L_z 𝔇ˡ_{m′m} = -m′ 𝔇ˡ_{m′m},                R_z 𝔇ˡ_{m′m} = m 𝔇ˡ_{m′m},
#     L_± 𝔇ˡ_{m′m} = -√((ℓ±m′)(ℓ∓m′+1)) 𝔇ˡ_{m′∓1,m},  R_± 𝔇ˡ_{m′m} = √((ℓ∓m)(ℓ±m+1)) 𝔇ˡ_{m′,m±1},
#     L_z ₛYₗₘ = m ₛYₗₘ,   R_z ₛYₗₘ = s ₛYₗₘ,
#     L_± ₛYₗₘ = √((ℓ∓m)(ℓ±m+1)) ₛYₗ,ₘ±₁,   R_± ₛYₗₘ = √((ℓ∓s)(ℓ±s+1)) ₛ±₁Yₗₘ,
#     ð = R₊,   ð̄ = -R₋.
#
# The explicit differential operators (`ExplicitOperators`) are applied, via automatic
# differentiation, to the package's own 𝔇 and ₛYₗₘ values, so these tests pin the sign
# conventions of the operators *and* of the Wigner/sYlm functions simultaneously.

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

@testitem "Operators: explicit definition on 𝔇" setup=[ExplicitOperators] begin
    import SphericalFunctions: D
    using Quaternionic
    using DoubleFloats
    using Random
    rng = Random.Xoshiro(123)
    const L = ExplicitOperators.L
    const R = ExplicitOperators.R
    for T ∈ [Float32, Float64, Double64, BigFloat]
        # Compare the explicit L and R operators, acting on 𝔇ˡ_{m′,m}, to the eigenvalue
        # and ladder relations of the conventions summary.
        ϵ = 100 * eps(T)
        for Q ∈ randn(rng, Rotor{T}, 6)
            for ℓ ∈ 0:4
                𝔇 = Q -> D(Q, ℓ)[ℓ]
                for m′ ∈ -ℓ:ℓ
                    for m ∈ -ℓ:ℓ
                        f(Q) = 𝔇(Q)[m′, m]

                        # L_z 𝔇_{m′m} = -m′ 𝔇_{m′m};  R_z 𝔇_{m′m} = m 𝔇_{m′m}
                        @test L(imz, f)(Q) ≈ -m′ * f(Q) atol=ϵ rtol=ϵ
                        @test R(imz, f)(Q) ≈ m * f(Q) atol=ϵ rtol=ϵ

                        # L₊ 𝔇_{m′m} = -√((ℓ+m′)(ℓ-m′+1)) 𝔇_{m′-1,m}
                        L₊f = L(imx, f)(Q) + im * L(imy, f)(Q)
                        if m′-1 ≥ -ℓ
                            @test L₊f ≈ -√T((ℓ+m′)*(ℓ-m′+1)) * 𝔇(Q)[m′-1, m] atol=ϵ rtol=ϵ
                        else
                            @test L₊f ≈ 0 atol=ϵ
                        end

                        # L₋ 𝔇_{m′m} = -√((ℓ-m′)(ℓ+m′+1)) 𝔇_{m′+1,m}
                        L₋f = L(imx, f)(Q) - im * L(imy, f)(Q)
                        if m′+1 ≤ ℓ
                            @test L₋f ≈ -√T((ℓ-m′)*(ℓ+m′+1)) * 𝔇(Q)[m′+1, m] atol=ϵ rtol=ϵ
                        else
                            @test L₋f ≈ 0 atol=ϵ
                        end

                        # R₊ 𝔇_{m′m} = √((ℓ-m)(ℓ+m+1)) 𝔇_{m′,m+1}
                        R₊f = R(imx, f)(Q) + im * R(imy, f)(Q)
                        if m+1 ≤ ℓ
                            @test R₊f ≈ √T((ℓ-m)*(ℓ+m+1)) * 𝔇(Q)[m′, m+1] atol=ϵ rtol=ϵ
                        else
                            @test R₊f ≈ 0 atol=ϵ
                        end

                        # R₋ 𝔇_{m′m} = √((ℓ+m)(ℓ-m+1)) 𝔇_{m′,m-1}
                        R₋f = R(imx, f)(Q) - im * R(imy, f)(Q)
                        if m-1 ≥ -ℓ
                            @test R₋f ≈ √T((ℓ+m)*(ℓ-m+1)) * 𝔇(Q)[m′, m-1] atol=ϵ rtol=ϵ
                        else
                            @test R₋f ≈ 0 atol=ϵ
                        end
                    end
                end
            end
        end
    end
end

@testitem "Operators: explicit definition on ₛYₗₘ" setup=[ExplicitOperators, Utilities] begin
    # The same, applied to the spin-weighted spherical harmonics, defined here by the
    # settled relation ₛYₗₘ(R) = (-1)^s √((2ℓ+1)/4π) conj(𝔇ˡₘ,₋ₛ(R)) in terms of the
    # package's own 𝔇.  This pins the R ladder operators' action on ₛYₗₘ: R₊ raises the
    # spin weight with coefficient √((ℓ-s)(ℓ+s+1)) and R₋ lowers it with √((ℓ+s)(ℓ-s+1)),
    # both positive; and it checks that the definition reproduces the closed-form ₛYₗₘ of
    # the `Utilities` snippet — the explicit sum over factorials given on the conventions
    # pages, which shares no code with the package.  The operators are applied to the whole
    # vector of ₛYₗₘ values at once (ForwardDiff differentiates vector-valued functions),
    # which keeps the runtime reasonable.
    #
    # The closed form is a function of the spherical coordinates (θ, ϕ) alone, so the natural
    # comparison points are the rotors `from_spherical_coordinates(θ, ϕ)`.  Any other rotor
    # has an extra phase: writing 𝐐 in terms of its Euler angles (α, β, γ) as
    # 𝐐 = from_spherical_coordinates(β, α) * exp(γ𝐤/2), the defining property of spin weight,
    # η(𝐐 exp(γ𝐤/2)) = exp(-isγ) η(𝐐) (conventions summary, "Spin-weighted functions"),
    # supplies it.  Both kinds of point are used — with that phase where it is needed — so
    # that the operator identities are still exercised at generic rotors.
    import SphericalFunctions: D
    using Quaternionic
    using Random
    rng = Random.Xoshiro(321)
    const L = ExplicitOperators.L
    const R = ExplicitOperators.R
    ℓₘₐₓ = 4
    T = Float64
    ϵ = 200 * eps(T)
    idx(ℓ, m) = ℓ*(ℓ+1) + m + 1  # ℓ-major, m increasing, from ℓ=0 (entries with ℓ<|s| are 0)
    function Ys(s)
        Q -> begin
            𝔇 = D(Q, ℓₘₐₓ)
            [
                ℓ < abs(s) ? zero(𝔇[0][0, 0]) : (-1)^s * √((2ℓ+1)/(4π)) * conj(𝔇[ℓ][m, -s])
                for ℓ ∈ 0:ℓₘₐₓ for m ∈ -ℓ:ℓ
            ]
        end
    end
    # The closed-form ₛYₗₘ of the `Utilities` snippet, evaluated at the rotor 𝐐 by way of its
    # Euler angles, with the spin-weight phase discussed above.
    function closed_form(s, ℓ, m, Q)
        α, β, γ = to_euler_angles(Q)
        sYlm(s, ℓ, m, β, α) * cis(-s * γ)
    end
    # γ = 0 for the first three, so the closed form applies to them with no phase at all.
    # The angles are kept away from θ ∈ {0, π} and ϕ ∈ πℤ: at those points the rotor's
    # α ± γ phases are exactly real or imaginary, and the package's `complex_powers!` (hence
    # `D`) is not differentiable there — its `√(-dc*(2+dc))` is evaluated at dc = 0, whose
    # derivative is infinite — so the ForwardDiff-based operators below would return NaN.
    Qs = [
        [from_spherical_coordinates(T(θ), T(ϕ)) for (θ, ϕ) ∈ ((0.4, 0.9), (1.0, 2.0), (2.5, -1.5))];
        randn(rng, Rotor{T}, 3)
    ]
    for Q ∈ Qs
        for s ∈ -2:2
            Y = Ys(s)(Q)
            # The definition agrees with the independent closed form.  The maximum error over
            # these points and spins measures 2.6e-15, well inside ϵ = 200eps(Float64) ≈
            # 4.4e-14.  (Accumulated into one assertion rather than one per mode.)
            @test maximum(
                abs(Y[idx(ℓ, m)] - closed_form(s, ℓ, m, Q))
                for ℓ ∈ abs(s):ℓₘₐₓ for m ∈ -ℓ:ℓ
            ) < ϵ
            Y₊ = abs(s+1) ≤ ℓₘₐₓ ? Ys(s+1)(Q) : zero(Y)
            Y₋ = abs(s-1) ≤ ℓₘₐₓ ? Ys(s-1)(Q) : zero(Y)
            LzY = L(imz, Ys(s))(Q)
            RzY = R(imz, Ys(s))(Q)
            L₊Y = L(imx, Ys(s))(Q) + im * L(imy, Ys(s))(Q)
            L₋Y = L(imx, Ys(s))(Q) - im * L(imy, Ys(s))(Q)
            R₊Y = R(imx, Ys(s))(Q) + im * R(imy, Ys(s))(Q)
            R₋Y = R(imx, Ys(s))(Q) - im * R(imy, Ys(s))(Q)
            for ℓ ∈ abs(s):ℓₘₐₓ
                for m ∈ -ℓ:ℓ
                    # L_z ₛYₗₘ = m ₛYₗₘ;  R_z ₛYₗₘ = s ₛYₗₘ
                    @test LzY[idx(ℓ, m)] ≈ m * Y[idx(ℓ, m)] atol=ϵ rtol=ϵ
                    @test RzY[idx(ℓ, m)] ≈ s * Y[idx(ℓ, m)] atol=ϵ rtol=ϵ
                    # L_± ₛYₗₘ = √((ℓ∓m)(ℓ±m+1)) ₛYₗ,ₘ±₁
                    @test L₊Y[idx(ℓ, m)] ≈ (m+1 ≤ ℓ ? √T((ℓ-m)*(ℓ+m+1)) * Y[idx(ℓ, m+1)] : 0) atol=ϵ rtol=ϵ
                    @test L₋Y[idx(ℓ, m)] ≈ (m-1 ≥ -ℓ ? √T((ℓ+m)*(ℓ-m+1)) * Y[idx(ℓ, m-1)] : 0) atol=ϵ rtol=ϵ
                    # R_± ₛYₗₘ = √((ℓ∓s)(ℓ±s+1)) ₛ±₁Yₗₘ
                    @test R₊Y[idx(ℓ, m)] ≈ (abs(s+1) ≤ ℓ ? √T((ℓ-s)*(ℓ+s+1)) * Y₊[idx(ℓ, m)] : 0) atol=ϵ rtol=ϵ
                    @test R₋Y[idx(ℓ, m)] ≈ (abs(s-1) ≤ ℓ ? √T((ℓ+s)*(ℓ-s+1)) * Y₋[idx(ℓ, m)] : 0) atol=ϵ rtol=ϵ
                end
            end
        end
    end
end

@testitem "Operators: composition" setup=[ExplicitOperators] begin
    # Test the order of operations:
    #   LₘLₙf(Q) = λ²∂ᵧ∂ᵨf(exp(ρn) exp(γm) Q)
    #   RₘRₙf(Q) = λ²∂ᵧ∂ᵨf(Q exp(γm) exp(ρn))
    import SphericalFunctions: D
    using Quaternionic
    using DoubleFloats
    import ForwardDiff
    using Random
    rng = Random.Xoshiro(123)

    const L = ExplicitOperators.L
    const R = ExplicitOperators.R

    for T ∈ [Float32, Float64, Double64]
        z = zero(T)
        function LL(m, n, f, Q)
            # L_m L_n f = (-i/2)² ∂ᵧ∂ᵨ f(e^{ρn} e^{γm} Q)
            - ForwardDiff.derivative(
                γ -> ForwardDiff.derivative(
                    ρ -> f(exp(ρ*n) * exp(γ*m) * Q),
                    z
                ),
                z
            ) / 4
        end
        function RR(m, n, f, Q)
            # R_m R_n f = (i/2)² ∂ᵧ∂ᵨ f(Q e^{γm} e^{ρn})
            - ForwardDiff.derivative(
                γ -> ForwardDiff.derivative(
                    ρ -> f(Q * exp(γ*m) * exp(ρ*n)),
                    z
                ),
                z
            ) / 4
        end

        ϵ = 100 * eps(T)
        M = randn(rng, QuatVec{T}, 3)
        N = randn(rng, QuatVec{T}, 3)
        for Q ∈ randn(rng, Rotor{T}, 5)
            for ℓ ∈ 0:3
                for m′ ∈ -ℓ:ℓ
                    for m ∈ -ℓ:ℓ
                        f(Q) = D(Q, ℓ)[ℓ][m′, m]
                        for n ∈ N
                            for mm ∈ M
                                @test L(mm, L(n, f))(Q) ≈ LL(mm, n, f, Q) atol=ϵ rtol=ϵ
                                @test R(mm, R(n, f))(Q) ≈ RR(mm, n, f, Q) atol=ϵ rtol=ϵ
                            end
                        end
                    end
                end
            end
        end
    end
end

@testitem "Operators: linearity" setup=[ExplicitOperators] begin
    import SphericalFunctions: D
    using Quaternionic
    using DoubleFloats
    using Random
    rng = Random.Xoshiro(123)
    const L = ExplicitOperators.L
    const R = ExplicitOperators.R
    for T ∈ [Float32, Float64, Double64]
        # Test L_{sg} = sL_{g}, R_{sg} = sR_{g}, L_{a+b} = L_{a}+L_{b}, and R_{a+b} = R_{a}+R_{b}
        ϵ = 100 * eps(T)
        Ss = randn(rng, T, 3)
        Gs = randn(rng, QuatVec{T}, 3)
        for Q ∈ randn(rng, Rotor{T}, 4)
            for ℓ ∈ 0:3
                for m′ ∈ -ℓ:ℓ
                    for m ∈ -ℓ:ℓ
                        f(Q) = D(Q, ℓ)[ℓ][m′, m]
                        for s ∈ Ss
                            for g ∈ Gs
                                @test L(s*g, f)(Q) ≈ s*L(g, f)(Q) atol=ϵ rtol=ϵ
                                @test R(s*g, f)(Q) ≈ s*R(g, f)(Q) atol=ϵ rtol=ϵ
                            end
                        end
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

@testitem "Operators: basis commutators" setup=[ExplicitOperators] begin
    # [L_𝐮, L_𝐯] = (i/2) L_{[𝐮,𝐯]},   [R_𝐮, R_𝐯] = (i/2) R_{[𝐮,𝐯]},   [L_𝐮, R_𝐯] = 0
    import SphericalFunctions: D
    using Quaternionic
    using DoubleFloats
    using Random
    rng = Random.Xoshiro(1234)

    const L = ExplicitOperators.L
    const R = ExplicitOperators.R

    for T ∈ [Float32, Float64, Double64]
        ϵ = 400 * eps(T)
        E = QuatVec{T}[imx, imy, imz]
        for Q ∈ randn(rng, Rotor{T}, 5)
            for ℓ ∈ 0:3
                for m′ ∈ -ℓ:ℓ
                    for m ∈ -ℓ:ℓ
                        f(Q) = D(Q, ℓ)[ℓ][m′, m]
                        for eⱼ ∈ E
                            for eₖ ∈ E
                                eⱼeₖ = QuatVec{T}(eⱼ * eₖ - eₖ * eⱼ) / 2
                                @test L(eⱼ, L(eₖ, f))(Q) - L(eₖ, L(eⱼ, f))(Q) ≈ im * L(eⱼeₖ, f)(Q) atol=ϵ rtol=ϵ
                                @test R(eⱼ, R(eₖ, f))(Q) - R(eₖ, R(eⱼ, f))(Q) ≈ im * R(eⱼeₖ, f)(Q) atol=ϵ rtol=ϵ
                                @test L(eⱼ, R(eₖ, f))(Q) - R(eₖ, L(eⱼ, f))(Q) ≈ zero(T) atol=4ϵ
                            end
                        end
                    end
                end
            end
        end
    end
end

@testitem "Operators: matrix commutators" begin
    import SphericalFunctions: L², Lz, L₊, L₋, R², Rz, R₊, R₋, ð, ð̄
    using DoubleFloats
    for T ∈ [Float32, Float64, Double64, BigFloat]
        # Test the following relations, as matrices on mode weights.  Note that the R
        # operators change the spin weight, so the operator for the appropriate input spin
        # weight must be used in each factor:
        # [L², Lz] = 0     [L², L₊] = 0     [L², L₋] = 0
        # [R², Rz] = 0     [R², R₊] = 0     [R², R₋] = 0
        # [Lz, L₊] = L₊    [Lz, L₋] = -L₋   [L₊, L₋] = 2Lz
        # [Rz, R₊] = R₊    [Rz, R₋] = -R₋   [R₊, R₋] = 2Rz
        # [Rz, ð] = ð      [Rz, ð̄] = -ð̄    [ð, ð̄] = -2Rz
        ϵ = 100 * eps(T)
        @testset "$ℓₘₐₓ" for ℓₘₐₓ ∈ 4:7
            for s in -3:3
                let ℓₘᵢₙ = 0
                    for Oᵢ ∈ [Lz, L₊, L₋, Rz]
                        for O² ∈ [L², R²]
                            let O²=O²(s, ℓₘᵢₙ, ℓₘₐₓ, T),
                                Oᵢ=Oᵢ(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                                # [O², Oᵢ] = 0
                                @test O²*Oᵢ-Oᵢ*O² ≈ 0*O² atol=ϵ rtol=ϵ
                            end
                        end
                    end
                    for O² ∈ [L², R²]
                        # [O², R₊] = 0 and [O², R₋] = 0, with the spin-weight shift
                        @test O²(s+1, ℓₘᵢₙ, ℓₘₐₓ, T)*R₊(s, ℓₘᵢₙ, ℓₘₐₓ, T) - R₊(s, ℓₘᵢₙ, ℓₘₐₓ, T)*O²(s, ℓₘᵢₙ, ℓₘₐₓ, T) ≈ 0*R₊(s, ℓₘᵢₙ, ℓₘₐₓ, T) atol=ϵ rtol=ϵ
                        @test O²(s-1, ℓₘᵢₙ, ℓₘₐₓ, T)*R₋(s, ℓₘᵢₙ, ℓₘₐₓ, T) - R₋(s, ℓₘᵢₙ, ℓₘₐₓ, T)*O²(s, ℓₘᵢₙ, ℓₘₐₓ, T) ≈ 0*R₋(s, ℓₘᵢₙ, ℓₘₐₓ, T) atol=ϵ rtol=ϵ
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
                        # [Rz, R₊] = R₊   (R₊ maps spin weight s to s+1)
                        @test (
                            Rz(s+1, ℓₘᵢₙ, ℓₘₐₓ, T)*R₊(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            - R₊(s, ℓₘᵢₙ, ℓₘₐₓ, T)*Rz(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            ≈ R₊(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                        ) atol=ϵ rtol=ϵ
                        # [Rz, R₋] = -R₋   (R₋ maps spin weight s to s-1)
                        @test (
                            Rz(s-1, ℓₘᵢₙ, ℓₘₐₓ, T)*R₋(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            - R₋(s, ℓₘᵢₙ, ℓₘₐₓ, T)*Rz(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            ≈ -R₋(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                        ) atol=ϵ rtol=ϵ
                        # [R₊, R₋] = 2Rz
                        @test (
                            R₊(s-1, ℓₘᵢₙ, ℓₘₐₓ, T)*R₋(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            - R₋(s+1, ℓₘᵢₙ, ℓₘₐₓ, T)*R₊(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            ≈ 2Rz(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                        ) atol=ϵ rtol=ϵ
                        # [Rz, ð] = ð
                        @test (
                            Rz(s+1, ℓₘᵢₙ, ℓₘₐₓ, T)*ð(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            - ð(s, ℓₘᵢₙ, ℓₘₐₓ, T)*Rz(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            ≈ ð(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                        ) atol=ϵ rtol=ϵ
                        # [Rz, ð̄] = -ð̄
                        @test (
                            Rz(s-1, ℓₘᵢₙ, ℓₘₐₓ, T)*ð̄(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            - ð̄(s, ℓₘᵢₙ, ℓₘₐₓ, T)*Rz(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            ≈ -ð̄(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                        ) atol=ϵ rtol=ϵ
                        # [ð, ð̄] = -2Rz (which, given the two identities just below, is
                        # [R₊, R₋] = 2Rz restated)
                        @test (
                            ð(s-1, ℓₘᵢₙ, ℓₘₐₓ, T)*ð̄(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            - ð̄(s+1, ℓₘᵢₙ, ℓₘₐₓ, T)*ð(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            ≈ -2Rz(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                        ) atol=ϵ rtol=ϵ
                        # ð = R₊ and ð̄ = -R₋.  These hold by construction — `ð` *is*
                        # defined as `R₊` and `ð̄` as `-R₋` in `src/utilities/operators.jl`
                        # — so they cannot fail; they are here to pin the aliasing itself,
                        # i.e. that a future definition of `ð` in its own right would still
                        # have to agree.  The thing that actually confirms the ð sign
                        # convention against something outside the package is the
                        # Newman–Penrose finite-difference item below.
                        @test ð(s, ℓₘᵢₙ, ℓₘₐₓ, T) == R₊(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                        @test ð̄(s, ℓₘᵢₙ, ℓₘₐₓ, T) == -R₋(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                    end
                end
            end
        end
    end
end

@testitem "Operators: Casimir" begin
    import SphericalFunctions: L², Lz, L₊, L₋, R², Rz, R₊, R₋
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
                        # R² = (R₊R₋ + R₋R₊ + 2Rz²)/2, with the spin-weight shifts
                        R1 = R²(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                        R2 = T.(Array(
                            R₊(s-1, ℓₘᵢₙ, ℓₘₐₓ, T) * R₋(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            .+ R₋(s+1, ℓₘᵢₙ, ℓₘₐₓ, T) * R₊(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                            .+ 2Rz(s, ℓₘᵢₙ, ℓₘₐₓ, T) * Rz(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                        ) / 2)
                        @test R1 ≈ R2 atol=ϵ rtol=ϵ
                    end
                end
            end
        end
    end
end

@testitem "Operators: Lx and Ly" begin
    import SphericalFunctions: L², Lz, L₊, L₋, Lx, Ly, ModeWeights, Ysize, spin
    using DoubleFloats: Double64
    using Random
    rng = Random.Xoshiro(90210)
    for T ∈ (Float32, Float64, Double64, BigFloat)
        ϵ = 100eps(T)
        for s ∈ -2:2, ℓₘₐₓ ∈ 3:5, ℓₘᵢₙ ∈ unique((abs(s), 0))
            lx, ly = Matrix(Lx(s, ℓₘᵢₙ, ℓₘₐₓ, T)), Matrix(Ly(s, ℓₘᵢₙ, ℓₘₐₓ, T))
            lp, lm = Matrix(L₊(s, ℓₘᵢₙ, ℓₘₐₓ, T)), Matrix(L₋(s, ℓₘᵢₙ, ℓₘₐₓ, T))
            lz, l² = Matrix(Lz(s, ℓₘᵢₙ, ℓₘₐₓ, T)), Matrix(L²(s, ℓₘᵢₙ, ℓₘₐₓ, T))
            ## The defining combinations, which hold exactly because no arithmetic is lost
            @test lx == (lp .+ lm) ./ 2
            @test 2im .* ly == lp .- lm  # (not ly == (lp .- lm) ./ 2im: complex division
            ##                             is inexact for some float types)
            ## Real / imaginary, and Hermitian, exactly
            @test eltype(lx) === T
            @test eltype(ly) === Complex{T}
            @test all(iszero, real(ly))
            @test lx == lx'
            @test ly == ly'
            ## Casimir and the su(2) commutator
            @test lx^2 + ly^2 + lz^2 ≈ l² atol=ϵ rtol=ϵ
            @test lx*ly - ly*lx ≈ im*lz atol=ϵ rtol=ϵ
            @test ly*lz - lz*ly ≈ im*lx atol=ϵ rtol=ϵ
            @test lz*lx - lx*lz ≈ im*ly atol=ϵ rtol=ϵ
        end
    end
    ## The `ModeWeights` methods preserve the spin weight and agree with the matrices
    for s ∈ -1:1, ℓₘₐₓ ∈ (3, 4)
        w = ModeWeights(randn(rng, ComplexF64, Ysize(abs(s), ℓₘₐₓ)), s)
        for (O, M) ∈ ((Lx, Lx(s, abs(s), ℓₘₐₓ)), (Ly, Ly(s, abs(s), ℓₘₐₓ)))
            @test spin(O(w)) == s
            @test parent(O(w)) == M * parent(w)
        end
    end
    ## `ℓₘᵢₙ` defaults to `abs(s)`
    for O ∈ (Lx, Ly), s ∈ -2:2, ℓₘₐₓ ∈ 2:4
        @test O(s, ℓₘₐₓ) == O(s, abs(s), ℓₘₐₓ, Float64)
        @test O(s, ℓₘₐₓ, Float32) == O(s, abs(s), ℓₘₐₓ, Float32)
    end
end

@testitem "Operators: default ℓₘᵢₙ" begin
    import SphericalFunctions: L², Lz, L₊, L₋, Lx, Ly, R², Rz, R₊, R₋, ð, ð̄
    for O ∈ (L², Lz, L₊, L₋, Lx, Ly, R², Rz, R₊, R₋, ð, ð̄)
        for s ∈ -2:2, ℓₘₐₓ ∈ 2:5
            @test O(s, ℓₘₐₓ) == O(s, abs(s), ℓₘₐₓ, Float64)
            @test O(s, ℓₘₐₓ, Float32) == O(s, abs(s), ℓₘₐₓ, Float32)
        end
    end
end

@testitem "Operators: applied to ₛYₗₘ values" setup=[Utilities] begin
    # `ð` and `ð̄` are matrices acting on mode weights; this item checks that they really do
    # implement the *differential* operators of the same name acting on the corresponding
    # functions on the sphere.  The independent reference is Newman and Penrose's coordinate
    # form of those operators (conventions summary, "Spin-weighted functions"),
    #
    #     ð η = -sinˢθ {∂_θ + (i/sinθ) ∂_ϕ} (sin⁻ˢθ η),
    #     ð̄ η = -sin⁻ˢθ {∂_θ - (i/sinθ) ∂_ϕ} (sinˢθ η),
    #
    # applied by fourth-order central differences to the closed-form ₛYₗₘ of the `Utilities`
    # snippet.  Neither the differential operator nor the harmonic comes from the package, so
    # this is a genuinely independent check of the matrices' entries.
    #
    # Concretely: for each single mode (ℓ, m) of spin weight `s`, the mode weights `ð * Y` are
    # synthesized with the closed-form harmonics of spin weight s+1 and compared with ð
    # applied to the closed-form ₛYₗₘ; and likewise for ð̄ with spin weight s-1.  Errors are
    # accumulated over each (T, ℓₘₐₓ, s) block and asserted once, so a broken operator
    # produces a handful of failures rather than thousands.
    import SphericalFunctions: ð, ð̄
    using DoubleFloats

    # Fourth-order central difference.  The differencing is done in BigFloat with h = 1e-15,
    # near the optimum at the default 256-bit precision (truncation ~ h⁴ ≈ 1e-60, roundoff ~
    # eps/h ≈ 1e-62).  Checked against the ladder relations ð ₛYₗₘ = √((ℓ-s)(ℓ+s+1)) ₛ₊₁Yₗₘ and
    # ð̄ ₛYₗₘ = -√((ℓ+s)(ℓ-s+1)) ₛ₋₁Yₗₘ for every mode used below, the reference reproduces them
    # to 5.3e-58, which is what limits the BigFloat tolerance chosen at the bottom.
    h = big"1e-15"
    ∂(f, x) = (-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)) / (12h)
    function ðNP(s, f, θ, ϕ)
        g(t, p) = sin(t)^(-s) * f(t, p)
        -sin(θ)^s * (∂(t -> g(t, ϕ), θ) + im * ∂(p -> g(θ, p), ϕ) / sin(θ))
    end
    function ð̄NP(s, f, θ, ϕ)
        g(t, p) = sin(t)^s * f(t, p)
        -sin(θ)^(-s) * (∂(t -> g(t, ϕ), θ) - im * ∂(p -> g(θ, p), ϕ) / sin(θ))
    end

    ℓₘᵢₙ = 0
    ℓₘₐₓs = 4:7
    # The mode-weight ordering, written out rather than taken from the package's `Yindex`
    ℓmpairs(ℓₘₐₓ) = [(ℓ, m) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ]
    allpairs = ℓmpairs(maximum(ℓₘₐₓs))

    # Two generic points, well away from the poles where the sinᵗθ factors blow up
    θϕs = [(big"0.7", big"1.3"), (big"2.4", big"-0.8")]

    # Closed-form harmonics at those points, for every spin weight that can appear
    Yvals = Dict(
        (σ, p) => [
            ℓ < abs(σ) ? zero(Complex{BigFloat}) : sYlm(σ, ℓ, m, θϕs[p]...)
            for (ℓ, m) ∈ allpairs
        ]
        for σ ∈ -4:4, p ∈ eachindex(θϕs)
    )
    # ð and ð̄ of each single closed-form harmonic, by the differential operators above
    refs = Dict{NTuple{4,Int}, NTuple{2,Complex{BigFloat}}}()
    for s ∈ -3:3, (ℓ, m) ∈ allpairs, p ∈ eachindex(θϕs)
        ℓ < abs(s) && continue
        f(t, q) = sYlm(s, ℓ, m, t, q)
        refs[(s, ℓ, m, p)] = (ðNP(s, f, θϕs[p]...), ð̄NP(s, f, θϕs[p]...))
    end

    for T ∈ [Float32, Float64, Double64, BigFloat]
        # Measured maximum errors over everything below: 1.1e-7 (Float32), 2.2e-16 (Float64),
        # 1.6e-32 (Double64) — all within one eps of the respective type — and 5.3e-58 for
        # BigFloat, where the finite-difference reference rather than the arithmetic sets the
        # floor, so the tolerance cannot be 100eps(BigFloat) ≈ 1e-75 there.
        ϵ = max(100 * eps(T), 1e-55)
        @testset "$ℓₘₐₓ" for ℓₘₐₓ ∈ ℓₘₐₓs
            prs = ℓmpairs(ℓₘₐₓ)
            n = length(prs)
            @test prs == allpairs[1:n]  # the synthesis below relies on this
            for s ∈ -3:3
                𝔡, 𝔡̄ = ð(s, ℓₘᵢₙ, ℓₘₐₓ, T), ð̄(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                Y = zeros(Complex{T}, n)
                maxerr = 0.0
                subthreshold_zero = true
                for (i, (ℓ, m)) ∈ enumerate(prs)
                    ℓ < abs(s) && continue
                    Y .= zero(Complex{T})
                    Y[i] = one(Complex{T})
                    cð, cð̄ = 𝔡 * Y, 𝔡̄ * Y
                    # Modes below the spin weight of the *result* must be exactly zero
                    subthreshold_zero &= all(
                        iszero(cð[j]) for (j, (ℓⱼ, _)) ∈ enumerate(prs) if ℓⱼ < abs(s+1)
                    )
                    subthreshold_zero &= all(
                        iszero(cð̄[j]) for (j, (ℓⱼ, _)) ∈ enumerate(prs) if ℓⱼ < abs(s-1)
                    )
                    for p ∈ eachindex(θϕs)
                        ðY, ð̄Y = refs[(s, ℓ, m, p)]
                        maxerr = max(maxerr, Float64(abs(
                            sum(cð[j] * Yvals[(s+1, p)][j] for j ∈ 1:n) - ðY
                        )))
                        maxerr = max(maxerr, Float64(abs(
                            sum(cð̄[j] * Yvals[(s-1, p)][j] for j ∈ 1:n) - ð̄Y
                        )))
                    end
                end
                @test subthreshold_zero
                @test maxerr < ϵ
            end
        end
    end
end


### Half-integer indices
#
# The operator matrices for half-odd-integer `s`, `ℓₘᵢₙ` and `ℓₘₐₓ`, spelled as `Rational`s
# with denominator 2.  The mode ordering is `Yrange`, and the reference values are the
# textbook eigenvalues and ladder coefficients, formed from `Rational`s and `Int`s in the
# tests themselves rather than through the package's numerator arithmetic.

@testitem "Operators: half-integer eigenvectors and ladders" begin
    import SphericalFunctions: L², Lz, L₊, L₋, R², Rz, R₊, R₋, ð, ð̄
    import SphericalFunctions: Ysize, Yindex, Yrange, HalfOddInteger
    using LinearAlgebra: diag
    using DoubleFloats
    h(x) = HalfOddInteger(x)

    for T ∈ (Float32, Float64, Double64, BigFloat)
        ϵ = 10eps(T)
        for s ∈ (-5//2, -3//2, -1//2, 1//2, 3//2, 5//2), ℓₘₐₓ ∈ (7//2, 11//2), ℓₘᵢₙ ∈ unique((1//2, abs(s)))
            n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
            pairs = Yrange(ℓₘᵢₙ, ℓₘₐₓ)
            @test pairs == [(ℓ, m) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ]
            sh = h(s)
            l², lz, lp, lm = (O(s, ℓₘᵢₙ, ℓₘₐₓ, T) for O ∈ (L², Lz, L₊, L₋))
            r², rz, rp, rm, d, d̄ = (O(s, ℓₘᵢₙ, ℓₘₐₓ, T) for O ∈ (R², Rz, R₊, R₋, ð, ð̄))
            for M ∈ (l², lz, lp, lm, r², rz, rp, rm, d, d̄)
                @test size(M) == (n, n)
                @test eltype(M) === T
            end
            for (i, (ℓ, m)) ∈ enumerate(pairs)
                e = zeros(T, n)
                e[i] = 1
                if ℓ < abs(sh)
                    for M ∈ (l², lz, lp, lm, r², rz, rp, rm, d, d̄)
                        @test iszero(M * e)
                    end
                    continue
                end
                ℓr, mr = Rational(ℓ), Rational(m)
                # L² e = ℓ(ℓ+1) e, Lz e = m e, Rz e = s e, exactly
                @test l² * e == T(ℓr * (ℓr + 1)) .* e
                @test r² * e == T(ℓr * (ℓr + 1)) .* e
                @test lz * e == T(mr) .* e
                @test rz * e == T(s) .* e
                # L₊ maps (ℓ, m) to (ℓ, m+1) with coefficient √((ℓ-m)(ℓ+m+1)), and L₋ maps it
                # to (ℓ, m-1) with √((ℓ+m)(ℓ-m+1)); (ℓ∓m) and (ℓ±m+1) are `Int`s
                expected = zeros(T, n)
                if m < ℓ
                    expected[Yindex(ℓ, m + 1, h(ℓₘᵢₙ))] = √(T((ℓ - m) * (ℓ + m + 1)))
                end
                @test lp * e ≈ expected atol=ϵ rtol=ϵ
                expected = zeros(T, n)
                if m > -ℓ
                    expected[Yindex(ℓ, m - 1, h(ℓₘᵢₙ))] = √(T((ℓ + m) * (ℓ - m + 1)))
                end
                @test lm * e ≈ expected atol=ϵ rtol=ϵ
                # ð raises the spin weight with √((ℓ-s)(ℓ+s+1)), ð̄ lowers it with
                # -√((ℓ+s)(ℓ-s+1)), and R₊ = ð, R₋ = -ð̄
                @test d * e ≈ √(T((ℓ - sh) * (ℓ + sh + 1))) .* e atol=ϵ rtol=ϵ
                @test d̄ * e ≈ -√(T((ℓ + sh) * (ℓ - sh + 1))) .* e atol=ϵ rtol=ϵ
                @test rp * e == d * e
                @test rm * e == -(d̄ * e)
                # Entries below the spin weight of the *result* are exactly zero
                if ℓ < abs(sh + 1)
                    @test iszero(d * e)
                end
                if ℓ < abs(sh - 1)
                    @test iszero(d̄ * e)
                end
            end
        end
    end
end

@testitem "Operators: half-integer commutators and Casimir" begin
    import SphericalFunctions: L², Lz, L₊, L₋, Lx, Ly, R², Rz, R₊, R₋, ð, ð̄
    using DoubleFloats
    for T ∈ (Float32, Float64, Double64, BigFloat)
        ϵ = 100eps(T)
        for s ∈ (-3//2, -1//2, 1//2, 3//2), ℓₘₐₓ ∈ (7//2, 9//2), ℓₘᵢₙ ∈ unique((1//2, abs(s)))
            lz, lp, lm, l² = (Array(O(s, ℓₘᵢₙ, ℓₘₐₓ, T)) for O ∈ (Lz, L₊, L₋, L²))
            # [Lz, L±] = ±L±, [L₊, L₋] = 2Lz, and L² commutes with all of them
            @test lz*lp - lp*lz ≈ lp atol=ϵ rtol=ϵ
            @test lz*lm - lm*lz ≈ -lm atol=ϵ rtol=ϵ
            @test lp*lm - lm*lp ≈ 2lz atol=ϵ rtol=ϵ
            for O ∈ (lz, lp, lm)
                @test l²*O - O*l² ≈ zero(l²) atol=ϵ rtol=ϵ
            end
            # L² = (L₊L₋ + L₋L₊)/2 + Lz²
            @test (lp*lm + lm*lp)/2 + lz*lz ≈ l² atol=ϵ rtol=ϵ
            # The same for the right operators, with the spin-weight shift in each factor
            rz(σ) = Rz(σ, ℓₘᵢₙ, ℓₘₐₓ, T)
            rp(σ) = R₊(σ, ℓₘᵢₙ, ℓₘₐₓ, T)
            rm(σ) = R₋(σ, ℓₘᵢₙ, ℓₘₐₓ, T)
            r²(σ) = R²(σ, ℓₘᵢₙ, ℓₘₐₓ, T)
            @test rz(s+1)*rp(s) - rp(s)*rz(s) ≈ rp(s) atol=ϵ rtol=ϵ
            @test rz(s-1)*rm(s) - rm(s)*rz(s) ≈ -rm(s) atol=ϵ rtol=ϵ
            @test rp(s-1)*rm(s) - rm(s+1)*rp(s) ≈ 2rz(s) atol=ϵ rtol=ϵ
            @test (rp(s-1)*rm(s) + rm(s+1)*rp(s))/2 + rz(s)*rz(s) ≈ r²(s) atol=ϵ rtol=ϵ
            @test r²(s) ≈ l² atol=ϵ rtol=ϵ
            @test r²(s+1)*rp(s) - rp(s)*r²(s) ≈ zero(rp(s)) atol=ϵ rtol=ϵ
            # ð and ð̄: [Rz, ð] = ð, [Rz, ð̄] = -ð̄, [ð, ð̄] = -2Rz
            dd(σ) = ð(σ, ℓₘᵢₙ, ℓₘₐₓ, T)
            dd̄(σ) = ð̄(σ, ℓₘᵢₙ, ℓₘₐₓ, T)
            @test rz(s+1)*dd(s) - dd(s)*rz(s) ≈ dd(s) atol=ϵ rtol=ϵ
            @test rz(s-1)*dd̄(s) - dd̄(s)*rz(s) ≈ -dd̄(s) atol=ϵ rtol=ϵ
            @test dd(s-1)*dd̄(s) - dd̄(s+1)*dd(s) ≈ -2rz(s) atol=ϵ rtol=ϵ
            @test dd(s) == rp(s)
            @test dd̄(s) == -rm(s)
            # Lx and Ly: the defining combinations exactly, and the su(2) relations
            lx, ly = Matrix(Lx(s, ℓₘᵢₙ, ℓₘₐₓ, T)), Matrix(Ly(s, ℓₘᵢₙ, ℓₘₐₓ, T))
            @test eltype(lx) === T
            @test eltype(ly) === Complex{T}
            @test lx == (lp .+ lm) ./ 2
            @test 2im .* ly == lp .- lm
            @test lx == lx'
            @test ly == ly'
            @test lx^2 + ly^2 + lz^2 ≈ l² atol=ϵ rtol=ϵ
            @test lx*ly - ly*lx ≈ im*lz atol=ϵ rtol=ϵ
            @test ly*lz - lz*ly ≈ im*lx atol=ϵ rtol=ϵ
            @test lz*lx - lx*lz ≈ im*ly atol=ϵ rtol=ϵ
        end
    end
end

@testitem "Operators: half-integer spellings, mixed kinds and the integer path" begin
    import SphericalFunctions: L², Lz, L₊, L₋, Lx, Ly, R², Rz, R₊, R₋, ð, ð̄, Ysize, HalfOddInteger
    using LinearAlgebra: Diagonal, Bidiagonal, Tridiagonal, diag
    h(x) = HalfOddInteger(x)
    ops = (L², Lz, L₊, L₋, Lx, Ly, R², Rz, R₊, R₋, ð, ð̄)

    # Every spelling gives the same matrix, `ℓₘᵢₙ` defaults to `abs(s)`, and the matrix has the
    # same structure as in the integer case
    for O ∈ ops, s ∈ (-3//2, 1//2, 3//2), ℓₘₐₓ ∈ (5//2, 7//2), T ∈ (Float64, Float32)
        M = O(s, abs(s), ℓₘₐₓ, T)
        @test O(h(s), h(abs(s)), h(ℓₘₐₓ), T) == M
        @test O(s, h(abs(s)), ℓₘₐₓ, T) == M
        @test O(s, ℓₘₐₓ, T) == M
        @test O(h(s), h(ℓₘₐₓ), T) == M
        if T === Float64
            @test O(s, ℓₘₐₓ) == M
            @test O(s, abs(s), ℓₘₐₓ) == M
        end
        @test size(M) == (Ysize(abs(s), ℓₘₐₓ), Ysize(abs(s), ℓₘₐₓ))
        @test typeof(M) === typeof(O(1, 1, 3, T))
        # A mixture of the two kinds of index is refused with a message naming both spellings,
        # and a `Rational` that is not a half-odd-integer with one naming the denominator
        mixed = "all be integers, like 3, or all be half-odd-integers, like 7//2"
        @test_throws mixed O(s, 0, ℓₘₐₓ, T)
        @test_throws mixed O(1, abs(s), ℓₘₐₓ, T)
        @test_throws mixed O(s, abs(s), 3, T)
        @test_throws mixed O(s, 3, T)
        @test_throws mixed O(h(s), 0, h(ℓₘₐₓ), T)
        @test_throws "denominator 2" O(1//3, 1//3, 7//3, T)
        @test_throws "denominator 2" O(s, 3//1, T)
    end
    # Each public function carries both layers: a boundary method typed `IndexSpelling`, and a
    # worker typed `where {IT<:IntegerHalf}`.  The worker is reached by dispatch rather than by
    # a separate name — one `IT` for all three indices is strictly more specific than three
    # independent `IndexSpelling`s — so anything that is not already three indices of one kind
    # lands on the boundary, which normalizes it or explains why it cannot.
    for O ∈ ops
        worker = which(O, (HalfOddInteger, HalfOddInteger, HalfOddInteger, Type{Float64}))
        boundary = which(O, (Rational{Int}, Rational{Int}, Rational{Int}, Type{Float64}))
        @test worker !== boundary
        # Three indices of one kind reach the worker, whichever kind
        @test which(O, (Int, Int, Int, Type{Float64})) === worker
        # A `Rational`, a mixture of kinds, or a mixture of integer widths does not
        @test which(O, (Int, HalfOddInteger, HalfOddInteger, Type{Float64})) === boundary
        @test which(O, (Int8, Int, Int, Type{Float64})) === boundary
        # ... and every spelling infers the same concrete result
        RT = Base.infer_return_type(O, (Rational{Int}, Rational{Int}, Rational{Int}, Type{Float64}))
        @test isconcretetype(RT)
        @test RT === Base.infer_return_type(O, (HalfOddInteger, HalfOddInteger, HalfOddInteger, Type{Float64}))
        @test RT === Base.infer_return_type(O, (Int, Int, Int, Type{Float64}))
        @test RT === Base.infer_return_type(O, (Rational{Int}, Rational{Int}, Type{Float64}))
    end
    # The half-integer eigenvalues of L², in every type, against the `Rational` arithmetic
    for T ∈ (Float32, Float64, BigFloat)
        @test diag(L²(1//2, 1//2, 21//2, T)) == [T(ℓ * (ℓ + 1)) for ℓ ∈ 1//2:21//2 for m ∈ -ℓ:ℓ]
        @test diag(Lz(1//2, 1//2, 21//2, T)) == [T(m) for ℓ ∈ 1//2:21//2 for m ∈ -ℓ:ℓ]
    end

    # The integer path is unchanged: small cases computed by hand, exactly ...
    @test diag(L²(0, 0, 2)) == [0, 2, 2, 2, 6, 6, 6, 6, 6]
    @test diag(L²(1, 0, 2)) == [0, 2, 2, 2, 6, 6, 6, 6, 6]
    @test diag(L²(2, 0, 2)) == [0, 0, 0, 0, 6, 6, 6, 6, 6]
    @test diag(Lz(0, 0, 1)) == [0, -1, 0, 1]
    @test diag(Lz(0, 1, 2)) == [-1, 0, 1, -2, -1, 0, 1, 2]
    @test L₊(0, 0, 1).ev == [0, √2, √2]
    @test L₋(0, 0, 1).ev == [0, √2, √2]
    @test L₊(0, 1, 2).ev == [√2, √2, 0, 2, √6, √6, 2]
    @test L₋(0, 1, 2).ev == [√2, √2, 0, 2, √6, √6, 2]
    @test diag(Rz(-2, 2, 3)) == fill(-2.0, 12)
    @test diag(Rz(1, 0, 1)) == [0, 1, 1, 1]
    @test diag(ð(1, 1, 2)) == [0, 0, 0, 2, 2, 2, 2, 2]
    @test diag(ð̄(1, 1, 2)) == [-√2, -√2, -√2, -√6, -√6, -√6, -√6, -√6]
    @test diag(R₊(-1, 0, 1)) == [0, √2, √2, √2]
    @test diag(R₋(1, 0, 1)) == [0, √2, √2, √2]
    @test diag(R₊(0, 0, 1)) == [0, √2, √2, √2]
    @test diag(R₋(0, 0, 1)) == [0, √2, √2, √2]
    @test Lx(0, 0, 1) == Tridiagonal([0, √2, √2] ./ 2, zeros(4), [0, √2, √2] ./ 2)
    @test Ly(0, 0, 1) == Tridiagonal(-im .* [0, √2, √2] ./ 2, zeros(4), im .* [0, √2, √2] ./ 2)
    # ... with the same matrix types and element types as before ...
    @test L²(0, 0, 2, Float32) isa Diagonal{Float32, Vector{Float32}}
    @test L₊(0, 0, 2) isa Bidiagonal{Float64, Vector{Float64}}
    @test L₋(0, 0, 2) isa Bidiagonal{Float64, Vector{Float64}}
    @test Lx(0, 0, 2) isa Tridiagonal{Float64, Vector{Float64}}
    @test Ly(0, 0, 2) isa Tridiagonal{ComplexF64, Vector{ComplexF64}}
    @test ð(0, 0, 2, BigFloat) isa Diagonal{BigFloat, Vector{BigFloat}}
    # ... narrower or mixed integer types still agree with `Int` ...
    @test L²(Int8(1), Int8(1), Int8(3)) == L²(1, 1, 3)
    @test L₊(Int8(-1), 1, 3) == L₊(-1, 1, 3)
    @test ð(Int32(1), Int8(1), 3, Float32) == ð(1, 1, 3, Float32)
    # ... and the size, now `Ysize`, is the former (ℓₘₐₓ+1)² - ℓₘᵢₙ² in every case
    for ℓₘᵢₙ ∈ 0:3, ℓₘₐₓ ∈ ℓₘᵢₙ:6, O ∈ (L₊, L₋, Lx, Ly)
        @test size(O(0, ℓₘᵢₙ, ℓₘₐₓ)) == ((ℓₘₐₓ+1)^2 - ℓₘᵢₙ^2, (ℓₘₐₓ+1)^2 - ℓₘᵢₙ^2)
    end
end

@testitem "DifferentialOperator: the operators are values" begin
    import SphericalFunctions: DifferentialOperator, Δspin, bandstructure, coefftype,
        DiagonalBand, SubdiagonalBand, SuperdiagonalBand, TridiagonalBand

    ops = (L², Lz, L₊, L₋, Lx, Ly, R², Rz, R₊, R₋, ð, ð̄)

    for op ∈ ops
        @test op isa DifferentialOperator
        # Zero-size singletons, so dispatching on one costs nothing and every trait folds away
        @test Base.issingletontype(typeof(op))
        @test sizeof(op) == 0
        # They say their own names, so `repr` is `ð` rather than `SpinRaising()`
        @test repr(op) == string(nameof(op))
        @test getfield(SphericalFunctions, nameof(op)) === op
    end
    # No two share a type, so no two can share a trait by accident
    @test length(unique(typeof.(ops))) == length(ops)

    # The spin weight each one moves
    @test map(Δspin, ops) == (0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 1, -1)
    # The band each occupies
    @test bandstructure(L²) isa DiagonalBand
    @test bandstructure(L₊) isa SubdiagonalBand
    @test bandstructure(L₋) isa SuperdiagonalBand
    @test bandstructure(Lx) isa TridiagonalBand && bandstructure(Ly) isa TridiagonalBand
    # Ly is the only one whose entries are complex
    @test coefftype(Ly, Float64) === ComplexF64
    @test all(coefftype(op, Float64) === Float64 for op ∈ ops if op !== Ly)
end

@testitem "DifferentialOperator: op * w matches the matrix, bit for bit" begin
    import SphericalFunctions: Δspin
    using DoubleFloats: Double64
    using Random

    rng = Random.Xoshiro(20260920)
    ops = (L², Lz, L₊, L₋, Lx, Ly, R², Rz, R₊, R₋, ð, ð̄)

    # The sweep deliberately includes the small and degenerate containers the older items
    # skip: a single ℓ block, and the one-mode ℓ = 0 case.
    for T ∈ (Float32, Float64, Double64), CT ∈ (T, Complex{T})
        for (s, ℓₘᵢₙ, ℓₘₐₓ) ∈ ((0,0,3), (-2,2,4), (1,0,3), (0,1,1), (0,0,0), (2,2,2))
            data = randn(rng, CT, Ysize(ℓₘᵢₙ, ℓₘₐₓ))
            w = ModeWeights(copy(data), s, ℓₘᵢₙ, ℓₘₐₓ)
            for op ∈ ops
                M = op(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                # `==`, not `≈`: the loop and the comprehension evaluate the same coefficient
                # functions, and the loop sums a tridiagonal row left to right exactly as
                # `LinearAlgebra` does, so the two agree to the last bit.
                @test parent(op * w) == M * data
                @test op * w == op(w)                      # the two spellings agree
                @test spin(op * w) == s + Δspin(op)        # ... and the label moves correctly
                @test SphericalFunctions.ℓₘᵢₙ(op * w) === ℓₘᵢₙ   # the ℓ range never does
                @test SphericalFunctions.ℓₘₐₓ(op * w) === ℓₘₐₓ
                @test eltype(op * w) === eltype(M * data)
                @test parent(w) == data                    # the input is untouched
            end
        end
    end
end

@testitem "DifferentialOperator: op * w allocates once, mul! allocates nothing" begin
    import SphericalFunctions: Δspin
    using LinearAlgebra: mul!
    using Random

    rng = Random.Xoshiro(4242)
    ℓₘₐₓ = 12
    w = ModeWeights(randn(rng, ComplexF64, Ysize(0, ℓₘₐₓ)), 0, 0, ℓₘₐₓ)

    # Measured inside a function, never at top level, where the result is meaningless
    apply(op, w) = op * w
    inplace(dst, op, w) = mul!(dst, op, w)

    for op ∈ (L², Lz, L₊, L₋, Lx, Ly, R₊, ð)
        dst = ModeWeights(similar(parent(w)), 0 + Δspin(op), 0, ℓₘₐₓ)
        apply(op, w); inplace(dst, op, w)                  # warm up
        @test (@allocated inplace(dst, op, w)) == 0
        # `op * w` allocates its result and nothing else — never the operator matrix, which
        # for the banded ones is several times larger
        @test (@allocated apply(op, w)) < (@allocated op(0, 0, ℓₘₐₓ) * parent(w))
        @test parent(dst) == parent(op * w)
    end
end

@testitem "DifferentialOperator: mul! refuses a bad destination" begin
    import SphericalFunctions: Δspin
    using LinearAlgebra: mul!
    using Random

    rng = Random.Xoshiro(99)
    w = ModeWeights(randn(rng, ComplexF64, Ysize(0, 3)), 0, 0, 3)

    # The destination must be labelled with what the operator actually produces
    @test_throws "gives s=1" mul!(similar(w), ð, w)
    @test_throws "ℓ ∈ 0:3" mul!(ModeWeights(zeros(ComplexF64, Ysize(0, 4)), 0, 0, 4), Lz, w)
    # ... and it may not be the input: the banded kernels read a neighbour
    @test_throws "aliases the input" mul!(w, Lz, w)
    # A correctly labelled, separate destination works
    dst = ModeWeights(similar(parent(w)), 1, 0, 3)
    @test parent(mul!(dst, ð, w)) == parent(ð * w)
end
