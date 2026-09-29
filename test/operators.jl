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
# differentiation, to functions built from the package's own values.  Applied to 𝔇 and to
# ₛYₗₘ, in the first two items, they pin the sign conventions of the Wigner and sYlm
# functions against the definitions above.  Applied to the function that a `ModeWeights`
# synthesizes, they pin the package's operators themselves: `op * w`, evaluated at a rotor,
# must be the derivative of `w` evaluated there.  The remaining items that use them check
# the helper and the automatic differentiation through 𝔇, which the others rely on, rather
# than any convention.

@testitem "Pretest ε and basis commutators" setup=[Utilities] begin
    import .Utilities: ε
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
    # both positive; and it checks that the definition reproduces `sYlm_closed_form` from
    # the `Utilities` module — the explicit sum over factorials given on the conventions
    # pages, which shares no code with the package.  The operators are applied to the whole
    # vector of ₛYₗₘ values at once (ForwardDiff differentiates vector-valued functions),
    # which keeps the runtime reasonable.
    #
    # The closed form is a function of the spherical coordinates (θ, ϕ) alone, so the
    # natural comparison points are the rotors `from_spherical_coordinates(θ, ϕ)`.  Any
    # other rotor has an extra phase: writing 𝐐 in terms of its Euler angles (α, β, γ) as
    # 𝐐 = from_spherical_coordinates(β, α) * exp(γ𝐤/2), the defining property of spin
    # weight, η(𝐐 exp(γ𝐤/2)) = exp(-isγ) η(𝐐) (conventions summary, "Spin-weighted
    # functions"), supplies it.  Both kinds of point are used — with that phase where it is
    # needed — so that the operator identities are still exercised at generic rotors.
    import SphericalFunctions: D
    import .Utilities: sYlm_closed_form
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
    # The closed-form ₛYₗₘ, evaluated at the rotor 𝐐 by way of its Euler angles, with the
    # spin-weight phase discussed above.
    function closed_form(s, ℓ, m, Q)
        α, β, γ = to_euler_angles(Q)
        sYlm_closed_form(s, ℓ, m, β, α) * cis(-s * γ)
    end
    # γ = 0 for the first four, so the closed form applies to them with no phase at all.
    # The fourth is exactly at the north pole, θ = 0, and the fifth exactly at the south
    # pole.  There the recurrence's split of the rotor into the half angles of β and the
    # phases of α ± γ is singular, and derivatives taken through it would be NaN, so the
    # calculators give the derivatives from their values instead, by the rules for automatic
    # differentiation (see `src/derivatives.jl`); these two check the ForwardDiff-based
    # operators through those rules.
    Qs = [
        [from_spherical_coordinates(T(θ), T(ϕ)) for (θ, ϕ) ∈ ((0.4, 0.9), (1.0, 2.0), (2.5, -1.5), (0.0, 0.9))];
        [Rotor{T}(zero(T), T(0.6), T(0.8), zero(T))];
        randn(rng, Rotor{T}, 3)
    ]
    for Q ∈ Qs
        for s ∈ -2:2
            Y = Ys(s)(Q)
            # The definition agrees with the independent closed form.  The maximum error
            # over these points and spins measures 2.6e-15, well inside ϵ = 200eps(Float64)
            # ≈ 4.4e-14.  (Accumulated into one assertion rather than one per mode.)
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

@testitem "Pretest: explicit operators nest as written" setup=[ExplicitOperators] begin
    # A check of the `ExplicitOperators` helper, on which the items that differentiate
    # functions twice rely: `L(m, L(n, f))` is L_m L_n f, so the derivative with respect to
    # the generator `n` is taken innermost,
    #   LₘLₙf(Q) = (-i/2)² ∂ᵧ∂ᵨf(exp(ρn) exp(γm) Q),
    #   RₘRₙf(Q) = (i/2)² ∂ᵧ∂ᵨf(Q exp(γm) exp(ρn)).
    # The function differentiated is a matrix element of 𝔇, but any smooth function would
    # do, so this tests the helper and not the package; one precision and a few elements
    # suffice.
    import SphericalFunctions: D
    using Quaternionic
    import ForwardDiff
    using Random
    rng = Random.Xoshiro(123)

    const L = ExplicitOperators.L
    const R = ExplicitOperators.R

    T = Float64
    z = zero(T)
    LL(m, n, f, Q) = -ForwardDiff.derivative(
        γ -> ForwardDiff.derivative(ρ -> f(exp(ρ*n) * exp(γ*m) * Q), z), z
    ) / 4
    RR(m, n, f, Q) = -ForwardDiff.derivative(
        γ -> ForwardDiff.derivative(ρ -> f(Q * exp(γ*m) * exp(ρ*n)), z), z
    ) / 4

    ϵ = 100 * eps(T)
    M = randn(rng, QuatVec{T}, 2)
    N = randn(rng, QuatVec{T}, 2)
    Q = randn(rng, Rotor{T})
    for (m′, m) ∈ ((0, 0), (1, -1), (-1, 0)), n ∈ N, mm ∈ M
        f(Q) = D(Q, 1)[1][m′, m]
        @test L(mm, L(n, f))(Q) ≈ LL(mm, n, f, Q) atol=ϵ rtol=ϵ
        @test R(mm, R(n, f))(Q) ≈ RR(mm, n, f, Q) atol=ϵ rtol=ϵ
        # ... and the two orders differ, so the check can tell them apart
        @test !isapprox(L(mm, L(n, f))(Q), LL(n, mm, f, Q); atol=ϵ, rtol=ϵ)
    end
end

@testitem "Operators: basis commutators" setup=[ExplicitOperators] begin
    # [L_𝐮, L_𝐯] = (i/2) L_{[𝐮,𝐯]},   [R_𝐮, R_𝐯] = (i/2) R_{[𝐮,𝐯]},   [L_𝐮, R_𝐯] = 0
    #
    # These identities hold for any smooth function of the rotor, so they pin no convention.
    # What they check is that second derivatives of 𝔇 taken by nested automatic
    # differentiation are consistent with its first derivatives, which the explicit
    # operators of the other items rely on.
    import SphericalFunctions: D
    using Quaternionic
    using DoubleFloats
    using Random
    rng = Random.Xoshiro(1234)

    const L = ExplicitOperators.L
    const R = ExplicitOperators.R

    for T ∈ [Float64, Double64]
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

@testitem "Operators: op * w is the derivative of the function w" setup=[ExplicitOperators] begin
    # The mode weights `w` define a function on the rotors, `Q -> w(Q)`, and applying an
    # operator to the weights must give the weights of that operator applied to the
    # function.  So `(op * w)(Q)` is compared with the explicit Lie derivatives of `Q ->
    # w(Q)` at `Q`, by automatic differentiation through the package's own evaluation, for
    # every operator: the left and right components, the ladder operators L± = Lx ± i Ly and
    # ð = R₊ = Rx + i Ry, ð̄ = -R₋ = -(Rx - i Ry), and the Casimirs L² = R² = Σₐ Lₐ Lₐ = Σₐ
    # Rₐ Rₐ.  Nothing is compared with a matrix, so this pins the operators' conventions and
    # signs, for integer and half-integer spin weights alike.
    using Quaternionic
    using Random
    rng = Random.Xoshiro(2718)
    const L = ExplicitOperators.L
    const R = ExplicitOperators.R

    # The conventions do not depend on the precision, so one type suffices.
    let T = Float64
        x, y, z = QuatVec{T}(imx), QuatVec{T}(imy), QuatVec{T}(imz)
        for (s, ℓₘᵢₙ, ℓₘₐₓ) ∈ ((0, 0, 4), (1, 1, 4), (-2, 2, 5), (1//2, 1//2, 7//2), (-3//2, 3//2, 9//2))
            w = ModeWeights(randn(rng, Complex{T}, Ysize(ℓₘᵢₙ, ℓₘₐₓ)), s, ℓₘᵢₙ, ℓₘₐₓ)
            f = Q -> w(Q)
            # The worst error measured over these rotors, relative to eps(T), is 190, at
            # ℓₘₐₓ = 5, most of it in the Casimirs, whose values are ℓ(ℓ+1) times larger
            # than the function's; the bound grows with ℓₘₐₓ² accordingly.
            ϵ = 50ℓₘₐₓ^2 * eps(T)
            worst = zero(T)
            # Generic rotors, which stay away from the poles, where the half-angle square
            # roots inside the evaluation have infinite derivatives
            for Q ∈ randn(rng, Rotor{T}, 4)
                Lx_, Ly_, Lz_ = L(x, f)(Q), L(y, f)(Q), L(z, f)(Q)
                Rx_, Ry_, Rz_ = R(x, f)(Q), R(y, f)(Q), R(z, f)(Q)
                L²_ = L(x, L(x, f))(Q) + L(y, L(y, f))(Q) + L(z, L(z, f))(Q)
                R²_ = R(x, R(x, f))(Q) + R(y, R(y, f))(Q) + R(z, R(z, f))(Q)
                expected = (
                    (Lz, Lz_), (Lx, Lx_), (Ly, Ly_), (L₊, Lx_ + im*Ly_), (L₋, Lx_ - im*Ly_),
                    (Rz, Rz_), (R₊, Rx_ + im*Ry_), (R₋, Rx_ - im*Ry_),
                    (ð, Rx_ + im*Ry_), (ð̄, -(Rx_ - im*Ry_)), (L², L²_), (R², R²_),
                )
                for (op, value) ∈ expected
                    worst = max(worst, abs((op * w)(Q) - value))
                end
            end
            @test worst < ϵ
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
    # applied by fourth-order central differences to the closed-form ₛYₗₘ,
    # `sYlm_closed_form` from the `Utilities` module.  Neither the differential operator nor
    # the harmonic comes from the package, so this is a truly independent check of the
    # matrices' entries.
    #
    # Concretely: for each single mode (ℓ, m) of spin weight `s`, the mode weights `ð * Y`
    # are synthesized with the closed-form harmonics of spin weight s+1 and compared with ð
    # applied to the closed-form ₛYₗₘ; and likewise for ð̄ with spin weight s-1.  Errors are
    # accumulated over each (T, ℓₘₐₓ, s) block and asserted once, so a broken operator
    # produces a handful of failures rather than thousands.
    import SphericalFunctions: ð, ð̄
    import .Utilities: sYlm_closed_form
    using DoubleFloats

    # Fourth-order central difference.  The differencing is done in BigFloat with h = 1e-15,
    # near the optimum at the default 256-bit precision (truncation ~ h⁴ ≈ 1e-60, roundoff ~
    # eps/h ≈ 1e-62).  Checked against the ladder relations ð ₛYₗₘ = √((ℓ-s)(ℓ+s+1)) ₛ₊₁Yₗₘ
    # and ð̄ ₛYₗₘ = -√((ℓ+s)(ℓ-s+1)) ₛ₋₁Yₗₘ for every mode used below, the reference
    # reproduces them to 5.3e-58, which is what limits the BigFloat tolerance chosen at the
    # bottom.
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
            ℓ < abs(σ) ? zero(Complex{BigFloat}) : sYlm_closed_form(σ, ℓ, m, θϕs[p]...)
            for (ℓ, m) ∈ allpairs
        ]
        for σ ∈ -4:4, p ∈ eachindex(θϕs)
    )
    # ð and ð̄ of each single closed-form harmonic, by the differential operators above
    refs = Dict{NTuple{4,Int}, NTuple{2,Complex{BigFloat}}}()
    for s ∈ -3:3, (ℓ, m) ∈ allpairs, p ∈ eachindex(θϕs)
        ℓ < abs(s) && continue
        f(t, q) = sYlm_closed_form(s, ℓ, m, t, q)
        refs[(s, ℓ, m, p)] = (ðNP(s, f, θϕs[p]...), ð̄NP(s, f, θϕs[p]...))
    end

    for T ∈ [Float32, Float64, Double64, BigFloat]
        # Measured maximum errors over everything below: 1.1e-7 (Float32), 2.2e-16
        # (Float64), 1.6e-32 (Double64) — all within one eps of the respective type — and
        # 5.3e-58 for BigFloat, where the finite-difference reference rather than the
        # arithmetic sets the floor, so the tolerance cannot be 100eps(BigFloat) ≈ 1e-75
        # there.
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
                # to (ℓ, m-1) with √((ℓ+m)(ℓ-m+1)), formed from the `Rational`s
                expected = zeros(T, n)
                if m < ℓ
                    expected[Yindex(ℓ, m + 1, h(ℓₘᵢₙ))] = √(T((ℓr - mr) * (ℓr + mr + 1)))
                end
                @test lp * e ≈ expected atol=ϵ rtol=ϵ
                expected = zeros(T, n)
                if m > -ℓ
                    expected[Yindex(ℓ, m - 1, h(ℓₘᵢₙ))] = √(T((ℓr + mr) * (ℓr - mr + 1)))
                end
                @test lm * e ≈ expected atol=ϵ rtol=ϵ
                # ð raises the spin weight with √((ℓ-s)(ℓ+s+1)), ð̄ lowers it with
                # -√((ℓ+s)(ℓ-s+1)), and R₊ = ð, R₋ = -ð̄
                @test d * e ≈ √(T((ℓr - s) * (ℓr + s + 1))) .* e atol=ϵ rtol=ϵ
                @test d̄ * e ≈ -√(T((ℓr + s) * (ℓr - s + 1))) .* e atol=ϵ rtol=ϵ
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


@testitem "Operators: half-integer spellings, mixed kinds and the integer path" setup=[InferenceChecks] begin
    import SphericalFunctions: L², Lz, L₊, L₋, Lx, Ly, R², Rz, R₊, R₋, ð, ð̄, Ysize, HalfOddInteger
    import SphericalFunctions: IndexType
    using LinearAlgebra: Diagonal, Bidiagonal, Tridiagonal, diag
    using .InferenceChecks: inferred_type
    h(x) = HalfOddInteger(x)
    ops = (L², Lz, L₊, L₋, Lx, Ly, R², Rz, R₊, R₋, ð, ð̄)

    # Every spelling gives the same matrix, `ℓₘᵢₙ` defaults to `abs(s)`, and the matrix has
    # the same structure as in the integer case
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
        # A mixture of the two kinds of index is refused with a message that describes both
        # kinds and says which index is of which, and a `Rational` that is not a
        # half-odd-integer with one that says what it is; the message names the operator as
        # it prints
        mixed = "must all be integers of type `Int`, like 3, or all be half-odd-integers"
        @test_throws ArgumentError O(s, 0, ℓₘₐₓ, T)
        @test_throws mixed O(s, 0, ℓₘₐₓ, T)
        @test_throws "mixes integers (ℓₘᵢₙ) with half-odd-integers (s, ℓₘₐₓ)" O(s, 0, ℓₘₐₓ, T)
        @test_throws mixed O(1, abs(s), ℓₘₐₓ, T)
        @test_throws mixed O(s, abs(s), 3, T)
        @test_throws mixed O(s, 3, T)
        @test_throws mixed O(h(s), 0, h(ℓₘₐₓ), T)
        @test_throws "The indices of one call to `$(nameof(O))`" O(s, 3, T)
        @test_throws ArgumentError O(1//3, 1//3, 7//3, T)
        @test_throws "1//3 is neither an integer nor a half-odd-integer" O(1//3, 1//3, 7//3, T)
        @test_throws "3//1 is a whole number; write it as the integer 3" O(s, 3//1, T)
    end
    # Each call shape of the operators is one `@index_methods` definition.  Three indices of
    # one kind reach the work method of that kind, `Rational`s reach the method that
    # converts them, and anything else (a mixture of kinds, or an integer of a type other
    # than `Int`) reaches the fallback, which explains the refusal.
    for O ∈ ops
        fallback = which(O, (IndexType, IndexType, IndexType, Type{Float64}))
        integer_work = which(O, (Int, Int, Int, Type{Float64}))
        half_work = which(O, (HalfOddInteger, HalfOddInteger, HalfOddInteger, Type{Float64}))
        conversion = which(O, (Rational{Int}, Rational{Int}, Rational{Int}, Type{Float64}))
        @test allunique((fallback, integer_work, half_work, conversion))
        @test which(O, (Int, HalfOddInteger, HalfOddInteger, Type{Float64})) === fallback
        @test which(O, (Int8, Int, Int, Type{Float64})) === fallback
        @test which(O, (UInt, Int, Int, Type{Float64})) === fallback
        @test which(O, (Bool, Int, Int, Type{Float64})) === fallback
        # ... and every spelling infers the same concrete result
        RT = inferred_type(O, (Rational{Int}, Rational{Int}, Rational{Int}, Type{Float64}))
        @test isconcretetype(RT)
        @test RT === inferred_type(O, (HalfOddInteger, HalfOddInteger, HalfOddInteger, Type{Float64}))
        @test RT === inferred_type(O, (Int, Int, Int, Type{Float64}))
        @test RT === inferred_type(O, (Rational{Int}, Rational{Int}, Type{Float64}))
    end
    # The half-integer eigenvalues of L², in every type, against the `Rational` arithmetic
    for T ∈ (Float32, Float64, BigFloat)
        @test diag(L²(1//2, 1//2, 21//2, T)) == [T(ℓ * (ℓ + 1)) for ℓ ∈ 1//2:21//2 for m ∈ -ℓ:ℓ]
        @test diag(Lz(1//2, 1//2, 21//2, T)) == [T(m) for ℓ ∈ 1//2:21//2 for m ∈ -ℓ:ℓ]
    end

    # The integer path: small cases computed by hand, exactly ...
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
    @test Ly(0, 0, 1) == Tridiagonal(-im .* [0, √2, √2] ./ 2, zeros(ComplexF64, 4), im .* [0, √2, √2] ./ 2)
    # ... with these matrix types and element types ...
    @test L²(0, 0, 2, Float32) isa Diagonal{Float32, Vector{Float32}}
    @test L₊(0, 0, 2) isa Bidiagonal{Float64, Vector{Float64}}
    @test L₋(0, 0, 2) isa Bidiagonal{Float64, Vector{Float64}}
    @test Lx(0, 0, 2) isa Tridiagonal{Float64, Vector{Float64}}
    @test Ly(0, 0, 2) isa Tridiagonal{ComplexF64, Vector{ComplexF64}}
    @test ð(0, 0, 2, BigFloat) isa Diagonal{BigFloat, Vector{BigFloat}}
    # ... while an index of an integer type other than `Int`, or a `Rational` of another
    # integer type, is refused with a message that says why and how to write it instead ...
    @test_throws ArgumentError L²(Int8(1), Int8(1), Int8(3))
    @test_throws "`Int8` is narrower than `Int`" L²(Int8(1), Int8(1), Int8(3))
    @test_throws "`Int8` is narrower than `Int`" L₊(Int8(-1), 1, 3)
    @test_throws "`Int32` is narrower than `Int`" ð(Int32(1), Int8(1), 3, Float32)
    @test_throws ArgumentError L₊(0, 0, UInt(3))
    @test_throws "`UInt64` is unsigned" L₊(0, 0, UInt(3))
    @test_throws "A `Bool` is not an index" Lz(true, 3)
    @test_throws "`BigInt` is wider than `Int`" ð̄(0, big(3))
    @test_throws "`Rational{Int8}` is not `Rational{Int}`" Lz(Int8(1)//Int8(2), 5//2)
    # ... and the size is `Ysize`, which is (ℓₘₐₓ+1)² - ℓₘᵢₙ² for integer indices
    for ℓₘᵢₙ ∈ 0:3, ℓₘₐₓ ∈ ℓₘᵢₙ:6, O ∈ (L₊, L₋, Lx, Ly)
        @test size(O(0, ℓₘᵢₙ, ℓₘₐₓ)) == ((ℓₘₐₓ+1)^2 - ℓₘᵢₙ^2, (ℓₘₐₓ+1)^2 - ℓₘᵢₙ^2)
    end
end

@testitem "DifferentialOperator: op * w matches the matrix, bit for bit" begin
    import SphericalFunctions: Δspin
    using DoubleFloats: Double64
    using Random

    rng = Random.Xoshiro(20260920)
    ops = (L², Lz, L₊, L₋, Lx, Ly, R², Rz, R₊, R₋, ð, ð̄)

    # The sweep deliberately includes the small and degenerate containers that the other
    # items skip, a single ℓ block and the one-mode ℓ = 0 case, and half-integer containers.
    for T ∈ (Float32, Float64, Double64), CT ∈ (T, Complex{T})
        for (s, ℓₘᵢₙ, ℓₘₐₓ) ∈ (
            (0,0,3), (-2,2,4), (1,0,3), (0,1,1), (0,0,0), (2,2,2),
            (1//2,1//2,7//2), (-3//2,3//2,5//2), (1//2,1//2,1//2),
        )
            data = randn(rng, CT, Ysize(ℓₘᵢₙ, ℓₘₐₓ))
            w = ModeWeights(copy(data), s, ℓₘᵢₙ, ℓₘₐₓ)
            for op ∈ ops
                M = op(s, ℓₘᵢₙ, ℓₘₐₓ, T)
                # `==`, not `≈`: the loop and the comprehension evaluate the same
                # coefficient functions, and the loop sums a tridiagonal row left to right
                # exactly as `LinearAlgebra` does, so for finite data the two agree to the
                # last bit (apart from the sign of a zero).  The exception is a matrix of
                # three rows or fewer, which `LinearAlgebra` before Julia 1.12 multiplies as
                # a dense one, through BLAS, whose rounding can differ in the last bit.
                if VERSION < v"1.12" && size(M, 1) ≤ 3
                    @test parent(op * w) ≈ M * data rtol=4eps(T)
                else
                    @test parent(op * w) == M * data
                end
                @test op * w == op(w)                      # the two spellings agree
                @test spin(op * w) == s + Δspin(op)        # ... and the label moves correctly
                @test SphericalFunctions.ℓₘᵢₙ(op * w) == ℓₘᵢₙ   # the ℓ range never does
                @test SphericalFunctions.ℓₘₐₓ(op * w) == ℓₘₐₓ
                @test eltype(op * w) === eltype(M * data)
                @test parent(w) == data                    # the input is untouched
            end
        end
    end

    # Where the data are not finite the two differ, as the docstring of
    # `DifferentialOperator` says: the matrix multiplies its stored zero diagonal by an
    # infinite weight, giving a NaN, whereas `op * w` reads only the band that the operator
    # occupies.
    let data = zeros(9)
        data[3] = Inf
        w = ModeWeights(data, 0, 0, 2)
        @test iszero(parent(L₊ * w)[3]) && parent(L₊ * w)[4] == Inf
        @test isnan((L₊(0, 0, 2) * data)[3])
    end
end

@testitem "DifferentialOperator: products and broadcasts" begin
    using Random

    rng = Random.Xoshiro(314)
    w = ModeWeights(randn(rng, ComplexF64, Ysize(0, 4)), 0, 0, 4)

    # `a * b * w` is `a * (b * w)`, for any number of operators ...
    @test ð̄ * ð * w == ð̄ * (ð * w)
    @test spin(ð̄ * ð * w) == 0
    @test L₊ * L₋ * Lz * w == L₊ * (L₋ * (Lz * w))
    @test ð * ð * ð̄ * w == ð * (ð * (ð̄ * w))
    @test spin(ð * ð * ð̄ * w) == 1
    # ... so that identities among the operators can be checked as they are written: the
    # Casimir L² = (L₊L₋ + L₋L₊)/2 + Lz², and ð̄ð = -(L² - Rz² - Rz), which follows from [ð,
    # ð̄] = -2Rz and L² = R²
    @test (L₊ * L₋ * w + L₋ * L₊ * w) / 2 + Lz * Lz * w ≈ L² * w
    @test ð̄ * ð * w ≈ -(L² * w - Rz * Rz * w - Rz * w)
    let wh = ModeWeights(randn(rng, ComplexF64, Ysize(1//2, 7//2)), 1//2, 1//2, 7//2)
        @test ð̄ * ð * wh ≈ -(L² * wh - Rz * Rz * wh - Rz * wh)
        @test spin(ð̄ * ð * wh) == 1//2
    end
    # A product of two operators on their own is not defined
    @test_throws MethodError ð̄ * ð

    # In a broadcast an operator is a scalar, so `op .* ws` applies it to each element, as
    # the call form `op.(ws)` does
    ws = [w, 2w, ModeWeights(randn(rng, ComplexF64, Ysize(1, 3)), 1, 1, 3)]
    @test ð .* ws == ð.(ws) == [ð * v for v ∈ ws]
    @test spin.(ð .* ws) == [1, 1, 2]
    @test Lz .* ws == [Lz * v for v ∈ ws]
end

@testitem "DifferentialOperator: op * w allocates once, mul! allocates nothing" begin
    import SphericalFunctions: Δspin
    using LinearAlgebra: mul!
    using Random

    rng = Random.Xoshiro(4242)
    ℓₘₐₓ = 12
    w = ModeWeights(randn(rng, ComplexF64, Ysize(0, ℓₘₐₓ)), 0, 0, ℓₘₐₓ)

    # Each measurement is made inside a function, in which the types of the arguments are
    # known; at top level it would include the dynamic dispatch of the call being measured,
    # which allocates on some versions of Julia.
    allocations_apply(op, w) = @allocated op * w
    allocations_inplace(dst, op, w) = @allocated mul!(dst, op, w)
    allocations_matrix(op, w, ℓₘₐₓ) = @allocated op(0, 0, ℓₘₐₓ) * parent(w)

    for op ∈ (L², Lz, L₊, L₋, Lx, Ly, R², Rz, R₊, R₋, ð, ð̄)
        dst = ModeWeights(similar(parent(w)), 0 + Δspin(op), 0, ℓₘₐₓ)
        # warm up
        allocations_apply(op, w); allocations_inplace(dst, op, w); allocations_matrix(op, w, ℓₘₐₓ)
        @test allocations_inplace(dst, op, w) == 0
        # `op * w` allocates its result and nothing else, never the operator matrix, which
        # the product with the matrix must allocate as well as its result
        @test allocations_apply(op, w) < allocations_matrix(op, w, ℓₘₐₓ)
        @test parent(dst) == parent(op * w)
    end
end

@testitem "DifferentialOperator: mul! refuses a bad destination" begin
    import SphericalFunctions: Δspin, spin
    using LinearAlgebra: mul!
    using Random

    rng = Random.Xoshiro(99)
    w = ModeWeights(randn(rng, ComplexF64, Ysize(0, 3)), 0, 0, 3)

    # The destination must be labelled with what the operator actually produces
    @test_throws "gives s=1" mul!(similar(w), ð, w)
    @test_throws "ℓ ∈ 0:3" mul!(ModeWeights(zeros(ComplexF64, Ysize(0, 4)), 0, 0, 4), Lz, w)
    # ... and it may not be the input: the banded kernels read a neighbor
    @test_throws "aliases the input" mul!(w, Lz, w)
    # A correctly labelled, separate destination works
    dst = ModeWeights(similar(parent(w)), 1, 0, 3)
    @test parent(mul!(dst, ð, w)) == parent(ð * w)
    # ... and so does a bare vector at least as long as the result, which comes back
    # labelled with what the operator produces (here s = 1), as a ModeWeights over its first
    # entries
    out = zeros(ComplexF64, Ysize(0, 3) + 1)
    w′ = mul!(out, ð, w)
    @test w′ isa ModeWeights && spin(w′) == 1 && parent(array_view(w′)) === out
    @test array_view(w′) == parent(ð * w)
    @test iszero(out[end])
    @test_throws "at least" mul!(zeros(ComplexF64, Ysize(0, 3) - 1), ð, w)
    @test_throws "aliases the input" mul!(parent(w), Lz, w)

    # Half-integer weights take the same route, with the destination labelled by `Δspin`
    wh = ModeWeights(randn(rng, ComplexF64, Ysize(1//2, 7//2)), 1//2, 1//2, 7//2)
    for op ∈ (ð, ð̄, Lx, L₋, Rz)
        dsth = ModeWeights(similar(parent(wh)), 1//2 + Δspin(op), 1//2, 7//2)
        @test parent(mul!(dsth, op, wh)) == parent(op * wh)
        @test spin(dsth) == spin(op * wh)
    end
    @test_throws "gives s=3//2" mul!(similar(wh), ð, wh)
end

@testitem "Operators: every band structure refuses an invalid range of ℓ" begin
    import SphericalFunctions: L², Lz, L₊, L₋, Lx, Ly, R², Rz, R₊, R₋, ð, ð̄

    # Each operator acts on the ordering that `Ysize` counts, so each refuses what `Ysize`
    # refuses — a negative ℓₘᵢₙ, or an ℓₘₐₓ below ℓₘᵢₙ-1 — whatever the shape of its matrix
    for op ∈ (L², Lz, L₊, L₋, Lx, Ly, R², Rz, R₊, R₋, ð, ð̄)
        @test_throws ArgumentError op(0, -2, 3)
        @test_throws "ℓₘᵢₙ=-2 must be non-negative" op(0, -2, 3)
        @test_throws "ℓₘᵢₙ=-1 must be non-negative" op(0, -1, 2)
        @test_throws ArgumentError op(0, 5, 2)
        @test_throws "ℓₘₐₓ=2 must be at least ℓₘᵢₙ-1=4" op(0, 5, 2)
        @test_throws "ℓₘₐₓ=1 must be at least ℓₘᵢₙ-1=2" op(0, 3, 1)
        @test_throws "ℓₘₐₓ=3 must be at least ℓₘᵢₙ-1=4" op(5, 3)  # ℓₘᵢₙ is |s| = 5 by default
        @test_throws ArgumentError op(-1//2, -3//2, 5//2)
        @test_throws "must be non-negative" op(-1//2, -3//2, 5//2)
        @test_throws "must be non-negative" op(1//2, -1//2, 5//2)
        @test_throws "must be at least" op(1//2, 9//2, 5//2)
        # The empty range ℓₘₐₓ = ℓₘᵢₙ-1 is legal, and gives an empty matrix
        @test size(op(5, 5, 4)) == (0, 0)
        @test size(op(1//2, 7//2, 5//2)) == (0, 0)
    end
end
