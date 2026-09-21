# Tests of the unbatched, single-matrix form of the `H` recurrence in `src/wigner/recurrence.jl`
# — `recurrence_step1!` … `recurrence_step6!`, `convert_H_to_d!` and `convert_H_to_D!`, the
# functions documented on `docs/src/40-api/01-internal.md`.
#
# Why this exists.  These functions take an `AbstractWignerMatrix` holding one whole `Hˡ` for
# one rotor; the engine that the package actually runs (`WignerHCalculator`) has its own,
# separate methods of the same names in `src/wigner/wigner_H_calculator.jl`, which work on a
# batched quarter-wedge.  Until 2026-09-11 the single-matrix path was driven by the internal
# `DenseWignerCalculator`, and the item "Wigner calculators vs DenseWignerCalculator"
# compared the two; that calculator was deleted along with `Deprecated`, leaving this path
# with no caller and no test.  This item restores the cross-check directly: it is a genuinely
# independent second implementation of the same recurrence (different loop structure,
# different storage, no batching), and it is what caught bug B1 of the v3 design memo.

@testitem "Dense H recurrence vs the batched engine" begin
    import SphericalFunctions as SF
    import SphericalFunctions: WignerMatrix, D, d
    import SphericalFunctions:
        recurrence_step1!, recurrence_step2!, recurrence_step3!,
        recurrence_step4!, recurrence_step5!, recurrence_step6!,
        convert_H_to_d!, convert_H_to_D!, spinor_phases
    using Quaternionic: Rotor, 𝐢, 𝐣, 𝐤
    using Random: Xoshiro

    # Drive the six documented steps by hand, exactly as the notes on the `H` recursion
    # describe them, and return the filled `Hˡ`.  `axes_[n+1]` holds the m′=0 axis of `Hⁿ`;
    # step 3 needs the axis one order *above* the target ℓ.
    function dense_H(::Type{NT}, ℓ, cosβ, sinβ; m′ₘₐₓ=ℓ) where {NT}
        axes_ = [WignerMatrix(zeros(NT, 2n + 1, 2n + 1), n) for n ∈ 0:ℓ+1]
        recurrence_step1!(axes_[1])                               # H⁰₀₀ = 1
        for n ∈ 1:ℓ+1
            recurrence_step2!(axes_[n+1], axes_[n], sinβ, cosβ)   # Hⁿ⁻¹₀ₘ -> Hⁿ₀ₘ
        end
        Hˡ = WignerMatrix(
            zeros(NT, 2m′ₘₐₓ + 1, 2ℓ + 1), ℓ; m′ₘₐₓ=m′ₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ
        )
        for m ∈ 0:ℓ
            Hˡ[0, m] = axes_[ℓ+1][0, m]
        end
        recurrence_step3!(Hˡ, axes_[ℓ+2], sinβ, cosβ)             # Hˡ⁺¹₀ₘ -> Hˡ₁ₘ
        recurrence_step4!(Hˡ, sinβ, cosβ)                         # ... -> Hˡₘ′₊₁ₘ
        recurrence_step5!(Hˡ, sinβ, cosβ)                         # ... -> Hˡₘ′₋₁ₘ
        recurrence_step6!(Hˡ)                                     # the symmetries
        Hˡ
    end

    # The poles (where the recurrence's special branches live), both signs of the double
    # cover, and a few seeded random rotors.
    rotors(::Type{T}) where {T} = [
        Rotor{T}(1); Rotor{T}(𝐢); Rotor{T}(𝐣); Rotor{T}(𝐤);
        -Rotor{T}(1); -Rotor{T}(𝐣);
        randn(Xoshiro(1729), Rotor{T}, 4)
    ]

    @testset "$T" for T ∈ (Float64, BigFloat)
        ℓₘₐₓ = 6
        # Measured worst case over everything below: 3.3e-16 (d) and 4.9e-16 (𝔇) in
        # Float64, i.e. about 2 eps; the two implementations differ only in rounding.
        atol = 20 * eps(T)
        for R ∈ rotors(T)
            eⁱᵝ, z₊, z₋, _, _ = spinor_phases(R)
            cosβ, sinβ = reim(eⁱᵝ)
            # e^{iα} = z₊ z₋ and e^{iγ} = z₊ conj(z₋), from z₊ = e^{i(α+γ)/2},
            # z₋ = e^{i(α-γ)/2}; taken this way, no Euler angle is ever extracted.
            eⁱᵅ, eⁱᵞ = z₊ * z₋, z₊ * conj(z₋)
            for ℓ ∈ 0:ℓₘₐₓ
                dref = d(eⁱᵝ, ℓ)[ℓ]
                𝔇ref = D(R, ℓ)[ℓ]
                # Errors are accumulated and asserted once per (rotor, ℓ): a per-element
                # `@test` would be tens of thousands of assertions.
                errd = zero(T)
                err𝔇 = zero(T)
                for m′ₘₐₓ ∈ 0:ℓ
                    Hd = dense_H(T, ℓ, cosβ, sinβ; m′ₘₐₓ)
                    convert_H_to_d!(Hd)
                    H𝔇 = dense_H(Complex{T}, ℓ, cosβ, sinβ; m′ₘₐₓ)
                    convert_H_to_D!(H𝔇, eⁱᵅ, eⁱᵞ)
                    for m′ ∈ -m′ₘₐₓ:m′ₘₐₓ, m ∈ -ℓ:ℓ
                        errd = max(errd, abs(Hd[m′, m] - dref[m′, m]))
                        err𝔇 = max(err𝔇, abs(H𝔇[m′, m] - 𝔇ref[m′, m]))
                    end
                end
                @test errd ≤ atol
                @test err𝔇 ≤ atol
            end
        end
    end

    # `H` itself, before the phases: the symmetries step 6 imposes must hold exactly, since
    # they are assignments from one stored element to another.
    let T = Float64, ℓ = 5
        eⁱᵝ, = spinor_phases(randn(Xoshiro(11), Rotor{T}))
        cosβ, sinβ = reim(eⁱᵝ)
        H = dense_H(T, ℓ, cosβ, sinβ)
        for m′ ∈ -ℓ:ℓ, m ∈ -ℓ:ℓ
            @test H[m′, m] == H[m, m′]
            @test H[m′, m] == H[-m′, -m]
        end
    end

    # Step 1 only initializes ℓ=0.
    @test_throws ErrorException recurrence_step1!(WignerMatrix(zeros(3, 3), 1))
end
