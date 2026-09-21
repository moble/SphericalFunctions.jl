# Tests of `WignerHCalculator`, the batched engine that runs the Wigner recurrences and
# produces the H wedge for several rotors at once.  The oracle for its values is a
# closed-form H, built from the definition of H in terms of the Wigner d function; see the
# first item below.

@testitem "WignerHCalculator vs closed-form H" setup=[HalfIntegerOracle] begin
    import SphericalFunctions: WignerHCalculator, recurrence!, wedge_value
    import .HalfIntegerOracle: d_oracle
    import DoubleFloats: Double64
    import Random

    # Reference, category 1 (a closed-form formula written out in the test).  H is not a
    # standard published quantity: it is *defined* by
    #
    #     d^ℓ_{m′,m}(β) = ϵ(m′) ϵ(-m) H^ℓ_{m′,m}(β),   ϵ(k) = (-1)^⌊k⌋ for k > 0, else 1
    #
    # (`docs/src/50-notes/01-H_recurrence.md`; design memo §5.2).  The ϵ are signs, so the
    # relation inverts to H = ϵ(m′) ϵ(-m) d, and `d` is taken from Varshalovich
    # Eq. 4.3.1(2) as transcribed in the `HalfIntegerOracle` setup module — a reference
    # that owes nothing to this package.
    #
    # Every index here is an integer, where the symmetry sign σ = sgn(m′) sgn(m) of
    # `H_recurrence.md` is identically +1, so H_{m′,m} = H_{m,m′} = H_{-m′,-m} and the one
    # expression above is the reference both for the elements stored in the wedge and for
    # the elements `wedge_value` reconstructs from them by symmetry.  (The half-integer
    # case, where σ is genuinely -1 whenever sgn(m′) ≠ sgn(m), reaches the same reads
    # through the dᴶ and 𝔇ᴶ comparisons of `test/wigner/half_integer.jl`.)
    ϵ(k) = ifelse(k > 0 && isodd(k), -1, 1)

    # The closed form is evaluated at four times the working precision and rounded to `T`,
    # so the reference has no more than half an ulp of its own error.  Evaluated *at*
    # the working precision it would be useless as an oracle: its alternating sum loses
    # about 270 eps at ℓ = 16 in Float64.
    function Href(::Type{T}, ℓ, m′, m, β) where {T<:Real}
        r = setprecision(BigFloat, 4 * precision(T) + 64) do
            ϵ(m′) * ϵ(-m) * d_oracle(ℓ, m′, m, BigFloat(β))
        end
        T(r)
    end

    ℓᵗᵒᵖ = 16  # the largest ℓ reached by any configuration below

    @testset "$T" for T in (Float64, Double64, BigFloat)
        rng = Random.Xoshiro(1234)
        # Both poles, their nearest neighbors, and four random angles — the same set as
        # `Utilities.βrange`, which cannot build its `T(0):eps:T(π)` range for Double64.
        β⃗ = T[0; nextfloat(T(0)); T.(rand(rng, 4)) .* T(π); prevfloat(T(π)); T(π)]
        @test length(β⃗) == 8
        @test iszero(first(β⃗)) && last(β⃗) == T(π)

        # H^ℓ_{m′,m}(β) depends on nothing but ℓ, m′, m and β, so one table of references
        # serves every (ℓₘₐₓ, m′ₘₐₓ, Nᵣ) configuration below: `ref[k][ℓ+1][m′+ℓ+1, m+ℓ+1]`
        # is H^ℓ_{m′,m}(β⃗[k]).
        ref = [
            [T[Href(T, ℓ, m′, m, β) for m′ in -ℓ:ℓ, m in -ℓ:ℓ] for ℓ in 0:ℓᵗᵒᵖ]
            for β in β⃗
        ]

        # Measured worst case over every configuration below, with no growth in ℓ:
        # 2.5 eps (Float64), 2.3 eps (Double64), 3.0 eps (BigFloat) at ℓ ≤ 16.
        atol = 10 * eps(T)

        @testset "ℓₘₐₓ=$ℓₘₐₓ" for ℓₘₐₓ in (0, 1, 2, 3, 7, 16)
            m′ₘₐₓs = ℓₘₐₓ ≤ 7 ? (0:ℓₘₐₓ) : (0, 1, 2, 5, 8, 16)
            @testset "m′ₘₐₓ=$m′ₘₐₓ" for m′ₘₐₓ in m′ₘₐₓs
                for Nᵣ in (1, 4)
                    calc = WignerHCalculator(β⃗[1:Nᵣ], ℓₘₐₓ; m′ₘₐₓ)
                    # The errors are accumulated and asserted once per configuration: a
                    # per-element `@test` here would be a million assertions, and a broken
                    # engine would spend minutes printing them.
                    stored = zero(T)   # worst error over the stored wedge
                    symmetric = zero(T)  # worst error over the reads via the symmetries
                    # Split the angles into batches of Nᵣ, so that both values of Nᵣ cover
                    # exactly the same set of angles.
                    for batch in Iterators.partition(eachindex(β⃗), Nᵣ)
                        for ℓ in 0:ℓₘₐₓ
                            if ℓ == 0
                                recurrence!(calc, β⃗[batch], ℓ)  # this batch's angles
                            else
                                recurrence!(calc, ℓ)  # the cheap sequential step
                            end
                            Hˡ = calc.Hˡ
                            @test Hˡ.ℓ == ℓ
                            m′ᵐᵃˣ = min(ℓ, m′ₘₐₓ)
                            for (iᵣ, k) in enumerate(batch)
                                Hᵏ = ref[k][ℓ+1]
                                # Every stored element of the wedge, read directly
                                for m′ in -m′ᵐᵃˣ:m′ᵐᵃˣ, m in abs(m′):ℓ
                                    stored = max(
                                        stored, abs(Hˡ[iᵣ, m′, m] - Hᵏ[m′+ℓ+1, m+ℓ+1])
                                    )
                                end
                                # Every element the wedge can supply, read via the symmetries
                                for m′ in -ℓ:ℓ, m in -ℓ:ℓ
                                    min(abs(m′), abs(m)) ≤ m′ₘₐₓ || continue
                                    symmetric = max(
                                        symmetric,
                                        abs(wedge_value(Hˡ, iᵣ, m′, m) - Hᵏ[m′+ℓ+1, m+ℓ+1])
                                    )
                                end
                            end
                        end
                    end
                    @test stored ≤ atol
                    @test symmetric ≤ atol
                end
            end
        end
    end
end


@testitem "WignerHCalculator ℓ ordering" begin
    import SphericalFunctions
    import SphericalFunctions: WignerHCalculator, HWedge, recurrence!
    import Random

    # Snapshot of the stored wedge for the current ℓ, in storage order
    function wedge(H::HWedge)
        [H[iᵣ, m′, m] for m′ in H.m′ₘᵢₙ:H.m′ₘₐₓ for m in abs(m′):H.ℓ for iᵣ in 1:H.Nᵣ]
    end

    T = Float64
    ℓₘₐₓ = 9
    rng = Random.Xoshiro(3)
    β⃗ = T[0; rand(rng, 2) .* T(π); T(π)]

    # Repeats, forward jumps, and backward moves, including returns to 0 and to ℓₘₐₓ
    order = (0, 3, 1, 9, 9, 5, 6, 0, 2, 8, 4, 9, 7, 7, 0, 9)

    for m′ₘₐₓ in (ℓₘₐₓ, 4)
        # The reference: every ℓ in increasing order, from a single set of rotor data
        sequential = WignerHCalculator(β⃗, ℓₘₐₓ; m′ₘₐₓ)
        reference = [wedge(recurrence!(sequential, ℓ).Hˡ) for ℓ in 0:ℓₘₐₓ]
        @test all(!isempty, reference)

        # Reusing the rotor data.  Every result must be identical to the sequential one:
        # a repeat recomputes from the same axes, a jump forward advances through the
        # intermediate ℓ values, and a move backward restarts from ℓ=0, so exactly the same
        # arithmetic is performed in every case.
        calc = WignerHCalculator(β⃗, ℓₘₐₓ; m′ₘₐₓ)
        for ℓ in order
            @test recurrence!(calc, ℓ) === calc
            @test SphericalFunctions.ℓ(calc) == ℓ
            @test calc.Hˡ.ℓ == ℓ
            @test wedge(calc.Hˡ) == reference[ℓ+1]
        end

        # Supplying the rotor data afresh at each ℓ (which restarts the recurrence)
        for ℓ in order
            recurrence!(calc, β⃗, ℓ)
            @test SphericalFunctions.ℓ(calc) == ℓ
            @test wedge(calc.Hˡ) == reference[ℓ+1]
        end
    end
end


@testitem "WignerHCalculator rotor inputs" begin
    import SphericalFunctions: WignerHCalculator, HWedge, recurrence!
    import Quaternionic: Quaternionic, Rotor, Quaternion
    import DoubleFloats: Double64
    import Random

    # The stored wedge of one rotor for the current ℓ
    wedge(H::HWedge, iᵣ) = [H[iᵣ, m′, m] for m′ in H.m′ₘᵢₙ:H.m′ₘₐₓ for m in abs(m′):H.ℓ]
    maxabsdiff(a, b) = maximum(abs.(a .- b))

    ℓₘₐₓ = 8
    @testset "$T" for T in (Float64, Double64, BigFloat)
        rng = Random.Xoshiro(2024)
        # Angles including both poles; α and γ are irrelevant to H and must not affect it
        β⃗ = T[0; T.(rand(rng, 4)) .* T(π); T(π)]
        Nᵣ = length(β⃗)
        α⃗ = T.(rand(rng, Nᵣ)) .* 2T(π)
        γ⃗ = T.(rand(rng, Nᵣ)) .* 2T(π)
        eⁱᵝ⃗ = cis.(β⃗)
        R⃗ = [Quaternionic.from_euler_angles(α, β, γ) for (α, β, γ) in zip(α⃗, β⃗, γ⃗)]
        # A quaternion that is not a `Rotor` is refused rather than normalized: `Rotor` is
        # what says a quaternion denotes a rotation.  (This used to feed `2 * Quaternion(R)`
        # in and check that the magnitude divided out, which it still does internally.)
        Q⃗ = [2 * Quaternion(R) for R in R⃗]
        @test_throws "Rotations are taken as" WignerHCalculator(Q⃗, ℓₘₐₓ)
        @test_throws "Rotations are taken as" WignerHCalculator(Q⃗[1], ℓₘₐₓ)
        @test eltype(β⃗) === T
        @test eltype(eⁱᵝ⃗) === Complex{T}
        @test eltype(R⃗) === Rotor{T}
        @test !(eltype(Q⃗) <: Rotor)

        # The rotor path computes cosβ and sinβ from the quaternion components rather than
        # from cis(β), so it agrees with the angle path only to a few eps (measured ≤ 3.3 eps
        # at ℓ ≤ 8).
        rotor_atol = 8eps(T)

        calcᵦ = WignerHCalculator(β⃗, ℓₘₐₓ)
        calcₑ = WignerHCalculator(eⁱᵝ⃗, ℓₘₐₓ)
        calcᵣ = WignerHCalculator(R⃗, ℓₘₐₓ)
        batched = Vector{Vector{Vector{T}}}(undef, ℓₘₐₓ + 1)  # batched[ℓ+1][iᵣ]
        for ℓ in 0:ℓₘₐₓ
            recurrence!(calcᵦ, ℓ)
            recurrence!(calcₑ, ℓ)
            recurrence!(calcᵣ, ℓ)
            batched[ℓ+1] = [wedge(calcᵦ.Hˡ, iᵣ) for iᵣ in 1:Nᵣ]
            for iᵣ in 1:Nᵣ
                # The angle is converted to the same phase, so the arithmetic is identical
                @test wedge(calcₑ.Hˡ, iᵣ) == batched[ℓ+1][iᵣ]
                @test maxabsdiff(wedge(calcᵣ.Hˡ, iᵣ), batched[ℓ+1][iᵣ]) ≤ rotor_atol
            end
        end

        # Single-element inputs for Nᵣ=1: a scalar β, a scalar eⁱᵝ, a single Rotor.  Each
        # rotor's slice of the batched result must equal the single-rotor result.
        for (iᵣ, (β, eⁱᵝ, R)) in enumerate(zip(β⃗, eⁱᵝ⃗, R⃗))
            calc₁ᵦ = WignerHCalculator(β, ℓₘₐₓ)
            calc₁ₑ = WignerHCalculator(eⁱᵝ, ℓₘₐₓ)
            calc₁ᵣ = WignerHCalculator(R, ℓₘₐₓ)
            for ℓ in 0:ℓₘₐₓ
                # Supplying the rotor data on every call is allowed; it restarts the recurrence
                recurrence!(calc₁ᵦ, β, ℓ)
                recurrence!(calc₁ₑ, eⁱᵝ, ℓ)
                recurrence!(calc₁ᵣ, R, ℓ)
                @test wedge(calc₁ᵦ.Hˡ, 1) == batched[ℓ+1][iᵣ]
                @test wedge(calc₁ₑ.Hˡ, 1) == batched[ℓ+1][iᵣ]
                @test maxabsdiff(wedge(calc₁ᵣ.Hˡ, 1), batched[ℓ+1][iᵣ]) ≤ rotor_atol
            end
        end

        # Length-1 vectors are also accepted for Nᵣ=1, at construction and later
        calc₁ = WignerHCalculator(β⃗[2:2], ℓₘₐₓ)
        recurrence!(calc₁, ℓₘₐₓ)
        @test wedge(calc₁.Hˡ, 1) == batched[ℓₘₐₓ+1][2]
        recurrence!(calc₁, eⁱᵝ⃗[3:3], ℓₘₐₓ)
        @test wedge(calc₁.Hˡ, 1) == batched[ℓₘₐₓ+1][3]
        recurrence!(calc₁, R⃗[4:4], ℓₘₐₓ)
        @test maxabsdiff(wedge(calc₁.Hˡ, 1), batched[ℓₘₐₓ+1][4]) ≤ rotor_atol
    end
end


@testitem "WignerHCalculator errors" begin
    import SphericalFunctions: WignerHCalculator, recurrence!, wedge_value
    import Quaternionic: Quaternionic

    ℓₘₐₓ = 4
    β⃗ = [0.1, 0.2, 0.3, 0.4]
    R = Quaternionic.from_euler_angles(0.1, 0.2, 0.3)

    # Invalid construction
    @test_throws ErrorException WignerHCalculator(0.3, ℓₘₐₓ; m′ₘₐₓ=ℓₘₐₓ+1)
    @test_throws ErrorException WignerHCalculator(0.3, ℓₘₐₓ; m′ₘₐₓ=-1)
    @test_throws ErrorException WignerHCalculator(0.3, -1)
    # Nᵣ is implied by the rotor data, so an empty batch is how one asks for no rotors
    @test_throws ErrorException WignerHCalculator(Float64[], ℓₘₐₓ)

    # Out-of-range ℓ, with and without fresh rotor data
    calc = WignerHCalculator(0.3, ℓₘₐₓ)
    @test_throws ErrorException recurrence!(calc, 0.3, -1)
    @test_throws ErrorException recurrence!(calc, 0.3, ℓₘₐₓ + 1)
    recurrence!(calc, 0.3, 2)
    @test_throws ErrorException recurrence!(calc, -1)
    @test_throws ErrorException recurrence!(calc, ℓₘₐₓ + 1)
    @test calc.Hˡ.ℓ == 2  # the rejected requests left the calculator where it was

    # Wrong number of rotors
    calc₄ = WignerHCalculator(β⃗, ℓₘₐₓ)
    @test_throws ErrorException recurrence!(calc₄, β⃗[1:3], 0)
    @test_throws ErrorException recurrence!(calc₄, [β⃗; 0.5], 0)
    @test_throws ErrorException recurrence!(calc₄, cis.(β⃗[1:2]), 0)
    @test_throws ErrorException recurrence!(calc₄, fill(R, 3), 0)
    @test_throws ErrorException recurrence!(calc, β⃗[1:2], 0)  # Nᵣ=1 given two

    # A single rotor for a calculator with Nᵣ>1
    @test_throws ErrorException recurrence!(calc₄, 0.3, 0)
    @test_throws ErrorException recurrence!(calc₄, cis(0.3), 0)
    @test_throws ErrorException recurrence!(calc₄, R, 0)

    # Reads outside the stored wedge
    recurrence!(calc, 0.3, 2)
    @test_throws BoundsError calc.Hˡ[1, 2, 1]  # m < |m′| is not stored ...
    @test wedge_value(calc.Hˡ, 1, 2, 1) == calc.Hˡ[1, 1, 2]  # ... but is available by symmetry
    @test_throws BoundsError calc.Hˡ[2, 0, 0]  # iᵣ out of range
    @test_throws BoundsError calc.Hˡ[0, 0, 0]
    @test_throws BoundsError wedge_value(calc.Hˡ, 2, 0, 0)
    @test_throws BoundsError calc.Hˡ[1, 0, 3]  # m > ℓ
    @test_throws BoundsError wedge_value(calc.Hˡ, 1, 0, 3)
    @test_throws BoundsError calc.Hˡ[1, -3, 3]  # |m′| > ℓ
    calc₁ = WignerHCalculator(0.3, ℓₘₐₓ; m′ₘₐₓ=1)
    recurrence!(calc₁, 0.3, ℓₘₐₓ)
    @test_throws BoundsError calc₁.Hˡ[1, 2, 3]  # |m′| > m′ₘₐₓ is not stored ...
    @test_throws ArgumentError wedge_value(calc₁.Hˡ, 1, 2, 3)  # ... nor obtainable by symmetry
    @test wedge_value(calc₁.Hˡ, 1, 3, 1) == calc₁.Hˡ[1, 1, 3]  # unlike |m| ≤ m′ₘₐₓ < |m′|
end


@testitem "WignerHCalculator similar and fill!" begin
    import SphericalFunctions
    import SphericalFunctions: WignerHCalculator, HWedge, recurrence!
    import SphericalFunctions: ℓₘₐₓ, ℓₘᵢₙ, m′ₘₐₓ, m′ₘᵢₙ, Nᵣ
    import DoubleFloats: Double64
    import MathChecker: checked, unchecked

    # Snapshot of the stored wedge for the current ℓ, in storage order
    function wedge(H::HWedge)
        [H[iᵣ, m′, m] for m′ in H.m′ₘᵢₙ:H.m′ₘₐₓ for m in abs(m′):H.ℓ for iᵣ in 1:H.Nᵣ]
    end

    @testset "$T" for T in (Float64, Double64, BigFloat)
        β⃗ = T[0, 0.1, 1.2, 2.9, π]
        calc = WignerHCalculator(β⃗, 6; m′ₘₐₓ=3)

        # `similar` gives the same sizes and types, with fresh storage and the same data
        s = similar(calc)
        @test typeof(s) === typeof(calc)
        @test ℓₘₐₓ(s) == ℓₘₐₓ(calc) == 6
        @test m′ₘₐₓ(s) == m′ₘₐₓ(calc) == 3
        @test m′ₘᵢₙ(s) == m′ₘᵢₙ(calc) == -3
        @test Nᵣ(s) == Nᵣ(calc) == length(β⃗)
        @test ℓₘᵢₙ(s) == ℓₘᵢₙ(calc) == 0
        @test eltype(parent(s.Hˡ)) === T
        @test parent(s.Hˡ) !== parent(calc.Hˡ)
        @test length(parent(s.Hˡ)) == length(parent(calc.Hˡ))

        # The two calculators do not share storage
        fill!(s, NaN)
        recurrence!(calc, 0)
        for ℓ in 1:6
            recurrence!(calc, ℓ)
        end
        @test all(isnan, parent(s.Hˡ))
        reference = copy(parent(calc.Hˡ))  # at ℓ=ℓₘₐₓ every stored entry is in use
        @test !any(isnan, reference)
        recurrence!(s, reverse(β⃗), 6)
        @test parent(calc.Hˡ) == reference
        @test parent(s.Hˡ) != reference

        # Wedges for every ℓ from a fresh calculator.  `similar` retains the rotor data,
        # so `fresh` needs none of its own — and its results must therefore be
        # `calc`'s own values exactly, which every comparison below relies on.
        fresh = similar(calc)
        wedges = [wedge(recurrence!(fresh, ℓ).Hˡ) for ℓ in 0:6]

        # `fill!(NaN)` poisons every buffer, so the recurrence has to rebuild everything it
        # reads; the results must be unchanged, and no NaN may survive in the used region.
        fill!(calc, NaN)
        @test all(isnan, parent(calc.Hˡ))
        @test all(isnan, parent(calc.h⃗ᵃ))
        @test all(isnan, parent(calc.h⃗ᵇ))
        for ℓ in 0:6  # sequential, reusing the rotor data
            recurrence!(calc, ℓ)
            @test wedge(calc.Hˡ) == wedges[ℓ+1]
        end
        @test parent(calc.Hˡ) == reference
        fill!(calc, NaN)
        recurrence!(calc, 6)  # direct jump to ℓₘₐₓ
        @test parent(calc.Hˡ) == reference
        fill!(calc, NaN)
        for ℓ in (4, 2, 6, 0, 5)  # jumps forward and restarts
            recurrence!(calc, ℓ)
            @test wedge(calc.Hˡ) == wedges[ℓ+1]
            @test !any(isnan, wedge(calc.Hˡ))
        end

        # Other fill values, and supplying the rotor data again after a fill
        @test fill!(calc, 7) === calc
        @test all(==(7), parent(calc.Hˡ))
        @test wedge(recurrence!(calc, 3).Hˡ) == wedges[4]
        fill!(calc, NaN)
        @test wedge(recurrence!(calc, β⃗, 5).Hˡ) == wedges[6]
    end

    # Signaling NaNs: with `MathChecker.Checked` storage any arithmetic on a poisoned entry
    # throws, so a run after `fill!(NaN)` proves that no stale or uninitialized value is
    # ever read.  Only the NaN check is enabled; the rest would flag things this test is
    # not about.
    @testset "Checked" begin
        T = Float64
        NC = checked(T; precision=false, nan=true, inf=false)
        β⃗ = T[0, 0.1, 1.2, 2.9, π]
        reference = WignerHCalculator(β⃗, 6; m′ₘₐₓ=3)
        wedges = [wedge(recurrence!(reference, ℓ).Hˡ) for ℓ in 0:6]

        # The calculator's element type is its angles' own, so the checked type is asked for
        # by giving the same angles as `NC` values.
        signaling = WignerHCalculator(NC.(β⃗), 6; m′ₘₐₓ=3)
        @test eltype(parent(signaling.Hˡ)) === NC
        for order in ((0, 1, 2, 3, 4, 5, 6), (6,), (3, 6, 1, 4, 0, 6))
            fill!(signaling, NaN)
            @test all(isnan, parent(signaling.Hˡ))
            for ℓ in order  # `fill!` keeps the rotor data, so none is supplied again here
                recurrence!(signaling, ℓ)  # throws NaNError if any poisoned entry is used
                values = unchecked.(wedge(signaling.Hˡ))
                # `@fastmath` in the Float64 path may contract or reorder operations, so the
                # agreement is not bitwise (measured ≤ 1 eps).
                @test maximum(abs.(values .- wedges[ℓ+1])) ≤ 8eps(T)
            end
        end
    end
end
