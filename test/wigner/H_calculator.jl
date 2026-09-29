# Tests of `HCalculator`, the batched engine that runs the Wigner recurrences and
# produces the H wedge for several rotors at once.  The oracle for its values is a
# closed-form H, built from the definition of H in terms of the Wigner d function; see the
# first item below.

@testitem "HCalculator vs closed-form H" setup=[HalfIntegerOracle, Utilities] begin
    import SphericalFunctions: HCalculator, recurrence!, wedge_value
    import .HalfIntegerOracle: d_oracle
    import .Utilities: βrange
    import DoubleFloats: Double64
    import Random

    # Reference, category 1 (a closed-form formula written out in the test).  H is not a
    # standard published quantity: it is *defined* by
    #
    #     d^ℓ_{m′,m}(β) = ϵ(m′) ϵ(-m) H^ℓ_{m′,m}(β),   ϵ(k) = (-1)^⌊k⌋ for k > 0, else 1
    #
    # (`docs/src/50-notes/01-H_recurrence.md`).  The ϵ are signs, so the relation inverts to
    # H = ϵ(m′) ϵ(-m) d, and `d` is taken from Varshalovich Eq. 4.3.1(2) as transcribed in the
    # `HalfIntegerOracle` setup module — a reference that owes nothing to this package.
    #
    # Every index here is an integer, where the symmetry sign σ = sgn(m′) sgn(m) of
    # `H_recurrence.md` is identically +1, so H_{m′,m} = H_{m,m′} = H_{-m′,-m} = H_{-m,-m′}.
    # The closed form is therefore evaluated on the stored wedge m ≥ |m′| only, and the
    # other elements of the reference are filled from those identities, which is exact for
    # integers.  (The half-integer case, where σ is -1 whenever sgn(m′) ≠ sgn(m), is the
    # next item, whose reference is evaluated element by element.)
    ϵ(k) = ifelse(k > 0 && isodd(k), -1, 1)

    # The closed form is evaluated at twice the working precision and more, and rounded to
    # `T`, so the reference has no more than half an ulp of its own error.  Evaluated *at*
    # the working precision it would be useless as an oracle: its alternating sum loses
    # about 270 eps (8 bits) at ℓ = 16 in Float64.
    function Href(::Type{T}, ℓ, m′, m, β) where {T<:Real}
        r = setprecision(BigFloat, 2 * precision(T) + 64) do
            ϵ(m′) * ϵ(-m) * d_oracle(ℓ, m′, m, BigFloat(β))
        end
        T(r)
    end
    # The whole reference block, `Hblock(T, ℓ, β)[m′+ℓ+1, m+ℓ+1]`
    function Hblock(::Type{T}, ℓ, β) where {T<:Real}
        B = Matrix{T}(undef, 2ℓ + 1, 2ℓ + 1)
        for m′ in -ℓ:ℓ, m in abs(m′):ℓ
            B[m′+ℓ+1, m+ℓ+1] = Href(T, ℓ, m′, m, β)
        end
        for m′ in -ℓ:ℓ, m in -ℓ:ℓ
            m ≥ abs(m′) && continue
            a, b = m′ ≥ abs(m) ? (m, m′) : -m′ ≥ abs(m) ? (-m, -m′) : (-m′, -m)
            B[m′+ℓ+1, m+ℓ+1] = B[a+ℓ+1, b+ℓ+1]
        end
        B
    end

    ℓᵗᵒᵖ = 16  # the largest ℓ reached by any configuration below

    @testset "$T" for T in (Float64, Double64, BigFloat)
        rng = Random.Xoshiro(1234)
        # Both poles, their nearest neighbors, and four random angles
        β⃗ = βrange(rng, T, 4)
        @test length(β⃗) == 8
        @test iszero(first(β⃗)) && last(β⃗) == T(π)

        # H^ℓ_{m′,m}(β) depends on nothing but ℓ, m′, m and β, so one table of references
        # serves every (ℓₘₐₓ, m′ₘₐₓ, Nᵣ) configuration below: `ref[k][ℓ+1][m′+ℓ+1, m+ℓ+1]`
        # is H^ℓ_{m′,m}(β⃗[k]).
        ref = [[Hblock(T, ℓ, β) for ℓ in 0:ℓᵗᵒᵖ] for β in β⃗]

        # Measured worst case over every configuration below, with no growth in ℓ:
        # 2.5 eps (Float64), 2.4 eps (Double64), 3.0 eps (BigFloat) at ℓ ≤ 16.
        atol = 10 * eps(T)

        @testset "ℓₘₐₓ=$ℓₘₐₓ" for ℓₘₐₓ in (0, 1, 2, 3, 7, 16)
            m′ₘₐₓs = ℓₘₐₓ ≤ 7 ? (0:ℓₘₐₓ) : (0, 1, 2, 5, 8, 16)
            @testset "m′ₘₐₓ=$m′ₘₐₓ" for m′ₘₐₓ in m′ₘₐₓs
                for Nᵣ in (1, 4)
                    calc = HCalculator(β⃗[1:Nᵣ], ℓₘₐₓ; m′ₘₐₓ)
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


@testitem "HCalculator: half-integer wedge_value vs closed-form H" setup=[HalfIntegerOracle] begin
    import SphericalFunctions: HCalculator, recurrence!, wedge_value, wedge_source,
        HalfOddInteger
    import .HalfIntegerOracle: d_oracle

    # `wedge_value` is the way the documentation gives to read the elements of H that the
    # wedge does not store, and for half-integer indices it has to apply the sign σ, which is
    # -1 for half of the elements that come from a transposition.  So every element that the
    # wedge can supply is compared here with the definition H = ϵ(m′) ϵ(-m) d, with d from
    # Varshalovich's closed form, evaluated for each element separately rather than filled in
    # by a symmetry, so that the reference owes nothing to σ.  ϵ(k) = (-1)^⌊k⌋ for k > 0.
    ϵ(k) = ifelse(k > 0 && isodd(floor(Int, k)), -1, 1)
    Href(ℓ, m′, m, β) = Float64(ϵ(m′) * ϵ(-m) * d_oracle(ℓ, m′, m, big(β)))

    # Measured worst case over everything below: 2.0 eps, over 12684 values of which 453
    # come from an image with σ = -1; 10 eps is asserted.
    atol = 10eps()
    βs = [0.0, 1.0e-3, 0.4, 1.1, 2.9, π - 1.0e-9, Float64(π)]
    worst, σ₋ = let worst = 0.0, σ₋ = 0
        for ℓₘₐₓ ∈ (1//2, 7//2, 15//2), m′ₘₐₓ ∈ unique((1//2, 3//2, ℓₘₐₓ))
            m′ₘₐₓ ≤ ℓₘₐₓ || continue
            calc = HCalculator(βs, ℓₘₐₓ; m′ₘₐₓ)
            for ℓ ∈ 1//2:1:ℓₘₐₓ
                H = recurrence!(calc, ℓ)
                for m′ ∈ -ℓ:ℓ, m ∈ -ℓ:ℓ
                    min(abs(m′), abs(m)) ≤ m′ₘₐₓ || continue
                    σ₋ += last(wedge_source(HalfOddInteger(m′), HalfOddInteger(m), H.m′ₘₐₓ)) == -1
                    for (iᵣ, β) ∈ enumerate(βs)
                        worst = max(worst, abs(wedge_value(H, iᵣ, m′, m) - Href(ℓ, m′, m, β)))
                    end
                end
            end
        end
        (worst, σ₋)
    end
    @test worst ≤ atol
    @test σ₋ == 453  # the comparison reaches the elements whose σ is -1
end


@testitem "HCalculator ℓ ordering" begin
    import SphericalFunctions
    import SphericalFunctions: HCalculator, HWedge, recurrence!
    import Random

    # Snapshot of the stored wedge for the current ℓ, in storage order
    function wedge(H::HWedge)
        [H[iᵣ, m′, m] for m′ in -H.m′ₘₐₓ:H.m′ₘₐₓ for m in abs(m′):H.ℓ for iᵣ in 1:H.Nᵣ]
    end

    T = Float64
    ℓₘₐₓ = 9
    rng = Random.Xoshiro(3)
    β⃗ = T[0; rand(rng, 2) .* T(π); T(π)]

    # Repeats, forward jumps, and backward moves, including returns to 0 and to ℓₘₐₓ
    order = (0, 3, 1, 9, 9, 5, 6, 0, 2, 8, 4, 9, 7, 7, 0, 9)

    for m′ₘₐₓ in (ℓₘₐₓ, 4)
        # The reference: every ℓ in increasing order, from a single set of rotor data
        sequential = HCalculator(β⃗, ℓₘₐₓ; m′ₘₐₓ)
        reference = [wedge(recurrence!(sequential, ℓ)) for ℓ in 0:ℓₘₐₓ]
        @test all(!isempty, reference)

        # A copy of the wedge is independent of the calculator, and survives its next step
        snapshot = copy(recurrence!(sequential, 3))
        recurrence!(sequential, 5)
        @test wedge(snapshot) == reference[4]

        # Reusing the rotor data.  Every result must be identical to the sequential one:
        # a repeat recomputes from the same axes, a jump forward advances through the
        # intermediate ℓ values, and a move backward restarts from ℓ=0, so exactly the same
        # arithmetic is performed in every case.
        calc = HCalculator(β⃗, ℓₘₐₓ; m′ₘₐₓ)
        for ℓ in order
            # The wedge comes back by identity — it is one mutable object, not a view
            @test recurrence!(calc, ℓ) === calc.Hˡ
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

    # The same for half-integer ℓ, whose axis runs at the integer order ℓ - 1/2, with both
    # poles among the angles
    ℓₘₐₓₕ = 19//2
    orderₕ = (1//2, 7//2, 3//2, 19//2, 19//2, 11//2, 13//2, 1//2, 5//2, 17//2)
    β⃗ₕ = T[0; rand(rng, 3) .* T(π); T(π)]
    for m′ₘₐₓ in (ℓₘₐₓₕ, 5//2, 1//2)
        sequential = HCalculator(β⃗ₕ, ℓₘₐₓₕ; m′ₘₐₓ)
        reference = Dict(ℓ => wedge(recurrence!(sequential, ℓ)) for ℓ in 1//2:1:ℓₘₐₓₕ)
        calc = HCalculator(β⃗ₕ, ℓₘₐₓₕ; m′ₘₐₓ)
        for ℓ in orderₕ
            @test recurrence!(calc, ℓ) === calc.Hˡ
            @test SphericalFunctions.ℓ(calc) == ℓ
            @test wedge(calc.Hˡ) == reference[ℓ]
        end
        for ℓ in orderₕ
            recurrence!(calc, β⃗ₕ, ℓ)
            @test wedge(calc.Hˡ) == reference[ℓ]
        end
    end
end


@testitem "HCalculator rotor inputs" begin
    import SphericalFunctions: HCalculator, HWedge, recurrence!
    import Quaternionic: Quaternionic, Rotor, Quaternion
    import DoubleFloats: Double64
    import Random

    # The stored wedge of one rotor for the current ℓ
    wedge(H::HWedge, iᵣ) = [H[iᵣ, m′, m] for m′ in -H.m′ₘₐₓ:H.m′ₘₐₓ for m in abs(m′):H.ℓ]
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
        # what says a quaternion denotes a rotation.  (A `Rotor` whose magnitude is not 1 is
        # accepted, and the magnitude divides out; see the item "Calculators accept
        # unnormalized rotors, and refuse phases that are not".)
        Q⃗ = [2 * Quaternion(R) for R in R⃗]
        @test_throws "Rotations are taken as" HCalculator(Q⃗, ℓₘₐₓ)
        @test_throws "Rotations are taken as" HCalculator(Q⃗[1], ℓₘₐₓ)
        @test eltype(β⃗) === T
        @test eltype(eⁱᵝ⃗) === Complex{T}
        @test eltype(R⃗) === Rotor{T}
        @test !(eltype(Q⃗) <: Rotor)

        # The rotor path computes cosβ and sinβ from the quaternion components rather than
        # from cis(β), so it agrees with the angle path only to a few eps (measured ≤ 3.3 eps
        # at ℓ ≤ 8).
        rotor_atol = 8eps(T)

        calcᵦ = HCalculator(β⃗, ℓₘₐₓ)
        calcₑ = HCalculator(eⁱᵝ⃗, ℓₘₐₓ)
        calcᵣ = HCalculator(R⃗, ℓₘₐₓ)
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
            calc₁ᵦ = HCalculator(β, ℓₘₐₓ)
            calc₁ₑ = HCalculator(eⁱᵝ, ℓₘₐₓ)
            calc₁ᵣ = HCalculator(R, ℓₘₐₓ)
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
        calc₁ = HCalculator(β⃗[2:2], ℓₘₐₓ)
        recurrence!(calc₁, ℓₘₐₓ)
        @test wedge(calc₁.Hˡ, 1) == batched[ℓₘₐₓ+1][2]
        recurrence!(calc₁, eⁱᵝ⃗[3:3], ℓₘₐₓ)
        @test wedge(calc₁.Hˡ, 1) == batched[ℓₘₐₓ+1][3]
        recurrence!(calc₁, R⃗[4:4], ℓₘₐₓ)
        @test maxabsdiff(wedge(calc₁.Hˡ, 1), batched[ℓₘₐₓ+1][4]) ≤ rotor_atol
    end
end


@testitem "HCalculator errors" setup=[RefusalChecks] begin
    import SphericalFunctions: HCalculator, recurrence!, wedge_value, HalfOddInteger
    import Quaternionic: Quaternionic

    ℓₘₐₓ = 4
    β⃗ = [0.1, 0.2, 0.3, 0.4]
    R = Quaternionic.from_euler_angles(0.1, 0.2, 0.3)

    # Invalid construction
    m′range = "must satisfy 0 ≤ m′ₘₐₓ ≤ ℓₘₐₓ"
    @test refuses(() -> HCalculator(0.3, ℓₘₐₓ; m′ₘₐₓ=ℓₘₐₓ+1), ArgumentError, m′range)
    @test refuses(() -> HCalculator(0.3, ℓₘₐₓ; m′ₘₐₓ=-1), ArgumentError, m′range)
    # A bad ℓₘₐₓ is named as such, rather than as the `m′ₘₐₓ` that defaults to it
    @test refuses(() -> HCalculator(0.3, -1), ArgumentError, "ℓₘₐₓ=-1 must be non-negative")
    @test refuses(
        () -> HCalculator(0.3, -1//2), ArgumentError, "ℓₘₐₓ=-1//2 must be non-negative"
    )
    @test refuses(
        () -> HCalculator(0.3, 7//2; m′ₘₐₓ=9//2), ArgumentError,
        "must satisfy 1//2 ≤ m′ₘₐₓ ≤ ℓₘₐₓ"
    )
    # Nᵣ is implied by the rotor data, and an empty batch would be a request for no rotors,
    # which is refused.  (`nrotors` owns this message and its exception type.)
    @test_throws "at least one rotor" HCalculator(Float64[], ℓₘₐₓ)

    # The indices are of one kind, and `Int` or half-odd-integers, in any spelling of either
    @test refuses(() -> HCalculator(0.3, 7//2; m′ₘₐₓ=1), ArgumentError, "keyword argument `m′ₘₐₓ`")
    @test refuses(() -> HCalculator(0.3, Int16(4)), ArgumentError, "narrower than `Int`")
    @test refuses(() -> HCalculator(0.3, 4; m′ₘₐₓ=Int8(2)), ArgumentError, "narrower than `Int`")

    # Out-of-range ℓ, with and without fresh rotor data
    calc = HCalculator(0.3, ℓₘₐₓ)
    @test refuses(() -> recurrence!(calc, 0.3, -1), ArgumentError, "out of bounds")
    @test refuses(() -> recurrence!(calc, 0.3, ℓₘₐₓ + 1), ArgumentError, "out of bounds")
    recurrence!(calc, 0.3, 2)
    @test refuses(() -> recurrence!(calc, -1), ArgumentError, "out of bounds")
    @test refuses(() -> recurrence!(calc, ℓₘₐₓ + 1), ArgumentError, "out of bounds")
    # ... and an ℓ that is not an index of this calculator's kind at all
    @test refuses(() -> recurrence!(calc, 2.0), ArgumentError, "so ℓ must be one too")
    @test refuses(() -> recurrence!(calc, 3//2), ArgumentError, "so ℓ must be one too")
    @test refuses(() -> recurrence!(calc, 2//1), ArgumentError, "so ℓ must be one too")
    @test refuses(() -> recurrence!(calc, true), ArgumentError, "so ℓ must be one too")
    @test calc.Hˡ.ℓ == 2  # the rejected requests left the calculator where it was
    # An integer of another type is the same index
    @test recurrence!(calc, Int8(3)) == recurrence!(HCalculator(0.3, ℓₘₐₓ), 3)

    # Wrong number of rotors
    calc₄ = HCalculator(β⃗, ℓₘₐₓ)
    @test refuses(() -> recurrence!(calc₄, β⃗[1:3], 0), DimensionMismatch, "Expected 4 rotors")
    @test refuses(() -> recurrence!(calc₄, [β⃗; 0.5], 0), DimensionMismatch, "Expected 4 rotors")
    @test refuses(
        () -> recurrence!(calc₄, cis.(β⃗[1:2]), 0), DimensionMismatch, "Expected 4 rotors"
    )
    @test refuses(() -> recurrence!(calc₄, fill(R, 3), 0), DimensionMismatch, "Expected 4 rotors")
    @test refuses(() -> recurrence!(calc, β⃗[1:2], 0), DimensionMismatch, "Expected 1 rotors")

    # A single rotor for a calculator with Nᵣ>1
    for single ∈ (0.3, cis(0.3), R)
        @test refuses(
            () -> recurrence!(calc₄, single, 0), DimensionMismatch, "A single rotor was given"
        )
    end

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
    calc₁ = HCalculator(0.3, ℓₘₐₓ; m′ₘₐₓ=1)
    recurrence!(calc₁, 0.3, ℓₘₐₓ)
    @test_throws BoundsError calc₁.Hˡ[1, 2, 3]  # |m′| > m′ₘₐₓ is not stored ...
    @test refuses(  # ... nor obtainable by symmetry
        () -> wedge_value(calc₁.Hˡ, 1, 2, 3), ArgumentError, "both |m′| and |m| exceed"
    )
    @test wedge_value(calc₁.Hˡ, 1, 3, 1) == calc₁.Hˡ[1, 1, 3]  # unlike |m| ≤ m′ₘₐₓ < |m′|
    # The indices of `wedge_value` must be of the wedge's kind, in any spelling of it
    @test wedge_value(calc₁.Hˡ, 1, Int8(3), 1) == wedge_value(calc₁.Hˡ, 1, 3, 1)
    @test refuses(() -> wedge_value(calc₁.Hˡ, 1, 1//2, 1), ArgumentError, "so m′ must be one too")
    Hₕ = recurrence!(HCalculator(0.3, 7//2), 7//2)
    @test wedge_value(Hₕ, 1, 1//2, -3//2) ==
        wedge_value(Hₕ, 1, HalfOddInteger(1//2), HalfOddInteger(-3//2))
    @test refuses(() -> wedge_value(Hₕ, 1, 1, 2), ArgumentError, "so m′ must be one too")
    @test refuses(() -> wedge_value(Hₕ, 1, 1//2, 1//1), ArgumentError, "so m must be one too")
end


@testitem "HCalculator: either spelling of the keyword, and of a half-integer" begin
    import SphericalFunctions: HCalculator, recurrence!, HalfOddInteger
    import SphericalFunctions

    # `mp_max` is the ASCII spelling of `m′ₘₐₓ`, and a half-integer keyword may be spelled as a
    # `Rational` whatever the spelling of ℓₘₐₓ
    for (ℓₘₐₓ, m′ₘₐₓ, ℓ) ∈ ((6, 2, 5), (HalfOddInteger(13//2), 3//2, 9//2), (13//2, HalfOddInteger(3//2), 9//2))
        reference = recurrence!(HCalculator([0.3, 1.1], ℓₘₐₓ; m′ₘₐₓ), ℓ)
        @test recurrence!(HCalculator([0.3, 1.1], ℓₘₐₓ; mp_max=m′ₘₐₓ), ℓ) == reference
        @test SphericalFunctions.m′ₘₐₓ(HCalculator(0.3, ℓₘₐₓ; mp_max=m′ₘₐₓ)) == m′ₘₐₓ
    end
end


@testitem "Calculators accept unnormalized rotors, and refuse phases that are not" setup=[RefusalChecks] begin
    import SphericalFunctions: D, d, sYlm, DCalculator, dCalculator, HCalculator, recurrence!,
        set_β!, ℓ, array_view
    import Quaternionic: Rotor, Quaternion, rotor

    # A `Rotor` whose magnitude is not 1 is representable, and denotes the same rotation as the
    # normalized one: `spinor_phases` divides the magnitude out.  Measured differences at most
    # 3.6e-16; 4 eps is asserted.
    Ru = Rotor{Float64}(0.6, 0.8, 0.4, 0.2)
    @test abs2(Quaternion(Ru)) ≈ 1.2
    Rn = rotor(Quaternion(Ru))
    ϵ = 4eps()
    for ℓₘₐₓ ∈ (3, 5//2)
        @test maximum(abs, array_view(D(Ru, ℓₘₐₓ)[ℓₘₐₓ]) - array_view(D(Rn, ℓₘₐₓ)[ℓₘₐₓ])) ≤ ϵ
        @test maximum(abs, array_view(d(Ru, ℓₘₐₓ)[ℓₘₐₓ]) - array_view(d(Rn, ℓₘₐₓ)[ℓₘₐₓ])) ≤ ϵ
        s = ℓₘₐₓ isa Integer ? 1 : 1//2
        @test maximum(abs, array_view(sYlm(Ru, ℓₘₐₓ, s)) - array_view(sYlm(Rn, ℓₘₐₓ, s))) ≤ ϵ
        calc = DCalculator(Rn, ℓₘₐₓ)
        @test maximum(abs, array_view(recurrence!(calc, Ru, ℓₘₐₓ)) -
            array_view(recurrence!(DCalculator(Rn, ℓₘₐₓ), ℓₘₐₓ))) ≤ ϵ
    end

    # A phase must have unit modulus.  The usual mistake is to pass β itself as a complex
    # number, which would otherwise give silently wrong values; a NaN is refused as well.
    for bad ∈ (complex(0.7), complex(NaN, NaN), complex(1.0, NaN))
        @test refuses(() -> d(bad, 3), DomainError, "must have unit modulus")
        @test refuses(() -> dCalculator(bad, 3), DomainError, "must have unit modulus")
        @test refuses(() -> HCalculator([cis(0.1), bad], 3), DomainError, "must have unit modulus")
    end
    # ... and a refused phase leaves the calculator, and what it had computed, as they were
    calc = dCalculator(0.3, 3)
    before = copy(recurrence!(calc, 2))
    @test refuses(() -> set_β!(calc, complex(0.7)), DomainError, "must have unit modulus")
    @test ℓ(calc) == 2
    @test recurrence!(calc, 2) == before
end



@testitem "HCalculator similar and fill!" begin
    import SphericalFunctions
    import SphericalFunctions: HCalculator, HWedge, recurrence!
    import SphericalFunctions: ℓₘₐₓ, ℓₘᵢₙ, m′ₘₐₓ, m′ₘᵢₙ, Nᵣ
    import DoubleFloats: Double64
    import MathChecker: checked, unchecked

    # Snapshot of the stored wedge for the current ℓ, in storage order
    function wedge(H::HWedge)
        [H[iᵣ, m′, m] for m′ in -H.m′ₘₐₓ:H.m′ₘₐₓ for m in abs(m′):H.ℓ for iᵣ in 1:H.Nᵣ]
    end

    @testset "$T" for T in (Float64, Double64, BigFloat)
        β⃗ = T[0, 0.1, 1.2, 2.9, π]
        calc = HCalculator(β⃗, 6; m′ₘₐₓ=3)

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
        wedges = [wedge(recurrence!(fresh, ℓ)) for ℓ in 0:6]

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
        @test wedge(recurrence!(calc, 3)) == wedges[4]
        fill!(calc, NaN)
        @test wedge(recurrence!(calc, β⃗, 5)) == wedges[6]
    end

    # Signaling NaNs: with `MathChecker.Checked` storage any arithmetic on a poisoned entry
    # throws, so a run after `fill!(NaN)` proves that no stale or uninitialized value is
    # ever read.  Only the NaN check is enabled; the rest would flag things this test is
    # not about.
    @testset "Checked" begin
        T = Float64
        NC = checked(T; precision=false, nan=true, inf=false)
        β⃗ = T[0, 0.1, 1.2, 2.9, π]
        reference = HCalculator(β⃗, 6; m′ₘₐₓ=3)
        wedges = [wedge(recurrence!(reference, ℓ)) for ℓ in 0:6]

        # The calculator's element type is its angles' own, so the checked type is asked for
        # by giving the same angles as `NC` values.
        signaling = HCalculator(NC.(β⃗), 6; m′ₘₐₓ=3)
        @test eltype(parent(signaling.Hˡ)) === NC
        for order in ((0, 1, 2, 3, 4, 5, 6), (6,), (3, 6, 1, 4, 0, 6))
            fill!(signaling, NaN)
            @test all(isnan, parent(signaling.Hˡ))
            for ℓ in order  # `fill!` keeps the rotor data, so none is supplied again here
                recurrence!(signaling, ℓ)  # throws NaNError if any poisoned entry is used
                values = unchecked.(wedge(signaling.Hˡ))
                # The `muladd` in `complex_powers!` is fused for `Float64` but not for the
                # wrapper type, so the agreement is not bitwise (measured ≤ 1 eps).
                @test maximum(abs.(values .- wedges[ℓ+1])) ≤ 8eps(T)
            end
        end
    end
end

@testitem "HCalculator: ℓ is ℓₘᵢₙ-1 until something is computed for the current data" begin
    import SphericalFunctions: HCalculator, dCalculator, recurrence!, set_β!, ℓ, ℓₘᵢₙ

    # As documented for every calculator, and as `dCalculator` does: `ℓ` is the order of the
    # block most recently computed, and `ℓₘᵢₙ - 1` when nothing has been computed since the
    # calculator was built or its data were last replaced.  The layout of the wedge would give
    # a different answer — ℓₘᵢₙ when fresh, and the old ℓ after `set_β!` or `fill!`, while
    # `show` says that nothing was computed and the wedge holds stale numbers.
    for (β, ℓₘₐₓ, ℓ₂) ∈ ((0.3, 4, 2), (0.3, 9//2, 5//2))
        for calc ∈ (HCalculator(β, ℓₘₐₓ), dCalculator(β, ℓₘₐₓ))
            @test ℓ(calc) == ℓₘᵢₙ(calc) - 1
            recurrence!(calc, ℓ₂)
            @test ℓ(calc) == ℓ₂
            set_β!(calc, 0.4)
            @test ℓ(calc) == ℓₘᵢₙ(calc) - 1
            recurrence!(calc, ℓ₂)
            @test ℓ(calc) == ℓ₂
            fill!(calc, 0)
            @test ℓ(calc) == ℓₘᵢₙ(calc) - 1
            recurrence!(calc, 0.5, ℓ₂)
            @test ℓ(calc) == ℓ₂
        end
    end
    @test occursin("nothing computed yet", sprint(show, HCalculator(0.3, 4)))
end

@testitem "Rotor-data setters: a refused angle leaves the calculator as it was" setup=[RefusalChecks] begin
    import SphericalFunctions
    import SphericalFunctions: HCalculator, dCalculator, sYlmCalculator, recurrence!
    import SphericalFunctions: set_β!, set_θ!

    # Every value is validated before anything is stored, so a setter that refuses its input
    # leaves the rotor data, and whatever the calculator had computed from them, exactly as
    # they were, and the next step continues from them.  The infinite angle is that of the
    # second rotor, so that the first would already have been stored if the angles were not
    # all checked first.
    β = [0.3, 0.5]
    for bad ∈ ([0.9, Inf], [0.9, -Inf], [0.9, NaN])
        for (ℓₘₐₓ, ℓ) ∈ ((6, 3), (13//2, 5//2))
            c = HCalculator(β, ℓₘₐₓ)
            recurrence!(c, ℓ)
            @test refuses(() -> set_β!(c, bad), DomainError, "so it has no phase")
            @test recurrence!(c, ℓ + 1) == recurrence!(HCalculator(β, ℓₘₐₓ), ℓ + 1)
        end

        # ... and the same through the calculators built on it
        c = dCalculator(β, 6)
        recurrence!(c, 3)
        @test refuses(() -> set_β!(c, bad), DomainError, "so it has no phase")
        @test copy(recurrence!(c, 4)) == copy(recurrence!(dCalculator(β, 6), 4))
        for Calculator ∈ (SphericalFunctions.sλlmCalculator, sYlmCalculator)
            c = Calculator(β, 6, 1)
            recurrence!(c, 3)
            @test refuses(() -> set_θ!(c, bad), DomainError, "so it has no phase")
            @test copy(recurrence!(c, 4)) == copy(recurrence!(Calculator(β, 6, 1), 4))
        end
    end
end

@testitem "HCalculator: axes labelled out of step are rebuilt rather than trusted" begin
    import SphericalFunctions as SF
    import SphericalFunctions: HCalculator, SSHT, recurrence!, array_view
    import SphericalFunctions: nmodes  # unexported
    using Random

    rng = Random.Xoshiro(1848)

    # The two axis buffers hold consecutive orders whenever the calculator's data are valid.
    # A transform used from several tasks at once, which the documentation of `SSHT` warns
    # against, can leave them labelled otherwise; the recurrence then starts again from the
    # beginning rather than stepping on from them, so that the misuse does not outlast itself.
    for (ℓₘₐₓ, ℓ) ∈ ((6, 5), (13//2, 9//2))
        c = HCalculator([0.3, 0.5], ℓₘₐₓ)
        recurrence!(c, ℓ)
        SF.h⃗ˡ⁺¹(c).ℓ = 1
        @test recurrence!(c, ℓ) == recurrence!(HCalculator([0.3, 0.5], ℓₘₐₓ), ℓ)
        @test recurrence!(c, ℓ + 1) == recurrence!(HCalculator([0.3, 0.5], ℓₘₐₓ), ℓ + 1)
    end

    # Each transform starts at ℓ = |s|, that is at the axis of order ⌊|s|⌋
    for (s, ℓₘₐₓ) ∈ ((2, 12), (1//2, 13//2))
        𝒯 = SSHT(s, ℓₘₐₓ)
        f̃ = randn(rng, ComplexF64, nmodes(𝒯))
        ϵ = 500 * eps()
        @test array_view(𝒯 \ (𝒯 * f̃)) ≈ f̃ atol=ϵ rtol=ϵ
        H = 𝒯.λ.H
        SF.h⃗ˡ(H).ℓ = SF.axis_ℓ(H, abs(𝒯.s))
        SF.h⃗ˡ⁺¹(H).ℓ = 5
        @test array_view(𝒯 \ (𝒯 * f̃)) ≈ f̃ atol=ϵ rtol=ϵ
        @test array_view(𝒯 \ (𝒯 * f̃)) ≈ f̃ atol=ϵ rtol=ϵ
    end
end


@testitem "HCalculator: the table of coefficients is refilled at every step" setup=[RefusalChecks] begin
    import SphericalFunctions: HCalculator, HAxis, HWedge, recurrence!, ℓₘᵢₙ, ℓₘₐₓ
    import SphericalFunctions: FixedSizeVector

    # Steps 4 and 5 read their m-side coefficients √δ²(ℓ, m) from a table, which every
    # `recurrence!` fills for its own ℓ before either step runs.  The table therefore has no
    # state that could go out of date: poisoning it before each step, whether the step
    # advances by one ℓ, jumps ahead, or restarts from a smaller ℓ, changes nothing.
    snapshot(H) = [H[iᵣ, m′, m] for iᵣ ∈ 1:H.Nᵣ for m′ ∈ -H.m′ₘₐₓ:H.m′ₘₐₓ for m ∈ abs(m′):H.ℓ]
    for (ℓmax, β) ∈ ((12, [0.3, 1.7, 2.9]), (23//2, [0.3, 1.7]), (12, 0.9), (23//2, 2.2))
        clean = HCalculator(β, ℓmax)
        poisoned = HCalculator(β, ℓmax)
        ℓs = collect(ℓₘᵢₙ(clean):ℓₘₐₓ(clean))
        for ℓ ∈ [ℓs; ℓs[end-3]; ℓs[2]; ℓs[end]]
            reference = snapshot(recurrence!(clean, ℓ))
            fill!(poisoned.d̄ₗ, NaN)
            @test isequal(snapshot(recurrence!(poisoned, ℓ)), reference)
        end
    end

    # The steps read the table under `@inbounds`, for every m below the largest ℓ of the
    # wedge, so the calculator's own constructor refuses one that is too short.
    Hˡ = HWedge(Float64, 1, 4)
    h⃗ᵃ, h⃗ᵇ = HAxis(Float64, 1, 5), HAxis(Float64, 1, 5)
    h⃗ᵇ.ℓ = 1
    eⁱᵝ = FixedSizeVector{ComplexF64}(undef, 1)
    eⁱᵝ[1] = cis(0.3)
    none = FixedSizeVector{Float64}(undef, 0)
    build(n) = HCalculator{Int, Float64, typeof(parent(Hˡ))}(
        h⃗ᵃ, h⃗ᵇ, Hˡ, eⁱᵝ, none, none, FixedSizeVector{Float64}(undef, n), 4, 4, Ref(false),
        Ref(false)
    )
    @test refuses(() -> build(3), DimensionMismatch, "the table of coefficients 4 entries")
    @test recurrence!(build(4), 4) === Hˡ
end
