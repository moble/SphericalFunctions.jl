# Tests for the Wigner-matrix calculators: `DCalculator`, `dCalculator`, the blocks
# `recurrence!` returns, the four block limits, batched rotors, and the convenience functions
# `D` and `d`.  The oracles are independent closed forms — Varshalovich Eq. 4.3.1(2) for `d`,
# the settled Euler factorization for `𝔇`, and the quaternionic form of Boyle (2016), all
# transcribed in the `HalfIntegerOracle` setup module — together with the explicit and
# formulaic matrices in `ExplicitWignerMatrices` and a set of metamorphic identities.

@testitem "DCalculator vs closed forms" setup=[HalfIntegerOracle, Utilities] begin
    import SphericalFunctions: DCalculator, recurrence!
    import .HalfIntegerOracle: d_oracle, D_oracle
    import .Utilities: Rrange
    using Quaternionic: Rotor, Quaternion, components, 𝐢, 𝐣, 𝐤
    using DoubleFloats: Double64
    using Random

    rng = Random.Xoshiro(1234)

    # Euler angles of a (possibly unnormalized) quaternion, in high precision.  Writing
    # R = exp(α𝐤/2) exp(β𝐣/2) exp(γ𝐤/2) and splitting it into the two complex parts
    # Rₛ = w + z𝑖 = ‖R‖ cos(β/2) e^{i(α+γ)/2} and Rₐ = y + x𝑖 = ‖R‖ sin(β/2) e^{i(γ-α)/2}
    # gives the angles below.  `Quaternionic.to_euler_angles` would do the same job, but it
    # reaches `β` through `acos`, which throws on a rotor that is normalized only to the
    # *working* precision; `atan(|Rₐ|, |Rₛ|)` is scale-free and needs no normalization.
    function euler_angles(R)
        w, x, y, z = BigFloat.(components(Quaternion(R)))
        ϕₛ, ϕₐ = angle(Complex(w, z)), angle(Complex(y, x))
        (ϕₛ - ϕₐ, 2 * atan(abs(Complex(y, x)), abs(Complex(w, z))), ϕₛ + ϕₐ)
    end

    # Reference 1, category 1 (a closed-form formula, from the conventions pages).  The
    # settled convention is 𝔇ˡ_{m′m}(R) = e^{-i m′ α} dˡ_{m′m}(β) e^{-i m γ}
    # (docs/src/30-conventions/01-summary.md), and `d` is Varshalovich Eq. 4.3.1(2) as
    # transcribed in `HalfIntegerOracle`.  It is evaluated at four times the working
    # precision and rounded, so the reference itself contributes at most half an ulp; at
    # the working precision the alternating sum in `d` would be worthless as an oracle.
    function Dref(::Type{T}, R, ℓₘₐₓ) where {T}
        blocks = setprecision(BigFloat, 4 * precision(T) + 64) do
            α, β, γ = euler_angles(R)
            [
                Complex{BigFloat}[
                    cis(-m′ * α) * d_oracle(ℓ, m′, m, β) * cis(-m * γ)
                    for m′ in -ℓ:ℓ, m in -ℓ:ℓ
                ]
                for ℓ in 0:ℓₘₐₓ
            ]
        end
        [Complex{T}.(b) for b in blocks]  # `blocks[ℓ+1][m′+ℓ+1, m+ℓ+1]`
    end

    # Reference 2, category 1 as well, and Float64 only: the quaternionic closed form of
    # Boyle (2016), which reaches 𝔇 straight from the components of R with no Euler
    # decomposition and no trigonometry, so it pins the convention independently of the
    # factorization above.  Its normalization constant is a ratio of `Int` factorials and
    # is therefore computed in Float64 whatever the element type, which is why it is used
    # only for Float64 and with a looser tolerance.
    @testset "$T" for T in (Float64, Double64)
        # Measured worst errors at ℓ ≤ 8: 7.1 eps (Float64) and 2.3 eps (Double64) against
        # reference 1, and 16.2 eps against reference 2.
        atol = 32 * eps(T)
        atolᵇ = 64 * eps(T)
        for ℓₘₐₓ in (0, 1, 2, 4, 8)
            rotors = Rrange(rng, T, 6)
            calc = DCalculator(first(rotors), ℓₘₐₓ)
            for R in rotors
                worst = zero(T)
                worstᵇ = zero(T)
                ref = Dref(T, R, ℓₘₐₓ)
                for ℓ in 0:ℓₘₐₓ
                    𝔇ˡ = if ℓ == 0
                        recurrence!(calc, R, ℓ)  # set the rotor and compute ℓ=0
                    else
                        recurrence!(calc, ℓ)  # reuse the rotor data
                    end
                    @test axes(𝔇ˡ) == (-ℓ:ℓ, -ℓ:ℓ)
                    @test eltype(𝔇ˡ) === Complex{T}
                    # The block is a view into the calculator's storage; stripping the
                    # offsets leaves a 1-based (2ℓ+1)×(2ℓ+1) matrix, and `collect` gives a
                    # plain `Matrix`
                    @test parent(𝔇ˡ) isa AbstractMatrix{Complex{T}}
                    @test axes(parent(𝔇ˡ)) == (1:2ℓ+1, 1:2ℓ+1)
                    @test collect(𝔇ˡ) isa Matrix{Complex{T}}
                    # The errors are accumulated and asserted once per rotor rather than
                    # element by element, so a broken engine reports a handful of failures
                    # instead of tens of thousands.
                    for m′ in -ℓ:ℓ, m in -ℓ:ℓ
                        worst = max(worst, abs(𝔇ˡ[m′, m] - ref[ℓ+1][m′+ℓ+1, m+ℓ+1]))
                        if T === Float64
                            worstᵇ = max(worstᵇ, abs(𝔇ˡ[m′, m] - D_oracle(Rotor(R), ℓ, m′, m)))
                        end
                    end
                end
                @test worst ≤ atol
                T === Float64 && @test worstᵇ ≤ atolᵇ
            end
        end
    end
end


@testitem "dCalculator vs closed form" setup=[HalfIntegerOracle, Utilities] begin
    import SphericalFunctions: dCalculator, recurrence!
    import .HalfIntegerOracle: d_oracle
    import .Utilities: βrange
    using Quaternionic: Rotor, from_euler_angles
    using DoubleFloats: Double64
    using Random

    rng = Random.Xoshiro(2345)

    # Reference, category 1 (a closed-form formula): Varshalovich Eq. 4.3.1(2), whose index
    # order is this package's, as transcribed in the `HalfIntegerOracle` setup module.  It
    # is evaluated at four times the working precision and rounded to `T`, so its own error
    # is at most half an ulp; evaluated in Float64 the alternating sum is already wrong in
    # the third-from-last digit by ℓ = 16.
    function dref(::Type{T}, β, ℓₘₐₓ) where {T}
        blocks = setprecision(BigFloat, 4 * precision(T) + 64) do
            βᵣ = BigFloat(β)
            [BigFloat[d_oracle(ℓ, m′, m, βᵣ) for m′ in -ℓ:ℓ, m in -ℓ:ℓ] for ℓ in 0:ℓₘₐₓ]
        end
        [T.(b) for b in blocks]  # `blocks[ℓ+1][m′+ℓ+1, m+ℓ+1]`
    end

    @testset "$T" for T in (Float64, Double64)
        # Measured worst error at ℓ ≤ 8: 3.8 eps (Float64) and 2.2 eps (Double64).
        atol = 20 * eps(T)
        for ℓₘₐₓ in (0, 1, 2, 4, 8)
            βs = βrange(rng, T, 6)
            calc = dCalculator(first(βs), ℓₘₐₓ)
            for β in βs
                eⁱᵝ = cis(β)
                ref = dref(T, β, ℓₘₐₓ)
                # The same β supplied as an angle, as a phase, and as a rotor whose α and γ
                # must be ignored
                R = from_euler_angles(T(2π * rand(rng)), β, T(2π * rand(rng)))
                for input in (β, eⁱᵝ, R)
                    worst = zero(T)
                    for ℓ in 0:ℓₘₐₓ
                        dˡ = if ℓ == 0
                            recurrence!(calc, input, ℓ)
                        else
                            recurrence!(calc, ℓ)
                        end
                        @test axes(dˡ) == (-ℓ:ℓ, -ℓ:ℓ)
                        @test eltype(dˡ) === T  # d is real
                        @test parent(dˡ) isa AbstractMatrix{T}
                        for m′ in -ℓ:ℓ, m in -ℓ:ℓ
                            worst = max(worst, abs(dˡ[m′, m] - ref[ℓ+1][m′+ℓ+1, m+ℓ+1]))
                        end
                    end
                    @test worst ≤ atol
                end
            end
        end
    end
end


@testitem "Wigner calculators vs explicit formulas" setup=[ExplicitWignerMatrices, Utilities] begin
    import SphericalFunctions: DCalculator, dCalculator, recurrence!
    import .Utilities: Rrange
    using Quaternionic: Rotor, 𝐢, 𝐣, 𝐤, to_euler_phases
    using Random

    rng = Random.Xoshiro(3456)
    ℓₘₐₓ = 3

    @testset "$T" for T in (Float64, BigFloat)
        # At ℓ ≤ 3 the closed forms lose nothing measurable: the worst errors measured here
        # are 4.7 eps for 𝔇, 3.0 eps for the d formula and 2.2 eps for the explicit d, in
        # both types.
        atol = 20 * eps(T)
        rotors = Rrange(rng, T, 6)
        calcD = DCalculator(first(rotors), ℓₘₐₓ)
        calcd = dCalculator(first(rotors), ℓₘₐₓ)
        for R in rotors
            eⁱᵅ, eⁱᵝ, eⁱᵞ = to_euler_phases(R)
            for ℓ in 0:ℓₘₐₓ
                𝔇ˡ = recurrence!(calcD, R, ℓ)
                dˡ = recurrence!(calcd, R, ℓ)
                for m′ in -ℓ:ℓ, m in -ℓ:ℓ
                    𝔇ᶠ = ExplicitWignerMatrices.D_formula(ℓ, m′, m, eⁱᵅ, eⁱᵝ, eⁱᵞ)
                    dᶠ = ExplicitWignerMatrices.d_formula(ℓ, m′, m, eⁱᵝ)
                    @test 𝔇ˡ[m′, m] ≈ 𝔇ᶠ atol=atol
                    @test dˡ[m′, m] ≈ dᶠ atol=atol
                    if ℓ ≤ 2  # the explicit expressions are coded only up to ℓ=2
                        dᵉ = ExplicitWignerMatrices.d_explicit(ℓ, m′, m, eⁱᵝ)
                        @test dˡ[m′, m] ≈ dᵉ atol=atol
                    end
                end
            end
        end
    end
end


@testitem "Wigner calculator block limits" setup=[RefusalChecks] begin
    import SphericalFunctions
    import SphericalFunctions: DCalculator, dCalculator, recurrence!
    using Quaternionic: Rotor
    using OffsetArrays: OffsetVector
    using Random

    rng = Random.Xoshiro(4567)
    ℓₘₐₓ = 5
    R = randn(rng, Rotor{Float64})
    R₂ = randn(rng, Rotor{Float64})

    # Full-range reference blocks for `R`, copied out for every ℓ
    function full_blocks(Ctor)
        calc = Ctor(R, ℓₘₐₓ)
        OffsetVector([copy(recurrence!(calc, ℓ)) for ℓ in 0:ℓₘₐₓ], 0:ℓₘₐₓ)
    end

    for (name, Ctor) in (("DCalculator", DCalculator), ("dCalculator", dCalculator))
        @testset "$name" begin
            full = full_blocks(Ctor)
            for m′ₘₐₓ in 0:ℓₘₐₓ, m′ₘᵢₙ in -ℓₘₐₓ:0, mₘₐₓ in (5, 3, 0), mₘᵢₙ in (-5, -2, 0)
                calc = Ctor(R, ℓₘₐₓ; m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
                @test SphericalFunctions.ℓₘₐₓ(calc) == ℓₘₐₓ
                @test SphericalFunctions.ℓₘᵢₙ(calc) == 0
                @test SphericalFunctions.m′ₘₐₓ(calc) == m′ₘₐₓ
                @test SphericalFunctions.m′ₘᵢₙ(calc) == m′ₘᵢₙ
                @test SphericalFunctions.mₘₐₓ(calc) == mₘₐₓ
                @test SphericalFunctions.mₘᵢₙ(calc) == mₘᵢₙ
                @test SphericalFunctions.Nᵣ(calc) == 1
                # `similar` reproduces the sizes and types, with its own storage
                twin = similar(calc)
                @test typeof(twin) === typeof(calc)
                @test (
                    SphericalFunctions.m′ₘₐₓ(twin), SphericalFunctions.m′ₘᵢₙ(twin),
                    SphericalFunctions.mₘₐₓ(twin), SphericalFunctions.mₘᵢₙ(twin)
                ) == (m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
                for ℓ in 0:ℓₘₐₓ
                    if ℓ == 0
                        blk = recurrence!(calc, R, ℓ)
                        blktwin = recurrence!(twin, R, ℓ)
                    else
                        blk = recurrence!(calc, ℓ)
                        blktwin = recurrence!(twin, ℓ)
                    end
                    m′r = max(-ℓ, m′ₘᵢₙ):min(ℓ, m′ₘₐₓ)
                    mr = max(-ℓ, mₘᵢₙ):min(ℓ, mₘₐₓ)
                    @test axes(blk) == (m′r, mr)
                    # Restricting the block never changes a value: the limited calculator
                    # runs exactly the same operations for the rows it keeps
                    @test all(blk[m′, m] == full[ℓ][m′, m] for m′ in m′r, m in mr)
                    @test blktwin == blk
                end
                # Storage is not shared between a calculator and its `similar`
                blktwin = recurrence!(twin, R₂, ℓₘₐₓ)
                blk = recurrence!(calc, ℓₘₐₓ)
                @test blktwin != blk
                @test all(
                    blk[m′, m] == full[ℓₘₐₓ][m′, m]
                    for m′ in axes(blk, 1), m in axes(blk, 2)
                )
            end

            # Invalid limits are rejected at construction
            limits(; kw...) = () -> Ctor(R, ℓₘₐₓ; kw...)
            small, large = "is too large for this index type", "is too large for ℓₘₐₓ"
            @test refuses(limits(m′ₘᵢₙ=1), ArgumentError, "m′ₘᵢₙ=1 $small")  # m′ₘᵢₙ > 0
            @test refuses(limits(mₘᵢₙ=1), ArgumentError, "mₘᵢₙ=1 $small")  # mₘᵢₙ > 0
            @test refuses(limits(m′ₘₐₓ=-1), ArgumentError, "m′ₘₐₓ=-1 is less than")
            @test refuses(limits(mₘₐₓ=-1), ArgumentError, "mₘₐₓ=-1 is less than")
            @test refuses(limits(m′ₘₐₓ=ℓₘₐₓ+1), ArgumentError, large)  # |limit| > ℓₘₐₓ
            @test refuses(limits(m′ₘₐₓ=ℓₘₐₓ+1, m′ₘᵢₙ=0), ArgumentError, large)
            @test refuses(limits(m′ₘₐₓ=0, m′ₘᵢₙ=-(ℓₘₐₓ+1)), ArgumentError, large)
            @test refuses(limits(mₘₐₓ=ℓₘₐₓ+1, mₘᵢₙ=0), ArgumentError, large)
            @test refuses(limits(mₘₐₓ=0, mₘᵢₙ=-(ℓₘₐₓ+1)), ArgumentError, large)
            @test refuses(limits(m′ₘₐₓ=1, m′ₘᵢₙ=2), ArgumentError, "is less than")  # max < min
            @test refuses(limits(mₘₐₓ=1, mₘᵢₙ=2), ArgumentError, "is less than")
            @test refuses(() -> Ctor(R, -1), ArgumentError, "must be non-negative")
            # The limits are indices of the kind of ℓₘₐₓ, and of type `Int`
            @test refuses(limits(m′ₘₐₓ=3//2), ArgumentError, "keyword argument `m′ₘₐₓ`")
            @test refuses(limits(mₘᵢₙ=Int8(-1)), ArgumentError, "narrower than `Int`")
            # The defaults themselves are valid at every ℓₘₐₓ, including 0
            @test SphericalFunctions.m′ₘₐₓ(Ctor(R, 0)) == 0

            # Each limit may also be spelled in ASCII, and the Unicode spelling wins where both
            # are given
            ascii = Ctor(R, ℓₘₐₓ; mp_max=3, mp_min=-2, m_max=4, m_min=-1)
            unicode = Ctor(R, ℓₘₐₓ; m′ₘₐₓ=3, m′ₘᵢₙ=-2, mₘₐₓ=4, mₘᵢₙ=-1)
            @test typeof(ascii) === typeof(unicode)
            @test copy(recurrence!(ascii, R, ℓₘₐₓ)) == copy(recurrence!(unicode, R, ℓₘₐₓ))
            @test SphericalFunctions.m′ₘₐₓ(Ctor(R, ℓₘₐₓ; mp_max=1, m′ₘₐₓ=2)) == 2
            # The default of each lower limit follows the upper limit, in either spelling
            @test SphericalFunctions.m′ₘᵢₙ(Ctor(R, ℓₘₐₓ; mp_max=2)) == -2
            @test SphericalFunctions.mₘᵢₙ(Ctor(R, ℓₘₐₓ; m_max=3)) == -3
        end
    end
end


@testitem "Wigner calculators batched rotors" setup=[RefusalChecks] begin
    import SphericalFunctions: DCalculator, dCalculator, recurrence!
    import SphericalFunctions: Nᵣ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ
    using Quaternionic: Rotor, to_euler_phases
    using OffsetArrays: OffsetVector
    using Random

    rng = Random.Xoshiro(5678)
    ℓₘₐₓ = 7
    N = 5
    rotors = randn(rng, Rotor{Float64}, N)
    # β ∈ [0, π] for each rotor, and the corresponding phases
    βs = [angle(to_euler_phases(R)[2]) for R in rotors]
    eⁱᵝs = cis.(βs)

    # Sequential ℓ on the batched calculator, comparing each rotor's block with the
    # single-rotor calculator's result for that rotor alone: identical operations, so
    # identical bits
    function check_batched(batched, single, data)
        @test Nᵣ(batched) == N
        @test Nᵣ(single) == 1
        for ℓ in 0:ℓₘₐₓ
            blk = if ℓ == 0
                recurrence!(batched, data, ℓ)
            else
                recurrence!(batched, ℓ)
            end
            m′r = max(-ℓ, m′ₘᵢₙ(batched)):min(ℓ, m′ₘₐₓ(batched))
            mr = max(-ℓ, mₘᵢₙ(batched)):min(ℓ, mₘₐₓ(batched))
            @test ndims(blk) == 3
            @test axes(blk) == (1:N, m′r, mr)
            for i in 1:N
                @test blk[i] == recurrence!(single, data[i], ℓ)
            end
        end
    end

    # Arbitrary ℓ order: backwards restarts, forwards advances, both reproduce the
    # sequential results exactly
    function check_order(batched, data)
        ref = OffsetVector([copy(recurrence!(batched, ℓ)) for ℓ in 0:ℓₘₐₓ], 0:ℓₘₐₓ)
        for ℓ in (4, 0, 7, 7, 2, 5, 1, 6, 3, 3, 0, 7)
            @test recurrence!(batched, ℓ) == ref[ℓ]
        end
        # Re-setting the rotor data and jumping straight to ℓ
        for ℓ in (7, 3, 0)
            @test recurrence!(batched, data, ℓ) == ref[ℓ]
        end
        # And from a fresh calculator of the same shape
        fresh = similar(batched)
        @test recurrence!(fresh, data, ℓₘₐₓ) == ref[ℓₘₐₓ]
    end

    @testset "DCalculator" begin
        check_batched(
            DCalculator(rotors, ℓₘₐₓ), DCalculator(first(rotors), ℓₘₐₓ), rotors
        )
        check_batched(
            DCalculator(rotors, ℓₘₐₓ; m′ₘₐₓ=2, mₘᵢₙ=-1),
            DCalculator(first(rotors), ℓₘₐₓ; m′ₘₐₓ=2, mₘᵢₙ=-1),
            rotors
        )
        check_order(DCalculator(rotors, ℓₘₐₓ), rotors)
        # Wrong number of rotors
        batched = DCalculator(rotors, ℓₘₐₓ)
        @test refuses(
            () -> recurrence!(batched, rotors[1:N-1], 0), DimensionMismatch,
            "This calculator handles Nᵣ=5 rotors, but got 4."
        )
        @test refuses(
            () -> recurrence!(batched, rotors[1], 0), DimensionMismatch,
            "This calculator handles Nᵣ=5 rotors, but a single rotor was given."
        )
        @test refuses(
            () -> recurrence!(DCalculator(first(rotors), ℓₘₐₓ), rotors, 0), DimensionMismatch,
            "This calculator handles Nᵣ=1 rotors, but got 5."
        )
    end

    @testset "dCalculator" begin
        for data in (βs, eⁱᵝs, rotors)
            check_batched(
                dCalculator(data, ℓₘₐₓ), dCalculator(first(data), ℓₘₐₓ), data
            )
            check_order(dCalculator(data, ℓₘₐₓ), data)
        end
        check_batched(
            dCalculator(βs, ℓₘₐₓ; m′ₘₐₓ=3, m′ₘᵢₙ=0, mₘₐₓ=4),
            dCalculator(first(βs), ℓₘₐₓ; m′ₘₐₓ=3, m′ₘᵢₙ=0, mₘₐₓ=4),
            βs
        )
        # Wrong number of angles
        batched = dCalculator(βs, ℓₘₐₓ)
        @test refuses(
            () -> recurrence!(batched, βs[1:2], 0), DimensionMismatch,
            "This calculator handles Nᵣ=5 rotors, but got 2."
        )
        @test refuses(
            () -> recurrence!(batched, βs[1], 0), DimensionMismatch,
            "This calculator handles Nᵣ=5 rotors, but a single rotor was given."
        )
        @test refuses(
            () -> recurrence!(dCalculator(first(βs), ℓₘₐₓ), βs, 0), DimensionMismatch,
            "This calculator handles Nᵣ=1 rotors, but got 5."
        )
    end
end


@testitem "Calculator block types are inferrable" begin
    import SphericalFunctions: DCalculator, dCalculator, sYlmCalculator
    import SphericalFunctions: recurrence!, isbatched, Nᵣ
    using Quaternionic: Rotor
    using Random

    # The payoff, measured the way a user's inner loop sees it: inside a function, where the
    # calculator's type is known, a step allocates nothing.  Measuring at top level instead
    # would report the boxing of a dynamically dispatched call and prove nothing.  This
    # covers the recurrence and the block together, since one call does both; a block
    # whose type were not concrete would show up here as the union split's allocation.
    blockallocs(c, ℓ) = (recurrence!(c, ℓ); @allocated recurrence!(c, ℓ))

    # Whether `recurrence!` returns a single block or a batch of them is a type parameter, not a
    # runtime test of `Nᵣ`, so the return type is concrete rather than a union of the two.
    # Guarding that here because nothing else would notice it silently regressing: the union
    # is split by the compiler, so the cost is one small allocation per call, not a failure.

    rng = Random.Xoshiro(2718)
    R = randn(rng, Rotor{Float64})
    rotors = randn(rng, Rotor{Float64}, 3)

    for (ℓₘₐₓ, ℓ) ∈ ((4, 3), (5//2, 3//2))
        @testset "ℓₘₐₓ = $ℓₘₐₓ" begin
            # A `Rotor` is acceptable rotor data for both calculators — `dCalculator`
            # takes the β Euler angle from it — so one set of inputs serves both here.
            for Ctor ∈ (DCalculator, dCalculator)
                single = Ctor(R, ℓₘₐₓ)
                @test !isbatched(single)
                @test isconcretetype(typeof(single))
                recurrence!(single, R, ℓ)
                @test isconcretetype(
                    Base.return_types(recurrence!, (typeof(single), typeof(ℓ)))[1]
                )
                @inferred recurrence!(single, ℓ)  # fails, without `@test`, if not inferred
                @test blockallocs(single, ℓ) == 0

                batched = Ctor(rotors, ℓₘₐₓ)
                @test isbatched(batched)
                @test isconcretetype(typeof(batched))
                recurrence!(batched, rotors, ℓ)
                @test isconcretetype(
                    Base.return_types(recurrence!, (typeof(batched), typeof(ℓ)))[1]
                )
                @inferred recurrence!(batched, ℓ)  # fails, without `@test`, if not inferred
                @test blockallocs(batched, ℓ) == 0
            end

            # ... and the same for the harmonics, where the spin-weight argument is a second
            # type parameter, so the block has four shapes to keep concrete rather than two
            s = ℓₘₐₓ isa Integer ? 2 : 3//2
            srange = ℓₘₐₓ isa Integer ? (-2:2) : (-3//2:3//2)
            for spec ∈ (s, srange)
                single = sYlmCalculator(R, ℓₘₐₓ, spec)
                @test !isbatched(single)
                recurrence!(single, R, ℓ)
                @test isconcretetype(
                    Base.return_types(recurrence!, (typeof(single), typeof(ℓ)))[1]
                )
                @inferred recurrence!(single, ℓ)  # fails, without `@test`, if not inferred
                @test blockallocs(single, ℓ) == 0

                batched = sYlmCalculator(rotors, ℓₘₐₓ, spec)
                @test isbatched(batched)
                recurrence!(batched, rotors, ℓ)
                @test isconcretetype(
                    Base.return_types(recurrence!, (typeof(batched), typeof(ℓ)))[1]
                )
                @inferred recurrence!(batched, ℓ)  # fails, without `@test`, if not inferred
                @test blockallocs(batched, ℓ) == 0

                # Slicing one spin weight out of a multi-spin block keeps that concreteness
                if spec isa AbstractRange
                    @test isconcretetype(typeof(recurrence!(single, ℓ)[s, :]))
                    @test isconcretetype(typeof(recurrence!(batched, ℓ)[:, s, :]))
                end
            end
        end
    end

    # `similar` must preserve the parameter, which it can only do by asserting it
    for c ∈ (DCalculator(R, 4), DCalculator(rotors, 4),
             sYlmCalculator(R, 4, 2), sYlmCalculator(rotors, 4, 2),
             sYlmCalculator(R, 4, -2:2), sYlmCalculator(rotors, 4, -2:2))
        @test typeof(@inferred similar(c)) === typeof(c)
        @test Nᵣ(similar(c)) == Nᵣ(c)
        @test isbatched(similar(c)) == isbatched(c)
    end
end

@testitem "Wigner calculators range errors" setup=[RefusalChecks] begin
    import SphericalFunctions
    import SphericalFunctions: DCalculator, dCalculator, recurrence!, WignerMatrix
    using Quaternionic: Rotor
    using Random

    rng = Random.Xoshiro(6789)
    ℓₘₐₓ = 4
    R = randn(rng, Rotor{Float64})

    for (name, Ctor) in (("DCalculator", DCalculator), ("dCalculator", dCalculator))
        @testset "$name" begin
            calc = Ctor(R, ℓₘₐₓ)
            # A calculator is not indexed at all; `recurrence!` is the only way in
            @test_throws MethodError calc[2]
            # Nothing has been computed yet, and `ℓ` says so
            @test SphericalFunctions.ℓ(calc) == SphericalFunctions.ℓₘᵢₙ(calc) - 1
            blk = recurrence!(calc, R, 2)
            @test SphericalFunctions.ℓ(calc) == 2
            @test size(blk) == (5, 5)
            # ℓ out of range for this calculator, or not an `Int`; the failed call leaves
            # the current block in place
            @test refuses(() -> recurrence!(calc, ℓₘₐₓ + 1), ArgumentError, "out of bounds")
            @test refuses(() -> recurrence!(calc, -1), ArgumentError, "out of bounds")
            kind = "The indices of this `$name` are integers of type `Int`, like 3; got ℓ ="
            @test refuses(() -> recurrence!(calc, 2.0), ArgumentError, "$kind 2.0::Float64")
            @test refuses(() -> recurrence!(calc, 5//2), ArgumentError, "$kind 5//2::")
            @test refuses(() -> recurrence!(calc, R, Int8(2)), ArgumentError, "$kind 2::")
            @test SphericalFunctions.ℓ(calc) == 2
            reference = copy(blk)
            # Recomputing from NaN-filled storage reproduces the result exactly, so no
            # uninitialized element is ever read
            fill!(calc, NaN)
            blk = recurrence!(calc, R, 2)
            @test blk == reference
            @test !any(isnan, blk)
            @test refuses(() -> recurrence!(calc, R, ℓₘₐₓ + 1), ArgumentError, "out of bounds")
        end
    end

    # A DCalculator needs the full rotor, not just β, and says what to use instead.  (Matching
    # "Rotor" alone would not do: a bare `MethodError` lists candidates that mention it.)
    calcD = DCalculator(R, ℓₘₐₓ)
    for β ∈ (0.3, cis(0.3), [0.3])
        @test refuses(
            () -> recurrence!(calcD, β, 0), ArgumentError,
            "use a dCalculator if only β is available"
        )
    end
    # ... whereas a dCalculator accepts any of the three forms
    calcd = dCalculator(R, ℓₘₐₓ)
    for input in (0.3, cis(0.3), R)
        # Deliberately *not* an `AbstractMatrix`: see the note on `AbstractBlock`.
        blk = recurrence!(calcd, input, 1)
        @test blk isa WignerMatrix && eltype(blk) === Float64
    end
end


@testitem "D and d convenience functions" setup=[RefusalChecks] begin
    import SphericalFunctions: DCalculator, dCalculator, recurrence!, D, d
    using Quaternionic: Rotor, to_euler_phases
    import SphericalFunctions: WignerSeries, WignerMatrix
    using Random

    rng = Random.Xoshiro(7890)
    ℓₘₐₓ = 5
    R = randn(rng, Rotor{Float64})
    R₂ = randn(rng, Rotor{Float64})
    # The phase and rotor forms differ from the angle form only through the rounding of β
    # itself, measured at 1.1 eps; Float32 results differ from Float64 by at most 1.0
    # eps(Float32)
    atol = 8eps(Float64)
    atol32 = 4eps(Float32)
    limits = (m′ₘₐₓ=2, m′ₘᵢₙ=-1, mₘₐₓ=3, mₘᵢₙ=0)

    @testset "D" begin
        𝔇 = D(R, ℓₘₐₓ)
        @test 𝔇 isa WignerSeries
        @test axes(𝔇) == (0:ℓₘₐₓ,)
        calc = DCalculator(R, ℓₘₐₓ)
        for ℓ in 0:ℓₘₐₓ
            @test 𝔇[ℓ] isa WignerMatrix && eltype(𝔇[ℓ]) === ComplexF64
            @test axes(𝔇[ℓ]) == (-ℓ:ℓ, -ℓ:ℓ)
            @test parent(𝔇[ℓ]) isa Matrix{ComplexF64}  # a copy, not a view
            @test 𝔇[ℓ] == recurrence!(calc, R, ℓ)
        end
        # Each block is an independent copy: computing more matrices, for a different rotor,
        # changes nothing
        snapshot = deepcopy(𝔇)
        𝔇₂ = D(R₂, ℓₘₐₓ)
        blk = recurrence!(calc, R₂, ℓₘₐₓ)
        @test all(𝔇[ℓ] == snapshot[ℓ] for ℓ in 0:ℓₘₐₓ)
        @test 𝔇₂[ℓₘₐₓ] == blk
        @test 𝔇₂[ℓₘₐₓ] != 𝔇[ℓₘₐₓ]
        # Block limits shrink the blocks without changing any value
        𝔇ₗ = D(R, ℓₘₐₓ; limits...)
        @test axes(𝔇ₗ) == (0:ℓₘₐₓ,)
        for ℓ in 0:ℓₘₐₓ
            m′r = max(-ℓ, limits.m′ₘᵢₙ):min(ℓ, limits.m′ₘₐₓ)
            mr = max(-ℓ, limits.mₘᵢₙ):min(ℓ, limits.mₘₐₓ)
            @test axes(𝔇ₗ[ℓ]) == (m′r, mr)
            @test all(𝔇ₗ[ℓ][m′, m] == 𝔇[ℓ][m′, m] for m′ in m′r, m in mr)
        end
        # ... in either spelling of the keywords
        𝔇ₐ = D(R, ℓₘₐₓ; mp_max=2, mp_min=-1, m_max=3, m_min=0)
        @test all(𝔇ₐ[ℓ] == 𝔇ₗ[ℓ] for ℓ in 0:ℓₘₐₓ)
        # `D` takes one rotor; a vector of them is what `DCalculator` is for
        @test refuses(() -> D([R, R₂], ℓₘₐₓ), ArgumentError, "DCalculator(R⃗, ℓₘₐₓ)")
        @test refuses(() -> D([R], ℓₘₐₓ; m′ₘₐₓ=1), ArgumentError, "takes a single rotor")
        # The element type follows the rotor's type
        @test eltype(D(Rotor{Float32}(R), 2)[2]) === ComplexF32
        @test eltype(D(R, 2)[2]) === ComplexF64
        @test eltype(D(Rotor{BigFloat}(R), 2)[2]) === Complex{BigFloat}
        # ... and lower precision only rounds the result
        𝔇₃₂ = D(Rotor{Float32}(R), 3)
        for ℓ in 0:3
            @test maximum(abs, 𝔇₃₂[ℓ] .- 𝔇[ℓ]) ≤ atol32
        end
    end

    @testset "d" begin
        eⁱᵝ = to_euler_phases(R)[2]
        β = angle(eⁱᵝ)  # ∈ [0, π]
        dβ = d(β, ℓₘₐₓ)
        @test dβ isa WignerSeries
        @test axes(dβ) == (0:ℓₘₐₓ,)
        calc = dCalculator(β, ℓₘₐₓ)
        for ℓ in 0:ℓₘₐₓ
            @test dβ[ℓ] isa WignerMatrix && eltype(dβ[ℓ]) === Float64
            @test axes(dβ[ℓ]) == (-ℓ:ℓ, -ℓ:ℓ)
            @test parent(dβ[ℓ]) isa Matrix{Float64}  # a copy, not a view
            @test dβ[ℓ] == recurrence!(calc, β, ℓ)
        end
        # Independent copies
        snapshot = deepcopy(dβ)
        d₂ = d(β / 3, ℓₘₐₓ)
        blk = recurrence!(calc, β / 3, ℓₘₐₓ)
        @test all(dβ[ℓ] == snapshot[ℓ] for ℓ in 0:ℓₘₐₓ)
        @test d₂[ℓₘₐₓ] == blk
        @test d₂[ℓₘₐₓ] != dβ[ℓₘₐₓ]
        # The phase and rotor forms give the same matrices (up to the rounding of β itself)
        for input in (eⁱᵝ, R)
            dᵢ = d(input, ℓₘₐₓ)
            @test dᵢ isa WignerSeries
            @test axes(dᵢ) == (0:ℓₘₐₓ,)
            for ℓ in 0:ℓₘₐₓ
                @test axes(dᵢ[ℓ]) == (-ℓ:ℓ, -ℓ:ℓ)
                @test maximum(abs, dᵢ[ℓ] .- dβ[ℓ]) ≤ atol
            end
        end
        # `d` takes one angle, phase or rotor; a vector of them is what `dCalculator` is for
        for input in ([β, β / 3], cis.([β]), [R, R₂])
            @test refuses(() -> d(input, ℓₘₐₓ), ArgumentError, "dCalculator(β⃗, ℓₘₐₓ)")
        end
        # Block limits, in either spelling
        @test all(
            d(β, ℓₘₐₓ; mp_max=2, mp_min=-1, m_max=3, m_min=0)[ℓ] == d(β, ℓₘₐₓ; limits...)[ℓ]
            for ℓ in 0:ℓₘₐₓ
        )
        for input in (β, eⁱᵝ, R)
            dₗ = d(input, ℓₘₐₓ; limits...)
            @test axes(dₗ) == (0:ℓₘₐₓ,)
            for ℓ in 0:ℓₘₐₓ
                m′r = max(-ℓ, limits.m′ₘᵢₙ):min(ℓ, limits.m′ₘₐₓ)
                mr = max(-ℓ, limits.mₘᵢₙ):min(ℓ, limits.mₘₐₓ)
                @test axes(dₗ[ℓ]) == (m′r, mr)
                if input === β  # the same computation, restricted
                    @test all(dₗ[ℓ][m′, m] == dβ[ℓ][m′, m] for m′ in m′r, m in mr)
                else
                    @test all(isapprox(dₗ[ℓ][m′, m], dβ[ℓ][m′, m]; atol) for m′ in m′r, m in mr)
                end
            end
        end
        # The element type follows the input
        @test eltype(d(Float32(β), 2)[2]) === Float32
        @test eltype(d(cis(Float32(β)), 2)[2]) === Float32
        @test eltype(d(Rotor{Float32}(R), 2)[2]) === Float32
        @test eltype(d(β, 2)[2]) === Float64
        @test eltype(d(eⁱᵝ, 2)[2]) === Float64
        @test eltype(d(R, 2)[2]) === Float64
        @test eltype(d(BigFloat(β), 2)[2]) === BigFloat
        @test eltype(d(Rotor{BigFloat}(R), 2)[2]) === BigFloat
        @test eltype(d(1, 2)[2]) === Float64  # an integer angle is promoted to Float64
        @test maximum(abs, d(Float32(β), 3)[3] .- dβ[3]) ≤ atol32
    end
end


@testitem "Wigner calculators vs independent references" setup=[HalfIntegerOracle] begin
    import SphericalFunctions
    import SphericalFunctions: DCalculator, dCalculator, recurrence!, D
    import .HalfIntegerOracle: d_oracle, D_oracle
    using Quaternionic: Rotor, Quaternion, components
    using LinearAlgebra: I, opnorm
    using Random

    # Every element of 𝔇ˡ, for several generic rotors, at every ℓ up to 6, is compared here
    # with two closed forms that owe nothing to the package, and every element of dˡ with
    # one of them, and both are checked against identities that the recurrence cannot
    # satisfy by accident.

    rng = Random.Xoshiro(8901)
    ℓₘₐₓ = 6
    T = Float64
    rotors = randn(rng, Rotor{T}, 5)

    # High-precision Euler angles of R, taken from the quaternion components (see the
    # "DCalculator vs closed forms" item for the derivation and for why
    # `to_euler_angles` is not used here).
    function euler_angles(R)
        w, x, y, z = BigFloat.(components(Quaternion(R)))
        ϕₛ, ϕₐ = angle(Complex(w, z)), angle(Complex(y, x))
        (ϕₛ - ϕₐ, 2 * atan(abs(Complex(y, x)), abs(Complex(w, z))), ϕₛ + ϕₐ)
    end

    # Reference 1 (category 1: closed forms).  𝔇ˡ_{m′m} = e^{-i m′ α} dˡ_{m′m}(β) e^{-i m γ}
    # from the conventions pages with Varshalovich Eq. 4.3.1(2) for d, and — for 𝔇 only —
    # the quaternionic form of Boyle (2016), which never decomposes R into Euler angles.
    # Both come from the `HalfIntegerOracle` setup module.  The first is evaluated at four
    # times Float64 precision and rounded, so it contributes at most half an ulp of its
    # own; the second is only Float64-accurate (its normalization is a ratio of `Int`
    # factorials), hence its looser tolerance.
    #
    # Measured worst errors at ℓ ≤ 6: 4.6 eps for 𝔇 and 2.6 eps for d against the
    # Varshalovich references, 10.4 eps for 𝔇 against Boyle (2016).
    atol = 20 * eps(T)
    atolᵇ = 64 * eps(T)
    atolᵘ = 32 * eps(T)  # the identities below are operator norms of (2ℓ+1)-square matrices
    atolʳ = 64 * eps(T)  # ... and the representation property multiplies two of them

    calcD = DCalculator(first(rotors), ℓₘₐₓ)
    calcd = dCalculator(first(rotors), ℓₘₐₓ)
    for R in rotors
        αᵣ, βᵣ, γᵣ = setprecision(BigFloat, 4 * precision(T) + 64) do
            euler_angles(R)
        end
        worstD = zero(T)
        worstd = zero(T)
        worstB = zero(T)
        worstF = zero(T)
        for ℓ in 0:ℓₘₐₓ
            𝔇ˡ, dˡ = if ℓ == 0
                recurrence!(calcD, R, ℓ), recurrence!(calcd, R, ℓ)
            else
                recurrence!(calcD, ℓ), recurrence!(calcd, ℓ)
            end
            @test size(𝔇ˡ) == size(dˡ) == (2ℓ + 1, 2ℓ + 1)
            setprecision(BigFloat, 4 * precision(T) + 64) do
                for m′ in -ℓ:ℓ, m in -ℓ:ℓ
                    dᵣ = d_oracle(ℓ, m′, m, βᵣ)
                    𝔇ᵣ = cis(-m′ * αᵣ) * dᵣ * cis(-m * γᵣ)
                    worstD = max(worstD, T(abs(Complex{BigFloat}(𝔇ˡ[m′, m]) - 𝔇ᵣ)))
                    worstd = max(worstd, T(abs(BigFloat(dˡ[m′, m]) - dᵣ)))
                    worstB = max(worstB, abs(𝔇ˡ[m′, m] - D_oracle(R, ℓ, m′, m)))
                    # Reference 3 (category 3: a metamorphic identity relating the two
                    # calculators, which run the same H recurrence but apply different
                    # phases to it).  The phases here come from the high-precision angles,
                    # so only the package's own error shows up.
                    𝔇ᶠ = Complex{T}(cis(-m′ * αᵣ)) * dˡ[m′, m] * Complex{T}(cis(-m * γᵣ))
                    worstF = max(worstF, abs(𝔇ˡ[m′, m] - 𝔇ᶠ))
                end
            end
            # Reference 2 (category 3: metamorphic identities).  𝔇ˡ is unitary and dˡ is
            # real orthogonal, measured to 9.0 and 8.4 eps respectively at ℓ ≤ 6.
            M = parent(𝔇ˡ)
            Md = parent(dˡ)
            @test opnorm(M * M' - I) ≤ atolᵘ
            @test opnorm(M' * M - I) ≤ atolᵘ
            @test opnorm(Md * Md' - I) ≤ atolᵘ
            @test opnorm(Md' * Md - I) ≤ atolᵘ
        end
        @test worstD ≤ atol
        @test worstd ≤ atol
        @test worstB ≤ atolᵇ
        @test worstF ≤ atol  # measured 4.5 eps
    end

    # The representation property 𝔇ˡ(R₁ R₂) = 𝔇ˡ(R₁) 𝔇ˡ(R₂) over every ordered pair of the
    # same rotors (category 3).  An operator norm again, and one that accumulates a whole
    # matrix product: measured 17.7 eps at ℓ ≤ 6.
    𝔇s = [D(R, ℓₘₐₓ) for R in rotors]
    for (i₁, R₁) in enumerate(rotors), (i₂, R₂) in enumerate(rotors)
        𝔇₁₂ = D(R₁ * R₂, ℓₘₐₓ)
        worst = zero(T)
        for ℓ in 0:ℓₘₐₓ
            worst = max(worst, opnorm(parent(𝔇s[i₁][ℓ]) * parent(𝔇s[i₂][ℓ]) - parent(𝔇₁₂[ℓ])))
        end
        @test worst ≤ atolʳ
    end
end


@testitem "Wigner calculators refuse integer indices of other types" setup=[RefusalChecks] begin
    import SphericalFunctions: DCalculator, dCalculator, recurrence!, D, d
    using Quaternionic: Rotor
    using Random

    # Index arithmetic is not closed under the other integer types — ℓ² overflows a narrow
    # one, and -m wraps around in an unsigned one — so an index of any integer type but `Int`
    # is refused, whether it is ℓₘₐₓ or a keyword limit, with a sentence saying why.
    rng = Random.Xoshiro(20260917)
    R = randn(rng, Rotor{Float64})
    for (IT, sentence) ∈ (
        (Int8, "narrower than `Int`"), (Int16, "narrower than `Int`"),
        (Int32, "narrower than `Int`"), (UInt, "is unsigned"), (Int128, "wider than `Int`"),
        (BigInt, "wider than `Int`"), (Bool, "A `Bool` is not an index"),
    )
        n = IT === Bool ? true : IT(3)
        @test refuses(() -> D(R, n), ArgumentError, sentence)
        @test refuses(() -> d(R, n), ArgumentError, sentence)
        for Ctor in (DCalculator, dCalculator)
            @test refuses(() -> Ctor(R, n), ArgumentError, sentence)
            @test refuses(() -> Ctor(R, 4; m′ₘₐₓ=n), ArgumentError, sentence)
            @test refuses(() -> Ctor(R, 4; mp_max=n), ArgumentError, sentence)
        end
    end
    # A `Rational` that is not a half-odd-integer of `Int`s is refused as well
    @test refuses(() -> D(R, 3//1), ArgumentError, "is a whole number")
    @test refuses(() -> D(R, big(7)//2), ArgumentError, "is not `Rational{Int}`")
    # ... and so is an ℓ of another type given to `recurrence!`, with the message that names
    # the calculator
    for (IT, sentence) ∈ (
        (Int8, "narrower than `Int`"), (UInt, "is unsigned"), (Int128, "wider than `Int`"),
        (BigInt, "wider than `Int`"), (Bool, "A `Bool` is not an index"),
    )
        n = IT === Bool ? true : IT(2)
        for (name, Ctor) ∈ (("DCalculator", DCalculator), ("dCalculator", dCalculator))
            calc = Ctor(R, 3)
            @test refuses(() -> recurrence!(calc, n), ArgumentError, sentence)
            @test refuses(
                () -> recurrence!(calc, n), ArgumentError,
                "The indices of this `$name` are integers of type `Int`, like 3; got ℓ = "
            )
        end
    end
end

@testitem "Calculators: a vector of rotor data is a batch, however short" begin
    import SphericalFunctions: DCalculator, dCalculator, sYlmCalculator, sλlmCalculator,
        recurrence!, isbatched, sYlm, WignerMatrix, WignerMatrixBatch, DegreeBlock,
        DegreeBlockBatch
    using Quaternionic: from_euler_angles
    using Test: @inferred

    # Whether the blocks have a rotor index is decided by whether the data is a vector, not by
    # its length, so that a loop written for a batch works for a batch of one, and `isbatched`
    # of an `sYlm` result agrees with the type of its blocks.
    R = from_euler_angles(0.3, 0.7, 1.1)
    @test !isbatched(DCalculator(R, 2)) && isbatched(DCalculator([R], 2))
    @test recurrence!(DCalculator(R, 2), 2) isa WignerMatrix
    𝔇ˡ = recurrence!(DCalculator([R], 2), 2)
    @test 𝔇ˡ isa WignerMatrixBatch
    @test 𝔇ˡ[1, 0, 0] == recurrence!(DCalculator(R, 2), 2)[0, 0]
    @test isbatched(dCalculator([0.7], 2)) && !isbatched(dCalculator(0.7, 2))
    @test isbatched(sYlmCalculator([R], 2, 0)) && !isbatched(sYlmCalculator(R, 2, 0))
    @test recurrence!(sYlmCalculator([R], 2, 0), 2) isa DegreeBlockBatch
    @test isbatched(sλlmCalculator([0.7], 2, 0))

    # ... and the values `sYlm` returns say the same as their blocks
    Y = sYlm([R], 2, 0)
    @test isbatched(Y) && Y[2] isa DegreeBlockBatch
    @test !isbatched(sYlm(R, 2, 0)) && sYlm(R, 2, 0)[2] isa DegreeBlock
    @test isbatched(sYlm([R], 2, -1:1)) && !isbatched(sYlm(R, 2, -1:1))

    # ... and a batched calculator says so when shown, since its blocks are indexed differently
    @test occursin("batched", sprint(show, DCalculator([R], 2)))
    @test !occursin("batched", sprint(show, DCalculator(R, 2)))
    @test occursin("batched", sprint(show, sYlmCalculator([R], 2, 0)))
    @test !occursin("batched", sprint(show, sYlmCalculator(R, 2, 0)))

    # The batchedness is known from the type of the argument, so construction is inferrable,
    # and so is the block
    @inferred DCalculator([R], 2)
    @inferred DCalculator(R, 2)
    @inferred recurrence!(DCalculator([R], 2), 2)
    @inferred DCalculator(R, 7//2; m′ₘₐₓ=1//2)
    @inferred DCalculator(R, 4; mp_max=2)
    @inferred dCalculator([0.7], 2)
    @inferred dCalculator(R, 7//2)
    @inferred sYlmCalculator([R], 2, -1:1)
    @inferred sYlmCalculator(R, 7//2, 1//2)
    @inferred sλlmCalculator([0.7], 2, 0)
    @inferred recurrence!(sYlmCalculator([R], 2, -1:1), 2)
end


@testitem "Wigner calculators: each block element is the wedge element `wedge_value` reads" begin
    import SphericalFunctions
    import SphericalFunctions: DCalculator, dCalculator, recurrence!, wedge_value, zpower, ϵ,
        HalfOddInteger, Nᵣ, isbatched, ℓₘᵢₙ, ℓₘₐₓ
    using Quaternionic: Rotor
    import Random

    # `materialize!` reads the wedge in runs along its rows, resolving the symmetries of H
    # once for each run rather than once for each element, as `wedge_source` does.  Each
    # element of every block is compared here, bit for bit, with the value built from
    # `wedge_value`, which applies the symmetries one element at a time, and from the same
    # power tables: the ϵ signs of d, and for 𝔇 the phase e^{-i(m′α+mγ)} = conj(z₊^(m′+m)
    # z₋^(m′-m)).  The restrictions include rows or columns narrower than the other range,
    # so that the wedge is narrowed too, and asymmetric ranges; the calculators are single
    # and batched.  The rotors and angles include some exactly at both poles, β = 0 and β =
    # π.
    rng = Random.Xoshiro(20260924)
    R⃗ = [randn(rng, Rotor{Float64}, 5); Rotor(1.0, 0.0, 0.0, 0.0); Rotor(0.0, 0.6, 0.8, 0.0)]
    β⃗ = [0.0, 0.4, 1.9, π, 2.7]
    function expected(calc, H, iᵣ, m′, m)
        RT = SphericalFunctions.floattype(calc)
        dᵐ′ᵐ = convert(RT, ϵ(m′) * ϵ(-m)) * wedge_value(H, iᵣ, m′, m)
        if eltype(calc.Wˡ) <: Complex
            dᵐ′ᵐ * conj(zpower(calc.engine.Z₊, iᵣ, m′ + m) * zpower(calc.engine.Z₋, iᵣ, m′ - m))
        else
            dᵐ′ᵐ
        end
    end
    integer_limits = [
        (;), (m′ₘₐₓ=2, m′ₘᵢₙ=-1, mₘₐₓ=3, mₘᵢₙ=-3), (mₘₐₓ=2, mₘᵢₙ=-2), (m′ₘₐₓ=2,),
        (mₘₐₓ=3, mₘᵢₙ=0, m′ₘₐₓ=7, m′ₘᵢₙ=-5), (mₘₐₓ=0, mₘᵢₙ=0), (m′ₘₐₓ=0, m′ₘᵢₙ=0),
        (m′ₘₐₓ=9, m′ₘᵢₙ=0, mₘₐₓ=1, mₘᵢₙ=-4), (mₘₐₓ=9, mₘᵢₙ=-1, m′ₘₐₓ=1, m′ₘᵢₙ=-2),
    ]
    half_limits = [
        (;), (m′ₘₐₓ=5//2, m′ₘᵢₙ=-3//2), (mₘₐₓ=3//2, mₘᵢₙ=-1//2), (mₘₐₓ=1//2,),
        (m′ₘₐₓ=17//2, m′ₘᵢₙ=-1//2, mₘₐₓ=3//2, mₘᵢₙ=-5//2), (m′ₘₐₓ=1//2, m′ₘᵢₙ=-1//2),
    ]
    count = Ref(0)
    for (ℓmax, limits) ∈ ((9, integer_limits), (17//2, half_limits)), lim ∈ limits,
            data ∈ (R⃗, R⃗[1], β⃗, β⃗[2])
        Ctor = eltype(data) <: Rotor ? DCalculator : dCalculator
        calc = Ctor(data, ℓmax; lim...)
        for ℓ ∈ ℓₘᵢₙ(calc):ℓₘₐₓ(calc)
            blk = recurrence!(calc, ℓ)
            H = calc.engine.H.Hˡ
            good = true
            for m′ ∈ axes(blk, isbatched(calc) ? 2 : 1), m ∈ axes(blk, isbatched(calc) ? 3 : 2),
                    iᵣ ∈ 1:Nᵣ(calc)
                value = isbatched(calc) ? blk[iᵣ, m′, m] : blk[m′, m]
                good &= isequal(value, expected(calc, H, iᵣ, m′, m))
                count[] += 1
            end
            @test good
        end
    end
    @test count[] > 50_000  # every element of every block of every case was compared
end


@testitem "Wigner calculators: restricting m narrows the wedge as restricting m′ does" begin
    import SphericalFunctions
    import SphericalFunctions: DCalculator, dCalculator, recurrence!, maxm′ₘₐₓ
    using Quaternionic: Rotor
    import Random

    # The wedge holds the rows |m′| ≤ W for every m, and an element of the block whose |m′|
    # exceeds W is read from its transpose, so W need only be the narrower of the widest m′
    # and the widest m.  A calculator restricted to a few columns then runs the recurrence of
    # one restricted to as many rows, rather than the whole of it, and its blocks are the
    # restriction of the full ones, bit for bit.
    rng = Random.Xoshiro(3)
    R⃗ = randn(rng, Rotor{Float64}, 3)
    for (ℓmax, narrow, wide) ∈ ((24, 2, 24), (47//2, 3//2, 47//2))
        full = DCalculator(R⃗, ℓmax)
        by_columns = DCalculator(R⃗, ℓmax; mₘₐₓ=narrow)
        by_rows = DCalculator(R⃗, ℓmax; m′ₘₐₓ=narrow)
        @test maxm′ₘₐₓ(by_columns.engine.H.Hˡ) == maxm′ₘₐₓ(by_rows.engine.H.Hˡ) == narrow
        @test maxm′ₘₐₓ(full.engine.H.Hˡ) == wide
        @test length(parent(by_columns.engine.H.Hˡ)) == length(parent(by_rows.engine.H.Hˡ))
        @test length(parent(by_columns.engine.H.Hˡ)) < length(parent(full.engine.H.Hˡ)) ÷ 5
        # `similar` builds the same narrow wedge
        @test maxm′ₘₐₓ(similar(by_columns).engine.H.Hˡ) == narrow
        # An asymmetric range is as wide as its larger end, and the narrower of the two
        # ranges decides; a lower limit alone narrows nothing
        @test maxm′ₘₐₓ(DCalculator(R⃗, ℓmax; mₘₐₓ=narrow, mₘᵢₙ=-narrow - 2).engine.H.Hˡ) == narrow + 2
        @test maxm′ₘₐₓ(DCalculator(R⃗, ℓmax; mₘₐₓ=narrow + 2, m′ₘₐₓ=narrow).engine.H.Hˡ) == narrow
        @test maxm′ₘₐₓ(DCalculator(R⃗, ℓmax; mₘᵢₙ=-narrow).engine.H.Hˡ) == wide
        for (ℓ, blk) ∈ by_columns
            ref = recurrence!(full, ℓ)
            @test all(
                isequal(blk[iᵣ, m′, m], ref[iᵣ, m′, m])
                for iᵣ ∈ 1:3, m′ ∈ axes(blk, 2), m ∈ axes(blk, 3)
            )
        end
        cβ = dCalculator([0.3, 2.1], ℓmax; mₘₐₓ=narrow, mₘᵢₙ=-narrow)
        fβ = dCalculator([0.3, 2.1], ℓmax)
        @test maxm′ₘₐₓ(cβ.engine.H.Hˡ) == narrow
        for (ℓ, blk) ∈ cβ
            ref = recurrence!(fβ, ℓ)
            @test all(
                isequal(blk[iᵣ, m′, m], ref[iᵣ, m′, m])
                for iᵣ ∈ 1:2, m′ ∈ axes(blk, 2), m ∈ axes(blk, 3)
            )
        end
    end
end


@testitem "Calculators: one rotor gives bit for bit the values it has in a batch" begin
    import SphericalFunctions: DCalculator, dCalculator, HCalculator, sYlmCalculator,
        sλlmCalculator, recurrence!, ℓₘₐₓ
    using Quaternionic: Rotor
    import Random

    # With one rotor there is nothing to vectorize, so the innermost loops of the recurrence
    # and of the assembly of the blocks write that case out as a single statement.  It is the
    # loop's own expression, so the values of a rotor are the same alone as in a batch.
    rng = Random.Xoshiro(11)
    R⃗ = randn(rng, Rotor{Float64}, 2)
    β⃗ = [0.4, 2.2]
    first_rotor(blk) = (a = Array(copy(blk)); a[1, ntuple(_ -> :, ndims(a) - 1)...])
    for ℓmax ∈ (13, 25//2)
        s = ℓmax isa Integer ? (-2:2) : (-3//2:1//2)
        narrow = ℓmax isa Integer ? 1 : 1//2
        for (batch, single) ∈ (
            (DCalculator(R⃗, ℓmax), DCalculator(R⃗[1], ℓmax)),
            (dCalculator(β⃗, ℓmax), dCalculator(β⃗[1], ℓmax)),
            (DCalculator(R⃗, ℓmax; mₘₐₓ=narrow), DCalculator(R⃗[1], ℓmax; mₘₐₓ=narrow)),
            (sYlmCalculator(R⃗, ℓmax, s), sYlmCalculator(R⃗[1], ℓmax, s)),
            (sYlmCalculator(R⃗, ℓmax, last(s)), sYlmCalculator(R⃗[1], ℓmax, last(s))),
            (sYlmCalculator(β⃗, ℓmax, s), sYlmCalculator(β⃗[1], ℓmax, s)),
            (sλlmCalculator(β⃗, ℓmax, s), sλlmCalculator(β⃗[1], ℓmax, s)),
        )
            @test all(
                isequal(first_rotor(b), Array(copy(a))) for ((_, b), (_, a)) ∈ zip(batch, single)
            )
        end
        batch, single = HCalculator(β⃗, ℓmax), HCalculator(β⃗[1], ℓmax)
        for ℓ ∈ (ℓₘₐₓ(batch) - 2):ℓₘₐₓ(batch)
            Hb, H₁ = recurrence!(batch, ℓ), recurrence!(single, ℓ)
            @test all(
                isequal(Hb[1, m′, m], H₁[1, m′, m])
                for m′ ∈ -Hb.m′ₘₐₓ:Hb.m′ₘₐₓ for m ∈ abs(m′):ℓ
            )
        end
    end
end
