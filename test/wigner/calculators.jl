# Tests for the v3 Wigner-matrix calculators: `DCalculator`, `dCalculator`, the
# `calc[ℓ]` views, the four block limits, batched rotors, and the convenience functions `D`
# and `d`.  The oracles are independent closed forms — Varshalovich Eq. 4.3.1(2) for `d`,
# the settled Euler factorization for `𝔇`, and the quaternionic form of Boyle (2016), all
# transcribed in the `HalfIntegerOracle` setup module — together with the explicit and
# formulaic matrices in `ExplicitWignerMatrices` and a set of metamorphic identities.

@testitem "DCalculator vs closed forms" setup=[HalfIntegerOracle, Utilities] begin
    import SphericalFunctions: DCalculator, recurrence!
    import .HalfIntegerOracle: d_oracle, D_oracle
    using Quaternionic: Rotor, Quaternion, components, 𝐢, 𝐣, 𝐤
    using DoubleFloats: Double64
    using Random

    Random.seed!(1234)  # `Rrange` draws from the default RNG

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
        # Measured worst errors at ℓ ≤ 8: 7.5 eps (Float64) and 2.1 eps (Double64) against
        # reference 1, and 15.6 eps against reference 2.
        atol = 32 * eps(T)
        atolᵇ = 64 * eps(T)
        for ℓₘₐₓ in (0, 1, 2, 4, 8)
            rotors = Rrange(T, 6)
            calc = DCalculator(first(rotors), ℓₘₐₓ)
            for R in rotors
                worst = zero(T)
                worstᵇ = zero(T)
                ref = Dref(T, R, ℓₘₐₓ)
                for ℓ in 0:ℓₘₐₓ
                    if ℓ == 0
                        recurrence!(calc, R, ℓ)  # set the rotor and compute ℓ=0
                    else
                        recurrence!(calc, ℓ)  # reuse the rotor data
                    end
                    𝔇ˡ = calc[ℓ]
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
    using Quaternionic: Rotor, from_euler_angles
    using DoubleFloats: Double64
    using Random

    Random.seed!(2345)  # `βrange` draws from the default RNG
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

    # `βrange(Double64, n)` overflows while building its step range, so for other types we
    # convert the Float64 values and add the exact endpoints of the type itself.
    βvalues(::Type{Float64}, n) = βrange(Float64, n)
    βvalues(::Type{T}, n) where {T} = T[T(0); T.(βrange(Float64, n)); T(π)]

    @testset "$T" for T in (Float64, Double64)
        # Measured worst error at ℓ ≤ 8: 3.0 eps (Float64) and 1.8 eps (Double64).
        atol = 20 * eps(T)
        for ℓₘₐₓ in (0, 1, 2, 4, 8)
            βs = βvalues(T, 6)
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
                        if ℓ == 0
                            recurrence!(calc, input, ℓ)
                        else
                            recurrence!(calc, ℓ)
                        end
                        dˡ = calc[ℓ]
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
    using Quaternionic: Rotor, 𝐢, 𝐣, 𝐤, to_euler_phases
    using Random

    Random.seed!(3456)  # `Rrange` draws from the default RNG
    ℓₘₐₓ = 3

    @testset "$T" for T in (Float64, BigFloat)
        # The naive closed-form sum loses digits in Float64; BigFloat has ~77 digits to spare
        atol = T === BigFloat ? big"1e-60" : 1e-13
        rotors = Rrange(T, 6)
        calcD = DCalculator(first(rotors), ℓₘₐₓ)
        calcd = dCalculator(first(rotors), ℓₘₐₓ)
        for R in rotors
            eⁱᵅ, eⁱᵝ, eⁱᵞ = to_euler_phases(R)
            for ℓ in 0:ℓₘₐₓ
                recurrence!(calcD, R, ℓ)
                recurrence!(calcd, R, ℓ)
                𝔇ˡ = calcD[ℓ]
                dˡ = calcd[ℓ]
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


@testitem "Wigner calculator block limits" begin
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
        OffsetVector([copy(recurrence!(calc, ℓ)[ℓ]) for ℓ in 0:ℓₘₐₓ], 0:ℓₘₐₓ)
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
                        recurrence!(calc, R, ℓ)
                        recurrence!(twin, R, ℓ)
                    else
                        recurrence!(calc, ℓ)
                        recurrence!(twin, ℓ)
                    end
                    blk = calc[ℓ]
                    m′r = max(-ℓ, m′ₘᵢₙ):min(ℓ, m′ₘₐₓ)
                    mr = max(-ℓ, mₘᵢₙ):min(ℓ, mₘₐₓ)
                    @test axes(blk) == (m′r, mr)
                    # Restricting the block never changes a value: the limited calculator
                    # runs exactly the same operations for the rows it keeps
                    @test all(blk[m′, m] == full[ℓ][m′, m] for m′ in m′r, m in mr)
                    @test twin[ℓ] == blk
                end
                # Storage is not shared between a calculator and its `similar`
                recurrence!(twin, R₂, ℓₘₐₓ)
                @test twin[ℓₘₐₓ] != calc[ℓₘₐₓ]
                @test all(calc[ℓₘₐₓ][m′, m] == full[ℓₘₐₓ][m′, m] for m′ in axes(calc[ℓₘₐₓ], 1), m in axes(calc[ℓₘₐₓ], 2))
            end

            # Invalid limits are rejected at construction
            @test_throws ErrorException Ctor(R, ℓₘₐₓ; m′ₘᵢₙ=1)  # m′ₘᵢₙ > 0
            @test_throws ErrorException Ctor(R, ℓₘₐₓ; mₘᵢₙ=1)  # mₘᵢₙ > 0
            @test_throws ErrorException Ctor(R, ℓₘₐₓ; m′ₘₐₓ=-1)  # m′ₘₐₓ < 0
            @test_throws ErrorException Ctor(R, ℓₘₐₓ; mₘₐₓ=-1)  # mₘₐₓ < 0
            @test_throws ErrorException Ctor(R, ℓₘₐₓ; m′ₘₐₓ=ℓₘₐₓ+1)  # |limit| > ℓₘₐₓ
            @test_throws ErrorException Ctor(R, ℓₘₐₓ; m′ₘₐₓ=ℓₘₐₓ+1, m′ₘᵢₙ=0)
            @test_throws ErrorException Ctor(R, ℓₘₐₓ; m′ₘₐₓ=0, m′ₘᵢₙ=-(ℓₘₐₓ+1))
            @test_throws ErrorException Ctor(R, ℓₘₐₓ; mₘₐₓ=ℓₘₐₓ+1, mₘᵢₙ=0)
            @test_throws ErrorException Ctor(R, ℓₘₐₓ; mₘₐₓ=0, mₘᵢₙ=-(ℓₘₐₓ+1))
            @test_throws ErrorException Ctor(R, ℓₘₐₓ; m′ₘₐₓ=1, m′ₘᵢₙ=2)  # max < min
            @test_throws ErrorException Ctor(R, ℓₘₐₓ; mₘₐₓ=1, mₘᵢₙ=2)
            @test_throws ErrorException Ctor(R, -1)
            # The defaults themselves are valid at every ℓₘₐₓ, including 0
            @test SphericalFunctions.m′ₘₐₓ(Ctor(R, 0)) == 0
        end
    end
end


@testitem "Wigner calculators batched rotors" begin
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
            if ℓ == 0
                recurrence!(batched, data, ℓ)
            else
                recurrence!(batched, ℓ)
            end
            blk = batched[ℓ]
            m′r = max(-ℓ, m′ₘᵢₙ(batched)):min(ℓ, m′ₘₐₓ(batched))
            mr = max(-ℓ, mₘᵢₙ(batched)):min(ℓ, mₘₐₓ(batched))
            @test ndims(blk) == 3
            @test axes(blk) == (1:N, m′r, mr)
            for i in 1:N
                recurrence!(single, data[i], ℓ)
                @test blk[i] == single[ℓ]
            end
        end
    end

    # Arbitrary ℓ order: backwards restarts, forwards advances, both reproduce the
    # sequential results exactly
    function check_order(batched, data)
        ref = OffsetVector([copy(recurrence!(batched, ℓ)[ℓ]) for ℓ in 0:ℓₘₐₓ], 0:ℓₘₐₓ)
        for ℓ in (4, 0, 7, 7, 2, 5, 1, 6, 3, 3, 0, 7)
            @test recurrence!(batched, ℓ)[ℓ] == ref[ℓ]
        end
        # Re-setting the rotor data and jumping straight to ℓ
        for ℓ in (7, 3, 0)
            @test recurrence!(batched, data, ℓ)[ℓ] == ref[ℓ]
        end
        # And from a fresh calculator of the same shape
        fresh = similar(batched)
        @test recurrence!(fresh, data, ℓₘₐₓ)[ℓₘₐₓ] == ref[ℓₘₐₓ]
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
        @test_throws ErrorException recurrence!(batched, rotors[1:N-1], 0)
        @test_throws ErrorException recurrence!(batched, rotors[1], 0)
        @test_throws ErrorException recurrence!(DCalculator(first(rotors), ℓₘₐₓ), rotors, 0)
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
        @test_throws ErrorException recurrence!(batched, βs[1:2], 0)
        @test_throws ErrorException recurrence!(batched, βs[1], 0)
        @test_throws ErrorException recurrence!(dCalculator(first(βs), ℓₘₐₓ), βs, 0)
    end
end


@testitem "Calculator block types are inferrable" begin
    import SphericalFunctions: DCalculator, dCalculator, sYlmCalculator
    import SphericalFunctions: recurrence!, isbatched, Nᵣ
    using Quaternionic: Rotor
    using Random

    # The payoff, measured the way a user's inner loop sees it: inside a function, where the
    # calculator's type is known, fetching a block allocates nothing.  Measuring at top level
    # instead would report the boxing of a dynamically dispatched call and prove nothing.
    blockallocs(c, ℓ) = (c[ℓ]; @allocated c[ℓ])
    blockallocs(c, ℓ, s) = (c[ℓ, s]; @allocated c[ℓ, s])

    # Whether `calc[ℓ]` returns a single block or a batch of them is a type parameter, not a
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
                @test isconcretetype(Base.return_types(getindex, (typeof(single), typeof(ℓ)))[1])
                @test (@inferred single[ℓ]) == single[ℓ]
                @test blockallocs(single, ℓ) == 0

                batched = Ctor(rotors, ℓₘₐₓ)
                @test isbatched(batched)
                @test isconcretetype(typeof(batched))
                recurrence!(batched, rotors, ℓ)
                @test isconcretetype(Base.return_types(getindex, (typeof(batched), typeof(ℓ)))[1])
                @test (@inferred batched[ℓ]) == batched[ℓ]
                @test blockallocs(batched, ℓ) == 0
            end

            # ... and the same for the harmonics, where the spin-weight argument is a second
            # type parameter, so `calc[ℓ]` has four shapes to keep concrete rather than two
            s = ℓₘₐₓ isa Integer ? 2 : 3//2
            srange = ℓₘₐₓ isa Integer ? (-2:2) : (-3//2:3//2)
            for spec ∈ (s, srange)
                single = sYlmCalculator(R, ℓₘₐₓ, spec)
                @test !isbatched(single)
                recurrence!(single, R, ℓ)
                @test isconcretetype(Base.return_types(getindex, (typeof(single), typeof(ℓ)))[1])
                @test (@inferred single[ℓ]) == single[ℓ]
                @test blockallocs(single, ℓ) == 0
                @test isconcretetype(Base.return_types(getindex, (typeof(single), typeof(ℓ), typeof(s)))[1])
                @test (@inferred single[ℓ, s]) == single[ℓ, s]
                @test blockallocs(single, ℓ, s) == 0

                batched = sYlmCalculator(rotors, ℓₘₐₓ, spec)
                @test isbatched(batched)
                recurrence!(batched, rotors, ℓ)
                @test isconcretetype(Base.return_types(getindex, (typeof(batched), typeof(ℓ)))[1])
                @test (@inferred batched[ℓ]) == batched[ℓ]
                @test blockallocs(batched, ℓ) == 0
                @test isconcretetype(Base.return_types(getindex, (typeof(batched), typeof(ℓ), typeof(s)))[1])
                @test (@inferred batched[ℓ, s]) == batched[ℓ, s]
                @test blockallocs(batched, ℓ, s) == 0
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

@testitem "Wigner calculators indexing errors" begin
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
            # Nothing has been computed yet, so no ℓ may be indexed — including values
            # outside 0:ℓₘₐₓ
            for ℓ in -1:ℓₘₐₓ+1
                @test_throws ErrorException calc[ℓ]
            end
            recurrence!(calc, R, 2)
            @test SphericalFunctions.ℓ(calc) == 2
            @test size(calc[2]) == (5, 5)
            # Only the current ℓ may be indexed, and the error says what to do about it
            for ℓ in (-1, 0, 1, 3, 4, 5)
                @test_throws "recurrence!" calc[ℓ]
            end
            # ℓ out of range for this calculator; the failed call leaves the current block
            # in place
            @test_throws ErrorException recurrence!(calc, ℓₘₐₓ + 1)
            @test_throws ErrorException recurrence!(calc, -1)
            @test SphericalFunctions.ℓ(calc) == 2
            reference = copy(calc[2])
            # `fill!` invalidates the current block ...
            fill!(calc, NaN)
            @test_throws ErrorException calc[2]
            # ... and recomputing from NaN-filled storage reproduces the result exactly, so
            # no uninitialized element is ever read
            recurrence!(calc, R, 2)
            @test calc[2] == reference
            @test !any(isnan, calc[2])
            @test_throws ErrorException recurrence!(calc, R, ℓₘₐₓ + 1)
        end
    end

    # A DCalculator needs the full rotor, not just β
    calcD = DCalculator(R, ℓₘₐₓ)
    @test_throws "Rotor" recurrence!(calcD, 0.3, 0)
    @test_throws "Rotor" recurrence!(calcD, cis(0.3), 0)
    @test_throws "Rotor" recurrence!(calcD, [0.3], 0)
    # ... whereas a dCalculator accepts any of the three forms
    calcd = dCalculator(R, ℓₘₐₓ)
    for input in (0.3, cis(0.3), R)
        # Deliberately *not* an `AbstractMatrix`: see the note on `AbstractWignerMatrix`.
        blk = recurrence!(calcd, input, 1)[1]
        @test blk isa WignerMatrix && eltype(blk) === Float64
    end
end


@testitem "D and d convenience functions" begin
    import SphericalFunctions: DCalculator, dCalculator, recurrence!, D, d
    using Quaternionic: Rotor, to_euler_phases
    import SphericalFunctions: WignerSeries, WignerMatrix
    using Random

    rng = Random.Xoshiro(7890)
    ℓₘₐₓ = 5
    R = randn(rng, Rotor{Float64})
    R₂ = randn(rng, Rotor{Float64})
    atol = 20 * ℓₘₐₓ * eps(Float64)
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
            @test 𝔇[ℓ] == recurrence!(calc, R, ℓ)[ℓ]
        end
        # Each block is an independent copy: computing more matrices, for a different rotor,
        # changes nothing
        snapshot = deepcopy(𝔇)
        𝔇₂ = D(R₂, ℓₘₐₓ)
        recurrence!(calc, R₂, ℓₘₐₓ)
        @test all(𝔇[ℓ] == snapshot[ℓ] for ℓ in 0:ℓₘₐₓ)
        @test 𝔇₂[ℓₘₐₓ] == calc[ℓₘₐₓ]
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
        # The element type follows the rotor's type
        @test eltype(D(Rotor{Float32}(R), 2)[2]) === ComplexF32
        @test eltype(D(R, 2)[2]) === ComplexF64
        @test eltype(D(Rotor{BigFloat}(R), 2)[2]) === Complex{BigFloat}
        # ... and lower precision only rounds the result
        𝔇₃₂ = D(Rotor{Float32}(R), 3)
        for ℓ in 0:3
            @test maximum(abs, 𝔇₃₂[ℓ] .- 𝔇[ℓ]) ≤ 20 * 3 * eps(Float32)
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
            @test dβ[ℓ] == recurrence!(calc, β, ℓ)[ℓ]
        end
        # Independent copies
        snapshot = deepcopy(dβ)
        d₂ = d(β / 3, ℓₘₐₓ)
        recurrence!(calc, β / 3, ℓₘₐₓ)
        @test all(dβ[ℓ] == snapshot[ℓ] for ℓ in 0:ℓₘₐₓ)
        @test d₂[ℓₘₐₓ] == calc[ℓₘₐₓ]
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
        # Block limits
        for input in (β, eⁱᵝ, R)
            dₗ = d(input, ℓₘₐₓ; limits...)
            @test axes(dₗ) == (0:ℓₘₐₓ,)
            for ℓ in 0:ℓₘₐₓ
                m′r = max(-ℓ, limits.m′ₘᵢₙ):min(ℓ, limits.m′ₘₐₓ)
                mr = max(-ℓ, limits.mₘᵢₙ):min(ℓ, limits.mₘₐₓ)
                @test axes(dₗ[ℓ]) == (m′r, mr)
                @test all(isapprox(dₗ[ℓ][m′, m], dβ[ℓ][m′, m]; atol) for m′ in m′r, m in mr)
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
        @test maximum(abs, d(Float32(β), 3)[3] .- dβ[3]) ≤ 20 * 3 * eps(Float32)
    end
end


@testitem "Wigner calculators vs independent references" setup=[HalfIntegerOracle] begin
    import SphericalFunctions
    import SphericalFunctions: DCalculator, dCalculator, recurrence!, D
    import .HalfIntegerOracle: d_oracle, D_oracle
    using Quaternionic: Rotor, Quaternion, components
    using LinearAlgebra: I, opnorm
    using Random

    # This item used to compare both calculators with `DenseWignerCalculator`, the older
    # single-rotor implementation kept as an internal oracle while the v3 engine was
    # written; that oracle has been deleted, and comparing the package to itself was never
    # an independent check anyway.  What it covered — every element of 𝔇ˡ and dˡ, for
    # several generic rotors, at every ℓ up to 6 — is covered here by two closed forms and
    # by identities the recurrence cannot satisfy by accident.

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
            if ℓ == 0
                recurrence!(calcD, R, ℓ)
                recurrence!(calcd, R, ℓ)
            else
                recurrence!(calcD, ℓ)
                recurrence!(calcd, ℓ)
            end
            𝔇ˡ = calcD[ℓ]
            dˡ = calcd[ℓ]
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


@testitem "Wigner calculators with narrow integer indices" begin
    import SphericalFunctions: DCalculator, dCalculator, recurrence!, D, d
    using Quaternionic: Rotor
    using Random

    # Indices all of one narrower integer type are kept as that type, all the way into the
    # calculator — whose `ℓ` reference must then be of that type too — and give the `Int`
    # result exactly.
    rng = Random.Xoshiro(20260917)
    R = randn(rng, Rotor{Float64})
    for IT in (Int8, Int16, Int32)
        @test D(R, IT(3)) == D(R, 3)
        @test d(R, IT(3)) == d(R, 3)
        for Ctor in (DCalculator, dCalculator)
            calc = Ctor(R, IT(3))
            @test calc.ℓ isa Base.RefValue{IT}
            @test recurrence!(calc, IT(2))[IT(2)] == recurrence!(Ctor(R, 3), 2)[2]
        end
    end
end
