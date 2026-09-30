# Tests for the iteration interface: `for (ℓ, 𝔇ˡ) ∈ calc`, `collect` on a calculator, the
# `set_R!`/`set_β!`/`set_θ!` family, the constructors that take the rotor data first —
# including the rule that the data alone fixes the element type, so that a mismatch is an
# error rather than a conversion — and `Ylm`.
#
# The oracle throughout is the package's own manual path — the `recurrence!` loop, and the
# convenience functions `D`, `d`, `sYlm` and `sYlm_matrix` that are built on it — because
# iteration is meant to do exactly the same arithmetic in exactly the same order and nothing
# else.  Every comparison is therefore bitwise, except in the two places where different
# code paths are being compared: a rotor's `β` reaches the recurrence through the
# quaternion's components rather than through `cis(β)`, and the `(θ, ϕ=0)` path skips the
# phase tables that a rotor at `ϕ = 0` still multiplies in.  Those two have a stated
# tolerance, with the measured error in a comment.

@testitem "Iteration reproduces D and d" begin
    import SphericalFunctions: DCalculator, dCalculator, D, d
    using Quaternionic: Rotor, to_euler_phases
    using Random

    rng = Random.Xoshiro(1414)
    N = 4
    rotors = randn(rng, Rotor{Float64}, N)
    # β ∈ [0, π] for each rotor, reached through the phase rather than through `acos`
    βs = [angle(to_euler_phases(R)[2]) for R ∈ rotors]

    @testset "ℓₘₐₓ = $ℓₘₐₓ" for ℓₘₐₓ ∈ (4, 7//2)
        # Nᵣ = 1.  `D` and `d` drive the same recurrence over the same ℓ in the same order,
        # so iteration has to reproduce them to the last bit.
        for (R, β) ∈ zip(rotors, βs)
            calc = DCalculator(R, ℓₘₐₓ)
            𝔇 = D(R, ℓₘₐₓ)
            visited = eltype(keys(calc))[]
            for (ℓ, 𝔇ˡ) ∈ calc
                @test all(𝔇ˡ[m′, m] == 𝔇[ℓ][m′, m] for m′ ∈ -ℓ:ℓ, m ∈ -ℓ:ℓ)
                push!(visited, ℓ)
            end
            @test visited == collect(keys(calc))

            calcd = dCalculator(β, ℓₘₐₓ)
            dm = d(β, ℓₘₐₓ)
            visitedd = eltype(keys(calcd))[]
            for (ℓ, dˡ) ∈ calcd
                @test all(dˡ[m′, m] == dm[ℓ][m′, m] for m′ ∈ -ℓ:ℓ, m ∈ -ℓ:ℓ)
                push!(visitedd, ℓ)
            end
            @test visitedd == collect(keys(calcd))
        end

        # Nᵣ > 1.  The batch reorders nothing: each rotor's slice of the batched block is
        # the single-rotor result exactly, so `D` and `d` are still the oracle.
        𝔇s = [D(R, ℓₘₐₓ) for R ∈ rotors]
        for (ℓ, 𝔇ˡ) ∈ DCalculator(rotors, ℓₘₐₓ)
            @test all(𝔇ˡ[i, m′, m] == 𝔇s[i][ℓ][m′, m] for i ∈ 1:N, m′ ∈ -ℓ:ℓ, m ∈ -ℓ:ℓ)
        end
        ds = [d(β, ℓₘₐₓ) for β ∈ βs]
        for (ℓ, dˡ) ∈ dCalculator(βs, ℓₘₐₓ)
            @test all(dˡ[i, m′, m] == ds[i][ℓ][m′, m] for i ∈ 1:N, m′ ∈ -ℓ:ℓ, m ∈ -ℓ:ℓ)
        end
    end
end


@testitem "Iteration reproduces sYlm, and the manual loop for half-integers" begin
    import SphericalFunctions: sYlmCalculator, sYlm, sYlm_matrix, recurrence!, Yindex
    using Quaternionic: Rotor
    using Random

    rng = Random.Xoshiro(2727)
    N = 4
    rotors = randn(rng, Rotor{Float64}, N)
    ℓₘₐₓ, sₘₐₓ = 5, 2

    # Integer ℓ.  `sYlm` and `sYlm_matrix` copy the very buffer these blocks view, so the
    # comparisons are bitwise, zeros for ℓ < |s| included.
    for R ∈ rotors
        calc = sYlmCalculator(R, ℓₘₐₓ, -sₘₐₓ:sₘₐₓ)
        for s ∈ -sₘₐₓ:sₘₐₓ
            Y = array_view(sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ=0))
            for (ℓ, ₛYₗ) ∈ calc
                @test all(ₛYₗ[s, m] == Y[Yindex(ℓ, m)] for m ∈ -ℓ:ℓ)
            end
        end
    end
    batched = sYlmCalculator(rotors, ℓₘₐₓ, -sₘₐₓ:sₘₐₓ)
    for s ∈ -sₘₐₓ:sₘₐₓ
        Y = sYlm_matrix(rotors, ℓₘₐₓ, s; ℓₘᵢₙ=0)
        for (ℓ, ₛYₗ) ∈ batched
            @test all(ₛYₗ[i, s, m] == Y[i, Yindex(ℓ, m)] for i ∈ 1:N, m ∈ -ℓ:ℓ)
        end
    end

    # Half-integer ℓ.  The flat interfaces are integer-only — the canonical `Yindex`
    # ordering is — so the oracle here is the manual `recurrence!` loop, and the batch is
    # compared with the single-rotor calculators.
    ℓₘₐₓₕ, sₘₐₓₕ = 5//2, 3//2
    for s ∈ (-3//2, -1//2, 1//2, 3//2)
        calc = sYlmCalculator(rotors[1], ℓₘₐₓₕ, -sₘₐₓₕ:sₘₐₓₕ)
        iterated = [ℓ => copy(ₛYₗ[s, :]) for (ℓ, ₛYₗ) ∈ calc]
        # Re-driving the same calculator by hand restarts the recurrence from ℓₘᵢₙ and runs
        # forward through the same ℓ, which is the same arithmetic again
        for (ℓ, ₛYₗ) ∈ iterated
            @test ₛYₗ == recurrence!(calc, ℓ)[s, :]
        end

        batchedₕ = sYlmCalculator(rotors, ℓₘₐₓₕ, -sₘₐₓₕ:sₘₐₓₕ)
        singles = [
            [copy(ₛYₗ[s, :]) for (_, ₛYₗ) ∈ sYlmCalculator(R, ℓₘₐₓₕ, -sₘₐₓₕ:sₘₐₓₕ)]
            for R ∈ rotors
        ]
        for (k, (ℓ, ₛYₗ)) ∈ enumerate(batchedₕ)
            row = ₛYₗ[:, s, :]
            @test all(row[i, m] == singles[i][k][m] for i ∈ 1:N, m ∈ -ℓ:ℓ)
        end
    end
end


@testitem "Iteration agrees with the manual loop" begin
    import SphericalFunctions: DCalculator, dCalculator, sYlmCalculator, recurrence!
    using Quaternionic: Rotor
    using Random

    rng = Random.Xoshiro(3141)
    rotors = randn(rng, Rotor{Float64}, 3)
    snapshot(iterable) = [ℓ => copy(block) for (ℓ, block) ∈ iterable]

    # Iteration is the manual `recurrence!` loop, one step per ℓ, and nothing else
    for calc ∈ (
        DCalculator(rotors[1], 4),
        dCalculator(0.7, 4),
        DCalculator(rotors, 4),
        DCalculator(rotors[1], 7//2),
        dCalculator(rotors[1], 7//2),
    )
        full = snapshot(calc)
        manual = [ℓ => copy(recurrence!(calc, ℓ)) for ℓ ∈ keys(calc)]
        @test manual == full
    end

    # One spin weight is the slice `ₛYₗ[s, :]`, whichever way the block was reached
    calc = sYlmCalculator(rotors[2], 5, -2:2)
    for s ∈ -2:2
        iterated = [ℓ => copy(ₛYₗ[s, :]) for (ℓ, ₛYₗ) ∈ calc]
        manual = [ℓ => copy(recurrence!(calc, ℓ)[s, :]) for ℓ ∈ keys(calc)]
        @test iterated == manual
        # ... and equals what a calculator built for that spin weight alone gives
        @test iterated == snapshot(sYlmCalculator(rotors[2], 5, s))
    end

    # A partial sweep is bit-for-bit the corresponding slice of a full one: the recurrence
    # runs through the intermediate ℓ either way
    calc = DCalculator(rotors[3], 5)
    full = snapshot(calc)
    partial(c, lo, hi) = [ℓ => copy(recurrence!(c, ℓ)) for ℓ ∈ lo:hi]
    @test partial(calc, 2, 4) == full[3:5]
    @test partial(calc, 2, 5) == full[3:end]
    @test partial(calc, 0, 1) == full[1:2]
    @test partial(calc, 5, 5) == full[6:6]
    calcₕ = DCalculator(rotors[3], 7//2)
    fullₕ = snapshot(calcₕ)
    @test partial(calcₕ, 3//2, 5//2) == fullₕ[2:3]
    calcY = sYlmCalculator(rotors[3], 5, -1:1)
    fullY = [ℓ => copy(ₛYₗ[-1, :]) for (ℓ, ₛYₗ) ∈ calcY]
    @test [ℓ => copy(recurrence!(calcY, ℓ)[-1, :]) for ℓ ∈ 1:3] == fullY[2:4]
end


@testitem "Calculator setters reach a freshly constructed state" setup=[RefusalChecks] begin
    import SphericalFunctions: DCalculator, dCalculator, HCalculator,
        sYlmCalculator, sλlmCalculator, HWedge, recurrence!, set_R!, set_β!, set_θ!,
        set_beta!, set_theta!
    using Quaternionic: Rotor, from_euler_angles
    using Random

    rng = Random.Xoshiro(4242)
    ℓₘₐₓ = 5
    rotors = randn(rng, Rotor{Float64}, 3)
    β = 0.7
    R = from_euler_angles(0.3, β, 1.1)
    θ = 1.3
    snapshot(iterable) = [ℓ => copy(block) for (ℓ, block) ∈ iterable]

    # `set_R!` on the calculators that need the whole rotor
    calc = DCalculator(rotors[1], ℓₘₐₓ)
    @test set_R!(calc, rotors[2]) === calc
    @test snapshot(calc) == snapshot(DCalculator(rotors[2], ℓₘₐₓ))
    calcₕ = DCalculator(rotors[1], 7//2)
    @test snapshot(set_R!(calcₕ, rotors[2])) == snapshot(DCalculator(rotors[2], 7//2))
    batched = DCalculator(rotors, ℓₘₐₓ)
    @test snapshot(set_R!(batched, reverse(rotors))) ==
        snapshot(DCalculator(reverse(rotors), ℓₘₐₓ))
    calcY = sYlmCalculator(rotors[1], ℓₘₐₓ, -2:2)
    @test set_R!(calcY, rotors[3]) === calcY
    @test snapshot(calcY) == snapshot(sYlmCalculator(rotors[3], ℓₘₐₓ, -2:2))

    # `set_β!`, in each of the three forms the angle may take
    calcd = dCalculator(0.25, ℓₘₐₓ)
    @test set_β!(calcd, β) === calcd
    @test snapshot(calcd) == snapshot(dCalculator(β, ℓₘₐₓ))
    @test snapshot(set_β!(calcd, cis(β))) == snapshot(dCalculator(cis(β), ℓₘₐₓ))
    @test snapshot(set_β!(dCalculator(0.25, 7//2), β)) ==
        snapshot(dCalculator(β, 7//2))
    # A rotor's β reaches the recurrence through the quaternion's components rather than
    # through `cis(β)`, so those two paths agree only to a few eps (measured 2.75 eps here)
    fromR = snapshot(set_β!(calcd, R))
    fromβ = snapshot(set_β!(calcd, β))
    @test maximum(maximum(abs.(a.second .- b.second)) for (a, b) ∈ zip(fromR, fromβ)) ≤ 8eps()

    # `set_β!` also serves the H calculator, which is not iterable: its wedge is read after
    # a manual step, entry by entry, in storage order
    wedge(H::HWedge) =
        [H[iᵣ, m′, m] for m′ ∈ -H.m′ₘₐₓ:H.m′ₘₐₓ for m ∈ abs(m′):H.ℓ for iᵣ ∈ 1:H.Nᵣ]
    H₁ = HCalculator(0.25, ℓₘₐₓ)
    @test set_β!(H₁, β) === H₁
    H₂ = HCalculator(β, ℓₘₐₓ)
    for ℓ ∈ 0:ℓₘₐₓ
        @test wedge(recurrence!(H₁, ℓ)) == wedge(recurrence!(H₂, ℓ))
    end
    wedgeᵣ = wedge(recurrence!(set_β!(H₁, R), ℓₘₐₓ))  # the rotor again, to 2.75 eps
    wedgeᵦ = wedge(recurrence!(H₂, ℓₘₐₓ))
    @test maximum(abs.(wedgeᵣ .- wedgeᵦ)) ≤ 8eps()

    # `set_θ!`, the entry point to the real functions ₛλₗₘ(θ) = ₛYₗₘ(θ, 0)
    calcθ = sYlmCalculator(0.25, ℓₘₐₓ, -2:2)
    @test set_θ!(calcθ, θ) === calcθ
    @test snapshot(calcθ) == snapshot(sYlmCalculator(θ, ℓₘₐₓ, -2:2))
    @test snapshot(set_θ!(calcθ, 0.0)) == snapshot(sYlmCalculator(0.0, ℓₘₐₓ, -2:2))

    # The ASCII spellings are the same functions
    @test set_beta! === set_β! && set_theta! === set_θ!

    # The wrong setter for a calculator is an error that names the right one
    @test refuses(() -> set_R!(dCalculator(β, ℓₘₐₓ), R), ArgumentError, "set_β!")
    @test refuses(() -> set_R!(HCalculator(β, ℓₘₐₓ), R), ArgumentError, "set_β!")
    @test refuses(() -> set_β!(DCalculator(R, ℓₘₐₓ), β), ArgumentError, "set_R!")
    @test refuses(() -> set_θ!(DCalculator(R, ℓₘₐₓ), θ), ArgumentError, "set_R!")
    @test refuses(() -> set_θ!(dCalculator(β, ℓₘₐₓ), θ), ArgumentError, "set_β!")
    @test refuses(() -> set_θ!(HCalculator(β, ℓₘₐₓ), θ), ArgumentError, "set_β!")
    # ... and says that both kinds of harmonic calculator evaluate at (θ, ϕ=0)
    @test refuses(() -> set_θ!(dCalculator(β, ℓₘₐₓ), θ), ArgumentError, "or an sλlmCalculator")
    # An sYlmCalculator takes rotors (`set_R!`) or angles (`set_θ!`), and neither setter
    # quietly accepts the other's data, singly or in a vector, so that `set_R!` never
    # switches a calculator to the (θ, ϕ=0) evaluation
    for calcY ∈ (sYlmCalculator(R, ℓₘₐₓ, -2:2), sYlmCalculator(rotors, ℓₘₐₓ, -2:2))
        @test refuses(() -> set_R!(calcY, θ), ArgumentError, "set_θ!")
        @test refuses(() -> set_R!(calcY, fill(θ, 3)), ArgumentError, "set_θ!")
        @test refuses(() -> set_θ!(calcY, R), ArgumentError, "set_R!")
        @test refuses(() -> set_θ!(calcY, fill(R, 3)), ArgumentError, "set_R!")
        @test calcY.phases[]  # still evaluating at the rotors
    end
    @test refuses(() -> set_θ!(sλlmCalculator(θ, ℓₘₐₓ, 1), R), ArgumentError, "nowhere to put")
    # β alone is not enough to place a point on the sphere, so an sYlmCalculator has no
    # `set_β!` at all
    @test_throws MethodError set_β!(sYlmCalculator(R, ℓₘₐₓ, -2:2), β)
end


@testitem "Iteration is restartable" begin
    import SphericalFunctions
    import SphericalFunctions: DCalculator, sYlmCalculator, recurrence!
    using Quaternionic: Rotor
    using Random

    rng = Random.Xoshiro(5353)
    R = randn(rng, Rotor{Float64})
    snapshot(iterable) = [ℓ => copy(block) for (ℓ, block) ∈ iterable]

    for calc ∈ (DCalculator(R, 5), DCalculator(R, 7//2))
        first_pass = snapshot(calc)
        # The iteration state is the next ℓ, not the calculator's internal position, so a
        # second pass starts over and gives the same values
        @test snapshot(calc) == first_pass
        # A pass cut short leaves nothing behind
        for (ℓ, _) ∈ calc
            ℓ ≥ SphericalFunctions.ℓₘᵢₙ(calc) + 1 && break
        end
        @test first(calc).first == SphericalFunctions.ℓₘᵢₙ(calc)
        @test snapshot(calc) == first_pass
        # Nor does a manual step in between
        recurrence!(calc, SphericalFunctions.ℓₘₐₓ(calc))
        @test snapshot(calc) == first_pass
    end

    # A partial sweep runs over whatever range it is given, wherever the calculator happens
    # to be standing, and leaves it fit for a later full pass
    calc = DCalculator(R, 5)
    full = snapshot(calc)
    partial = [ℓ => copy(recurrence!(calc, ℓ)) for ℓ ∈ 2:4]
    @test partial == full[3:5]
    @test snapshot(calc) == full
    @test [ℓ => copy(recurrence!(calc, ℓ)) for ℓ ∈ 2:4] == partial

    # The same for a spin-weighted calculator, whose block holds every spin weight at once,
    # so that two spin weights of one ℓ come from the one block rather than from two passes
    calcY = sYlmCalculator(R, 4, -1:1)
    first_pass = [ℓ => copy(ₛYₗ[1, :]) for (ℓ, ₛYₗ) ∈ calcY]
    for (ℓ, _) ∈ calcY
        ℓ == 2 && break
    end
    @test first(calcY).first == 0
    @test [ℓ => copy(ₛYₗ[1, :]) for (ℓ, ₛYₗ) ∈ calcY] == first_pass
    @test [ℓ => copy(ₛYₗ[-1, :]) for (ℓ, ₛYₗ) ∈ calcY] != first_pass
end


@testitem "One calculator over several rotors" begin
    import SphericalFunctions
    import SphericalFunctions: DCalculator, dCalculator, sYlmCalculator,
        recurrence!, set_R!, set_β!
    using Quaternionic: Rotor, to_euler_phases
    using Random

    rng = Random.Xoshiro(6464)
    rotors = randn(rng, Rotor{Float64}, 4)
    βs = [angle(to_euler_phases(R)[2]) for R ∈ rotors]
    snapshot(iterable) = [ℓ => copy(block) for (ℓ, block) ∈ iterable]

    # One workspace walked over many rotors is the whole point of the setters, and it must
    # give exactly what a fresh calculator per rotor would
    for ℓₘₐₓ ∈ (4, 7//2)
        calc = DCalculator(rotors[1], ℓₘₐₓ)
        for R ∈ rotors
            set_R!(calc, R)
            @test snapshot(calc) == snapshot(DCalculator(R, ℓₘₐₓ))
        end
        calcd = dCalculator(βs[1], ℓₘₐₓ)
        for β ∈ βs
            set_β!(calcd, β)
            @test snapshot(calcd) == snapshot(dCalculator(β, ℓₘₐₓ))
        end
        # The three-argument `recurrence!` replaces the data in the same way
        for R ∈ rotors
            recurrence!(calc, R, SphericalFunctions.ℓₘᵢₙ(calc))
            @test snapshot(calc) == snapshot(DCalculator(R, ℓₘₐₓ))
        end
    end

    calcY = sYlmCalculator(rotors[1], 4, -2:2)
    for R ∈ rotors
        set_R!(calcY, R)
        @test snapshot(calcY) == snapshot(sYlmCalculator(R, 4, -2:2))
    end

    # Batched, with the rotors rotated through the batch
    batched = DCalculator(rotors, 4)
    for k ∈ 0:length(rotors)-1
        R⃗ = circshift(rotors, k)
        set_R!(batched, R⃗)
        @test snapshot(batched) == snapshot(DCalculator(R⃗, 4))
    end
end


@testitem "Iteration allocates nothing and is inferrable" begin
    import SphericalFunctions
    import SphericalFunctions: DCalculator, dCalculator, sYlmCalculator, recurrence!
    using Quaternionic: Rotor
    using Random

    # Measured inside functions, as a user's inner loop sees it: at top level the boxing of
    # a dynamically dispatched call would be counted and would prove nothing.
    function trace(calc)  # Nᵣ = 1
        t = 0.0
        for (ℓ, block) ∈ calc
            t += abs(block[ℓ, ℓ])
        end
        t
    end
    function trace_batched(calc)  # Nᵣ > 1
        t = 0.0
        for (ℓ, block) ∈ calc
            t += abs(block[1, ℓ, ℓ])
        end
        t
    end
    function trace_spin(calc, s)  # Nᵣ = 1
        t = 0.0
        for (ℓ, block) ∈ calc
            t += abs(block[s, ℓ])
        end
        t
    end
    function trace_spin_batched(calc, s)  # Nᵣ > 1
        t = 0.0
        for (ℓ, block) ∈ calc
            t += abs(block[1, s, ℓ])
        end
        t
    end
    function trace_bare_spin(calc)  # a single-spin sYlmCalculator, Nᵣ = 1
        t = 0.0
        for (ℓ, block) ∈ calc
            t += abs(block[ℓ])
        end
        t
    end
    function trace_bare_spins(calc, s)  # a multi-spin sYlmCalculator, Nᵣ = 1
        t = 0.0
        for (ℓ, block) ∈ calc
            t += abs(block[s, ℓ])
        end
        t
    end
    function trace_range(calc, ℓₘᵢₙ, ℓₘₐₓ)  # a partial sweep, driven by hand
        t = 0.0
        for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ
            block = recurrence!(calc, ℓ)
            t += abs(block[ℓ, ℓ])
        end
        t
    end

    rng = Random.Xoshiro(7575)
    R = randn(rng, Rotor{Float64})
    rotors = randn(rng, Rotor{Float64}, 8)

    # Before Julia 1.12 a loop that reads the blocks costs a few small allocations in all,
    # rather than any per ℓ (measured at most 80 bytes on 1.10, with half-integer ℓ); from
    # 1.12 on it costs nothing.
    budget = VERSION ≥ v"1.12" ? 0 : 128
    for (ℓₘₐₓ, lo, hi) ∈ ((16, 2, 12), (25//2, 3//2, 21//2))
        single = DCalculator(R, ℓₘₐₓ)
        batched = DCalculator(rotors, ℓₘₐₓ)
        singled = dCalculator(0.7, ℓₘₐₓ)
        trace(single); trace_batched(batched); trace(singled)  # warm-up
        trace_range(single, lo, hi)
        @test (@allocated trace(single)) ≤ budget
        @test (@allocated trace_batched(batched)) ≤ budget
        @test (@allocated trace(singled)) ≤ budget
        @test (@allocated trace_range(single, lo, hi)) ≤ budget
    end

    calcY = sYlmCalculator(R, 16, -2:2)
    batchedY = sYlmCalculator(rotors, 16, -2:2)
    trace_spin(calcY, 1)
    @test (@allocated trace_spin(calcY, 1)) ≤ budget
    trace_spin_batched(batchedY, 1)
    @test (@allocated trace_spin_batched(batchedY, 1)) ≤ budget
    # Bare iteration of an sYlmCalculator, for both shapes of block
    calcY1 = sYlmCalculator(R, 16, 1)
    trace_bare_spin(calcY1)
    @test (@allocated trace_bare_spin(calcY1)) ≤ budget
    trace_bare_spins(calcY, 1)
    @test (@allocated trace_bare_spins(calcY, 1)) ≤ budget

    # `iterate` returns `nothing` at the end, so its type is a `Union` by construction: the
    # one-argument `@inferred` is *supposed* to fail on it, and the two-argument form —
    # naming the exact union — is what pins the type down.
    for calc ∈ (
        DCalculator(R, 4), DCalculator(rotors, 4), dCalculator(0.7, 4),
        DCalculator(R, 7//2), DCalculator(rotors, 7//2),
        sYlmCalculator(R, 4, 1), sYlmCalculator(rotors, 4, 1),
        sYlmCalculator(R, 4, -2:2), sYlmCalculator(rotors, 4, -2:2),
        sYlmCalculator(R, 7//2, 1//2), sYlmCalculator(R, 7//2, -3//2:3//2),
    )
        IT = typeof(SphericalFunctions.ℓₘᵢₙ(calc))
        @test isconcretetype(eltype(calc))
        @test (@inferred Union{Nothing, Tuple{eltype(calc), IT}} iterate(calc)) isa Tuple
        @test Base.return_types(iterate, (typeof(calc),))[1] ===
            Union{Nothing, Tuple{eltype(calc), IT}}
        @test Base.return_types(iterate, (typeof(calc), IT))[1] ===
            Union{Nothing, Tuple{eltype(calc), IT}}
    end
    # A hand-driven step is inferrable too, which is what makes the partial sweep above
    # allocation-free
    for calc ∈ (DCalculator(R, 4), sYlmCalculator(R, 4, -2:2))
        @test isconcretetype(Base.promote_op(recurrence!, typeof(calc), Int))
        @test (@inferred recurrence!(calc, 2)) !== nothing
    end
end


@testitem "collect copies every block" begin
    import SphericalFunctions: DCalculator, sYlmCalculator, recurrence!
    using Quaternionic: Rotor
    using Random

    rng = Random.Xoshiro(8686)
    R = randn(rng, Rotor{Float64})
    rotors = randn(rng, Rotor{Float64}, 3)

    calc = DCalculator(R, 4)
    reference = [ℓ => copy(block) for (ℓ, block) ∈ calc]
    v = collect(calc)
    @test v isa Vector
    @test v == reference
    # 1-based on the outside, natural indices on the inside — unlike `D` and `d`, which are
    # indexed by ℓ
    @test axes(v) == (1:length(calc),)
    @test [p.first for p ∈ v] == collect(keys(calc))
    @test axes(v[3].second) == (-2:2, -2:2)
    # The copies are independent of each other and of the calculator's storage
    v[2].second[1, -1] = 17
    @test v[2].second[1, -1] == 17
    @test all(v[k] == reference[k] for k ∈ eachindex(v) if k != 2)
    recurrence!(calc, 4)
    @test all(v[k] == reference[k] for k ∈ eachindex(v) if k != 2)
    # Copying is what changes the element type: iteration yields views, `collect` arrays
    @test eltype(v) !== eltype(calc)
    @test eltype(v) === typeof(v[1])

    # The same for a batch, for half-integer ℓ, and for the spin-weighted iterator
    vb = collect(DCalculator(rotors, 4))
    @test axes(vb) == (1:5,)
    @test axes(vb[3].second) == (1:3, -2:2, -2:2)
    vₕ = collect(DCalculator(R, 7//2))
    @test axes(vₕ) == (1:4,)
    @test [p.first for p ∈ vₕ] == [1//2, 3//2, 5//2, 7//2]
    @test vₕ[2].second == [ℓ => copy(b) for (ℓ, b) ∈ DCalculator(R, 7//2)][2].second
    calcY = sYlmCalculator(R, 4, -2:2)
    vY = collect(calcY)
    @test axes(vY) == (1:5,)
    @test axes(vY[3].second) == (-2:2, -2:2)
    @test vY == [ℓ => copy(block) for (ℓ, block) ∈ calcY]

    # `values` yields the blocks alone, as it does for a `WignerSeries`, and `collect` of it
    # copies them for the same reason
    vals = values(calc)
    @test length(vals) == length(calc)
    @test eltype(vals) === eltype(calc).parameters[2]
    @test [copy(b) for b ∈ vals] == last.(reference)
    cv = collect(vals)
    @test cv == last.(reference)
    @test all(parent(cv[k]) !== parent(cv[k+1]) for k ∈ 1:length(cv)-1)
    recurrence!(calc, 1)
    @test cv == last.(reference)
    @test collect(values(DCalculator(R, 7//2))) == last.(vₕ)
end


@testitem "Calculator container interface" setup=[RefusalChecks] begin
    import SphericalFunctions
    import SphericalFunctions: DCalculator, dCalculator, sYlmCalculator, recurrence!
    using Quaternionic: Rotor
    using Random

    rng = Random.Xoshiro(9797)
    R = randn(rng, Rotor{Float64})
    rotors = randn(rng, Rotor{Float64}, 3)

    for calc ∈ (
        DCalculator(R, 4), dCalculator(0.7, 4), DCalculator(rotors, 4),
        DCalculator(R, 7//2), dCalculator(0.7, 7//2),
    )
        ℓs = SphericalFunctions.ℓₘᵢₙ(calc):SphericalFunctions.ℓₘₐₓ(calc)
        @test keys(calc) == ℓs
        @test keys(calc) isa UnitRange  # never a `WignerRange`, which cannot be collected
        @test length(calc) == length(ℓs)
        @test length([ℓ for (ℓ, _) ∈ calc]) == length(calc)
        @test eltype(calc) === typeof(first(calc))
        @test isconcretetype(eltype(calc))
        @test pairs(calc) === calc
        @test Base.IteratorSize(typeof(calc)) === Base.HasLength()
        @test Base.IteratorEltype(typeof(calc)) === Base.HasEltype()
    end

    # `keys` is what a partial sweep is written against, and a sweep over any part of it
    # gives the blocks of those ℓ
    calc = DCalculator(R, 5)
    full = [copy(block) for (_, block) ∈ calc]
    for (lo, hi) ∈ ((0, 5), (2, 4), (3, 3), (0, 0))
        @test lo:hi ⊆ keys(calc)
        @test [copy(recurrence!(calc, ℓ)) for ℓ ∈ lo:hi] == full[(lo:hi) .+ 1]
    end

    # A type that does not fix the index type has no known element type, rather than one
    # that is looked for without end
    @test eltype(SphericalFunctions.DCalculator) === Any
    @test eltype(SphericalFunctions.dCalculator) === Any
    @test eltype(SphericalFunctions.sYlmCalculator) === Any
    @test eltype(SphericalFunctions.WignerCalculator) === Any

    # A multi-spin calculator iterates over its own whole block, and reports the same keys
    # as a single-spin one over the same ℓ
    calcY = sYlmCalculator(R, 4, -2:2)
    @test keys(calcY) == 0:4
    @test length(calcY) == 5
    @test isconcretetype(eltype(calcY))
    @test pairs(calcY) === calcY
    @test keys(sYlmCalculator(R, 4, -1)) == keys(calcY)

    # An ℓ outside the calculator's own is refused rather than silently clamped
    @test refuses(() -> recurrence!(calc, -1), ArgumentError, "out of bounds")
    @test refuses(() -> recurrence!(calc, 6), ArgumentError, "out of bounds")
    @test refuses(() -> recurrence!(calcY, 5), ArgumentError, "out of bounds")
end


@testitem "Calculators take their rotor data first" begin
    import SphericalFunctions
    import SphericalFunctions: DCalculator, dCalculator, HCalculator,
        sYlmCalculator, floattype
    using Quaternionic: Rotor, from_spherical_coordinates
    using Random

    rng = Random.Xoshiro(1212)
    R64 = randn(rng, Rotor{Float64})
    R32 = Rotor{Float32}(R64)
    rotors = randn(rng, Rotor{Float64}, 4)
    βs = [0.0, 0.3, 1.7, 2.9]
    blocktype(calc) = eltype(first(calc).second)

    # The element type follows the rotor data ...
    @test blocktype(DCalculator(R32, 3)) === ComplexF32
    @test blocktype(dCalculator(0.5f0, 3)) === Float32
    @test blocktype(dCalculator(Rotor{Float32}(R64), 3)) === Float32
    @test blocktype(sYlmCalculator(R32, 3, -1:1)) === ComplexF32
    @test blocktype(DCalculator(R64, 3)) === ComplexF64
    @test eltype(HCalculator(0.5f0, 3).eⁱᵝ) === ComplexF32
    # ... which is what `floattype` reports, for every kind of calculator
    @test floattype(DCalculator(R32, 3)) === Float32
    @test floattype(dCalculator(0.5f0, 3)) === Float32
    @test floattype(HCalculator(0.5f0, 3)) === Float32
    @test floattype(sYlmCalculator(R32, 3, -1:1)) === Float32
    @test floattype(DCalculator(R64, 3)) === Float64

    # ... and nothing else does: there is no element-type argument to override it, in either
    # direction, so a call that passes one has no method at all
    @test_throws MethodError DCalculator(R32, 3, Float64)
    @test_throws MethodError DCalculator(R64, 3, Float32)
    @test_throws MethodError dCalculator(0.5f0, 3, BigFloat)
    @test_throws MethodError HCalculator(0.5f0, 3, Float64)
    @test_throws MethodError sYlmCalculator(R32, 3, 1, Float64)
    # ... nor a constructor of the parametric type, which `show` prints, and which names the
    # element type
    @test_throws MethodError SphericalFunctions.DCalculator{Int, Float64}(R32, 3)
    @test_throws MethodError SphericalFunctions.WignerCalculator{Int, Float32, ComplexF32}(
        Rotor{BigFloat}(R64), 3
    )
    @test_throws MethodError SphericalFunctions.dCalculator{Int, Float64}(0.5f0, 3)
    # To compute in another type, build the rotor data in that type — which is also the
    # honest way to say it, since the type of the data is the claim being made about it
    @test blocktype(DCalculator(Rotor{BigFloat}(R64), 3)) === Complex{BigFloat}
    @test blocktype(dCalculator(big(0.5), 3)) === BigFloat
    @test blocktype(sYlmCalculator(Rotor{BigFloat}(R64), 3, -1:1)) === Complex{BigFloat}
    @test blocktype(DCalculator(Rotor{Float32}(R64), 3)) === ComplexF32

    # A vector argument gives a batch of exactly that length; `Nᵣ` is implied by it, and is
    # not a keyword argument anywhere
    @test SphericalFunctions.Nᵣ(DCalculator(rotors, 3)) == 4
    @test SphericalFunctions.Nᵣ(dCalculator(βs, 3)) == 4
    @test SphericalFunctions.Nᵣ(HCalculator(βs, 3)) == 4
    @test SphericalFunctions.Nᵣ(sYlmCalculator(rotors, 3, -1:1)) == 4
    @test SphericalFunctions.Nᵣ(DCalculator(rotors[1:1], 3)) == 1
    @test SphericalFunctions.Nᵣ(DCalculator(R64, 3)) == 1
    # A vector of one rotor is still a batch: batchedness follows the type of the data
    @test SphericalFunctions.isbatched(DCalculator(rotors[1:1], 3)) == true
    @test SphericalFunctions.isbatched(DCalculator(R64, 3)) == false

    # An empty vector describes no rotors at all, which is not a calculator
    @test_throws "at least one rotor" DCalculator(Rotor{Float64}[], 3)
    @test_throws "at least one rotor" dCalculator(Float64[], 3)
    @test_throws "at least one rotor" HCalculator(Float64[], 3)
    @test_throws "at least one rotor" sYlmCalculator(Rotor{Float64}[], 3, -1:1)

    # An ℓₘₐₓ-first spelling has no method at all, so such a call fails at the call site
    # rather than dispatching with the element type in the rotor's place
    @test_throws MethodError DCalculator(3)
    @test_throws MethodError DCalculator(3, Float64)
    @test_throws MethodError dCalculator(3, Float64)
    @test_throws MethodError HCalculator(3, Float64)
    @test_throws MethodError sYlmCalculator(3, 1, Float64)
    # A `Rational` ℓₘₐₓ selects the half-integer path
    @test SphericalFunctions.ℓₘₐₓ(DCalculator(R64, 7//2)) == 7//2
    @test first(DCalculator(R64, 7//2)).first == 1//2

    # `sYlmCalculator(θ, …)` is the (θ, ϕ=0) path: the real functions ₛλₗₘ(θ).  A rotor at
    # ϕ = 0 gives the same values to a few eps (measured ≤ 2.5 eps over these angles, and
    # exactly equal at several of them, because the phase it multiplies in is exactly 1
    # there), while a rotor at ϕ ≠ 0 does not give real values at all.  The `phases` flag is
    # what records which of the two kinds of data a calculator holds.
    for θ ∈ (0.0, 0.3, 1.0, 1.3, 2.9, Float64(π))
        calcθ = sYlmCalculator(θ, 4, -2:2)
        calcR = sYlmCalculator(Rotor(from_spherical_coordinates(θ, 0.0)), 4, -2:2)
        calcϕ = sYlmCalculator(Rotor(from_spherical_coordinates(θ, 0.9)), 4, -2:2)
        @test calcθ.phases[] == false
        @test calcR.phases[] == true
        for s ∈ -2:2
            fromθ = [copy(block[s, :]) for (_, block) ∈ calcθ]
            fromR = [copy(block[s, :]) for (_, block) ∈ calcR]
            @test all(all(iszero, imag.(block)) for block ∈ fromθ)
            @test maximum(maximum(abs.(a .- b)) for (a, b) ∈ zip(fromθ, fromR)) ≤ 8eps()
        end
        # Away from ϕ = 0 the harmonics are complex, so the angle path is not merely a
        # different spelling of a rotor
        if 0 < θ < π
            fromϕ = [copy(block[1, :]) for (_, block) ∈ calcϕ]
            @test any(any(!iszero, imag.(block)) for block ∈ fromϕ)
        end
    end
    # A vector of angles is a batch of them, just as a vector of rotors is
    @test SphericalFunctions.Nᵣ(sYlmCalculator([0.3, 1.3, 2.1], 4, -2:2)) == 3
    @test sYlmCalculator([0.3, 1.3, 2.1], 4, -2:2).phases[] == false
    let
        batched = sYlmCalculator([0.3, 1.3, 2.1], 4, -2:2)  # 1.3 is the second of the three
        single = [copy(block[1, :]) for (_, block) ∈ sYlmCalculator(1.3, 4, -2:2)]
        for (k, (ℓ, block)) ∈ enumerate(batched)
            @test all(block[2, 1, m] == single[k][m] for m ∈ -ℓ:ℓ)
        end
    end
end


@testitem "Angles in place of a rotor give exactly the rotor's values" begin
    import SphericalFunctions: D, DCalculator, sYlm, sYlmCalculator, Ylm, YlmCalculator
    using Quaternionic: from_euler_angles, from_spherical_coordinates
    using DoubleFloats: Double64

    # `===` on every element, so that the angle forms are seen to be the rotor forms rather
    # than merely close to them; a signed zero would be caught too.
    function identical(a, b)
        pairs_a, pairs_b = collect(a), collect(b)
        length(pairs_a) == length(pairs_b) && all(
            ℓa == ℓb && all(x === y for (x, y) ∈ zip(collect(blka), collect(blkb)))
            for ((ℓa, blka), (ℓb, blkb)) ∈ zip(pairs_a, pairs_b)
        )
    end

    for T ∈ (Float64, Float32, Double64)
        α, β, γ = T(3)/10, T(11)/10, T(24)/10
        θ, ϕ = T(7)/10, T(4)
        R = from_euler_angles(α, β, γ)
        Rθϕ = from_spherical_coordinates(θ, ϕ)
        for (ℓₘₐₓ, s, spins) ∈ ((4, 2, -2:2), (7//2, 1//2, -3//2:3//2))
            @test identical(D(α, β, γ, ℓₘₐₓ), D(R, ℓₘₐₓ))
            @test identical(DCalculator(α, β, γ, ℓₘₐₓ), DCalculator(R, ℓₘₐₓ))
            @test identical(sYlm(θ, ϕ, ℓₘₐₓ, s), sYlm(Rθϕ, ℓₘₐₓ, s))
            @test identical(sYlm(θ, ϕ, ℓₘₐₓ, spins), sYlm(Rθϕ, ℓₘₐₓ, spins))
            @test identical(sYlmCalculator(θ, ϕ, ℓₘₐₓ, s), sYlmCalculator(Rθϕ, ℓₘₐₓ, s))
        end
        # The keywords pass through
        @test identical(D(α, β, γ, 3; mₘₐₓ=1), D(R, 3; mₘₐₓ=1))
        @test identical(
            DCalculator(α, β, γ, 3; m′ₘₐₓ=2, mₘᵢₙ=-1), DCalculator(R, 3; m′ₘₐₓ=2, mₘᵢₙ=-1)
        )
        @test identical(sYlm(θ, ϕ, 5, 1; ℓₘᵢₙ=3), sYlm(Rθϕ, 5, 1; ℓₘᵢₙ=3))
        @test identical(Ylm(θ, ϕ, 5), Ylm(Rθϕ, 5))
        @test identical(Ylm(θ, ϕ, 5; ℓₘᵢₙ=2), Ylm(Rθϕ, 5; ℓₘᵢₙ=2))
        @test identical(YlmCalculator(θ, ϕ, 4), YlmCalculator(Rθϕ, 4))
    end

    # The indices obey the same rules as ever, and a forgotten spin weight has no method
    @test_throws ArgumentError Ylm(0.7, 4.0, 1//2)
    @test_throws ArgumentError YlmCalculator(0.7, 4.0, 3//2)
    @test_throws ArgumentError D(0.3, 1.1, 2.4, Int8(3))
    @test_throws ArgumentError sYlm(0.7, 4.0, 2, 3)
    @test_throws MethodError sYlmCalculator(0.7, 4.0, 3)
end


@testitem "A calculator's element type is fixed by its data" begin
    import SphericalFunctions
    import SphericalFunctions: DCalculator, dCalculator, HCalculator,
        sYlmCalculator, D, d, sYlm, sYlm_matrix, Ylm, floattype, set_R!, set_β!, set_θ!,
        array_view, recurrence!
    using Quaternionic: Rotor, Quaternion, QuatVec, rotor
    using StaticArrays: SVector
    using Random

    rng = Random.Xoshiro(1515)
    rotors = randn(rng, Rotor{Float64}, 2)
    rotors32 = Rotor{Float32}.(rotors)
    rotorsb = Rotor{BigFloat}.(rotors)

    # The setters replace a calculator's data, not the type it works in, so the new data must
    # give that same type.  A mismatch is an error rather than a silent conversion: narrowing
    # a BigFloat rotor into a Float64 calculator would throw away precision that nobody chose
    # to throw away, and widening a Float32 one would claim an accuracy that is not there.
    @test_throws "works in Float64" set_R!(DCalculator(rotors[1], 3), rotors32[1])
    @test_throws "works in Float32" set_R!(DCalculator(rotors32[1], 3), rotors[1])
    @test_throws "works in Float64" set_R!(DCalculator(rotors[1], 3), rotorsb[1])
    @test_throws "works in Float64" set_R!(DCalculator(rotors, 3), rotors32)
    @test_throws "works in Float64" set_R!(sYlmCalculator(rotors[1], 3, -1:1), rotors32[1])
    @test_throws "works in Float64" set_R!(sYlmCalculator(rotors, 3, -1:1), rotors32)
    # `set_β!` accepts the angle, the phase e^{iβ} or a rotor, and the rule reaches all three
    @test_throws "works in Float64" set_β!(dCalculator(0.25, 3), 0.5f0)
    @test_throws "works in Float64" set_β!(dCalculator(0.25, 3), cis(0.5f0))
    @test_throws "works in Float64" set_β!(dCalculator(0.25, 3), rotors32[1])
    @test_throws "works in Float32" set_β!(dCalculator(0.25f0, 3), 0.5)
    @test_throws "works in Float64" set_β!(HCalculator(0.25, 3), 0.5f0)
    @test_throws "works in Float64" set_β!(HCalculator(0.25, 3), cis(0.5f0))
    @test_throws "works in Float64" set_β!(HCalculator(0.25, 3), rotors32[1])
    @test_throws "works in Float64" set_θ!(sYlmCalculator(0.25, 3, -1:1), 0.5f0)
    @test_throws "works in Float32" set_θ!(sYlmCalculator(0.25f0, 3, -1:1), 0.5)
    # `recurrence!(calc, data, ℓ)` replaces the data too, and follows the same rule rather than
    # silently converting, which would return a BigFloat rotor's block as a ComplexF64 one.  The
    # check comes before anything is written, so a refused call leaves the calculator as it was.
    import SphericalFunctions: recurrence!
    rotorsb1 = Rotor{BigFloat}(rotors[1])
    @test_throws "works in Float64" recurrence!(DCalculator(rotors[1], 3), rotorsb1, 2)
    @test_throws "works in Float64" recurrence!(DCalculator(rotors[1], 3), rotors32[1], 2)
    @test_throws "works in Float64" recurrence!(dCalculator(0.25, 3), 0.5f0, 2)
    @test_throws "works in Float64" recurrence!(dCalculator(0.25, 3), rotorsb1, 2)
    @test_throws "works in Float64" recurrence!(HCalculator(0.25, 3), 0.5f0, 2)
    @test_throws "works in Float64" recurrence!(sYlmCalculator(rotors[1], 3, 1), rotorsb1, 2)
    @test_throws "works in Float64" recurrence!(sYlmCalculator(0.25, 3, 1), 0.5f0, 2)
    let c = DCalculator(rotors[1], 3)
        before = copy(recurrence!(c, 2))
        @test_throws "works in Float64" recurrence!(c, rotorsb1, 2)
        @test recurrence!(c, 2) == before
        @test recurrence!(c, rotors[2], 2) == D(rotors[2], 3)[2]  # the matching type works
    end
    # `similar(calc, data)` builds a second workspace of exactly the calculator's type, so it
    # is just as strict; the Nᵣ check it has always had is tested with the rest of `similar`
    @test_throws "works in Float64" similar(DCalculator(rotors[1], 3), rotors32[1])
    @test_throws "works in Float64" similar(sYlmCalculator(rotors[1], 3, -1:1), rotors32[1])
    @test_throws "works in Float64" similar(HCalculator(0.25, 3), 0.5f0)
    # Data of the calculator's own type is accepted, in every one of these forms
    @test floattype(set_R!(DCalculator(rotors32[1], 3), rotors32[2])) === Float32
    @test floattype(set_β!(dCalculator(0.25f0, 3), 0.5f0)) === Float32
    @test floattype(set_β!(HCalculator(0.25, 3), cis(0.5))) === Float64
    @test floattype(set_θ!(sYlmCalculator(0.25, 3, -1:1), 0.5)) === Float64
    @test floattype(similar(DCalculator(rotorsb[1], 3), rotorsb[2])) === BigFloat

    # A calculator of `Float32` cannot be pointed at a `Float64` rotor, either, and one
    # pointed at a rotor of its own type gives what `sYlm` gives, in that type
    @test_throws "works in Float32" set_R!(sYlmCalculator(rotors32[1], 4, -1:1), rotors[1])
    for (c, Y) ∈ (
        (set_R!(sYlmCalculator(rotors[2], 4, 1), rotors[1]), sYlm(rotors[1], 4, 1)),
        (set_R!(sYlmCalculator(rotors32[2], 4, 1), rotors32[1]), sYlm(rotors32[1], 4, 1)),
    )
        @test all(
            eltype(array_view(b)) === eltype(array_view(Y)) && b == Y[ℓ]
            for (ℓ, b) ∈ c if ℓ ∈ keys(Y)
        )
    end

    # A vector of rotor data must say what it holds, so an abstract or ambiguous element type
    # is refused instead of guessed at, as it would be by re-boxing such a vector into a
    # `Vector{AbstractQuaternion}` and falling back on Float64.
    anyvector = Any[rotors[1], rotors[2]]
    abstractvector = Rotor[rotors[1], rotors[2]]
    mixedvector = Union{Rotor{Float64}, Rotor{Float32}}[rotors[1], rotors32[2]]
    for bad ∈ (anyvector, abstractvector, mixedvector)
        @test_throws "Cannot build a calculator" DCalculator(bad, 3)
        @test_throws "Cannot build a calculator" sYlmCalculator(bad, 3, -1:1)
        @test_throws "Cannot build a calculator" set_R!(DCalculator(rotors, 3), bad)
        @test_throws "Cannot build a calculator" set_R!(sYlmCalculator(rotors, 3, -1:1), bad)
    end
    @test_throws "Cannot build a calculator" dCalculator(Any[0.3, 0.5], 3)
    @test_throws "Cannot build a calculator" HCalculator(Any[0.3, 0.5], 3)
    @test_throws "Cannot build a calculator" dCalculator(Number[0.3, 0.5], 3)
    # The error names the offending type and says what to do about it
    err = try DCalculator(anyvector, 3) catch e; e end
    @test err isa ArgumentError
    @test occursin("Vector{Any}", err.msg)
    @test occursin("should be converted", err.msg)

    # Rotations are taken as `Rotor`s, or as `Quaternion`s, which denote the rotations of
    # their normalizations.  A `QuatVec` is a vector rather than a rotation at all, and is
    # not silently reinterpreted: every entry point — the convenience functions, `w(R)`, the
    # transforms, the calculators and the setters — refuses it with a message that names
    # `exp(v/2)`.
    q = Quaternion(0.3, 0.5, 0.7, 0.11)
    qv = QuatVec(0.0, 0.0, 1.0)
    w = SphericalFunctions.ModeWeights(randn(ComplexF64, 9), 0)
    for bad ∈ (qv,)
        @test_throws "Rotations are taken as" D(bad, 2)
        @test_throws "Rotations are taken as" D(bad, 3//2)
        @test_throws "Rotations are taken as" d(bad, 2)
        @test_throws "Rotations are taken as" sYlm(bad, 2, 0)
        @test_throws "Rotations are taken as" sYlm([bad, bad], 2, 0)
        @test_throws "Rotations are taken as" Ylm(bad, 2)
        @test_throws "Rotations are taken as" Ylm([bad], 2)
        @test_throws "Rotations are taken as" sYlm_matrix([bad, bad], 2, 0)
        @test_throws "Rotations are taken as" w(bad)
        @test_throws "Rotations are taken as" w([bad])
        @test_throws "Rotations are taken as" SphericalFunctions.SSHTMatrix(0, 0; Rθϕ=[bad])
        @test_throws "Rotations are taken as" DCalculator(bad, 2)
        @test_throws "Rotations are taken as" dCalculator(bad, 2)
        @test_throws "Rotations are taken as" HCalculator(bad, 2)
        @test_throws "Rotations are taken as" sYlmCalculator(bad, 2, -0:0)
        @test_throws "Rotations are taken as" set_R!(DCalculator(rotors[1], 2), bad)
        @test_throws "Rotations are taken as" DCalculator([bad, bad], 2)
    end
    # The message names what to write instead, and that spelling works
    @test floattype(DCalculator(exp(qv/2), 2)) === Float64
    # Every entry point takes a `Quaternion` as the rotation of its normalization
    R = rotor(q)
    @test D(q, 2) ≈ D(R, 2)
    @test D(q, 3//2) ≈ D(R, 3//2)
    @test d(q, 2) ≈ d(R, 2)
    @test array_view(sYlm(q, 2, 0)) ≈ array_view(sYlm(R, 2, 0))
    @test array_view(sYlm([q, q], 2, 0)) ≈ array_view(sYlm([R, R], 2, 0))
    @test array_view(Ylm(q, 2)) ≈ array_view(Ylm(R, 2))
    @test sYlm_matrix([q, q], 2, 0) ≈ sYlm_matrix([R, R], 2, 0)
    @test w(q) ≈ w(R)
    @test w([q]) ≈ w([R])
    @test SphericalFunctions.rotors(SphericalFunctions.SSHTMatrix(0, 0; Rθϕ=[q])) ≈ [R]
    @test recurrence!(DCalculator(q, 2), 2) ≈ recurrence!(DCalculator(R, 2), 2)
    @test floattype(DCalculator(q, 2)) === floattype(dCalculator(q, 2)) === Float64
    @test floattype(HCalculator(q, 2)) === floattype(sYlmCalculator(q, 2, -0:0)) === Float64
    @test recurrence!(set_R!(DCalculator(rotors[1], 2), q), 2) ≈ recurrence!(DCalculator(R, 2), 2)
    @test SphericalFunctions.Nᵣ(DCalculator([q, q], 2)) == 2

    # Concretely typed data is untouched by any of this, in each of its forms
    @test SphericalFunctions.Nᵣ(DCalculator(rotors, 3)) == 2
    # (A `Vector{Quaternion}` is one of these forms too; it is covered by the tests above.)
    @test SphericalFunctions.Nᵣ(DCalculator(SVector{2}(rotors[1], rotors[2]), 3)) == 2
    @test SphericalFunctions.Nᵣ(dCalculator([0.3, 0.5], 3)) == 2
    @test SphericalFunctions.Nᵣ(dCalculator(cis.([0.3, 0.5]), 3)) == 2
    @test SphericalFunctions.Nᵣ(sYlmCalculator([0.3, 0.5], 3, -1:1)) == 2
    @test floattype(DCalculator(SVector{2}(rotors32[1], rotors32[2]), 3)) === Float32
end


@testitem "similar retains the rotor data" begin
    import SphericalFunctions
    import SphericalFunctions: DCalculator, dCalculator, sYlmCalculator
    using Quaternionic: Rotor
    using Random

    rng = Random.Xoshiro(2323)
    rotors = randn(rng, Rotor{Float64}, 3)
    snapshot(iterable) = [ℓ => copy(block) for (ℓ, block) ∈ iterable]

    # A calculator without rotor data is not representable, so `similar` has to retain the
    # rotor data; only the computed results are absent, and iterating recomputes them
    for calc ∈ (
        DCalculator(rotors[1], 4),
        DCalculator(rotors, 4; m′ₘₐₓ=2, mₘᵢₙ=-3),
        dCalculator(0.7, 4),
        DCalculator(rotors[1], 7//2),
        dCalculator(rotors, 7//2),
    )
        reference = snapshot(calc)
        s = similar(calc)
        @test typeof(s) === typeof(calc)
        @test s !== calc
        @test SphericalFunctions.Nᵣ(s) == SphericalFunctions.Nᵣ(calc)
        @test snapshot(s) == reference
        @test snapshot(calc) == reference  # and the original is undisturbed
    end

    # `similar(calc, R)` keeps the parameters and replaces the data, and needs the same Nᵣ
    calc = DCalculator(rotors[1], 4)
    s = similar(calc, rotors[2])
    @test typeof(s) === typeof(calc)
    @test snapshot(s) == snapshot(DCalculator(rotors[2], 4))
    @test_throws "Nᵣ" similar(calc, rotors)
    batched = DCalculator(rotors, 4)
    @test snapshot(similar(batched, reverse(rotors))) ==
        snapshot(DCalculator(reverse(rotors), 4))
    @test_throws "Nᵣ" similar(batched, rotors[1])

    # For an sYlmCalculator the data includes the `phases` flag, which is what distinguishes
    # the (θ, ϕ=0) path; `similar` of a θ-constructed calculator must still give θ values
    calcθ = sYlmCalculator(1.3, 4, -2:2)
    referenceθ = snapshot(calcθ)
    sθ = similar(calcθ)
    @test typeof(sθ) === typeof(calcθ)
    @test sθ.phases[] == calcθ.phases[] == false
    @test snapshot(sθ) == referenceθ
    @test snapshot(similar(calcθ, 0.9)) == snapshot(sYlmCalculator(0.9, 4, -2:2))
    calcR = sYlmCalculator(rotors[1], 4, -2:2)
    sR = similar(calcR)
    @test sR.phases[] == calcR.phases[] == true
    @test snapshot(sR) == snapshot(calcR)
    @test snapshot(similar(calcR, rotors[2])) == snapshot(sYlmCalculator(rotors[2], 4, -2:2))
end


@testitem "Ylm is the spin-zero case of sYlm" setup=[RefusalChecks] begin
    import SphericalFunctions: Ylm, sYlm
    using Quaternionic: Rotor
    using Random

    rng = Random.Xoshiro(3434)
    rotors = randn(rng, Rotor{Float64}, 3)

    for R ∈ rotors
        @test array_view(Ylm(R, 5)) == array_view(sYlm(R, 5, 0))
        @test array_view(Ylm(R, 5; ℓₘᵢₙ=0)) == array_view(sYlm(R, 5, 0; ℓₘᵢₙ=0))
        @test array_view(Ylm(R, 5; ℓₘᵢₙ=2)) == array_view(sYlm(R, 5, 0; ℓₘᵢₙ=2))
        @test array_view(Ylm(R, 5; ell_min=2)) == array_view(sYlm(R, 5, 0; ℓₘᵢₙ=2))
        @test array_view(Ylm(R, 0)) == array_view(sYlm(R, 0, 0))
        # The two agree in whatever type the rotor is given in, which is the only thing
        # that decides it
        @test array_view(Ylm(Rotor{Float32}(R), 5)) == array_view(sYlm(Rotor{Float32}(R), 5, 0))
        @test array_view(Ylm(Rotor{BigFloat}(R), 5)) == array_view(sYlm(Rotor{BigFloat}(R), 5, 0))
    end
    @test eltype(array_view(Ylm(Rotor{Float32}(rotors[1]), 3))) === ComplexF32
    @test eltype(array_view(Ylm(Rotor{BigFloat}(rotors[1]), 3))) === Complex{BigFloat}
    # Neither function has an element-type keyword argument
    @test_throws MethodError array_view(Ylm(rotors[1], 3; T=BigFloat))
    @test_throws MethodError array_view(sYlm(rotors[1], 3, 0; T=BigFloat))

    # Half-integer ℓ goes with half-integer spin weight, so these functions have no
    # half-integer analogue, and a half-integer index is refused with a message that says so
    # and names the functions that do take one
    for f ∈ (() -> Ylm(rotors[1], 7//2), () -> Ylm(rotors, 7//2))
        @test refuses(f, ArgumentError, "does not accept half-odd-integers")
    end
    @test refuses(() -> Ylm(rotors[1], 5; ℓₘᵢₙ=1//2), ArgumentError, "keyword argument `ℓₘᵢₙ`")
    @test refuses(() -> Ylm(rotors[1], 7//2), ArgumentError, "`sYlm` and `sYlmCalculator`")
    @test refuses(() -> YlmCalculator(rotors[1], 7//2), ArgumentError, "`sYlm` and `sYlmCalculator`")
    # An integer of another type than `Int` is refused as everywhere else, and a
    # floating-point ℓ is not an index at all
    @test refuses(() -> Ylm(rotors[1], Int32(3)), ArgumentError, "narrower than `Int`")
    @test refuses(() -> YlmCalculator(rotors[1], Int32(3)), ArgumentError, "narrower than `Int`")
    @test_throws MethodError array_view(Ylm(rotors[1], 3.0))
end


@testitem "Iteration's deliberate refusal" begin
    import SphericalFunctions: DCalculator, HCalculator, sYlmCalculator,
        recurrence!
    using Quaternionic: Rotor
    using Random

    rng = Random.Xoshiro(4545)
    R = randn(rng, Rotor{Float64})

    # An sYlmCalculator is built for the spin weights it serves, so bare iteration has
    # something to yield; what is refused is a spin weight it was not built for, which is out
    # of bounds of the block.
    calcY = sYlmCalculator(R, 3, -1:1)
    @test first(calcY).first == 0
    @test length(collect(calcY)) == 4
    @test_throws BoundsError recurrence!(calcY, 0)[2, :]

    # An HCalculator's only block is the wedge itself — one mutable object handed back
    # by identity, which `copy` cannot preserve — so it is not iterable at all.  The error
    # names the manual loop it has always had.
    calcH = HCalculator(0.7, 3)
    err = try iterate(calcH) catch e; e end
    @test err isa ArgumentError
    @test occursin("not iterable", err.msg)
    @test occursin("recurrence!", err.msg)
    @test occursin("DCalculator", err.msg)
    @test_throws "not iterable" [x for x ∈ calcH]
    # And that manual loop is unaffected
    @test recurrence!(calcH, 2).ℓ == 2

    # Neither refusal touches the calculators that are iterable
    @test length(collect(DCalculator(R, 3))) == 4
end
