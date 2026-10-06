# Tests of the real flavor of the harmonics, `sλlmCalculator`, which shares its struct, its
# recurrence, and its blocks with the complex `sYlmCalculator` and differs only in the
# number type it stores.
#
# The load-bearing item is the first: ₛλₗₘ must be *bit-for-bit* the real part of the
# complex harmonic at (θ, 0) for an integer spin weight, and ± its imaginary part for a
# half-odd one, because that is the whole claim.  The relation is written out here rather
# than taken from the package, which is what makes this a regression test rather than a
# tautology:
#
#     λreal(x, ::Integer) = real(x)
#     λreal(x, s::HalfOddInteger) = ifelse(mod(2s, 4) == 1, imag(x), -imag(x))
#
# The complex reference must be given the *angle*, at construction or with `set_θ!`.
# Against a rotor built from the same θ the values agree only to about 5 ulps, because the
# rotor path multiplies by a unit phase that the angle path never forms — a rounding
# difference, not a discrepancy, and asserted as such at the end of the first item.

@testitem "sλlmCalculator is bit-for-bit the real part of the complex harmonics" begin
    import SphericalFunctions: sλlmCalculator, sYlmCalculator, sYlm, set_θ!, recurrence!,
        floattype, array_view
    using DoubleFloats: Double64
    using Quaternionic: from_spherical_coordinates

    # The relation of the header, written out
    λref(x, s) = isinteger(s) ? real(x) : (mod(Int(2s), 4) == 1 ? imag(x) : -imag(x))

    for T ∈ (Float64, Double64)
        for (ℓmax, spins) ∈ ((4, (0, 1, -2, 2)), (7//2, (1//2, -1//2, 3//2, -3//2)))
            for s ∈ spins
                θ = T(7) / 10
                # One calculator of each flavor is built from the angle, and another is
                # built from a different angle and then given this one with `set_θ!`.
                # Since `set_θ!` replaces all of the rotor data, the two calculators of one
                # flavor must agree exactly.
                cY = sYlmCalculator(θ, ℓmax, s)          # complex, angle mode
                cλ = sλlmCalculator(θ, ℓmax, s)
                cY′ = set_θ!(sYlmCalculator(T(1) / 10, ℓmax, s), θ)
                cλ′ = set_θ!(sλlmCalculator(T(1) / 10, ℓmax, s), θ)
                @test floattype(cλ) === T && eltype(cλ.Yˡ) === T
                for ℓ ∈ abs(s):ℓmax
                    blkY = recurrence!(cY, ℓ)
                    blkλ = recurrence!(cλ, ℓ)
                    blkY′ = recurrence!(cY′, ℓ)
                    blkλ′ = recurrence!(cλ′, ℓ)
                    for m ∈ -ℓ:ℓ
                        # `===` rather than `==`, so that a signed zero would be caught too.
                        @test blkλ[m] === λref(blkY[m], s)
                        @test blkλ′[m] === λref(blkY′[m], s)
                        @test blkλ′[m] === blkλ[m]
                    end
                end
            end
        end
    end

    # The angle path and the rotor path differ only by rounding: measured worst case over
    # every s and ℓ below was 5 ulps of the largest value, and 20 is asserted.
    for s ∈ (0, -2, 1//2, -3//2)
        ℓmax = s isa Rational ? 7//2 : 4
        θ = 0.7
        a = array_view(sYlm(from_spherical_coordinates(θ, 0.0), ℓmax, s))
        cλ = sλlmCalculator(θ, ℓmax, s)
        b = reduce(vcat, [Array(recurrence!(cλ, ℓ)) for ℓ ∈ abs(s):ℓmax])
        scale = maximum(abs, b)
        @test maximum(abs, λref.(a, s) .- b) < 20eps(scale)
    end
end

@testitem "sλlmCalculator allocates no phase tables" begin
    import SphericalFunctions: sλlmCalculator, number_type

    for (ℓmax, s) ∈ ((6, -2), (11//2, 3//2))
        cλ = sλlmCalculator(0.3, ℓmax, s)
        cY = sYlmCalculator(0.3, ℓmax, s)
        # The tables are empty rather than merely unread, which is the point of the `K`
        # trick in `allocate_Y`: they cost nothing at all for the real flavor.
        @test isempty(cλ.engine.Z₊) && isempty(cλ.engine.Z₋)
        @test !isempty(cY.engine.Z₊) && !isempty(cY.engine.Z₋)
        # The rotor axis, which comes first, is still there, and only the powers are missing
        @test size(cλ.engine.Z₊) == (size(cY.engine.Z₊, 1), 0)
        @test size(cY.engine.Z₊) == (1, 2SphericalFunctions.power_extent(ℓmax, abs(s)) + 1)
        @test number_type(cλ) === Float64
        @test number_type(cY) === ComplexF64
        @test eltype(cλ.Yˡ) === Float64
        # The shared H recurrence is the same size either way — that is what makes this a
        # saving rather than a trade.
        @test size(parent(cλ.engine.H.Hˡ)) == size(parent(cY.engine.H.Hˡ))
    end
end

@testitem "sλlmCalculator refuses rotors" setup=[RefusalChecks] begin
    import SphericalFunctions: sλlmCalculator
    using Quaternionic: Rotor, from_spherical_coordinates
    using Random

    rng = Random.Xoshiro(20260920)
    R = randn(rng, Rotor{Float64})
    cλ = sλlmCalculator(0.3, 4, -2)

    # A rotor specifies α and γ, whose phases a real calculator has nowhere to put.  It says
    # so rather than silently dropping them.
    @test refuses(() -> set_R!(cλ, R), ArgumentError, "nowhere to put the α and γ phases")
    @test refuses(() -> set_R!(cλ, [R]), ArgumentError, "nowhere to put the α and γ phases")
    @test refuses(() -> set_θ!(cλ, R), ArgumentError, "nowhere to put the α and γ phases")
    @test refuses(() -> sλlmCalculator(R, 4, -2), ArgumentError, "cannot take a rotor")
    @test refuses(
        () -> SphericalFunctions.set_rotors!(cλ, "nonsense"), ArgumentError, "needs angles θ"
    )

    # `set_θ!` is what serves it, and re-setting works as for the complex flavor
    set_θ!(cλ, 0.9)
    reference = sλlmCalculator(0.9, 4, -2)
    @test collect(recurrence!(cλ, 3)) == collect(recurrence!(reference, 3))

    # The complex flavor still takes both
    cY = sYlmCalculator(0.3, 4, -2)
    @test set_R!(cY, R) === cY
    @test set_θ!(cY, 0.3) === cY
end

@testitem "sλlmCalculator: blocks, batches, spin ranges, and half-integers" begin
    import SphericalFunctions: sλlmCalculator, slambdalmCalculator, set_θ!, recurrence!,
        array_view
    import SphericalFunctions: DegreeBlock, DegreeBlockBatch, SpinMatrix, SpinMatrixBatch
    import SphericalFunctions: ℓₘᵢₙ, ℓₘₐₓ, spins, spin, Nᵣ, isbatched

    θ⃗ = [0.3, 0.7, 1.1]
    N = length(θ⃗)

    # The same four shapes of block as an `sYlmCalculator`, holding reals
    one_one   = sλlmCalculator(θ⃗[1], 4, -2)
    many_one  = sλlmCalculator(θ⃗,    4, -2)
    one_many  = sλlmCalculator(θ⃗[1], 4, -2:2)
    many_many = sλlmCalculator(θ⃗,    4, -2:2)
    @test recurrence!(one_one, 3)   isa DegreeBlock
    @test recurrence!(many_one, 3)  isa DegreeBlockBatch
    @test recurrence!(one_many, 3)  isa SpinMatrix
    @test recurrence!(many_many, 3) isa SpinMatrixBatch
    @test all(
        eltype(array_view(recurrence!(c, 3))) === Float64
        for c ∈ (one_one, many_one, one_many, many_many)
    )
    @test axes(recurrence!(many_many, 3)) == (1:N, -2:2, -3:3)
    @test Nᵣ(many_one) == N && !isbatched(one_one) && isbatched(many_one)
    @test spins(one_many) == -2:2 && spin(one_one) == -2

    # Iteration runs over every ℓ from 0, as for the complex flavor
    seen = Int[]
    for (ℓ, block) ∈ one_one
        push!(seen, ℓ)
        @test axes(block) == (-ℓ:ℓ,)
    end
    @test seen == collect(0:4)

    # Half-odd indices
    H = sλlmCalculator(θ⃗[1], 7//2, 1//2)
    @test ℓₘᵢₙ(H) == 1//2 && ℓₘₐₓ(H) == 7//2
    @test axes(recurrence!(H, 3//2)) == (-3//2:3//2,)
    @test eltype(array_view(recurrence!(H, 3//2))) === Float64

    # A calculator can be reused across angles with `set_θ!`, which is the reason it exists,
    # for one angle or for several, and for integer or half-odd indices
    for (calc, θ) ∈ (
        (sλlmCalculator(θ⃗[1], 4, -1:1), θ⃗[3]),
        (sλlmCalculator(θ⃗, 4, -2), reverse(θ⃗)),
        (sλlmCalculator(θ⃗[1], 7//2, -3//2:3//2), θ⃗[2]),
    )
        @test set_θ!(calc, θ) === calc
        fresh = sλlmCalculator(θ, ℓₘₐₓ(calc), calc.s)
        @test all(b == f for ((_, b), (_, f)) ∈ zip(calc, fresh))
    end

    # The ASCII spelling is the same constructor
    @test slambdalmCalculator === sλlmCalculator
end

@testitem "The ring transforms are unchanged by the real tables" begin
    import SphericalFunctions: SSHT, SSHTRS, SSHTMinimal, rotors, Ysize, sλlmCalculator, Nᵣ
    using Random

    # `SSHTRS` and `SSHTMinimal` build their Λ tables with an `sλlmCalculator`.  The precise
    # claim is that the tables are those of the complex harmonics at (θ, 0), read through
    # the relation of the header, so it is asserted with `===` on the table itself; the
    # round trips below are the looser end-to-end confirmation.
    λref(x, s) = isinteger(s) ? real(x) : (mod(Int(2s), 4) == 1 ? imag(x) : -imag(x))

    rng = Random.Xoshiro(20260920)
    for (ℓmax, s) ∈ ((6, 0), (6, -2), (11//2, 1//2), (11//2, -3//2))
        𝒯rs = SSHTRS(s, ℓmax)
        @test 𝒯rs.λ isa sλlmCalculator
        # The table the innermost loops read, against what the complex calculator gave
        cY = sYlmCalculator(collect(𝒯rs.θ), ℓmax, s)
        for ℓ ∈ abs(s):ℓmax
            blkλ = recurrence!(𝒯rs.λ, ℓ)
            blkY = recurrence!(cY, ℓ)
            for iᵣ ∈ 1:Nᵣ(𝒯rs.λ), m ∈ -ℓ:ℓ
                @test blkλ[iᵣ, m] === λref(blkY[iᵣ, m], s)
            end
        end

        # End to end, against the package's synthesis matrix `sYlm_matrix`, which is
        # computed separately from these ring tables.  Measured worst case over
        # the cases here was 11 eps, and 200 is asserted.
        f̃ = randn(rng, ComplexF64, Ysize(abs(s), ℓmax))
        f = 𝒯rs * copy(f̃)
        reference = sYlm_matrix(rotors(𝒯rs), ℓmax, s) * f̃
        @test f ≈ reference atol=200eps(Float64) * maximum(abs, reference)
        # ... and the algorithm still inverts itself: worst case 3 eps, 200 asserted.
        @test 𝒯rs \ copy(f) ≈ f̃ atol=200eps(Float64) * maximum(abs, f̃)

        # `SSHTMinimal` is defined only for integer spin weights.  Against `sYlm_matrix`
        # its synthesis measured 12 eps at worst, and its round trip 11 eps; 200 is asserted
        # for both.  (On the rings of `sorted_rings` the round trip at s = -2 would measure
        # 26_000 eps, which is a property of those rings and not of these tables.)
        if isinteger(s)
            𝒯min = SSHTMinimal(s, ℓmax)
            g = 𝒯min * copy(f̃)
            reference = sYlm_matrix(rotors(𝒯min), ℓmax, s) * f̃
            @test g ≈ reference atol=200eps(Float64) * maximum(abs, reference)
            @test 𝒯min \ copy(g) ≈ f̃ atol=200eps(Float64) * maximum(abs, f̃)
        end
    end
end
