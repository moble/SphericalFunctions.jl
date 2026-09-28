# Tests of the real flavor of the harmonics — `sλlmCalculator`, `sλlm`, `sλlm!` and
# `sλlm_matrix` — which share their struct, their recurrence and their containers with the
# complex `sYlmCalculator` and differ only in the number type they store.
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
# The complex reference must be built in *angle* mode.  Against a rotor built from the same
# θ the values agree only to about 5 ulps, because the rotor path multiplies by a unit phase
# that the angle path never forms — a rounding difference, not a discrepancy, and asserted
# as such at the end of the first item.

@testitem "sλlm is bit-for-bit the real part of the complex harmonics" begin
    import SphericalFunctions: sλlm, sλlmCalculator, ℓₘᵢₙ, ℓₘₐₓ
    using DoubleFloats: Double64
    using Quaternionic: from_spherical_coordinates

    # The relation of the header, written out
    λref(x, s) = isinteger(s) ? real(x) : (mod(Int(2s), 4) == 1 ? imag(x) : -imag(x))

    for T ∈ (Float64, Double64)
        for (ℓmax, spins) ∈ ((4, (0, 1, -2, 2)), (7//2, (1//2, -1//2, 3//2, -3//2)))
            for s ∈ spins
                θ = T(7) / 10
                cY = sYlmCalculator(θ, ℓmax, s)          # complex, angle mode
                cλ = sλlmCalculator(θ, ℓmax, s)
                # The flat form runs the same recurrence to completion, so it agrees with the
                # calculator exactly rather than approximately.
                flat = sλlm(θ, ℓmax, s)
                @test eltype(array_view(flat)) === T
                for ℓ ∈ ℓₘᵢₙ(flat):ℓₘₐₓ(flat)
                    blkY = recurrence!(cY, ℓ)
                    blkλ = recurrence!(cλ, ℓ)
                    for m ∈ -ℓ:ℓ
                        # `===` rather than `==`, so that a signed zero would be caught too.
                        @test blkλ[m] === λref(blkY[m], s)
                        @test flat[ℓ][m] === blkλ[m]
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
        b = array_view(sλlm(θ, ℓmax, s))
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
        @test isempty(cλ.Z₊) && isempty(cλ.Z₋)
        @test !isempty(cY.Z₊) && !isempty(cY.Z₋)
        # The rotor axis, which comes first, is still there, and only the powers are missing
        @test size(cλ.Z₊) == (size(cY.Z₊, 1), 0)
        @test size(cY.Z₊) == (1, 2ℓmax + 1)
        @test number_type(cλ) === Float64
        @test number_type(cY) === ComplexF64
        @test eltype(cλ.Yˡ) === Float64
        # The shared H recurrence is the same size either way — that is what makes this a
        # saving rather than a trade.
        @test size(parent(cλ.H.Hˡ)) == size(parent(cY.H.Hˡ))
    end
end

@testitem "sλlmCalculator refuses rotors" setup=[RefusalChecks] begin
    import SphericalFunctions: sλlmCalculator, sλlm, sλlm_matrix
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
    # ... and the flat forms take an angle, so a rotor is not even a method
    @test_throws MethodError sλlm(R, 4, -2)
    @test_throws MethodError sλlm_matrix([R], 4, -2)

    # `set_θ!` is what serves it, and re-setting works as for the complex flavor
    set_θ!(cλ, 0.9)
    reference = sλlmCalculator(0.9, 4, -2)
    @test collect(recurrence!(cλ, 3)) == collect(recurrence!(reference, 3))

    # The complex flavor still takes both
    cY = sYlmCalculator(0.3, 4, -2)
    @test set_R!(cY, R) === cY
    @test set_θ!(cY, 0.3) === cY
end

@testitem "sλlm: containers, batches, spin ranges and half-integers" setup=[RefusalChecks] begin
    import SphericalFunctions: sλlm, sλlm!, sλlm_matrix, sλlmCalculator, slambdalm,
        slambdalm!, slambdalm_matrix, slambdalmCalculator
    import SphericalFunctions: DegreeBlock, DegreeBlockBatch, SpinMatrix, SpinMatrixBatch
    import SphericalFunctions: ℓₘᵢₙ, ℓₘₐₓ, spins, spin, Nᵣ, isbatched, Yindex

    θ⃗ = [0.3, 0.7, 1.1]
    N = length(θ⃗)

    # The same four shapes as `sYlm`, holding reals
    one_one   = sλlm(θ⃗[1], 4, -2)
    many_one  = sλlm(θ⃗,    4, -2)
    one_many  = sλlm(θ⃗[1], 4, -2:2)
    many_many = sλlm(θ⃗,    4, -2:2)
    @test one_one[3]   isa DegreeBlock
    @test many_one[3]  isa DegreeBlockBatch
    @test one_many[3]  isa SpinMatrix
    @test many_many[3] isa SpinMatrixBatch
    @test all(eltype(Y[3]) === Float64 for Y ∈ (one_one, many_one, one_many, many_many))
    @test axes(many_many[3]) == (1:N, -2:2, -3:3)
    @test Nᵣ(many_one) == N && !isbatched(one_one)
    @test spins(one_many) == -2:2 && spin(one_one) == -2

    # The flat batched form is exactly `sλlm_matrix`, as for the complex pair
    @test array_view(many_one) == sλlm_matrix(θ⃗, 4, -2)
    @test array_view(many_many) == sλlm_matrix(θ⃗, 4, -2:2)
    @test array_view(many_one)[2, :] == array_view(sλlm(θ⃗[2], 4, -2))

    # Iteration reads the same as a calculator's
    seen = Int[]
    for (ℓ, block) ∈ one_one
        push!(seen, ℓ)
        @test block == one_one[ℓ]
    end
    @test seen == collect(2:4)

    # Half-odd indices, including ℓₘᵢₙ = 1/2
    H = sλlm(θ⃗[1], 7//2, 1//2)
    @test ℓₘᵢₙ(H) == 1//2 && ℓₘₐₓ(H) == 7//2
    @test eltype(array_view(H)) === Float64
    for ℓ ∈ 1//2:7//2, m ∈ -ℓ:ℓ
        @test H[ℓ][m] == array_view(H)[Yindex(ℓ, m, 1//2)]
    end

    # `sλlm!` writes into existing storage, including through the container
    Y = zeros(Float64, length(array_view(one_one)))
    sλlm!(Y, θ⃗[1], 4, -2)
    @test Y == array_view(one_one)
    container = sλlm(θ⃗[2], 4, -2)
    sλlm!(container, θ⃗[1], 4, -2)
    @test array_view(container) == array_view(one_one)
    # The element type must match the calculator's, and says so
    @test refuses(
        () -> sλlm!(zeros(Float32, length(Y)), θ⃗[1], 4, -2), ArgumentError,
        "element type must be Float64"
    )
    @test_throws MethodError sλlm!(zeros(ComplexF64, length(Y)), θ⃗[1], 4, -2)

    # A calculator can be reused across angles, which is the reason it exists, in each of
    # its forms: for everything it serves, for one spin weight of a range, and into a
    # container
    calc = sλlmCalculator(θ⃗[1], 4, -2)
    Y2 = similar(Y)
    sλlm!(Y2, calc, θ⃗[3])
    @test Y2 == array_view(sλlm(θ⃗[3], 4, -2))
    cr = sλlmCalculator(θ⃗[1], 4, -1:1)
    v = zeros(Ysize(1, 4))
    @test array_view(sλlm!(v, cr, 0.7, 1)) == array_view(sλlm(0.7, 4, 1))
    @test parent(array_view(sλlm!(v, cr, 0.7, 1))) === v
    @test array_view(sλlm!(v, cr, 0.9, 1; ell_min=1)) == array_view(sλlm(0.9, 4, 1))
    Λr = sλlm(0.1, 4, -1:1)
    @test sλlm!(Λr, cr, 0.7) === Λr
    @test array_view(Λr) == array_view(sλlm(0.7, 4, -1:1))
    Λ1 = sλlm(0.1, 4, 1)
    @test sλlm!(Λ1, cr, 0.9, 1) === Λ1
    @test array_view(Λ1) == array_view(sλlm(0.9, 4, 1))
    # ... and for half-odd indices, with ℓₘᵢₙ spelled as a `Rational`
    crₕ = sλlmCalculator(θ⃗[1], 7//2, -3//2:3//2)
    vₕ = zeros(Ysize(1//2, 7//2))
    @test array_view(sλlm!(vₕ, crₕ, 0.7, 3//2; ℓₘᵢₙ=1//2)) ==
        array_view(sλlm(0.7, 7//2, 3//2; ℓₘᵢₙ=1//2))
    # A calculator built from a vector, even of one angle, cannot fill one angle's values
    @test refuses(
        () -> sλlm!(zeros(Ysize(1, 4)), sλlmCalculator([0.1], 4, 1), 0.7), ArgumentError,
        "`sλlm!` needs a calculator built for a single rotor"
    )

    # The ASCII spellings are the same functions, and the keyword has one too
    @test slambdalm === sλlm && slambdalm! === sλlm! && slambdalm_matrix === sλlm_matrix
    @test slambdalmCalculator === sλlmCalculator
    @test array_view(slambdalm(0.7, 4, 1; ell_min=2)) == array_view(sλlm(0.7, 4, 1; ℓₘᵢₙ=2))
end

@testitem "The ring transforms are unchanged by the real tables" begin
    import SphericalFunctions: SSHT, SSHTRS, SSHTMinimal, rotors, Ysize, sλlmCalculator
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
            for iᵣ ∈ 1:size(𝒯rs.λ.Yˡ, 1), m ∈ -ℓ:ℓ
                @test blkλ[iᵣ, m] === λref(blkY[iᵣ, m], s)
            end
        end

        # End to end, against the closed-form synthesis matrix.  Measured worst case over
        # the cases here was 11 eps, and 200 is asserted.
        f̃ = randn(rng, ComplexF64, Ysize(abs(s), ℓmax))
        f = 𝒯rs * copy(f̃)
        reference = sYlm_matrix(rotors(𝒯rs), ℓmax, s) * f̃
        @test f ≈ reference atol=200eps(Float64) * maximum(abs, reference)
        # ... and the algorithm still inverts itself: worst case 3 eps, 200 asserted.
        @test 𝒯rs \ copy(f) ≈ f̃ atol=200eps(Float64) * maximum(abs, f̃)

        # `SSHTMinimal` is defined only for integer spin weights.  Against the closed form
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
