# Tests of `HarmonicValues`, the container `sYlm` returns — `src/mode_weights/containers.jl`.
#
# The point of the container is that one loop reads the same whichever of the four shapes it
# was handed, and that the flat array the transforms want is still one `array_view` call away.

@testitem "HarmonicValues: the four shapes" begin
    using Quaternionic: Rotor
    import SphericalFunctions: DegreeBlock, DegreeBlockBatch, SpinMatrix, SpinMatrixBatch
    import SphericalFunctions: ℓₘᵢₙ, ℓₘₐₓ, spins, spin, Nᵣ, isbatched
    using Random

    rng = Random.Xoshiro(2026)
    ℓmax, N = 5, 3
    R = randn(rng, Rotor{Float64})
    Rs = randn(rng, Rotor{Float64}, N)
    s, sr = -2, -2:2

    one_one  = sYlm(R,  ℓmax, s)     # one rotor,  one spin
    many_one = sYlm(Rs, ℓmax, s)     # many rotors, one spin
    one_many = sYlm(R,  ℓmax, sr)    # one rotor,  a range
    many_many = sYlm(Rs, ℓmax, sr)   # many rotors, a range

    @test one_one[3]   isa DegreeBlock
    @test many_one[3]  isa DegreeBlockBatch
    @test one_many[3]  isa SpinMatrix
    @test many_many[3] isa SpinMatrixBatch

    @test axes(one_one[3])   == (-3:3,)
    @test axes(many_one[3])  == (1:N, -3:3)
    @test axes(one_many[3])  == (-2:2, -3:3)
    @test axes(many_many[3]) == (1:N, -2:2, -3:3)

    @test Nᵣ(one_one) == 1 && Nᵣ(many_one) == N
    @test !isbatched(one_one) && isbatched(many_many)
    @test spins(one_one) == s:s && spin(one_one) == s
    @test spins(one_many) == sr
    @test_throws MethodError spin(one_many)
    @test ℓₘᵢₙ(one_one) == abs(s) && ℓₘₐₓ(one_one) == ℓmax

    # Every shape agrees with the flat storage at the canonical index, which is the property
    # that lets `array_view` be handed to a transform
    for ℓ ∈ abs(s):ℓmax, m ∈ -ℓ:ℓ
        i = Yindex(ℓ, m, abs(s))
        @test one_one[ℓ][m]  == array_view(one_one)[i]
        @test many_one[ℓ][2, m] == array_view(many_one)[2, i]
    end
    for ℓ ∈ 0:ℓmax, m ∈ -ℓ:ℓ, (j, σ) ∈ enumerate(sr)
        i = Yindex(ℓ, m, 0)
        @test one_many[ℓ][σ, m]     == array_view(one_many)[j, i]
        @test many_many[ℓ][2, σ, m] == array_view(many_many)[2, j, i]
    end

    # The batched flat forms are exactly what `sYlm_matrix` gives
    @test array_view(many_one) == sYlm_matrix(Rs, ℓmax, s)
    @test array_view(many_many) == sYlm_matrix(Rs, ℓmax, sr)
    # ... and one rotor's row of the batch is the single-rotor result
    @test array_view(many_one)[2, :] == array_view(sYlm(Rs[2], ℓmax, s))
end

@testitem "HarmonicValues: iteration and blocks are views" begin
    using Quaternionic: Rotor
    using Random

    rng = Random.Xoshiro(7)
    R = randn(rng, Rotor{Float64})
    Y = sYlm(R, 4, -1)

    # Iterating gives `ℓ => block`, as a calculator does, so a loop reads the same either way
    seen = Int[]
    for (ℓ, block) ∈ Y
        push!(seen, ℓ)
        @test block == Y[ℓ]
    end
    @test seen == collect(1:4)
    @test length(Y) == 4              # the number of blocks
    @test collect(keys(Y)) == 1:4
    @test firstindex(Y) == 1 && lastindex(Y) == 4

    # A block is a view, so writing through it writes into the container
    Y[3][0] = 17
    @test array_view(Y)[Yindex(3, 0, 1)] == 17

    # An ℓ the container does not hold says so
    @test_throws ArgumentError Y[1//2]
    @test_throws BoundsError Y[5]
    @test_throws BoundsError Y[0]
end

@testitem "HarmonicValues: half-integer indices" begin
    using Quaternionic: Rotor
    import SphericalFunctions: DegreeBlock, SpinMatrix, ℓₘᵢₙ, ℓₘₐₓ, ishalfinteger
    using Random

    rng = Random.Xoshiro(8)
    R = randn(rng, Rotor{Float64})
    Y = sYlm(R, 7//2, 1//2)

    @test ishalfinteger(Y)
    @test ℓₘᵢₙ(Y) == 1//2 && ℓₘₐₓ(Y) == 7//2
    @test Y[3//2] isa DegreeBlock
    @test axes(Y[3//2]) == (-3//2:3//2,)
    for ℓ ∈ 1//2:7//2, m ∈ -ℓ:ℓ
        @test Y[ℓ][m] == array_view(Y)[Yindex(ℓ, m, 1//2)]
    end
    # A whole number asked of a half-integer container is told what the container holds
    @test_throws ArgumentError Y[2]

    Ys = sYlm(R, 7//2, -3//2:3//2)
    @test Ys[3//2] isa SpinMatrix
    @test axes(Ys[3//2]) == (-3//2:3//2, -3//2:3//2)
end

@testitem "HarmonicValues: sYlm! writes through the container" begin
    using Quaternionic: Rotor
    using Random

    rng = Random.Xoshiro(9)
    R₁ = randn(rng, Rotor{Float64})
    R₂ = randn(rng, Rotor{Float64})

    Y = sYlm(R₁, 4, -2)
    before = copy(array_view(Y))
    sYlm!(Y, R₂, 4, -2)
    @test array_view(Y) == array_view(sYlm(R₂, 4, -2))
    @test array_view(Y) != before
end

# The items above check the four shapes and the indexing.  This one covers the rest of the
# container interface — the type-level queries, the copy/equality pair, and display — none of
# which the transform tests reach.

@testitem "HarmonicValues: type queries, copying, equality and display" begin
    using Quaternionic: Rotor
    import SphericalFunctions: HarmonicValues, AbstractModeContainer
    import SphericalFunctions: ishalfinteger, ℓₘᵢₙ, ℓₘₐₓ, Nᵣ, isbatched, spins, half_integer
    using Random

    rng = Random.Xoshiro(2026)
    R = randn(rng, Rotor{Float64})
    Rs = randn(rng, Rotor{Float64}, 3)
    Y = sYlm(R, 4, -2)

    @test Y isa AbstractModeContainer
    @test eltype(Y) == ComplexF64
    @test eltype(typeof(Y)) == ComplexF64        # the type-level method, used by generic code
    @test ℓₘᵢₙ(Y) == 2 && ℓₘₐₓ(Y) == 4
    @test !ishalfinteger(Y)
    @test Nᵣ(Y) == 1 && !isbatched(Y)
    @test spins(Y) == -2:-2

    # A half-integer container answers `ishalfinteger` the other way
    Yh = sYlm(R, half_integer(7//2), half_integer(1//2))
    @test ishalfinteger(Yh)
    @test ℓₘᵢₙ(Yh) == half_integer(1//2) && ℓₘₐₓ(Yh) == half_integer(7//2)

    # `copy` is independent of the original
    c = copy(Y)
    @test c == Y
    @test c !== Y
    original = Y[3][0]
    c[3][0] = 12345.0 + 0.0im
    @test c[3][0] == 12345.0 + 0.0im
    @test Y[3][0] == original            # writing through the copy does not reach the original
    @test c != Y

    # Equality compares the spin, the ℓ range, the rotor count and the data
    @test sYlm(R, 4, -2) == Y
    @test sYlm(R, 3, -2) != Y                    # a different ℓₘₐₓ
    @test sYlm(R, 4, -1) != Y                    # a different spin
    @test sYlm(Rs, 4, -2) != Y                   # a different rotor count

    # Iteration is by `ℓ => block` pairs, and the container is its own `pairs`
    @test Base.IteratorSize(typeof(Y)) == Base.HasLength()
    @test pairs(Y) === Y
    @test [ℓ for (ℓ, _) ∈ Y] == collect(2:4)

    # `show` names the type, the ℓ range and the spin; the batched form also says how many
    # rotors, and a spin range prints as a range
    s = sprint(show, Y)
    @test occursin("HarmonicValues", s)
    @test occursin("ℓ ∈ 2:4", s)
    @test occursin("s = -2", s)
    @test !occursin("rotors", s)                 # a single rotor is not mentioned

    sb = sprint(show, sYlm(Rs, 4, -2))
    @test occursin("3 rotors", sb)

    sr = sprint(show, sYlm(R, 4, -2:2))
    @test occursin("s ∈ -2:2", sr)

    # The three-argument `show` lists each ℓ block in turn
    s3 = sprint(show, MIME("text/plain"), Y)
    @test occursin("HarmonicValues", s3)
    for ℓ ∈ 2:4
        @test occursin("ℓ = $ℓ", s3)
    end

    # The mode axis must be exactly the length the ℓ range calls for
    @test_throws ArgumentError HarmonicValues(zeros(ComplexF64, 5), -2, 2, 4, 1)
    @test_throws "mode axis has length" HarmonicValues(zeros(ComplexF64, 5), -2, 2, 4, 1)
    @test_throws "Ysize" HarmonicValues(zeros(ComplexF64, 5), -2, 2, 4, 1)
end
