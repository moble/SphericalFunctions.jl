# Tests of `HarmonicValues`, the container `sYlm` returns — `src/mode_weights/containers.jl`.
#
# The point of the container is that one loop reads the same whichever of the four shapes it
# was handed, and that the flat array the transforms want is still one `strided` call away.

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
    # that lets `strided` be handed to a transform
    for ℓ ∈ abs(s):ℓmax, m ∈ -ℓ:ℓ
        i = Yindex(ℓ, m, abs(s))
        @test one_one[ℓ][m]  == strided(one_one)[i]
        @test many_one[ℓ][2, m] == strided(many_one)[2, i]
    end
    for ℓ ∈ 0:ℓmax, m ∈ -ℓ:ℓ, (j, σ) ∈ enumerate(sr)
        i = Yindex(ℓ, m, 0)
        @test one_many[ℓ][σ, m]     == strided(one_many)[j, i]
        @test many_many[ℓ][2, σ, m] == strided(many_many)[2, j, i]
    end

    # The batched flat forms are exactly what `sYlm_matrix` gives
    @test strided(many_one) == sYlm_matrix(Rs, ℓmax, s)
    @test strided(many_many) == sYlm_matrix(Rs, ℓmax, sr)
    # ... and one rotor's row of the batch is the single-rotor result
    @test strided(many_one)[2, :] == strided(sYlm(Rs[2], ℓmax, s))
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
    @test strided(Y)[Yindex(3, 0, 1)] == 17

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
        @test Y[ℓ][m] == strided(Y)[Yindex(ℓ, m, 1//2)]
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
    before = copy(strided(Y))
    sYlm!(Y, R₂, 4, -2)
    @test strided(Y) == strided(sYlm(R₂, 4, -2))
    @test strided(Y) != before
end
