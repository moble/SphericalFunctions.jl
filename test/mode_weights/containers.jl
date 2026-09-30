# Tests of `HarmonicValues`, the container `sYlm` returns —
# `src/containers/series.jl`.
#
# The point of the container is that one loop reads the same whichever of the four shapes it
# was handed, and that the flat array the transforms use can be obtained with `array_view`.

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

    # Iterating gives `ℓ => block`, as a calculator does, so a loop reads the same either
    # way
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

# The items above check the four shapes and the indexing.  This one covers the rest of the
# container interface — the type-level queries, the copy/equality pair, and display — none
# of which the transform tests reach.

@testitem "HarmonicValues: type queries, copying, equality and display" begin
    using Quaternionic: Rotor
    import SphericalFunctions: HarmonicValues, AbstractDegreeSeries
    import SphericalFunctions: ishalfinteger, ℓₘᵢₙ, ℓₘₐₓ, Nᵣ, isbatched, spins, HalfOddInteger
    import SphericalFunctions: array_view, relabel
    using Random

    rng = Random.Xoshiro(2026)
    R = randn(rng, Rotor{Float64})
    Rs = randn(rng, Rotor{Float64}, 3)
    Y = sYlm(R, 4, -2)

    @test Y isa AbstractDegreeSeries{Int}
    # The element type is that of the iteration, `ℓ => block`, as for a calculator, so that
    # `collect` works; the number type is that of the storage
    @test eltype(Y) === eltype(typeof(Y)) === typeof(first(iterate(Y)))
    @test eltype(array_view(Y)) === ComplexF64
    @test ℓₘᵢₙ(Y) == 2 && ℓₘₐₓ(Y) == 4
    @test !ishalfinteger(Y)
    @test Nᵣ(Y) == 1 && !isbatched(Y)
    @test spins(Y) == -2:-2

    # A half-integer container answers `ishalfinteger` the other way
    Yh = sYlm(R, HalfOddInteger(7//2), HalfOddInteger(1//2))
    @test ishalfinteger(Yh)
    @test ℓₘᵢₙ(Yh) == HalfOddInteger(1//2) && ℓₘₐₓ(Yh) == HalfOddInteger(7//2)

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

    # Against a plain array, `==`, `isequal` and `≈` compare the numbers of `array_view`, in
    # either order; between two containers they compare the labels as well, as for
    # `ModeWeights`, so the same numbers labelled with the opposite spin are neither equal
    # nor approximately equal
    @test Y == array_view(Y) && array_view(Y) == Y
    @test isequal(Y, array_view(Y)) && isequal(array_view(Y), Y)
    @test Y ≈ array_view(Y) && array_view(Y) ≈ Y
    @test Y ≈ array_view(Y) .+ 1e-14 && !(Y ≈ array_view(Y) .+ 1e-3)
    @test Y ≈ copy(Y) && isequal(Y, copy(Y))
    relabelled = HarmonicValues(copy(array_view(Y)), 2, 2, 4, 1)
    @test relabelled == array_view(Y)
    @test relabelled != Y && !(relabelled ≈ Y) && !isequal(relabelled, Y)

    # Broadcasting reads the numbers of `array_view`, and gives a plain array; a broadcast
    # that combines two containers, or writes one into another, pairs their numbers by
    # position, and so requires their labels to agree, while a plain array is paired as it
    # is
    @test Y .+ 1 isa Vector{ComplexF64} && Y .+ 1 == array_view(Y) .+ 1
    @test 2 .* Y == 2 .* array_view(Y) && Y .- copy(Y) == zeros(21)
    @test (similar(Y) .= Y) == Y
    labels = "must have the same labels"
    @test_throws ArgumentError Y .+ relabelled
    @test_throws labels Y .+ relabelled
    @test_throws labels similar(Y) .= relabelled
    @test_throws labels similar(Y) .= 2 .* relabelled .+ Y
    @test Y .+ array_view(relabelled) == 2 .* array_view(Y)
    @test (similar(Y) .= array_view(relabelled)) == array_view(Y)

    # `similar` keeps the labels and the shape, optionally with another number type, and
    # `.=` writes into the storage, for every shape
    Ys = (
        Y, sYlm(Rs, 4, -2), sYlm(R, HalfOddInteger(7//2), HalfOddInteger(1//2)),
        HarmonicValues(randn(rng, ComplexF64, 3, 24), -1:1, 1, 4, 1),
        HarmonicValues(randn(rng, ComplexF64, 2, 3, 24), -1:1, 1, 4, 2),
    )
    for X ∈ Ys
        Z = similar(X)
        @test Z isa typeof(X) && size(array_view(Z)) == size(array_view(X))
        @test (spins(Z), ℓₘᵢₙ(Z), ℓₘₐₓ(Z), Nᵣ(Z)) == (spins(X), ℓₘᵢₙ(X), ℓₘₐₓ(X), Nᵣ(X))
        @test array_view(Z) !== array_view(X)
        Z32 = similar(X, ComplexF32)
        @test eltype(array_view(Z32)) === ComplexF32 && ℓₘₐₓ(Z32) == ℓₘₐₓ(X)
        @test (Z .= 2 .* array_view(X)) === Z
        @test Z == 2 .* array_view(X)
        # `relabel` wraps an array of the same shape without copying it, and refuses another
        A = 3 .* array_view(X)
        W = relabel(X, A)
        @test array_view(W) === A && W ≈ 3 .* array_view(X)
        @test (spins(W), ℓₘᵢₙ(W), ℓₘₐₓ(W), Nᵣ(W)) == (spins(X), ℓₘᵢₙ(X), ℓₘₐₓ(X), Nᵣ(X))
        @test_throws DimensionMismatch relabel(X, zeros(ComplexF64, 2, size(array_view(X))...))
        # `collect` gives the pairs, as a comprehension over the container does
        @test collect(X) == [ℓ => X[ℓ] for ℓ ∈ keys(X)]
        @test eltype(collect(X)) === eltype(X)
    end

    # Iteration is by `ℓ => block` pairs, and the container is its own `pairs`
    @test Base.IteratorSize(typeof(Y)) == Base.HasLength()
    @test pairs(Y) === Y
    @test [ℓ for (ℓ, _) ∈ Y] == collect(2:4)

    # Indexing, `Y[ℓ, :]`, `first` and `last`, with or without a count, and `only` give
    # blocks, as for a `WignerSeries`, whatever the kind of index
    @test Y[3, :] == Y[3] && Y[3, :] isa typeof(Y[3])
    @test first(Y) == Y[2] && last(Y) == Y[4]
    @test first(Y) isa typeof(Y[2]) && last(Y) isa typeof(Y[4])
    @test first(Y, 2) == [Y[2], Y[3]] && last(Y, 2) == [Y[3], Y[4]]
    @test first(Y, 10) == [Y[2], Y[3], Y[4]] && last(Y, 10) == first(Y, 3)
    @test isempty(first(Y, 0)) && isempty(last(Y, 0))
    @test_throws ArgumentError first(Y, -1)
    @test only(sYlm(R, 2, 2)) == sYlm(R, 2, 2)[2]
    @test_throws ArgumentError only(Y)
    @test_throws "hold 3 blocks" only(Y)
    h = HalfOddInteger
    @test first(Yh) == Yh[1//2] && last(Yh) == Yh[7//2] && Yh[h(3//2), :] == Yh[3//2]
    @test first(Yh, 2) == [Yh[1//2], Yh[3//2]] && last(Yh, 1) == [Yh[7//2]]
    @test first(sYlm(Rs, 4, -2)) == sYlm(Rs, 4, -2)[2]
    @test last(sYlm(R, 4, -2:2), 1) == [sYlm(R, 4, -2:2)[4]]

    # `hash` agrees with `isequal`, which, against a plain array, compares the numbers alone,
    # so the hash is that of the numbers; `isequal` between two containers counts the labels,
    # and follows `isequal` of the numbers for NaN and signed zeros
    @test hash(Y) == hash(copy(Y)) == hash(array_view(Y))
    @test length(Set([Y, copy(Y)])) == 1 && length(unique([Y, copy(Y), sYlm(R, 4, -1)])) == 2
    Yn = copy(Y)
    array_view(Yn)[1] = NaN
    @test isequal(Yn, Yn) && isequal(Yn, copy(Yn)) && Yn != Yn
    @test hash(Yn) == hash(copy(Yn))
    Yz, Yz′ = copy(Y), copy(Y)
    array_view(Yz)[1], array_view(Yz′)[1] = 0.0, -0.0
    @test Yz == Yz′ && !isequal(Yz, Yz′)

    # `show` names the type, the ℓ range and the spin; the batched form also says how many
    # rotors, even for a batch of one, and a spin range prints as a range
    s = sprint(show, Y)
    @test s == "HarmonicValues{ComplexF64} for ℓ ∈ 2:4, s=-2"
    @test !occursin("rotor", s)                  # a single rotor is not mentioned

    sb = sprint(show, sYlm(Rs, 4, -2))
    @test occursin("3 rotors", sb)
    @test endswith(sprint(show, sYlm([R], 4, -2)), ", 1 rotor")

    sr = sprint(show, sYlm(R, 4, -2:2))
    @test occursin("s ∈ -2:2", sr)

    # The three-argument `show` lists each ℓ block in turn, but where the output is limited,
    # as at the REPL, it lists only the first two and the last two of more than four
    s3 = sprint(show, MIME("text/plain"), Y)
    @test occursin("HarmonicValues", s3)
    for ℓ ∈ 2:4
        @test occursin("ℓ = $ℓ", s3)
    end
    Y200 = sYlm(R, 200, 0)
    full = sprint(show, MIME("text/plain"), Y200)
    limited = sprint(
        show, MIME("text/plain"), Y200; context=(:limit => true, :displaysize => (24, 80))
    )
    @test all(ℓ -> occursin(" ℓ = $ℓ:", full), 0:200)
    @test all(ℓ -> occursin(" ℓ = $ℓ:", limited), (0, 1, 199, 200))
    @test !occursin(" ℓ = 2:", limited) && occursin("⋮", limited)
    @test count('\n', limited) < 100 < count('\n', full)

    # The mode axis must be exactly the length the ℓ range calls for
    @test_throws ArgumentError HarmonicValues(zeros(ComplexF64, 5), -2, 2, 4, 1)
    @test_throws "mode axis has length" HarmonicValues(zeros(ComplexF64, 5), -2, 2, 4, 1)
    @test_throws "Ysize" HarmonicValues(zeros(ComplexF64, 5), -2, 2, 4, 1)
end

@testitem "HarmonicValues: the labels must describe the storage" begin
    import SphericalFunctions: HarmonicValues, Ysize, spins, Nᵣ, isbatched, HalfOddInteger
    using Random

    rng = Random.Xoshiro(2027)
    n = Ysize(2, 4)

    # The four shapes, each with its labels, are accepted, and report what was given
    for (dims, s, N) in (((n,), 2, 1), ((3, n), 2, 3), ((5, n), -2:2, 1), ((3, 5, n), -2:2, 3), ((1, n), 2, 1))
        Y = HarmonicValues(randn(rng, ComplexF64, dims...), s, 2, 4, N)
        @test Nᵣ(Y) == N && spins(Y) == (s isa Integer ? (s:s) : s)
        @test isbatched(Y) == (length(dims) == (s isa Integer ? 2 : 3))
    end

    # A rotor axis that `Nᵣ` does not describe, in either direction, is refused ...
    @test_throws DimensionMismatch HarmonicValues(randn(rng, ComplexF64, 3, n), 2, 2, 4, 1)
    @test_throws "rotor axis of length 3, but Nᵣ=1" HarmonicValues(randn(rng, ComplexF64, 3, n), 2, 2, 4, 1)
    @test_throws "rotor axis of length 1, but Nᵣ=2" HarmonicValues(randn(rng, ComplexF64, 1, n), 2, 2, 4, 2)
    @test_throws "no rotor axis" HarmonicValues(randn(rng, ComplexF64, n), 2, 2, 4, 3)
    @test_throws "no rotor axis" HarmonicValues(randn(rng, ComplexF64, 5, n), -2:2, 2, 4, 2)
    # ... as are a spin axis of the wrong length, an empty range of spin weights, and a rank
    # that fits none of the shapes
    @test_throws "spin axis of the storage has length 3" HarmonicValues(randn(rng, ComplexF64, 3, n), -2:2, 2, 4, 1)
    @test_throws "spin axis of the storage has length 5" HarmonicValues(randn(rng, ComplexF64, 2, 5, n), 1:3, 2, 4, 2)
    @test_throws "is empty" HarmonicValues(randn(rng, ComplexF64, 0, n), 2:1, 2, 4, 1)
    @test_throws "has 3 dimensions" HarmonicValues(randn(rng, ComplexF64, 1, 1, n), 2, 2, 4, 1)
    @test_throws "has 1 dimensions" HarmonicValues(randn(rng, ComplexF64, n), -2:2, 2, 4, 1)
    @test_throws DimensionMismatch HarmonicValues(randn(rng, ComplexF64, 2, 2, 5, n), -2:2, 2, 4, 2)

    # The indices are all integers of type `Int` or all half-odd-integers, the spin weights
    # included, so a half-integer spin weight is refused for an integer ℓ range; a `Rational`
    # is converted, alone or as the endpoints of a range
    @test_throws "and so mixes integers (ℓₘᵢₙ, ℓₘₐₓ) with half-odd-integers (s)" HarmonicValues(randn(rng, ComplexF64, n), 1//2, 2, 4, 1)
    @test_throws "`Int32` is narrower than `Int`" HarmonicValues(randn(rng, ComplexF64, n), Int32(2), 2, 4, 1)
    @test_throws "must be a unit range" HarmonicValues(randn(rng, ComplexF64, 3, n), 1:1:3, 2, 4, 1)
    m = Ysize(1//2, 7//2)
    Yh = HarmonicValues(randn(rng, ComplexF64, 2, m), -1//2:1//2, 1//2, 7//2, 1)
    @test spins(Yh) isa UnitRange{HalfOddInteger} && spins(Yh) == -1//2:1//2
    @test HarmonicValues(randn(rng, ComplexF64, m), 3//2, 1//2, 7//2, 1) isa
        HarmonicValues{ComplexF64, HalfOddInteger, HalfOddInteger}
end
