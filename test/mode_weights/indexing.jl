# Tests for the closed-form indexing of the canonical mode-weight ordering
#
#     [ f(ℓ, m) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ ],
#
# provided by `Ysize`, `Yindex` and `Yrange`.  The oracle is that ordering itself, enumerated
# directly by the local helper `ℓmpairs(ℓₘᵢₙ, ℓₘₐₓ)`: `Ysize` must be its length, `Yrange`
# must be the list, and `Yindex(ℓ, m, ℓₘᵢₙ)` must be the 1-based position of `(ℓ, m)` in it.
# The enumeration shares no code with the closed forms under test, so it is an independent
# reference; the structural invariants of the ordering (first pair at index 1, consecutive
# pairs at consecutive indices, last pair at index `Ysize`) are checked as well.
#
# The item names use an "Indexing: " prefix so that they are distinguishable — to a human
# reading a results list, and to `runtests.jl`'s `occursin` filters — from the v2 items,
# which were also named `Ysize`, `Yindex` and `Yrange`.  Those items are gone with the
# `Deprecated` module, but the prefix is kept: the names are what people filter on.

@testitem "Indexing: Ysize" begin
    import SphericalFunctions: Ysize

    # The ordering, built by brute force: `Ysize` is supposed to count these pairs
    ℓmpairs(ℓₘᵢₙ, ℓₘₐₓ) = [(ℓ, m) for ℓ in ℓₘᵢₙ:ℓₘₐₓ for m in -ℓ:ℓ]

    # Agreement with the enumeration over every valid pair, including the empty range
    # ℓₘₐₓ = ℓₘᵢₙ - 1, whose size is 0
    for ℓₘᵢₙ in 0:6, ℓₘₐₓ in ℓₘᵢₙ-1:12
        n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
        @test n == length(ℓmpairs(ℓₘᵢₙ, ℓₘₐₓ))
        @test n == (ℓₘₐₓ + 1)^2 - ℓₘᵢₙ^2
        # An independent count of the pairs in the ordering
        @test n == sum(2ℓ + 1 for ℓ in ℓₘᵢₙ:ℓₘₐₓ; init=0)
        @test n ≥ 0
    end
    for ℓₘᵢₙ in 0:6
        @test Ysize(ℓₘᵢₙ, ℓₘᵢₙ - 1) == 0
        @test Ysize(ℓₘᵢₙ, ℓₘᵢₙ) == 2ℓₘᵢₙ + 1
    end

    # The one-argument form is ℓₘᵢₙ = 0
    for ℓₘₐₓ in -1:12
        @test Ysize(ℓₘₐₓ) == Ysize(0, ℓₘₐₓ)
        @test Ysize(ℓₘₐₓ) == (ℓₘₐₓ + 1)^2
        @test Ysize(ℓₘₐₓ) == length(ℓmpairs(0, ℓₘₐₓ))
    end

    # A size can never be negative (GitHub issue #52).  The bare formula (ℓₘₐₓ+1)² - ℓₘᵢₙ²,
    # which is what v2 returned, goes negative whenever ℓₘᵢₙ > ℓₘₐₓ + 1, even though the
    # ordering it is meant to count is simply empty there; v3 throws instead ...
    @test (0 + 1)^2 - 3^2 < 0
    @test isempty(ℓmpairs(3, 0))
    @test_throws ArgumentError Ysize(3, 0)
    for ℓₘᵢₙ in 0:6, ℓₘₐₓ in -4:ℓₘᵢₙ-2
        @test isempty(ℓmpairs(ℓₘᵢₙ, ℓₘₐₓ))
        @test_throws ArgumentError Ysize(ℓₘᵢₙ, ℓₘₐₓ)
    end
    for ℓₘₐₓ in -4:-2
        @test_throws ArgumentError Ysize(ℓₘₐₓ)
    end
    # ... and also rejects a negative ℓₘᵢₙ, even where the formula would give a positive number
    for ℓₘᵢₙ in -3:-1, ℓₘₐₓ in -1:6
        @test_throws ArgumentError Ysize(ℓₘᵢₙ, ℓₘₐₓ)
    end
    # The error message names the offending argument
    @test_throws "ℓₘₐₓ" Ysize(3, 0)
    @test_throws "ℓₘᵢₙ" Ysize(-1, 3)

    # Narrower integer types, alone or mixed, give the same numbers; (10+1)² = 121 fits in Int8
    for IT in (Int8, Int32)
        for ℓₘᵢₙ in 0:6, ℓₘₐₓ in ℓₘᵢₙ-1:10
            @test Ysize(IT(ℓₘᵢₙ), IT(ℓₘₐₓ)) == length(ℓmpairs(ℓₘᵢₙ, ℓₘₐₓ))
            @test Ysize(IT(ℓₘᵢₙ), IT(ℓₘₐₓ)) == Ysize(ℓₘᵢₙ, ℓₘₐₓ)
            @test Ysize(IT(ℓₘᵢₙ), IT(ℓₘₐₓ)) isa Integer
            @test Ysize(IT(ℓₘᵢₙ), ℓₘₐₓ) == Ysize(ℓₘᵢₙ, ℓₘₐₓ)
            @test Ysize(ℓₘᵢₙ, IT(ℓₘₐₓ)) == Ysize(ℓₘᵢₙ, ℓₘₐₓ)
        end
        for ℓₘₐₓ in -1:10
            @test Ysize(IT(ℓₘₐₓ)) == Ysize(ℓₘₐₓ)
        end
        @test_throws ArgumentError Ysize(IT(3), IT(1))
        @test_throws ArgumentError Ysize(IT(-1), IT(3))
        @test_throws ArgumentError Ysize(IT(-2))
    end
end


@testitem "Indexing: Yindex" begin
    import SphericalFunctions: Yindex, Ysize

    # The ordering, built by brute force: `Yindex(ℓ, m, ℓₘᵢₙ)` is supposed to be the 1-based
    # position of `(ℓ, m)` in this list
    ℓmpairs(ℓₘᵢₙ, ℓₘₐₓ) = [(ℓ, m) for ℓ in ℓₘᵢₙ:ℓₘₐₓ for m in -ℓ:ℓ]
    positions(ℓₘᵢₙ, ℓₘₐₓ) = Dict(p => i for (i, p) in enumerate(ℓmpairs(ℓₘᵢₙ, ℓₘₐₓ)))

    ℓₘₐₓ = 12

    # Agreement with the enumerated position for every (ℓ, m, ℓₘᵢₙ) with ℓₘᵢₙ ≤ ℓ ≤ ℓₘₐₓ
    for ℓₘᵢₙ in 0:ℓₘₐₓ
        position = positions(ℓₘᵢₙ, ℓₘₐₓ)
        for ℓ in ℓₘᵢₙ:ℓₘₐₓ, m in -ℓ:ℓ
            i = Yindex(ℓ, m, ℓₘᵢₙ)
            @test i == position[(ℓ, m)]
            @test i == ℓ*(ℓ+1) - ℓₘᵢₙ^2 + m + 1
            @test i isa Integer
        end
    end
    # ℓₘᵢₙ defaults to 0
    let position = positions(0, ℓₘₐₓ)
        for ℓ in 0:ℓₘₐₓ, m in -ℓ:ℓ
            @test Yindex(ℓ, m) == Yindex(ℓ, m, 0)
            @test Yindex(ℓ, m) == position[(ℓ, m)]
        end
    end

    # The ordering is ℓ-major with m increasing: the first pair (ℓₘᵢₙ, -ℓₘᵢₙ) is at index 1,
    # stepping m by one steps the index by one, the step from (ℓ, ℓ) to (ℓ+1, -ℓ-1) is also
    # one, and the last pair (ℓₘₐₓ, ℓₘₐₓ) is at index Ysize(ℓₘᵢₙ, ℓₘₐₓ)
    for ℓₘᵢₙ in 0:ℓₘₐₓ
        @test Yindex(ℓₘᵢₙ, -ℓₘᵢₙ, ℓₘᵢₙ) == 1
        for ℓ in ℓₘᵢₙ:ℓₘₐₓ
            for m in -ℓ:ℓ-1
                @test Yindex(ℓ, m+1, ℓₘᵢₙ) == Yindex(ℓ, m, ℓₘᵢₙ) + 1
            end
            @test Yindex(ℓ+1, -ℓ-1, ℓₘᵢₙ) == Yindex(ℓ, ℓ, ℓₘᵢₙ) + 1
            # Every ℓ in the range serves as an ℓₘₐₓ here: the last index of the ℓ block is
            # the size up to and including ℓ, and its first index follows the size up to ℓ-1
            @test Yindex(ℓ, ℓ, ℓₘᵢₙ) == Ysize(ℓₘᵢₙ, ℓ)
            @test Yindex(ℓ, -ℓ, ℓₘᵢₙ) == Ysize(ℓₘᵢₙ, ℓ - 1) + 1
        end
        # Altogether, the indices enumerate 1:Ysize exactly once each, in order
        indices = [Yindex(ℓ, m, ℓₘᵢₙ) for ℓ in ℓₘᵢₙ:ℓₘₐₓ for m in -ℓ:ℓ]
        @test indices == 1:Ysize(ℓₘᵢₙ, ℓₘₐₓ)
        @test indices == 1:length(ℓmpairs(ℓₘᵢₙ, ℓₘₐₓ))
    end

    # Narrower integer types, alone or mixed; ℓ ≤ 10 keeps every index within Int8
    for IT in (Int8, Int32)
        for ℓₘᵢₙ in 0:4
            position = positions(ℓₘᵢₙ, 10)
            for ℓ in ℓₘᵢₙ:10, m in -ℓ:ℓ
                @test Yindex(IT(ℓ), IT(m), IT(ℓₘᵢₙ)) == position[(ℓ, m)]
                @test Yindex(IT(ℓ), IT(m), IT(ℓₘᵢₙ)) == Yindex(ℓ, m, ℓₘᵢₙ)
                @test Yindex(IT(ℓ), IT(m)) == Yindex(ℓ, m)
                @test Yindex(IT(ℓ), m, ℓₘᵢₙ) == Yindex(ℓ, m, ℓₘᵢₙ)
            end
        end
    end
end


@testitem "Indexing: Yrange" begin
    import SphericalFunctions: Yrange, Yindex, Ysize

    # The ordering, built by brute force: `Yrange` is supposed to return exactly this list
    ℓmpairs(ℓₘᵢₙ, ℓₘₐₓ) = [(ℓ, m) for ℓ in ℓₘᵢₙ:ℓₘₐₓ for m in -ℓ:ℓ]

    for ℓₘᵢₙ in 0:7, ℓₘₐₓ in ℓₘᵢₙ-1:12
        r = Yrange(ℓₘᵢₙ, ℓₘₐₓ)
        @test r isa AbstractVector
        @test eltype(r) <: Tuple{Integer, Integer}
        @test length(r) == Ysize(ℓₘᵢₙ, ℓₘₐₓ)
        # The list itself, pair by pair
        @test r == ℓmpairs(ℓₘᵢₙ, ℓₘₐₓ)
        @test length(r) == length(ℓmpairs(ℓₘᵢₙ, ℓₘₐₓ))
        # `Yrange` and `Yindex` are inverse to each other
        for ℓ in ℓₘᵢₙ:ℓₘₐₓ, m in -ℓ:ℓ
            @test r[Yindex(ℓ, m, ℓₘᵢₙ)] == (ℓ, m)
        end
        for (i, (ℓ, m)) in enumerate(r)
            @test Yindex(ℓ, m, ℓₘᵢₙ) == i
            @test ℓₘᵢₙ ≤ ℓ ≤ ℓₘₐₓ
            @test -ℓ ≤ m ≤ ℓ
        end
        # ℓ-major, m increasing, no repeats
        @test issorted(r)
        @test allunique(r)
        if !isempty(r)
            @test first(r) == (ℓₘᵢₙ, -ℓₘᵢₙ)
            @test last(r) == (ℓₘₐₓ, ℓₘₐₓ)
        end
    end

    # The one-argument form is ℓₘᵢₙ = 0
    for ℓₘₐₓ in -1:12
        @test Yrange(ℓₘₐₓ) == Yrange(0, ℓₘₐₓ)
        @test Yrange(ℓₘₐₓ) == ℓmpairs(0, ℓₘₐₓ)
    end

    # Narrower integer types
    for IT in (Int8, Int32)
        for ℓₘᵢₙ in 0:3, ℓₘₐₓ in ℓₘᵢₙ-1:10
            r = Yrange(IT(ℓₘᵢₙ), IT(ℓₘₐₓ))
            @test r == Yrange(ℓₘᵢₙ, ℓₘₐₓ)
            @test r == ℓmpairs(ℓₘᵢₙ, ℓₘₐₓ)
            @test length(r) == Ysize(IT(ℓₘᵢₙ), IT(ℓₘₐₓ))
        end
        @test Yrange(IT(5)) == Yrange(5)
    end
end


# The half-integer items below repeat the checks above for half-odd indices.  The oracle is
# again the ordering itself, but enumerated in `Rational` arithmetic — a `UnitRange` of
# `Rational`s steps by one — which shares nothing with the `Int` closed forms under test.
# Half-odd-integers are spelled as `Rational`s at the public boundary, as a user would write
# them, and the `HalfOddInteger` spelling is checked against that.

@testitem "Indexing: half-integer Ysize" begin
    import SphericalFunctions: Ysize, HalfOddInteger

    # The ordering, built by brute force in `Rational`s
    ℓmpairs(ℓₘᵢₙ, ℓₘₐₓ) = [(ℓ, m) for ℓ in ℓₘᵢₙ:ℓₘₐₓ for m in -ℓ:ℓ]

    # Agreement with the enumeration over every valid pair, including the empty range
    # ℓₘₐₓ = ℓₘᵢₙ - 1, whose size is 0
    for ℓₘᵢₙ in 1//2:11//2, ℓₘₐₓ in ℓₘᵢₙ-1:25//2
        n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
        @test n == length(ℓmpairs(ℓₘᵢₙ, ℓₘₐₓ))
        # The integer formula, evaluated in `Rational`s, is a whole number here too
        @test n == (ℓₘₐₓ + 1)^2 - ℓₘᵢₙ^2
        # An independent count of the pairs in the ordering
        @test n == sum(2ℓ + 1 for ℓ in ℓₘᵢₙ:ℓₘₐₓ; init=0)
        @test n isa Int
        @test n ≥ 0
        # The `HalfOddInteger` spelling gives the same value
        @test Ysize(HalfOddInteger(ℓₘᵢₙ), HalfOddInteger(ℓₘₐₓ)) === n
    end
    for ℓₘᵢₙ in 1//2:11//2
        @test Ysize(ℓₘᵢₙ, ℓₘᵢₙ - 1) == 0
        @test Ysize(ℓₘᵢₙ, ℓₘᵢₙ) == 2ℓₘᵢₙ + 1
    end

    # The one-argument form starts at ℓₘᵢₙ = 1/2, the smallest half-odd ℓ
    for ℓₘₐₓ in -1//2:25//2
        @test Ysize(ℓₘₐₓ) == Ysize(1//2, ℓₘₐₓ)
        @test Ysize(ℓₘₐₓ) == (ℓₘₐₓ + 1)^2 - 1//4
        @test Ysize(ℓₘₐₓ) == length(ℓmpairs(1//2, ℓₘₐₓ))
        @test Ysize(HalfOddInteger(ℓₘₐₓ)) === Ysize(ℓₘₐₓ)
    end
    @test Ysize(-1//2) == 0

    # Validation is as for integers: ℓₘₐₓ below ℓₘᵢₙ - 1 is an error rather than a negative
    # size, and so is a negative ℓₘᵢₙ
    for ℓₘᵢₙ in 1//2:11//2, ℓₘₐₓ in -7//2:ℓₘᵢₙ-2
        @test isempty(ℓmpairs(ℓₘᵢₙ, ℓₘₐₓ))
        @test_throws ArgumentError Ysize(ℓₘᵢₙ, ℓₘₐₓ)
    end
    for ℓₘₐₓ in -7//2:-3//2
        @test_throws ArgumentError Ysize(ℓₘₐₓ)
    end
    for ℓₘᵢₙ in -5//2:-1//2, ℓₘₐₓ in -1//2:11//2
        @test_throws ArgumentError Ysize(ℓₘᵢₙ, ℓₘₐₓ)
    end
    # The error message names the offending argument
    @test_throws "ℓₘₐₓ" Ysize(1//2, -3//2)
    @test_throws "ℓₘᵢₙ" Ysize(-1//2, 3//2)
end


@testitem "Indexing: half-integer Yindex" begin
    import SphericalFunctions: Yindex, Ysize, HalfOddInteger

    # The ordering, built by brute force in `Rational`s: `Yindex(ℓ, m, ℓₘᵢₙ)` is supposed to
    # be the 1-based position of `(ℓ, m)` in this list
    ℓmpairs(ℓₘᵢₙ, ℓₘₐₓ) = [(ℓ, m) for ℓ in ℓₘᵢₙ:ℓₘₐₓ for m in -ℓ:ℓ]
    positions(ℓₘᵢₙ, ℓₘₐₓ) = Dict(p => i for (i, p) in enumerate(ℓmpairs(ℓₘᵢₙ, ℓₘₐₓ)))

    ℓₘₐₓ = 25//2

    # Agreement with the enumerated position for every (ℓ, m, ℓₘᵢₙ) with ℓₘᵢₙ ≤ ℓ ≤ ℓₘₐₓ, in
    # both spellings
    for ℓₘᵢₙ in 1//2:ℓₘₐₓ
        position = positions(ℓₘᵢₙ, ℓₘₐₓ)
        for ℓ in ℓₘᵢₙ:ℓₘₐₓ, m in -ℓ:ℓ
            i = Yindex(ℓ, m, ℓₘᵢₙ)
            @test i == position[(ℓ, m)]
            # The integer formula, evaluated in `Rational`s, is a whole number here too
            @test i == ℓ*(ℓ+1) - ℓₘᵢₙ^2 + m + 1
            @test i isa Int
            @test Yindex(HalfOddInteger(ℓ), HalfOddInteger(m), HalfOddInteger(ℓₘᵢₙ)) === i
        end
    end
    # ℓₘᵢₙ defaults to 1/2
    let position = positions(1//2, ℓₘₐₓ)
        for ℓ in 1//2:ℓₘₐₓ, m in -ℓ:ℓ
            @test Yindex(ℓ, m) == Yindex(ℓ, m, 1//2)
            @test Yindex(ℓ, m) == position[(ℓ, m)]
            @test Yindex(HalfOddInteger(ℓ), HalfOddInteger(m)) === Yindex(ℓ, m)
        end
    end

    # The structural invariants of the ordering, as in the integer item: first pair at
    # index 1, consecutive pairs at consecutive indices, the ℓ block ending at Ysize(ℓₘᵢₙ, ℓ)
    for ℓₘᵢₙ in 1//2:ℓₘₐₓ
        @test Yindex(ℓₘᵢₙ, -ℓₘᵢₙ, ℓₘᵢₙ) == 1
        for ℓ in ℓₘᵢₙ:ℓₘₐₓ
            for m in -ℓ:ℓ-1
                @test Yindex(ℓ, m+1, ℓₘᵢₙ) == Yindex(ℓ, m, ℓₘᵢₙ) + 1
            end
            @test Yindex(ℓ+1, -ℓ-1, ℓₘᵢₙ) == Yindex(ℓ, ℓ, ℓₘᵢₙ) + 1
            @test Yindex(ℓ, ℓ, ℓₘᵢₙ) == Ysize(ℓₘᵢₙ, ℓ)
            @test Yindex(ℓ, -ℓ, ℓₘᵢₙ) == Ysize(ℓₘᵢₙ, ℓ - 1) + 1
        end
        # Altogether, the indices enumerate 1:Ysize exactly once each, in order
        indices = [Yindex(ℓ, m, ℓₘᵢₙ) for ℓ in ℓₘᵢₙ:ℓₘₐₓ for m in -ℓ:ℓ]
        @test indices == 1:Ysize(ℓₘᵢₙ, ℓₘₐₓ)
        @test indices == 1:length(ℓmpairs(ℓₘᵢₙ, ℓₘₐₓ))
    end
end


@testitem "Indexing: half-integer Yrange" begin
    import SphericalFunctions: Yrange, Yindex, Ysize, HalfOddInteger

    # The ordering, built by brute force in `Rational`s: `Yrange` is supposed to return
    # exactly this list, as `HalfOddInteger` pairs
    ℓmpairs(ℓₘᵢₙ, ℓₘₐₓ) = [(ℓ, m) for ℓ in ℓₘᵢₙ:ℓₘₐₓ for m in -ℓ:ℓ]

    for ℓₘᵢₙ in 1//2:13//2, ℓₘₐₓ in ℓₘᵢₙ-1:25//2
        r = Yrange(ℓₘᵢₙ, ℓₘₐₓ)
        @test r isa AbstractVector
        # Decision E1: the pairs are of `HalfOddInteger`s, however the call was spelled, and
        # the element type is the same for the empty range
        @test eltype(r) === Tuple{HalfOddInteger, HalfOddInteger}
        @test r == Yrange(HalfOddInteger(ℓₘᵢₙ), HalfOddInteger(ℓₘₐₓ))
        @test length(r) == Ysize(ℓₘᵢₙ, ℓₘₐₓ)
        # The list itself, pair by pair, compared through `HalfOddInteger == Rational` and
        # again after converting back to the `Rational` spelling
        @test r == ℓmpairs(ℓₘᵢₙ, ℓₘₐₓ)
        rational_pairs = [(Rational(ℓ), Rational(m)) for (ℓ, m) in r]
        @test rational_pairs == ℓmpairs(ℓₘᵢₙ, ℓₘₐₓ)
        # `Yrange` and `Yindex` are inverse to each other
        for ℓ in ℓₘᵢₙ:ℓₘₐₓ, m in -ℓ:ℓ
            @test r[Yindex(ℓ, m, ℓₘᵢₙ)] == (ℓ, m)
        end
        for (i, (ℓ, m)) in enumerate(r)
            @test Yindex(ℓ, m, HalfOddInteger(ℓₘᵢₙ)) == i
            @test ℓₘᵢₙ ≤ Rational(ℓ) ≤ ℓₘₐₓ
            @test -ℓ ≤ m ≤ ℓ
        end
        # ℓ-major, m increasing, no repeats
        @test issorted(r)
        @test allunique(rational_pairs)
        if !isempty(r)
            @test first(r) == (ℓₘᵢₙ, -ℓₘᵢₙ)
            @test last(r) == (ℓₘₐₓ, ℓₘₐₓ)
        end
    end

    # The empty range
    @test isempty(Yrange(1//2, -1//2))
    @test isempty(Yrange(-1//2))
    @test eltype(Yrange(1//2, -1//2)) === Tuple{HalfOddInteger, HalfOddInteger}

    # The one-argument form is ℓₘᵢₙ = 1/2
    for ℓₘₐₓ in -1//2:25//2
        @test Yrange(ℓₘₐₓ) == Yrange(1//2, ℓₘₐₓ)
        @test Yrange(ℓₘₐₓ) == ℓmpairs(1//2, ℓₘₐₓ)
        @test Yrange(HalfOddInteger(ℓₘₐₓ)) == Yrange(ℓₘₐₓ)
    end
end


@testitem "Indexing: half-integer spellings and mixed kinds" begin
    import SphericalFunctions: Ysize, Yindex, Yrange, HalfOddInteger

    # The `Rational` and `HalfOddInteger` spellings of a half-odd index are interchangeable,
    # argument by argument
    h = HalfOddInteger
    @test Ysize(1//2, 7//2) === Ysize(h(1//2), h(7//2)) === Ysize(1//2, h(7//2)) === Ysize(h(1//2), 7//2)
    @test Ysize(7//2) === Ysize(h(7//2))
    @test Yindex(5//2, -3//2, 1//2) === Yindex(h(5//2), h(-3//2), h(1//2)) === Yindex(5//2, h(-3//2), 1//2)
    @test Yindex(5//2, -3//2) === Yindex(h(5//2), h(-3//2)) === Yindex(h(5//2), -3//2)
    @test Yrange(1//2, 5//2) == Yrange(h(1//2), h(5//2)) == Yrange(1//2, h(5//2))
    @test Yrange(5//2) == Yrange(h(5//2))

    # Mixing an integer index with a half-odd one is refused, with a message that names both
    # spellings and the offending values, rather than with a bare `MethodError`
    msg = "all be integers, like 3, or all be half-odd-integers, like 7//2"
    @test_throws ArgumentError Ysize(0, 7//2)
    @test_throws msg Ysize(0, 7//2)
    @test_throws msg Ysize(7//2, 0)
    @test_throws msg Ysize(0, h(7//2))
    @test_throws msg Yindex(3//2, 1//2, 0)
    @test_throws msg Yindex(3//2, 0)
    @test_throws msg Yindex(1, h(1//2))
    @test_throws msg Yrange(0, 7//2)
    @test_throws msg Yrange(1//2, 7)
    @test_throws "got 0, 7//2" Ysize(0, 7//2)

    # A `Rational` with denominator 1 is not accepted as a spelling of an integer index.
    # `half_integer` rejects every `Rational` whose denominator is not 2, as it does at the
    # Wigner constructors, so an integer index is spelled as an `Integer`; this is
    # deliberate, and the message says what was expected
    @test_throws "must have denominator 2" Ysize(3//1)
    @test_throws "must have denominator 2" Ysize(0//1, 3//1)
    @test_throws "must have denominator 2" Yindex(3//1, 1//1)
    @test_throws "must have denominator 2" Yindex(3//1, 1//1, 0//1)
    @test_throws "must have denominator 2" Yrange(3//1)
    @test_throws "must have denominator 2" Yrange(0//1, 3//1)
    # ... as are the other denominators
    @test_throws "must have denominator 2" Ysize(1//2, 5//3)
    @test_throws "must have denominator 2" Yindex(1//2, 1//4)
end


@testitem "Indexing: half-integer arithmetic stays in `Int`" begin
    import SphericalFunctions: Ysize, Yindex, Yrange, HalfOddInteger
    using Test: @inferred

    # The closed forms on `HalfOddInteger` arguments infer to `Int`, through every method
    # including the `Rational` boundary, and no `Rational` appears in the compiled code.  A
    # correctness test could not see this: the answers would be the same either way.
    ℓ, m, ℓₘᵢₙ = HalfOddInteger(7//2), HalfOddInteger(-3//2), HalfOddInteger(1//2)
    @test @inferred(Ysize(ℓₘᵢₙ, ℓ)) === Int((7//2 + 1)^2 - (1//2)^2)
    @test @inferred(Ysize(ℓ)) === Ysize(ℓₘᵢₙ, ℓ)
    @test @inferred(Yindex(ℓ, m, ℓₘᵢₙ)) === Int((7//2)*(7//2 + 1) - (1//2)^2 + (-3//2) + 1)
    @test @inferred(Yindex(ℓ, m)) === Yindex(ℓ, m, ℓₘᵢₙ)
    @test @inferred(Yrange(ℓₘᵢₙ, ℓ)) isa Vector{Tuple{HalfOddInteger, HalfOddInteger}}

    H, Q = HalfOddInteger, Rational{Int}
    for sig in ((H, H), (H,), (Q, Q), (Q,), (H, Q), (Q, H))
        @test Base.return_types(Ysize, sig) == [Int]
    end
    for sig in ((H, H, H), (H, H), (Q, Q, Q), (Q, Q), (H, Q, H), (Q, H))
        @test Base.return_types(Yindex, sig) == [Int]
    end
    for sig in ((H, H), (H,), (Q, Q), (Q,))
        @test Base.return_types(Yrange, sig) == [Vector{Tuple{H, H}}]
    end
    for (f, sig) in ((Ysize, (H, H)), (Ysize, (H,)), (Yindex, (H, H, H)), (Yindex, (H, H)))
        ir = string(Base.code_typed(f, sig; optimize=true)[1][1])
        @test !occursin("Rational", ir)
    end

    # The integer path is untouched
    @test @inferred(Ysize(2, 5)) === 32
    @test @inferred(Yindex(3, -2)) === 11
    @test Base.return_types(Ysize, (Int, Int)) == [Int]
    @test Base.return_types(Yindex, (Int, Int, Int)) == [Int]
end
