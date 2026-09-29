# Tests of `HalfOddInteger` itself.
#
# The package's own test suite would pass against a badly broken version of this type — it
# exercises only the handful of operations the recurrences perform, and would not notice, for
# example, that `step` on a range of these had silently become `0`.  The failure modes here
# are mostly silent rather than loud, so they need testing directly.

@testitem "HalfOddInteger: construction and conversion" begin
    using SphericalFunctions: HalfOddInteger

    # The public constructor takes a *value*, so that it agrees with `convert` — which is
    # what the calculators call on a user-supplied `ℓ`.
    @test HalfOddInteger(5//2) === convert(HalfOddInteger, 5//2)
    @test HalfOddInteger(HalfOddInteger(5//2)) === HalfOddInteger(5//2)
    @test numerator(HalfOddInteger(5//2)) == 5
    @test denominator(HalfOddInteger(5//2)) == 2
    @test numerator(HalfOddInteger(-5//2)) == -5

    # A whole number is not a half-odd-integer, however spelled.
    @test_throws InexactError HalfOddInteger(2)
    @test_throws InexactError convert(HalfOddInteger, 2)
    @test_throws "must have denominator 2" HalfOddInteger(2//1)
    @test_throws "must have denominator 2" HalfOddInteger(4//2)
    @test_throws "must have denominator 2" HalfOddInteger(5//3)
    @test_throws DomainError HalfOddInteger(5//3)

    # The numerator of a `Rational` of any integer type is converted to the stored `Int`;
    # only a numerator that does not fit an `Int` is refused, with the ordinary
    # `InexactError`
    for x in (Int8(7)//Int8(2), Int16(7)//Int16(2), Int32(7)//Int32(2), UInt8(7)//UInt8(2), big(7)//2)
        @test HalfOddInteger(x) === HalfOddInteger(7//2)
    end
    @test HalfOddInteger(Int8(-7)//Int8(2)) === HalfOddInteger(-7//2)
    @test HalfOddInteger(big(-1)//2) === HalfOddInteger(-1//2)
    @test_throws InexactError HalfOddInteger((big(2)^70 + 1)//2)

    # Round trips, through `Rational` and through every float type of `Base`: exact wherever
    # the type can represent n/2, and otherwise correctly rounded, as for 2^40 + 1 in
    # `Float32` (where `Float32(n)/2` is the rounding of n/2, since halving is exact)
    for n in (-7, -1, 1, 3, 9, 2^40 + 1)
        x = HalfOddInteger(n//2)
        @test Rational(x) === n//2
        @test Rational{Int}(x) === n//2
        @test Rational{BigInt}(x) == n//2 && Rational{BigInt}(x) isa Rational{BigInt}
        @test Float64(x) === n/2
        @test float(x) === n/2
        @test AbstractFloat(x) === n/2
        @test convert(Float64, x) === n/2
        @test Float32(x) === Float32(n)/2
        @test BigFloat(x) == big(n)/2 && BigFloat(x) isa BigFloat
        abs(n) < 2048 && @test Float16(x) === Float16(n)/2
        @test x == n//2
        @test n//2 == x
        @test x != n
        @test string(x) == "$n//2"   # displayed as the caller would write it
    end
end


@testitem "HalfOddInteger: conversion to the DoubleFloats types" begin
    using SphericalFunctions: HalfOddInteger, index_value
    using DoubleFloats: Double16, Double32, Double64

    # The package extension gives each of these a conversion of its own, which is exact.  A
    # conversion defined once for every `AbstractFloat` would instead be ambiguous with
    # `DoubleFloats`' own `Double64(x::T) where {T<:Real}`.
    h(x) = HalfOddInteger(x)
    for T in (Double16, Double32, Double64), a in (-7//2, -1//2, 1//2, 3//2, 21//2)
        @test T(h(a)) isa T
        @test T(h(a)) == T(a)
        @test convert(T, h(a)) == T(a)
        @test h(a) == T(a) && T(a) == h(a)
    end
    # ... even for a numerator wider than the significand of the high word, which
    # `DoubleFloats`' conversion of an `Int` would round
    @test Double64(h((2^61 + 1)//2)) == Double64(big(2)^61 + 1) / 2
    @test Double32(h((2^45 + 1)//2)) == Double32(big(2)^45 + 1) / 2
    @test Double16(h((2^15 + 1)//2)) == Double16(big(2)^15 + 1) / 2

    # `index_value` gives an index as an element of a float type.  An integer is returned as
    # it is, for the typed comprehension that receives it to convert; a half-odd-integer
    # becomes the float exactly, through its numerator, in every float type
    @test index_value(Float64, 3) === 3
    @test index_value(Float32, -2) === -2
    for T in (Float16, Float32, Float64, BigFloat, Double64), a in (-7//2, -1//2, 1//2, 3//2, 21//2)
        v = index_value(T, h(a))
        @test v isa T
        @test v == a
    end
end


@testitem "HalfOddInteger: arithmetic" begin
    using SphericalFunctions: HalfOddInteger

    a, b = HalfOddInteger(7//2), HalfOddInteger(3//2)
    h(x) = HalfOddInteger(x)

    # Sums and differences of two of them leave the type; adding an integer stays in it.
    @test a + b === 5
    @test a - b === 2
    @test b - a === -2
    @test a + 1 === HalfOddInteger(9//2)
    @test 1 + a === HalfOddInteger(9//2)
    @test a - 1 === HalfOddInteger(5//2)
    @test 1 - a === HalfOddInteger(-5//2)
    @test -a === HalfOddInteger(-7//2)
    @test abs(HalfOddInteger(-7//2)) === a

    # An integer operand of any type is converted to `Int` first, so that the result is the
    # same as for an `Int`, or an `InexactError` if the value does not fit
    for n in (Int8(3), Int16(3), Int32(3), UInt8(3), UInt(3), big(3), Int128(3))
        @test a + n === h(13//2)
        @test n + a === h(13//2)
        @test a - n === h(1//2)
        @test n - a === h(-1//2)
    end
    @test a + true === h(9//2) && a - false === a
    @test_throws InexactError a + typemax(UInt)
    @test_throws InexactError a - big(2)^70
    @test_throws InexactError big(2)^70 * a

    # Multiplication is defined only by an even integer, and gives an `Int`, whatever the
    # type of the multiplier.  An unsigned multiplier gives the signed product.
    @test 2a === 7
    @test a * 2 === 7
    @test 4a === 14
    @test -2a === -7
    @test 0a === 0
    @test UInt(2) * h(-1//2) === -1
    @test big(2) * h(-1//2) === -1
    @test Int8(4) * h(-7//2) === -14
    # An odd multiplier is outside the domain of the operation
    @test_throws DomainError 3a
    @test_throws DomainError a * 3
    @test_throws DomainError UInt(3) * a
    @test_throws "may only be multiplied by an even integer" 3a
    let e = try -5a catch e; e end
        @test e isa DomainError && e.val === -5
    end

    # Comparison, against itself and against integers.
    @test b < a && a > b && b ≤ a && a ≥ b
    @test HalfOddInteger(-1//2) < 0 < HalfOddInteger(1//2)
    @test HalfOddInteger(1//2) ≤ 1 && !(HalfOddInteger(3//2) ≤ 1)
    @test sort([a, b, -a]) == [-a, b, a]
    # A half-odd-integer lies strictly between two integers, and the comparisons agree with
    # `Rational` arithmetic at every sign, including against integers of other types
    for x in -9//2:9//2, n in -5:5
        @test (h(x) < n) === (x < n) === (h(x) ≤ n)
        @test (n < h(x)) === (n < x) === (n ≤ h(x))
        @test (h(x) > n) === (x > n) && (h(x) ≥ n) === (x ≥ n)
        @test (h(x) < Int8(n)) === (x < n) && (UInt8(n + 5) ≤ h(x + 5)) === (n ≤ x)
        @test (h(x) < big(n)) === (x < n) && (big(n) < h(x)) === (n < x)
    end
    # ... and at the extremes of `Int`, where doubling the integer would overflow
    for x in (h(1//2), h(-1//2), h(typemax(Int)//2), h(-typemax(Int)//2))
        @test x < typemax(Int) && x ≤ typemax(Int)
        @test typemax(Int) > x && typemax(Int) ≥ x
        @test x > typemin(Int) && x ≥ typemin(Int)
        @test typemin(Int) < x && typemin(Int) ≤ x
        @test x < typemax(UInt) && x < big(2)^100 && x > -big(2)^100
    end
    @test !(h(typemax(Int)//2) < typemax(Int) ÷ 2) && h(typemax(Int)//2) < typemax(Int) ÷ 2 + 1
    @test h(-typemax(Int)//2) < -(typemax(Int) ÷ 2) && !(h(-typemax(Int)//2) < -(typemax(Int) ÷ 2) - 1)
end


@testitem "HalfOddInteger: equality with Rationals and floats, and isinteger" begin
    using SphericalFunctions: HalfOddInteger
    using DoubleFloats: Double64

    h(x) = HalfOddInteger(x)
    for n in (-2001, -7, -1, 1, 3, 101)
        x = h(n//2)
        # Equal to the `Rational` and to the float of the same value, in either order, and
        # `isequal` to both, as `Base`'s own numbers are
        for y in (n//2, big(n)//2, n/2, Float32(n/2), big(n)/2, Double64(n)/2)
            @test x == y && y == x
            @test isequal(x, y) && isequal(y, x)
        end
        # Never equal to a whole number, however it is written, nor to a neighbouring value
        for y in (n ÷ 2, n//1, float(n), n/2 + 1/4, n/2 + 1, nextfloat(n/2), prevfloat(n/2), NaN, Inf)
            @test x != y && y != x
            @test !isequal(x, y) && !isequal(y, x)
        end
        @test !isinteger(x)
    end
    # Doubling a float near the top of its range overflows to `Inf`, which is equal to nothing
    @test h(1//2) != floatmax(Float64) && h(1//2) != floatmax(Float16)
    @test h(typemax(Int)//2) != floatmax(Float64)
end


@testitem "HalfOddInteger: values that do not exist" begin
    using SphericalFunctions: HalfOddInteger

    # `Base`'s `Number` fall-backs would form these by conversion — `one(::Type{T})` is
    # `convert(T, 1)` — and refuse with an `InexactError` about a conversion the caller
    # never asked for.  These say what is wrong instead.  A value silently returned from
    # `one` would make `step` on a range of these wrong, which is not detectable downstream.
    a = HalfOddInteger(3//2)
    @test_throws ArgumentError zero(HalfOddInteger)
    @test_throws ArgumentError one(HalfOddInteger)
    @test_throws ArgumentError oneunit(HalfOddInteger)
    @test_throws ArgumentError zero(a)
    @test_throws ArgumentError one(a)
    @test_throws ArgumentError oneunit(a)
    @test_throws "has no zero" zero(a)
    @test_throws "has no one" one(a)
    @test_throws "has no oneunit" oneunit(a)
    @test_throws InexactError Int(a)
end


@testitem "HalfOddInteger: ranges" begin
    using SphericalFunctions: HalfOddInteger

    h(x) = HalfOddInteger(x)
    lo, hi = h(-5//2), h(5//2)
    r = lo:hi
    @test r isa UnitRange{HalfOddInteger}
    @test step(r) === 1
    @test length(r) === 6
    @test first(r) === lo && last(r) === hi
    @test collect(r) == HalfOddInteger.([-5//2, -3//2, -1//2, 1//2, 3//2, 5//2])
    @test h(1//2) ∈ r
    @test h(7//2) ∉ r

    # An empty range must report length 0 rather than a negative number.
    @test length(h(5//2):h(1//2)) == 0
    @test isempty(h(5//2):h(1//2))
    @test isempty(collect(h(5//2):h(1//2)))

    # The descending form, which step 5 of the recurrence uses, and `reverse` gives.
    rd = hi:-1:lo
    @test rd isa StepRange{HalfOddInteger, Int}
    @test reverse(r) == rd && reverse(r) isa StepRange{HalfOddInteger, Int}
    @test collect(rd) == reverse(collect(r))

    # Membership means `==` to some element, as for `Base`'s ranges, whether the range runs
    # upward or downward, and however the value is written
    for range in (r, rd, reverse(r))
        for x in (-5//2, -1//2, 5//2, h(3//2), -1.5, 2.5, big(1)//2)
            @test x ∈ range
        end
        for x in (-7//2, 7//2, h(-7//2), 0, 2, -3, 0.0, 1.0, 0.25, 3.5, 1//3, 3//1, NaN, Inf)
            @test x ∉ range
        end
    end
    # ... including a range whose step is not 1, whose members are a whole number of steps
    # apart, and an empty one
    for (range, members) in (
        (h(-5//2):2:h(5//2), (-5//2, -1//2, 3//2)),
        (h(5//2):-2:h(-5//2), (5//2, 1//2, -3//2)),
        (h(-7//2):3:h(7//2), (-7//2, -1//2, 5//2)),
        (h(1//2):1:h(-1//2), ()),
        (h(1//2):-1:h(3//2), ()),
    )
        @test collect(range) == collect(members)
        for n in -9:2:9
            @test (n//2 ∈ range) === (n//2 ∈ members)
            @test (h(n//2) ∈ range) === (n//2 ∈ members)
            @test (n/2 ∈ range) === (n//2 ∈ members)
        end
        @test 0 ∉ range && 1 ∉ range && 0.25 ∉ range
    end
    # A half-odd-integer is not a member of a range of integers
    @test h(1//2) ∉ 0:3
    @test h(1//2) ∉ 3:-1:0
    @test h(-1//2) ∉ -3:3

    # Iterating, and testing membership, allocate nothing.
    f(r) = (s = 0; for m in r; s += (m - first(r)); end; s)
    g(r, xs) = count(x -> x ∈ r, xs)
    xs = [-5//2, 1//2, 2//1, 7//2]
    f(r); f(rd); g(r, xs); g(rd, xs)
    @test (@allocated f(r)) == 0
    @test (@allocated f(rd)) == 0
    @test g(r, xs) == g(rd, xs) == 2
    @test (@allocated g(r, xs)) == 0
    @test (@allocated g(rd, xs)) == 0
end


@testitem "HalfOddInteger: the recurrence arithmetic stays in `Int`" setup=[InferenceChecks] begin
    using SphericalFunctions: HalfOddInteger, HCalculator, δ², sgn, ϵ
    import .InferenceChecks: inferred_type, dynamic_calls

    # This is the whole point of the type, and it is exactly the property that no
    # correctness test can see: if these inferred `Rational` instead, every answer would
    # still be right and the half-integer path would run roughly thirty times slower.
    @test inferred_type(ℓ -> 2ℓ + 1, (HalfOddInteger,)) === Int
    @test inferred_type(δ², (HalfOddInteger, HalfOddInteger)) === Int
    @test inferred_type(δ², (Int, Int)) === Int
    @test inferred_type(ϵ, (HalfOddInteger,)) === Int
    @test inferred_type(sgn, (HalfOddInteger,)) === Int
    # ... and `2ℓ` compiles to nothing but a read of the numerator: the check that refuses an
    # odd multiplier folds away for the literal 2, so that nothing is called
    @test dynamic_calls(ℓ -> 2ℓ, (HalfOddInteger,)) == 0
    let code = only(Base.code_typed(ℓ -> 2ℓ, (HalfOddInteger,); optimize=true)).first.code
        @test !any(ex -> Meta.isexpr(ex, :invoke), code)
        @test !any(ex -> ex isa Core.GotoIfNot, code)
    end

    # ... and no `Rational` survives anywhere in the hot loops themselves.  The types
    # inspected are those of real calculators, single and batched:
    # `HCalculator{HalfOddInteger, Float64}` is missing its storage parameter, and code
    # typed for that `UnionAll` reads the storage as `Any` and dispatches dynamically, which
    # would hide a `Rational` that appears only in the specialization that runs.  (The
    # absence of dynamic calls shows that the whole loop was inferred, so that the absence
    # of `Rational` means something.)
    for H ∈ (HCalculator(0.3, 7//2), HCalculator([0.3, 0.4], 7//2))
        @test isconcretetype(typeof(H))
        @test typeof(H) <: HCalculator{HalfOddInteger, Float64}
        for step! in (SphericalFunctions.recurrence_step4!, SphericalFunctions.recurrence_step5!,
                      SphericalFunctions.recurrence_seed!)
            ir = string(Base.code_typed(step!, (typeof(H),); optimize=true)[1][1])
            @test !occursin("Rational", ir)
            @test dynamic_calls(step!, (typeof(H),)) == 0
        end
    end

    # `ϵ` agrees with the closed form for both index types, from the one definition.
    @test all(ϵ(HalfOddInteger(n//2)) == (n > 0 && isodd((n-1)÷2) ? -1 : 1) for n in -21:2:21)
    @test all(ϵ(m) == (m > 0 && isodd(m) ? -1 : 1) for m in -10:10)
end


@testitem "HalfOddInteger: the floor of an index" setup=[InferenceChecks] begin
    using SphericalFunctions: HalfOddInteger, floor_int
    import .InferenceChecks: inferred_type

    # It gives an `Int`, and the one expression is right at both signs.  The floating-point
    # `floor` is the reference.
    for n in -21:2:21
        x = n//2
        a = HalfOddInteger(x)
        @test floor_int(a) === Int(floor(x)) === Int(floor(float(x)))
        @test floor_int(a) < a < floor_int(a) + 1
    end
    @test floor_int(HalfOddInteger(7//2)) === 3
    @test floor_int(HalfOddInteger(-1//2)) === -1
    @test floor_int(HalfOddInteger(-7//2)) === -4
    # It does not overflow at the extremes of the numerator
    @test floor_int(HalfOddInteger(typemax(Int)//2)) === typemax(Int) ÷ 2
    @test floor_int(HalfOddInteger((typemin(Int) + 1)//2)) === typemin(Int) ÷ 2

    # For an integer it is the identity, as an `Int`
    for n in (-3, 0, 2, Int8(5), UInt8(7), big(4), true)
        @test floor_int(n) === Int(n)
    end
    @test_throws InexactError floor_int(big(2)^70)
    @test inferred_type(floor_int, (HalfOddInteger,)) === Int
    @test inferred_type(floor_int, (Int8,)) === Int

    # `Base.floor` and `Base.ceil` have no methods for the type, because a method of either
    # for a new argument type would invalidate compiled code of `Base`
    @test_throws MethodError floor(Int, HalfOddInteger(1//2))
    @test_throws MethodError ceil(HalfOddInteger(1//2))
end


@testitem "HalfOddInteger: a Rational with denominator 1 is not an index" begin
    using SphericalFunctions: ℓ
    using SphericalFunctions: Ylm, D, ModeWeights
    using Quaternionic: Rotor

    # An integer index is an `Integer`, and a half-odd one a `HalfOddInteger` or a
    # `Rational` with denominator 2, so every container refuses `2//1` as an ℓ, as a
    # `ModeWeights` does
    R = Rotor(1.0)
    Y = Ylm(R, 3)
    𝔇 = D(R, 3)
    w = ModeWeights(zeros(16))
    @test ℓ(Y[2]) == ℓ(𝔇[2]) == ℓ(w[2, :]) == 2
    @test_throws ArgumentError Y[2//1]
    @test_throws ArgumentError 𝔇[2//1]
    @test_throws ArgumentError w[2//1, :]
    @test_throws "2//1" Y[2//1]
    @test_throws "2//1" 𝔇[2//1]
    @test_throws "2//1" w[2//1, :]
end


@testitem "HalfOddInteger: hashing agrees with Rational" begin
    using SphericalFunctions: HalfOddInteger, Yrange, Ysize, Yindex

    h = HalfOddInteger

    # A half-odd-integer is `isequal` to the `Rational` of the same value, so the two must
    # hash alike, and so must the `Float64`, which is `isequal` to the `Rational`
    for x ∈ (1//2, -1//2, 7//2, -7//2, 101//2, -2001//2)
        @test isequal(h(x), x)
        @test isequal(h(x), float(x))
        @test hash(h(x)) == hash(x) == hash(float(x)) == hash(Float32(x)) == hash(big(x))
        @test hash(h(x), UInt(1729)) == hash(x, UInt(1729)) == hash(float(x), UInt(1729))
    end
    @test hash(h(1//2)) != hash(h(3//2))
    @test hash(h(1//2)) != hash(h(-1//2))

    # The index values the package hands back can therefore be collected, deduplicated and
    # used as keys, as integer indices can
    r = Yrange(1//2, 25//2)
    @test allunique(r)
    @test length(Set(r)) == length(r) == Ysize(1//2, 25//2)
    @test length(unique(vcat(r, r))) == length(r)
    @test Rational.(unique(first.(r))) == [(2k + 1)//2 for k ∈ 0:12]
    @test unique([h(1//2), h(1//2), h(3//2)]) == [h(1//2), h(3//2)]
    positions = Dict(zip(r, eachindex(r)))
    @test positions[(h(3//2), h(-1//2))] == Yindex(3//2, -1//2)
    # ... and a key finds the equal value of the other type
    @test positions[(3//2, -1//2)] == Yindex(3//2, -1//2)
    @test positions[(1.5, -0.5)] == Yindex(3//2, -1//2)
    @test Dict(7//2 => 1)[h(7//2)] == 1
    @test Dict(3.5 => 1)[h(7//2)] == 1
    @test Dict(h(7//2) => 1)[3.5] == 1
    @test h(7//2) ∈ Set([7//2])
    @test h(7//2) ∈ Set([3.5])
end
