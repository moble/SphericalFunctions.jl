# Tests of `HalfOddInteger` itself.
#
# The package's own test suite would pass against a badly broken version of this type — it
# exercises only the handful of operations the recurrences perform, and would not notice, for
# example, that `step` on a range of these had silently become `0`.  The failure modes here
# are mostly silent rather than loud, so they need testing directly.

@testitem "HalfOddInteger: construction and conversion" begin
    using SphericalFunctions: HalfOddInteger, half_integer

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

    # `half_integer` normalizes whichever spelling the caller used.
    @test half_integer(3) === 3
    @test half_integer(3//2) === HalfOddInteger(3//2)
    @test half_integer(HalfOddInteger(3//2)) === HalfOddInteger(3//2)
    @test_throws "must have denominator 2" half_integer(3//4)

    # Round trips.
    for n in (-7, -1, 1, 3, 9)
        x = HalfOddInteger(n//2)
        @test Rational(x) == n//2
        @test Rational{Int}(x) == n//2
        @test Float64(x) == n/2
        @test float(x) == n/2
        @test x == n//2
        @test n//2 == x
        @test x != n
        @test string(x) == "$n//2"   # displayed as the caller would write it
    end
end


@testitem "HalfOddInteger: arithmetic" begin
    using SphericalFunctions: HalfOddInteger

    a, b = HalfOddInteger(7//2), HalfOddInteger(3//2)

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

    # Multiplication is defined only by an even integer, and gives an `Integer`.
    @test 2a === 7
    @test a * 2 === 7
    @test 4a === 14
    @test_throws AssertionError 3a
    @test_throws AssertionError a * 3

    # Comparison, against itself and against integers.
    @test b < a && a > b && b ≤ a && a ≥ b
    @test HalfOddInteger(-1//2) < 0 < HalfOddInteger(1//2)
    @test HalfOddInteger(1//2) ≤ 1 && !(HalfOddInteger(3//2) ≤ 1)
    @test sort([a, b, -a]) == [-a, b, a]

    # `floor` is what ε(m) = (-1)^⌊m⌋ needs.
    @test floor(Int, HalfOddInteger(7//2)) == 3
    @test floor(Int, HalfOddInteger(-1//2)) == -1
    @test floor(HalfOddInteger(-7//2)) == -4
end


@testitem "HalfOddInteger: values that do not exist" begin
    using SphericalFunctions: HalfOddInteger

    # `Base`'s `Number` fall-backs would otherwise manufacture these from the positional
    # constructor: `one(::Type{T})` is `convert(T, 1)`.  Silently returning 1/2 from `one`
    # would then make `step` on a range of these `0`, which is not detectable downstream.
    a = HalfOddInteger(3//2)
    @test_throws ArgumentError zero(HalfOddInteger)
    @test_throws ArgumentError one(HalfOddInteger)
    @test_throws ArgumentError oneunit(HalfOddInteger)
    @test_throws ArgumentError zero(a)
    @test_throws ArgumentError one(a)
    @test_throws ArgumentError oneunit(a)
    @test_throws InexactError Int(a)
end


@testitem "HalfOddInteger: ranges" begin
    using SphericalFunctions: HalfOddInteger

    lo, hi = HalfOddInteger(-5//2), HalfOddInteger(5//2)
    r = lo:hi
    @test r isa UnitRange{HalfOddInteger}
    @test step(r) == 1
    @test length(r) == 6
    @test first(r) === lo && last(r) === hi
    @test collect(r) == HalfOddInteger.([-5//2, -3//2, -1//2, 1//2, 3//2, 5//2])
    @test HalfOddInteger(1//2) ∈ r
    @test HalfOddInteger(7//2) ∉ r

    # An empty range must report length 0 rather than a negative number.
    @test length(HalfOddInteger(5//2):HalfOddInteger(1//2)) == 0
    @test isempty(HalfOddInteger(5//2):HalfOddInteger(1//2))

    # The descending form, which step 5 of the recurrence uses.
    rd = hi:-1:lo
    @test rd isa StepRange{HalfOddInteger, Int}
    @test collect(rd) == reverse(collect(r))

    # Iterating allocates nothing.
    f(r) = (s = 0; for m in r; s += (m - first(r)); end; s)
    f(r)
    @test (@allocated f(r)) == 0
end


@testitem "HalfOddInteger: the recurrence arithmetic stays in `Int`" begin
    using SphericalFunctions: HalfOddInteger, WignerHCalculator, δ², sgn, ϵ

    # This is the whole point of the type, and it is exactly the property that no correctness
    # test can see: if these inferred `Rational` instead, every answer would still be right
    # and the half-integer path would run roughly thirty times slower.
    @test Base.infer_return_type(ℓ -> 2ℓ + 1, (HalfOddInteger,)) === Int
    @test Base.infer_return_type(δ², (HalfOddInteger, HalfOddInteger)) === Int
    @test Base.infer_return_type(δ², (Int, Int)) === Int
    @test Base.infer_return_type(ϵ, (HalfOddInteger,)) === Int
    @test Base.infer_return_type(sgn, (HalfOddInteger,)) === Int

    # ... and no `Rational` survives anywhere in the hot loops themselves.
    for step! in (SphericalFunctions.recurrence_step4!, SphericalFunctions.recurrence_step5!,
                  SphericalFunctions.recurrence_seed!)
        ir = string(Base.code_typed(
            step!, (WignerHCalculator{HalfOddInteger, Float64},); optimize=true
        )[1][1])
        @test !occursin("Rational", ir)
    end

    # `ϵ` agrees with the closed form for both index types, from the one definition.
    @test all(ϵ(HalfOddInteger(n//2)) == (n > 0 && isodd((n-1)÷2) ? -1 : 1) for n in -21:2:21)
    @test all(ϵ(m) == (m > 0 && isodd(m) ? -1 : 1) for m in -10:10)
end


@testitem "HalfOddInteger: half-integer ceil" begin
    using SphericalFunctions: HalfOddInteger

    # `ceil`, like `floor`, gives a whole number and so returns an `Int`, and the one
    # expression is right at both signs.  The floating-point `ceil` is the reference.
    for x in (-5//2, -3//2, -1//2, 1//2, 3//2, 5//2)
        a = HalfOddInteger(x)
        @test ceil(a) == ceil(float(x))
        @test ceil(a) isa Int
        @test ceil(Int, a) === ceil(a)
        @test ceil(Int8, a) === Int8(ceil(float(x)))
        # The floor and the ceiling bracket a half-odd-integer
        @test ceil(a) == floor(a) + 1
        @test floor(a) < x < ceil(a)
    end
    @test ceil(HalfOddInteger(-1//2)) === 0
    @test ceil(HalfOddInteger(1//2)) === 1

    # `ceil` of an `Integer` is the identity, and is unchanged by the new methods
    for n in (-3, 0, 2, Int8(5))
        @test ceil(n) === n
        @test ceil(Int, n) === Int(n)
    end
end


@testitem "HalfOddInteger: half-integer indices agree in kind" begin
    using SphericalFunctions: HalfOddInteger, half_integer, half_integers

    # `half_integers` normalizes each argument as `half_integer` does, and returns the tuple
    # when the results are all of one kind — integers of different types included
    @test half_integers(3, 4, 5) === (3, 4, 5)
    @test half_integers(Int8(3), 4) === (Int8(3), 4)
    @test half_integers(3//2, 5//2) === (HalfOddInteger(3//2), HalfOddInteger(5//2))
    @test half_integers(HalfOddInteger(3//2), 5//2) === (HalfOddInteger(3//2), HalfOddInteger(5//2))
    @test half_integers(7//2) === (HalfOddInteger(7//2),)
    @test half_integers(7) === (7,)

    # A mixture of the two kinds is refused, with a message that names both spellings and
    # the offending values
    msg = "all be integers, like 3, or all be half-odd-integers, like 7//2"
    @test_throws ArgumentError half_integers(0, 7//2)
    @test_throws msg half_integers(0, 7//2)
    @test_throws msg half_integers(7//2, 0)
    @test_throws msg half_integers(HalfOddInteger(7//2), 0)
    @test_throws msg half_integers(1//2, 3//2, 2)
    @test_throws "got 0, 7//2" half_integers(0, 7//2)

    # Whatever `half_integer` rejects is rejected here as well.  In particular a `Rational`
    # with denominator 1 is not a spelling of an integer index — the same rule the Wigner
    # constructors apply — and this is deliberate
    @test_throws "must have denominator 2" half_integer(3//1)
    @test_throws "must have denominator 2" half_integers(3//1)
    @test_throws "must have denominator 2" half_integers(3, 4//1)
    @test_throws "must have denominator 2" half_integers(1//2, 5//3)
    @test_throws ArgumentError half_integers(1.5)

    # No `Rational` survives normalization, and the kind check is settled by the argument
    # types alone: a mixed signature is inferred never to return at all
    @test Base.infer_return_type(half_integers, (Rational{Int}, Rational{Int})) === Tuple{HalfOddInteger, HalfOddInteger}
    @test Base.infer_return_type(half_integers, (Int, Int, Int)) === Tuple{Int, Int, Int}
    @test Base.infer_return_type(half_integers, (Int, Rational{Int}, Int)) === Union{}
end


@testitem "HalfOddInteger: the boundary helpers" begin
    using SphericalFunctions: HalfOddInteger, HalfInteger, IndexSpelling, unify_indices, index_value
    using SphericalFunctions: half_integer, Ysize, Yindex, L², sorted_rings, ModeWeights
    using DoubleFloats: Double64
    h(x) = HalfOddInteger(x)

    # `IndexSpelling` is exactly the set of spellings a boundary method accepts: the two kinds
    # of index, and the `Rational` in which users write a half-odd-integer
    @test IndexSpelling === Union{HalfInteger, Rational}
    for x in (3, Int8(3), big(3), h(7//2), 7//2, 7//1, big(7)//2)
        @test x isa IndexSpelling
    end
    for x in (3.5, 3.0, 7//2 + 0im, "3", nothing)
        @test !(x isa IndexSpelling)
    end

    # `unify_indices` is `half_integers` followed by `promote`: integers of different types
    # come back as one type, half-odd-integers come back as they are, and a mixture is refused
    # with the message naming both spellings
    @test unify_indices(3, 4, 5) === (3, 4, 5)
    @test unify_indices(Int8(3), 4) === (3, 4)
    @test unify_indices(Int8(-1), Int16(2), 3) === (-1, 2, 3)
    @test unify_indices(1//2, 7//2) === (h(1//2), h(7//2))
    @test unify_indices(h(1//2), 7//2, h(3//2)) === (h(1//2), h(7//2), h(3//2))
    @test unify_indices(7//2) === (h(7//2),)
    @test unify_indices(7) === (7,)
    mixed = "all be integers, like 3, or all be half-odd-integers, like 7//2"
    @test_throws mixed unify_indices(0, 7//2)
    @test_throws mixed unify_indices(h(7//2), Int8(0))
    @test_throws "must have denominator 2" unify_indices(1//2, 5//3)
    @test_throws "must have denominator 2" unify_indices(3, 4//1)
    @test Base.infer_return_type(unify_indices, (Rational{Int}, Rational{Int})) === Tuple{HalfOddInteger, HalfOddInteger}
    @test Base.infer_return_type(unify_indices, (Int8, Int)) === Tuple{Int, Int}
    @test Base.infer_return_type(unify_indices, (Int, Rational{Int})) === Union{}

    # Every `Rational` in `IndexSpelling` converts, whatever its integer type: the numerator is
    # brought to the stored `Int`, and only a numerator that does not fit an `Int` is refused,
    # with the ordinary `InexactError`
    for x in (Int8(7)//Int8(2), Int16(7)//Int16(2), Int32(7)//Int32(2), UInt8(7)//UInt8(2), big(7)//2)
        @test half_integer(x) === h(7//2)
        @test unify_indices(x, x) === (h(7//2), h(7//2))
    end
    @test half_integer(Int8(-7)//Int8(2)) === h(-7//2)
    @test half_integer(big(-1)//2) === h(-1//2)
    @test_throws InexactError half_integer((big(2)^70 + 1)//2)
    # ... and one boundary of each family accepts such a spelling and gives the `Int` result
    @test Ysize(Int8(1)//Int8(2), Int8(7)//Int8(2)) == Ysize(1//2, 7//2)
    @test Yindex(Int8(3)//Int8(2), Int8(1)//Int8(2), Int8(1)//Int8(2)) == Yindex(3//2, 1//2, 1//2)
    @test L²(big(1)//2, big(7)//2) == L²(1//2, 7//2)
    @test sorted_rings(Int32(1)//Int32(2), Int32(7)//Int32(2)) == sorted_rings(1//2, 7//2)
    @test ModeWeights(zeros(20), big(1)//2) isa ModeWeights{Float64, HalfOddInteger}  # Ysize(1//2, 7//2)
    @test ModeWeights(zeros(20), Int8(1)//Int8(2), Int8(1)//Int8(2), Int8(7)//Int8(2)) == ModeWeights(zeros(20), 1//2, 1//2, 7//2)

    # `index_value` gives an index as an element of a float type.  An integer is returned as it
    # is, for the typed comprehension that receives it to convert; a half-odd-integer becomes
    # the float exactly, through its numerator, in every float type — `Double64` included,
    # whose own constructor is ambiguous with the package's `T(::HalfOddInteger)`
    @test index_value(Float64, 3) === 3
    @test index_value(Float32, Int8(-2)) === Int8(-2)
    for T in (Float16, Float32, Float64, BigFloat, Double64), a in (-7//2, -1//2, 1//2, 3//2, 21//2)
        v = index_value(T, h(a))
        @test v isa T
        @test v == a
    end
    @test_throws MethodError Double64(h(1//2))  # the ambiguity `index_value` exists to avoid
end
