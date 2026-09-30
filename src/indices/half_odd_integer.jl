"""
    HalfOddInteger(x)

A half-odd-integer — one of ``…, -3//2, -1//2, 1//2, 3//2, …`` — stored as its odd
`numerator`, so that the represented value is `numerator//2`.  Together with the `Integer`s
these make up the [`IntegerHalf`](@ref)s, the index types over which the Wigner recurrences
are defined.

The point of the type is that the arithmetic which actually occurs in those recurrences
lands back in `Int`:

```jldoctest
julia> using SphericalFunctions: HalfOddInteger

julia> ℓ = HalfOddInteger(7//2)
7//2

julia> 2ℓ + 1  # an `Int`; compiles to `ℓ.numerator + 1`
8

julia> m = HalfOddInteger(3//2);

julia> (ℓ - m) * (ℓ + m + 1)  # an `Int`
12
```

Sums and differences of two `HalfOddInteger`s are `Int`s; adding or subtracting an integer
gives back a `HalfOddInteger`.  Multiplication is defined only by an *even* integer, and
returns an `Int`; `2ℓ` is the case that matters and it costs nothing.  Multiplying by an odd
integer throws a `DomainError`.  An integer operand of any type is converted to `Int` first,
so that one too large for an `Int` throws an `InexactError` rather than wrapping around.

A `HalfOddInteger` is ordered against integers and against other `HalfOddInteger`s, and is
`==` to the `Rational` or the floating-point number of the same value, with which it also
hashes alike.  It is never `==` to an integer, and `isinteger` is `false` for it.

The type is deliberately spare: it defines only the operations the recurrences and the
public interface actually need, and in particular it defines **no** `promote_rule`.  An
operation that has not been anticipated is therefore a loud `MethodError` rather than a
silent fall-back to `Rational` arithmetic, which is roughly thirty times slower and would
otherwise be invisible.  If you need an operation that is missing, add it deliberately.

Because a `HalfOddInteger` is never zero and never one, `zero`, `one` and `oneunit` all
throw rather than returning a value of some other type.

Users need not construct these directly: every public entry point in this package that
accepts a half-integer index also accepts a `Rational{Int}` with denominator 2, and converts
it automatically.  Use `Rational(x)` to convert back, or `Float64(x)` for a floating-point
value.  The conversion to a float type is defined for `Float16`, `Float32`, `Float64` and
`BigFloat`, and for `Double16`, `Double32` and `Double64` when DoubleFloats is loaded; any
other float type `T` converts through `T(Rational(x))`.

See also [`IntegerHalf`](@ref).
"""
struct HalfOddInteger <: Real
    numerator::Int

    # Construct from the *numerator*, unchecked.  Internal: the parity is known by
    # construction at every call site, and the public constructor below is value-semantic,
    # as Julia's numeric conventions require.
    global unsafe_half_odd_integer(numerator::Int) = new(numerator)
end

# These constructors take a *value*, not a numerator, so that `HalfOddInteger(x)` and
# `convert(HalfOddInteger, x)` agree, as Julia's `Number` machinery assumes they do.  Build
# one from a numerator with `unsafe_half_odd_integer`, or from a value with
# `HalfOddInteger(7//2)`.
HalfOddInteger(x::HalfOddInteger) = x

# The numerator is converted to the stored `Int`, so that a `Rational` of any integer type —
# `Int8(7)//Int8(2)` or `big(7)//2` — denotes the same half-odd-integer as `7//2` does.  A
# numerator too large for an `Int` fails with the ordinary `InexactError`.
#
# Any other denominator is refused with a `DomainError`, which stores the offending value
# and formats it only when the error is shown.  A message interpolated here would put
# `print_to_string` into the error path, and since effects are inferred through every
# branch, it would taint this method's effects — and those of every caller — so that no
# call could be evaluated at compile time.  As with `*` below, the error is built by a
# separate, non-inlined function, so that the inlined body is just the test and the `new`.
@inline function HalfOddInteger(x::Rational)
    denominator(x) == 2 || throw(denominator_error(x))
    unsafe_half_odd_integer(Int(numerator(x)))
end
@noinline denominator_error(x) =
    DomainError(x, "A `HalfOddInteger` must have denominator 2.")

# A whole number is not a half-odd-integer, so `HalfOddInteger(2)`, and with it
# `convert(HalfOddInteger, 2)`, is refused with an `InexactError`, as a conversion that
# cannot represent its argument is.  The public entry points check the kind of an index
# before they convert it, and refuse a whole number with a message of their own.
HalfOddInteger(x::Integer) = throw(InexactError(:HalfOddInteger, HalfOddInteger, x))

"""
    IntegerHalf

`Union{Integer, HalfOddInteger}` — the index types over which the Wigner recurrences are
defined.  A container or calculator is parameterized by one of these, so that whether its
indices are integers or half-odd-integers is a property of its *type* rather than of its
values, and the two cases dispatch to separate code.

See also [`HalfOddInteger`](@ref).
"""
const IntegerHalf = Union{Integer, HalfOddInteger}

# The smallest degree of an index type: 0 for integers and 1//2 for half-odd-integers.  The
# accessor `ℓₘᵢₙ` gives this for a block or a calculator, whose degrees are all of one index
# type, but it has no method for an index type or an index, and much of the code that needs
# this has an argument or a keyword named `ℓₘᵢₙ`, which would shadow the accessor anyway.
# The result depends on the type alone, so every call folds to a constant.
lowest_index(::Type{IT}) where {IT<:Integer} = zero(IT)
lowest_index(::Type{IT}) where {IT<:HalfOddInteger} = unsafe_half_odd_integer(1)


### Arithmetic.
#
# This is the whole point of the type: `ℓ ± m` is an `Int`, `ℓ ± 1` is a `HalfOddInteger`,
# and `2ℓ` is an `Int`, so an expression like `√((ℓ-m) * (ℓ+m+1))` — copied straight from
# the reference — is integer arithmetic throughout, with the shifts hidden in `+` and `-`
# rather than written out as `>>1` on twice-indices.  Every method is `@inline`: without
# that the check in `*` is enough to keep `2ℓ` from folding.
#
# An integer operand is converted to `Int`, the type of the stored numerator, before it is
# doubled.  For an `Int` the conversion is the identity and costs nothing; for another type
# it keeps the arithmetic in `Int`, rather than handing an unsigned or a `BigInt` sum to
# `unsafe_half_odd_integer`, and a value too large for an `Int` is an `InexactError` rather
# than a silent wrap-around.

@inline Base.:+(a::HalfOddInteger, b::HalfOddInteger) = (a.numerator + b.numerator) >> 1
@inline Base.:-(a::HalfOddInteger, b::HalfOddInteger) = (a.numerator - b.numerator) >> 1
@inline Base.:+(a::HalfOddInteger, n::Integer) = unsafe_half_odd_integer(a.numerator + 2Int(n))
@inline Base.:+(n::Integer, a::HalfOddInteger) = unsafe_half_odd_integer(a.numerator + 2Int(n))
@inline Base.:-(a::HalfOddInteger, n::Integer) = unsafe_half_odd_integer(a.numerator - 2Int(n))
@inline Base.:-(n::Integer, a::HalfOddInteger) = unsafe_half_odd_integer(2Int(n) - a.numerator)
@inline Base.:-(a::HalfOddInteger) = unsafe_half_odd_integer(-a.numerator)

# Multiply by an *even* integer, giving an `Int`.  Half-odd-integers are not closed under
# multiplication — (1/2)(1/2) = 1/4 is not one — so rather than return a type-unstable
# `Union`, this refuses an odd multiplier, which no multiplication in the recurrences has:
# their multiplier is 2, and the result is 2ℓ.  The refusal is a `DomainError` raised by a
# separate, non-inlined function, so that the inlined check is a single test of the low bit,
# which folds away for the literal 2.  One consequence is that `Base`'s `sum` of a range of
# `HalfOddInteger`s of odd length, which multiplies the first element by the length, throws;
# no code here forms such a sum.  (A docstring here would have nowhere to live in the
# manual, so this is a comment; the rule is stated in `HalfOddInteger`'s own docstring.)
@inline function Base.:*(n::Integer, a::HalfOddInteger)
    iseven(n) || throw(odd_multiplier_error(n))
    (Int(n) >> 1) * a.numerator
end
@inline Base.:*(a::HalfOddInteger, n::Integer) = n * a
@noinline odd_multiplier_error(n) =
    DomainError(n, "A `HalfOddInteger` may only be multiplied by an even integer.")

@inline Base.abs(a::HalfOddInteger) = unsafe_half_odd_integer(abs(a.numerator))

# `numerator` and `denominator` agree with what they would give for the equivalent
# `Rational`, so code that takes a half-integer apart does not have to care which of the two
# types it was handed.
@inline Base.numerator(a::HalfOddInteger) = a.numerator
@inline Base.denominator(::HalfOddInteger) = 2


### The floor of an index.
#
# `ϵ(m) = (-1)^⌊m⌋` needs the floor of an index, and the Fourier index of a half-odd `m` is
# its floor.  `floor_int` is the one expression for it that serves both index types, so that
# the code using it needs no `IT`-dependent branch; for an integer it is the identity,
# converted to `Int`.  It is an internal function, rather than a method of `Base.floor`,
# because a method of `floor` for a new argument type invalidates much of the compiled code of
# `Base` that calls it with an argument whose type is not known, including the printing of
# `@time`.
#
# The floor of a half-odd-integer is a whole number, which the type has no way to hold, so it
# is an `Int`.  For an odd numerator `a`, the arithmetic shift `a >> 1` rounds `a/2` toward -∞
# at either sign — `7 >> 1 == 3` and `-1 >> 1 == -1` — and so is the floor.  It cannot
# overflow.
@inline floor_int(n::Integer) = Int(n)
@inline floor_int(a::HalfOddInteger) = a.numerator >> 1


### Comparison.
#
# Defined against `Integer` as well as against itself, because loop bounds and validation
# compare indices to `0` and to `ℓₘᵢₙ`.  With no `promote_rule` these would otherwise fail.
# A half-odd `a` is never equal to an integer `n`, so `a < n` and `a ≤ n` are both
# `⌊a⌋ < n`, and `n < a` and `n ≤ a` are both `n ≤ ⌊a⌋`.  The floor is formed in `Int` and
# then compared with `n`, which Julia does exactly for every integer type; doubling `n`
# instead, to compare it with the numerator, would overflow for |n| > typemax(Int) ÷ 2, and
# so misorder `typemax(Int)`.

@inline Base.:(<)(a::HalfOddInteger, b::HalfOddInteger) = a.numerator < b.numerator
@inline Base.:(<)(a::HalfOddInteger, n::Integer) = floor_int(a) < n
@inline Base.:(<)(n::Integer, a::HalfOddInteger) = n ≤ floor_int(a)
@inline Base.:(<=)(a::HalfOddInteger, b::HalfOddInteger) = a.numerator <= b.numerator
@inline Base.:(<=)(a::HalfOddInteger, n::Integer) = floor_int(a) < n
@inline Base.:(<=)(n::Integer, a::HalfOddInteger) = n ≤ floor_int(a)
@inline Base.:(==)(a::HalfOddInteger, b::HalfOddInteger) = a.numerator == b.numerator
# A half-odd-integer is never equal to a whole number, and is equal to a `Rational` only if
# that `Rational` has denominator 2 — and is in reduced form, which Julia's `Rational`
# always is.  It is equal to a float `x` when `2x` is its numerator: doubling a binary float
# is exact short of overflow, which gives `Inf`, and Julia compares a float with an integer
# exactly.
@inline Base.:(==)(::HalfOddInteger, ::Integer) = false
@inline Base.:(==)(::Integer, ::HalfOddInteger) = false
@inline Base.:(==)(a::HalfOddInteger, x::Rational) =
    denominator(x) == 2 && a.numerator == numerator(x)
@inline Base.:(==)(x::Rational, a::HalfOddInteger) = a == x
@inline Base.:(==)(a::HalfOddInteger, x::AbstractFloat) = 2x == a.numerator
@inline Base.:(==)(x::AbstractFloat, a::HalfOddInteger) = a == x
@inline Base.isinteger(::HalfOddInteger) = false

# `isequal` follows `==`, so a half-odd-integer that is `isequal` to a `Rational` or a float
# must hash as that value does, or a `Dict` or `Set` keyed by one would not find the other.
# Hashing the `Rational` itself makes the two agree by construction, and agree with the hash
# of the equal float too, because `Base` makes the hashes of `Rational`s and floats agree.
Base.hash(a::HalfOddInteger, h::UInt) = hash(Rational(a), h)


### Values that do not exist.
#
# `Base`'s `Number` fall-backs would form these by conversion — `one(::Type{T})` is
# `convert(T, 1)` — which the value-semantic constructor refuses with an `InexactError` about
# a conversion the caller never asked for.  These methods say what is wrong instead.  They
# must throw, whatever they say: a zero or a one of some other type would propagate silently
# into range machinery, which forms a step as `oneunit(T) - zero(T)`.

Base.zero(::Type{HalfOddInteger}) =
    throw(ArgumentError("`HalfOddInteger` has no zero: 0 is not a half-odd-integer."))
Base.one(::Type{HalfOddInteger}) =
    throw(ArgumentError("`HalfOddInteger` has no one: 1 is not a half-odd-integer."))
Base.oneunit(::Type{HalfOddInteger}) =
    throw(ArgumentError("`HalfOddInteger` has no oneunit: 1 is not a half-odd-integer."))
Base.zero(::HalfOddInteger) = zero(HalfOddInteger)
Base.one(::HalfOddInteger) = one(HalfOddInteger)
Base.oneunit(::HalfOddInteger) = oneunit(HalfOddInteger)

# More informative than the `MethodError` that would otherwise arise, and it is what the
# `Integer` path's `Int(m - mₘᵢₙ)` offsets would hit if one were ever reached with the wrong
# index type.
Base.Int(a::HalfOddInteger) = throw(InexactError(:Int, Int, a))


### Ranges.
#
# `m′ₘᵢₙ(w):m′ₘₐₓ(w)` builds a `UnitRange{HalfOddInteger}` without help — `unitrange_last`
# forms only `HalfOddInteger ± Int` expressions, from the `Int` difference of the endpoints.
# But `Base`'s `length` for a non-`Integer` unit range reaches for `zero(T)`, and `step` is
# defined as `oneunit(T) - zero(T)`, so both must be given directly.  With these two, the
# ranges iterate, `collect`, `sort` and `in` all work, allocation-free.  A descending
# `m′ₘₐₓ:-1:m′ₘᵢₙ`, which `reverse` also gives, is a `StepRange{HalfOddInteger, Int}`; it
# iterates without help, and its membership test is given below.

Base.length(r::UnitRange{HalfOddInteger}) = max(0, (last(r) - first(r)) + 1)
Base.step(::UnitRange{HalfOddInteger}) = 1

# Membership means `==` to some element, as it does for `Base`'s ranges: `3//2 ∈ 1//2:5//2`
# is true however the `3//2` is spelled.  `Base`'s `in(::Real, ::AbstractRange)` would
# promote the value to the range's type, which `HalfOddInteger` deliberately cannot do, so
# the test is written out: a value is a member when twice it is an odd integer between the
# numerators of the endpoints, and, for a range whose step is not 1, when it is a whole number
# of steps from the first element.  (A whole number never is, nor anything that is not a
# multiple of 1/2.)  The methods whose value is a `HalfOddInteger` settle the tie between the
# ones whose value is any `Real` and `Base`'s `in(::T, ::AbstractRange{T})`.
@inline Base.in(x::HalfOddInteger, r::UnitRange{HalfOddInteger}) = first(r) ≤ x ≤ last(r)
@inline Base.in(x::Real, r::UnitRange{HalfOddInteger}) = in_half_odd_range(x, first(r), last(r))
@inline Base.in(x::HalfOddInteger, r::StepRange{HalfOddInteger, <:Integer}) =
    in_half_odd_range(x, r)
@inline Base.in(x::Real, r::StepRange{HalfOddInteger, <:Integer}) = in_half_odd_range(x, r)
@inline function in_half_odd_range(x::Real, lo::HalfOddInteger, hi::HalfOddInteger)
    t = 2x
    isinteger(t) && !isinteger(x) && numerator(lo) ≤ t ≤ numerator(hi)
end
@inline function in_half_odd_range(x::Real, r::StepRange{HalfOddInteger, <:Integer})
    isempty(r) && return false
    t = 2x
    a, b = numerator(first(r)), numerator(last(r))
    lo, hi = minmax(a, b)
    # `t` is a whole number between two `Int`s once the tests before the last have passed, so
    # `Int(t)` is exact.
    isinteger(t) && !isinteger(x) && lo ≤ t ≤ hi && iszero(rem(Int(t) - a, 2step(r)))
end


### Conversion and display.

Base.Rational{T}(a::HalfOddInteger) where {T<:Integer} = T(a.numerator) // T(2)
Base.Rational(a::HalfOddInteger) = a.numerator // 2

# The conversion to a float is defined type by type, rather than once for every
# `T<:AbstractFloat`, because such a method is exactly as specific as a float package's own
# `T(x::Real)` constructor, and the two are then ambiguous for every half-odd-integer: this
# is so for `DoubleFloats`' `Double64(x::T) where {T<:Real}`.  The extension
# `SphericalFunctionsDoubleFloatsExt` defines the conversions to `Double16`, `Double32` and
# `Double64`.
for F ∈ (Float16, Float32, Float64, BigFloat)
    @eval (::Type{$F})(a::HalfOddInteger) = $F(a.numerator) / 2
end
Base.AbstractFloat(a::HalfOddInteger) = Float64(a)
Base.float(a::HalfOddInteger) = Float64(a)

# An index as an element of a floating-point matrix or vector of type `T`.  An integer is
# returned as is, for the typed comprehension that receives it to convert as it always has.
# A half-odd-integer is converted through its numerator: the `Int` 2x becomes a `T` and is
# halved, which is exact in every binary floating-point type.  This exists because the
# direct conversion `T(x)` is defined only for the float types named above, whereas the
# numerator route depends on nothing but `T(::Int)`, which every number type has.
@inline index_value(::Type{T}, x::Integer) where {T} = x
@inline index_value(::Type{T}, x::HalfOddInteger) where {T} = T(2x) / 2

# Displayed as `5//2` rather than `5/2`: that is the notation users write at every entry
# point, it round-trips through `HalfOddInteger(5//2)`, and it cannot be misread as a
# floating-point division.  The type itself is named in `summary`, so nothing is hidden.
Base.show(io::IO, a::HalfOddInteger) = print(io, a.numerator, "//2")
