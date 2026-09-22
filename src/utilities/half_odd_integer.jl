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

Sums and differences of two `HalfOddInteger`s are `Integer`s; adding or subtracting an
`Integer` gives back a `HalfOddInteger`.  Multiplication is defined only by an *even*
integer, and returns an `Integer`; `2ℓ` is the case that matters and it costs nothing.

The type is deliberately spare: it defines only the operations the recurrences and the
public interface actually need, and in particular it defines **no** `promote_rule`.  An
operation that has not been anticipated is therefore a loud `MethodError` rather than a
silent fall-back to `Rational` arithmetic, which is roughly thirty times slower and would
otherwise be invisible.  If you need an operation that is missing, add it deliberately.

Because a `HalfOddInteger` is never zero and never one, `zero`, `one` and `oneunit` all
throw rather than returning a value of some other type.

Users need not construct these directly: every public entry point that accepts a
half-integer index also accepts a `Rational` with denominator 2, and converts.  Use
`Rational(x)` to convert back.

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
# `convert(HalfOddInteger, x)` agree — Julia's `Number` machinery assumes they do, and
# `convert` is what the calculators call on a user-supplied `ℓ`.  Build one from a numerator
# with `unsafe_half_odd_integer`, or from a value with `HalfOddInteger(7//2)`.
HalfOddInteger(x::HalfOddInteger) = x

# The numerator is converted to the stored `Int`, so that a `Rational` of any integer type —
# `Int8(7)//Int8(2)` or `big(7)//2` — denotes the same half-odd-integer as `7//2` does.  A
# numerator too large for an `Int` fails with the ordinary `InexactError`.
function HalfOddInteger(x::Rational)
    denominator(x) == 2 || throw(ArgumentError(
        "A `HalfOddInteger` must have denominator 2; got $x."
    ))
    unsafe_half_odd_integer(Int(numerator(x)))
end

# A whole number is not a half-odd-integer.  This is the error a user sees on, for example,
# `recurrence!(calc, R, 2)` for a half-integer calculator.
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

# The types in which an index may be passed to a boundary method: an `Integer`, a
# `HalfOddInteger`, or the `Rational` with denominator 2 in which users write a
# half-odd-integer.  A method whose index arguments are typed with this normalizes them,
# through `half_integers` or `unify_indices` below, may check their kind against a
# container's, and then re-dispatches to a method that sees only `Integer` or
# `HalfOddInteger` values, so that no method body ever sees a `Rational`.  (A docstring here
# would have to be placed in the manual, where the alias would mean nothing to a reader; the
# rule it encodes is stated under `half_integers`.)
const IndexArgument = Union{IntegerHalf, Rational}


### Arithmetic.
#
# This is the whole point of the type: `ℓ ± m` is an `Int`, `ℓ ± 1` is a `HalfOddInteger`,
# and `2ℓ` is an `Int`, so an expression like `√((ℓ-m) * (ℓ+m+1))` — copied straight from
# the reference — is integer arithmetic throughout, with the shifts hidden in `+` and `-`
# rather than written out as `>>1` on twice-indices.  Every method is `@inline`: without
# that the assertion in `*` is enough to keep `2ℓ` from folding.

@inline Base.:+(a::HalfOddInteger, b::HalfOddInteger) = (a.numerator + b.numerator) >> 1
@inline Base.:-(a::HalfOddInteger, b::HalfOddInteger) = (a.numerator - b.numerator) >> 1
@inline Base.:+(a::HalfOddInteger, n::Integer) = unsafe_half_odd_integer(a.numerator + 2n)
@inline Base.:+(n::Integer, a::HalfOddInteger) = unsafe_half_odd_integer(a.numerator + 2n)
@inline Base.:-(a::HalfOddInteger, n::Integer) = unsafe_half_odd_integer(a.numerator - 2n)
@inline Base.:-(n::Integer, a::HalfOddInteger) = unsafe_half_odd_integer(2n - a.numerator)
@inline Base.:-(a::HalfOddInteger) = unsafe_half_odd_integer(-a.numerator)

# Multiply by an *even* integer, giving an `Integer`.  Half-odd-integers are not closed
# under multiplication — (1/2)(1/2) = 1/4 is not one — so rather than return a type-unstable
# `Union`, this asserts what is true of every multiplication the recurrences actually
# perform: the multiplier is 2, and the result is 2ℓ.  (A docstring here would have nowhere
# to live in the manual, so this is a comment; the rule is stated in `HalfOddInteger`'s own
# docstring.)
@inline function Base.:*(n::Integer, a::HalfOddInteger)
    @assert iseven(n) "A `HalfOddInteger` may only be multiplied by an even integer; got $n."
    (n >> 1) * a.numerator
end
@inline Base.:*(a::HalfOddInteger, n::Integer) = n * a

@inline Base.abs(a::HalfOddInteger) = unsafe_half_odd_integer(abs(a.numerator))

# `numerator` and `denominator` agree with what they would give for the equivalent
# `Rational`, so code that takes a half-integer apart does not have to care which of the two
# types it was handed.
@inline Base.numerator(a::HalfOddInteger) = a.numerator
@inline Base.denominator(::HalfOddInteger) = 2


### Comparison.
#
# Defined against `Integer` as well as against itself, because loop bounds and validation
# compare indices to `0` and to `ℓₘᵢₙ`.  With no `promote_rule` these would otherwise fail.

@inline Base.:(<)(a::HalfOddInteger, b::HalfOddInteger) = a.numerator < b.numerator
@inline Base.:(<)(a::HalfOddInteger, n::Integer) = a.numerator < 2n
@inline Base.:(<)(n::Integer, a::HalfOddInteger) = 2n < a.numerator
@inline Base.:(<=)(a::HalfOddInteger, b::HalfOddInteger) = a.numerator <= b.numerator
@inline Base.:(<=)(a::HalfOddInteger, n::Integer) = a.numerator <= 2n
@inline Base.:(<=)(n::Integer, a::HalfOddInteger) = 2n <= a.numerator
@inline Base.:(==)(a::HalfOddInteger, b::HalfOddInteger) = a.numerator == b.numerator
# A half-odd-integer is never equal to a whole number, and is equal to a `Rational` only if
# that `Rational` has denominator 2 — and is in reduced form, which Julia's `Rational`
# always is.
@inline Base.:(==)(::HalfOddInteger, ::Integer) = false
@inline Base.:(==)(::Integer, ::HalfOddInteger) = false
@inline Base.:(==)(a::HalfOddInteger, x::Rational) = denominator(x) == 2 && a.numerator == numerator(x)
@inline Base.:(==)(x::Rational, a::HalfOddInteger) = a == x

# `ε(m) = (-1)^⌊m⌋` needs the floor, and this is the one expression that serves both index
# types, so that the recurrences need no `IT`-dependent branch for it.
@inline Base.floor(::Type{T}, a::HalfOddInteger) where {T<:Integer} = T((a.numerator - 1) >> 1)
@inline Base.floor(a::HalfOddInteger) = floor(Int, a)

# `sorted_rings` uses the spin weight as a count of ulps to break the ties in its sort, and
# `ceil(s)` is the one expression for that count which serves both index types, since for an
# `Integer` it is the identity.  Like `floor`, this returns an `Int` rather than a
# `HalfOddInteger`: the ceiling of a half-odd-integer is a whole number, and the type has no
# way to hold one.  For an odd numerator `a`, the ceiling of `a/2` is exactly `(a+1)/2` at
# either sign — ⌈1/2⌉ = 1 = (1+1)/2 and ⌈-1/2⌉ = 0 = (-1+1)/2 — and because `a+1` is even
# the shift halves it exactly, so no rounding direction enters and the one expression is
# right for negative and positive values alike.
@inline Base.ceil(::Type{T}, a::HalfOddInteger) where {T<:Integer} = T((a.numerator + 1) >> 1)
@inline Base.ceil(a::HalfOddInteger) = ceil(Int, a)


### Values that do not exist.
#
# `Base`'s `Number` fall-backs would otherwise manufacture these from the positional
# constructor — `one(::Type{T})` is `convert(T, 1)`, which would return 1/2 — and the wrong
# values then propagate silently into range machinery as a zero step.  These must throw.

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
# only ever forms `start + floor(stop - start)`, which is `HalfOddInteger + Int`.  But
# `Base`'s `length` for a non-`Integer` unit range reaches for `zero(T)`, and `step` is
# defined as `oneunit(T) - zero(T)`, so both must be given directly.  With these two, the
# ranges iterate, `collect`, `sort` and `in` all work, allocation-free.  A descending
# `m′ₘₐₓ:-1:m′ₘᵢₙ` is a `StepRange{HalfOddInteger, Int}` and needs nothing.

Base.length(r::UnitRange{HalfOddInteger}) = max(0, (last(r) - first(r)) + 1)
Base.step(::UnitRange{HalfOddInteger}) = 1
Base.step(::Type{UnitRange{HalfOddInteger}}) = 1


### Conversion and display.

Base.Rational{T}(a::HalfOddInteger) where {T<:Integer} = T(a.numerator) // T(2)
Base.Rational(a::HalfOddInteger) = a.numerator // 2
(::Type{T})(a::HalfOddInteger) where {T<:AbstractFloat} = T(a.numerator) / 2
Base.AbstractFloat(a::HalfOddInteger) = Float64(a)
Base.float(a::HalfOddInteger) = Float64(a)

# An index as an element of a floating-point matrix or vector of type `T`.  An integer is
# returned as is, for the typed comprehension that receives it to convert as it always has.
# A half-odd-integer is converted through its numerator: the `Int` 2x becomes a `T` and is
# halved, which is exact in every binary floating-point type.  This exists because the
# direct conversion `T(x)` is not resolved by every float type — `DoubleFloats` defines
# `Double64(x::T) where {T<:Real}`, which is exactly as specific as the
# `(::Type{T})(a::HalfOddInteger) where {T<:AbstractFloat}` above, so `Double64(x)` is an
# ambiguity error — whereas the numerator route depends on nothing but `T(::Int)`.
@inline index_value(::Type{T}, x::Integer) where {T} = x
@inline index_value(::Type{T}, x::HalfOddInteger) where {T} = T(2x) / 2

# Displayed as `5//2` rather than `5/2`: that is the notation users write at every entry
# point, it round-trips through `HalfOddInteger(5//2)`, and it cannot be misread as a
# floating-point division.  The type itself is named in `summary`, so nothing is hidden.
Base.show(io::IO, a::HalfOddInteger) = print(io, a.numerator, "//2")


"""
    half_integer(x)

Normalize a user-supplied index to the package's internal index type: an `Integer` is
returned unchanged, and a `Rational` with denominator 2 becomes a [`HalfOddInteger`](@ref).
Anything else is an error.

Every public entry point that takes an index calls this, so that callers may keep writing
`3//2` while the recurrences see a `HalfOddInteger`.  Use `Rational(x)` to convert back.
"""
@inline half_integer(x::Integer) = x
@inline half_integer(x::HalfOddInteger) = x
@inline function half_integer(x::Rational)
    # The refusal names both types, which the constructor's own message — written for a
    # caller who asked for a `HalfOddInteger` by name — does not; the constructor then does
    # the conversion, including the `InexactError` for a numerator too large for an `Int`.
    denominator(x) == 2 || throw(ArgumentError(
        "A `Rational` index must have denominator 2, such as 7//2;\n"
        * "an integer index is given as an `Integer`, such as 3.  Got $x."
    ))
    HalfOddInteger(x)
end
half_integer(x) = throw(ArgumentError(
    "An index must be an `Integer`, a `HalfOddInteger`, or a `Rational` with denominator 2;"
    * " got $x of type $(typeof(x))."
))

"""
    half_integers(xs...)

Normalize several user-supplied indices at once, returning them as a tuple.  Each is passed
through [`half_integer`](@ref), and the results are then required to be all of one kind —
all `Integer`s or all [`HalfOddInteger`](@ref)s — because the package never mixes the two
kinds of index within a single call: an integer ``ℓ`` goes with an integer ``m`` and
``ℓₘᵢₙ``, and a half-odd ``ℓ`` with half-odd ones.  A call that mixes them, such as
`Ysize(0, 7//2)`, is refused with an `ArgumentError` that names both kinds, in place of
the bare `MethodError` that dispatch alone would produce.

This is the tool of the boundary methods: a method that accepts `Rational` arguments calls
this on its index arguments and re-dispatches on the result, so that no method body ever
sees a `Rational`.  A `Rational` with denominator 1, such as `3//1`, is refused by
`half_integer` like any other `Rational` whose denominator is not 2; an integer index is
passed as an `Integer`.
"""
@inline function half_integers(xs...)
    ys = map(half_integer, xs)
    if !(all(y -> y isa Integer, ys) || all(y -> y isa HalfOddInteger, ys))
        throw(ArgumentError(
            "The indices in one call must all be integers, like 3, or all be "
            * "half-odd-integers, like 7//2; got " * join(xs, ", ") * "."
        ))
    end
    ys
end

# Normalize the indices that describe a set of mode weights — the spin weight, ℓₘᵢₙ and ℓₘₐₓ
# — to one concrete index type.  `half_integers` turns each `Rational` into a
# `HalfOddInteger` and refuses a mixture of integers and half-odd-integers with an
# explanation; `promote` then unifies integers of different concrete types, as the callers'
# own arithmetic would have done before half-integer indices were admitted, and returns
# `HalfOddInteger`s unchanged, since they are already of one type (the type has no
# `promote_rule`, and none is needed here).  This is what the boundary methods of the flat
# `sYlm` functions, the `ModeWeights` constructors, the operators and the pixelizations call
# before re-dispatching to their workers.  (A docstring here would have to be placed in the
# manual; the rule is stated under `half_integers`.)
@inline unify_indices(xs...) = promote(half_integers(xs...)...)

# Normalize the pair of indices an `sYlmCalculator` is built from, where the spin weight may
# be either a single value or an ascending range of them.  The two forms are deliberately
# kept apart: a scalar argument normalizes to a scalar and a range to a `UnitRange`, because
# the calculator stores whichever it was given and that is what settles the shape of the
# block handed back by `recurrence!`.
#
# A range is normalized from its two endpoints rather than element by element, and its step
# is compared by value rather than converted.  The reason is that `-3//2:3//2` is a
# `UnitRange{Rational{Int}}` whose step is `1//1`, which `half_integer` rightly refuses; the
# comparison `step(s) == 1` is true for that type and for the `UnitRange{Int}` and
# `UnitRange{HalfOddInteger}` types alike.  The result is rebuilt with the colon, which
# forms only `start + floor(stop - start)` and so asks `HalfOddInteger` for nothing it
# lacks.
@inline spin_indices(ℓₘₐₓ, s) = unify_indices(ℓₘₐₓ, s)
function spin_indices(ℓₘₐₓ, s::AbstractRange)
    # A descending range is answered before the step is complained about in general, because
    # the two ways to produce one — `3//2:-1:-3//2` and the empty `3//2:-3//2` — are
    # a single mistake with a single remedy, and naming that remedy is more use than naming
    # the step.
    if step(s) < 0 || isempty(s)
        throw(ArgumentError(
            "The range of spin weights $s runs downward or is empty.  A range runs from its "
            * "lower limit to its upper one, so 3//2:-1:-3//2 is written -3//2:3//2."
        ))
    end
    if step(s) != 1
        throw(ArgumentError(
            "The spin weights of one calculator are consecutive, so the range's step must be "
            * "1; got $s, whose step is $(step(s))."
        ))
    end
    ℓₘₐₓ, lo, hi = unify_indices(ℓₘₐₓ, first(s), last(s))
    (ℓₘₐₓ, lo:hi)
end

# Convert the index-valued keyword arguments of a public entry point.  Only the names that
# are indices are touched; anything else (`Nᵣ`, say) is passed through untouched.
const INDEX_KEYWORDS = (
    :m′ₘₐₓ, :m′ₘᵢₙ, :mₘₐₓ, :mₘᵢₙ, :mp_max, :mp_min, :m_max, :m_min,
    # `SpinMatrix` and `SpinMatrixBatch` take their spin bounds by keyword, and those
    # are indices of the same kind as the rest; without them the `Rational`-ℓ
    # constructors hand an unconverted `Rational` to a method typed on `IT`.
    :sₘₐₓ, :sₘᵢₙ, :s_max, :s_min,
)
function half_integer_kwargs(kwargs)
    pairs(NamedTuple(
        k => (k in INDEX_KEYWORDS ? half_integer(v) : v) for (k, v) in pairs(kwargs)
    ))
end

"""
    isindex(IT, x)

Whether `x` denotes a legal index of type `IT` — an `Integer` for an integer `IT`, or a
half-odd-integer (however written) for `IT === HalfOddInteger`.

This exists so that a container can give a helpful message about its own index set before
`convert` throws a bare `InexactError` about the type.
"""
@inline isindex(::Type{<:Integer}, x::Integer) = true
@inline isindex(::Type{<:Integer}, x::Rational) = isinteger(x)
@inline isindex(::Type{<:Integer}, x) = false
@inline isindex(::Type{HalfOddInteger}, x::HalfOddInteger) = true
@inline isindex(::Type{HalfOddInteger}, x::Rational) = denominator(x) == 2
@inline isindex(::Type{HalfOddInteger}, x) = false
