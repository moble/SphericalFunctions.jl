# Closed-form indexing into the canonical ordering of mode weights,
#
#     [ f(ℓ, m) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ ],
#
# i.e., ℓ-major with m increasing.  These formulas (and their SymPy derivations) are inherited
# from the v2 code; they are the only place in the package that knows this layout.
#
# Each function has three kinds of method.  The `Integer` methods are the v2 formulas, and are
# generic over every integer type.  The `HalfOddInteger` methods are the same formulas written
# on the doubled indices 2ℓ, 2m and 2ℓₘᵢₙ, which are odd `Int`s, so that the arithmetic never
# leaves `Int`: the quarter-integers that would appear in ℓ(ℓ+1) and ℓₘᵢₙ² cancel against each
# other, and the final division by 4 is exact.  The boundary methods accept `Integer`,
# `HalfOddInteger` and `Rational` indices alike, normalize them with `half_integers`, and
# re-dispatch; they are also what refuses a call such as `Ysize(0, 7//2)`, which mixes the
# two kinds of index, with an explanation rather than a bare `MethodError`.

"""
    Ysize(ℓₘₐₓ)
    Ysize(ℓₘᵢₙ, ℓₘₐₓ)

Total number of mode weights ``(ℓ, m)`` with ``ℓₘᵢₙ ≤ ℓ ≤ ℓₘₐₓ`` and ``-ℓ ≤ m ≤ ℓ``, which is
``(ℓₘₐₓ+1)^2 - ℓₘᵢₙ^2``.

The indices may be integers or half-odd-integers, the latter passed either as `Rational`s
with denominator 2 (`7//2`) or as [`HalfOddInteger`](@ref)s.  Both indices in one call must be
of the same kind; a call that mixes them, such as `Ysize(0, 7//2)`, throws an `ArgumentError`
saying so.  The formula holds for either kind — for half-odd ``ℓₘᵢₙ`` and ``ℓₘₐₓ`` the quarters
in the two squares cancel — and the result is always a whole number.  `ℓₘᵢₙ` defaults to the
smallest ``ℓ`` of the given kind, which is 0 for integers and 1/2 for half-odd-integers.  An
`ArgumentError` is thrown if ``ℓₘᵢₙ < 0`` or if ``ℓₘₐₓ < ℓₘᵢₙ - 1`` (the empty range
``ℓₘₐₓ = ℓₘᵢₙ - 1`` has size 0).

See also [`Yindex`](@ref) and [`Yrange`](@ref).
"""
function Ysize(ℓₘᵢₙ::Integer, ℓₘₐₓ::Integer)
    if ℓₘᵢₙ < 0
        throw(ArgumentError("ℓₘᵢₙ=$ℓₘᵢₙ must be non-negative."))
    end
    if ℓₘₐₓ < ℓₘᵢₙ - 1
        throw(ArgumentError("ℓₘₐₓ=$ℓₘₐₓ must be at least ℓₘᵢₙ-1=$(ℓₘᵢₙ-1)."))
    end
    (ℓₘₐₓ + 1)^2 - ℓₘᵢₙ^2
end
function Ysize(ℓₘᵢₙ::HalfOddInteger, ℓₘₐₓ::HalfOddInteger)
    if ℓₘᵢₙ < 0
        throw(ArgumentError("ℓₘᵢₙ=$ℓₘᵢₙ must be non-negative."))
    end
    if ℓₘₐₓ < ℓₘᵢₙ - 1
        throw(ArgumentError("ℓₘₐₓ=$ℓₘₐₓ must be at least ℓₘᵢₙ-1=$(ℓₘᵢₙ-1)."))
    end
    # (ℓₘₐₓ+1)² - ℓₘᵢₙ² on the doubled indices, which are odd `Int`s.  The difference of the
    # two squares is a multiple of 4, so the division is exact and the result is an `Int`.
    ((2ℓₘₐₓ + 2)^2 - (2ℓₘᵢₙ)^2) ÷ 4
end
# The one-argument form starts the ordering at the floor of the index type: 0 for an integer
# ℓₘₐₓ, and 1/2 for a half-odd one.
Ysize(ℓₘₐₓ::IT) where {IT<:IntegerHalf} = Ysize(ℓₘᵢₙ(IT), ℓₘₐₓ)
# The boundary methods, which are the only ones that see a `Rational`.  `IndexArgument` is
# defined beside `IntegerHalf` in `half_odd_integer.jl`.
Ysize(ℓₘₐₓ::Rational) = Ysize(half_integer(ℓₘₐₓ))
function Ysize(ℓₘᵢₙ::IndexArgument, ℓₘₐₓ::IndexArgument)
    Ysize(half_integers(ℓₘᵢₙ, ℓₘₐₓ)...)
end

"""
    Yindex(ℓ, m)
    Yindex(ℓ, m, ℓₘᵢₙ)

Index of the mode weight ``(ℓ, m)`` in the canonical ordering
`[f(ℓ, m) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ]`, which is ``ℓ(ℓ+1) - ℓₘᵢₙ^2 + m + 1``.

As for [`Ysize`](@ref), the indices may be integers or half-odd-integers, passed as
`Rational`s with denominator 2 or as [`HalfOddInteger`](@ref)s, and all of the indices in one
call must be of the same kind.  The formula gives a whole number for either kind, because for
half-odd indices the quarters in ``ℓ(ℓ+1)`` and ``ℓₘᵢₙ^2`` cancel.  `ℓₘᵢₙ` defaults to the
smallest ``ℓ`` of the given kind, which is 0 for integers and 1/2 for half-odd-integers.  No
bounds are checked.

See also [`Ysize`](@ref) and [`Yrange`](@ref).
"""
@inline Yindex(ℓ::Integer, m::Integer, ℓₘᵢₙ::Integer) = ℓ*(ℓ+1) - ℓₘᵢₙ^2 + m + 1
@inline Yindex(ℓ::Integer, m::Integer) = Yindex(ℓ, m, ℓₘᵢₙ(typeof(ℓ)))
# ℓ(ℓ+1) - ℓₘᵢₙ² + m + 1 on the doubled indices a = 2ℓ, b = 2m and c = 2ℓₘᵢₙ, all odd `Int`s.
# Four times the index is a(a+2) - c² + 2b + 4; for odd a, b and c this is a multiple of 4, so
# the division is exact and the result is an `Int`.
@inline function Yindex(ℓ::HalfOddInteger, m::HalfOddInteger, ℓₘᵢₙ::HalfOddInteger)
    a, b, c = 2ℓ, 2m, 2ℓₘᵢₙ
    (a*(a+2) - c^2 + 2b + 4) ÷ 4
end
@inline Yindex(ℓ::HalfOddInteger, m::HalfOddInteger) = Yindex(ℓ, m, ℓₘᵢₙ(HalfOddInteger))
# The boundary methods, which are the only ones that see a `Rational`.
@inline function Yindex(ℓ::IndexArgument, m::IndexArgument)
    Yindex(half_integers(ℓ, m)...)
end
@inline function Yindex(ℓ::IndexArgument, m::IndexArgument, ℓₘᵢₙ::IndexArgument)
    Yindex(half_integers(ℓ, m, ℓₘᵢₙ)...)
end

"""
    Yrange(ℓₘₐₓ)
    Yrange(ℓₘᵢₙ, ℓₘₐₓ)

Vector of the ``(ℓ, m)`` pairs in the canonical ordering, so that `Yrange(ℓₘᵢₙ, ℓₘₐₓ)[i]` is
the pair stored at index `i`.

As for [`Ysize`](@ref), the indices may be integers or half-odd-integers, passed as
`Rational`s with denominator 2 or as [`HalfOddInteger`](@ref)s, and both indices in one call
must be of the same kind.  For half-odd-integers the pairs are of `HalfOddInteger`s —
whatever the type of the arguments — as is every other index value the package hands back;
use `Rational(x)` to convert an element to the more familiar `Rational`.  `ℓₘᵢₙ` defaults to
the smallest ``ℓ`` of the given kind, which is 0 for integers and 1/2 for half-odd-integers.

See also [`Ysize`](@ref) and [`Yindex`](@ref).
"""
function Yrange(ℓₘᵢₙ::Integer, ℓₘₐₓ::Integer)
    [(ℓ, m) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ]
end
# `UnitRange{HalfOddInteger}` iterates, so the comprehension is the same one.
function Yrange(ℓₘᵢₙ::HalfOddInteger, ℓₘₐₓ::HalfOddInteger)
    [(ℓ, m) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ]
end
Yrange(ℓₘₐₓ::IT) where {IT<:IntegerHalf} = Yrange(ℓₘᵢₙ(IT), ℓₘₐₓ)
# The boundary methods, which are the only ones that see a `Rational`.
Yrange(ℓₘₐₓ::Rational) = Yrange(half_integer(ℓₘₐₓ))
function Yrange(ℓₘᵢₙ::IndexArgument, ℓₘₐₓ::IndexArgument)
    Yrange(half_integers(ℓₘᵢₙ, ℓₘₐₓ)...)
end
