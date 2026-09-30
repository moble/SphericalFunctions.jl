# Closed-form indexing into the canonical ordering of mode weights,
#
#     [ f(ℓ, m) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ ],
#
# i.e., ℓ-major with m increasing.  These formulas are the only place in the package that
# relies on this layout.
#
# Each function is defined with `@index_methods`, so that its body sees indices that are all
# `Int`s or all `HalfOddInteger`s, and a call that mixes the two kinds, such as `Ysize(0,
# 7//2)`, or that passes an integer of another type, such as an `Int8`, is rejected.  The
# body branches on the index type, which is settled when the method is compiled.  For `Int`s
# it is the formula of the docstring, evaluated as written.  For `HalfOddInteger`s it is the
# same formula written on the doubled indices 2ℓ, 2m and 2ℓₘᵢₙ, which are odd `Int`s, so
# that the arithmetic never leaves `Int`: the quarter-integers that would appear in ℓ(ℓ+1)
# and ℓₘᵢₙ² cancel against each other, and the final division by 4 is exact.  (The doubled
# form is exact for integers too, but the branch keeps the integer arithmetic as short as it
# can be.)

"""
    Ysize(ℓₘₐₓ)
    Ysize(ℓₘᵢₙ, ℓₘₐₓ)

Total number of mode weights ``(ℓ, m)`` with ``ℓₘᵢₙ ≤ ℓ ≤ ℓₘₐₓ`` and ``-ℓ ≤ m ≤ ℓ``, which
is ``(ℓₘₐₓ+1)^2 - ℓₘᵢₙ^2``.

The indices may be integers of type `Int`, or half-odd-integers passed either as
`Rational{Int}`s with denominator 2 (`7//2`) or as [`HalfOddInteger`](@ref)s.  Both indices
in one call must be of the same kind; a call that mixes them, such as `Ysize(0, 7//2)`, or
that passes an integer of another type, such as an `Int8`, throws an `ArgumentError` saying
so.  The formula holds for either kind — for half-odd ``ℓₘᵢₙ`` and ``ℓₘₐₓ`` the quarters in
the two squares cancel — and the result is always a whole number.  `ℓₘᵢₙ` defaults to the
smallest ``ℓ`` of the given kind, which is 0 for integers and 1/2 for half-odd-integers.  An
`ArgumentError` is thrown if ``ℓₘᵢₙ < 0`` or if ``ℓₘₐₓ < ℓₘᵢₙ - 1`` (the empty range ``ℓₘₐₓ
= ℓₘᵢₙ - 1`` has size 0).

See also [`Yindex`](@ref) and [`Yrange`](@ref).
"""
function Ysize end

@index_methods function Ysize(ℓₘᵢₙ::IT, ℓₘₐₓ::IT) where {IT<:IndexType}
    if ℓₘᵢₙ < 0
        throw(ArgumentError("ℓₘᵢₙ=$ℓₘᵢₙ must be non-negative."))
    end
    if ℓₘₐₓ < ℓₘᵢₙ - 1
        throw(ArgumentError("ℓₘₐₓ=$ℓₘₐₓ must be at least ℓₘᵢₙ-1=$(ℓₘᵢₙ-1)."))
    end
    if IT === Int
        (ℓₘₐₓ + 1)^2 - ℓₘᵢₙ^2
    else
        # (ℓₘₐₓ+1)² - ℓₘᵢₙ² on the doubled indices, which are odd `Int`s.  The difference of
        # the two squares is a multiple of 4, so the division is exact.
        ((2ℓₘₐₓ + 2)^2 - (2ℓₘᵢₙ)^2) ÷ 4
    end
end
# The one-argument form starts the ordering at the floor of the index type: 0 for an integer
# ℓₘₐₓ, and 1/2 for a half-odd one.
@index_methods Ysize(ℓₘₐₓ::IT) where {IT<:IndexType} = Ysize(ℓₘᵢₙ(IT), ℓₘₐₓ)

"""
    Yindex(ℓ, m)
    Yindex(ℓ, m, ℓₘᵢₙ)

Index of the mode weight ``(ℓ, m)`` in the canonical ordering `[f(ℓ, m) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ
for m ∈ -ℓ:ℓ]`, which is ``ℓ(ℓ+1) - ℓₘᵢₙ^2 + m + 1``.

As for [`Ysize`](@ref), the indices may be integers of type `Int`, or half-odd-integers
passed as `Rational{Int}`s with denominator 2 or as [`HalfOddInteger`](@ref)s, and all of
the indices in one call must be of the same kind.  The formula gives a whole number for
either kind, because for half-odd indices the quarters in ``ℓ(ℓ+1)`` and ``ℓₘᵢₙ^2`` cancel.
`ℓₘᵢₙ` defaults to the smallest ``ℓ`` of the given kind, which is 0 for integers and 1/2 for
half-odd-integers.  No bounds are checked.

See also [`Ysize`](@ref) and [`Yrange`](@ref).
"""
function Yindex end

# `Yindex` is called in the inner loops of the package with `Int` or `HalfOddInteger`
# indices, which reach the work methods directly, and `@inline` applies to every generated
# method.  The default of `ℓₘᵢₙ` is evaluated after the indices have been converted, so it may
# use `IT`; the call without `ℓₘᵢₙ` has its own methods, generated from the default.
@index_methods @inline function Yindex(ℓ::IT, m::IT, ℓₘᵢₙ::IT=ℓₘᵢₙ(IT)) where {IT<:IndexType}
    if IT === Int
        ℓ*(ℓ+1) - ℓₘᵢₙ^2 + m + 1
    else
        # ℓ(ℓ+1) - ℓₘᵢₙ² + m + 1 on the doubled indices a = 2ℓ, b = 2m and c = 2ℓₘᵢₙ, all odd
        # `Int`s.  Four times the index is a(a+2) - c² + 2b + 4; for odd a, b and c this is a
        # multiple of 4, so the division is exact.
        a, b, c = 2ℓ, 2m, 2ℓₘᵢₙ
        (a*(a+2) - c^2 + 2b + 4) ÷ 4
    end
end

"""
    Yrange(ℓₘₐₓ)
    Yrange(ℓₘᵢₙ, ℓₘₐₓ)

Vector of the ``(ℓ, m)`` pairs in the canonical ordering, so that `Yrange(ℓₘᵢₙ, ℓₘₐₓ)[i]` is
the pair stored at index `i`.

As for [`Ysize`](@ref), the indices may be integers of type `Int`, or half-odd-integers
passed as `Rational{Int}`s with denominator 2 or as [`HalfOddInteger`](@ref)s, and both
indices in one call must be of the same kind.  For half-odd-integers the pairs are of
`HalfOddInteger`s — whatever the type of the arguments — as is every other index value the
package hands back; use `Rational(x)` to convert an element to the more familiar `Rational`.
`ℓₘᵢₙ` defaults to the smallest ``ℓ`` of the given kind, which is 0 for integers and 1/2 for
half-odd-integers.

See also [`Ysize`](@ref) and [`Yindex`](@ref).
"""
function Yrange end

# The length of the ordering is `Ysize`, so the vector is allocated once at that length and
# then filled, which also gives `Yrange` exactly the validation and the messages of `Ysize`.
@index_methods function Yrange(ℓₘᵢₙ::IT, ℓₘₐₓ::IT) where {IT<:IndexType}
    mode_pairs!(Vector{Tuple{IT, IT}}(undef, Ysize(ℓₘᵢₙ, ℓₘₐₓ)), ℓₘᵢₙ, ℓₘₐₓ)
end
# `UnitRange{HalfOddInteger}` iterates, so this is the same loop for either kind of index.
function mode_pairs!(ordering, ℓₘᵢₙ, ℓₘₐₓ)
    i = 0
    for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ, m ∈ -ℓ:ℓ
        ordering[i += 1] = (ℓ, m)
    end
    ordering
end
@index_methods Yrange(ℓₘₐₓ::IT) where {IT<:IndexType} = Yrange(ℓₘᵢₙ(IT), ℓₘₐₓ)
