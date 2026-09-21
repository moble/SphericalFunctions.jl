### Accessors
#
# These names are shared by every container and calculator in the package, so their
# docstrings live here — ahead of the files that add methods to them — rather than being
# attached to whichever method happens to be defined first.  Each has an ASCII alias for
# use where the subscripted Unicode names are inconvenient.

"""
    ℓ(x)
    ell(x)

The value of ``ℓ`` that `x` currently holds: the order of a block ([`WignerMatrix`](@ref),
[`DegreeBlock`](@ref), [`HWedge`](@ref), …), or of the block most recently computed by a
calculator.  For a calculator on which [`recurrence!`](@ref) has not yet been called, this
is `ℓₘᵢₙ(x) - 1`.
"""
function ℓ end

"""
    ℓₘᵢₙ(x)
    ellmin(x)

The smallest ``ℓ`` that `x` can hold: `0` for an integer index type and `1//2` for a
half-integer one.  Also defined on the index type itself, as `ℓₘᵢₙ(Int)` and
`ℓₘᵢₙ(HalfOddInteger)`; a `Rational` is a spelling that users write at the entry points rather
than an index type, so `ℓₘᵢₙ(Rational{Int})` is not defined.
"""
function ℓₘᵢₙ end

"""
    ℓₘₐₓ(x)
    ellmax(x)

The largest ``ℓ`` that `x` can hold — the value its storage was sized for.
"""
function ℓₘₐₓ end

"""
    m′ₘₐₓ(x)
    mpmax(x)

The largest ``m'`` in the block `x` holds or returns.  For a calculator this is the keyword
argument of the same name, which restricts the rows that are computed.
"""
function m′ₘₐₓ end

"""
    m′ₘᵢₙ(x)
    mpmin(x)

The smallest ``m'`` in the block `x` holds or returns.  Note that both ``±ℓₘᵢₙ`` must lie in
`m′ₘᵢₙ(x):m′ₘₐₓ(x)`, because the recurrence seeds the half-integer ladder from the pair of
rows ``m' = ±1/2``; for integer indices that reduces to the familiar `m′ₘᵢₙ ≤ 0 ≤ m′ₘₐₓ`.
"""
function m′ₘᵢₙ end

"""
    mₘₐₓ(x)
    mmax(x)

The largest ``m`` in the block `x` holds or returns.
"""
function mₘₐₓ end

"""
    mₘᵢₙ(x)
    mmin(x)

The smallest ``m`` in the block `x` holds or returns.  As for [`m′ₘᵢₙ`](@ref), the range
must bracket ``±ℓₘᵢₙ``.
"""
function mₘᵢₙ end

"""
    sₘₐₓ(x)
    smax(x)

The largest ``s`` in the block `x` holds.  Unlike [`m′ₘₐₓ`](@ref), this is under no
obligation to bracket zero or to stay within ``±ℓ``: a block of spin-weighted harmonics may
hold any consecutive run of spin weights, and those with ``|s| > ℓ`` are zero.
"""
function sₘₐₓ end

"""
    sₘᵢₙ(x)
    smin(x)

The smallest ``s`` in the block `x` holds; see [`sₘₐₓ`](@ref).
"""
function sₘᵢₙ end

"""
    spins(c)

The spin weights an [`sYlmCalculator`](@ref) serves, as a range.  A calculator built for a
single spin weight reports the one-element range containing it, so that this accessor answers
in the same currency whichever way the calculator was built; [`spin`](@ref) gives the value
itself, and has no method for a calculator built for several.
"""
function spins end

"""
    Nᵣ(x)
    Nr(x)

The number of rotors `x` handles at once.  A calculator with `Nᵣ > 1` evaluates a whole
batch in one pass, and its blocks have a leading rotor index: `calc[ℓ][iᵣ, m′, m]`.
"""
function Nᵣ end

"""
    ishalfinteger(w)

Whether the natural indices of `w` are half-odd-integers (`true`) or integers (`false`).
This is a property of the index type, so it is known at compile time.

(Before version 3 this was called `isrational`, because half-integer indices were
represented as `Rational`s with denominator 2; they are now [`HalfOddInteger`](@ref)s.)
"""
function ishalfinteger end

"""
    isbatched(c)

Whether the calculator `c` was built for a batch of rotors (`true`, when `Nᵣ > 1`) or for a
single one (`false`).  Like [`ishalfinteger`](@ref), this is a property of the type rather
than of the data, so it is known at compile time — which is what lets `calc[ℓ]` have a
single, inferrable return type rather than a union of the batched and unbatched ones.
"""
function isbatched end
