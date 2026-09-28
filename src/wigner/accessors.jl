### Accessors
#
# These names are shared by every container and calculator in the package, so their
# docstrings live here — ahead of the files that add methods to them — rather than being
# attached to whichever method happens to be defined first.  Each has an ASCII alias for use
# where the subscripted Unicode names are inconvenient, passed as the ASCII aliases of the
# keyword arguments are: `ell_min` for `ℓₘᵢₙ`, `mp_max` for `m′ₘₐₓ`, and so on.

"""
    ℓ(x)
    ell(x)

The value of ``ℓ`` that `x` currently holds: the degree of a block ([`WignerMatrix`](@ref),
[`DegreeBlock`](@ref), [`HWedge`](@ref), …), or of the block most recently computed by a
calculator.  For a calculator on which [`recurrence!`](@ref) has not been called since it
was built, or since its data were last replaced (by [`set_R!`](@ref), [`set_β!`](@ref),
[`set_θ!`](@ref) or `fill!`), this is `ℓₘᵢₙ(x) - 1`.
"""
function ℓ end

"""
    ℓₘᵢₙ(x)
    ell_min(x)

The smallest ``ℓ`` that `x` holds or can hold.

For an index, a block and a calculator this is the smallest degree of the index type: `0`
for integer indices and `1//2` for half-integer ones.  It is also defined on the index type
itself, as `ℓₘᵢₙ(Int)` and `ℓₘᵢₙ(HalfOddInteger)`; a `Rational` is something users pass at
the entry points rather than an index type, so `ℓₘᵢₙ(Rational{Int})` is not defined.  For a
[`ModeWeights`](@ref), a [`HarmonicValues`](@ref), a [`WignerSeries`](@ref) and a transform
it is the smallest ``ℓ`` stored, which may be larger, such as the ``|s|`` of the weights
that an analysis returns.
"""
function ℓₘᵢₙ end

"""
    ℓₘₐₓ(x)
    ell_max(x)

The largest ``ℓ`` that `x` can hold — the value its storage was sized for.
"""
function ℓₘₐₓ end

"""
    m′ₘₐₓ(x)
    mp_max(x)

The largest ``m'`` in the block `x` holds or returns.  For a calculator this is the keyword
argument of the same name, which restricts the rows that are computed.
"""
function m′ₘₐₓ end

"""
    m′ₘᵢₙ(x)
    mp_min(x)

The smallest ``m'`` in the block `x` holds or returns.  Note that both ``±ℓₘᵢₙ`` must lie in
`m′ₘᵢₙ(x):m′ₘₐₓ(x)`, because the recurrence seeds the half-integer ladder from the pair of
rows ``m' = ±1/2``; for integer indices that reduces to the familiar `m′ₘᵢₙ ≤ 0 ≤ m′ₘₐₓ`.

The blocks with a single ``m`` axis, [`DegreeBlock`](@ref), [`SpinMatrix`](@ref) and their
batches, have no ``m'``, and this and [`m′ₘₐₓ`](@ref) refuse them with an explanation.
"""
function m′ₘᵢₙ end

"""
    mₘₐₓ(x)
    m_max(x)

The largest ``m`` in the block `x` holds or returns.
"""
function mₘₐₓ end

"""
    mₘᵢₙ(x)
    m_min(x)

The smallest ``m`` in the block `x` holds or returns.  For a [`WignerMatrix`](@ref), a
[`WignerMatrixBatch`](@ref) or a calculator, the range must bracket ``±ℓₘᵢₙ``, as for
[`m′ₘᵢₙ`](@ref), because the recurrence computes those blocks.  The blocks that no
recurrence fills, [`DegreeBlock`](@ref), [`SpinMatrix`](@ref) and their batches, require
only `-ℓ ≤ mₘᵢₙ(x) ≤ mₘₐₓ(x) ≤ ℓ`.
"""
function mₘᵢₙ end

"""
    sₘₐₓ(x)
    s_max(x)

The largest ``s`` in the block `x` holds.  Unlike [`m′ₘₐₓ`](@ref), this is under no
obligation to bracket zero or to stay within ``±ℓ``: a block of spin-weighted harmonics may
hold any consecutive run of spin weights, and those with ``|s| > ℓ`` are zero.
"""
function sₘₐₓ end

"""
    sₘᵢₙ(x)
    s_min(x)

The smallest ``s`` in the block `x` holds; see [`sₘₐₓ`](@ref).
"""
function sₘᵢₙ end

"""
    spins(c)

The spin weights an [`sYlmCalculator`](@ref) serves, as a range.  A calculator built for a
single spin weight reports the one-element range containing it, so that this accessor
answers in the same currency whichever way the calculator was built; [`spin`](@ref) gives
the value itself, and has no method for a calculator built for several.

For a [`SpinMatrix`](@ref) or a [`SpinMatrixBatch`](@ref), which is what such a calculator
yields when it serves several spin weights, this is the range of spin weights the block
holds, `sₘᵢₙ(b):sₘₐₓ(b)`.  For a [`HarmonicValues`](@ref) it is the range of spin weights
the values were computed for, again a one-element range when there is only one.
"""
function spins end

"""
    Nᵣ(x)
    Nr(x)

The number of rotors `x` handles at once.  A calculator built from a vector of rotor data
evaluates the whole batch in one pass, and its blocks have a leading rotor index, `𝔇ˡ[iᵣ,
m′, m]`, even when the vector has only one element; [`isbatched`](@ref) says which kind a
calculator is, since `Nᵣ == 1` does not.
"""
function Nᵣ end

"""
    ishalfinteger(w)

Whether the natural indices of `w` are half-odd-integers (`true`) or integers (`false`).
This is a property of the index type — half-odd-integer indices are represented as
[`HalfOddInteger`](@ref)s — so it is known at compile time.
"""
function ishalfinteger end

"""
    isbatched(c)

Whether the calculator `c` was built for a batch of rotors (`true`, when it was built from a
vector of rotor data, of any length) or for a single one (`false`).  Like
[`ishalfinteger`](@ref), this is a property of the type rather than of the data, so it is
known at compile time — which is what lets the block that [`recurrence!`](@ref) returns have
a single, inferrable type rather than a union of the batched and unbatched ones.  For a
[`HarmonicValues`](@ref) it says whether the storage has a rotor axis, which again it has
exactly when the values were computed for a vector of rotors.  For a block it says whether
the block has a leading rotor axis: it is `true` for a [`WignerMatrixBatch`](@ref), a
[`DegreeBlockBatch`](@ref) and a [`SpinMatrixBatch`](@ref), which is what a batched
calculator returns, and `false` for the other blocks.
"""
function isbatched end


# The ASCII aliases of the accessors, which are the same functions under other names.
const ell = ℓ
const ell_min = ℓₘᵢₙ
const ell_max = ℓₘₐₓ
const mp_max = m′ₘₐₓ
const mp_min = m′ₘᵢₙ
const m_max = mₘₐₓ
const m_min = mₘᵢₙ
const s_max = sₘₐₓ
const s_min = sₘᵢₙ
const Nr = Nᵣ
