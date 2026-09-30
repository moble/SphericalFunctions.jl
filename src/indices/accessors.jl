### Accessors
#
# These names are shared by every container and calculator in the package, so their
# docstrings live here — ahead of the files that add methods to them — rather than being
# attached to whichever method happens to be defined first.  Those whose names are not ASCII
# have ASCII aliases, defined in `aliases.jl` and spelled as the ASCII aliases of the
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

For a block or a calculator, this is the smallest degree of the index type: `0` for integer
indices and `1//2` for half-integer ones.  For a [`ModeWeights`](@ref), a
[`HarmonicValues`](@ref), a [`WignerSeries`](@ref), or a transform, it is the smallest
``ℓ`` stored, which may be larger, such as the ``|s|`` of the weights that an analysis
returns.
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
    spin(w)

The spin weight of a [`ModeWeights`](@ref) vector, of an [`SSHT`](@ref) transform, or of an
[`sYlmCalculator`](@ref) built for a single one.  A calculator built for a range of spin
weights has no single value to report, so it has no method here; ask it for [`spins`](@ref
SphericalFunctions.spins) instead, which answers for either kind.

A function of spin weight ``s`` has ``R_z f = s f``, and is expanded in the harmonics
``{}_{s}Y_{ℓ,m}`` with ``ℓ ≥ |s|``.  The spin weight is kept alongside the numbers because
nothing about the numbers themselves reveals it.

```jldoctest
julia> using SphericalFunctions

julia> spin(ModeWeights(zeros(ComplexF64, 21), -2))
-2

julia> spin(SSHT(1, 4))
1
```

See also [`modes`](@ref), [`ModeWeights`](@ref), [`SSHT`](@ref), and
[`spins`](@ref SphericalFunctions.spins).
"""
function spin end

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
    floattype(x)
    floattype(T)

The floating-point type in which `x`, or anything of type `T`, is computed.  For a
calculator or a transform, this is the type it works in, and its results are that type, or
`Complex` of it.  For rotor data — a `Rotor` or other `Quaternion`, an angle `β::Real`, a
phase `e^{iβ}::Complex`, or an `AbstractVector` of any one of those — it is the type in
which a calculator built from that data works: the `float` of the data's component type, so
that a `Rotor{Float32}` gives a `Float32` calculator, a `Rotor{Double64}` a `Double64` one,
and the angle `3` a `Float64` one.  As with `eltype`, the method for a value is that for its
type.

The type of the rotor data is the *only* thing that decides the type in which a calculator
works.  There is no way to set it independently of the data: to compute in some other type,
convert the data — which is also the honest way to say it, since the type of the data is
the claim being made about the points.  The methods that replace a calculator's data —
[`set_R!`](@ref), [`set_β!`](@ref), [`set_θ!`](@ref) — require the new data to agree with
what this reports.  A `Quaternion` counts as a rotor here; we always divide out by the norm,
and accept `Quaternion` for compatibility with various automatic-differentiation packages.

The rotor data must therefore commit to a type to go on: data whose component type is
abstract, such as a `Complex{Real}` or a `Vector{Rotor{Real}}`, or a vector whose element
type is abstract or a `Union`, is refused with an `ArgumentError` rather than guessed at, as
is a `QuatVec`, which does not denote a rotation (see [`not_a_rotor`](@ref)).
"""
function floattype end

"""
    nmodes(𝒯)

Number of mode weights the transform `𝒯` works with, `Ysize(abs(s), ℓₘₐₓ)`.
"""
function nmodes end

"""
    npixels(𝒯)

Number of sample points (function values) the transform `𝒯` works with.
"""
function npixels end

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

