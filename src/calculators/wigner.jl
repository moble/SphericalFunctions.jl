# The last parameter, `B`, records whether the calculator was built from a *vector* of rotor
# data — batched, with blocks indexed `[iᵣ, m′, m]` — or from a single rotor.  It is decided
# by the type of that data rather than by its length, so that a one-element vector gives
# batched blocks like any other vector, and so that `B` is known at compile time: the block
# then has one type rather than a union of the batched and unbatched ones, because the
# branch in `block` below is resolved at compile time.  Only the predicate is lifted into
# the type, not `Nᵣ` itself: `SSHTRS` sets `Nᵣ = Nθ`, which grows with ℓₘₐₓ, and
# parameterizing on the count would recompile the recurrence for every resolution.

"""
    WignerCalculator{IT, RT, NT, ST, B, FT, L}

Calculator producing Wigner's ``𝔇`` matrices (when `NT` is `Complex{RT}`) or ``d`` matrices
(when `NT` is `RT`) for `Nᵣ` rotors at a time, one ``ℓ`` at a time.  Use the constructors
[`DCalculator`](@ref) and [`dCalculator`](@ref).
The type parameters are as follows:
- `IT` is the index type, `Int` or [`HalfOddInteger`](@ref).
- `RT` is the real type of the rotor data, and of the elements.
- `NT` is the number type of the matrix elements, `Complex{RT}` or `RT`.
- `ST` is the storage type of the ``H`` wedge.
- `B` is `true` exactly when the calculator was built from a vector of rotor data, of any
  length, and is what [`isbatched`](@ref) reads.
- `FT` is the real type in which the recurrence runs, which is `RT` itself unless `RT`
  holds derivatives, as a dual number does.
- `L` is `Nothing`, unless the calculator lifts the blocks of a calculator of the values of
  its rotors into blocks that hold derivatives, when it is the type of the data for that.

Internally this holds an [`HCalculator`](@ref), which does the actual recurrence, and the
tables of the powers of the rotors' phases, in the same form as the calculators of the
harmonics, plus a buffer into which the requested block of the matrix is written for the
current ``ℓ``; that block is what [`recurrence!`](@ref) returns, as a
[`WignerMatrix`](@ref) indexed naturally by `[m′, m]`, or a [`WignerMatrixBatch`](@ref)
indexed by `[iᵣ, m′, m]` when the calculator was built from a vector of rotor data.

A calculator of ``𝔇`` whose rotors hold derivatives, such as dual numbers, runs the
recurrence on the values of those rotors, in `FT`, and gives each block its derivatives
from the angular-momentum operators, as described in `src/derivatives/kernels.jl`; the
recurrence itself is never differentiated.

Because `B` is a type parameter, the type of the block is known at compile time, and a loop
over the blocks is inferrable.
"""
struct WignerCalculator{
    IT, RT<:Real, NT<:Union{RT, Complex{RT}}, ST, B, FT<:Real, L
} <: AbstractCalculator{IT}
    engine::SphericalFunctionsEngine{IT, FT, ST}  # the recurrence and the power tables
    Wˡ::Array{NT, 3}  # [iᵣ, m′, m] stored rows of the block for the current ℓ, from the first
    rotors::Vector{Quaternion{RT}}  # the rotors themselves; empty when NT is real
    m′ₘₐₓ::IT
    m′ₘᵢₙ::IT
    mₘₐₓ::IT
    mₘᵢₙ::IT
    m′ₘₐₓˢ::IT  # the limits of the rows and columns stored in Wˡ, which extend one beyond
    m′ₘᵢₙˢ::IT  # the limits of the block along one axis when the block is restricted in
    mₘₐₓˢ::IT   # both m′ and m; see `stored_limits`
    mₘᵢₙˢ::IT
    ℓ::Base.RefValue{IT}  # ℓ of the block currently in Wˡ; ℓₘᵢₙ-1 if none
    lift::L  # `nothing`, or the calculator of the rotors' values and their generators
    # `materialize!` writes `Wˡ` and reads the power tables under `@inbounds`, for every
    # rotor of the engine and every (m′, m) within the limits, so the buffers must be large
    # enough for those; as for `HCalculator`, this checks them once, as they are brought
    # together.
    function WignerCalculator{IT, RT, NT, ST, B, FT, L}(
        engine, Wˡ, rotors, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, m′ₘₐₓˢ, m′ₘᵢₙˢ, mₘₐₓˢ, mₘᵢₙˢ, ℓ, lift
    ) where {IT, RT<:Real, NT<:Union{RT, Complex{RT}}, ST, B, FT<:Real, L}
        let n = Nᵣ(engine), K = NT <: Complex ? 2ℓₘₐₓ(engine) + 1 : 0,
                Z₊ = engine.Z₊, Z₋ = engine.Z₋
            if !(
                size(Wˡ, 1) ≥ n
                && size(Wˡ, 2) ≥ Int(m′ₘₐₓˢ - m′ₘᵢₙˢ) + 1 && size(Wˡ, 3) ≥ Int(mₘₐₓˢ - mₘᵢₙˢ) + 1
                && size(Z₊, 1) ≥ n && size(Z₊, 2) ≥ K && size(Z₋, 1) ≥ n && size(Z₋, 2) ≥ K
                && length(rotors) == (NT <: Complex ? n : 0)
                && m′ₘₐₓˢ ≥ m′ₘₐₓ && m′ₘᵢₙˢ ≤ m′ₘᵢₙ && mₘₐₓˢ ≥ mₘₐₓ && mₘᵢₙˢ ≤ mₘᵢₙ
            )
                throw(DimensionMismatch(
                    "The buffers of a WignerCalculator for Nᵣ=$n rotors, m′ ∈ $m′ₘᵢₙ:$m′ₘₐₓ "
                    * "(stored as $m′ₘᵢₙˢ:$m′ₘₐₓˢ) and m ∈ $mₘᵢₙ:$mₘₐₓ (stored as "
                    * "$mₘᵢₙˢ:$mₘₐₓˢ) are too small: the "
                    * "block has size $(size(Wˡ)), the power tables $(size(Z₊)) and "
                    * "$(size(Z₋)), which need a row for each rotor and at least $K columns, "
                    * "and there are $(length(rotors)) rotors."
                ))
            end
        end
        new{IT, RT, NT, ST, B, FT, L}(
            engine, Wˡ, rotors, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, m′ₘₐₓˢ, m′ₘᵢₙˢ, mₘₐₓˢ, mₘᵢₙˢ,
            ℓ, lift
        )
    end
end

# The limits of the rows and columns that a calculator of 𝔇 stores for the block it
# returns.  The derivatives of the block, which the rules for automatic differentiation give
# in terms of its values (see `src/derivatives/kernels.jl`), couple each element either to
# its neighbors in the same column (from the left) or to those in the same row (from the
# right).  A block whose rows are all of -ℓ:ℓ is differentiated from the left, and one whose
# columns are, from the right, so that every neighbor needed is in the block already.  Only
# a block restricted in both needs values beyond its limits, and for that one, one row or
# column more is stored on each side, within ±ℓₘₐₓ, along whichever axis is the wider, so
# that the wedge, whose width is that of the narrower (see `allocate_W`), is widened only
# when the two are equally wide.  The widened limits still bracket ±ℓₘᵢₙ, as the originals
# do.  A calculator of `d` stores exactly its block.
function stored_limits(ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT) where {IT}
    full(hi, lo) = hi == ℓₘₐₓ && lo == -ℓₘₐₓ
    widened(hi, lo) = (min(hi + 1, ℓₘₐₓ), max(lo - 1, -ℓₘₐₓ))
    if full(m′ₘₐₓ, m′ₘᵢₙ) || full(mₘₐₓ, mₘᵢₙ)
        (m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    elseif max(m′ₘₐₓ, -m′ₘᵢₙ) ≥ max(mₘₐₓ, -mₘᵢₙ)
        (widened(m′ₘₐₓ, m′ₘᵢₙ)..., mₘₐₓ, mₘᵢₙ)
    else
        (m′ₘₐₓ, m′ₘᵢₙ, widened(mₘₐₓ, mₘᵢₙ)...)
    end
end

# Whether the derivatives of a calculator's blocks are taken from the left, coupling
# neighbors in a column, or from the right, coupling neighbors in a row, as described under
# `stored_limits`.
function derivatives_from_left(c::WignerCalculator)
    if c.m′ₘₐₓ == ℓₘₐₓ(c) && c.m′ₘᵢₙ == -ℓₘₐₓ(c)
        true
    elseif c.mₘₐₓ == ℓₘₐₓ(c) && c.mₘᵢₙ == -ℓₘₐₓ(c)
        false
    else
        c.m′ₘₐₓˢ != c.m′ₘₐₓ || c.m′ₘᵢₙˢ != c.m′ₘᵢₙ
    end
end

# Allocate the buffers without touching them.  PRIVATE: see the note on `allocate_H`.  The
# rotor data here is the engine's (see `copy_rotor_data!`) *plus* the rotors themselves, so
# a caller that copies rather than sets must copy both, as `copy_rotor_state!` does.  A
# calculator of 𝔇 whose real type holds derivatives is given a calculator of the rotors'
# values, whose engine it holds as its own, and whose blocks it lifts; see `allocate_lift`.
function allocate_W(
    ::Type{IT}, ::Type{RT}, ::Type{NT}, ℓₘₐₓ::IT,
    m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT, Nᵣ::Int, ::Val{B}
) where {IT<:IntegerHalf, RT<:Real, NT<:Union{RT, Complex{RT}}, B}
    validate_degree(ℓₘₐₓ)
    validate_axis(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, "m′")
    validate_axis(ℓₘₐₓ, mₘₐₓ, mₘᵢₙ, "m")
    m′ₘₐₓˢ, m′ₘᵢₙˢ, mₘₐₓˢ, mₘᵢₙˢ = if NT <: Complex
        stored_limits(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    else
        (m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    end
    Wˡ = Array{NT, 3}(undef, Nᵣ, Int(m′ₘₐₓˢ - m′ₘᵢₙˢ) + 1, Int(mₘₐₓˢ - mₘᵢₙˢ) + 1)
    rotors = Vector{Quaternion{RT}}(undef, NT <: Complex ? Nᵣ : 0)
    # The field is a `RefValue{IT}`, so the type is given explicitly, as in `allocate_Y`.
    ℓ = Ref{IT}(lowest_index(IT) - 1)
    if NT <: Complex && value_type(RT) !== RT
        let inner = allocate_W(
            IT, value_type(RT), Complex{value_type(RT)}, ℓₘₐₓ,
            m′ₘₐₓˢ, m′ₘᵢₙˢ, mₘₐₓˢ, mₘᵢₙˢ, Nᵣ, Val(B)
        )
            lift = allocate_lift(RT, inner, Nᵣ)
            WignerCalculator{IT, RT, NT, typeof(parent(inner.engine.H.Hˡ)), B, recurrence_type(RT), typeof(lift)}(
                inner.engine, Wˡ, rotors,
                m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, m′ₘₐₓˢ, m′ₘᵢₙˢ, mₘₐₓˢ, mₘᵢₙˢ, ℓ, lift
            )
        end
    else
        # The wedge holds the rows |m′| ≤ W of H, for every m, and `materialize!` reads an
        # element whose |m′| exceeds W from its transpose, whose |m| does not.  So W need
        # only be the narrower of the widest m′ and the widest m of the stored block, and
        # restricting either the rows or the columns of the block to a narrow band narrows
        # the recurrence with it.  Both limits bracket ±ℓₘᵢₙ, so W does too, as the
        # recurrence requires.
        W = min(max(m′ₘₐₓˢ, -m′ₘᵢₙˢ), max(mₘₐₓˢ, -mₘᵢₙˢ))
        engine = allocate_engine(IT, RT, ℓₘₐₓ, W, Nᵣ, NT <: Complex)
        WignerCalculator{IT, RT, NT, typeof(parent(engine.H.Hˡ)), B, RT, Nothing}(
            engine, Wˡ, rotors,
            m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, m′ₘₐₓˢ, m′ₘᵢₙˢ, mₘₐₓˢ, mₘᵢₙˢ, ℓ, nothing
        )
    end
end

"""
    DCalculator(R, ℓₘₐₓ; m′ₘₐₓ=ℓₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ, mₘₐₓ=ℓₘₐₓ, mₘᵢₙ=-mₘₐₓ)
    DCalculator(α, β, γ, ℓₘₐₓ; m′ₘₐₓ=ℓₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ, mₘₐₓ=ℓₘₐₓ, mₘᵢₙ=-mₘₐₓ)

Calculator for Wigner's ``𝔇^{(ℓ)}_{m′,m}(R)`` matrices, for ``ℓ ≤ ℓₘₐₓ``, with elements of
type `Complex{RT}`.  The keyword arguments restrict the block of each matrix that is
computed and returned; they may also be spelled `mp_max`, `mp_min`, `m_max` and `m_min`.
Restricting either the rows ``m′`` or the columns ``m`` to a narrow band also narrows the
recurrence itself, which at each ``ℓ`` costs in proportion to the narrower of the two
ranges.

The rotor is given at construction, so that the calculator is usable the moment it exists,
and so that the element type is the rotor's own: a `Rotor{Float32}` gives a `Float32`
calculator, and `floattype(calc)` reports it.  There is no argument to override that — to
compute in another type, convert the rotor, which is the honest way to say that those are
the values to be treated as exact.  Give `R` as an `AbstractVector` of rotors to get a
calculator that handles all `Nᵣ = length(R)` of them at once — which is substantially faster
per rotor, for the reason described under [Reusing the workspace](@ref
interface_wigner_matrices).  Later rotors are supplied with [`set_R!`](@ref).  A single
rotor may instead be given by its Euler angles: `DCalculator(α, β, γ, ℓₘₐₓ)` is
`DCalculator(from_euler_angles(α, β, γ), ℓₘₐₓ)`.

The calculator is iterable, yielding one ``ℓ`` at a time:

```julia
calc = DCalculator(R, ℓₘₐₓ)
for (ℓ, 𝔇ˡ) ∈ calc
    # 𝔇ˡ[m′, m] with m′, m ∈ -ℓ:ℓ
end
```

For a calculator built from a vector of rotors — of any length, even one — each block is
indexed as `[iᵣ, m′, m]` instead.  The block is a view into the calculator's storage and is
overwritten by the next step, so `copy` it if it must survive (the copy keeps the natural
indices), or `collect` it to get an ordinary 1-based `Array`: a `Matrix`, or a
three-dimensional `Array` for a block with a rotor index.  `collect(calc)` copies every
block, so it is safe.  To sweep part of the range, or to take the values of ``ℓ`` in some
other order, loop over [`recurrence!`](@ref) yourself.

A calculator is a mutable workspace, which every step and every setter overwrites, so one
calculator must not be used by two tasks at once; `similar(calc)` gives each task a
calculator of its own, holding the same rotor data.

The convention is ``𝔇^{(ℓ)}_{m′,m}(𝐑_{α,β,γ}) = e^{-im′α}\\, d^{(ℓ)}_{m′,m}(β)\\,
e^{-imγ}``; see the "Conventions" section of the documentation.

# Half-integer indices

`DCalculator(R, 7//2)` — a `Rational` `ℓₘₐₓ` with denominator 2, or a
[`HalfOddInteger`](@ref) — gives a calculator for half-integer ``ℓ, m′, m``.  All four
keyword limits must then be half-integers too, in either spelling; `recurrence!` accepts
only half-integer `ℓ` and returns a [`WignerMatrix`](@ref) (or a [`WignerMatrixBatch`](@ref)
when built from a vector) whose indices are half-odd-integers; it is indexed the same way.
The double cover is respected exactly: ``𝔇(-R) = -𝔇(R)``.

See also [`D`](@ref) for a simpler interface when the matrices for only one rotor are
needed, [`dCalculator`](@ref) for the real ``d`` matrices, and [`recurrence!`](@ref).
"""
const DCalculator{IT, RT, ST, B} = WignerCalculator{IT, RT, Complex{RT}, ST, B} where {IT, RT<:Real, ST, B}

"""
    dCalculator(β, ℓₘₐₓ; m′ₘₐₓ=ℓₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ, mₘₐₓ=ℓₘₐₓ, mₘᵢₙ=-mₘₐₓ)

Calculator for Wigner's real ``d^{(ℓ)}_{m′,m}(β)`` matrices, for ``ℓ ≤ ℓₘₐₓ``, with elements
of the angle's own floating-point type.  The first argument may be the angle ``β``, the
phase ``e^{iβ}``, or a `Rotor` or other `Quaternion`, or an `AbstractVector` of `Nᵣ` of any
one of those; later values are supplied with [`set_β!`](@ref).  Otherwise this behaves
exactly like [`DCalculator`](@ref) — including iteration, the ASCII spellings of the
keywords, and half-integer ``ℓ`` for a `Rational` or [`HalfOddInteger`](@ref) `ℓₘₐₓ`.

A `Rotor` contributes the ``β ∈ [0, π]`` of its canonical Euler decomposition.  If the rotor
was built from a ``β`` outside that range, that ``β`` is folded into its ``α`` and ``γ``,
which a ``d`` matrix does not see, so `d` of such a rotor equals `d` of its angle only for
``β ∈ [0, π]``; this is consistent with [`D`](@ref) of the same rotor.

Half-integer ``d`` has period ``4π`` in ``β``, so an angle or a `Rotor` determines it
unambiguously, while a bare phase ``e^{iβ}`` determines it only up to the double-cover sign
``(-1)^{2ℓ}``; the branch ``β ∈ (-π, π]`` is used in that case.

See also [`d`](@ref).
"""
const dCalculator{IT, RT, ST, B} = WignerCalculator{IT, RT, RT, ST, B, RT, Nothing} where {IT, RT<:Real, ST, B}
# As for `sλlmCalculator`, a calculator of `d` never lifts the blocks of another, so the last
# two parameters above are fixed.

# The keyword limits are normalized by `@index_methods` against the kind of `ℓₘₐₓ`, so that
# either spelling of each, and either spelling of a half-integer, reaches the work method as
# an `Int` or a `HalfOddInteger`.  `R` is untyped, but `ℓₘₐₓ` is an index, so a call with
# ℓₘₐₓ first — `DCalculator(ℓₘₐₓ, Float64)` — is an immediate `MethodError` at the call site
# rather than something that dispatches with the element type in the rotor's place.  The
# element type is `floattype(R)`, which depends on the type of `R` alone, so that the
# compiler settles the concrete return type, including the `B` parameter that the block's
# type depends on.
#
# No method of `DCalculator` or `dCalculator` takes the element type as an argument, since
# that would be a public way to override it, and there is none by design — the type of the
# rotor data is the only thing that decides it.  For the same reason there is no
# constructor `WignerCalculator{IT, RT, NT}(R, ℓₘₐₓ)`: `DCalculator{Int, Float64}` names
# that very type, and is what `show` prints.  `test/wigner/iteration.jl` asserts exactly
# that, with `@test_throws MethodError`.
@index_methods function DCalculator(
    R, ℓₘₐₓ::IT;
    mp_max::IndexType=ℓₘₐₓ, m′ₘₐₓ::IndexType=mp_max,
    mp_min::IndexType=-m′ₘₐₓ, m′ₘᵢₙ::IndexType=mp_min,
    m_max::IndexType=ℓₘₐₓ, mₘₐₓ::IndexType=m_max,
    m_min::IndexType=-mₘₐₓ, mₘᵢₙ::IndexType=m_min
) where {IT<:IndexType}
    RT = floattype(R)
    set_rotors!(
        allocate_W(
            IT, RT, Complex{RT}, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, nrotors(R), batched_data(R)
        ),
        R
    )
end
@index_methods function DCalculator(
    α::Real, β::Real, γ::Real, ℓₘₐₓ::IndexType;
    mp_max::IndexType=ℓₘₐₓ, m′ₘₐₓ::IndexType=mp_max,
    mp_min::IndexType=-m′ₘₐₓ, m′ₘᵢₙ::IndexType=mp_min,
    m_max::IndexType=ℓₘₐₓ, mₘₐₓ::IndexType=m_max,
    m_min::IndexType=-mₘₐₓ, mₘᵢₙ::IndexType=m_min
)
    DCalculator(Quaternionic.from_euler_angles(α, β, γ), ℓₘₐₓ; m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
end
@index_methods function dCalculator(
    β, ℓₘₐₓ::IT;
    mp_max::IndexType=ℓₘₐₓ, m′ₘₐₓ::IndexType=mp_max,
    mp_min::IndexType=-m′ₘₐₓ, m′ₘᵢₙ::IndexType=mp_min,
    m_max::IndexType=ℓₘₐₓ, mₘₐₓ::IndexType=m_max,
    m_min::IndexType=-mₘₐₓ, mₘᵢₙ::IndexType=m_min
) where {IT<:IndexType}
    RT = floattype(β)
    set_rotors!(
        allocate_W(IT, RT, RT, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, nrotors(β), batched_data(β)),
        β
    )
end

"""
    similar(calc)
    similar(calc, R)

A second calculator of exactly the type of `calc` — the same limits, number of rotors and
element type — with storage of its own and nothing computed.  With one argument it holds a
copy of the rotor data of `calc`, so that stepping it gives the values `calc` would give;
with two it holds the rotor data `R` instead, which must describe the same number of rotors,
in the same floating-point type.  This applies to every calculator: [`DCalculator`](@ref),
[`dCalculator`](@ref), [`HCalculator`](@ref), [`sYlmCalculator`](@ref) and
[`sλlmCalculator`](@ref).

Unlike `similar` of an array, this does copy the rotor data, because a calculator without
rotor data is not representable, and the rotor data cannot be recovered from what a
calculator stores.  The second argument is rotor data rather than an element type; to
compute in another type, build a calculator from rotor data of that type.

A calculator is a mutable workspace, which every step and every setter overwrites, so one
calculator must not be used by two tasks at once.  `similar` is the way to give each task a
calculator of its own; the two share nothing that either one changes.
"""
function Base.similar(
    c::WignerCalculator{IT, RT, NT, ST, B, FT, L}
) where {IT, RT, NT, ST, B, FT<:Real, L}
    # A second workspace holding the same rotor data, with nothing computed.  The assertion
    # is what keeps this inferrable: `Nᵣ(c)` is a field lookup, so the constructor cannot
    # know `B`, but the copy necessarily has the same parameters as the original.  The data
    # is copied buffer by buffer rather than derived again; see the comment on
    # `copy_rotor_data!` for why it cannot be derived again from the `H` wedge.
    c′ = allocate_W(
        IT, RT, NT, ℓₘₐₓ(c), c.m′ₘₐₓ, c.m′ₘᵢₙ, c.mₘₐₓ, c.mₘᵢₙ, Nᵣ(c), Val(B)
    )::WignerCalculator{IT, RT, NT, ST, B, FT, L}
    copy_rotor_state!(c′, c)
end
function Base.similar(
    c::WignerCalculator{IT, RT, NT, ST, B, FT, L}, R
) where {IT, RT, NT, ST, B, FT<:Real, L}
    check_rotor_count(c, R)
    check_rotor_type(c, R)
    set_rotors!(
        allocate_W(
            IT, RT, NT, ℓₘₐₓ(c), c.m′ₘₐₓ, c.m′ₘᵢₙ, c.mₘₐₓ, c.mₘᵢₙ, Nᵣ(c), Val(B)
        )::WignerCalculator{IT, RT, NT, ST, B, FT, L},
        R
    )
end

ℓ(c::WignerCalculator) = c.ℓ[]
ℓₘₐₓ(c::WignerCalculator) = ℓₘₐₓ(c.engine)
m′ₘₐₓ(c::WignerCalculator) = c.m′ₘₐₓ
m′ₘᵢₙ(c::WignerCalculator) = c.m′ₘᵢₙ
mₘₐₓ(c::WignerCalculator) = c.mₘₐₓ
mₘᵢₙ(c::WignerCalculator) = c.mₘᵢₙ
Nᵣ(c::WignerCalculator) = Nᵣ(c.engine)
floattype(::Type{<:WignerCalculator{IT, RT}}) where {IT, RT} = RT
isbatched(::WignerCalculator{IT, RT, NT, ST, B}) where {IT, RT, NT, ST, B} = B

# A batched calculator says so, because its blocks have a rotor index that those of an
# unbatched one with the same `Nᵣ = 1` do not.
function Base.show(io::IO, c::WignerCalculator{IT, RT, NT}) where {IT, RT, NT}
    print(
        io,
        NT <: Complex ? "D" : "d", "Calculator{$IT, $RT} for ",
        "ℓₘₐₓ=$(ℓₘₐₓ(c)), m′=$(c.m′ₘᵢₙ):$(c.m′ₘₐₓ), m=$(c.mₘᵢₙ):$(c.mₘₐₓ), Nᵣ=$(Nᵣ(c))",
        isbatched(c) ? ", batched" : "",
        c.ℓ[] < ℓₘᵢₙ(c) ? " (nothing computed yet)" : ", currently at ℓ=$(c.ℓ[])"
    )
end
function Base.show(io::IO, ::MIME"text/plain", c::WignerCalculator)
    show(io, c)
end
# The constructor that built a calculator, which is what a message calls it.
container_name(::WignerCalculator{IT, RT, NT}) where {IT, RT, NT} =
    NT <: Complex ? "DCalculator" : "dCalculator"

"""
    fill!(c::WignerCalculator, v)

Fill the axis, wedge and output buffers of `c` with the value `v` and mark the current
results as invalid.  The stored rotor data — `e^{iβ}`, the half angles, the phase powers
`Z₊` and `Z₋` — is deliberately *not* touched, so `recurrence!(c, ℓ)` still has everything
it needs, exactly as for [`HCalculator`](@ref).  Useful for testing that no uninitialized
storage is ever read: everything the recurrence is responsible for writing is poisoned,
while everything `set_rotors!` is responsible for writing is left alone.
"""
function Base.fill!(c::WignerCalculator{IT, RT, NT}, v::Number) where {IT, RT, NT}
    fill!(c.engine, real(v))
    fill!(c.Wˡ, convert(NT, v))
    c.ℓ[] = lowest_index(IT) - 1
    c
end


### Rotor data

# For 𝔇 we need the full rotor: eⁱᵝ for the recurrence, and the powers of z₊ and z₋ for the
# phases.  The rotors are kept too, as quaternions, which is what the rules for automatic
# differentiation read (see `src/derivatives/kernels.jl`).  A quaternion that is not a unit
# quaternion denotes the rotation of its normalization, as it does throughout.
#
# Setting the rotors is two steps: `store_rotors!` copies them, and `set_rotor_data!`
# computes everything else from them.  The extensions for Enzyme and Mooncake declare the
# second step to have no derivatives, so that those tools differentiate only the copy, and
# follow the rotors' tangents into the calculator's copy of them, which is what their rules
# for `compute_block!` read.
function set_rotors!(
    c::WignerCalculator{IT, RT, Complex{RT}, ST, B, FT, Nothing}, R::AbstractVector{<:RotorLike}
) where {IT, RT<:Real, ST, B, FT<:Real}
    # The loops write the calculator's 1-based buffers at the input's own indices, under
    # `@inbounds`, so an offset vector would write outside them.
    Base.require_one_based_indexing(R)
    check_rotor_count(c, R)
    store_rotors!(c.rotors, R)
    set_rotor_data!(c, R)
    c
end
function set_rotor_data!(c::WignerCalculator{IT, RT, Complex{RT}}, R::AbstractVector) where {IT, RT<:Real}
    # The results are marked invalid before the first rotor is replaced.
    c.ℓ[] = lowest_index(IT) - 1
    set_rotor_data!(c.engine, R)
end
# The rotor data of a calculator of 𝔇 or of `d` is its engine's.
function copy_rotor_data!(c′::WignerCalculator, c::WignerCalculator)
    copy_rotor_data!(c′.engine, c.engine)
    c′
end
# A calculator that lifts the blocks of the calculator of its rotors' values gives that
# calculator the values, and computes the generators of its own rotors' derivatives, which
# are all that the lifting needs.
function set_rotors!(
    c::WignerCalculator{IT, RT, Complex{RT}, ST, B, FT, <:Lift}, R::AbstractVector{<:RotorLike}
) where {IT, RT<:Real, ST, B, FT<:Real}
    Base.require_one_based_indexing(R)
    check_rotor_count(c, R)
    c.ℓ[] = lowest_index(IT) - 1
    store_rotors!(c.rotors, R)
    set_rotors!(c.lift.inner, LiftedValues(c.rotors))
    set_generators!(c.lift, derivatives_from_left(c), c.rotors)
    c
end
function set_rotors!(c::WignerCalculator{IT, RT, Complex{RT}}, R::RotorLike) where {IT, RT<:Real}
    check_rotor_count(c, R)
    set_rotors!(c, @SVector [R])
end
function set_rotors!(c::WignerCalculator{IT, RT, Complex{RT}}, R) where {IT, RT<:Real}
    throw(ArgumentError(
        "A DCalculator needs rotors, given as `Rotor`s or `Quaternion`s — one, or an "
        * "AbstractVector of $(Nᵣ(c)) of them — not $(typeof(R)); use a dCalculator if only "
        * "β is available."
    ))
end

# For d we need only eⁱᵝ, and accept anything the H calculator accepts.  That calculator
# validates everything before it replaces anything, so if it refuses the data this
# calculator is left exactly as it was, and its own `ℓ` is reset only once the new data are
# in place.
function set_rotors!(c::WignerCalculator{IT, RT, RT}, R) where {IT, RT<:Real}
    set_rotors!(c.engine.H, R)
    c.ℓ[] = lowest_index(IT) - 1
    c
end


### Driver

function recurrence!(c::WignerCalculator, R, ℓ)
    check_ℓ(c.engine.H, ℓ, c)
    check_rotor_type(c, R)  # as `set_R!` and `set_β!` do, rather than silently converting
    set_rotors!(c, R)
    recurrence!(c, ℓ)
end
function recurrence!(c::WignerCalculator{IT}, ℓ) where {IT}
    let ℓ = checked_index(IT, ℓ, c, "ℓ")
        compute_block!(c, ℓ)
        current_block(c, ℓ)
    end
end

# Compute the stored rows of the block of degree ℓ into `Wˡ`: one step of the recurrence and
# the materialization of its result, or, for a calculator that lifts the blocks of the
# calculator of its rotors' values, that calculator's step and the lifting.  Every sweep of
# a calculator is made of these steps, and they are what the extensions for Enzyme and
# Mooncake attach their rules to, so that the recurrence itself is never differentiated.
function compute_block!(c::WignerCalculator{IT, RT, NT, ST, B, FT, Nothing}, ℓ::IT) where {IT, RT, NT, ST, B, FT<:Real}
    recurrence!(c.engine.H, ℓ)
    materialize!(c, ℓ)
    c
end
function compute_block!(c::WignerCalculator{IT, RT, NT, ST, B, FT, <:Lift}, ℓ::IT) where {IT, RT, NT, ST, B, FT<:Real}
    compute_block!(c.lift.inner, ℓ)
    lift!(c, ℓ)
    c.ℓ[] = ℓ
    c
end

# The block for the ``ℓ`` just computed, restricted to the m′ and m limits the calculator
# was built with.
current_block(c::WignerCalculator{IT}, ℓ::IT) where {IT} =
    block(c, ℓ, m′range(c, ℓ), mrange(c, ℓ))

# Ranges of m′ and m in the block for a given ℓ, which is what `current_block` labels, and
# the ranges of the rows and columns stored, which is what `materialize!` writes.
m′range(c::WignerCalculator, ℓ) = max(-ℓ, c.m′ₘᵢₙ):min(ℓ, c.m′ₘₐₓ)
stored_m′range(c::WignerCalculator, ℓ) = max(-ℓ, c.m′ₘᵢₙˢ):min(ℓ, c.m′ₘₐₓˢ)
stored_mrange(c::WignerCalculator, ℓ) = max(-ℓ, c.mₘᵢₙˢ):min(ℓ, c.mₘₐₓˢ)
mrange(c::WignerCalculator, ℓ) = max(-ℓ, c.mₘᵢₙ):min(ℓ, c.mₘₐₓ)

# Write the block of the d matrix (real NT) or 𝔇 matrix (complex NT) for the current ℓ into
# c.Wˡ, applying the ϵ signs relating H to d, and for 𝔇 the phases e^{-i(m′α+mγ)}.
#
# The block is written a column at a time, in the order of its storage, and the source of
# each element in the wedge is resolved by the column's runs rather than element by element
# through `wedge_source`, whose index arithmetic would otherwise cost more than the
# recurrence that computed the wedge.  For a column m, the symmetries that `wedge_source`
# applies divide the m′ axis into three parts:
#
#     -m′ > |m|     H[-m, -m′]                     row -m, one run, descending in m′
#     |m′| ≤ |m|    H[m′, m] or σ H[-m′, -m]       rows m′ or -m′, by the sign of m
#     m′ > |m|      σ H[m, m′]                     row m, one run, ascending in m′
#
# with σ = `transpose_sign(m′, m)`.  These are exactly the branches that `wedge_source`
# takes; in the middle part the element m′ = m = 0 belongs with m ≥ 0, as it does there.
# The two runs read rows |m| ≤ W, where W is the number of rows of the wedge above the axis,
# and the middle part reads rows |m′| ≤ W; every element of the block is of one kind or the
# other, because the wedge is as wide as the narrower of the m′ and m ranges (see
# `allocate_W`), so the runs are needed only for the columns |m| ≤ W, and the middle part
# never reaches beyond |m′| = W.  The result is compared with `wedge_value`, element by
# element, in `test/wigner/calculators.jl`.
#
# Nothing here depends on whether the indices are integers or half-odd-integers: the
# exponents m′±m of z₊ and z₋ are `Integer`s either way (see "Step 7" of
# `docs/src/50-notes/01-H_recurrence.md`), and ϵ and σ are already general.
function materialize!(c::WignerCalculator{IT, RT, NT}, ℓ::IT) where {IT, RT, NT}
    # A calculator of 𝔇 always holds the phases of its rotors, and one of `d` never does.
    let H = c.engine.H.Hˡ, Wˡ = c.Wˡ, Z₊ = c.engine.Z₊, Z₋ = c.engine.Z₋, Nᵣ = Nᵣ(c),
            Hp = parent(H), phases = NT <: Complex, conjugate = Val(true)
        if H.ℓ != ℓ
            error("The H wedge holds ℓ=$(H.ℓ), but ℓ=$ℓ was requested.")
        end
        W = m′ₘₐₓ(H)
        m′ₘᵢₙw = m′ₘᵢₙ(H)
        ri = row_index(H)
        m′lo, m′hi = first(stored_m′range(c, ℓ)), last(stored_m′range(c, ℓ))
        mlo, mhi = first(stored_mrange(c, ℓ)), last(stored_mrange(c, ℓ))
        # `r` is the offset of the first element of a row, H[±m, |m|]; element iᵣ of H[a, b]
        # is at `ri[(a - m′ₘᵢₙw) + 1] - 1 + Nᵣ * (b - |a|) + iᵣ`, as in `wedge_offset`.
        @inbounds for m ∈ mlo:mhi
            a = abs(m)
            j = Int(m - mlo) + 1
            if a ≤ W
                # -m′ > |m|: d = ϵ(m′) ϵ(-m) H[-m, -m′], with σ = 1
                r = ri[(-m - m′ₘᵢₙw) + 1] - 1
                for m′ ∈ m′lo:min(m′hi, -a - 1)
                    coefficient = convert(RT, ϵ(m′) * ϵ(-m))
                    materialize_element!(
                        Wˡ, Hp, Z₊, Z₋, Nᵣ, Int(m′ - m′lo) + 1, j, r + Nᵣ * Int(-m′ - a),
                        coefficient, m′ + m, m′ - m, phases, conjugate
                    )
                end
            end
            b = min(a, W)
            if m ≥ 0
                # m ≥ |m′|: the element is stored, with σ = 1
                for m′ ∈ max(m′lo, -b):min(m′hi, b)
                    coefficient = convert(RT, ϵ(m′) * ϵ(-m))
                    materialize_element!(
                        Wˡ, Hp, Z₊, Z₋, Nᵣ, Int(m′ - m′lo) + 1, j,
                        ri[(m′ - m′ₘᵢₙw) + 1] - 1 + Nᵣ * Int(m - abs(m′)),
                        coefficient, m′ + m, m′ - m, phases, conjugate
                    )
                end
            else
                # -m ≥ |m′|: σ H[-m′, -m]
                for m′ ∈ max(m′lo, -b):min(m′hi, b)
                    coefficient = convert(RT, transpose_sign(m′, m) * ϵ(m′) * ϵ(-m))
                    materialize_element!(
                        Wˡ, Hp, Z₊, Z₋, Nᵣ, Int(m′ - m′lo) + 1, j,
                        ri[(-m′ - m′ₘᵢₙw) + 1] - 1 + Nᵣ * Int(-m - abs(m′)),
                        coefficient, m′ + m, m′ - m, phases, conjugate
                    )
                end
            end
            if a ≤ W
                # m′ > |m|: σ H[m, m′]
                r = ri[(m - m′ₘᵢₙw) + 1] - 1
                for m′ ∈ max(m′lo, a + 1):m′hi
                    coefficient = convert(RT, transpose_sign(m′, m) * ϵ(m′) * ϵ(-m))
                    materialize_element!(
                        Wˡ, Hp, Z₊, Z₋, Nᵣ, Int(m′ - m′lo) + 1, j, r + Nᵣ * Int(m′ - a),
                        coefficient, m′ + m, m′ - m, phases, conjugate
                    )
                end
            end
        end
    end
    c.ℓ[] = ℓ
    c
end

# `isbatched(c)` reads the type parameter, so this branch is resolved at compile time and
# the method has a single concrete return type.  The same containers are returned for
# integer and half-odd-integer indices alike; see the note on `AbstractBlock` for why
# they are not `OffsetArray`s even where an `OffsetArray` could represent them.
#
# The blocks are built with the inner constructors, as `D_series` builds its own: the
# limits are the calculator's, which were validated when it was built, so that the public
# constructors would only repeat the normalization and validation of the indices at every ℓ.
function block(c::WignerCalculator{IT, RT, NT}, ℓ::IT, m′r, mr) where {IT<:IntegerHalf, RT, NT}
    # The block begins `o′` rows and `o` columns into what is stored; see `stored_limits`.
    o′ = Int(first(m′r) - first(stored_m′range(c, ℓ)))
    o = Int(first(mr) - first(stored_mrange(c, ℓ)))
    rows, cols = (o′ + 1):(o′ + length(m′r)), (o + 1):(o + length(mr))
    if isbatched(c)
        let p = view(c.Wˡ, :, rows, cols)
            WignerMatrixBatch{IT, NT, typeof(p)}(
                p, ℓ, last(m′r), first(m′r), last(mr), first(mr), size(p, 1)
            )
        end
    else
        let p = view(c.Wˡ, 1, rows, cols)
            WignerMatrix{IT, NT, typeof(p)}(p, ℓ, last(m′r), first(m′r), last(mr), first(mr))
        end
    end
end


### Convenience functions

"""
    D(R, ℓₘₐₓ; m′ₘₐₓ=ℓₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ, mₘₐₓ=ℓₘₐₓ, mₘᵢₙ=-mₘₐₓ)
    D(α, β, γ, ℓₘₐₓ; m′ₘₐₓ=ℓₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ, mₘₐₓ=ℓₘₐₓ, mₘᵢₙ=-mₘₐₓ)

Wigner's ``𝔇^{(ℓ)}_{m′,m}(R)`` matrices for all ``ℓ ≤ ℓₘₐₓ``, for the single rotor `R`, or
for the rotor `from_euler_angles(α, β, γ)` of the Euler angles ``(α, β, γ)``.

The result is a [`WignerSeries`](@ref), indexed by ``ℓ`` and then naturally by ``(m′, m)``:
`D(R, ℓₘₐₓ)[ℓ][m′, m]`.  Each block is a [`WignerMatrix`](@ref); `Matrix` (or `collect`)
gives the ordinary ``(2ℓ+1)×(2ℓ+1)`` `Matrix` with rows and columns in order of increasing
``m′`` and ``m``.  The keyword arguments restrict the block of each matrix that is computed,
and may also be spelled `mp_max`, `mp_min`, `m_max` and `m_min`.  The convention is
``𝔇^{(ℓ)}_{m′,m}(𝐑_{α,β,γ}) = e^{-im′α}\\, d^{(ℓ)}_{m′,m}(β)\\, e^{-imγ}``; see the
"Conventions" section of the documentation.

`ℓₘₐₓ` may be a half-integer, as a `Rational` with denominator 2 or a
[`HalfOddInteger`](@ref), as may the keyword limits; then ``ℓ`` runs over ``1/2, 3/2, …,
ℓₘₐₓ`` rather than ``0, 1, …, ℓₘₐₓ``.  The container types are the same either way, and
`D(R, ℓₘₐₓ)[ℓ][m′, m]` reads the same.

This function allocates all of its output on every call, and takes a single rotor.  To
evaluate the matrices for many rotors, or to avoid holding every ``ℓ`` at once, use a
[`DCalculator`](@ref) instead, which allocates once, computes one ``ℓ`` at a time, and takes
a vector of rotors as a batch.

See also [`d`](@ref) and [`sYlm`](@ref).
"""
@index_methods function D(
    R::RotorLike, ℓₘₐₓ::IT;
    mp_max::IndexType=ℓₘₐₓ, m′ₘₐₓ::IndexType=mp_max,
    mp_min::IndexType=-m′ₘₐₓ, m′ₘᵢₙ::IndexType=mp_min,
    m_max::IndexType=ℓₘₐₓ, mₘₐₓ::IndexType=m_max,
    m_min::IndexType=-mₘₐₓ, mₘᵢₙ::IndexType=m_min
) where {IT<:IndexType}
    D_series(D_array(R, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ), ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
end
@index_methods function D(
    α::Real, β::Real, γ::Real, ℓₘₐₓ::IndexType;
    mp_max::IndexType=ℓₘₐₓ, m′ₘₐₓ::IndexType=mp_max,
    mp_min::IndexType=-m′ₘₐₓ, m′ₘᵢₙ::IndexType=mp_min,
    m_max::IndexType=ℓₘₐₓ, mₘₐₓ::IndexType=m_max,
    m_min::IndexType=-mₘₐₓ, mₘᵢₙ::IndexType=m_min
)
    D(Quaternionic.from_euler_angles(α, β, γ), ℓₘₐₓ; m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
end

"""
    d(β, ℓₘₐₓ; m′ₘₐₓ=ℓₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ, mₘₐₓ=ℓₘₐₓ, mₘᵢₙ=-mₘₐₓ)

Wigner's real ``d^{(ℓ)}_{m′,m}(β)`` matrices for all ``ℓ ≤ ℓₘₐₓ``, for the single angle `β`,
which may also be given as the phase ``e^{iβ}`` or as a `Rotor` or other `Quaternion`.  The
result is indexed as `d(β, ℓₘₐₓ)[ℓ][m′, m]`.  See [`D`](@ref) for details, including the
keyword arguments and their ASCII spellings; this function is the real, ``β``-only analogue.

A `Rotor` contributes the ``β ∈ [0, π]`` of its canonical Euler decomposition, as described
under [`dCalculator`](@ref), so `d(R, ℓₘₐₓ)` equals `d(β, ℓₘₐₓ)` for the angle ``β`` used to
build `R` only when ``β ∈ [0, π]``.

`ℓₘₐₓ` may be a half-integer, as a `Rational` with denominator 2 or a
[`HalfOddInteger`](@ref).  Half-integer ``d`` has period ``4π`` in ``β``, so an angle or a
`Rotor` determines it unambiguously, but a bare phase ``e^{iβ}`` determines it only up to
the sign ``(-1)^{2ℓ}`` (the branch ``β ∈ (-π, π]`` is used).
"""
@index_methods function d(
    β::Union{Real, Complex, RotorLike}, ℓₘₐₓ::IT;
    mp_max::IndexType=ℓₘₐₓ, m′ₘₐₓ::IndexType=mp_max,
    mp_min::IndexType=-m′ₘₐₓ, m′ₘᵢₙ::IndexType=mp_min,
    m_max::IndexType=ℓₘₐₓ, mₘₐₓ::IndexType=m_max,
    m_min::IndexType=-mₘₐₓ, mₘᵢₙ::IndexType=m_min
) where {IT<:IndexType}
    calc = dCalculator(β, ℓₘₐₓ; m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    WignerSeries(
        [copy(recurrence!(calc, ℓ)) for ℓ ∈ lowest_index(IT):ℓₘₐₓ], lowest_index(IT), ℓₘₐₓ
    )
end

# The values of `D(R, ℓₘₐₓ)` are computed as one matrix for each ℓ, from ℓₘᵢₙ up, which are
# then labelled without being copied.  This split is for automatic differentiation:
# `D_array` takes the rotor and returns plain arrays, which every tool can handle, so it is
# the function to which the extensions for ChainRulesCore and ReverseDiff attach their rules
# (see `src/derivatives/kernels.jl`); the other tools differentiate it through the
# calculator.
function D_array(
    R::RotorLike, ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT
) where {IT<:IntegerHalf}
    calc = DCalculator(R, ℓₘₐₓ; m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    [copy(parent(recurrence!(calc, ℓ))) for ℓ ∈ lowest_index(IT):ℓₘₐₓ]
end

function D_series(
    blocks::AbstractVector{<:AbstractMatrix{NT}}, ℓₘₐₓ::IT,
    m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT
) where {NT, IT<:IntegerHalf}
    series = map(enumerate(lowest_index(IT):ℓₘₐₓ)) do (i, ℓ)
        p = blocks[i]
        m′r, mr = max(-ℓ, m′ₘᵢₙ):min(ℓ, m′ₘₐₓ), max(-ℓ, mₘᵢₙ):min(ℓ, mₘₐₓ)
        WignerMatrix{IT, NT, typeof(p)}(p, ℓ, last(m′r), first(m′r), last(mr), first(mr))
    end
    WignerSeries(series, lowest_index(IT), ℓₘₐₓ)
end

# `D` and `d` take one rotor.  A vector of them is what a calculator is for, and is refused
# with that advice rather than with a bare `MethodError`; so is a `QuatVec`, which denotes
# a vector rather than a rotation, with the advice of `not_a_rotor`.
function D(R⃗::AbstractVector{<:RotorLike}, ℓₘₐₓ; kwargs...)
    throw(ArgumentError(
        "`D` takes a single rotor; for a vector of $(length(R⃗)) "
        * "rotor$(length(R⃗) == 1 ? "" : "s") use "
        * "`DCalculator(R⃗, ℓₘₐₓ)`, which computes all of them at once, one ℓ at a time, "
        * "in blocks indexed [iᵣ, m′, m]."
    ))
end
function d(β⃗::AbstractVector{<:Union{Real, Complex, RotorLike}}, ℓₘₐₓ; kwargs...)
    throw(ArgumentError(
        "`d` takes a single angle, phase or rotor; for a vector of $(length(β⃗)) of them use "
        * "`dCalculator(β⃗, ℓₘₐₓ)`, which computes all of them at once, one ℓ at a time, in "
        * "blocks indexed [iᵣ, m′, m]."
    ))
end
D(R::NonRotorData, ℓₘₐₓ; kwargs...) = throw(ArgumentError(not_a_rotor(R)))
d(R::NonRotorData, ℓₘₐₓ; kwargs...) = throw(ArgumentError(not_a_rotor(R)))
