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

A calculator whose rotor data hold derivatives, such as dual numbers, runs the recurrence on
the values of those data, in `FT`, and gives each block its derivatives from the
angular-momentum operators (for ``d``, from the one about the ``y`` axis), as described in
`src/derivatives/kernels.jl`; the recurrence itself is never differentiated.

Because `B` is a type parameter, the type of the block is known at compile time, and a loop
over the blocks is inferrable.
"""
struct WignerCalculator{
    IT, RT<:Real, NT<:Union{RT, Complex{RT}}, ST, B, FT<:Real, L
} <: AbstractCalculator{IT}
    engine::SphericalFunctionsEngine{IT, FT, ST}  # the recurrence and the power tables
    Wˡ::Vector{NT}  # the block for the current ℓ, dense as [iᵣ, m′, m], in its leading entries
    rotors::Vector{Quaternion{RT}}  # the rotors themselves; empty when NT is real
    angles::Vector{RT}  # the angles β of the rotor data; empty when NT is complex, unused for floats
    m′ₘₐₓ::IT
    m′ₘᵢₙ::IT
    mₘₐₓ::IT
    mₘᵢₙ::IT
    m′ₘₐₓᵈ::IT  # the limits of the rows and columns from which the derivatives of a block
    m′ₘᵢₙᵈ::IT  # are computed, which extend one beyond those of the block along one axis
    mₘₐₓᵈ::IT   # when the block is restricted in both m′ and m; see `derivative_limits`
    mₘᵢₙᵈ::IT
    ℓ::Base.RefValue{IT}  # ℓ of the block currently in Wˡ; ℓₘᵢₙ-1 if none
    lift::L  # `nothing`, or the calculator of the rotors' values and their generators
    # `materialize!` writes `Wˡ` and reads the power tables under `@inbounds`, for every
    # rotor of the engine and every (m′, m) within the limits, so the buffers must be large
    # enough for those; as for `HCalculator`, this checks them once, as they are brought
    # together.  The block of the largest ℓ is the largest, and `Wˡ` must hold it, unless it
    # is empty: a calculator that writes every block into a destination given to it holds no
    # block of its own (see `allocate_W`).  A step of such a calculator into `Wˡ` is
    # refused, by `check_destination` in `materialize!`, or, where the calculator lifts the
    # blocks of another, by the bounds check of the view of `Wˡ` into which `lift!` writes.
    function WignerCalculator{IT, RT, NT, ST, B, FT, L}(
        engine, Wˡ, rotors, angles, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, m′ₘₐₓᵈ, m′ₘᵢₙᵈ, mₘₐₓᵈ, mₘᵢₙᵈ, ℓ,
        lift
    ) where {IT, RT<:Real, NT<:Union{RT, Complex{RT}}, ST, B, FT<:Real, L}
        let n = Nᵣ(engine),
                K = NT <: Complex ? 2power_extent(ℓₘₐₓ(engine), engine.H.m′ₘₐₓ) + 1 : 0,
                Z₊ = engine.Z₊, Z₋ = engine.Z₋
            if !(
                (isempty(Wˡ)
                 || length(Wˡ) ≥ n * (Int(m′ₘₐₓ - m′ₘᵢₙ) + 1) * (Int(mₘₐₓ - mₘᵢₙ) + 1))
                && size(Z₊, 1) ≥ n && size(Z₊, 2) ≥ K && size(Z₋, 1) ≥ n && size(Z₋, 2) ≥ K
                && length(rotors) == (NT <: Complex ? n : 0)
                && length(angles) == (NT <: Complex ? 0 : n)
                && m′ₘₐₓᵈ ≥ m′ₘₐₓ && m′ₘᵢₙᵈ ≤ m′ₘᵢₙ && mₘₐₓᵈ ≥ mₘₐₓ && mₘᵢₙᵈ ≤ mₘᵢₙ
            )
                throw(DimensionMismatch(
                    "The buffers of a WignerCalculator for Nᵣ=$n rotors, m′ ∈ $m′ₘᵢₙ:$m′ₘₐₓ "
                    * "(differentiated from $m′ₘᵢₙᵈ:$m′ₘₐₓᵈ) and m ∈ $mₘᵢₙ:$mₘₐₓ "
                    * "(differentiated from $mₘᵢₙᵈ:$mₘₐₓᵈ) are too small: the "
                    * "block has length $(length(Wˡ)), the power tables $(size(Z₊)) and "
                    * "$(size(Z₋)), which need a row for each rotor and at least $K columns, "
                    * "and there are $(length(rotors)) rotors and $(length(angles)) angles."
                ))
            end
        end
        new{IT, RT, NT, ST, B, FT, L}(
            engine, Wˡ, rotors, angles, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, m′ₘₐₓᵈ, m′ₘᵢₙᵈ, mₘₐₓᵈ,
            mₘᵢₙᵈ, ℓ, lift
        )
    end
end

# The limits of the rows and columns from which the derivatives of the block of a calculator
# of 𝔇 or of `d` are computed.  The derivatives, which the rules for automatic
# differentiation give in terms of the values (see `src/derivatives/kernels.jl`), couple
# each element either to its neighbors in the same column (from the left) or to those in the
# same row (from the right).  A block whose rows are all of -ℓ:ℓ is differentiated from the
# left, and one whose columns are, from the right, so that every neighbor needed is in the
# block already.  Only a block restricted in both needs values beyond its limits: one row or
# column more on each side, within ±ℓₘₐₓ, along whichever axis is the wider, so that the
# wedge, whose width is that of the narrower (see `allocate_W`), is widened only when the
# two are equally wide.  The widened limits still bracket ±ℓₘᵢₙ, as the originals do.  Those
# values are not kept in the calculator's block: a calculator that lifts the blocks of
# another is given one whose own limits are these (see `allocate_W`), and the rules for the
# other tools materialize them from the wedge (see `derivative_values`).
function derivative_limits(ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT) where {IT}
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
# `derivative_limits`.  A block of `d` whose columns are all of -ℓ:ℓ is differentiated from
# the right even when its rows are too: the coefficients of the derivative from the right
# are the same down a column, the order in which a block is stored, along which the kernels
# then run fastest (see `src/derivatives/kernels.jl`).
function derivatives_from_left(c::WignerCalculator{IT, RT, NT}) where {IT, RT, NT}
    full_rows = c.m′ₘₐₓ == ℓₘₐₓ(c) && c.m′ₘᵢₙ == -ℓₘₐₓ(c)
    full_columns = c.mₘₐₓ == ℓₘₐₓ(c) && c.mₘᵢₙ == -ℓₘₐₓ(c)
    if full_columns && (NT <: Real || !full_rows)
        false
    elseif full_rows
        true
    else
        c.m′ₘₐₓᵈ != c.m′ₘₐₓ || c.m′ₘᵢₙᵈ != c.m′ₘᵢₙ
    end
end

# Allocate the buffers without touching them.  PRIVATE: see the note on `allocate_H`.  The
# rotor data here is the engine's (see `copy_rotor_data!`) *plus* the rotors or angles
# themselves, so a caller that copies rather than sets must copy both, as
# `copy_rotor_state!` does.  A calculator whose real type holds derivatives is given a
# calculator of the rotors' values, whose engine it holds as its own, and whose blocks it
# lifts; see `allocate_lift`.
#
# The block buffer `Wˡ` has `block_length` entries, by default as many as the block of the
# largest ℓ has, so that every block fits in it.  The calculator of `D` and `d` is given 0
# (see `series_calculator`), since `wigner_arrays` writes each of its blocks into one new
# vector, and nothing would read `Wˡ`.  The calculator of values inside a calculator that
# lifts its blocks is always given the whole buffer, because `lift!` reads the blocks there.
function allocate_W(
    ::Type{IT}, ::Type{RT}, ::Type{NT}, ℓₘₐₓ::IT,
    m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT, Nᵣ::Int, ::Val{B},
    block_length::Int=Nᵣ * (Int(m′ₘₐₓ - m′ₘᵢₙ) + 1) * (Int(mₘₐₓ - mₘᵢₙ) + 1)
) where {IT<:IntegerHalf, RT<:Real, NT<:Union{RT, Complex{RT}}, B}
    validate_degree(ℓₘₐₓ)
    validate_axis(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, "m′")
    validate_axis(ℓₘₐₓ, mₘₐₓ, mₘᵢₙ, "m")
    m′ₘₐₓᵈ, m′ₘᵢₙᵈ, mₘₐₓᵈ, mₘᵢₙᵈ = derivative_limits(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    Wˡ = Vector{NT}(undef, block_length)
    rotors = Vector{Quaternion{RT}}(undef, NT <: Complex ? Nᵣ : 0)
    angles = Vector{RT}(undef, NT <: Complex ? 0 : Nᵣ)
    # The field is a `RefValue{IT}`, so the type is given explicitly, as in `allocate_Y`.
    ℓ = Ref{IT}(lowest_index(IT) - 1)
    if value_type(RT) !== RT
        let inner = allocate_W(
            IT, value_type(RT), NT <: Complex ? Complex{value_type(RT)} : value_type(RT),
            ℓₘₐₓ, m′ₘₐₓᵈ, m′ₘᵢₙᵈ, mₘₐₓᵈ, mₘᵢₙᵈ, Nᵣ, Val(B)
        )
            lift = allocate_lift(RT, NT, inner, Nᵣ)
            WignerCalculator{IT, RT, NT, typeof(parent(inner.engine.H.Hˡ)), B, recurrence_type(RT), typeof(lift)}(
                inner.engine, Wˡ, rotors, angles,
                m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, m′ₘₐₓᵈ, m′ₘᵢₙᵈ, mₘₐₓᵈ, mₘᵢₙᵈ, ℓ, lift
            )
        end
    else
        # The wedge holds the rows |m′| ≤ W of H, for every m, and `materialize!` reads an
        # element whose |m′| exceeds W from its transpose, whose |m| does not.  So W need
        # only be the narrower of the widest m′ and the widest m within the derivative
        # limits, and restricting either the rows or the columns of the block to a narrow
        # band narrows the recurrence with it.  Both limits bracket ±ℓₘᵢₙ, so W does too, as
        # the recurrence requires.
        W = min(max(m′ₘₐₓᵈ, -m′ₘᵢₙᵈ), max(mₘₐₓᵈ, -mₘᵢₙᵈ))
        engine = allocate_engine(IT, RT, ℓₘₐₓ, W, Nᵣ, NT <: Complex)
        WignerCalculator{IT, RT, NT, typeof(parent(engine.H.Hˡ)), B, RT, Nothing}(
            engine, Wˡ, rotors, angles,
            m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, m′ₘₐₓᵈ, m′ₘᵢₙᵈ, mₘₐₓᵈ, mₘᵢₙᵈ, ℓ, nothing
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
const dCalculator{IT, RT, ST, B} = WignerCalculator{IT, RT, RT, ST, B} where {IT, RT<:Real, ST, B}

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
# in place.  The angles are kept too (see `store_angles!`), which is what the rules for
# automatic differentiation read; as for 𝔇, the extensions for Enzyme and Mooncake declare
# `set_rotor_data!` to have no derivatives, so that those tools differentiate only the copy
# of the angles.  A calculator that lifts the blocks of another gives that calculator the
# values of the rotor data, and computes the generators of its angles' derivatives.
function set_rotors!(c::WignerCalculator{IT, RT, RT, ST, B, FT, Nothing}, R) where {IT, RT<:Real, ST, B, FT<:Real}
    set_rotor_data!(c, R)
    store_angles!(c.angles, R)
    c
end
function set_rotor_data!(c::WignerCalculator{IT, RT, RT}, R) where {IT, RT<:Real}
    set_rotors!(c.engine.H, R)
    c.ℓ[] = lowest_index(IT) - 1
    nothing
end
function set_rotors!(c::WignerCalculator{IT, RT, RT, ST, B, FT, <:Lift}, R) where {IT, RT<:Real, ST, B, FT<:Real}
    set_rotors!(c.lift.inner, LiftedValues(R))
    c.ℓ[] = lowest_index(IT) - 1
    store_angles!(c.angles, R)
    set_generators!(c.lift, derivatives_from_left(c), c.angles)
    c
end
function set_rotors!(
    c::WignerCalculator{IT, RT, RT, ST, B, FT, <:Lift}, R::Union{Real, Complex, RotorLike}
) where {IT, RT<:Real, ST, B, FT<:Real}
    check_rotor_count(c, R)
    set_rotors!(c, @SVector [R])
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

# Compute the block of degree ℓ into `A`, dense as [iᵣ, m′, m] after its first `o` entries:
# one step of the recurrence and the materialization of its result, or, for a calculator
# that lifts the blocks of the calculator of its rotors' values, that calculator's step and
# the lifting.  `A` is the calculator's own `Wˡ` with `o = 0`, or another array with linear
# indexing, such as the vector into which `wigner_arrays` writes every block of `D` and `d`,
# so that the values reach their destination without passing through `Wˡ`; the calculator
# then no longer holds a block.  Every sweep of a calculator is made of these steps, and
# they are what the extensions for Enzyme and Mooncake attach their rules to, so that the
# recurrence itself is never differentiated.
compute_block!(c::WignerCalculator, ℓ) = compute_block!(c, ℓ, c.Wˡ, 0)
function compute_block!(
    c::WignerCalculator{IT, RT, NT, ST, B, FT, Nothing}, ℓ::IT, A, o::Int
) where {IT, RT, NT, ST, B, FT<:Real}
    recurrence!(c.engine.H, ℓ)
    materialize!(c, ℓ, m′range(c, ℓ), mrange(c, ℓ), A, o)
    c.ℓ[] = holds_block(c, A, o) ? ℓ : lowest_index(IT) - 1
    c
end
function compute_block!(
    c::WignerCalculator{IT, RT, NT, ST, B, FT, <:Lift}, ℓ::IT, A, o::Int
) where {IT, RT, NT, ST, B, FT<:Real}
    compute_block!(c.lift.inner, ℓ)
    lift!(c, ℓ, A, o)
    c.ℓ[] = holds_block(c, A, o) ? ℓ : lowest_index(IT) - 1
    c
end
holds_block(c::WignerCalculator, A, o) = A === c.Wˡ && o == 0

# Ranges of m′ and m in the block for a given ℓ, which is what `current_block` labels, and
# the ranges from which the derivatives of that block are computed (see
# `derivative_limits`).
m′range(c::WignerCalculator, ℓ) = max(-ℓ, c.m′ₘᵢₙ):min(ℓ, c.m′ₘₐₓ)
derivative_m′range(c::WignerCalculator, ℓ) = max(-ℓ, c.m′ₘᵢₙᵈ):min(ℓ, c.m′ₘₐₓᵈ)
derivative_mrange(c::WignerCalculator, ℓ) = max(-ℓ, c.mₘᵢₙᵈ):min(ℓ, c.mₘₐₓᵈ)
mrange(c::WignerCalculator, ℓ) = max(-ℓ, c.mₘᵢₙ):min(ℓ, c.mₘₐₓ)

# The block of degree ℓ after the first `o` entries of the destination `A` of
# `compute_block!`, as a 3-dimensional array [iᵣ, m′, m], batched or not.  The lifting of a
# block and the rules for automatic differentiation read and write blocks through this.
function block_array(c::WignerCalculator, A::AbstractArray, ℓ, o::Int=0)
    dims = (Nᵣ(c), length(m′range(c, ℓ)), length(mrange(c, ℓ)))
    reshape(view(A, (o + 1):(o + prod(dims))), dims)
end

# The values from which the derivatives of the block of degree ℓ are computed, over
# `derivative_m′range` and `derivative_mrange`, as a 3-dimensional array [iᵣ, m′, m], just
# after the step that computed the block into `A` after its first `o` entries: the block
# itself, where those are the ranges of the block, and otherwise a new array, into which the
# wedge of ℓ, which the engine still holds, is materialized over the wider ranges.
function derivative_values(
    c::WignerCalculator{IT, RT, NT}, ℓ::IT, A::AbstractArray=c.Wˡ, o::Int=0
) where {IT, RT, NT}
    rows, cols = derivative_m′range(c, ℓ), derivative_mrange(c, ℓ)
    if rows == m′range(c, ℓ) && cols == mrange(c, ℓ)
        block_array(c, A, ℓ, o)
    else
        A = Vector{NT}(undef, Nᵣ(c) * length(rows) * length(cols))
        materialize!(c, ℓ, rows, cols, A, 0)
        reshape(view(A, 1:length(A)), Nᵣ(c), length(rows), length(cols))
    end
end

# Write the rows `m′r` and columns `mr` of the block of the d matrix (real NT) or 𝔇 matrix
# (complex NT) for the current ℓ into `A`, dense as [iᵣ, m′, m] after its first `o` entries,
# applying the ϵ signs relating H to d, and for 𝔇 the phases e^{-i(m′α+mγ)}.  The ranges
# are those of the block, or, for the rules for automatic differentiation, the wider ones of
# `derivative_limits`, which the wedge reaches; either way they include ±ℓₘᵢₙ.
#
# The block is written as runs down its columns and along its rows, each of which reads one
# row of the wedge in order, rather than element by element through `wedge_source`, whose
# index arithmetic would otherwise cost more than the recurrence that computed the wedge.
# The symmetries that `wedge_source` applies divide the block into four parts:
#
#     -m′ > |m|           H[-m, -m′]      down the column m, along the row -m, descending
#     m′ > |m|            σ H[m, m′]      down the column m, along the row m, ascending
#     m ≥ |m′|            H[m′, m]        along the row m′, along the row m′, ascending
#     -m ≥ |m′|, m < 0    σ H[-m′, -m]    along the row m′, along the row -m′, descending
#
# with σ = `transpose_sign(m′, m)`.  These are exactly the branches that `wedge_source`
# takes, and the element m′ = m = 0 belongs with m ≥ 0, as it does there.  The runs down the
# columns read the rows |m| ≤ W of the wedge, where W is the number of its rows above the
# axis, and those along the rows read the rows |m′| ≤ W; every element of the block is in
# one part or another, because the wedge is as wide as the narrower of the m′ and m ranges
# (see `allocate_W`).  Within a run the sign ϵ(m′) ϵ(-m) σ is constant or alternates, and
# the exponents m′ ± m of the phases each step by ±1, so every run is a loop with no branch
# in it (see `materialize_run!`).  The runs along the rows store to the block at a stride;
# they were measured to be faster than filling each column in the order of its storage.
#
# The element (-m′, -m) is read from the same element of the wedge as (m′, m), from the
# first part if (m′, m) is in the second, and from the fourth if it is in the third, and its
# phase is the conjugate of that of (m′, m).  So the part of a block of 𝔇 within the
# symmetric limits -M′:M′ by -M:M is written as runs of the first and third parts, each of
# which also writes the partners of its elements, with the one phase product of each element
# serving both; for integer indices the element (0, 0) is its own partner, and is written
# alone.  The rest of the block, beyond ±M′ or ±M, is written as runs of each of the four
# parts.  The conjugate of a product and the product of the conjugates agree in value but
# not always in the signs of their zeros, so a partner may hold -0.0 where its own product
# would give 0.0.  A block of `d` is written the same way, its partners with their
# coefficients and no phases, except for one rotor at large ℓ (see below).  The result is
# compared with `wedge_value`, element by element, in `test/wigner/calculators.jl`.
#
# Nothing here depends on whether the indices are integers or half-odd-integers: the
# exponents m′±m of z₊ and z₋ are `Integer`s either way (see "Step 7" of
# `docs/src/50-notes/01-H_recurrence.md`), and ϵ and σ are already general.
function materialize!(
    c::WignerCalculator{IT, RT, NT}, ℓ::IT, m′r::AbstractUnitRange{IT},
    mr::AbstractUnitRange{IT}, A::AbstractArray, o::Int
) where {IT, RT, NT}
    if c.engine.H.Hˡ.ℓ != ℓ
        error("The H wedge holds ℓ=$(c.engine.H.Hˡ.ℓ), but ℓ=$ℓ was requested.")
    end
    check_destination(c, ℓ, m′r, mr, A, o)
    # One rotor is written by a method in which `Nᵣ` is `Val(1)`, so that the compiler sees
    # the unit strides of its runs, and several by one in which it is the `Int` it is.
    if Nᵣ(c) == 1
        materialize!(c, ℓ, m′r, mr, A, o, Val(1))
    else
        materialize!(c, ℓ, m′r, mr, A, o, Nᵣ(c))
    end
    nothing
end
function materialize!(
    c::WignerCalculator{IT, RT, NT}, ℓ::IT, m′r::AbstractUnitRange{IT},
    mr::AbstractUnitRange{IT}, A::AbstractArray, o::Int, Nᵣ::Union{Val{1}, Int}
) where {IT, RT, NT}
    # A calculator of 𝔇 always holds the phases of its rotors, and one of `d` never does.
    let H = c.engine.H.Hˡ, Hp = parent(H), Z₊ = c.engine.Z₊, Z₋ = c.engine.Z₋,
            n = rotor_count(Nᵣ), phases = Val(NT <: Complex), conjugate = Val(true),
            W = m′ₘₐₓ(H), m′ₘᵢₙw = m′ₘᵢₙ(H),
            m′lo = first(m′r), m′hi = last(m′r), mlo = first(mr), mhi = last(mr)
        # The offset, before the first rotor, of the element (m′, m) in `A`; the steps from
        # one element of a run to the next, down a column and along a row of the block; and
        # the coefficient σ ϵ(m′) ϵ(-m) of an element.
        inblock(m′, m) = o + n * (Int(m′ - m′lo) + length(m′r) * Int(m - mlo))
        down, along = n, n * length(m′r)
        coefficient(σ, m′, m) = convert(RT, σ * ϵ(m′) * ϵ(-m))
        # A run of `len` elements from (m′, m), at the step Δ in the block, reading the
        # wedge from H[a, b] at the step Δh, with the exponents m′ + m and m′ - m of the
        # phases stepping by Δ₊ and Δ₋, and the partners `partners` (see
        # `materialize_run!`).
        @inline run!(len, (m′, m), Δ, partners, (a, b), Δh, (Δ₊, Δ₋), coefficients) =
            materialize_run!(
                A, Hp, Z₊, Z₋, Nᵣ, Int(len), (inblock(m′, m), Δ), partners,
                (wedge_offset(H, a, b, m′ₘᵢₙw), n * Δh),
                (power_offset(Z₊, m′ + m, n), n * Δ₊), (power_offset(Z₋, m′ - m, n), n * Δ₋),
                coefficients, phases, conjugate
            )
        # The four parts, without partners, within the rows m′s and the columns ms.  The
        # part -m ≥ |m′| is written before the runs down the columns, and the part m ≥ |m′|
        # after them, so that the elements that read any one element of the wedge are
        # written in the order of their storage in the block.
        @inline function unpartnered!(m′s, ms)
            (isempty(m′s) || isempty(ms)) && return nothing
            for m′ ∈ max(first(m′s), -W):min(last(m′s), W)
                a′ = abs(m′)
                # -m ≥ |m′| with m < 0: σ H[-m′, -m], where ϵ(-m) alternates.  For integer
                # indices -|m′| is 0 when m′ = 0, and that element belongs with the part
                # m ≥ |m′|.
                lo, hi = first(ms), min(last(ms), -a′)
                if IT <: Integer && iszero(a′)
                    hi = min(hi, -one(hi))
                end
                if lo ≤ hi
                    run!(
                        hi - lo + 1, (m′, lo), along, nothing, (-m′, -lo), -1, (1, -1),
                        (coefficient(transpose_sign(m′, lo), m′, lo),
                         coefficient(transpose_sign(m′, lo + 1), m′, lo + 1))
                    )
                end
            end
            for m ∈ max(first(ms), -W):min(last(ms), W)
                a = abs(m)
                # -m′ > |m|: ϵ(m′) ϵ(-m) H[-m, -m′], with ϵ(m′) = 1 and σ = 1
                lo, hi = first(m′s), min(last(m′s), -a - 1)
                if lo ≤ hi
                    cₗ = coefficient(1, lo, m)
                    run!(hi - lo + 1, (lo, m), down, nothing, (-m, -lo), -1, (1, 1), (cₗ, cₗ))
                end
                # m′ > |m|: σ H[m, m′], where ϵ(m′) alternates
                lo, hi = max(first(m′s), a + 1), last(m′s)
                if lo ≤ hi
                    run!(
                        hi - lo + 1, (lo, m), down, nothing, (m, lo), 1, (1, 1),
                        (coefficient(transpose_sign(lo, m), lo, m),
                         coefficient(transpose_sign(lo + 1, m), lo + 1, m))
                    )
                end
            end
            for m′ ∈ max(first(m′s), -W):min(last(m′s), W)
                a′ = abs(m′)
                # m ≥ |m′|: H[m′, m], with ϵ(-m) = 1 and σ = 1
                lo, hi = max(first(ms), a′), last(ms)
                if lo ≤ hi
                    cₗ = coefficient(1, m′, lo)
                    run!(hi - lo + 1, (m′, lo), along, nothing, (m′, lo), 1, (1, -1), (cₗ, cₗ))
                end
            end
            nothing
        end
        # A block of `d` of one rotor at ℓ ≥ 48 is written without partners, which was
        # measured to be faster there (by 9% at ℓₘₐₓ = 64 and 12% at 128), although partners
        # are faster below that and for every batch.
        NT <: Complex || Nᵣ isa Int || ℓ < 48 || return unpartnered!(m′r, mr)
        # The symmetric part, -M′:M′ by -M:M.  Down the column m, the part -m′ > |m|, whose
        # partners are the part m′ > |m| of the column -m, σ H[-m, -m′] with ϵ(-m′)
        # alternating.
        M′, M = min(m′hi, -m′lo), min(mhi, -mlo)
        for m ∈ max(-M, -W):min(M, W)
            a = abs(m)
            lo, hi = -M′, -a - 1
            if lo ≤ hi
                cₗ = coefficient(1, lo, m)
                partners = (
                    inblock(-lo, -m), -down,
                    (coefficient(transpose_sign(-lo, -m), -lo, -m),
                     coefficient(transpose_sign(-lo - 1, -m), -lo - 1, -m))
                )
                run!(hi - lo + 1, (lo, m), down, partners, (-m, -lo), -1, (1, 1), (cₗ, cₗ))
            end
        end
        # Along the row m′, the part m ≥ |m′|, whose partners are the part -m ≥ |m′| of the
        # row -m′, σ H[m′, m] with ϵ(m) alternating.
        for m′ ∈ max(-M′, -W):min(M′, W)
            a′ = abs(m′)
            lo, hi = a′, M
            if IT <: Integer && iszero(a′) && lo ≤ hi
                cₗ = coefficient(1, 0, 0)
                run!(1, (m′, lo), along, nothing, (m′, lo), 1, (1, -1), (cₗ, cₗ))
                lo += 1
            end
            if lo ≤ hi
                cₗ = coefficient(1, m′, lo)
                partners = (
                    inblock(-m′, -lo), -along,
                    (coefficient(transpose_sign(-m′, -lo), -m′, -lo),
                     coefficient(transpose_sign(-m′, -lo - 1), -m′, -lo - 1))
                )
                run!(hi - lo + 1, (m′, lo), along, partners, (m′, lo), 1, (1, -1), (cₗ, cₗ))
            end
        end
        # The rest of the block: the rows beyond ±M′, then the columns beyond ±M.  These are
        # a loop rather than four calls, so that `unpartnered!` is inlined once: four copies
        # of it took 0.6 s to compile for each new type, against 0.25 s, and ran no faster.
        margins = (
            (m′lo:(-M′ - 1), mr), ((M′ + 1):m′hi, mr), (-M′:M′, mlo:(-M - 1)), (-M′:M′, (M + 1):mhi)
        )
        for (m′s, ms) ∈ margins
            unpartnered!(m′s, ms)
        end
    end
    nothing
end

# `materialize!` writes its destination under `@inbounds`, and reads the wedge only within
# the derivative limits, which it is as wide as.  Both ranges must also include ±ℓₘᵢₙ, as
# the ranges of every block and of the derivative limits do, so that the kernel may rely on
# it.
function check_destination(c::WignerCalculator{IT}, ℓ, m′r, mr, A, o) where {IT}
    Base.require_one_based_indexing(A)
    n = Nᵣ(c) * length(m′r) * length(mr)
    brackets(r) = first(r) ≤ -lowest_index(IT) && lowest_index(IT) ≤ last(r)
    if !(
        o ≥ 0 && o + n ≤ length(A) && brackets(m′r) && brackets(mr)
        && first(m′r) ≥ first(derivative_m′range(c, ℓ)) && last(m′r) ≤ last(derivative_m′range(c, ℓ))
        && first(mr) ≥ first(derivative_mrange(c, ℓ)) && last(mr) ≤ last(derivative_mrange(c, ℓ))
    )
        error(
            "The rows $m′r and columns $mr of the block of ℓ=$ℓ for Nᵣ=$(Nᵣ(c)) were "
            * "requested after $o entries of an array of length $(length(A)), from a "
            * "calculator that computes the rows $(derivative_m′range(c, ℓ)) and columns "
            * "$(derivative_mrange(c, ℓ)); each range must include ±$(lowest_index(IT))."
        )
    end
    nothing
end

# The block for the ``ℓ`` just computed, restricted to the m′ and m limits the calculator
# was built with, over the calculator's buffer, which holds it in its leading entries.
# `isbatched(c)` reads the type parameter, so this branch is resolved at compile time and
# the method has a single concrete return type.  The same containers are returned for
# integer and half-odd-integer indices alike; see the note on `AbstractBlock` for why they
# are not `OffsetArray`s even where an `OffsetArray` could represent them.
#
# The blocks are built with the inner constructors, as `D_series` builds its own: the
# limits are the calculator's, which were validated when it was built, so that the public
# constructors would only repeat the normalization and validation of the indices at every ℓ.
function current_block(c::WignerCalculator{IT, RT, NT}, ℓ::IT) where {IT<:IntegerHalf, RT, NT}
    m′r, mr = m′range(c, ℓ), mrange(c, ℓ)
    if isbatched(c)
        WignerMatrixBatch{IT, NT, Vector{NT}}(
            c.Wˡ, ℓ, last(m′r), first(m′r), last(mr), first(mr), Nᵣ(c)
        )
    else
        WignerMatrix{IT, NT, Vector{NT}}(c.Wˡ, ℓ, last(m′r), first(m′r), last(mr), first(mr))
    end
end


### Convenience functions

"""
    D(R, ℓₘₐₓ; m′ₘₐₓ=ℓₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ, mₘₐₓ=ℓₘₐₓ, mₘᵢₙ=-mₘₐₓ)
    D(α, β, γ, ℓₘₐₓ; m′ₘₐₓ=ℓₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ, mₘₐₓ=ℓₘₐₓ, mₘᵢₙ=-mₘₐₓ)

Wigner's ``𝔇^{(ℓ)}_{m',m}(R)`` matrices for all ``ℓ ≤ ℓₘₐₓ``, for the single rotor `R`, or
for the rotor `from_euler_angles(α, β, γ)` of the Euler angles ``(α, β, γ)``.

The result is a [`WignerSeries`](@ref), indexed by ``ℓ`` and then naturally by ``(m', m)``:
`D(R, ℓₘₐₓ)[ℓ][m′, m]`.  Each block is a [`WignerMatrix`](@ref); `Matrix` (or `collect`)
gives the ordinary ``(2ℓ+1)×(2ℓ+1)`` `Matrix` with rows and columns in order of increasing
``m'`` and ``m``.  The keyword arguments restrict the block of each matrix that is computed,
and may also be spelled `mp_max`, `mp_min`, `m_max` and `m_min`.  The convention is
``𝔇^{(ℓ)}_{m',m}(𝐑_{α,β,γ}) = e^{-im'α}\\, d^{(ℓ)}_{m',m}(β)\\, e^{-imγ}``; see the
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
    D_series(d_array(β, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ), ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
end

# The values of `D(R, ℓₘₐₓ)` are computed as one matrix for each ℓ, from ℓₘᵢₙ up, over
# consecutive parts of one vector, and are then labelled without being copied.  This split
# is for automatic differentiation: `D_array` takes the rotor and returns plain arrays,
# which every tool can handle, so it is the function to which the extensions for
# ChainRulesCore and ReverseDiff attach their rules (see `src/derivatives/kernels.jl`); the
# other tools differentiate it through the calculator.
function D_array(
    R::RotorLike, ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT
) where {IT<:IntegerHalf}
    wigner_arrays(
        series_calculator(Complex{floattype(R)}, R, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    )
end

# The same for `d` of an angle, a phase, or a rotor.
function d_array(
    β::Union{Real, Complex, RotorLike}, ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT
) where {IT<:IntegerHalf}
    wigner_arrays(series_calculator(floattype(β), β, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ))
end

# The calculator with which `D` and `d` compute their blocks, as do the rules for `D_array`:
# the one that `DCalculator` builds for the single rotor `x` when `NT` is complex, or that
# `dCalculator` builds for the single angle, phase, or rotor `x` when `NT` is real, but
# without a block buffer of its own (see `allocate_W`), since `wigner_arrays` writes every
# block into one new vector.
function series_calculator(
    ::Type{NT}, x, ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT
) where {NT, IT}
    set_rotors!(
        allocate_W(IT, real(NT), NT, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, 1, Val(false), 0), x
    )
end

# The blocks of every ℓ of a calculator of a single rotor, from its ℓₘᵢₙ up, as matrices
# over consecutive parts of one new vector, into which the calculator writes them directly.
# After each step, `f(ℓ, A, o)` is called with that vector `A` and the offset `o` of the
# block of ℓ in it, which is how the rules for automatic differentiation keep what they need
# of each step (see `wigner_arrays_with_derivative_values`).
wigner_arrays(calc::WignerCalculator) = wigner_arrays(Returns(nothing), calc)
function wigner_arrays(f::F, calc::WignerCalculator{IT, RT, NT}) where {F, IT, RT, NT}
    ℓs = lowest_index(IT):ℓₘₐₓ(calc)
    dims(ℓ) = (length(m′range(calc, ℓ)), length(mrange(calc, ℓ)))
    A = Vector{NT}(undef, sum(prod ∘ dims, ℓs))
    block(o, (n′, n)) = reshape(view(A, (o + 1):(o + n′ * n)), n′, n)
    blocks = Vector{typeof(block(0, (0, 0)))}(undef, length(ℓs))
    o = 0
    for (i, ℓ) ∈ enumerate(ℓs)
        compute_block!(calc, ℓ, A, o)
        f(ℓ, A, o)
        blocks[i] = block(o, dims(ℓ))
        o += length(blocks[i])
    end
    blocks
end

function D_series(
    blocks::AbstractVector{<:AbstractMatrix{NT}}, ℓₘₐₓ::IT,
    m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT
) where {NT, IT<:IntegerHalf}
    series = map(enumerate(lowest_index(IT):ℓₘₐₓ)) do (i, ℓ)
        p = vec(blocks[i])
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
