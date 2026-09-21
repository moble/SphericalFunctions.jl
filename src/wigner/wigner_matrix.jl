import Base: @propagate_inbounds

"""
    AbstractWignerMatrix{IT, NT, ST}

Abstract base type for Wigner rotation‐matrix objects of a specific ``ℓ`` value.
- `IT` is the index type (an `Integer` or a `HalfOddInteger`), governing the allowed
  ranges of `m′` and `m`.
- `NT` is the number type (e.g., `ComplexF64` for D-matrices or `Float64` for d-matrices).
- `ST` is the storage type (typically `Matrix{NT}`, but other `AbstractMatrix{NT}` storage
  can be used).

The basic concrete subtypes ([`WignerMatrix`](@ref), and its aliases
[`WignerDMatrix`](@ref) and [`WignerdMatrix`](@ref)) default to storing their data in a
`Matrix{NT}` and implement `getindex` and `setindex!` so that one can use `w[m′, m]`, plus
`Matrix(w)` to materialize the block as an ordinary matrix.  The indices may be negative or
positive, and must lie in `m′ₘᵢₙ(w):m′ₘₐₓ(w)` and `mₘᵢₙ(w):mₘₐₓ(w)` respectively.

These types are *not* `AbstractMatrix`es, because half-integer indices cannot satisfy that
interface: `axes` must be integer ranges there, and `-3//2:3//2` is not one.  Consequently
linear algebra (`*`, `'`, `\\`, `lu`) does not apply to them; call `Matrix(w)` first.
Broadcasting works but returns an ordinary 1-based `Array`, dropping the natural indices.

# Methods

Methods defined for `AbstractWignerMatrix` objects include:
- `parent(w)`: the underlying data array.
- `ℓ(w)` or `ell(w)`: the value of ``ℓ``.
- `m′ₘₐₓ(w)` or `mpmax(w)`, `m′ₘᵢₙ(w)` or `mpmin(w)`: the range of ``m′``.
- `mₘₐₓ(w)` or `mmax(w)`, `mₘᵢₙ(w)` or `mmin(w)`: the range of ``m``.
- `ℓₘᵢₙ(w)` or `ellmin(w)`: the minimum value of ``ℓ``, which is either 0 or 1//2.
- `ishalfinteger(w)`: whether the indices are half-odd-integers.
- `size(w)`, `size(w, d)`: the size of the block represented, which may be smaller than
  `size(parent(w))`.
- `length(w)`: the number of elements of that block.
- `getindex(w, m′, m)`: get the value at index `(m′, m)`.
- `setindex!(w, v, m′, m)`: set the value at index `(m′, m)`.
- `axes(w)`, `axes(w, d)`: the axes of the matrix, which are 2-tuples of ranges for the `m′`
  and `m` indices.
- `Matrix(w)`, `Array(w)`, `collect(w)`: the block as an ordinary 1-based array.
- `copy(w)`, `similar(w)`, `==`, `ndims(w)`, iteration.

Note that there is deliberately *no* linear indexing (`w[i]`), for the reason given above.

# Implementation

Any new subtypes of `AbstractWignerMatrix` should inherit from this type and re-implement any of the
methods mentioned above that are not appropriate for the new type.  Specifically, the
default implementations assume that subtypes store the fields
- `parent::ST`: the underlying storage type.
- `ℓ::IT`: the value of ``ℓ``.
- `m′ₘₐₓ::IT`: the maximum value of ``m′``.

plus `m′ₘᵢₙ`, `mₘₐₓ` and `mₘᵢₙ` for the types that hold a restricted block.  For example,
if the parent array is not stored as the `parent` field, then the `parent(w)` method should
be re-implemented to return the correct parent object; likewise the default `getindex` and
`setindex!` assume that `parent(w)` is a 1-based matrix whose `[1, 1]` element is
`w[m′ₘᵢₙ(w), mₘᵢₙ(w)]`.
"""
abstract type AbstractWignerMatrix{IT<:IntegerHalf, NT, ST<:AbstractArray{NT}} end
# Note that this is deliberately *not* a subtype of `AbstractMatrix`: the natural indices
# `(m′, m)` may be half-odd-integers, which cannot satisfy the `AbstractArray` interface
# (integer `axes`).  The array-like methods that make sense are defined explicitly below.

### General methods for all AbstractWignerMatrix types

Base.parent(w::AbstractWignerMatrix) = w.parent

ℓ(w::AbstractWignerMatrix{IT}) where {IT} = w.ℓ
ℓₘᵢₙ(::IT) where {IT} = ℓₘᵢₙ(IT)
ℓₘᵢₙ(::Type{IT}) where {IT} = error("No method defined for `ℓₘᵢₙ(::Type{$IT})`.")
ℓₘᵢₙ(::Type{IT}) where {IT<:Integer} = zero(IT)
ℓₘᵢₙ(::Type{IT}) where {IT<:HalfOddInteger} = unsafe_half_odd_integer(1)
ℓₘᵢₙ(::AbstractWignerMatrix{IT}) where {IT} = ℓₘᵢₙ(IT)

m′ₘₐₓ(w::AbstractWignerMatrix{IT}) where {IT} = w.m′ₘₐₓ
m′ₘᵢₙ(w::AbstractWignerMatrix{IT}) where {IT} = w.m′ₘᵢₙ
mₘₐₓ(w::AbstractWignerMatrix{IT}) where {IT} = w.mₘₐₓ
mₘᵢₙ(w::AbstractWignerMatrix{IT}) where {IT} = w.mₘᵢₙ

const ell = ℓ
const ellmin = ℓₘᵢₙ
const mpmax = m′ₘₐₓ
const mpmin = m′ₘᵢₙ
const mmax = mₘₐₓ
const mmin = mₘᵢₙ
const smax = sₘₐₓ
const smin = sₘᵢₙ

ishalfinteger(::AbstractWignerMatrix{IT}) where {IT<:Integer} = false
ishalfinteger(::AbstractWignerMatrix{IT}) where {IT<:HalfOddInteger} = true

Base.eltype(::AbstractWignerMatrix{IT, NT, ST}) where {IT, NT, ST} = NT
Base.eltype(::Type{<:AbstractWignerMatrix{IT, NT, ST}}) where {IT, NT, ST} = NT
# Generic fallbacks refer to the underlying storage; `WignerMatrix` overrides these with the
# size of the block it represents (see below).
Base.size(w::AbstractWignerMatrix{IT, NT, ST}) where {IT, NT, ST} = size(parent(w))
Base.length(w::AbstractWignerMatrix{IT, NT, ST}) where {IT, NT, ST} = length(parent(w))

struct WignerRange{T<:IntegerHalf} <: AbstractUnitRange{T}
    start::T
    stop::T

    WignerRange(r::UnitRange{T}) where {T} = new{T}(r.start, r.stop)
end
# A `WignerRange` is indexed by position, as `Base`'s ranges are — `r[1]` is `first(r)`, and
# `firstindex` and `lastindex` below say so — so its own axes are 1-based.  An earlier version
# made `axes(r)` the range itself, after the manner of an identity-offset range, which left
# `axes` and `getindex` disagreeing: broadcasting a function over the range then read positions
# that were really values, silently wrong for an integer range and an error for a half-odd one.
# The natural, possibly half-odd, bounds are what `axes` of a *container* reports, and those
# are what `inds2string` below displays.
@inline Base.axes(r::WignerRange) = (axes1(r),)
@inline axes1(r::WignerRange) = Base.OneTo(length(r))
if VERSION < v"1.8.2"
    Base.axes1(r::WignerRange) = axes1(r)
end
Base.inds2string(inds::NTuple{2, WignerRange}) =
    string(
        "(", inds[1].start, ":", inds[1].stop, ")",
        "×",
        "(", inds[2].start, ":", inds[2].stop, ")"
    )
Base.inds2string(inds::Tuple{UnitRange, WignerRange, WignerRange}) =
    string(
        "(", inds[1].start, ":", inds[1].stop, ")",
        "×",
        "(", inds[2].start, ":", inds[2].stop, ")",
        "×",
        "(", inds[3].start, ":", inds[3].stop, ")"
    )
Base.firstindex(r::WignerRange) = 1
Base.lastindex(r::WignerRange) = length(r)
# As for `UnitRange{HalfOddInteger}` in `half_odd_integer.jl`: `Base`'s `step` for an
# `AbstractUnitRange{T}` is `oneunit(T) - zero(T)`, and its `length` independently forms
# `oneunit(zero(stop) - zero(start))`; both reach for values that `HalfOddInteger` deliberately
# lacks, so both must be given directly.  Without these two, `show` of a half-integer axis —
# `axes(w[ℓ, :])` for a half-integer `ModeWeights`, say — and `lastindex` above both throw.
# The integer case is left to `Base`.
Base.step(::WignerRange{HalfOddInteger}) = 1
Base.length(r::WignerRange{HalfOddInteger}) = max(0, (last(r) - first(r)) + 1)
function Base.getindex(v::WignerRange, i::Bool)
    throw(ArgumentError("invalid index: $i of type Bool"))
end
@propagate_inbounds function Base.getindex(v::WignerRange{T}, i::Integer) where {T}
    val = convert(T, v.start + (i - oneunit(i)))
    @boundscheck (i>0 && val <= v.stop && val >= v.start) || throw(BoundsError(v, i))
    val
end
# A `WignerRange{T}` has unit step and holds every value of type `T` between its endpoints,
# so membership of a `T` is just the bracket.  (The generic
# `in(::Real, ::AbstractRange{<:Real})` would additionally test integrality of the offset,
# which is automatic here: `x - first(r)` is an `Integer` by construction for both index
# types.)
@inline Base.in(x::T, r::WignerRange{T}) where {T<:IntegerHalf} = first(r) ≤ x ≤ last(r)
# A value of the *other* index type is never a member, and neither is anything else.
@inline Base.in(::Real, ::WignerRange) = false
# `Base` has `in(::Integer, ::AbstractUnitRange{<:Integer})`, which is neither more nor less
# specific than either method above, so without this one `1 ∈ axes(w, 1)` on an
# integer-indexed container is an ambiguity error rather than an answer.
@inline Base.in(x::Integer, r::WignerRange{<:Integer}) = first(r) ≤ x ≤ last(r)


### Bounds checking
#
# `m ∈ axes(w, 2)` is the obvious spelling of the checks below, but it builds a fresh
# `WignerRange` and then calls the generic `in`, which together cost roughly ten times the
# load they guard: 145 ns per element for a half-integer block, against 0.9 ns for the
# `OffsetMatrix` that earlier versions returned.  The helpers here test the stored limits
# directly instead, and give identical answers.

# `lo ≤ m ≤ hi`.  There is no parity test to perform: `m` is of the container's own index
# type, so a whole number cannot reach a half-integer container in the first place.  (Before
# `HalfOddInteger` existed, indices were `Rational`s and `0//1` had to be rejected here to
# keep it from being floored onto a neighbouring element.)
@inline inrange(::Type{IT}, m, lo, hi) where {IT} = lo ≤ m ≤ hi

function Base.axes(w::AbstractWignerMatrix{IT}) where {IT}
    (WignerRange(m′ₘᵢₙ(w):m′ₘₐₓ(w)), WignerRange(mₘᵢₙ(w):mₘₐₓ(w)))
end
# Trailing dimensions behave as they do for `AbstractArray`, so that generic code written
# against a plain array works unchanged here: `axes(w, d)` is `OneTo(1)`
# and `size(w, d)` is `1` for `d > ndims(w)`.
Base.axes(w::AbstractWignerMatrix, d::Integer) = d ≤ ndims(w) ? axes(w)[d] : Base.OneTo(1)
Base.size(w::AbstractWignerMatrix, d::Integer) = d ≤ ndims(w) ? size(w)[d] : 1
Base.ndims(::AbstractWignerMatrix) = 2
Base.ndims(::Type{<:AbstractWignerMatrix}) = 2


"""
    Matrix(w::AbstractWignerMatrix)

Materialize the block of the Wigner matrix represented by `w` as an ordinary `Matrix`, with
rows and columns in order of increasing `m′` and `m` (so the element `w[m′, m]` is at
`[Int(m′-m′ₘᵢₙ)+1, Int(m-mₘᵢₙ)+1]`).
"""
function Base.Matrix(w::AbstractWignerMatrix{IT, NT}) where {IT, NT}
    [w[m′, m] for m′ ∈ m′ₘᵢₙ(w):m′ₘₐₓ(w), m ∈ mₘᵢₙ(w):mₘₐₓ(w)]
end
Base.Array(w::AbstractWignerMatrix) = Matrix(w)
Base.collect(w::AbstractWignerMatrix) = Matrix(w)

function Base.:(==)(w1::AbstractWignerMatrix, w2::AbstractWignerMatrix)
    ℓ(w1) == ℓ(w2) && axes(w1) == axes(w2) && Matrix(w1) == Matrix(w2)
end

function Base.summary(io::IO, w::AbstractWignerMatrix{IT, NT}) where {IT, NT}
    print(io, Base.inds2string(axes(w)), " ", nameof(typeof(w)), "{", IT, ", ", NT, "} for ℓ=", ℓ(w))
end
Base.show(io::IO, w::AbstractWignerMatrix) = summary(io, w)
function Base.show(io::IO, ::MIME"text/plain", w::AbstractWignerMatrix)
    summary(io, w)
    println(io, ":")
    Base.print_array(io, Matrix(w))
end

@propagate_inbounds function Base.getindex(w::AbstractWignerMatrix{IT}, m′::IT, m::IT) where {IT}
    @boundscheck if !(
        inrange(IT, m′, m′ₘᵢₙ(w), m′ₘₐₓ(w)) && inrange(IT, m, mₘᵢₙ(w), mₘₐₓ(w))
    )
        throw(BoundsError(w, (m′, m)))
    end
    @inbounds Base.parent(w)[(m′-m′ₘᵢₙ(w))+1, (m-mₘᵢₙ(w))+1]
end

@propagate_inbounds function Base.setindex!(w::AbstractWignerMatrix{IT}, v, m′::IT, m::IT) where {IT}
    @boundscheck if !(
        inrange(IT, m′, m′ₘᵢₙ(w), m′ₘₐₓ(w)) && inrange(IT, m, mₘᵢₙ(w), mₘₐₓ(w))
    )
        throw(BoundsError(w, (m′, m)))
    end
    @inbounds Base.parent(w)[(m′-m′ₘᵢₙ(w))+1, (m-mₘᵢₙ(w))+1] = v
end


### Indexing with `Rational`s.
#
# A half-integer container is indexed by `HalfOddInteger`s, but `w[1//2, -3//2]` is what a
# caller naturally writes (and is what earlier versions of this package required).  These
# methods convert and re-dispatch.  There is no corresponding method for integer containers:
# `w[1//1, 0//1]` was never accepted and still is not.

@propagate_inbounds Base.getindex(w::AbstractWignerMatrix{IT}, m′::Rational, m::Rational) where
    {IT<:HalfOddInteger} = w[HalfOddInteger(m′), HalfOddInteger(m)]
@propagate_inbounds Base.setindex!(w::AbstractWignerMatrix{IT}, v, m′::Rational, m::Rational) where
    {IT<:HalfOddInteger} = (w[HalfOddInteger(m′), HalfOddInteger(m)] = v)


function validate_index_ranges(ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT) where
    {IT<:Union{Signed, HalfOddInteger}}
    # ℓₘₐₓ must be at least as big as ℓₘᵢₙ(ℓₘₐₓ)
    if ℓₘₐₓ < ℓₘᵢₙ(ℓₘₐₓ)
        error("ℓₘₐₓ=$ℓₘₐₓ must be non-negative.")
    end

    # The m′ and m ranges must be ordered correctly
    if m′ₘₐₓ < m′ₘᵢₙ
        error("m′ₘₐₓ=$m′ₘₐₓ is less than m′ₘᵢₙ=$m′ₘᵢₙ.")
    end
    if mₘₐₓ < mₘᵢₙ
        error("mₘₐₓ=$mₘₐₓ is less than mₘᵢₙ=$mₘᵢₙ.")
    end

    # The m′ and m ranges must bracket ±ℓₘᵢₙ (i.e., include 0 for integers, or ±1/2 for
    # half-integers; the recurrence needs both rows ±1/2 in the latter case).
    if m′ₘₐₓ < ℓₘᵢₙ(ℓₘₐₓ)
        error("m′ₘₐₓ=$m′ₘₐₓ is too small for this index type, $IT.")
    end
    if m′ₘᵢₙ > -ℓₘᵢₙ(ℓₘₐₓ)
        error("m′ₘᵢₙ=$m′ₘᵢₙ is too large for this index type, $IT.")
    end
    if mₘₐₓ < ℓₘᵢₙ(ℓₘₐₓ)
        error("mₘₐₓ=$mₘₐₓ is too small for this index type, $IT.")
    end
    if mₘᵢₙ > -ℓₘᵢₙ(ℓₘₐₓ)
        error("mₘᵢₙ=$mₘᵢₙ is too large for this index type, $IT.")
    end

    # The m′ and m values must be in range for ℓₘₐₓ
    if abs(m′ₘₐₓ) > ℓₘₐₓ
        error("|m′ₘₐₓ|=|$m′ₘₐₓ| is too large for ℓₘₐₓ=$ℓₘₐₓ.")
    end
    if abs(m′ₘᵢₙ) > ℓₘₐₓ
        error("|m′ₘᵢₙ|=|$m′ₘᵢₙ| is too large for ℓₘₐₓ=$ℓₘₐₓ.")
    end
    if abs(mₘₐₓ) > ℓₘₐₓ
        error("|mₘₐₓ|=|$mₘₐₓ| is too large for ℓₘₐₓ=$ℓₘₐₓ.")
    end
    if abs(mₘᵢₙ) > ℓₘₐₓ
        error("|mₘᵢₙ|=|$mₘᵢₙ| is too large for ℓₘₐₓ=$ℓₘₐₓ.")
    end

end

function validate_index_ranges(ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT) where
    {IT<:Union{Signed, HalfOddInteger}}
    # ℓₘₐₓ must be at least as big as ℓₘᵢₙ(ℓₘₐₓ)
    if ℓₘₐₓ < ℓₘᵢₙ(ℓₘₐₓ)
        error("ℓₘₐₓ=$ℓₘₐₓ must be non-negative.")
    end

    # The m′ range must be ordered correctly
    if m′ₘₐₓ < m′ₘᵢₙ
        error("m′ₘₐₓ=$m′ₘₐₓ is less than m′ₘᵢₙ=$m′ₘᵢₙ.")
    end

    # The m′ range must bracket ±ℓₘᵢₙ
    if m′ₘₐₓ < ℓₘᵢₙ(ℓₘₐₓ)
        error("m′ₘₐₓ=$m′ₘₐₓ is too small for this index type, $IT.")
    end
    if m′ₘᵢₙ > -ℓₘᵢₙ(ℓₘₐₓ)
        error("m′ₘᵢₙ=$m′ₘᵢₙ is too large for this index type, $IT.")
    end

    # The m′ values must be in range for ℓₘₐₓ
    if abs(m′ₘₐₓ) > ℓₘₐₓ
        error("|m′ₘₐₓ|=|$m′ₘₐₓ| is too large for ℓₘₐₓ=$ℓₘₐₓ.")
    end
    if abs(m′ₘᵢₙ) > ℓₘₐₓ
        error("|m′ₘᵢₙ|=|$m′ₘᵢₙ| is too large for ℓₘₐₓ=$ℓₘₐₓ.")
    end

end


@doc raw"""
    WignerMatrix{IT, NT, ST} <: AbstractWignerMatrix{IT, NT, ST}

General concrete subtype of [`AbstractWignerMatrix`](@ref) for Wigner rotation matrices,
which can include D-matrices (when `NT` is complex) or d-matrices (when `NT` is real).

In general, the storage type `ST` can be any `AbstractMatrix{NT}`, but should be 1-based.
That is, the storage should generally be either a `Matrix` or a view.  That matrix will
represent a rectangular array of values representing some or all of the Wigner matrix for a
specific ``ℓ`` value.  The first dimension corresponds to the `m′` index, and the second
dimension corresponds to the `m` index.  The allowed ranges of `m′` and `m` are governed by
the fields `m′ₘₐₓ`, `m′ₘᵢₙ`, `mₘₐₓ`, and `mₘᵢₙ`, which must satisfy
```math
\begin{aligned}
-ℓₘₐₓ &≤ m′ₘᵢₙ ≤ -ℓₘᵢₙ ≤ ℓₘᵢₙ ≤ m′ₘₐₓ ≤ ℓₘₐₓ, \\
-ℓₘₐₓ &≤ mₘᵢₙ ≤ -ℓₘᵢₙ ≤ ℓₘᵢₙ ≤ mₘₐₓ ≤ ℓₘₐₓ,
\end{aligned}
```
where `ℓₘᵢₙ` is either 0 or 1//2 depending on whether `IT` is an integer or rational type.
Both rows `±ℓₘᵢₙ` must be included because the recurrence seeds the half-integer ladder from
the pair of rows `m′ = ±1/2`; for integers this reduces to the familiar
`m′ₘᵢₙ ≤ 0 ≤ m′ₘₐₓ`.

"""
struct WignerMatrix{IT, NT, ST} <: AbstractWignerMatrix{IT, NT, ST}
    parent::ST
    ℓ::IT
    m′ₘₐₓ::IT
    m′ₘᵢₙ::IT
    mₘₐₓ::IT
    mₘᵢₙ::IT
end

# The size of the *block* represented (the parent storage may be larger)
function Base.iterate(w::WignerMatrix, state=1)
    n₁, n₂ = size(w)
    state > n₁ * n₂ && return nothing
    i, j = (state - 1) % n₁, (state - 1) ÷ n₁  # column-major over the block
    (w[w.m′ₘᵢₙ + i, w.mₘᵢₙ + j], state + 1)
end

function WignerMatrix(parent::AbstractMatrix, ℓ::Rational; kwargs...)
    WignerMatrix(parent, half_integer(ℓ); half_integer_kwargs(kwargs)...)
end
function WignerMatrix(
    parent::ST, ℓ::IT;
    mp_max::IT=ℓ, mp_min::IT=-ℓ, m_max::IT=ℓ, m_min::IT=-ℓ,
    m′ₘₐₓ::IT=mp_max, m′ₘᵢₙ::IT=mp_min, mₘₐₓ::IT=m_max, mₘᵢₙ::IT=m_min
) where {IT<:IntegerHalf, NT, ST<:AbstractMatrix{NT}}
    validate_index_ranges(ℓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    s₁, s₂ = size(parent)
    if s₁ < Int(m′ₘₐₓ - m′ₘᵢₙ + 1)
        error(
            "The extent of the first dimension in the input data must be at least "
            * "m′ₘₐₓ-m′ₘᵢₙ+1=$m′ₘₐₓ-$m′ₘᵢₙ+1=$(Int(m′ₘₐₓ - m′ₘᵢₙ + 1)); it is $s₁."
        )
    end
    if s₂ < Int(mₘₐₓ - mₘᵢₙ + 1)
        error(
            "The extent of the second dimension in the input data must be at least "
            * "mₘₐₓ-mₘᵢₙ+1=$mₘₐₓ-$mₘᵢₙ+1=$(Int(mₘₐₓ - mₘᵢₙ + 1)); it is $s₂."
        )
    end
    WignerMatrix{IT, NT, ST}(parent, ℓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
end


function Base.copy(w::WignerMatrix{IT, NT}) where {IT, NT}
    let p = copy(parent(w))
        WignerMatrix{IT, NT, typeof(p)}(p, w.ℓ, w.m′ₘₐₓ, w.m′ₘᵢₙ, w.mₘₐₓ, w.mₘᵢₙ)
    end
end

"""
    similar(w::WignerMatrix, [T=eltype(w)])

A new `WignerMatrix` with the same ℓ and the same natural `(m′, m)` axes as `w`, with
uninitialized storage of element type `T`.  The storage is a plain `Matrix` sized exactly to
the block, even when `parent(w)` is a larger array or a view.
"""
Base.similar(w::WignerMatrix) = similar(w, eltype(w))
function Base.similar(w::WignerMatrix{IT}, ::Type{T}) where {IT, T}
    let p = Matrix{T}(undef, size(w))
        WignerMatrix{IT, T, typeof(p)}(p, w.ℓ, w.m′ₘₐₓ, w.m′ₘᵢₙ, w.mₘₐₓ, w.mₘᵢₙ)
    end
end


"""
    WignerMatrixBatch{IT, NT, ST} <: AbstractWignerMatrix{IT, NT, ST}

`Nᵣ` Wigner matrices of one ``ℓ``, stored together and indexed as `w[iᵣ, m′, m]`.  This is
what [`recurrence!`](@ref) returns when `Nᵣ > 1`, for either kind of index.

The storage `parent(w)` is 1-based and 3-dimensional, ordered `[iᵣ, m′, m]`, exactly as in
the calculator.  `w[iᵣ]` gives the [`WignerMatrix`](@ref) view of one rotor's matrix, which
is then indexed naturally as `w[iᵣ][m′, m]`.

See also [`WignerMatrix`](@ref) and [`WignerSeries`](@ref).
"""
struct WignerMatrixBatch{IT, NT, ST} <: AbstractWignerMatrix{IT, NT, ST}
    parent::ST
    ℓ::IT
    m′ₘₐₓ::IT
    m′ₘᵢₙ::IT
    mₘₐₓ::IT
    mₘᵢₙ::IT
    Nᵣ::Int
end

function WignerMatrixBatch(
    parent::ST, ℓ::IT;
    m′ₘₐₓ::IT=ℓ, m′ₘᵢₙ::IT=-ℓ, mₘₐₓ::IT=ℓ, mₘᵢₙ::IT=-ℓ
) where {IT, NT, ST<:AbstractArray{NT, 3}}
    validate_index_ranges(ℓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    s₀, s₁, s₂ = size(parent)
    if s₁ < Int(m′ₘₐₓ - m′ₘᵢₙ + 1)
        error(
            "The extent of the second dimension in the input data must be at least "
            * "m′ₘₐₓ-m′ₘᵢₙ+1=$m′ₘₐₓ-$m′ₘᵢₙ+1=$(Int(m′ₘₐₓ - m′ₘᵢₙ + 1)); it is $s₁."
        )
    end
    if s₂ < Int(mₘₐₓ - mₘᵢₙ + 1)
        error(
            "The extent of the third dimension in the input data must be at least "
            * "mₘₐₓ-mₘᵢₙ+1=$mₘₐₓ-$mₘᵢₙ+1=$(Int(mₘₐₓ - mₘᵢₙ + 1)); it is $s₂."
        )
    end
    WignerMatrixBatch{IT, NT, ST}(parent, ℓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, s₀)
end

Nᵣ(w::WignerMatrixBatch) = w.Nᵣ

@propagate_inbounds function Base.getindex(
    w::WignerMatrixBatch{IT}, iᵣ::Integer, m′::IT, m::IT
) where {IT}
    @boundscheck if !(
        1 ≤ iᵣ ≤ w.Nᵣ
        && inrange(IT, m′, w.m′ₘᵢₙ, w.m′ₘₐₓ) && inrange(IT, m, w.mₘᵢₙ, w.mₘₐₓ)
    )
        throw(BoundsError(w, (iᵣ, m′, m)))
    end
    @inbounds parent(w)[iᵣ, Int(m′ - w.m′ₘᵢₙ) + 1, Int(m - w.mₘᵢₙ) + 1]
end

@propagate_inbounds function Base.setindex!(
    w::WignerMatrixBatch{IT}, v, iᵣ::Integer, m′::IT, m::IT
) where {IT}
    @boundscheck if !(
        1 ≤ iᵣ ≤ w.Nᵣ
        && inrange(IT, m′, w.m′ₘᵢₙ, w.m′ₘₐₓ) && inrange(IT, m, w.mₘᵢₙ, w.mₘₐₓ)
    )
        throw(BoundsError(w, (iᵣ, m′, m)))
    end
    @inbounds parent(w)[iᵣ, Int(m′ - w.m′ₘᵢₙ) + 1, Int(m - w.mₘᵢₙ) + 1] = v
end

# See the note on `Rational` indexing above.
@propagate_inbounds Base.getindex(w::WignerMatrixBatch{IT}, iᵣ::Integer, m′::Rational, m::Rational) where
    {IT<:HalfOddInteger} = w[iᵣ, HalfOddInteger(m′), HalfOddInteger(m)]
@propagate_inbounds Base.setindex!(w::WignerMatrixBatch{IT}, v, iᵣ::Integer, m′::Rational, m::Rational) where
    {IT<:HalfOddInteger} = (w[iᵣ, HalfOddInteger(m′), HalfOddInteger(m)] = v)

"""
    w[iᵣ]

The [`WignerMatrix`](@ref) of rotor `iᵣ` in a [`WignerMatrixBatch`](@ref), as a view; index
it naturally as `w[iᵣ][m′, m]`.
"""
@propagate_inbounds function Base.getindex(w::WignerMatrixBatch{IT, NT}, iᵣ::Integer) where {IT, NT}
    @boundscheck if !(1 ≤ iᵣ ≤ w.Nᵣ)
        throw(BoundsError(w, (iᵣ,)))
    end
    let p = view(parent(w), iᵣ, :, :)
        WignerMatrix{IT, NT, typeof(p)}(p, w.ℓ, w.m′ₘₐₓ, w.m′ₘᵢₙ, w.mₘₐₓ, w.mₘᵢₙ)
    end
end

function Base.copy(w::WignerMatrixBatch{IT, NT}) where {IT, NT}
    let p = copy(parent(w))
        WignerMatrixBatch{IT, NT, typeof(p)}(
            p, w.ℓ, w.m′ₘₐₓ, w.m′ₘᵢₙ, w.mₘₐₓ, w.mₘᵢₙ, w.Nᵣ
        )
    end
end

"""
    similar(w::WignerMatrixBatch, [T=eltype(w)])

A new `WignerMatrixBatch` with the same ℓ, `Nᵣ`, and natural axes as `w`, with uninitialized
storage of element type `T`.
"""
Base.similar(w::WignerMatrixBatch) = similar(w, eltype(w))
function Base.similar(w::WignerMatrixBatch{IT}, ::Type{T}) where {IT, T}
    let p = Array{T, 3}(undef, size(w))
        WignerMatrixBatch{IT, T, typeof(p)}(
            p, w.ℓ, w.m′ₘₐₓ, w.m′ₘᵢₙ, w.mₘₐₓ, w.mₘᵢₙ, w.Nᵣ
        )
    end
end

# Iteration in the same `[iᵣ, m′, m]` order as `Array(w)`, so that `sum`, `maximum` and the
# other reducers work on a batch exactly as they do on an ordinary 3-d array.
function Base.iterate(w::WignerMatrixBatch, state=1)
    state > length(w) && return nothing
    n₀, n₁, _ = size(w)
    iᵣ = (state - 1) % n₀
    i = ((state - 1) ÷ n₀) % n₁
    j = (state - 1) ÷ (n₀ * n₁)
    (@inbounds w[1 + iᵣ, w.m′ₘᵢₙ + i, w.mₘᵢₙ + j], state + 1)
end

Base.Array(w::WignerMatrixBatch) =
    [w[iᵣ, m′, m] for iᵣ ∈ 1:w.Nᵣ, m′ ∈ w.m′ₘᵢₙ:w.m′ₘₐₓ, m ∈ w.mₘᵢₙ:w.mₘₐₓ]
Base.collect(w::WignerMatrixBatch) = Array(w)
Base.Matrix(w::WignerMatrixBatch) =
    error("A WignerMatrixBatch is 3-dimensional; use `Array(w)`, or `Matrix(w[iᵣ])`.")

function Base.:(==)(w1::WignerMatrixBatch, w2::WignerMatrixBatch)
    ℓ(w1) == ℓ(w2) && axes(w1) == axes(w2) && Array(w1) == Array(w2)
end

function Base.show(io::IO, ::MIME"text/plain", w::WignerMatrixBatch)
    summary(io, w)
    println(io, ":")
    Base.print_array(io, Array(w))
end


"""
    DegreeBlock{IT, NT, ST}

One harmonic degree's worth of values, indexed naturally by the order ``m``: `v[m]` for
``mₘᵢₙ ≤ m ≤ mₘₐₓ``.  This is the 1-dimensional sibling of [`WignerMatrix`](@ref).

The name is deliberately neutral, because the same shape serves three things: a block of
spin-weighted harmonics at one ``ℓ``, from [`sYlm`](@ref) or [`sYlmCalculator`](@ref)'s
`ₛYₗ[s, :]`; one ``ℓ``'s worth of mode weights, from [`ModeWeights`](@ref)'s `w[ℓ, :]`; and
one spin weight's row of a [`SpinMatrix`](@ref), from `b[s, :]`.  All three are values at a
fixed degree indexed by order, whatever they mean.

The storage `parent(v)` is 1-based, in order of increasing ``m``.

See also [`DegreeBlockBatch`](@ref).
"""
struct DegreeBlock{IT, NT, ST<:AbstractVector{NT}} <: AbstractWignerMatrix{IT, NT, ST}
    parent::ST
    ℓ::IT
    mₘₐₓ::IT
    mₘᵢₙ::IT
end

function DegreeBlock(parent::AbstractVector, ℓ::Rational; kwargs...)
    DegreeBlock(parent, half_integer(ℓ); half_integer_kwargs(kwargs)...)
end
function DegreeBlock(parent::ST, ℓ::IT; mₘₐₓ::IT=ℓ, mₘᵢₙ::IT=-ℓ) where {IT<:IntegerHalf, NT, ST<:AbstractVector{NT}}
    if length(parent) < Int(mₘₐₓ - mₘᵢₙ) + 1
        error(
            "The input data must have length at least mₘₐₓ-mₘᵢₙ+1="
            * "$mₘₐₓ-$mₘᵢₙ+1=$(Int(mₘₐₓ - mₘᵢₙ) + 1); it is $(length(parent))."
        )
    end
    DegreeBlock{IT, NT, ST}(parent, ℓ, mₘₐₓ, mₘᵢₙ)
end

Base.parent(v::DegreeBlock) = v.parent
ℓ(v::DegreeBlock) = v.ℓ
ℓₘᵢₙ(::DegreeBlock{IT}) where {IT} = ℓₘᵢₙ(IT)
mₘₐₓ(v::DegreeBlock) = v.mₘₐₓ
mₘᵢₙ(v::DegreeBlock) = v.mₘᵢₙ
Base.firstindex(v::DegreeBlock) = v.mₘᵢₙ
Base.lastindex(v::DegreeBlock) = v.mₘₐₓ
Base.keys(v::DegreeBlock) = v.mₘᵢₙ:v.mₘₐₓ

@propagate_inbounds function Base.getindex(v::DegreeBlock{IT}, m::IT) where {IT}
    @boundscheck if !inrange(IT, m, v.mₘᵢₙ, v.mₘₐₓ)
        throw(BoundsError(v, m))
    end
    @inbounds parent(v)[Int(m - v.mₘᵢₙ) + 1]
end
@propagate_inbounds function Base.setindex!(v::DegreeBlock{IT}, x, m::IT) where {IT}
    @boundscheck if !inrange(IT, m, v.mₘᵢₙ, v.mₘₐₓ)
        throw(BoundsError(v, m))
    end
    @inbounds parent(v)[Int(m - v.mₘᵢₙ) + 1] = x
end

# See the note on `Rational` indexing above.
@propagate_inbounds Base.getindex(v::DegreeBlock{IT}, m::Rational) where
    {IT<:HalfOddInteger} = v[HalfOddInteger(m)]
@propagate_inbounds Base.setindex!(v::DegreeBlock{IT}, x, m::Rational) where
    {IT<:HalfOddInteger} = (v[HalfOddInteger(m)] = x)

function Base.iterate(v::DegreeBlock, state=1)
    state > length(v) && return nothing
    (@inbounds v[v.mₘᵢₙ + (state - 1)], state + 1)
end
Base.Vector(v::DegreeBlock) = [v[m] for m ∈ v.mₘᵢₙ:v.mₘₐₓ]
Base.Array(v::DegreeBlock) = Vector(v)
Base.collect(v::DegreeBlock) = Vector(v)
function Base.copy(v::DegreeBlock{IT, NT}) where {IT, NT}
    let p = copy(parent(v))
        DegreeBlock{IT, NT, typeof(p)}(p, v.ℓ, v.mₘₐₓ, v.mₘᵢₙ)
    end
end

"""
    similar(v::DegreeBlock, [T=eltype(v)])

A new `DegreeBlock` with the same ℓ and natural `m` axis as `v`, with uninitialized storage
of element type `T`.
"""
Base.similar(v::DegreeBlock) = similar(v, eltype(v))
function Base.similar(v::DegreeBlock{IT}, ::Type{T}) where {IT, T}
    let p = Vector{T}(undef, length(v))
        DegreeBlock{IT, T, typeof(p)}(p, v.ℓ, v.mₘₐₓ, v.mₘᵢₙ)
    end
end
function Base.:(==)(v1::DegreeBlock, v2::DegreeBlock)
    ℓ(v1) == ℓ(v2) && axes(v1) == axes(v2) && Vector(v1) == Vector(v2)
end
function Base.summary(io::IO, v::DegreeBlock{IT, NT}) where {IT, NT}
    print(io, "(", v.mₘᵢₙ, ":", v.mₘₐₓ, ") DegreeBlock{", IT, ", ", NT, "} for ℓ=", v.ℓ)
end
Base.show(io::IO, v::DegreeBlock) = summary(io, v)
function Base.show(io::IO, ::MIME"text/plain", v::DegreeBlock)
    summary(io, v)
    println(io, ":")
    Base.print_array(io, Vector(v))
end


"""
    DegreeBlockBatch{IT, NT, ST}

`Nᵣ` rows of values of a single ``ℓ``, stored together and indexed as `v[iᵣ, m]`.  This is
the 1-dimensional sibling of [`WignerMatrixBatch`](@ref), and is what
a batched [`sYlmCalculator`](@ref) yields for one spin weight with half-integer indices.
`v[iᵣ]` gives the [`DegreeBlock`](@ref) view of one rotor's row.

The storage `parent(v)` is 1-based and 2-dimensional, ordered `[iᵣ, m]`.
"""
struct DegreeBlockBatch{IT, NT, ST<:AbstractMatrix{NT}} <: AbstractWignerMatrix{IT, NT, ST}
    parent::ST
    ℓ::IT
    mₘₐₓ::IT
    mₘᵢₙ::IT
    Nᵣ::Int
end

function DegreeBlockBatch(
    parent::ST, ℓ::IT; mₘₐₓ::IT=ℓ, mₘᵢₙ::IT=-ℓ
) where {IT, NT, ST<:AbstractMatrix{NT}}
    s₀, s₁ = size(parent)
    if s₁ < Int(mₘₐₓ - mₘᵢₙ) + 1
        error(
            "The extent of the second dimension in the input data must be at least "
            * "mₘₐₓ-mₘᵢₙ+1=$mₘₐₓ-$mₘᵢₙ+1=$(Int(mₘₐₓ - mₘᵢₙ) + 1); it is $s₁."
        )
    end
    DegreeBlockBatch{IT, NT, ST}(parent, ℓ, mₘₐₓ, mₘᵢₙ, s₀)
end

Base.parent(v::DegreeBlockBatch) = v.parent
ℓ(v::DegreeBlockBatch) = v.ℓ
ℓₘᵢₙ(::DegreeBlockBatch{IT}) where {IT} = ℓₘᵢₙ(IT)
mₘₐₓ(v::DegreeBlockBatch) = v.mₘₐₓ
mₘᵢₙ(v::DegreeBlockBatch) = v.mₘᵢₙ
Nᵣ(v::DegreeBlockBatch) = v.Nᵣ

@propagate_inbounds function Base.getindex(
    v::DegreeBlockBatch{IT}, iᵣ::Integer, m::IT
) where {IT}
    @boundscheck if !(1 ≤ iᵣ ≤ v.Nᵣ && inrange(IT, m, v.mₘᵢₙ, v.mₘₐₓ))
        throw(BoundsError(v, (iᵣ, m)))
    end
    @inbounds parent(v)[iᵣ, Int(m - v.mₘᵢₙ) + 1]
end
@propagate_inbounds function Base.setindex!(
    v::DegreeBlockBatch{IT}, x, iᵣ::Integer, m::IT
) where {IT}
    @boundscheck if !(1 ≤ iᵣ ≤ v.Nᵣ && inrange(IT, m, v.mₘᵢₙ, v.mₘₐₓ))
        throw(BoundsError(v, (iᵣ, m)))
    end
    @inbounds parent(v)[iᵣ, Int(m - v.mₘᵢₙ) + 1] = x
end

# See the note on `Rational` indexing above.
@propagate_inbounds Base.getindex(v::DegreeBlockBatch{IT}, iᵣ::Integer, m::Rational) where
    {IT<:HalfOddInteger} = v[iᵣ, HalfOddInteger(m)]
@propagate_inbounds Base.setindex!(v::DegreeBlockBatch{IT}, x, iᵣ::Integer, m::Rational) where
    {IT<:HalfOddInteger} = (v[iᵣ, HalfOddInteger(m)] = x)

@propagate_inbounds function Base.getindex(v::DegreeBlockBatch{IT, NT}, iᵣ::Integer) where {IT, NT}
    @boundscheck if !(1 ≤ iᵣ ≤ v.Nᵣ)
        throw(BoundsError(v, (iᵣ,)))
    end
    let p = view(parent(v), iᵣ, :)
        DegreeBlock{IT, NT, typeof(p)}(p, v.ℓ, v.mₘₐₓ, v.mₘᵢₙ)
    end
end

Base.Matrix(v::DegreeBlockBatch) = [v[iᵣ, m] for iᵣ ∈ 1:v.Nᵣ, m ∈ v.mₘᵢₙ:v.mₘₐₓ]
Base.Array(v::DegreeBlockBatch) = Matrix(v)
Base.collect(v::DegreeBlockBatch) = Matrix(v)
function Base.copy(v::DegreeBlockBatch{IT, NT}) where {IT, NT}
    let p = copy(parent(v))
        DegreeBlockBatch{IT, NT, typeof(p)}(p, v.ℓ, v.mₘₐₓ, v.mₘᵢₙ, v.Nᵣ)
    end
end

"""
    similar(v::DegreeBlockBatch, [T=eltype(v)])

A new `DegreeBlockBatch` with the same ℓ, `Nᵣ`, and natural `m` axis as `v`, with
uninitialized storage of element type `T`.
"""
Base.similar(v::DegreeBlockBatch) = similar(v, eltype(v))
function Base.similar(v::DegreeBlockBatch{IT}, ::Type{T}) where {IT, T}
    let p = Matrix{T}(undef, size(v))
        DegreeBlockBatch{IT, T, typeof(p)}(p, v.ℓ, v.mₘₐₓ, v.mₘᵢₙ, v.Nᵣ)
    end
end

# Iteration in the same `[iᵣ, m]` order as `Matrix(v)`.
function Base.iterate(v::DegreeBlockBatch, state=1)
    state > length(v) && return nothing
    n₀, _ = size(v)
    iᵣ = (state - 1) % n₀
    i = (state - 1) ÷ n₀
    (@inbounds v[1 + iᵣ, v.mₘᵢₙ + i], state + 1)
end
function Base.:(==)(v1::DegreeBlockBatch, v2::DegreeBlockBatch)
    ℓ(v1) == ℓ(v2) && axes(v1) == axes(v2) && Matrix(v1) == Matrix(v2)
end
function Base.summary(io::IO, v::DegreeBlockBatch{IT, NT}) where {IT, NT}
    print(
        io, "(1:", v.Nᵣ, ")×(", v.mₘᵢₙ, ":", v.mₘₐₓ, ") ",
        "DegreeBlockBatch{", IT, ", ", NT, "} for ℓ=", v.ℓ
    )
end
Base.show(io::IO, v::DegreeBlockBatch) = summary(io, v)
function Base.show(io::IO, ::MIME"text/plain", v::DegreeBlockBatch)
    summary(io, v)
    println(io, ":")
    Base.print_array(io, Matrix(v))
end


"""
    SpinMatrix{IT, NT, ST}

The values of a single ``ℓ`` for a range of spin weights, indexed naturally by ``(s, m)``:
`b[s, m]` for ``sₘᵢₙ ≤ s ≤ sₘₐₓ`` and ``mₘᵢₙ ≤ m ≤ mₘₐₓ``, and `b[s, :]` for one whole row as
a [`DegreeBlock`](@ref).  This is what an [`sYlmCalculator`](@ref) built for a range of spin
weights yields for each ``ℓ``.

The storage `parent(b)` is 1-based and 2-dimensional, ordered `[s, m]`.

The spin axis is under none of the restrictions a [`WignerMatrix`](@ref) places on its ``m′``:
it is whatever range of spin weights was asked for, so it may lie wholly on one side of zero,
and it may reach beyond ``ℓ`` — the harmonics with ``|s| > ℓ`` simply vanish, and a calculator
stores them as zeros.

See also [`SpinMatrixBatch`](@ref) and [`DegreeBlock`](@ref).
"""
struct SpinMatrix{IT, NT, ST<:AbstractMatrix{NT}} <: AbstractWignerMatrix{IT, NT, ST}
    parent::ST
    ℓ::IT
    sₘₐₓ::IT
    sₘᵢₙ::IT
    mₘₐₓ::IT
    mₘᵢₙ::IT
end

function SpinMatrix(parent::AbstractMatrix, ℓ::Rational; kwargs...)
    SpinMatrix(parent, half_integer(ℓ); half_integer_kwargs(kwargs)...)
end
function SpinMatrix(
    parent::ST, ℓ::IT; sₘₐₓ::IT, sₘᵢₙ::IT, mₘₐₓ::IT=ℓ, mₘᵢₙ::IT=-ℓ
) where {IT<:IntegerHalf, NT, ST<:AbstractMatrix{NT}}
    s₁, s₂ = size(parent)
    if s₁ < Int(sₘₐₓ - sₘᵢₙ) + 1
        error(
            "The extent of the first dimension in the input data must be at least "
            * "sₘₐₓ-sₘᵢₙ+1=$sₘₐₓ-$sₘᵢₙ+1=$(Int(sₘₐₓ - sₘᵢₙ) + 1); it is $s₁."
        )
    end
    if s₂ < Int(mₘₐₓ - mₘᵢₙ) + 1
        error(
            "The extent of the second dimension in the input data must be at least "
            * "mₘₐₓ-mₘᵢₙ+1=$mₘₐₓ-$mₘᵢₙ+1=$(Int(mₘₐₓ - mₘᵢₙ) + 1); it is $s₂."
        )
    end
    SpinMatrix{IT, NT, ST}(parent, ℓ, sₘₐₓ, sₘᵢₙ, mₘₐₓ, mₘᵢₙ)
end

Base.parent(b::SpinMatrix) = b.parent
ℓ(b::SpinMatrix) = b.ℓ
ℓₘᵢₙ(::SpinMatrix{IT}) where {IT} = ℓₘᵢₙ(IT)
sₘₐₓ(b::SpinMatrix) = b.sₘₐₓ
sₘᵢₙ(b::SpinMatrix) = b.sₘᵢₙ
mₘₐₓ(b::SpinMatrix) = b.mₘₐₓ
mₘᵢₙ(b::SpinMatrix) = b.mₘᵢₙ
ishalfinteger(::SpinMatrix{IT}) where {IT<:Integer} = false
ishalfinteger(::SpinMatrix{IT}) where {IT<:HalfOddInteger} = true
Base.keys(b::SpinMatrix) = b.sₘᵢₙ:b.sₘₐₓ

@propagate_inbounds function Base.getindex(b::SpinMatrix{IT}, s::IT, m::IT) where {IT}
    @boundscheck if !(inrange(IT, s, b.sₘᵢₙ, b.sₘₐₓ) && inrange(IT, m, b.mₘᵢₙ, b.mₘₐₓ))
        throw(BoundsError(b, (s, m)))
    end
    @inbounds parent(b)[Int(s - b.sₘᵢₙ) + 1, Int(m - b.mₘᵢₙ) + 1]
end
@propagate_inbounds function Base.setindex!(b::SpinMatrix{IT}, x, s::IT, m::IT) where {IT}
    @boundscheck if !(inrange(IT, s, b.sₘᵢₙ, b.sₘₐₓ) && inrange(IT, m, b.mₘᵢₙ, b.mₘₐₓ))
        throw(BoundsError(b, (s, m)))
    end
    @inbounds parent(b)[Int(s - b.sₘᵢₙ) + 1, Int(m - b.mₘᵢₙ) + 1] = x
end

# See the note on `Rational` indexing above.
@propagate_inbounds Base.getindex(b::SpinMatrix{IT}, s::Rational, m::Rational) where
    {IT<:HalfOddInteger} = b[HalfOddInteger(s), HalfOddInteger(m)]
@propagate_inbounds Base.setindex!(b::SpinMatrix{IT}, x, s::Rational, m::Rational) where
    {IT<:HalfOddInteger} = (b[HalfOddInteger(s), HalfOddInteger(m)] = x)

"""
    b[s, :]

The row of one spin weight of a [`SpinMatrix`](@ref), as a [`DegreeBlock`](@ref) view; it is
then indexed naturally as `b[s, :][m]`.  The spelling follows [`ModeWeights`](@ref)'s
`w[ℓ, :]`, so that a loop over spin weights reads the same in both places.
"""
@propagate_inbounds function Base.getindex(b::SpinMatrix{IT, NT}, s::IT, ::Colon) where {IT, NT}
    @boundscheck if !inrange(IT, s, b.sₘᵢₙ, b.sₘₐₓ)
        throw(BoundsError(b, (s, :)))
    end
    let p = view(parent(b), Int(s - b.sₘᵢₙ) + 1, :)
        DegreeBlock{IT, NT, typeof(p)}(p, b.ℓ, b.mₘₐₓ, b.mₘᵢₙ)
    end
end
@propagate_inbounds Base.getindex(b::SpinMatrix{IT}, s::Rational, ::Colon) where
    {IT<:HalfOddInteger} = b[HalfOddInteger(s), :]

Base.Matrix(b::SpinMatrix) = [b[s, m] for s ∈ b.sₘᵢₙ:b.sₘₐₓ, m ∈ b.mₘᵢₙ:b.mₘₐₓ]
Base.Array(b::SpinMatrix) = Matrix(b)
Base.collect(b::SpinMatrix) = Matrix(b)
function Base.copy(b::SpinMatrix{IT, NT}) where {IT, NT}
    let p = copy(parent(b))
        SpinMatrix{IT, NT, typeof(p)}(p, b.ℓ, b.sₘₐₓ, b.sₘᵢₙ, b.mₘₐₓ, b.mₘᵢₙ)
    end
end

"""
    similar(b::SpinMatrix, [T=eltype(b)])

A new `SpinMatrix` with the same ℓ and the same natural `(s, m)` axes as `b`, with
uninitialized storage of element type `T`.  The storage is a plain `Matrix` sized exactly to
the block, even when `parent(b)` is a larger array or a view.
"""
Base.similar(b::SpinMatrix) = similar(b, eltype(b))
function Base.similar(b::SpinMatrix{IT}, ::Type{T}) where {IT, T}
    let p = Matrix{T}(undef, size(b))
        SpinMatrix{IT, T, typeof(p)}(p, b.ℓ, b.sₘₐₓ, b.sₘᵢₙ, b.mₘₐₓ, b.mₘᵢₙ)
    end
end

# Iteration in the same `[s, m]` column-major order as `Matrix(b)`.
function Base.iterate(b::SpinMatrix, state=1)
    n₁, n₂ = size(b)
    state > n₁ * n₂ && return nothing
    i, j = (state - 1) % n₁, (state - 1) ÷ n₁
    (@inbounds b[b.sₘᵢₙ + i, b.mₘᵢₙ + j], state + 1)
end
function Base.:(==)(b1::SpinMatrix, b2::SpinMatrix)
    ℓ(b1) == ℓ(b2) && axes(b1) == axes(b2) && Matrix(b1) == Matrix(b2)
end
function Base.summary(io::IO, b::SpinMatrix{IT, NT}) where {IT, NT}
    print(
        io, "(", b.sₘᵢₙ, ":", b.sₘₐₓ, ")×(", b.mₘᵢₙ, ":", b.mₘₐₓ, ") ",
        "SpinMatrix{", IT, ", ", NT, "} for ℓ=", b.ℓ
    )
end
Base.show(io::IO, b::SpinMatrix) = summary(io, b)
function Base.show(io::IO, ::MIME"text/plain", b::SpinMatrix)
    summary(io, b)
    println(io, ":")
    Base.print_array(io, Matrix(b))
end


"""
    SpinMatrixBatch{IT, NT, ST}

`Nᵣ` [`SpinMatrix`](@ref) blocks of one ``ℓ``, stored together and indexed as `b[iᵣ, s, m]`.
This is what a batched [`sYlmCalculator`](@ref) yields when it was built for a range
of spin weights.  `b[iᵣ]` gives the `SpinMatrix` view of one rotor's block, which is then indexed
naturally as `b[iᵣ][s, m]`.

The storage `parent(b)` is 1-based and 3-dimensional, ordered `[iᵣ, s, m]`, exactly as in the
calculator.

See also [`SpinMatrix`](@ref) and [`DegreeBlockBatch`](@ref).
"""
struct SpinMatrixBatch{IT, NT, ST<:AbstractArray{NT, 3}} <: AbstractWignerMatrix{IT, NT, ST}
    parent::ST
    ℓ::IT
    sₘₐₓ::IT
    sₘᵢₙ::IT
    mₘₐₓ::IT
    mₘᵢₙ::IT
    Nᵣ::Int
end

function SpinMatrixBatch(parent::AbstractArray{<:Any, 3}, ℓ::Rational; kwargs...)
    SpinMatrixBatch(parent, half_integer(ℓ); half_integer_kwargs(kwargs)...)
end
function SpinMatrixBatch(
    parent::ST, ℓ::IT; sₘₐₓ::IT, sₘᵢₙ::IT, mₘₐₓ::IT=ℓ, mₘᵢₙ::IT=-ℓ
) where {IT<:IntegerHalf, NT, ST<:AbstractArray{NT, 3}}
    s₀, s₁, s₂ = size(parent)
    if s₁ < Int(sₘₐₓ - sₘᵢₙ) + 1
        error(
            "The extent of the second dimension in the input data must be at least "
            * "sₘₐₓ-sₘᵢₙ+1=$sₘₐₓ-$sₘᵢₙ+1=$(Int(sₘₐₓ - sₘᵢₙ) + 1); it is $s₁."
        )
    end
    if s₂ < Int(mₘₐₓ - mₘᵢₙ) + 1
        error(
            "The extent of the third dimension in the input data must be at least "
            * "mₘₐₓ-mₘᵢₙ+1=$mₘₐₓ-$mₘᵢₙ+1=$(Int(mₘₐₓ - mₘᵢₙ) + 1); it is $s₂."
        )
    end
    SpinMatrixBatch{IT, NT, ST}(parent, ℓ, sₘₐₓ, sₘᵢₙ, mₘₐₓ, mₘᵢₙ, s₀)
end

Base.parent(b::SpinMatrixBatch) = b.parent
ℓ(b::SpinMatrixBatch) = b.ℓ
ℓₘᵢₙ(::SpinMatrixBatch{IT}) where {IT} = ℓₘᵢₙ(IT)
sₘₐₓ(b::SpinMatrixBatch) = b.sₘₐₓ
sₘᵢₙ(b::SpinMatrixBatch) = b.sₘᵢₙ
mₘₐₓ(b::SpinMatrixBatch) = b.mₘₐₓ
mₘᵢₙ(b::SpinMatrixBatch) = b.mₘᵢₙ
Nᵣ(b::SpinMatrixBatch) = b.Nᵣ
ishalfinteger(::SpinMatrixBatch{IT}) where {IT<:Integer} = false
ishalfinteger(::SpinMatrixBatch{IT}) where {IT<:HalfOddInteger} = true

@propagate_inbounds function Base.getindex(
    b::SpinMatrixBatch{IT}, iᵣ::Integer, s::IT, m::IT
) where {IT}
    @boundscheck if !(
        1 ≤ iᵣ ≤ b.Nᵣ
        && inrange(IT, s, b.sₘᵢₙ, b.sₘₐₓ) && inrange(IT, m, b.mₘᵢₙ, b.mₘₐₓ)
    )
        throw(BoundsError(b, (iᵣ, s, m)))
    end
    @inbounds parent(b)[iᵣ, Int(s - b.sₘᵢₙ) + 1, Int(m - b.mₘᵢₙ) + 1]
end
@propagate_inbounds function Base.setindex!(
    b::SpinMatrixBatch{IT}, x, iᵣ::Integer, s::IT, m::IT
) where {IT}
    @boundscheck if !(
        1 ≤ iᵣ ≤ b.Nᵣ
        && inrange(IT, s, b.sₘᵢₙ, b.sₘₐₓ) && inrange(IT, m, b.mₘᵢₙ, b.mₘₐₓ)
    )
        throw(BoundsError(b, (iᵣ, s, m)))
    end
    @inbounds parent(b)[iᵣ, Int(s - b.sₘᵢₙ) + 1, Int(m - b.mₘᵢₙ) + 1] = x
end

# See the note on `Rational` indexing above.
@propagate_inbounds Base.getindex(
    b::SpinMatrixBatch{IT}, iᵣ::Integer, s::Rational, m::Rational
) where {IT<:HalfOddInteger} = b[iᵣ, HalfOddInteger(s), HalfOddInteger(m)]
@propagate_inbounds Base.setindex!(
    b::SpinMatrixBatch{IT}, x, iᵣ::Integer, s::Rational, m::Rational
) where {IT<:HalfOddInteger} = (b[iᵣ, HalfOddInteger(s), HalfOddInteger(m)] = x)

"""
    b[:, s, :]

The rows of one spin weight of a [`SpinMatrixBatch`](@ref), over every rotor, as a
[`DegreeBlockBatch`](@ref) view; it is then indexed naturally as `b[:, s, :][iᵣ, m]`.  The
spelling follows [`SpinMatrix`](@ref)'s `b[s, :]`.
"""
@propagate_inbounds function Base.getindex(
    b::SpinMatrixBatch{IT, NT}, ::Colon, s::IT, ::Colon
) where {IT, NT}
    @boundscheck if !inrange(IT, s, b.sₘᵢₙ, b.sₘₐₓ)
        throw(BoundsError(b, (:, s, :)))
    end
    let p = view(parent(b), :, Int(s - b.sₘᵢₙ) + 1, :)
        DegreeBlockBatch{IT, NT, typeof(p)}(p, b.ℓ, b.mₘₐₓ, b.mₘᵢₙ, b.Nᵣ)
    end
end
@propagate_inbounds Base.getindex(b::SpinMatrixBatch{IT}, ::Colon, s::Rational, ::Colon) where
    {IT<:HalfOddInteger} = b[:, HalfOddInteger(s), :]

"""
    b[iᵣ]

The [`SpinMatrix`](@ref) of rotor `iᵣ` in a [`SpinMatrixBatch`](@ref), as a view; index it
naturally as `b[iᵣ][s, m]`.
"""
@propagate_inbounds function Base.getindex(b::SpinMatrixBatch{IT, NT}, iᵣ::Integer) where {IT, NT}
    @boundscheck if !(1 ≤ iᵣ ≤ b.Nᵣ)
        throw(BoundsError(b, (iᵣ,)))
    end
    let p = view(parent(b), iᵣ, :, :)
        SpinMatrix{IT, NT, typeof(p)}(p, b.ℓ, b.sₘₐₓ, b.sₘᵢₙ, b.mₘₐₓ, b.mₘᵢₙ)
    end
end

function Base.Array(b::SpinMatrixBatch)
    [b[iᵣ, s, m] for iᵣ ∈ 1:b.Nᵣ, s ∈ b.sₘᵢₙ:b.sₘₐₓ, m ∈ b.mₘᵢₙ:b.mₘₐₓ]
end
Base.collect(b::SpinMatrixBatch) = Array(b)
function Base.copy(b::SpinMatrixBatch{IT, NT}) where {IT, NT}
    let p = copy(parent(b))
        SpinMatrixBatch{IT, NT, typeof(p)}(p, b.ℓ, b.sₘₐₓ, b.sₘᵢₙ, b.mₘₐₓ, b.mₘᵢₙ, b.Nᵣ)
    end
end

"""
    similar(b::SpinMatrixBatch, [T=eltype(b)])

A new `SpinMatrixBatch` with the same ℓ, `Nᵣ`, and natural `(s, m)` axes as `b`, with
uninitialized storage of element type `T`.
"""
Base.similar(b::SpinMatrixBatch) = similar(b, eltype(b))
function Base.similar(b::SpinMatrixBatch{IT}, ::Type{T}) where {IT, T}
    let p = Array{T, 3}(undef, size(b))
        SpinMatrixBatch{IT, T, typeof(p)}(p, b.ℓ, b.sₘₐₓ, b.sₘᵢₙ, b.mₘₐₓ, b.mₘᵢₙ, b.Nᵣ)
    end
end

# Iteration in the same `[iᵣ, s, m]` order as `Array(b)`.
function Base.iterate(b::SpinMatrixBatch, state=1)
    state > length(b) && return nothing
    n₀, n₁, _ = size(b)
    iᵣ = (state - 1) % n₀
    i = ((state - 1) ÷ n₀) % n₁
    j = (state - 1) ÷ (n₀ * n₁)
    (@inbounds b[1 + iᵣ, b.sₘᵢₙ + i, b.mₘᵢₙ + j], state + 1)
end
function Base.:(==)(b1::SpinMatrixBatch, b2::SpinMatrixBatch)
    ℓ(b1) == ℓ(b2) && axes(b1) == axes(b2) && Array(b1) == Array(b2)
end
function Base.summary(io::IO, b::SpinMatrixBatch{IT, NT}) where {IT, NT}
    print(
        io, "(1:", b.Nᵣ, ")×(", b.sₘᵢₙ, ":", b.sₘₐₓ, ")×(", b.mₘᵢₙ, ":", b.mₘₐₓ, ") ",
        "SpinMatrixBatch{", IT, ", ", NT, "} for ℓ=", b.ℓ
    )
end
Base.show(io::IO, b::SpinMatrixBatch) = summary(io, b)
function Base.show(io::IO, ::MIME"text/plain", b::SpinMatrixBatch)
    summary(io, b)
    println(io, ":")
    Base.print_array(io, Array(b))
end


"""
    WignerSeries{IT, VT}

The blocks of a Wigner matrix for every ``ℓ`` from `ℓₘᵢₙ` to `ℓₘₐₓ`, indexed by ``ℓ``:
`s[ℓ]` is the block of order `ℓ`, and `s[ℓ][m′, m]` an element of it.

This is what [`D`](@ref) and [`d`](@ref) return, for either kind of index.  It is iterable
and has `length`, `first`, and `last`.

See also [`WignerMatrix`](@ref) and [`WignerMatrixBatch`](@ref).
"""
struct WignerSeries{IT, VT<:AbstractVector}
    blocks::VT  # 1-based; blocks[i] is the block for ℓ = ℓₘᵢₙ + (i-1)
    ℓₘᵢₙ::IT
    ℓₘₐₓ::IT
    function WignerSeries(blocks::VT, ℓₘᵢₙ::IT, ℓₘₐₓ::IT) where {IT, VT<:AbstractVector}
        if length(blocks) != Int(ℓₘₐₓ - ℓₘᵢₙ) + 1
            error(
                "Got $(length(blocks)) blocks, but ℓ ∈ $ℓₘᵢₙ:$ℓₘₐₓ needs "
                * "$(Int(ℓₘₐₓ - ℓₘᵢₙ) + 1)."
            )
        end
        new{IT, VT}(blocks, ℓₘᵢₙ, ℓₘₐₓ)
    end
end

ℓₘᵢₙ(s::WignerSeries) = s.ℓₘᵢₙ
ℓₘₐₓ(s::WignerSeries) = s.ℓₘₐₓ
Base.parent(s::WignerSeries) = s.blocks
Base.length(s::WignerSeries) = length(s.blocks)
Base.eltype(s::WignerSeries) = eltype(s.blocks)
Base.eltype(::Type{<:WignerSeries{IT, VT}}) where {IT, VT} = eltype(VT)
Base.axes(s::WignerSeries) = (WignerRange(s.ℓₘᵢₙ:s.ℓₘₐₓ),)
Base.axes(s::WignerSeries, d::Integer) = d ≤ 1 ? axes(s)[d] : Base.OneTo(1)
Base.ndims(::WignerSeries) = 1
Base.ndims(::Type{<:WignerSeries}) = 1
Base.size(s::WignerSeries) = (length(s),)
Base.size(s::WignerSeries, d::Integer) = d ≤ 1 ? length(s) : 1
Base.keys(s::WignerSeries) = s.ℓₘᵢₙ:s.ℓₘₐₓ
Base.iterate(s::WignerSeries, state=1) = iterate(s.blocks, state)
Base.firstindex(s::WignerSeries) = s.ℓₘᵢₙ
Base.lastindex(s::WignerSeries) = s.ℓₘₐₓ

@propagate_inbounds function Base.getindex(s::WignerSeries{IT}, ℓ) where {IT}
    # Deliberately *not* inside `@boundscheck`, and deliberately before the `convert`: a
    # whole number asked of a half-integer series should be told what this series holds,
    # rather than shown the bare `InexactError` that `convert` would raise.
    if !isindex(IT, ℓ)
        throw(ArgumentError(
            "ℓ=$ℓ is not one of the ℓ values of this series, which runs over "
            * "ℓ ∈ $(s.ℓₘᵢₙ):$(s.ℓₘₐₓ) in steps of 1."
        ))
    end
    let ℓ = convert(IT, ℓ)
        @boundscheck if ℓ < s.ℓₘᵢₙ || ℓ > s.ℓₘₐₓ
            throw(BoundsError(s, ℓ))
        end
        @inbounds s.blocks[Int(ℓ - s.ℓₘᵢₙ) + 1]
    end
end

Base.copy(s::WignerSeries) = WignerSeries(map(copy, s.blocks), s.ℓₘᵢₙ, s.ℓₘₐₓ)
Base.similar(s::WignerSeries) = WignerSeries(map(similar, s.blocks), s.ℓₘᵢₙ, s.ℓₘₐₓ)
function Base.:(==)(s1::WignerSeries, s2::WignerSeries)
    ℓₘᵢₙ(s1) == ℓₘᵢₙ(s2) && ℓₘₐₓ(s1) == ℓₘₐₓ(s2) &&
        all(b1 == b2 for (b1, b2) ∈ zip(s1.blocks, s2.blocks))
end

function Base.show(io::IO, s::WignerSeries{IT}) where {IT}
    print(io, "WignerSeries{$IT} for ℓ ∈ $(s.ℓₘᵢₙ):$(s.ℓₘₐₓ)")
end
function Base.show(io::IO, ::MIME"text/plain", s::WignerSeries)
    show(io, s)
    println(io, ":")
    for ℓ ∈ s.ℓₘᵢₙ:s.ℓₘₐₓ
        println(io, " ℓ = ", ℓ, ":")
        show(io, MIME("text/plain"), s[ℓ])
        println(io)
    end
end


"""
    WignerDMatrix{IT, RT, ST}
    WignerDMatrix(parent, ℓ; m′ₘₐₓ=ℓ, m′ₘᵢₙ=-ℓ, mₘₐₓ=ℓ, mₘᵢₙ=-ℓ)
    WignerDMatrix(Complex{RT}, ℓ, m′ₘₐₓ=ℓ)

One block of Wigner's ``𝔇^{(ℓ)}_{m′,m}`` matrix: a [`WignerMatrix`](@ref) whose number type
is complex.  This is the type of each element of the [`WignerSeries`](@ref) that [`D`](@ref)
returns.

The first constructor wraps an existing `AbstractMatrix{Complex{RT}}`, which must be at
least `(m′ₘₐₓ-m′ₘᵢₙ+1) × (mₘₐₓ-mₘᵢₙ+1)`; the second allocates uninitialized storage for the
whole block.  `ℓ` may be an `Integer` or a half-integer `Rational` (denominator 2), and the
limits must then match it.

# Example
```julia
w = WignerDMatrix(ComplexF64, 2)   # uninitialized 5×5 block, m′, m ∈ -2:2
w[1, -2] = 3.0 + 0im
w[1, -2]                           # 3.0 + 0.0im
```

This is a type *alias*, not a distinct type, so `show` and `summary` name the underlying
`WignerMatrix`.  See also [`WignerdMatrix`](@ref) for the real ``d`` matrices.
"""
const WignerDMatrix{IT, RT, ST} = WignerMatrix{IT, Complex{RT}, ST} where {IT, RT<:Real, ST<:AbstractMatrix{Complex{RT}}}

"""
    WignerdMatrix{IT, RT, ST}
    WignerdMatrix(parent, ℓ; m′ₘₐₓ=ℓ, m′ₘᵢₙ=-ℓ, mₘₐₓ=ℓ, mₘᵢₙ=-ℓ)
    WignerdMatrix(RT, ℓ, m′ₘₐₓ=ℓ)

One block of Wigner's real ``d^{(ℓ)}_{m′,m}(β)`` matrix: a [`WignerMatrix`](@ref) whose
number type is real.  This is the type of each element of the [`WignerSeries`](@ref) that
[`d`](@ref) returns for half-integer `ℓₘₐₓ`.  It is the real analogue of
[`WignerDMatrix`](@ref) in every other respect, including being an alias rather than a
distinct type.

# Example
```julia
w = WignerdMatrix(Float64, 3//2)   # uninitialized 4×4 block, m′, m ∈ -3//2:3//2
w[1//2, -1//2] = 0.25
```
"""
const WignerdMatrix{IT, RT, ST} = WignerMatrix{IT, RT, ST} where {IT, RT<:Real, ST<:AbstractMatrix{RT}}

# Constructors for WignerDMatrix (complex)
function WignerDMatrix(parent::AbstractMatrix{Complex{RT}}, ℓ::Rational; kwargs...) where {RT<:Real}
    WignerDMatrix(parent, half_integer(ℓ); half_integer_kwargs(kwargs)...)
end
function WignerDMatrix(parent::ST, ℓ::IT; kwargs...) where {IT<:IntegerHalf, RT<:Real, ST<:AbstractMatrix{Complex{RT}}}
    WignerMatrix(parent, ℓ; kwargs...)
end
function WignerDMatrix(parent::ST, ℓ::IT; kwargs...) where {IT, RT<:Real, ST<:AbstractMatrix{RT}}
    error(
        "WignerDMatrix only supports complex types; the input type is $RT.\n"
        * "Perhaps you meant to use WignerdMatrix?\n"
    )
end
function WignerDMatrix(::Type{Complex{RT}}, ℓ::Rational, m′ₘₐₓ=ℓ; kwargs...) where {RT<:Real}
    WignerDMatrix(Complex{RT}, half_integer(ℓ), half_integer(m′ₘₐₓ); half_integer_kwargs(kwargs)...)
end
function WignerDMatrix(::Type{Complex{RT}}, ℓ::IT, m′ₘₐₓ::IT=ℓ; kwargs...) where {RT<:Real, IT<:IntegerHalf}
    # Validate before sizing the storage: `Int(2ℓ + 1)` on, say, ℓ = 5//3 would otherwise
    # throw a bare `InexactError` instead of the message explaining the denominator rule.
    validate_index_ranges(ℓ, m′ₘₐₓ, -m′ₘₐₓ)
    parent = Matrix{Complex{RT}}(undef, Int(m′ₘₐₓ - (-m′ₘₐₓ) + 1), Int(2ℓ + 1))
    WignerMatrix(parent, ℓ; m′ₘₐₓ=m′ₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ, kwargs...)
end

# Constructors for WignerdMatrix (real)
function WignerdMatrix(parent::AbstractMatrix{RT}, ℓ::Rational; kwargs...) where {RT<:Real}
    WignerdMatrix(parent, half_integer(ℓ); half_integer_kwargs(kwargs)...)
end
function WignerdMatrix(parent::ST, ℓ::IT; kwargs...) where {IT<:IntegerHalf, RT<:Real, ST<:AbstractMatrix{RT}}
    WignerMatrix(parent, ℓ; kwargs...)
end
function WignerdMatrix(parent::ST, ℓ::IT; kwargs...) where {IT, RT<:Real, ST<:AbstractMatrix{Complex{RT}}}
    error(
        "WignerdMatrix only supports real types; the input type is Complex{$RT}.\n"
        * "Perhaps you meant to use WignerDMatrix?"
    )
end
function WignerdMatrix(::Type{RT}, ℓ::Rational, m′ₘₐₓ=ℓ; kwargs...) where {RT<:Real}
    WignerdMatrix(RT, half_integer(ℓ), half_integer(m′ₘₐₓ); half_integer_kwargs(kwargs)...)
end
function WignerdMatrix(::Type{RT}, ℓ::IT, m′ₘₐₓ::IT=ℓ; kwargs...) where {RT<:Real, IT<:IntegerHalf}
    validate_index_ranges(ℓ, m′ₘₐₓ, -m′ₘₐₓ)
    parent = Matrix{RT}(undef, Int(m′ₘₐₓ - (-m′ₘₐₓ) + 1), Int(2ℓ + 1))
    WignerMatrix(parent, ℓ; m′ₘₐₓ=m′ₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ, kwargs...)
end


@testitem "WignerMatrix" begin
    import SphericalFunctions: WignerDMatrix, WignerdMatrix,
        parent, ell, mpmax, mpmin, mmax, mmin, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, ℓₘᵢₙ,
        HalfOddInteger, half_integer

    # Check that mixed-up types throw an error
    @test_throws "WignerDMatrix only supports complex types" WignerDMatrix(rand(Float64, 3, 3), 1)
    @test_throws "WignerdMatrix only supports real types" WignerdMatrix(rand(ComplexF64, 3, 3), 1)
    @test_throws "WignerDMatrix only supports complex types" WignerDMatrix(rand(Float64, 2, 2), 1//2)
    @test_throws "WignerdMatrix only supports real types" WignerdMatrix(rand(ComplexF64, 2, 2), 1//2)

    # Check that a negative ℓ value throws an error
    @test_throws "ℓₘₐₓ=-1 must be non-negative." WignerDMatrix(rand(ComplexF64, 3, 3), -1)
    @test_throws "ℓₘₐₓ=-1 must be non-negative." WignerdMatrix(rand(Float64, 3, 3), -1)
    @test_throws "ℓₘₐₓ=-1//2 must be non-negative." WignerDMatrix(rand(ComplexF64, 2, 2), -1//2)
    @test_throws "ℓₘₐₓ=-1//2 must be non-negative." WignerdMatrix(rand(Float64, 2, 2), -1//2)

    # Check that a `Rational` ℓ which is not a half-odd-integer throws an error.  `1//3` has
    # the wrong denominator; `2//2` normalizes to the whole number `1`, which is not a
    # half-odd-integer at all, so it is an `InexactError` rather than a denominator complaint.
    @test_throws "must have denominator 2" WignerDMatrix(rand(ComplexF64, 3, 3), 1//3)
    @test_throws "must have denominator 2" WignerdMatrix(rand(Float64, 3, 3), 1//3)
    @test_throws "must have denominator 2" WignerDMatrix(rand(ComplexF64, 2, 2), 1//3)
    @test_throws "must have denominator 2" WignerdMatrix(rand(Float64, 2, 2), 1//3)
    @test_throws "must have denominator 2" WignerDMatrix(rand(ComplexF64, 3, 3), 2//2)
    @test_throws "must have denominator 2" WignerdMatrix(rand(Float64, 3, 3), 2//2)
    @test_throws "must have denominator 2" WignerDMatrix(rand(ComplexF64, 2, 2), 2//2)
    @test_throws "must have denominator 2" WignerdMatrix(rand(Float64, 2, 2), 2//2)

    #for ℓ ∈ Any[collect(0:8); collect(1//2:15//2)]
    ℓₘₐₓ = 2
    # Encode on twice-indices, so that the arithmetic is `Int` for both index types (a
    # `HalfOddInteger` may only be multiplied by an even integer).  `2x + 6` is in `1:11`
    # for every index used here, so base 25 keeps the encoding injective.
    code(x) = 2x + 6
    encode(ℓ, m′, m) = code(ℓ) + code(m′)*25 + code(m)*625
    for ℓ ∈ Any[collect(0:ℓₘₐₓ); half_integer.(collect(1//2:(ℓₘₐₓ+1//2)))]
        mₘ = ℓ

        # These tests are old; the input array can be larger than necessary now.
        # # Check that ℓ < m′ₘₐₓ and ℓ ≠ mₘₐₓ throw errors
        # @test_throws "greater than 0 and less than or equal to 2ℓ+1=" WignerDMatrix(Array{ComplexF64}(undef, Int(2ℓ)+2, Int(2ℓ)+1), ℓ)
        # @test_throws "greater than 0 and less than or equal to 2ℓ+1=" WignerdMatrix(Array{Float64}(undef, Int(2ℓ)+2, Int(2ℓ)+1), ℓ)
        # @test_throws "in the input data must be 2ℓ+1=" WignerDMatrix(Array{ComplexF64}(undef, Int(2ℓ)+1, Int(2ℓ)+2), ℓ)
        # @test_throws "in the input data must be 2ℓ+1=" WignerdMatrix(Array{Float64}(undef, Int(2ℓ)+1, Int(2ℓ)+2), ℓ)

        # # Check that the input is at least as big as needed for the given ℓ
        @test_throws "The extent of the first dimension" WignerDMatrix(Array{ComplexF64}(undef, Int(2ℓ)+0, Int(2ℓ)+1), ℓ)
        @test_throws "The extent of the first dimension" WignerdMatrix(Array{Float64}(undef, Int(2ℓ)+0, Int(2ℓ)+1), ℓ)
        @test_throws "The extent of the second dimension" WignerDMatrix(Array{ComplexF64}(undef, Int(2ℓ)+1, Int(2ℓ)+0), ℓ)
        @test_throws "The extent of the second dimension" WignerdMatrix(Array{Float64}(undef, Int(2ℓ)+1, Int(2ℓ)+0), ℓ)

        # Check that a mismatch between integer/half-integer throws an error
        if ℓ>0 && ℓ isa Int
            @test_throws "The extent of the first dimension" WignerDMatrix(rand(ComplexF64, 2ℓ, 2ℓ+1), ℓ)
            @test_throws "The extent of the first dimension" WignerdMatrix(rand(Float64, 2ℓ, 2ℓ+1), ℓ)
        elseif ℓ isa HalfOddInteger
            @test_throws "The extent of the first dimension" WignerDMatrix(rand(ComplexF64, Int(2ℓ), Int(2ℓ+1)), ℓ)
            @test_throws "The extent of the first dimension" WignerdMatrix(rand(Float64, Int(2ℓ), Int(2ℓ+1)), ℓ)
        end
        @test_throws "The extent of the second dimension" WignerDMatrix(rand(ComplexF64, Int(2ℓ+1), Int(2ℓ)), ℓ)
        @test_throws "The extent of the second dimension" WignerdMatrix(rand(Float64, Int(2ℓ+1), Int(2ℓ)), ℓ)

        # Check that a data array with a dimension of 0 extent throws an error.
        @test_throws r"The extent of the second dimension.*; it is 0." WignerDMatrix(Array{ComplexF64}(undef, Int(2ℓ)+1, 0), ℓ)
        @test_throws r"The extent of the first dimension.*; it is 0." WignerDMatrix(Array{ComplexF64}(undef, 0, Int(2ℓ)+1), ℓ)
        @test_throws r"The extent of the second dimension.*; it is 0." WignerdMatrix(Array{Float64}(undef, Int(2ℓ)+1, 0), ℓ)
        @test_throws r"The extent of the first dimension.*; it is 0." WignerdMatrix(Array{Float64}(undef, 0, Int(2ℓ)+1), ℓ)

        for m′ₘ ∈ ℓₘᵢₙ(ℓ):ℓ
            # Make a big, dumb array full of the explicit indices.
            data = [
                encode(ℓ, m′, m)
                for m′ ∈ -m′ₘ:m′ₘ, m ∈ -mₘ:mₘ
            ]
            # Check that indexing works as expected.
            for (WignerMatrixType, NT) ∈ ((WignerDMatrix, ComplexF64), (WignerdMatrix, Float64))
                w = WignerMatrixType(NT.(data), ℓ; m′ₘₐₓ=m′ₘ, m′ₘᵢₙ=-m′ₘ, mₘₐₓ=mₘ, mₘᵢₙ=-mₘ)
                @test Base.parent(w) == data
                @test ell(w) == ℓ
                @test mpmax(w) == m′ₘ
                @test mmax(w) == ℓ
                @test mpmin(w) == -mpmax(w)
                @test mmin(w) == -mmax(w)
                for m ∈ -mₘ:mₘ
                    for m′ ∈ -m′ₘ:m′ₘ
                        @test w[m′, m] == encode(ℓ, m′, m)
                    end
                end
            end
        end

        for m′ₘ ∈ ℓₘᵢₙ(ℓ):ℓ
            for WignerMatrixType ∈ (WignerDMatrix, WignerdMatrix)
                data = rand(
                    WignerMatrixType<:WignerDMatrix ? ComplexF64 : Float64,
                    Int(2m′ₘ)+1, Int(2mₘ)+1
                )
                w = WignerMatrixType(data, ℓ; m′ₘₐₓ=m′ₘ, m′ₘᵢₙ=-m′ₘ, mₘₐₓ=mₘ, mₘᵢₙ=-mₘ)

                # Check that the data array is stored correctly.
                @test Base.parent(w) == data
                @test ell(w) == ℓ
                @test m′ₘₐₓ(w) == m′ₘ
                @test mₘₐₓ(w) == ℓ
                @test m′ₘᵢₙ(w) == -m′ₘₐₓ(w)
                @test mₘᵢₙ(w) == -mₘₐₓ(w)

                # These containers are deliberately not `AbstractArray`s, and their axes
                # deliberately do not meet the array interface's demand for an
                # `AbstractUnitRange{<:Integer}` that is its own axis — a half-odd axis
                # cannot be an index set at all.  A `WignerRange` is instead an ordinary
                # range of index values, indexed by position as `Base`'s ranges are, so that
                # its values can be collected and broadcast over.
                @test typeof(axes(w)) <: NTuple{2, AbstractUnitRange}
                @test axes.(axes(w), 1) == map(a -> Base.OneTo(length(a)), axes(w))
                @test all(collect(a) == [first(a) + k for k ∈ 0:length(a)-1] for a ∈ axes(w))
            end
        end
    end
end


### The shared array interface for the block containers.
#
# Every block container answers two questions about itself, and the whole array interface is
# written once in terms of the answers: `axis_roles` names its axes, in the order they are
# indexed, and `natural_axes` gives the matching ranges.  Before version 3 each container spelled
# out `size`, `length`, `ndims`, `axes` and their trailing-dimension forms for itself, which came
# to some forty near-identical methods differing only in rank and in which fields held the bounds.
#
# `getindex` and `setindex!` are deliberately *not* written this way, and neither are the
# accessors `ℓ`, `mₘₐₓ`, `sₘᵢₙ` and the rest.  Those are the hot path: the bounds check there
# reads the stored limits directly rather than building a range to test membership in, which is
# the difference the note above `inrange` measures at 145 ns against 0.9 ns per element.
#
# `HWedge` and `HAxis` are outside this union on purpose.  They are workspaces for the
# recursion rather than blocks handed to a caller, their storage is triangular, and `size` of
# one is the size of that storage rather than of any block; they keep the generic
# `AbstractWignerMatrix` methods above.

const BlockContainer = Union{
    WignerMatrix, WignerMatrixBatch,
    DegreeBlock, DegreeBlockBatch,
    SpinMatrix, SpinMatrixBatch,
}

# The roles are a property of the type, so these fold away at compile time.  A batched
# container's leading axis is an ordinary 1-based rotor position, not a natural index, which is
# why it is the one role whose range is a plain `UnitRange` rather than a `WignerRange`.
axis_roles(::Type{<:WignerMatrix}) = (:m′, :m)
axis_roles(::Type{<:WignerMatrixBatch}) = (:iᵣ, :m′, :m)
axis_roles(::Type{<:DegreeBlock}) = (:m,)
axis_roles(::Type{<:DegreeBlockBatch}) = (:iᵣ, :m)
axis_roles(::Type{<:SpinMatrix}) = (:s, :m)
axis_roles(::Type{<:SpinMatrixBatch}) = (:iᵣ, :s, :m)
axis_roles(w::BlockContainer) = axis_roles(typeof(w))

@inline natural_axes(w::WignerMatrix) =
    (WignerRange(w.m′ₘᵢₙ:w.m′ₘₐₓ), WignerRange(w.mₘᵢₙ:w.mₘₐₓ))
@inline natural_axes(w::WignerMatrixBatch) =
    (1:w.Nᵣ, WignerRange(w.m′ₘᵢₙ:w.m′ₘₐₓ), WignerRange(w.mₘᵢₙ:w.mₘₐₓ))
@inline natural_axes(v::DegreeBlock) = (WignerRange(v.mₘᵢₙ:v.mₘₐₓ),)
@inline natural_axes(v::DegreeBlockBatch) = (1:v.Nᵣ, WignerRange(v.mₘᵢₙ:v.mₘₐₓ))
@inline natural_axes(b::SpinMatrix) =
    (WignerRange(b.sₘᵢₙ:b.sₘₐₓ), WignerRange(b.mₘᵢₙ:b.mₘₐₓ))
@inline natural_axes(b::SpinMatrixBatch) =
    (1:b.Nᵣ, WignerRange(b.sₘᵢₙ:b.sₘₐₓ), WignerRange(b.mₘᵢₙ:b.mₘₐₓ))

# The extent is taken as the difference of the endpoints rather than as `length` of the range,
# because that is an `Int` by construction for a half-odd-integer axis, where `Int(ℓ)` throws.
@inline axis_extent(r) = Int(last(r) - first(r)) + 1

@inline Base.axes(w::BlockContainer) = natural_axes(w)
@inline Base.size(w::BlockContainer) = map(axis_extent, natural_axes(w))
Base.length(w::BlockContainer) = prod(size(w))
Base.ndims(w::BlockContainer) = length(axis_roles(w))
Base.ndims(::Type{T}) where {T<:BlockContainer} = length(axis_roles(T))
# Trailing dimensions behave as they do for `AbstractArray`, so that generic code written
# against a plain array works unchanged: `axes(w, d)` is `OneTo(1)` and `size(w, d)` is `1`.
Base.axes(w::BlockContainer, d::Integer) = d ≤ ndims(w) ? axes(w)[d] : Base.OneTo(1)
Base.size(w::BlockContainer, d::Integer) = d ≤ ndims(w) ? size(w)[d] : 1
