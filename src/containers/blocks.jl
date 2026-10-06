"""
    AbstractBlock{IT, NT, ST}

Abstract supertype of the blocks, the containers that hold the values of one ``ℓ``, indexed
by their natural indices ``m′``, ``m`` and ``s`` rather than by position.  There are six:
[`WignerMatrix`](@ref), indexed `w[m′, m]`, [`WignerMatrixBatch`](@ref), `w[iᵣ, m′, m]`,
[`DegreeBlock`](@ref), `v[m]`, [`DegreeBlockBatch`](@ref), `v[iᵣ, m]`, [`SpinMatrix`](@ref),
`b[s, m]`, and [`SpinMatrixBatch`](@ref), `b[iᵣ, s, m]`.
- `IT` is the index type, `Int` or [`HalfOddInteger`](@ref), which is the type of ``ℓ`` and
  of every natural index.
- `NT` is the number type, such as `ComplexF64` for a block of ``𝔇`` or `Float64` for one
  of ``d``.
- `ST` is the type of the storage, `parent(w)`, which is one-dimensional and 1-based (see
  below).

A block's storage is one-dimensional, and element k of the block, in the order of `Array(w)`
(the column-major order of its axes), is entry k of the storage; entries beyond `length(w)`
are not read.  The constructors take a 1-based array of any shape with at least `length(w)`
entries and use it in its linear order: the storage is `vec` of that array, which shares its
memory, and an array of exactly the block's shape is the block.  A calculator's blocks are
the leading entries of the calculator's one buffer, so `parent(w)` of such a block is that
buffer.  The blocks of [`D`](@ref) and [`d`](@ref) are consecutive parts of one vector,
which holds them all, so `parent(w)` of each is a view of its own part.
[`array_view`](@ref)`(w)` is the block's entries of the storage, reshaped to the block's
shape, without a copy: contiguous for a calculator's block, for the blocks of `D` and `d`,
and for the results of `copy`, `similar`, and [`relabel`](@ref).  A slice of a batch,
`w[iᵣ]` or `b[iᵣ]`, and the row `b[s, :]` of a `SpinMatrix`, are views of an arithmetic
progression of the storage, with the step `Nᵣ` or the number of spin weights.  The one
exception is the slice `b[:, s, :]` of a [`SpinMatrixBatch`](@ref), whose elements are not
one such progression: its storage is a strided matrix view of exactly its shape, and its
`array_view` is that matrix.

These types are *not* `AbstractArray`s, because half-integer indices cannot satisfy that
interface: `axes` must be integer ranges there, and `-3//2:3//2` is not one.  Consequently
linear algebra (`*`, `'`, `\\`, `lu`) does not apply to them; [`array_view`](@ref) gives the
numbers as an ordinary 1-based array, and [`relabel`](@ref) puts the labels back.

# Methods

Every block has `parent(w)`, the storage; `ℓ(w)`, the degree; `ℓₘᵢₙ(w)`, the smallest
degree of the index type, 0 or 1//2; `eltype(w)`, the number type; [`ishalfinteger`](@ref);
and `summary` and `show`.  Each block defines its own indexing, since the natural indices
differ from one to the next, and all six share an array-like interface as well:
- `axes(w)` are the ranges of the natural indices, one per axis, with `1:Nᵣ` for the rotor
  axis of a batch, and `size(w)`, `length(w)` and `ndims(w)` describe the block, whose
  storage may be larger, as a calculator's buffer is;
- the limits are read with [`m′ₘₐₓ`](@ref), [`m′ₘᵢₙ`](@ref), [`mₘₐₓ`](@ref), [`mₘᵢₙ`](@ref),
  [`sₘₐₓ`](@ref) and [`sₘᵢₙ`](@ref), for the axes a block has, and [`isbatched`](@ref) says
  whether it has a rotor axis;
- `checkbounds(Bool, w, inds...)` tests natural indices;
- `Array(w)` and `collect(w)` copy the block into an ordinary 1-based array, as do
  `Matrix(w)` for a two-axis block and `Vector(v)` for a `DegreeBlock`;
- `copy(w)` and `similar(w)` keep the labels, and iteration visits the elements in the order
  of `Array(w)`;
- `==`, `isequal`, `hash` and `≈` compare the labels — the kind of block, ``ℓ`` and the axes
  — as well as the numbers, so the same numbers under different labels are not equal;
- broadcasting gives an ordinary 1-based `Array`, and a broadcast that combines blocks, or
  writes into one with `.=`, requires their labels to agree, since it pairs the elements by
  position.

Note that there is deliberately *no* linear indexing (`w[i]`), for the reason given above.
"""
abstract type AbstractBlock{IT<:IntegerHalf, NT, ST<:AbstractArray{NT}} end
# Note that this is deliberately *not* a subtype of `AbstractMatrix`: the natural indices
# `(m′, m)` may be half-odd-integers, which cannot satisfy the `AbstractArray` interface
# (integer `axes`).  The array-like methods that make sense are defined explicitly below.

### General methods for all AbstractBlock types

Base.parent(w::AbstractBlock) = w.parent

ℓ(w::AbstractBlock{IT}) where {IT} = w.ℓ
# The smallest degree of the index type of the block.
ℓₘᵢₙ(::AbstractBlock{IT}) where {IT} = lowest_index(IT)

# The limits of the axes are fields of the blocks that have those axes.  Only a
# `WignerMatrix` and a `WignerMatrixBatch` have an m′ axis, and the other blocks refuse the
# m′ accessors with an explanation, at the end of this file, rather than with an error about
# a missing field.
mₘₐₓ(w::AbstractBlock{IT}) where {IT} = w.mₘₐₓ
mₘᵢₙ(w::AbstractBlock{IT}) where {IT} = w.mₘᵢₙ

ishalfinteger(::AbstractBlock{IT}) where {IT<:Integer} = false
ishalfinteger(::AbstractBlock{IT}) where {IT<:HalfOddInteger} = true

Base.eltype(::AbstractBlock{IT, NT, ST}) where {IT, NT, ST} = NT
Base.eltype(::Type{<:AbstractBlock{IT, NT, ST}}) where {IT, NT, ST} = NT

# The range of one natural index, as `axes` of a block reports it.  `T` is the index type,
# `Int` or `HalfOddInteger`.
struct WignerRange{T<:IntegerHalf} <: AbstractUnitRange{T}
    start::T
    stop::T

    WignerRange(r::UnitRange{T}) where {T} = new{T}(r.start, r.stop)
end
# A `WignerRange` is indexed by position, as `Base`'s ranges are — `r[1]` is `first(r)`, and
# `firstindex` and `lastindex` below say so — so its own axes are 1-based.  Making `axes(r)`
# the range itself, after the manner of an identity-offset range, would leave `axes` and
# `getindex` disagreeing: broadcasting a function over the range would then read positions
# that were really values, silently wrong for an integer range and an error for a half-odd
# one.  The natural, possibly half-odd, bounds are what `axes` of a *container* reports, and
# those are what `axes_string` below displays.
@inline Base.axes(r::WignerRange) = (axes1(r),)
@inline axes1(r::WignerRange) = Base.OneTo(length(r))
# Every `WignerRange` has unit step, so it is shown as a `UnitRange` is, without the step,
# which is how the summaries of the containers and the documentation write the same ranges.
Base.show(io::IO, r::WignerRange) = print(io, first(r), ":", last(r))
# The axes of a container, as its summary writes them: `(-2:2)×(-2:2)`.
axes_string(axes::Tuple) = join(("($(first(a)):$(last(a)))" for a ∈ axes), "×")
Base.firstindex(r::WignerRange) = 1
Base.lastindex(r::WignerRange) = length(r)
# As for `UnitRange{HalfOddInteger}` in `half_odd_integer.jl`: `Base`'s `step` for an
# `AbstractUnitRange{T}` is `oneunit(T) - zero(T)`, and its `length` independently forms
# `oneunit(zero(stop) - zero(start))`; both reach for values that `HalfOddInteger`
# deliberately lacks, so both must be given directly.  Without these two, `show` of a
# half-integer axis — `axes(w[ℓ, :])` for a half-integer `ModeWeights`, say — and
# `lastindex` above both throw.  The integer case is left to `Base`.
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
# so membership of a `T` is just the bracket.  (The generic `in(::Real,
# ::AbstractRange{<:Real})` would additionally test integrality of the offset, which is
# automatic here: `x - first(r)` is an `Integer` by construction for both index types.)
@inline Base.in(x::T, r::WignerRange{T}) where {T<:IntegerHalf} = first(r) ≤ x ≤ last(r)
# Any other value is a member when it is `==` to one, as for `Base`'s ranges: a `Rational`
# or a float equal to a half-odd member — `3//2 ∈ axes(w, 1)` — or a float equal to an
# integer member — `1.0 ∈ axes(w, 1)`.  A value of the other kind of index never is, which
# for a `HalfOddInteger` in an integer range follows from `isinteger`.
@inline Base.in(x::Real, r::WignerRange{HalfOddInteger}) = in_half_odd_range(x, first(r), last(r))
@inline Base.in(x::Real, r::WignerRange{<:Integer}) = isinteger(x) && first(r) ≤ x ≤ last(r)
# The general method just above ties with the `IntegerHalf` one when the value and the range
# are both half-odd, and this settles it, as the last method below does for integers.
@inline Base.in(x::HalfOddInteger, r::WignerRange{HalfOddInteger}) = first(r) ≤ x ≤ last(r)
# `Base` has `in(::Integer, ::AbstractUnitRange{<:Integer})`, which is neither more nor less
# specific than either method above, so without this one `1 ∈ axes(w, 1)` on an
# integer-indexed container is an ambiguity error rather than an answer.  Aqua's ambiguity
# check catches its absence.
@inline Base.in(x::Integer, r::WignerRange{<:Integer}) = first(r) ≤ x ≤ last(r)
# That one is deliberately loose — it has to cover `in(::Int8, ::WignerRange{Int})` and the
# rest of the intersection with `Base`'s method — which leaves it ambiguous in turn with the
# `IntegerHalf` method above whenever the value and the range share one integer type.  This
# third method ties the two together and so is more specific than both, which settles it.
# Aqua does *not* catch its absence: `Test.detect_ambiguities` reports nothing for that pair
# on either Julia 1.12 or 1.13, even though `0 ∈ WignerRange(-2:3)` throws without it.  The
# direct calls in `test/wigner/wigner_matrix.jl` are the only guard.
@inline Base.in(x::T, r::WignerRange{T}) where {T<:Integer} = first(r) ≤ x ≤ last(r)


### Bounds checking
#
# `m ∈ axes(w, 2)` is the obvious way to write the checks below, but it builds a fresh
# `WignerRange` and then calls the generic `in`, which together cost about 160 times the
# load they guard: 145 ns per element for a half-integer block, against 0.9 ns for a load
# from an ordinary matrix.  The helper here tests the stored limits directly instead, and
# gives identical answers.

# `lo ≤ m ≤ hi`.  There is no parity test to perform: `m` is of the container's own index
# type, so a whole number cannot reach a half-integer container in the first place.
@inline inrange(m, lo, hi) = lo ≤ m ≤ hi


### Validation of the limits
#
# A block with an m′ axis — a `WignerMatrix`, a `WignerMatrixBatch`, and the calculators and
# wedges that fill them — must have its limits in order, within ±ℓₘₐₓ, and bracketing ±ℓₘᵢₙ:
# the recurrence starts from the row m′ = 0 for integer indices and from the pair of rows m′
# = ±1/2 for half-odd-integers, and computes the rest outward from there, so a range that
# misses them cannot be filled.  The same holds for the m axis of those blocks.  Their
# degree is checked by `validate_degree`, and then each axis by `validate_axis`.  The
# one-axis blocks — `DegreeBlock`, `SpinMatrix` and their batches — label storage that no
# recurrence fills, so they check only that their m range is in order and within ±ℓ, in
# `validate_m_range`; their spin axis need only be in order.  These are all checks of the
# caller's arguments, and the messages are built only on the branches that throw them.

# The part of the bracketing rule that a refused limit broke, for the messages below.
bracket_rule(::Type{<:Integer}, name) =
    "the range of $name must include 0, where the recurrence starts"
bracket_rule(::Type{HalfOddInteger}, name) =
    "the range of $name must include both -1//2 and 1//2, where the recurrence starts"

# The largest degree must be non-negative.  Every caller checks it before the limits of any
# axis, because the limits default to ±ℓₘₐₓ, and a bad ℓₘₐₓ would otherwise be reported as a
# bad limit that the caller never gave.
function validate_degree(ℓₘₐₓ::IT) where {IT<:IntegerHalf}
    if ℓₘₐₓ < lowest_index(IT)
        throw(ArgumentError("ℓₘₐₓ=$ℓₘₐₓ must be non-negative."))
    end
    nothing
end

# The limits `hi` and `lo` of the axis `name`, "m′" or "m", of the blocks of degree up to
# `ℓₘₐₓ` must be in order, must bracket ±ℓₘᵢₙ (see above), and must lie within ±ℓₘₐₓ; they
# are checked in that order.
function validate_axis(ℓₘₐₓ::IT, hi::IT, lo::IT, name::String) where {IT<:IntegerHalf}
    if hi < lo
        throw(ArgumentError("$(name)ₘₐₓ=$hi is less than $(name)ₘᵢₙ=$lo."))
    end
    if hi < lowest_index(IT)
        throw(ArgumentError(
            "$(name)ₘₐₓ=$hi is too small for this index type, $IT: "
            * "$(bracket_rule(IT, name))."
        ))
    end
    if lo > -lowest_index(IT)
        throw(ArgumentError(
            "$(name)ₘᵢₙ=$lo is too large for this index type, $IT: "
            * "$(bracket_rule(IT, name))."
        ))
    end
    if abs(hi) > ℓₘₐₓ
        throw(ArgumentError("|$(name)ₘₐₓ|=|$hi| is too large for ℓₘₐₓ=$ℓₘₐₓ."))
    end
    if abs(lo) > ℓₘₐₓ
        throw(ArgumentError("|$(name)ₘᵢₙ|=|$lo| is too large for ℓₘₐₓ=$ℓₘₐₓ."))
    end
    nothing
end

# The one-axis rule for the m axis of a `DegreeBlock`, a `SpinMatrix` or one of their
# batches.
function validate_m_range(ℓ::IT, mₘₐₓ::IT, mₘᵢₙ::IT) where {IT<:IntegerHalf}
    if ℓ < lowest_index(IT)
        throw(ArgumentError("ℓ=$ℓ must be non-negative."))
    end
    if mₘₐₓ < mₘᵢₙ
        throw(ArgumentError("mₘₐₓ=$mₘₐₓ is less than mₘᵢₙ=$mₘᵢₙ."))
    end
    if mₘᵢₙ < -ℓ || mₘₐₓ > ℓ
        throw(ArgumentError(
            "The range of m, $mₘᵢₙ:$mₘₐₓ, must lie within -ℓ:ℓ, which is $(-ℓ):$ℓ for ℓ=$ℓ."
        ))
    end
    nothing
end

# The spin axis of a `SpinMatrix` or a `SpinMatrixBatch` may be any run of spin weights (see
# the `SpinMatrix` docstring), so it is checked only for its order.
function validate_s_range(sₘₐₓ::IT, sₘᵢₙ::IT) where {IT<:IntegerHalf}
    if sₘₐₓ < sₘᵢₙ
        throw(ArgumentError("sₘₐₓ=$sₘₐₓ is less than sₘᵢₙ=$sₘᵢₙ."))
    end
    nothing
end


### Storage
#
# A block's storage is one-dimensional: element k of the block, in the order of `Array(w)`,
# is entry k of `parent(w)`, and the entries beyond `length(w)` are not read.  The only
# other storage is the matrix of a slice `b[:, s, :]` of a `SpinMatrixBatch` (see the
# `DegreeBlockBatch` that it builds), whose elements are not one arithmetic progression
# through the storage of `b`; it is a strided view of exactly the slice's shape, and is read
# as a matrix.  Every method that reads the storage goes through the four functions below,
# which are the one place that tells the two apart: `storage_index` for one element,
# `elements` for all of them as an array of the block's shape (which `array_view`, the
# reductions, broadcasting, the copying forms, the comparisons, and `show` read),
# `leading_slice` for the slices of the leading axis, and `holds`, with `check_storage`, for
# the comparison of the storage with the block.  The one exception is iteration, which reads
# entry k of the storage by its linear index, since that is element k in either form (see
# `iterate`).
#
# The accessors index the storage itself, at the position `storage_index` gives, rather than
# `elements`: building the reshaped view at each access is free for a `Vector`, but for
# other storage, such as `vec` of a view of a matrix, it was measured to allocate at every
# access on Julia 1.10, and for a strided view on 1.13 as well.
#
# The inner constructors check that the storage holds every element the limits describe,
# because the accessors and iteration read it under `@inbounds`; they are also what `copy`,
# `similar`, `relabel`, and the slices `w[iᵣ]`, `b[s, :]`, `b[:, s, :]`, and `b[iᵣ]` are
# built with.  The ordering of the limits and their relation to ℓ are left to the outer
# constructors, with one exception: an axis of negative extent is refused here, because two
# of them would multiply to a positive `length`, which iteration would then read.  A vector
# can be resized after the constructor has compared it with the limits, so every access to a
# block over a vector compares its length with the block again, and any access to a block
# whose vector has become shorter than the block is refused; see `check_storage`.
#
# These run on every view that `w[iᵣ]` or `b[s, :]` builds, so the messages are formatted
# only on the way to an error.

# The axes of a block as its inner constructor describes them: `(hi, lo, name)` for a
# natural index, whose extent is `hi - lo + 1`, and `(Nᵣ,)` for the rotor axis of a batch.
@inline described_extent(a::Tuple{Any, Any, String}) = Int(a[1] - a[2]) + 1
@inline described_extent(a::Tuple{Int}) = a[1]

@inline function check_extent(parent::AbstractArray, axes::Tuple)
    dims = map(described_extent, axes)
    (all(≥(0), dims) && holds(parent, dims)) || extent_error(parent, axes)
    nothing
end

# Whether `parent` holds a block whose axes have the extents `dims`: a vector must be at
# least as long as the block, and the matrix of a slice `b[:, s, :]` must have its shape.
@inline holds(parent::AbstractVector, dims::Tuple) = length(parent) ≥ prod(dims)
@inline holds(parent::AbstractMatrix, dims::Tuple{Int, Int}) = size(parent) == dims

@noinline function extent_error(parent, axes::Tuple)
    for a ∈ axes
        n = described_extent(a)
        n ≥ 0 && continue
        if a isa Tuple{Int}
            throw(ArgumentError("The number of rotors Nᵣ=$n must not be negative."))
        end
        hi, lo, name = a
        throw(ArgumentError(
            "$(name)ₘₐₓ=$hi is less than $(name)ₘᵢₙ=$lo by more than one, which would give "
            * "the block an axis of extent $n."
        ))
    end
    dims = map(described_extent, axes)
    if parent isa AbstractVector && length(axes) == 1
        throw(DimensionMismatch(
            "The input data must have at least $(extent_formula(only(axes))) entries; it "
            * "has $(length(parent))."
        ))
    elseif parent isa AbstractVector
        throw(DimensionMismatch(
            "The input data must have at least $(prod(dims)) = $(join(dims, "×")) entries, "
            * "the number of elements of a block with "
            * "$(serial_list(map(axis_description, axes))); it has $(length(parent))."
        ))
    else
        throw(DimensionMismatch(
            "The storage of a slice must have the shape $dims of the block; it has the "
            * "shape $(size(parent))."
        ))
    end
end
# The extent of the one axis of a `DegreeBlock`, as the message above writes it:
# `mₘₐₓ-mₘᵢₙ+1=2-(-2)+1=5`, with a negative lower limit parenthesized.
function extent_formula((hi, lo, name)::Tuple{Any, Any, String})
    lo_text = lo < 0 ? "($lo)" : "$lo"
    "$(name)ₘₐₓ-$(name)ₘᵢₙ+1=$hi-$lo_text+1=$(Int(hi - lo) + 1)"
end
# The axes of a block of two or three axes, as the message above lists them: `Nᵣ=3 rotors,
# m′ ∈ -2:2, and m ∈ -2:2`.
axis_description(a::Tuple{Int}) = "Nᵣ=$(a[1]) rotors"
axis_description((hi, lo, name)::Tuple{Any, Any, String}) = "$name ∈ $lo:$hi"
serial_list(items) = length(items) == 2 ? join(items, " and ") : join(items, ", ", ", and ")

# The position in the storage of the element at the 1-based positions along the block's
# axes: its linear position in storage of one dimension, for a block of one, two, or three
# axes, and the positions themselves in the matrix of a slice `b[:, s, :]`.  Every accessor
# has compared the positions with the block's limits before it asks.  The linear position is
# computed here rather than by `LinearIndices`, whose own bounds check, which remains when
# bounds checking is forced on, as `Pkg.test` forces it, can crash Enzyme's reverse mode as
# it compiles a loop over the blocks of a calculator.  Only a `DegreeBlockBatch` may have a
# matrix as storage (see its type bound below); the methods for that storage, here and in
# the functions that follow, are written for any block over a matrix, because the block
# types are defined after them.
@inline storage_index(::AbstractBlock{IT, NT, <:AbstractVector}, i::Int) where {IT, NT} = i
@inline storage_index(w::AbstractBlock{IT, NT, <:AbstractVector}, i::Int, j::Int) where {IT, NT} =
    i + size(w)[1] * (j - 1)
@inline function storage_index(
    w::AbstractBlock{IT, NT, <:AbstractVector}, i::Int, j::Int, k::Int
) where {IT, NT}
    n₁, n₂ = size(w)
    i + n₁ * ((j - 1) + n₂ * (k - 1))
end
@inline storage_index(v::AbstractBlock{IT, NT, <:AbstractMatrix}, iᵣ::Int, j::Int) where {IT, NT} =
    CartesianIndex(iᵣ, j)

# The elements of `w`, 1-based, as an array of its shape; `array_view` is this with the
# storage compared with the block first.
@inline elements(w::AbstractBlock{IT, NT, <:AbstractVector}) where {IT<:IntegerHalf, NT} =
    reshape(leading_view(parent(w), length(w)), size(w))
@inline elements(v::AbstractBlock{IT, NT, <:AbstractMatrix}) where {IT<:IntegerHalf, NT} =
    parent(v)

# The elements of `w` whose first index is `i`, in order, as a vector: an arithmetic
# progression through storage of one dimension, whose step is the extent of the first axis,
# or a row of the matrix of a slice `b[:, s, :]`.  These are the storage of the slices
# `w[iᵣ]`, `v[iᵣ]`, and `b[iᵣ]` of a batch, and `b[s, :]` of a `SpinMatrix`.
@inline function leading_slice(w::AbstractBlock{IT, NT, <:AbstractVector}, i::Int) where {IT, NT}
    N = size(w, 1)
    view(parent(w), i:N:(i + N * (length(w) ÷ N - 1)))
end
@inline leading_slice(v::AbstractBlock{IT, NT, <:AbstractMatrix}, i::Int) where {IT, NT} =
    view(parent(v), i, :)

# The leading `n` entries of `p`, as a view, with no bounds check of its own: every caller
# has compared the storage with the block first, in `check_storage`, which the extension for
# Enzyme declares inactive (see the note there).  The extension for ChainRules
# differentiates it as `view(p, Base.OneTo(n))`.
@inline leading_view(p::AbstractVector, n::Int) = Base.unsafe_view(p, Base.OneTo(n))

# Only a vector can be resized after the constructor has compared its length with the limits
# — a caller's `Vector`, or a calculator's buffer, which `parent(w)` hands out — so every
# access compares its length with the block again, and refuses any access once it is shorter
# than the block.
@inline function check_storage(w::AbstractBlock{IT, NT, <:AbstractVector}) where {IT<:IntegerHalf, NT}
    if length(parent(w)) < length(w)
        throw(storage_error(w))
    end
    nothing
end
@inline check_storage(v::AbstractBlock{IT, NT, <:AbstractMatrix}) where {IT<:IntegerHalf, NT} =
    nothing
@noinline function storage_error(w::AbstractBlock)
    DimensionMismatch(
        "The storage of this $(sprint(summary, w)) has length $(length(parent(w))), but "
        * "the block has $(length(w)) elements, which are its leading entries.  A block "
        * "uses its storage without copying it, so the vector must not be resized."
    )
end


@doc raw"""
    WignerMatrix{IT, NT, ST} <: AbstractBlock{IT, NT, ST}
    WignerMatrix(parent, ℓ; m′ₘₐₓ=ℓ, m′ₘᵢₙ=-m′ₘₐₓ, mₘₐₓ=ℓ, mₘᵢₙ=-mₘₐₓ)

General concrete subtype of [`AbstractBlock`](@ref) for Wigner rotation matrices,
which can include D-matrices (when `NT` is complex) or d-matrices (when `NT` is real).
- `IT` is the index type, `Int` or [`HalfOddInteger`](@ref).
- `NT` is the number type of the elements.
- `ST` is the type of the storage, an `AbstractVector{NT}` (see [`AbstractBlock`](@ref)).

The block holds some or all of the Wigner matrix for a specific ``ℓ``, as a rectangular
array whose first dimension corresponds to the `m′` index and whose second corresponds to
the `m` index; how it is stored is described under [`AbstractBlock`](@ref).  The allowed
ranges of `m′` and `m` are governed by the fields `m′ₘₐₓ`, `m′ₘᵢₙ`, `mₘₐₓ`, and `mₘᵢₙ`,
which must satisfy
```math
\begin{aligned}
-ℓₘₐₓ &≤ m′ₘᵢₙ ≤ -ℓₘᵢₙ ≤ ℓₘᵢₙ ≤ m′ₘₐₓ ≤ ℓₘₐₓ, \\
-ℓₘₐₓ &≤ mₘᵢₙ ≤ -ℓₘᵢₙ ≤ ℓₘᵢₙ ≤ mₘₐₓ ≤ ℓₘₐₓ,
\end{aligned}
```
where `ℓₘᵢₙ` is either 0 or 1//2 depending on whether `IT` is an integer or half-integer
type.  Both rows `±ℓₘᵢₙ` must be included because the recurrence seeds the half-integer
ladder from the pair of rows `m′ = ±1/2`; for integers this reduces to the familiar `m′ₘᵢₙ ≤
0 ≤ m′ₘₐₓ`.

The constructor wraps `parent` without copying it.  `parent` is a 1-based array of any shape
with at least `(m′ₘₐₓ-m′ₘᵢₙ+1) × (mₘₐₓ-mₘᵢₙ+1)` entries, read in its linear order, so that a
matrix of exactly that size is the block.  `ℓ` and the limits must all be integers of type
`Int`, or all half-odd-integers, each a [`HalfOddInteger`](@ref) or a `Rational{Int}` with
denominator 2; the limits default to the whole block, and each lower limit defaults to minus
the corresponding upper one.  The keywords may also be spelled `mp_max`, `mp_min`, `m_max`,
and `m_min`; where both spellings of one are given, the Unicode one is used.
"""
struct WignerMatrix{IT, NT, ST<:AbstractVector{NT}} <: AbstractBlock{IT, NT, ST}
    parent::ST
    ℓ::IT
    m′ₘₐₓ::IT
    m′ₘᵢₙ::IT
    mₘₐₓ::IT
    mₘᵢₙ::IT
    # The parent is indexed as 1-based throughout, much of it under `@inbounds`, so an
    # offset array would be read and written outside its storage; `ModeWeights` and
    # `HarmonicValues` refuse one in the same way.  For the same reason the parent must
    # reach every element of the block (see "Storage" above).
    function WignerMatrix{IT, NT, ST}(parent, ℓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ) where {IT, NT, ST}
        Base.require_one_based_indexing(parent)
        check_extent(parent, ((m′ₘₐₓ, m′ₘᵢₙ, "m′"), (mₘₐₓ, mₘᵢₙ, "m")))
        new{IT, NT, ST}(parent, ℓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    end
end

"""
    Matrix(w::WignerMatrix)

Materialize the block of the Wigner matrix represented by `w` as an ordinary `Matrix`.  The
rows and columns are in order of increasing `m′` and `m`, so that the element `w[m′, m]` is
at `[Int(m′-m′ₘᵢₙ)+1, Int(m-mₘᵢₙ)+1]`.
"""
Base.Matrix(w::WignerMatrix) = Matrix(array_view(w))

# Indexing by `(m′, m)` is what a `WignerMatrix` means by two indices, and only that.  The
# other blocks index differently — `[iᵣ, m′, m]` for `WignerMatrixBatch`, `[s, m]` for
# `SpinMatrix`, and `[iᵣ, m]` for `DegreeBlockBatch` — and define their own methods; were
# these methods on `AbstractBlock`, a two-index call on a three-index block would silently
# read the wrong element instead of being an error.
Base.checkbounds(::Type{Bool}, w::WignerMatrix{IT}, m′::IT, m::IT) where {IT} =
    inrange(m′, w.m′ₘᵢₙ, w.m′ₘₐₓ) && inrange(m, w.mₘᵢₙ, w.mₘₐₓ)

@propagate_inbounds function Base.getindex(
    w::WignerMatrix{IT}, m′::IT, m::IT
) where {IT<:IntegerHalf}
    @boundscheck checkbounds(Bool, w, m′, m) || throw(BoundsError(w, (m′, m)))
    @boundscheck check_storage(w)
    @inbounds parent(w)[storage_index(w, Int(m′ - w.m′ₘᵢₙ) + 1, Int(m - w.mₘᵢₙ) + 1)]
end

@propagate_inbounds function Base.setindex!(
    w::WignerMatrix{IT}, v, m′::IT, m::IT
) where {IT<:IntegerHalf}
    @boundscheck checkbounds(Bool, w, m′, m) || throw(BoundsError(w, (m′, m)))
    @boundscheck check_storage(w)
    @inbounds parent(w)[storage_index(w, Int(m′ - w.m′ₘᵢₙ) + 1, Int(m - w.mₘᵢₙ) + 1)] = v
end


### Indexing with indices of other types.
#
# The accessors above take indices of the container's own type, `Int` or `HalfOddInteger`.
# A half-integer container is indexed by `HalfOddInteger`s, but `w[1//2, -3//2]` is what a
# caller naturally writes, so each block container also has methods that take any
# `IndexType`, convert each natural index with `checked_index`, and re-dispatch; either
# natural index may be written either way.  An index that is not of the container's kind, or
# that the index methods would refuse — a narrow integer, or a `Rational` whose denominator
# is not 2 or whose integer type is not `Int` — is refused there with the reason, rather
# than with the bare `MethodError` that dispatch alone would give.  These methods are less
# specific than the accessors of the container's own index type, so an index of that type,
# which is what every loop over the elements passes, never reaches them.  That needs the
# bound `IT<:IntegerHalf` on the accessors: without it, their signature would admit index
# types outside `IndexType`, and would not be the more specific of the two.

@propagate_inbounds Base.getindex(
    w::WignerMatrix{IT}, m′::IndexType, m::IndexType
) where {IT} = w[checked_index(IT, m′, w, "m′"), checked_index(IT, m, w, "m")]
@propagate_inbounds Base.setindex!(
    w::WignerMatrix{IT}, v, m′::IndexType, m::IndexType
) where {IT} = (w[checked_index(IT, m′, w, "m′"), checked_index(IT, m, w, "m")] = v)

@index_methods function WignerMatrix(
    parent::ST, ℓ::IT;
    mp_max::IndexType=ℓ, m′ₘₐₓ::IndexType=mp_max,
    mp_min::IndexType=-m′ₘₐₓ, m′ₘᵢₙ::IndexType=mp_min,
    m_max::IndexType=ℓ, mₘₐₓ::IndexType=m_max,
    m_min::IndexType=-mₘₐₓ, mₘᵢₙ::IndexType=m_min
) where {IT<:IndexType, NT, ST<:AbstractArray{NT}}
    validate_degree(ℓ)
    validate_axis(ℓ, m′ₘₐₓ, m′ₘᵢₙ, "m′")
    validate_axis(ℓ, mₘₐₓ, mₘᵢₙ, "m")
    Base.require_one_based_indexing(parent)
    p = vec(parent)
    WignerMatrix{IT, NT, typeof(p)}(p, ℓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
end

# A block with the labels of `w` over the storage `p`, through the inner constructor, which
# checks that the storage reaches every element; `copy`, `similar` and `relabel` of every
# block are written with this (see "The shared array interface" at the end of this file).
rewrap(w::WignerMatrix{IT}, p) where {IT} =
    WignerMatrix{IT, eltype(p), typeof(p)}(p, w.ℓ, w.m′ₘₐₓ, w.m′ₘᵢₙ, w.mₘₐₓ, w.mₘᵢₙ)


"""
    WignerMatrixBatch{IT, NT, ST} <: AbstractBlock{IT, NT, ST}
    WignerMatrixBatch(parent, ℓ; m′ₘₐₓ=ℓ, m′ₘᵢₙ=-m′ₘₐₓ, mₘₐₓ=ℓ, mₘᵢₙ=-mₘₐₓ)

`Nᵣ` Wigner matrices of one ``ℓ``, stored together and indexed as `w[iᵣ, m′, m]`.  This is
what [`recurrence!`](@ref) returns for a calculator built from a vector of rotor data, of
any length, for either kind of index.
- `IT` is the index type, `Int` or [`HalfOddInteger`](@ref).
- `NT` is the number type of the elements.
- `ST` is the type of the storage, an `AbstractVector{NT}` (see [`AbstractBlock`](@ref)).

The elements are ordered `[iᵣ, m′, m]`, exactly as in the calculator.  The constructor takes
a 1-based array of any shape, whose first extent is `Nᵣ` and which has at least as many
entries as the block, read in its linear order, so that an array of exactly the block's
shape is the block; it takes the indices and the limits, with their ASCII spellings, exactly
as [`WignerMatrix`](@ref) does, and under the same rules.  `w[iᵣ]` gives the
[`WignerMatrix`](@ref) view of one rotor's matrix, which is then indexed naturally as
`w[iᵣ][m′, m]`.

See also [`WignerMatrix`](@ref) and [`WignerSeries`](@ref).
"""
struct WignerMatrixBatch{IT, NT, ST<:AbstractVector{NT}} <: AbstractBlock{IT, NT, ST}
    parent::ST
    ℓ::IT
    m′ₘₐₓ::IT
    m′ₘᵢₙ::IT
    mₘₐₓ::IT
    mₘᵢₙ::IT
    Nᵣ::Int
    # As for `WignerMatrix`: the parent must be 1-based, and must reach every element.
    function WignerMatrixBatch{IT, NT, ST}(parent, ℓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, Nᵣ) where {IT, NT, ST}
        Base.require_one_based_indexing(parent)
        check_extent(parent, ((Nᵣ,), (m′ₘₐₓ, m′ₘᵢₙ, "m′"), (mₘₐₓ, mₘᵢₙ, "m")))
        new{IT, NT, ST}(parent, ℓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, Nᵣ)
    end
end

@index_methods function WignerMatrixBatch(
    parent::ST, ℓ::IT;
    mp_max::IndexType=ℓ, m′ₘₐₓ::IndexType=mp_max,
    mp_min::IndexType=-m′ₘₐₓ, m′ₘᵢₙ::IndexType=mp_min,
    m_max::IndexType=ℓ, mₘₐₓ::IndexType=m_max,
    m_min::IndexType=-mₘₐₓ, mₘᵢₙ::IndexType=m_min
) where {IT<:IndexType, NT, ST<:AbstractArray{NT}}
    validate_degree(ℓ)
    validate_axis(ℓ, m′ₘₐₓ, m′ₘᵢₙ, "m′")
    validate_axis(ℓ, mₘₐₓ, mₘᵢₙ, "m")
    Base.require_one_based_indexing(parent)
    Nᵣ = size(parent, 1)
    p = vec(parent)
    WignerMatrixBatch{IT, NT, typeof(p)}(p, ℓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, Nᵣ)
end

Nᵣ(w::WignerMatrixBatch) = w.Nᵣ

Base.checkbounds(::Type{Bool}, w::WignerMatrixBatch{IT}, iᵣ::Integer, m′::IT, m::IT) where {IT} =
    1 ≤ iᵣ ≤ w.Nᵣ && inrange(m′, w.m′ₘᵢₙ, w.m′ₘₐₓ) && inrange(m, w.mₘᵢₙ, w.mₘₐₓ)

@propagate_inbounds function Base.getindex(
    w::WignerMatrixBatch{IT}, iᵣ::Integer, m′::IT, m::IT
) where {IT<:IntegerHalf}
    @boundscheck checkbounds(Bool, w, iᵣ, m′, m) || throw(BoundsError(w, (iᵣ, m′, m)))
    @boundscheck check_storage(w)
    @inbounds parent(w)[
        storage_index(w, Int(iᵣ), Int(m′ - w.m′ₘᵢₙ) + 1, Int(m - w.mₘᵢₙ) + 1)
    ]
end

@propagate_inbounds function Base.setindex!(
    w::WignerMatrixBatch{IT}, v, iᵣ::Integer, m′::IT, m::IT
) where {IT<:IntegerHalf}
    @boundscheck checkbounds(Bool, w, iᵣ, m′, m) || throw(BoundsError(w, (iᵣ, m′, m)))
    @boundscheck check_storage(w)
    @inbounds parent(w)[
        storage_index(w, Int(iᵣ), Int(m′ - w.m′ₘᵢₙ) + 1, Int(m - w.mₘᵢₙ) + 1)
    ] = v
end

# See the note on indexing with indices of other types above.
@propagate_inbounds Base.getindex(
    w::WignerMatrixBatch{IT}, iᵣ::Integer, m′::IndexType, m::IndexType
) where {IT} = w[iᵣ, checked_index(IT, m′, w, "m′"), checked_index(IT, m, w, "m")]
@propagate_inbounds Base.setindex!(
    w::WignerMatrixBatch{IT}, v, iᵣ::Integer, m′::IndexType, m::IndexType
) where {IT} = (w[iᵣ, checked_index(IT, m′, w, "m′"), checked_index(IT, m, w, "m")] = v)

"""
    w[iᵣ]

The [`WignerMatrix`](@ref) of rotor `iᵣ` in a [`WignerMatrixBatch`](@ref), as a view; index
it naturally as `w[iᵣ][m′, m]`.
"""
@propagate_inbounds function Base.getindex(w::WignerMatrixBatch{IT, NT}, iᵣ::Integer) where {IT, NT}
    @boundscheck if !(1 ≤ iᵣ ≤ w.Nᵣ)
        throw(BoundsError(w, (iᵣ,)))
    end
    @boundscheck check_storage(w)
    let p = leading_slice(w, Int(iᵣ))
        WignerMatrix{IT, NT, typeof(p)}(p, w.ℓ, w.m′ₘₐₓ, w.m′ₘᵢₙ, w.mₘₐₓ, w.mₘᵢₙ)
    end
end

rewrap(w::WignerMatrixBatch{IT}, p) where {IT} =
    WignerMatrixBatch{IT, eltype(p), typeof(p)}(
        p, w.ℓ, w.m′ₘₐₓ, w.m′ₘᵢₙ, w.mₘₐₓ, w.mₘᵢₙ, w.Nᵣ
    )

Base.Matrix(w::WignerMatrixBatch) = throw(ArgumentError(
    "A WignerMatrixBatch is 3-dimensional; use `Array(w)`, or `Matrix(w[iᵣ])`."
))


"""
    DegreeBlock{IT, NT, ST}

One harmonic degree's worth of values, indexed naturally by the order ``m``: `v[m]` for
``mₘᵢₙ ≤ m ≤ mₘₐₓ``.  This is the 1-dimensional sibling of [`WignerMatrix`](@ref).
- `IT` is the index type, `Int` or [`HalfOddInteger`](@ref).
- `NT` is the number type of the elements.
- `ST` is the type of the storage, an `AbstractVector{NT}` (see [`AbstractBlock`](@ref)).

The name is deliberately neutral, because the same shape serves three things: a block of
spin-weighted harmonics at one ``ℓ``, from [`sYlm`](@ref) or [`sYlmCalculator`](@ref)'s
`ₛYₗ[s, :]`; one ``ℓ``'s worth of mode weights, from [`ModeWeights`](@ref)'s `w[ℓ, :]`; and
one spin weight's row of a [`SpinMatrix`](@ref), from `b[s, :]`.  All three are values at a
fixed degree indexed by order, whatever they mean.

The elements are in order of increasing ``m``.  The constructor

    DegreeBlock(parent, ℓ; mₘₐₓ=ℓ, mₘᵢₙ=-mₘₐₓ)

wraps `parent` without copying it.  `parent` is a 1-based array of any shape with at least
`mₘₐₓ-mₘᵢₙ+1` entries, read in its linear order, so that a vector of exactly that length is
the block.  The limits must satisfy `-ℓ ≤ mₘᵢₙ ≤ mₘₐₓ ≤ ℓ`; unlike the limits of a
[`WignerMatrix`](@ref), these need not bracket zero, because no recurrence fills the block.
`ℓ` and the limits must all be integers of type `Int`, or all half-odd-integers, each a
[`HalfOddInteger`](@ref) or a `Rational{Int}` with denominator 2.  The keywords may also be
spelled `m_max` and `m_min`.

See also [`DegreeBlockBatch`](@ref).
"""
struct DegreeBlock{IT, NT, ST<:AbstractVector{NT}} <: AbstractBlock{IT, NT, ST}
    parent::ST
    ℓ::IT
    mₘₐₓ::IT
    mₘᵢₙ::IT
    # As for `WignerMatrix`: the parent must be 1-based, and must reach every element.
    function DegreeBlock{IT, NT, ST}(parent, ℓ, mₘₐₓ, mₘᵢₙ) where {IT, NT, ST}
        Base.require_one_based_indexing(parent)
        check_extent(parent, ((mₘₐₓ, mₘᵢₙ, "m"),))
        new{IT, NT, ST}(parent, ℓ, mₘₐₓ, mₘᵢₙ)
    end
end

@index_methods function DegreeBlock(
    parent::ST, ℓ::IT;
    m_max::IndexType=ℓ, mₘₐₓ::IndexType=m_max,
    m_min::IndexType=-mₘₐₓ, mₘᵢₙ::IndexType=m_min
) where {IT<:IndexType, NT, ST<:AbstractArray{NT}}
    validate_m_range(ℓ, mₘₐₓ, mₘᵢₙ)
    Base.require_one_based_indexing(parent)
    p = vec(parent)
    DegreeBlock{IT, NT, typeof(p)}(p, ℓ, mₘₐₓ, mₘᵢₙ)
end

Base.firstindex(v::DegreeBlock) = v.mₘᵢₙ
Base.lastindex(v::DegreeBlock) = v.mₘₐₓ
Base.keys(v::DegreeBlock) = v.mₘᵢₙ:v.mₘₐₓ

Base.checkbounds(::Type{Bool}, v::DegreeBlock{IT}, m::IT) where {IT} =
    inrange(m, v.mₘᵢₙ, v.mₘₐₓ)

@propagate_inbounds function Base.getindex(v::DegreeBlock{IT}, m::IT) where {IT<:IntegerHalf}
    @boundscheck checkbounds(Bool, v, m) || throw(BoundsError(v, m))
    @boundscheck check_storage(v)
    @inbounds parent(v)[storage_index(v, Int(m - v.mₘᵢₙ) + 1)]
end
@propagate_inbounds function Base.setindex!(
    v::DegreeBlock{IT}, x, m::IT
) where {IT<:IntegerHalf}
    @boundscheck checkbounds(Bool, v, m) || throw(BoundsError(v, m))
    @boundscheck check_storage(v)
    @inbounds parent(v)[storage_index(v, Int(m - v.mₘᵢₙ) + 1)] = x
end

# See the note on indexing with indices of other types above.
@propagate_inbounds Base.getindex(v::DegreeBlock{IT}, m::IndexType) where {IT} =
    v[checked_index(IT, m, v, "m")]
@propagate_inbounds Base.setindex!(v::DegreeBlock{IT}, x, m::IndexType) where {IT} =
    (v[checked_index(IT, m, v, "m")] = x)

Base.Vector(v::DegreeBlock) = Vector(array_view(v))
rewrap(v::DegreeBlock{IT}, p) where {IT} =
    DegreeBlock{IT, eltype(p), typeof(p)}(p, v.ℓ, v.mₘₐₓ, v.mₘᵢₙ)


"""
    DegreeBlockBatch{IT, NT, ST}

`Nᵣ` rows of values of a single ``ℓ``, stored together and indexed as `v[iᵣ, m]`.  This is
the 1-dimensional sibling of [`WignerMatrixBatch`](@ref), and is what a batched
[`sYlmCalculator`](@ref) built for one spin weight yields, for either kind of index.
`v[iᵣ]` gives the [`DegreeBlock`](@ref) view of one rotor's row.
- `IT` is the index type, `Int` or [`HalfOddInteger`](@ref).
- `NT` is the number type of the elements.
- `ST` is the type of the storage, an `AbstractVector{NT}`, or an `AbstractMatrix{NT}` for
  the slice `b[:, s, :]` of a [`SpinMatrixBatch`](@ref) (see [`AbstractBlock`](@ref)).

The elements are ordered `[iᵣ, m]`.  The constructor, `DegreeBlockBatch(parent, ℓ; mₘₐₓ=ℓ,
mₘᵢₙ=-mₘₐₓ)`, takes a 1-based array of any shape, whose first extent is `Nᵣ` and which has
at least as many entries as the block, read in its linear order, so that a matrix of exactly
the block's shape is the block; it takes the indices and the limits exactly as
[`DegreeBlock`](@ref) does, and under the same rules.
"""
struct DegreeBlockBatch{IT, NT, ST<:Union{AbstractVector{NT}, AbstractMatrix{NT}}} <: AbstractBlock{IT, NT, ST}
    # The storage is a matrix only for a slice `b[:, s, :]` of a `SpinMatrixBatch`; see
    # "Storage" above.
    parent::ST
    ℓ::IT
    mₘₐₓ::IT
    mₘᵢₙ::IT
    Nᵣ::Int
    # As for `WignerMatrix`: the parent must be 1-based, and must reach every element.
    function DegreeBlockBatch{IT, NT, ST}(parent, ℓ, mₘₐₓ, mₘᵢₙ, Nᵣ) where {IT, NT, ST}
        Base.require_one_based_indexing(parent)
        check_extent(parent, ((Nᵣ,), (mₘₐₓ, mₘᵢₙ, "m")))
        new{IT, NT, ST}(parent, ℓ, mₘₐₓ, mₘᵢₙ, Nᵣ)
    end
end

@index_methods function DegreeBlockBatch(
    parent::ST, ℓ::IT;
    m_max::IndexType=ℓ, mₘₐₓ::IndexType=m_max,
    m_min::IndexType=-mₘₐₓ, mₘᵢₙ::IndexType=m_min
) where {IT<:IndexType, NT, ST<:AbstractArray{NT}}
    validate_m_range(ℓ, mₘₐₓ, mₘᵢₙ)
    Base.require_one_based_indexing(parent)
    Nᵣ = size(parent, 1)
    p = vec(parent)
    DegreeBlockBatch{IT, NT, typeof(p)}(p, ℓ, mₘₐₓ, mₘᵢₙ, Nᵣ)
end

Nᵣ(v::DegreeBlockBatch) = v.Nᵣ

Base.checkbounds(::Type{Bool}, v::DegreeBlockBatch{IT}, iᵣ::Integer, m::IT) where {IT} =
    1 ≤ iᵣ ≤ v.Nᵣ && inrange(m, v.mₘᵢₙ, v.mₘₐₓ)

@propagate_inbounds function Base.getindex(
    v::DegreeBlockBatch{IT}, iᵣ::Integer, m::IT
) where {IT<:IntegerHalf}
    @boundscheck checkbounds(Bool, v, iᵣ, m) || throw(BoundsError(v, (iᵣ, m)))
    @boundscheck check_storage(v)
    @inbounds parent(v)[storage_index(v, Int(iᵣ), Int(m - v.mₘᵢₙ) + 1)]
end
@propagate_inbounds function Base.setindex!(
    v::DegreeBlockBatch{IT}, x, iᵣ::Integer, m::IT
) where {IT<:IntegerHalf}
    @boundscheck checkbounds(Bool, v, iᵣ, m) || throw(BoundsError(v, (iᵣ, m)))
    @boundscheck check_storage(v)
    @inbounds parent(v)[storage_index(v, Int(iᵣ), Int(m - v.mₘᵢₙ) + 1)] = x
end

# See the note on indexing with indices of other types above.
@propagate_inbounds Base.getindex(v::DegreeBlockBatch{IT}, iᵣ::Integer, m::IndexType) where {IT} =
    v[iᵣ, checked_index(IT, m, v, "m")]
@propagate_inbounds Base.setindex!(
    v::DegreeBlockBatch{IT}, x, iᵣ::Integer, m::IndexType
) where {IT} = (v[iᵣ, checked_index(IT, m, v, "m")] = x)

@propagate_inbounds function Base.getindex(v::DegreeBlockBatch{IT, NT}, iᵣ::Integer) where {IT, NT}
    @boundscheck if !(1 ≤ iᵣ ≤ v.Nᵣ)
        throw(BoundsError(v, (iᵣ,)))
    end
    @boundscheck check_storage(v)
    let p = leading_slice(v, Int(iᵣ))
        DegreeBlock{IT, NT, typeof(p)}(p, v.ℓ, v.mₘₐₓ, v.mₘᵢₙ)
    end
end

Base.Matrix(v::DegreeBlockBatch) = Matrix(array_view(v))
rewrap(v::DegreeBlockBatch{IT}, p) where {IT} =
    DegreeBlockBatch{IT, eltype(p), typeof(p)}(p, v.ℓ, v.mₘₐₓ, v.mₘᵢₙ, v.Nᵣ)


"""
    SpinMatrix{IT, NT, ST}

The values of a single ``ℓ`` for a range of spin weights, indexed naturally by ``(s, m)``:
`b[s, m]` for ``sₘᵢₙ ≤ s ≤ sₘₐₓ`` and ``mₘᵢₙ ≤ m ≤ mₘₐₓ``, and `b[s, :]` for one whole row
as a [`DegreeBlock`](@ref).  This is what an [`sYlmCalculator`](@ref) built for a range of
spin weights yields for each ``ℓ``.
- `IT` is the index type, `Int` or [`HalfOddInteger`](@ref).
- `NT` is the number type of the elements.
- `ST` is the type of the storage, an `AbstractVector{NT}` (see [`AbstractBlock`](@ref)).

The elements are ordered `[s, m]`.  The constructor

    SpinMatrix(parent, ℓ; sₘₐₓ, sₘᵢₙ, mₘₐₓ=ℓ, mₘᵢₙ=-mₘₐₓ)

wraps `parent` without copying it.  `parent` is a 1-based array of any shape with at least
as many entries as the block, read in its linear order, so that a matrix of exactly the
block's shape is the block.  The spin limits are required, and the keywords may also be
spelled `s_max`, `s_min`, `m_max`, and `m_min`.  `ℓ` and every limit must all be integers of
type `Int`, or all half-odd-integers, each a [`HalfOddInteger`](@ref) or a `Rational{Int}`
with denominator 2.  The ``m`` axis obeys the rule of [`DegreeBlock`](@ref), `-ℓ ≤ mₘᵢₙ ≤
mₘₐₓ ≤ ℓ`.

The spin axis, [`spins`](@ref)`(b)`, is under none of the restrictions a
[`WignerMatrix`](@ref) places on its ``m′``: it is whatever range of spin weights was asked
for, so it may lie wholly on one side of zero, and it may reach beyond ``ℓ`` — the harmonics
with ``|s| > ℓ`` simply vanish, and a calculator stores them as zeros.  It must only be in
order, `sₘᵢₙ ≤ sₘₐₓ`.

See also [`SpinMatrixBatch`](@ref) and [`DegreeBlock`](@ref).
"""
struct SpinMatrix{IT, NT, ST<:AbstractVector{NT}} <: AbstractBlock{IT, NT, ST}
    parent::ST
    ℓ::IT
    sₘₐₓ::IT
    sₘᵢₙ::IT
    mₘₐₓ::IT
    mₘᵢₙ::IT
    # As for `WignerMatrix`: the parent must be 1-based, and must reach every element.
    function SpinMatrix{IT, NT, ST}(parent, ℓ, sₘₐₓ, sₘᵢₙ, mₘₐₓ, mₘᵢₙ) where {IT, NT, ST}
        Base.require_one_based_indexing(parent)
        check_extent(parent, ((sₘₐₓ, sₘᵢₙ, "s"), (mₘₐₓ, mₘᵢₙ, "m")))
        new{IT, NT, ST}(parent, ℓ, sₘₐₓ, sₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    end
end

# The spin limits are required, but each has two spellings, so neither spelling can be a
# required keyword; the absence of both is refused as Julia would refuse a missing keyword.
@index_methods function SpinMatrix(
    parent::ST, ℓ::IT;
    s_max::Union{Nothing, IndexType}=nothing, sₘₐₓ::Union{Nothing, IndexType}=s_max,
    s_min::Union{Nothing, IndexType}=nothing, sₘᵢₙ::Union{Nothing, IndexType}=s_min,
    m_max::IndexType=ℓ, mₘₐₓ::IndexType=m_max,
    m_min::IndexType=-mₘₐₓ, mₘᵢₙ::IndexType=m_min
) where {IT<:IndexType, NT, ST<:AbstractArray{NT}}
    sₘₐₓ === nothing && throw(UndefKeywordError(:sₘₐₓ))
    sₘᵢₙ === nothing && throw(UndefKeywordError(:sₘᵢₙ))
    validate_s_range(sₘₐₓ, sₘᵢₙ)
    validate_m_range(ℓ, mₘₐₓ, mₘᵢₙ)
    Base.require_one_based_indexing(parent)
    p = vec(parent)
    SpinMatrix{IT, NT, typeof(p)}(p, ℓ, sₘₐₓ, sₘᵢₙ, mₘₐₓ, mₘᵢₙ)
end

sₘₐₓ(b::SpinMatrix) = b.sₘₐₓ
sₘᵢₙ(b::SpinMatrix) = b.sₘᵢₙ
# The spin axis is `spins(b)`, as it is for the calculator that produced the block.  It is
# deliberately not `keys(b)`: iteration visits every (s, m) element, and `Base` builds
# `pairs`, and with it `findmax`, `argmax`, `findall` and `findfirst`, by zipping `keys`
# with the iteration, so a `keys` that named only the spin weights would pair each of them
# with an element and report it as that element's position.  Without `keys` those functions
# are a `MethodError`, as they are for a `WignerMatrix`.
spins(b::SpinMatrix) = b.sₘᵢₙ:b.sₘₐₓ

Base.checkbounds(::Type{Bool}, b::SpinMatrix{IT}, s::IT, m::IT) where {IT} =
    inrange(s, b.sₘᵢₙ, b.sₘₐₓ) && inrange(m, b.mₘᵢₙ, b.mₘₐₓ)

@propagate_inbounds function Base.getindex(
    b::SpinMatrix{IT}, s::IT, m::IT
) where {IT<:IntegerHalf}
    @boundscheck checkbounds(Bool, b, s, m) || throw(BoundsError(b, (s, m)))
    @boundscheck check_storage(b)
    @inbounds parent(b)[storage_index(b, Int(s - b.sₘᵢₙ) + 1, Int(m - b.mₘᵢₙ) + 1)]
end
@propagate_inbounds function Base.setindex!(
    b::SpinMatrix{IT}, x, s::IT, m::IT
) where {IT<:IntegerHalf}
    @boundscheck checkbounds(Bool, b, s, m) || throw(BoundsError(b, (s, m)))
    @boundscheck check_storage(b)
    @inbounds parent(b)[storage_index(b, Int(s - b.sₘᵢₙ) + 1, Int(m - b.mₘᵢₙ) + 1)] = x
end

# See the note on indexing with indices of other types above.
@propagate_inbounds Base.getindex(
    b::SpinMatrix{IT}, s::IndexType, m::IndexType
) where {IT} = b[checked_index(IT, s, b, "s"), checked_index(IT, m, b, "m")]
@propagate_inbounds Base.setindex!(
    b::SpinMatrix{IT}, x, s::IndexType, m::IndexType
) where {IT} = (b[checked_index(IT, s, b, "s"), checked_index(IT, m, b, "m")] = x)

"""
    b[s, :]

The row of one spin weight of a [`SpinMatrix`](@ref), as a [`DegreeBlock`](@ref) view; it is
then indexed naturally as `b[s, :][m]`.  The notation follows [`ModeWeights`](@ref)'s `w[ℓ,
:]`, so that a loop over spin weights reads the same in both places.
"""
@propagate_inbounds function Base.getindex(
    b::SpinMatrix{IT, NT}, s::IT, ::Colon
) where {IT<:IntegerHalf, NT}
    @boundscheck inrange(s, b.sₘᵢₙ, b.sₘₐₓ) || throw(BoundsError(b, (s, :)))
    @boundscheck check_storage(b)
    let p = leading_slice(b, Int(s - b.sₘᵢₙ) + 1)
        DegreeBlock{IT, NT, typeof(p)}(p, b.ℓ, b.mₘₐₓ, b.mₘᵢₙ)
    end
end
@propagate_inbounds Base.getindex(b::SpinMatrix{IT}, s::IndexType, ::Colon) where {IT} =
    b[checked_index(IT, s, b, "s"), :]

Base.Matrix(b::SpinMatrix) = Matrix(array_view(b))
rewrap(b::SpinMatrix{IT}, p) where {IT} =
    SpinMatrix{IT, eltype(p), typeof(p)}(p, b.ℓ, b.sₘₐₓ, b.sₘᵢₙ, b.mₘₐₓ, b.mₘᵢₙ)


"""
    SpinMatrixBatch{IT, NT, ST}

`Nᵣ` [`SpinMatrix`](@ref) blocks of one ``ℓ``, stored together and indexed as `b[iᵣ, s, m]`.
This is what a batched [`sYlmCalculator`](@ref) yields when it was built for a range of spin
weights.  `b[iᵣ]` gives the `SpinMatrix` view of one rotor's block, which is then indexed
naturally as `b[iᵣ][s, m]`.
- `IT` is the index type, `Int` or [`HalfOddInteger`](@ref).
- `NT` is the number type of the elements.
- `ST` is the type of the storage, an `AbstractVector{NT}` (see [`AbstractBlock`](@ref)).

The elements are ordered `[iᵣ, s, m]`, exactly as in the calculator.  The constructor,
`SpinMatrixBatch(parent, ℓ; sₘₐₓ, sₘᵢₙ, mₘₐₓ=ℓ, mₘᵢₙ=-mₘₐₓ)`, takes a 1-based array of any
shape, whose first extent is `Nᵣ` and which has at least as many entries as the block, read
in its linear order, so that an array of exactly the block's shape is the block; it takes
the indices and the limits exactly as [`SpinMatrix`](@ref) does, and under the same rules.

See also [`SpinMatrix`](@ref) and [`DegreeBlockBatch`](@ref).
"""
struct SpinMatrixBatch{IT, NT, ST<:AbstractVector{NT}} <: AbstractBlock{IT, NT, ST}
    parent::ST
    ℓ::IT
    sₘₐₓ::IT
    sₘᵢₙ::IT
    mₘₐₓ::IT
    mₘᵢₙ::IT
    Nᵣ::Int
    # As for `WignerMatrix`: the parent must be 1-based, and must reach every element.
    function SpinMatrixBatch{IT, NT, ST}(parent, ℓ, sₘₐₓ, sₘᵢₙ, mₘₐₓ, mₘᵢₙ, Nᵣ) where {IT, NT, ST}
        Base.require_one_based_indexing(parent)
        check_extent(parent, ((Nᵣ,), (sₘₐₓ, sₘᵢₙ, "s"), (mₘₐₓ, mₘᵢₙ, "m")))
        new{IT, NT, ST}(parent, ℓ, sₘₐₓ, sₘᵢₙ, mₘₐₓ, mₘᵢₙ, Nᵣ)
    end
end

# As for `SpinMatrix`, the spin limits are required in either spelling.
@index_methods function SpinMatrixBatch(
    parent::ST, ℓ::IT;
    s_max::Union{Nothing, IndexType}=nothing, sₘₐₓ::Union{Nothing, IndexType}=s_max,
    s_min::Union{Nothing, IndexType}=nothing, sₘᵢₙ::Union{Nothing, IndexType}=s_min,
    m_max::IndexType=ℓ, mₘₐₓ::IndexType=m_max,
    m_min::IndexType=-mₘₐₓ, mₘᵢₙ::IndexType=m_min
) where {IT<:IndexType, NT, ST<:AbstractArray{NT}}
    sₘₐₓ === nothing && throw(UndefKeywordError(:sₘₐₓ))
    sₘᵢₙ === nothing && throw(UndefKeywordError(:sₘᵢₙ))
    validate_s_range(sₘₐₓ, sₘᵢₙ)
    validate_m_range(ℓ, mₘₐₓ, mₘᵢₙ)
    Base.require_one_based_indexing(parent)
    Nᵣ = size(parent, 1)
    p = vec(parent)
    SpinMatrixBatch{IT, NT, typeof(p)}(p, ℓ, sₘₐₓ, sₘᵢₙ, mₘₐₓ, mₘᵢₙ, Nᵣ)
end

sₘₐₓ(b::SpinMatrixBatch) = b.sₘₐₓ
sₘᵢₙ(b::SpinMatrixBatch) = b.sₘᵢₙ
Nᵣ(b::SpinMatrixBatch) = b.Nᵣ
spins(b::SpinMatrixBatch) = b.sₘᵢₙ:b.sₘₐₓ

Base.checkbounds(::Type{Bool}, b::SpinMatrixBatch{IT}, iᵣ::Integer, s::IT, m::IT) where {IT} =
    1 ≤ iᵣ ≤ b.Nᵣ && inrange(s, b.sₘᵢₙ, b.sₘₐₓ) && inrange(m, b.mₘᵢₙ, b.mₘₐₓ)

@propagate_inbounds function Base.getindex(
    b::SpinMatrixBatch{IT}, iᵣ::Integer, s::IT, m::IT
) where {IT<:IntegerHalf}
    @boundscheck checkbounds(Bool, b, iᵣ, s, m) || throw(BoundsError(b, (iᵣ, s, m)))
    @boundscheck check_storage(b)
    @inbounds parent(b)[
        storage_index(b, Int(iᵣ), Int(s - b.sₘᵢₙ) + 1, Int(m - b.mₘᵢₙ) + 1)
    ]
end
@propagate_inbounds function Base.setindex!(
    b::SpinMatrixBatch{IT}, x, iᵣ::Integer, s::IT, m::IT
) where {IT<:IntegerHalf}
    @boundscheck checkbounds(Bool, b, iᵣ, s, m) || throw(BoundsError(b, (iᵣ, s, m)))
    @boundscheck check_storage(b)
    @inbounds parent(b)[
        storage_index(b, Int(iᵣ), Int(s - b.sₘᵢₙ) + 1, Int(m - b.mₘᵢₙ) + 1)
    ] = x
end

# See the note on indexing with indices of other types above.
@propagate_inbounds Base.getindex(
    b::SpinMatrixBatch{IT}, iᵣ::Integer, s::IndexType, m::IndexType
) where {IT} = b[iᵣ, checked_index(IT, s, b, "s"), checked_index(IT, m, b, "m")]
@propagate_inbounds Base.setindex!(
    b::SpinMatrixBatch{IT}, x, iᵣ::Integer, s::IndexType, m::IndexType
) where {IT} = (b[iᵣ, checked_index(IT, s, b, "s"), checked_index(IT, m, b, "m")] = x)

"""
    b[:, s, :]

The rows of one spin weight of a [`SpinMatrixBatch`](@ref), over every rotor, as a
[`DegreeBlockBatch`](@ref) view; it is then indexed naturally as `b[:, s, :][iᵣ, m]`.  The
notation follows [`SpinMatrix`](@ref)'s `b[s, :]`.
"""
@propagate_inbounds function Base.getindex(
    b::SpinMatrixBatch{IT, NT}, ::Colon, s::IT, ::Colon
) where {IT<:IntegerHalf, NT}
    @boundscheck inrange(s, b.sₘᵢₙ, b.sₘₐₓ) || throw(BoundsError(b, (:, s, :)))
    check_storage(b)
    let p = view(elements(b), :, Int(s - b.sₘᵢₙ) + 1, :)
        DegreeBlockBatch{IT, NT, typeof(p)}(p, b.ℓ, b.mₘₐₓ, b.mₘᵢₙ, b.Nᵣ)
    end
end
@propagate_inbounds Base.getindex(
    b::SpinMatrixBatch{IT}, ::Colon, s::IndexType, ::Colon
) where {IT} = b[:, checked_index(IT, s, b, "s"), :]

"""
    b[iᵣ]

The [`SpinMatrix`](@ref) of rotor `iᵣ` in a [`SpinMatrixBatch`](@ref), as a view; index it
naturally as `b[iᵣ][s, m]`.
"""
@propagate_inbounds function Base.getindex(b::SpinMatrixBatch{IT, NT}, iᵣ::Integer) where {IT, NT}
    @boundscheck if !(1 ≤ iᵣ ≤ b.Nᵣ)
        throw(BoundsError(b, (iᵣ,)))
    end
    @boundscheck check_storage(b)
    let p = leading_slice(b, Int(iᵣ))
        SpinMatrix{IT, NT, typeof(p)}(p, b.ℓ, b.sₘₐₓ, b.sₘᵢₙ, b.mₘₐₓ, b.mₘᵢₙ)
    end
end

rewrap(b::SpinMatrixBatch{IT}, p) where {IT} =
    SpinMatrixBatch{IT, eltype(p), typeof(p)}(p, b.ℓ, b.sₘₐₓ, b.sₘᵢₙ, b.mₘₐₓ, b.mₘᵢₙ, b.Nᵣ)


### The shared array interface of the blocks.
#
# Every block answers two questions about itself, and the whole array interface is written
# once in terms of the answers: `axis_roles` names its axes, in the order they are indexed,
# and `natural_axes` gives the matching ranges.  A third, `rewrap`, defined beside each
# block, puts its labels on other storage, which is what `copy`, `similar` and `relabel`
# need.
#
# `getindex` and `setindex!` are deliberately *not* written this way, and neither are the
# accessors `ℓ`, `mₘₐₓ`, `sₘᵢₙ` and the rest.  Those are the hot path: the bounds check
# there reads the stored limits directly rather than building a range to test membership in,
# which is the difference the note above `inrange` measures at 145 ns against 0.9 ns per
# element.

# The roles are a property of the type, so these fold away at compile time.  A batched
# block's leading axis is an ordinary 1-based rotor position, not a natural index, which is
# why it is the one role whose range is a plain `UnitRange` rather than a `WignerRange`.
axis_roles(::Type{<:WignerMatrix}) = (:m′, :m)
axis_roles(::Type{<:WignerMatrixBatch}) = (:iᵣ, :m′, :m)
axis_roles(::Type{<:DegreeBlock}) = (:m,)
axis_roles(::Type{<:DegreeBlockBatch}) = (:iᵣ, :m)
axis_roles(::Type{<:SpinMatrix}) = (:s, :m)
axis_roles(::Type{<:SpinMatrixBatch}) = (:iᵣ, :s, :m)
axis_roles(w::AbstractBlock) = axis_roles(typeof(w))

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

# The extent is taken as the difference of the endpoints rather than as `length` of the
# range, because that is an `Int` by construction for a half-odd-integer axis, where
# `Int(ℓ)` throws.
@inline axis_extent(r) = Int(last(r) - first(r)) + 1

@inline Base.axes(w::AbstractBlock) = natural_axes(w)
@inline Base.size(w::AbstractBlock) = map(axis_extent, natural_axes(w))
Base.length(w::AbstractBlock) = prod(size(w))
Base.ndims(w::AbstractBlock) = length(axis_roles(w))
# The form for a type names the six, rather than `AbstractBlock`, because a block type with
# free parameters, such as `WignerMatrix` or `WignerMatrix{Int}`, is not a subtype of
# `AbstractBlock`: its parameters are bounded only by that supertype, not by the type
# itself.
Base.ndims(::Type{T}) where {T<:Union{
    WignerMatrix, WignerMatrixBatch,
    DegreeBlock, DegreeBlockBatch,
    SpinMatrix, SpinMatrixBatch,
}} = length(axis_roles(T))
# Trailing dimensions behave as they do for `AbstractArray`, so that generic code written
# against a plain array works unchanged: `axes(w, d)` is `OneTo(1)` and `size(w, d)` is `1`.
Base.axes(w::AbstractBlock, d::Integer) = d ≤ ndims(w) ? axes(w)[d] : Base.OneTo(1)
Base.size(w::AbstractBlock, d::Integer) = d ≤ ndims(w) ? size(w)[d] : 1
# The throwing form, as for an array.
function Base.checkbounds(w::AbstractBlock, I...)
    checkbounds(Bool, w, I...) || throw(BoundsError(w, I))
    nothing
end

# A block is a batch when its leading axis is the rotor axis; the roles are constants of the
# type, so this is one too.
isbatched(w::AbstractBlock) = first(axis_roles(w)) === :iᵣ

# Only a `WignerMatrix` and a `WignerMatrixBatch` have an m′ axis.  The other blocks are
# told so, rather than being shown an error about a missing field.
m′ₘₐₓ(w::Union{WignerMatrix, WignerMatrixBatch}) = w.m′ₘₐₓ
m′ₘᵢₙ(w::Union{WignerMatrix, WignerMatrixBatch}) = w.m′ₘᵢₙ
m′ₘₐₓ(w::AbstractBlock) = throw(no_m′_axis(w))
m′ₘᵢₙ(w::AbstractBlock) = throw(no_m′_axis(w))
@noinline no_m′_axis(w) = ArgumentError(
    "A `$(nameof(typeof(w)))` has no m′ axis; its axes are "
    * join(("`$r`" for r ∈ axis_roles(w)), ", ") * ", and the limits of the natural ones "
    * "are read with `mₘₐₓ` and `mₘᵢₙ`, and for a spin axis `sₘₐₓ` and `sₘᵢₙ`."
)

# Everything below reads the elements through `array_view(w)`, which is an array of the
# block's shape (see "Storage" above), in the order of `Array(w)`.
Base.Array(w::AbstractBlock) = Array(array_view(w))
Base.collect(w::AbstractBlock) = Array(w)

# Iteration visits the elements in the order of `Array(w)`, entry by entry of the storage,
# so that a loop over a block, and the functions that fold over its elements one by one,
# such as `foldl`, `any`, and `all`, see them as they would see the elements of an ordinary
# array.  (The matrix of a slice `b[:, s, :]` is read in its linear order, which is that of
# `Array(w)`.)
@inline function Base.iterate(w::AbstractBlock, i::Int=1)
    check_storage(w)
    i > length(w) && return nothing
    (@inbounds parent(w)[i], i + 1)
end

# Reductions over a block — `sum`, `prod`, `maximum`, `count`, and the others that `Base`
# builds on `mapreduce` — run over its array, `array_view(w)`, rather than through the
# iteration, and give exactly what they give for that array.  `Base` reduces a linearly
# indexed array pairwise and vectorized, so a block over contiguous storage, such as a
# calculator's block or one built on an array, reduces exactly as `Array(w)` does; over the
# strided storage of a slice, the terms may be associated differently.  Only `init` is
# passed on: `dims` is refused, as it is for any collection that is not an array.  `f` and
# `op` are type parameters because Julia does not otherwise specialize a method on functions
# that it only passes on, and the call to the array's reduction would then be dynamic and
# would allocate.
Base.mapreduce(f::F, op::OP, w::AbstractBlock; init=nothing) where {F, OP} =
    init === nothing ? mapreduce(f, op, array_view(w)) : mapreduce(f, op, array_view(w); init)

# The copy is of the block alone, as a `Vector` of its elements.
Base.copy(w::AbstractBlock) = rewrap(w, collect(vec(array_view(w))))

"""
    similar(w::AbstractBlock, [T=eltype(w)])

A new block of the same kind as `w`, with the same ℓ and the same natural axes (and so, for
a batch, the same `Nᵣ`), over uninitialized storage of element type `T`.  The storage is a
`Vector` of exactly `length(w)` entries, even when `parent(w)` is a calculator's buffer or a
view.
"""
Base.similar(w::AbstractBlock, ::Type{T}=eltype(w)) where {T} =
    rewrap(w, Vector{T}(undef, length(w)))

function Base.summary(io::IO, w::AbstractBlock{IT, NT}) where {IT, NT}
    print(io, axes_string(axes(w)), " ", nameof(typeof(w)), "{", IT, ", ", NT, "} for ℓ=", ℓ(w))
end
Base.show(io::IO, w::AbstractBlock) = summary(io, w)
# The elements are printed from `array_view(w)` rather than from `Array(w)`, because a view
# does not read them: uninitialized storage of a non-bits type such as `BigFloat` then
# prints as `#undef`, where `Array(w)` would throw an `UndefRefError`.
function Base.show(io::IO, ::MIME"text/plain", w::AbstractBlock)
    summary(io, w)
    println(io, ":")
    Base.print_array(io, array_view(w))
end

# The labels of a block are its kind — the roles of its axes — its ℓ, and the limits of its
# axes, including the number of rotors of a batch.  The same numbers under different labels
# are not the same block, so `==`, `isequal`, `≈`, and `hash` count the labels as well as
# the numbers, which they read through `array_view` without copying them.  Blocks of the
# same labels but different number types compare their numbers as arrays do, and hash as
# arrays do, so that `hash` agrees with `isequal` between them too.
block_labels(w::AbstractBlock) =
    (axis_roles(w), ℓ(w), map(a -> (first(a), last(a)), natural_axes(w)))
Base.:(==)(a::AbstractBlock, b::AbstractBlock) =
    block_labels(a) == block_labels(b) && array_view(a) == array_view(b)
Base.isequal(a::AbstractBlock, b::AbstractBlock) =
    isequal(block_labels(a), block_labels(b)) && isequal(array_view(a), array_view(b))
Base.isapprox(a::AbstractBlock, b::AbstractBlock; kwargs...) =
    block_labels(a) == block_labels(b) && isapprox(array_view(a), array_view(b); kwargs...)
Base.hash(w::AbstractBlock, h::UInt) =
    hash(array_view(w), hash(block_labels(w), hash(:AbstractBlock, h)))
