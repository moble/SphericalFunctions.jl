import Base: @propagate_inbounds

"""
    AbstractWignerMatrix{IT, NT, ST}

Abstract supertype of the containers that hold the values of one ``ℓ``, indexed by their
natural indices ``m′``, ``m`` and ``s`` rather than by position.
- `IT` is the index type, `Int` or [`HalfOddInteger`](@ref), which is the type of ``ℓ`` and
  of every natural index.
- `NT` is the number type, such as `ComplexF64` for a block of ``𝔇`` or `Float64` for one
  of ``d``.
- `ST` is the type of the storage, which is 1-based.

These types are *not* `AbstractArray`s, because half-integer indices cannot satisfy that
interface: `axes` must be integer ranges there, and `-3//2:3//2` is not one.  Consequently
linear algebra (`*`, `'`, `\\`, `lu`) does not apply to them; [`array_view`](@ref) gives the
numbers as an ordinary 1-based array, and [`relabel`](@ref) puts the labels back.

# Methods

Every subtype has `parent(w)`, the storage; `ℓ(w)`, the degree; `ℓₘᵢₙ(w)`, the smallest
degree of the index type, 0 or 1//2; `eltype(w)`, the number type; [`ishalfinteger`](@ref);
and `summary` and `show`.  Each subtype defines its own indexing, since the natural indices
differ from one to the next.

The six blocks — [`WignerMatrix`](@ref), indexed `w[m′, m]`, [`WignerMatrixBatch`](@ref),
`w[iᵣ, m′, m]`, [`DegreeBlock`](@ref), `v[m]`, [`DegreeBlockBatch`](@ref), `v[iᵣ, m]`,
[`SpinMatrix`](@ref), `b[s, m]`, and [`SpinMatrixBatch`](@ref), `b[iᵣ, s, m]` — share an
array-like interface as well:
- `axes(w)` are the ranges of the natural indices, one per axis, with `1:Nᵣ` for the rotor
  axis of a batch, and `size(w)`, `length(w)` and `ndims(w)` describe the block, whose
  storage may be larger;
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

The recursion workspaces [`HWedge`](@ref) and [`HAxis`](@ref) are subtypes too, but not
blocks: their storage is triangular, `size` and `length` describe that storage, and they
have none of the array-like interface above beyond `axes`, `ndims` and `size(w, d)`.
"""
abstract type AbstractWignerMatrix{IT<:IntegerHalf, NT, ST<:AbstractArray{NT}} end
# Note that this is deliberately *not* a subtype of `AbstractMatrix`: the natural indices
# `(m′, m)` may be half-odd-integers, which cannot satisfy the `AbstractArray` interface
# (integer `axes`).  The array-like methods that make sense are defined explicitly below.

### General methods for all AbstractWignerMatrix types

Base.parent(w::AbstractWignerMatrix) = w.parent

ℓ(w::AbstractWignerMatrix{IT}) where {IT} = w.ℓ
# The smallest degree of an index type, of an index, and of the blocks, whose degrees are
# all of one index type.  Other arguments have no method, so that `applicable(ℓₘᵢₙ, x)` says
# whether `x` has one.
ℓₘᵢₙ(x::IntegerHalf) = ℓₘᵢₙ(typeof(x))
ℓₘᵢₙ(::Type{IT}) where {IT<:Integer} = zero(IT)
ℓₘᵢₙ(::Type{IT}) where {IT<:HalfOddInteger} = unsafe_half_odd_integer(1)
ℓₘᵢₙ(::AbstractWignerMatrix{IT}) where {IT} = ℓₘᵢₙ(IT)

# The limits of the axes are fields of the types that have those axes; the blocks with a
# single m axis refuse the m′ accessors with an explanation, below, rather than with an
# error about a missing field.
m′ₘₐₓ(w::AbstractWignerMatrix{IT}) where {IT} = w.m′ₘₐₓ
m′ₘᵢₙ(w::AbstractWignerMatrix{IT}) where {IT} = w.m′ₘᵢₙ
mₘₐₓ(w::AbstractWignerMatrix{IT}) where {IT} = w.mₘₐₓ
mₘᵢₙ(w::AbstractWignerMatrix{IT}) where {IT} = w.mₘᵢₙ

ishalfinteger(::AbstractWignerMatrix{IT}) where {IT<:Integer} = false
ishalfinteger(::AbstractWignerMatrix{IT}) where {IT<:HalfOddInteger} = true

Base.eltype(::AbstractWignerMatrix{IT, NT, ST}) where {IT, NT, ST} = NT
Base.eltype(::Type{<:AbstractWignerMatrix{IT, NT, ST}}) where {IT, NT, ST} = NT
# These describe the storage, which is what they mean for the recursion workspaces `HWedge`
# and `HAxis`; the blocks, whose storage may be larger than the block, have their own
# methods in terms of their axes (see "The shared array interface" at the end of this file).
Base.size(w::AbstractWignerMatrix{IT, NT, ST}) where {IT, NT, ST} = size(parent(w))
Base.length(w::AbstractWignerMatrix{IT, NT, ST}) where {IT, NT, ST} = length(parent(w))

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

function Base.axes(w::AbstractWignerMatrix{IT}) where {IT}
    (WignerRange(m′ₘᵢₙ(w):m′ₘₐₓ(w)), WignerRange(mₘᵢₙ(w):mₘₐₓ(w)))
end
# Trailing dimensions behave as they do for `AbstractArray`, so that generic code written
# against a plain array works unchanged here: `axes(w, d)` is `OneTo(1)` and `size(w, d)` is
# `1` for `d > ndims(w)`.  The rank is that of the axes, which is 2 for the generic `axes`
# above; a subtype with axes of another rank, such as `HWedge`, says so for its type too.
# `size(w, d)` compares `d` with the length of `size(w)` itself, because the generic `size`
# describes the storage, whose rank may differ from that of the axes.
Base.axes(w::AbstractWignerMatrix, d::Integer) = d ≤ ndims(w) ? axes(w)[d] : Base.OneTo(1)
Base.size(w::AbstractWignerMatrix, d::Integer) = d ≤ length(size(w)) ? size(w)[d] : 1
Base.ndims(w::AbstractWignerMatrix) = length(axes(w))
Base.ndims(::Type{<:AbstractWignerMatrix}) = 2

function Base.summary(io::IO, w::AbstractWignerMatrix{IT, NT}) where {IT, NT}
    print(io, axes_string(axes(w)), " ", nameof(typeof(w)), "{", IT, ", ", NT, "} for ℓ=", ℓ(w))
end
Base.show(io::IO, w::AbstractWignerMatrix) = summary(io, w)
function Base.show(io::IO, ::MIME"text/plain", w::AbstractWignerMatrix)
    summary(io, w)
    println(io, ":")
    Base.print_array(io, stored_elements(w))
end

# Every container keeps its elements in the leading block of its 1-based storage, in the
# order of `Array(w)`, so this view equals `Array(w)`.  The `show` methods print it instead
# because it does not read the elements: uninitialized storage of a non-bits type such as
# `BigFloat` then prints as `#undef`, where `Array(w)` would throw an `UndefRefError`.
stored_elements(w) = view(parent(w), map(n -> 1:n, size(w))...)



### Validation of the limits
#
# A block with an m′ axis — a `WignerMatrix`, a `WignerMatrixBatch`, and the calculators and
# wedges that fill them — must have its limits in order, within ±ℓₘₐₓ, and bracketing ±ℓₘᵢₙ:
# the recurrence starts from the row m′ = 0 for integer indices and from the pair of rows m′
# = ±1/2 for half-odd-integers, and computes the rest outward from there, so a range that
# misses them cannot be filled.  The same holds for the m axis of those blocks.  The
# one-axis blocks — `DegreeBlock`, `SpinMatrix` and their batches — label storage that no
# recurrence fills, so they check only that their m range is in order and within ±ℓ, in
# `validate_m_range`; their spin axis need only be in order.  These are all checks of the
# caller's arguments, and the messages are built only on the branches that throw them.

# The part of the bracketing rule that a refused limit broke, for the messages below.
bracket_rule(::Type{<:Integer}, name) =
    "the range of $name must include 0, where the recurrence starts"
bracket_rule(::Type{HalfOddInteger}, name) =
    "the range of $name must include both -1//2 and 1//2, where the recurrence starts"

function validate_index_ranges(ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT) where
    {IT<:Union{Signed, HalfOddInteger}}
    # ℓₘₐₓ must be at least as big as ℓₘᵢₙ(ℓₘₐₓ)
    if ℓₘₐₓ < ℓₘᵢₙ(ℓₘₐₓ)
        throw(ArgumentError("ℓₘₐₓ=$ℓₘₐₓ must be non-negative."))
    end

    # The m′ and m ranges must be ordered correctly
    if m′ₘₐₓ < m′ₘᵢₙ
        throw(ArgumentError("m′ₘₐₓ=$m′ₘₐₓ is less than m′ₘᵢₙ=$m′ₘᵢₙ."))
    end
    if mₘₐₓ < mₘᵢₙ
        throw(ArgumentError("mₘₐₓ=$mₘₐₓ is less than mₘᵢₙ=$mₘᵢₙ."))
    end

    # The m′ and m ranges must bracket ±ℓₘᵢₙ (see above)
    if m′ₘₐₓ < ℓₘᵢₙ(ℓₘₐₓ)
        throw(ArgumentError(
            "m′ₘₐₓ=$m′ₘₐₓ is too small for this index type, $IT: $(bracket_rule(IT, "m′"))."
        ))
    end
    if m′ₘᵢₙ > -ℓₘᵢₙ(ℓₘₐₓ)
        throw(ArgumentError(
            "m′ₘᵢₙ=$m′ₘᵢₙ is too large for this index type, $IT: $(bracket_rule(IT, "m′"))."
        ))
    end
    if mₘₐₓ < ℓₘᵢₙ(ℓₘₐₓ)
        throw(ArgumentError(
            "mₘₐₓ=$mₘₐₓ is too small for this index type, $IT: $(bracket_rule(IT, "m"))."
        ))
    end
    if mₘᵢₙ > -ℓₘᵢₙ(ℓₘₐₓ)
        throw(ArgumentError(
            "mₘᵢₙ=$mₘᵢₙ is too large for this index type, $IT: $(bracket_rule(IT, "m"))."
        ))
    end

    # The m′ and m values must be in range for ℓₘₐₓ
    if abs(m′ₘₐₓ) > ℓₘₐₓ
        throw(ArgumentError("|m′ₘₐₓ|=|$m′ₘₐₓ| is too large for ℓₘₐₓ=$ℓₘₐₓ."))
    end
    if abs(m′ₘᵢₙ) > ℓₘₐₓ
        throw(ArgumentError("|m′ₘᵢₙ|=|$m′ₘᵢₙ| is too large for ℓₘₐₓ=$ℓₘₐₓ."))
    end
    if abs(mₘₐₓ) > ℓₘₐₓ
        throw(ArgumentError("|mₘₐₓ|=|$mₘₐₓ| is too large for ℓₘₐₓ=$ℓₘₐₓ."))
    end
    if abs(mₘᵢₙ) > ℓₘₐₓ
        throw(ArgumentError("|mₘᵢₙ|=|$mₘᵢₙ| is too large for ℓₘₐₓ=$ℓₘₐₓ."))
    end
    nothing
end

function validate_index_ranges(ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT) where
    {IT<:Union{Signed, HalfOddInteger}}
    # ℓₘₐₓ must be at least as big as ℓₘᵢₙ(ℓₘₐₓ)
    if ℓₘₐₓ < ℓₘᵢₙ(ℓₘₐₓ)
        throw(ArgumentError("ℓₘₐₓ=$ℓₘₐₓ must be non-negative."))
    end

    # The m′ range must be ordered correctly
    if m′ₘₐₓ < m′ₘᵢₙ
        throw(ArgumentError("m′ₘₐₓ=$m′ₘₐₓ is less than m′ₘᵢₙ=$m′ₘᵢₙ."))
    end

    # The m′ range must bracket ±ℓₘᵢₙ
    if m′ₘₐₓ < ℓₘᵢₙ(ℓₘₐₓ)
        throw(ArgumentError(
            "m′ₘₐₓ=$m′ₘₐₓ is too small for this index type, $IT: $(bracket_rule(IT, "m′"))."
        ))
    end
    if m′ₘᵢₙ > -ℓₘᵢₙ(ℓₘₐₓ)
        throw(ArgumentError(
            "m′ₘᵢₙ=$m′ₘᵢₙ is too large for this index type, $IT: $(bracket_rule(IT, "m′"))."
        ))
    end

    # The m′ values must be in range for ℓₘₐₓ
    if abs(m′ₘₐₓ) > ℓₘₐₓ
        throw(ArgumentError("|m′ₘₐₓ|=|$m′ₘₐₓ| is too large for ℓₘₐₓ=$ℓₘₐₓ."))
    end
    if abs(m′ₘᵢₙ) > ℓₘₐₓ
        throw(ArgumentError("|m′ₘᵢₙ|=|$m′ₘᵢₙ| is too large for ℓₘₐₓ=$ℓₘₐₓ."))
    end
    nothing
end

# The one-axis rule for the m axis of a `DegreeBlock`, a `SpinMatrix` or one of their
# batches.
function validate_m_range(ℓ::IT, mₘₐₓ::IT, mₘᵢₙ::IT) where {IT<:Union{Signed, HalfOddInteger}}
    if ℓ < ℓₘᵢₙ(ℓ)
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
function validate_s_range(sₘₐₓ::IT, sₘᵢₙ::IT) where {IT<:Union{Signed, HalfOddInteger}}
    if sₘₐₓ < sₘᵢₙ
        throw(ArgumentError("sₘₐₓ=$sₘₐₓ is less than sₘᵢₙ=$sₘᵢₙ."))
    end
    nothing
end


### Storage extents
#
# The natural-index accessors compare an index with a container's limits and then index the
# storage under `@inbounds`, and iteration reads every element the same way, so the storage
# must reach every element the limits describe.  That is checked in each inner constructor,
# because the inner constructors are also what `copy`, `similar` and the views `w[iᵣ]`,
# `b[s, :]`, `b[:, s, :]` and `b[iᵣ]` are built with.  Storage larger than the block is
# legitimate: a calculator's blocks sit in storage sized for its largest ℓ.  The ordering of
# the limits and their relation to ℓ are left to the outer constructors, with one exception:
# an axis of negative extent is refused here, because two of them would multiply to a
# positive `length`, which iteration would then read.  A `DegreeBlock` is checked again at
# each access, because its storage may be a caller's `Vector`, which can be resized after
# construction; see `check_storage` below.
#
# These run on every view that `w[iᵣ]` or `b[s, :]` builds, so the messages are formatted
# only on the way to an error.

@inline function check_extent(parent, d::Int, hi, lo, name::String)
    (0 ≤ Int(hi - lo) + 1 ≤ size(parent, d)) || extent_error(parent, d, hi, lo, name)
    nothing
end
@inline function check_extent(parent, Nᵣ::Int)
    (0 ≤ Nᵣ ≤ size(parent, 1)) || extent_error(parent, Nᵣ)
    nothing
end

@noinline function extent_error(parent, d::Int, hi, lo, name::String)
    n = Int(hi - lo) + 1
    # A negative lower limit is parenthesized, so that the difference reads `2-(-2)+1`.
    lo_text = lo < 0 ? "($lo)" : "$lo"
    if n < 0
        throw(ArgumentError(
            "$(name)ₘₐₓ=$hi is less than $(name)ₘᵢₙ=$lo by more than one, which would give "
            * "the block an axis of extent $n."
        ))
    elseif parent isa AbstractVector
        throw(DimensionMismatch(
            "The input data must have length at least "
            * "$(name)ₘₐₓ-$(name)ₘᵢₙ+1=$hi-$lo_text+1=$n; it is $(length(parent))."
        ))
    else
        throw(DimensionMismatch(
            "The extent of the $(("first", "second", "third")[d]) dimension in the input data "
            * "must be at least $(name)ₘₐₓ-$(name)ₘᵢₙ+1=$hi-$lo_text+1=$n; it is "
            * "$(size(parent, d))."
        ))
    end
end
@noinline function extent_error(parent, Nᵣ::Int)
    if Nᵣ < 0
        throw(ArgumentError("The number of rotors Nᵣ=$Nᵣ must not be negative."))
    else
        throw(DimensionMismatch(
            "The extent of the first dimension in the input data must be at least the number "
            * "of rotors Nᵣ=$Nᵣ; it is $(size(parent, 1))."
        ))
    end
end


@doc raw"""
    WignerMatrix{IT, NT, ST} <: AbstractWignerMatrix{IT, NT, ST}
    WignerMatrix(parent, ℓ; m′ₘₐₓ=ℓ, m′ₘᵢₙ=-m′ₘₐₓ, mₘₐₓ=ℓ, mₘᵢₙ=-mₘₐₓ)

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
where `ℓₘᵢₙ` is either 0 or 1//2 depending on whether `IT` is an integer or half-integer
type.  Both rows `±ℓₘᵢₙ` must be included because the recurrence seeds the half-integer
ladder from the pair of rows `m′ = ±1/2`; for integers this reduces to the familiar `m′ₘᵢₙ ≤
0 ≤ m′ₘₐₓ`.

The constructor wraps `parent`, which must be 1-based and at least `(m′ₘₐₓ-m′ₘᵢₙ+1) ×
(mₘₐₓ-mₘᵢₙ+1)`, without copying it.  `ℓ` and the limits must all be integers of type `Int`,
or all half-odd-integers, each a [`HalfOddInteger`](@ref) or a `Rational{Int}` with
denominator 2; the limits default to the whole block, and each lower limit defaults to minus
the corresponding upper one.  The keywords may also be spelled `mp_max`, `mp_min`, `m_max`
and `m_min`; where both spellings of one are given, the Unicode one is used.
"""
struct WignerMatrix{IT, NT, ST} <: AbstractWignerMatrix{IT, NT, ST}
    parent::ST
    ℓ::IT
    m′ₘₐₓ::IT
    m′ₘᵢₙ::IT
    mₘₐₓ::IT
    mₘᵢₙ::IT
    # The parent is indexed as 1-based throughout, much of it under `@inbounds`, so an
    # offset array would be read and written outside its storage; `ModeWeights` and
    # `HarmonicValues` refuse one in the same way.  For the same reason the parent must
    # reach every element of the block (see "Storage extents" above).
    function WignerMatrix{IT, NT, ST}(parent, ℓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ) where {IT, NT, ST}
        Base.require_one_based_indexing(parent)
        check_extent(parent, 1, m′ₘₐₓ, m′ₘᵢₙ, "m′")
        check_extent(parent, 2, mₘₐₓ, mₘᵢₙ, "m")
        new{IT, NT, ST}(parent, ℓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    end
end

"""
    Matrix(w::WignerMatrix)

Materialize the block of the Wigner matrix represented by `w` as an ordinary `Matrix`.  The
rows and columns are in order of increasing `m′` and `m`, so that the element `w[m′, m]` is
at `[Int(m′-m′ₘᵢₙ)+1, Int(m-mₘᵢₙ)+1]`.
"""
Base.Matrix(w::WignerMatrix) = Matrix(stored_elements(w))

# Indexing by `(m′, m)` is what a `WignerMatrix` means by two indices, and only that.  The
# other containers index differently — `[iᵣ, m′, m]` for `WignerMatrixBatch` and `HWedge`,
# `[s, m]` for `SpinMatrix`, `[iᵣ, m]` for `DegreeBlockBatch` and `HAxis` — and define their
# own methods; were these methods on `AbstractWignerMatrix`, a two-index call on a
# three-index container would silently read the wrong element instead of being an error.
Base.checkbounds(::Type{Bool}, w::WignerMatrix{IT}, m′::IT, m::IT) where {IT} =
    inrange(m′, w.m′ₘᵢₙ, w.m′ₘₐₓ) && inrange(m, w.mₘᵢₙ, w.mₘₐₓ)

@propagate_inbounds function Base.getindex(
    w::WignerMatrix{IT}, m′::IT, m::IT
) where {IT<:IntegerHalf}
    @boundscheck checkbounds(Bool, w, m′, m) || throw(BoundsError(w, (m′, m)))
    @inbounds Base.parent(w)[(m′-m′ₘᵢₙ(w))+1, (m-mₘᵢₙ(w))+1]
end

@propagate_inbounds function Base.setindex!(
    w::WignerMatrix{IT}, v, m′::IT, m::IT
) where {IT<:IntegerHalf}
    @boundscheck checkbounds(Bool, w, m′, m) || throw(BoundsError(w, (m′, m)))
    @inbounds Base.parent(w)[(m′-m′ₘᵢₙ(w))+1, (m-mₘᵢₙ(w))+1] = v
end


### Indexing with indices of other types.
#
# The accessors above take indices of the container's own type, `Int` or `HalfOddInteger`.
# A half-integer container is indexed by `HalfOddInteger`s, but `w[1//2, -3//2]` is what a
# caller naturally writes, so each block container also has methods that take any
# `IndexType`, convert each natural index with `container_index` below, and re-dispatch;
# either natural index may be written either way.  An index that is not of the container's
# kind, or that the index methods would refuse — a narrow integer, or a `Rational` whose
# denominator is not 2 or whose integer type is not `Int` — is refused there with the
# reason, rather than with the bare `MethodError` that dispatch alone would give.  These
# methods are less specific than the accessors of the container's own index type, so an
# index of that type, which is what every loop over the elements passes, never reaches them.
# That needs the bound `IT<:IntegerHalf` on the accessors: without it, their signature would
# admit index types outside `IndexType`, and would not be the more specific of the two.

# The index of a container's own kind, from an index as a caller wrote it.  The accessors of
# the containers are written by hand, because the kind they accept is that of the container
# rather than one that dispatch on the arguments could choose, but they accept what the
# index methods accept: an `Int` for an integer container, and a `HalfOddInteger` or a
# `Rational{Int}` with denominator 2 for a half-integer one.  Anything else is refused with
# the reason the index methods give; a narrow integer, in particular, is told to be
# converted, because the arithmetic that finds an element's position is not closed under it.
@inline container_index(::Type{Int}, x::Int, c, name) = x
@inline container_index(::Type{HalfOddInteger}, x::HalfOddInteger, c, name) = x
@inline function container_index(::Type{HalfOddInteger}, x::Rational{Int}, c, name)
    is_half_odd_index(x) || throw(container_index_error(HalfOddInteger, x, c, name))
    half_odd_index(x)
end
container_index(::Type{IT}, x, c, name) where {IT} =
    throw(container_index_error(IT, x, c, name))
@noinline function container_index_error(::Type{IT}, x, c, name) where {IT}
    kind = IT === HalfOddInteger ? (
        "half-odd-integers, each a `HalfOddInteger` or a `Rational{Int}` with denominator 2, "
        * "like 7//2"
    ) : "integers of type `Int`, like 3"
    message = (
        "The indices of this `$(container_name(c))` are $kind; got $name = $(typed_repr(x))."
    )
    index_kind(x) === nothing && (message *= "  " * index_problem(x))
    ArgumentError(message)
end
# The name of a container, or of a calculator whose index is checked in the same way, as a
# message gives it.
container_name(c) = nameof(typeof(c))

@propagate_inbounds Base.getindex(
    w::WignerMatrix{IT}, m′::IndexType, m::IndexType
) where {IT} = w[container_index(IT, m′, w, "m′"), container_index(IT, m, w, "m")]
@propagate_inbounds Base.setindex!(
    w::WignerMatrix{IT}, v, m′::IndexType, m::IndexType
) where {IT} = (w[container_index(IT, m′, w, "m′"), container_index(IT, m, w, "m")] = v)

# Column-major over the *block* represented, whose elements are in range by construction
# (the parent storage may be larger).
function Base.iterate(w::WignerMatrix, state=1)
    n₁, n₂ = size(w)
    state > n₁ * n₂ && return nothing
    i, j = (state - 1) % n₁, (state - 1) ÷ n₁
    (@inbounds(w[w.m′ₘᵢₙ + i, w.mₘᵢₙ + j]), state + 1)
end

@index_methods function WignerMatrix(
    parent::ST, ℓ::IT;
    mp_max::IndexType=ℓ, m′ₘₐₓ::IndexType=mp_max,
    mp_min::IndexType=-m′ₘₐₓ, m′ₘᵢₙ::IndexType=mp_min,
    m_max::IndexType=ℓ, mₘₐₓ::IndexType=m_max,
    m_min::IndexType=-mₘₐₓ, mₘᵢₙ::IndexType=m_min
) where {IT<:IndexType, NT, ST<:AbstractMatrix{NT}}
    validate_index_ranges(ℓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
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
function Base.similar(w::WignerMatrix{IT}, ::Type{T}=eltype(w)) where {IT, T}
    let p = Matrix{T}(undef, size(w))
        WignerMatrix{IT, T, typeof(p)}(p, w.ℓ, w.m′ₘₐₓ, w.m′ₘᵢₙ, w.mₘₐₓ, w.mₘᵢₙ)
    end
end


"""
    WignerMatrixBatch{IT, NT, ST} <: AbstractWignerMatrix{IT, NT, ST}
    WignerMatrixBatch(parent, ℓ; m′ₘₐₓ=ℓ, m′ₘᵢₙ=-m′ₘₐₓ, mₘₐₓ=ℓ, mₘᵢₙ=-mₘₐₓ)

`Nᵣ` Wigner matrices of one ``ℓ``, stored together and indexed as `w[iᵣ, m′, m]`.  This is
what [`recurrence!`](@ref) returns for a calculator built from a vector of rotor data, of
any length, for either kind of index.

The storage `parent(w)` is 1-based and 3-dimensional, ordered `[iᵣ, m′, m]`, exactly as in
the calculator, and `Nᵣ` is its first extent.  The constructor takes the indices and the
limits, with their ASCII spellings, exactly as [`WignerMatrix`](@ref) does, and under the
same rules.  `w[iᵣ]` gives the [`WignerMatrix`](@ref) view of one rotor's matrix, which is
then indexed naturally as `w[iᵣ][m′, m]`.

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
    # As for `WignerMatrix`: the parent must be 1-based, and must reach every element.
    function WignerMatrixBatch{IT, NT, ST}(parent, ℓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, Nᵣ) where {IT, NT, ST}
        Base.require_one_based_indexing(parent)
        check_extent(parent, Nᵣ)
        check_extent(parent, 2, m′ₘₐₓ, m′ₘᵢₙ, "m′")
        check_extent(parent, 3, mₘₐₓ, mₘᵢₙ, "m")
        new{IT, NT, ST}(parent, ℓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, Nᵣ)
    end
end

@index_methods function WignerMatrixBatch(
    parent::ST, ℓ::IT;
    mp_max::IndexType=ℓ, m′ₘₐₓ::IndexType=mp_max,
    mp_min::IndexType=-m′ₘₐₓ, m′ₘᵢₙ::IndexType=mp_min,
    m_max::IndexType=ℓ, mₘₐₓ::IndexType=m_max,
    m_min::IndexType=-mₘₐₓ, mₘᵢₙ::IndexType=m_min
) where {IT<:IndexType, NT, ST<:AbstractArray{NT, 3}}
    validate_index_ranges(ℓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    WignerMatrixBatch{IT, NT, ST}(parent, ℓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, size(parent, 1))
end

Nᵣ(w::WignerMatrixBatch) = w.Nᵣ

Base.checkbounds(::Type{Bool}, w::WignerMatrixBatch{IT}, iᵣ::Integer, m′::IT, m::IT) where {IT} =
    1 ≤ iᵣ ≤ w.Nᵣ && inrange(m′, w.m′ₘᵢₙ, w.m′ₘₐₓ) && inrange(m, w.mₘᵢₙ, w.mₘₐₓ)

@propagate_inbounds function Base.getindex(
    w::WignerMatrixBatch{IT}, iᵣ::Integer, m′::IT, m::IT
) where {IT<:IntegerHalf}
    @boundscheck checkbounds(Bool, w, iᵣ, m′, m) || throw(BoundsError(w, (iᵣ, m′, m)))
    @inbounds parent(w)[iᵣ, Int(m′ - w.m′ₘᵢₙ) + 1, Int(m - w.mₘᵢₙ) + 1]
end

@propagate_inbounds function Base.setindex!(
    w::WignerMatrixBatch{IT}, v, iᵣ::Integer, m′::IT, m::IT
) where {IT<:IntegerHalf}
    @boundscheck checkbounds(Bool, w, iᵣ, m′, m) || throw(BoundsError(w, (iᵣ, m′, m)))
    @inbounds parent(w)[iᵣ, Int(m′ - w.m′ₘᵢₙ) + 1, Int(m - w.mₘᵢₙ) + 1] = v
end

# See the note on indexing with indices of other types above.
@propagate_inbounds Base.getindex(
    w::WignerMatrixBatch{IT}, iᵣ::Integer, m′::IndexType, m::IndexType
) where {IT} = w[iᵣ, container_index(IT, m′, w, "m′"), container_index(IT, m, w, "m")]
@propagate_inbounds Base.setindex!(
    w::WignerMatrixBatch{IT}, v, iᵣ::Integer, m′::IndexType, m::IndexType
) where {IT} = (w[iᵣ, container_index(IT, m′, w, "m′"), container_index(IT, m, w, "m")] = v)

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
function Base.similar(w::WignerMatrixBatch{IT}, ::Type{T}=eltype(w)) where {IT, T}
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

Base.Matrix(w::WignerMatrixBatch) = throw(ArgumentError(
    "A WignerMatrixBatch is 3-dimensional; use `Array(w)`, or `Matrix(w[iᵣ])`."
))


"""
    DegreeBlock{IT, NT, ST}

One harmonic degree's worth of values, indexed naturally by the order ``m``: `v[m]` for
``mₘᵢₙ ≤ m ≤ mₘₐₓ``.  This is the 1-dimensional sibling of [`WignerMatrix`](@ref).

The name is deliberately neutral, because the same shape serves three things: a block of
spin-weighted harmonics at one ``ℓ``, from [`sYlm`](@ref) or [`sYlmCalculator`](@ref)'s
`ₛYₗ[s, :]`; one ``ℓ``'s worth of mode weights, from [`ModeWeights`](@ref)'s `w[ℓ, :]`; and
one spin weight's row of a [`SpinMatrix`](@ref), from `b[s, :]`.  All three are values at a
fixed degree indexed by order, whatever they mean.

The storage `parent(v)` is 1-based, in order of increasing ``m``.  The constructor

    DegreeBlock(parent, ℓ; mₘₐₓ=ℓ, mₘᵢₙ=-mₘₐₓ)

wraps it without copying, and requires `-ℓ ≤ mₘᵢₙ ≤ mₘₐₓ ≤ ℓ`; unlike the limits of a
[`WignerMatrix`](@ref), these need not bracket zero, because no recurrence fills the block.
`ℓ` and the limits must all be integers of type `Int`, or all half-odd-integers, each a
[`HalfOddInteger`](@ref) or a `Rational{Int}` with denominator 2.  The keywords may also be
spelled `m_max` and `m_min`.

See also [`DegreeBlockBatch`](@ref).
"""
struct DegreeBlock{IT, NT, ST<:AbstractVector{NT}} <: AbstractWignerMatrix{IT, NT, ST}
    parent::ST
    ℓ::IT
    mₘₐₓ::IT
    mₘᵢₙ::IT
    # As for `WignerMatrix`: the parent must be 1-based, and must reach every element.
    function DegreeBlock{IT, NT, ST}(parent, ℓ, mₘₐₓ, mₘᵢₙ) where {IT, NT, ST}
        Base.require_one_based_indexing(parent)
        check_extent(parent, 1, mₘₐₓ, mₘᵢₙ, "m")
        new{IT, NT, ST}(parent, ℓ, mₘₐₓ, mₘᵢₙ)
    end
end

@index_methods function DegreeBlock(
    parent::ST, ℓ::IT;
    m_max::IndexType=ℓ, mₘₐₓ::IndexType=m_max,
    m_min::IndexType=-mₘₐₓ, mₘᵢₙ::IndexType=m_min
) where {IT<:IndexType, NT, ST<:AbstractVector{NT}}
    validate_m_range(ℓ, mₘₐₓ, mₘᵢₙ)
    DegreeBlock{IT, NT, ST}(parent, ℓ, mₘₐₓ, mₘᵢₙ)
end

Base.firstindex(v::DegreeBlock) = v.mₘᵢₙ
Base.lastindex(v::DegreeBlock) = v.mₘₐₓ
Base.keys(v::DegreeBlock) = v.mₘᵢₙ:v.mₘₐₓ

# Of all the blocks, only a `DegreeBlock` may have a caller's `Vector` as its storage — from
# `DegreeBlock(v, ℓ)` or `relabel`, for example — which can be resized after the constructor
# has compared its length with the limits.  So the accessors compare the position of the
# element with the length of the storage as well as with the limits, and iteration compares
# each position before it reads.
@inline function check_storage(v::DegreeBlock, i)
    if i > length(parent(v))
        throw(storage_error(v, i))
    end
    nothing
end
@noinline function storage_error(v::DegreeBlock, i)
    DimensionMismatch(
        "The storage of this DegreeBlock for ℓ=$(v.ℓ), with m ∈ $(v.mₘᵢₙ):$(v.mₘₐₓ), has "
        * "length $(length(parent(v))), but the limits need an entry at position $i.  A "
        * "`DegreeBlock` uses its vector as storage without copying it, so the vector must not "
        * "be resized."
    )
end
# The elements are the leading `length(v)` entries of the storage, which must still be
# there.
function stored_elements(v::DegreeBlock)
    check_storage(v, length(v))
    view(parent(v), 1:length(v))
end

Base.checkbounds(::Type{Bool}, v::DegreeBlock{IT}, m::IT) where {IT} =
    inrange(m, v.mₘᵢₙ, v.mₘₐₓ)

@propagate_inbounds function Base.getindex(v::DegreeBlock{IT}, m::IT) where {IT<:IntegerHalf}
    i = Int(m - v.mₘᵢₙ) + 1
    @boundscheck checkbounds(Bool, v, m) || throw(BoundsError(v, m))
    @boundscheck check_storage(v, i)
    @inbounds parent(v)[i]
end
@propagate_inbounds function Base.setindex!(
    v::DegreeBlock{IT}, x, m::IT
) where {IT<:IntegerHalf}
    i = Int(m - v.mₘᵢₙ) + 1
    @boundscheck checkbounds(Bool, v, m) || throw(BoundsError(v, m))
    @boundscheck check_storage(v, i)
    @inbounds parent(v)[i] = x
end

# See the note on indexing with indices of other types above.
@propagate_inbounds Base.getindex(v::DegreeBlock{IT}, m::IndexType) where {IT} =
    v[container_index(IT, m, v, "m")]
@propagate_inbounds Base.setindex!(v::DegreeBlock{IT}, x, m::IndexType) where {IT} =
    (v[container_index(IT, m, v, "m")] = x)

# Element `state` of the block is at position `state` of its storage.
function Base.iterate(v::DegreeBlock, state=1)
    state > length(v) && return nothing
    check_storage(v, state)
    (@inbounds parent(v)[state], state + 1)
end
Base.Vector(v::DegreeBlock) = Vector(stored_elements(v))
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
function Base.similar(v::DegreeBlock{IT}, ::Type{T}=eltype(v)) where {IT, T}
    let p = Vector{T}(undef, length(v))
        DegreeBlock{IT, T, typeof(p)}(p, v.ℓ, v.mₘₐₓ, v.mₘᵢₙ)
    end
end
function Base.summary(io::IO, v::DegreeBlock{IT, NT}) where {IT, NT}
    print(io, "(", v.mₘᵢₙ, ":", v.mₘₐₓ, ") DegreeBlock{", IT, ", ", NT, "} for ℓ=", v.ℓ)
end
Base.show(io::IO, v::DegreeBlock) = summary(io, v)
function Base.show(io::IO, ::MIME"text/plain", v::DegreeBlock)
    summary(io, v)
    println(io, ":")
    Base.print_array(io, stored_elements(v))
end


"""
    DegreeBlockBatch{IT, NT, ST}

`Nᵣ` rows of values of a single ``ℓ``, stored together and indexed as `v[iᵣ, m]`.  This is
the 1-dimensional sibling of [`WignerMatrixBatch`](@ref), and is what a batched
[`sYlmCalculator`](@ref) built for one spin weight yields, for either kind of index.
`v[iᵣ]` gives the [`DegreeBlock`](@ref) view of one rotor's row.

The storage `parent(v)` is 1-based and 2-dimensional, ordered `[iᵣ, m]`, and `Nᵣ` is its
first extent.  The constructor, `DegreeBlockBatch(parent, ℓ; mₘₐₓ=ℓ, mₘᵢₙ=-mₘₐₓ)`, takes the
indices and the limits exactly as [`DegreeBlock`](@ref) does, and under the same rules.
"""
struct DegreeBlockBatch{IT, NT, ST<:AbstractMatrix{NT}} <: AbstractWignerMatrix{IT, NT, ST}
    parent::ST
    ℓ::IT
    mₘₐₓ::IT
    mₘᵢₙ::IT
    Nᵣ::Int
    # As for `WignerMatrix`: the parent must be 1-based, and must reach every element.
    function DegreeBlockBatch{IT, NT, ST}(parent, ℓ, mₘₐₓ, mₘᵢₙ, Nᵣ) where {IT, NT, ST}
        Base.require_one_based_indexing(parent)
        check_extent(parent, Nᵣ)
        check_extent(parent, 2, mₘₐₓ, mₘᵢₙ, "m")
        new{IT, NT, ST}(parent, ℓ, mₘₐₓ, mₘᵢₙ, Nᵣ)
    end
end

@index_methods function DegreeBlockBatch(
    parent::ST, ℓ::IT;
    m_max::IndexType=ℓ, mₘₐₓ::IndexType=m_max,
    m_min::IndexType=-mₘₐₓ, mₘᵢₙ::IndexType=m_min
) where {IT<:IndexType, NT, ST<:AbstractMatrix{NT}}
    validate_m_range(ℓ, mₘₐₓ, mₘᵢₙ)
    DegreeBlockBatch{IT, NT, ST}(parent, ℓ, mₘₐₓ, mₘᵢₙ, size(parent, 1))
end

Nᵣ(v::DegreeBlockBatch) = v.Nᵣ

Base.checkbounds(::Type{Bool}, v::DegreeBlockBatch{IT}, iᵣ::Integer, m::IT) where {IT} =
    1 ≤ iᵣ ≤ v.Nᵣ && inrange(m, v.mₘᵢₙ, v.mₘₐₓ)

@propagate_inbounds function Base.getindex(
    v::DegreeBlockBatch{IT}, iᵣ::Integer, m::IT
) where {IT<:IntegerHalf}
    @boundscheck checkbounds(Bool, v, iᵣ, m) || throw(BoundsError(v, (iᵣ, m)))
    @inbounds parent(v)[iᵣ, Int(m - v.mₘᵢₙ) + 1]
end
@propagate_inbounds function Base.setindex!(
    v::DegreeBlockBatch{IT}, x, iᵣ::Integer, m::IT
) where {IT<:IntegerHalf}
    @boundscheck checkbounds(Bool, v, iᵣ, m) || throw(BoundsError(v, (iᵣ, m)))
    @inbounds parent(v)[iᵣ, Int(m - v.mₘᵢₙ) + 1] = x
end

# See the note on indexing with indices of other types above.
@propagate_inbounds Base.getindex(v::DegreeBlockBatch{IT}, iᵣ::Integer, m::IndexType) where {IT} =
    v[iᵣ, container_index(IT, m, v, "m")]
@propagate_inbounds Base.setindex!(
    v::DegreeBlockBatch{IT}, x, iᵣ::Integer, m::IndexType
) where {IT} = (v[iᵣ, container_index(IT, m, v, "m")] = x)

@propagate_inbounds function Base.getindex(v::DegreeBlockBatch{IT, NT}, iᵣ::Integer) where {IT, NT}
    @boundscheck if !(1 ≤ iᵣ ≤ v.Nᵣ)
        throw(BoundsError(v, (iᵣ,)))
    end
    let p = view(parent(v), iᵣ, :)
        DegreeBlock{IT, NT, typeof(p)}(p, v.ℓ, v.mₘₐₓ, v.mₘᵢₙ)
    end
end

Base.Matrix(v::DegreeBlockBatch) = Matrix(stored_elements(v))
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
function Base.similar(v::DegreeBlockBatch{IT}, ::Type{T}=eltype(v)) where {IT, T}
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
    Base.print_array(io, stored_elements(v))
end


"""
    SpinMatrix{IT, NT, ST}

The values of a single ``ℓ`` for a range of spin weights, indexed naturally by ``(s, m)``:
`b[s, m]` for ``sₘᵢₙ ≤ s ≤ sₘₐₓ`` and ``mₘᵢₙ ≤ m ≤ mₘₐₓ``, and `b[s, :]` for one whole row
as a [`DegreeBlock`](@ref).  This is what an [`sYlmCalculator`](@ref) built for a range of
spin weights yields for each ``ℓ``.

The storage `parent(b)` is 1-based and 2-dimensional, ordered `[s, m]`.  The constructor

    SpinMatrix(parent, ℓ; sₘₐₓ, sₘᵢₙ, mₘₐₓ=ℓ, mₘᵢₙ=-mₘₐₓ)

wraps it without copying.  The spin limits are required, and the keywords may also be
spelled `s_max`, `s_min`, `m_max` and `m_min`.  `ℓ` and every limit must all be integers of
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
struct SpinMatrix{IT, NT, ST<:AbstractMatrix{NT}} <: AbstractWignerMatrix{IT, NT, ST}
    parent::ST
    ℓ::IT
    sₘₐₓ::IT
    sₘᵢₙ::IT
    mₘₐₓ::IT
    mₘᵢₙ::IT
    # As for `WignerMatrix`: the parent must be 1-based, and must reach every element.
    function SpinMatrix{IT, NT, ST}(parent, ℓ, sₘₐₓ, sₘᵢₙ, mₘₐₓ, mₘᵢₙ) where {IT, NT, ST}
        Base.require_one_based_indexing(parent)
        check_extent(parent, 1, sₘₐₓ, sₘᵢₙ, "s")
        check_extent(parent, 2, mₘₐₓ, mₘᵢₙ, "m")
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
) where {IT<:IndexType, NT, ST<:AbstractMatrix{NT}}
    sₘₐₓ === nothing && throw(UndefKeywordError(:sₘₐₓ))
    sₘᵢₙ === nothing && throw(UndefKeywordError(:sₘᵢₙ))
    validate_s_range(sₘₐₓ, sₘᵢₙ)
    validate_m_range(ℓ, mₘₐₓ, mₘᵢₙ)
    SpinMatrix{IT, NT, ST}(parent, ℓ, sₘₐₓ, sₘᵢₙ, mₘₐₓ, mₘᵢₙ)
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
    @inbounds parent(b)[Int(s - b.sₘᵢₙ) + 1, Int(m - b.mₘᵢₙ) + 1]
end
@propagate_inbounds function Base.setindex!(
    b::SpinMatrix{IT}, x, s::IT, m::IT
) where {IT<:IntegerHalf}
    @boundscheck checkbounds(Bool, b, s, m) || throw(BoundsError(b, (s, m)))
    @inbounds parent(b)[Int(s - b.sₘᵢₙ) + 1, Int(m - b.mₘᵢₙ) + 1] = x
end

# See the note on indexing with indices of other types above.
@propagate_inbounds Base.getindex(
    b::SpinMatrix{IT}, s::IndexType, m::IndexType
) where {IT} = b[container_index(IT, s, b, "s"), container_index(IT, m, b, "m")]
@propagate_inbounds Base.setindex!(
    b::SpinMatrix{IT}, x, s::IndexType, m::IndexType
) where {IT} = (b[container_index(IT, s, b, "s"), container_index(IT, m, b, "m")] = x)

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
    let p = view(parent(b), Int(s - b.sₘᵢₙ) + 1, :)
        DegreeBlock{IT, NT, typeof(p)}(p, b.ℓ, b.mₘₐₓ, b.mₘᵢₙ)
    end
end
@propagate_inbounds Base.getindex(b::SpinMatrix{IT}, s::IndexType, ::Colon) where {IT} =
    b[container_index(IT, s, b, "s"), :]

Base.Matrix(b::SpinMatrix) = Matrix(stored_elements(b))
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
function Base.similar(b::SpinMatrix{IT}, ::Type{T}=eltype(b)) where {IT, T}
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
    Base.print_array(io, stored_elements(b))
end


"""
    SpinMatrixBatch{IT, NT, ST}

`Nᵣ` [`SpinMatrix`](@ref) blocks of one ``ℓ``, stored together and indexed as `b[iᵣ, s, m]`.
This is what a batched [`sYlmCalculator`](@ref) yields when it was built for a range of spin
weights.  `b[iᵣ]` gives the `SpinMatrix` view of one rotor's block, which is then indexed
naturally as `b[iᵣ][s, m]`.

The storage `parent(b)` is 1-based and 3-dimensional, ordered `[iᵣ, s, m]`, exactly as in
the calculator, and `Nᵣ` is its first extent.  The constructor, `SpinMatrixBatch(parent, ℓ;
sₘₐₓ, sₘᵢₙ, mₘₐₓ=ℓ, mₘᵢₙ=-mₘₐₓ)`, takes the indices and the limits exactly as
[`SpinMatrix`](@ref) does, and under the same rules.

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
    # As for `WignerMatrix`: the parent must be 1-based, and must reach every element.
    function SpinMatrixBatch{IT, NT, ST}(parent, ℓ, sₘₐₓ, sₘᵢₙ, mₘₐₓ, mₘᵢₙ, Nᵣ) where {IT, NT, ST}
        Base.require_one_based_indexing(parent)
        check_extent(parent, Nᵣ)
        check_extent(parent, 2, sₘₐₓ, sₘᵢₙ, "s")
        check_extent(parent, 3, mₘₐₓ, mₘᵢₙ, "m")
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
) where {IT<:IndexType, NT, ST<:AbstractArray{NT, 3}}
    sₘₐₓ === nothing && throw(UndefKeywordError(:sₘₐₓ))
    sₘᵢₙ === nothing && throw(UndefKeywordError(:sₘᵢₙ))
    validate_s_range(sₘₐₓ, sₘᵢₙ)
    validate_m_range(ℓ, mₘₐₓ, mₘᵢₙ)
    SpinMatrixBatch{IT, NT, ST}(parent, ℓ, sₘₐₓ, sₘᵢₙ, mₘₐₓ, mₘᵢₙ, size(parent, 1))
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
    @inbounds parent(b)[iᵣ, Int(s - b.sₘᵢₙ) + 1, Int(m - b.mₘᵢₙ) + 1]
end
@propagate_inbounds function Base.setindex!(
    b::SpinMatrixBatch{IT}, x, iᵣ::Integer, s::IT, m::IT
) where {IT<:IntegerHalf}
    @boundscheck checkbounds(Bool, b, iᵣ, s, m) || throw(BoundsError(b, (iᵣ, s, m)))
    @inbounds parent(b)[iᵣ, Int(s - b.sₘᵢₙ) + 1, Int(m - b.mₘᵢₙ) + 1] = x
end

# See the note on indexing with indices of other types above.
@propagate_inbounds Base.getindex(
    b::SpinMatrixBatch{IT}, iᵣ::Integer, s::IndexType, m::IndexType
) where {IT} = b[iᵣ, container_index(IT, s, b, "s"), container_index(IT, m, b, "m")]
@propagate_inbounds Base.setindex!(
    b::SpinMatrixBatch{IT}, x, iᵣ::Integer, s::IndexType, m::IndexType
) where {IT} = (b[iᵣ, container_index(IT, s, b, "s"), container_index(IT, m, b, "m")] = x)

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
    let p = view(parent(b), :, Int(s - b.sₘᵢₙ) + 1, :)
        DegreeBlockBatch{IT, NT, typeof(p)}(p, b.ℓ, b.mₘₐₓ, b.mₘᵢₙ, b.Nᵣ)
    end
end
@propagate_inbounds Base.getindex(
    b::SpinMatrixBatch{IT}, ::Colon, s::IndexType, ::Colon
) where {IT} = b[:, container_index(IT, s, b, "s"), :]

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
function Base.similar(b::SpinMatrixBatch{IT}, ::Type{T}=eltype(b)) where {IT, T}
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
    Base.print_array(io, stored_elements(b))
end


"""
    WignerSeries{IT, VT}
    WignerSeries(blocks, ℓₘᵢₙ, ℓₘₐₓ)

The blocks of a Wigner matrix for every ``ℓ`` from `ℓₘᵢₙ` to `ℓₘₐₓ`, indexed by ``ℓ``:
`s[ℓ]` is the block of degree `ℓ`, and `s[ℓ][m′, m]` an element of it.  For a half-integer
series `ℓ` may be written as a [`HalfOddInteger`](@ref) or as a `Rational{Int}` with
denominator 2, and for an integer series it is an `Int`.

This is what [`D`](@ref) and [`d`](@ref) return, for either kind of index.  Like a
calculator, it iterates as `ℓ => block` pairs, so that `for (ℓ, 𝔇ˡ) ∈ D(R, ℓₘₐₓ)` reads
exactly as the same loop over a [`DCalculator`](@ref); `keys` is the range of ``ℓ``,
`values` gives the blocks alone, and `length` counts them.  `first` and `last` give the
first and last blocks, as indexing does, and so do `first(s, n)` and `last(s, n)`, which
give vectors of the first or last `n` blocks, and `only(s)`, the one block of a series that
has only one.

Two series are `==`, `isequal` or `≈` when they have the same range of ``ℓ`` and their
blocks are, block by block; `≈` applies its tolerances to each block separately.

The constructor takes a 1-based vector of blocks without copying it, and requires block `i`
to have ``ℓ = ℓₘᵢₙ + i - 1``, with the index type of the bounds.  The bounds must both be
integers of type `Int`, or both half-odd-integers, each a [`HalfOddInteger`](@ref) or a
`Rational{Int}` with denominator 2.

See also [`WignerMatrix`](@ref) and [`WignerMatrixBatch`](@ref).
"""
struct WignerSeries{IT<:IntegerHalf, VT<:AbstractVector}
    blocks::VT  # 1-based; blocks[i] is the block for ℓ = ℓₘᵢₙ + (i-1)
    ℓₘᵢₙ::IT
    ℓₘₐₓ::IT
    # Indexing finds the block of each ℓ at the position its label gives, so each block is
    # compared with its position once, here, rather than on every access.
    @index_methods function WignerSeries(
        blocks::VT, ℓₘᵢₙ::IT, ℓₘₐₓ::IT
    ) where {IT<:IndexType, VT<:AbstractVector}
        Base.require_one_based_indexing(blocks)  # `blocks[i]` is the block for ℓₘᵢₙ + (i-1)
        if length(blocks) != (ℓₘₐₓ - ℓₘᵢₙ) + 1
            throw(DimensionMismatch(
                "Got $(length(blocks)) blocks, but ℓ ∈ $ℓₘᵢₙ:$ℓₘₐₓ needs "
                * "$(max(0, (ℓₘₐₓ - ℓₘᵢₙ) + 1))."
            ))
        end
        for (i, block) ∈ enumerate(blocks)
            check_series_block(block, i, ℓₘᵢₙ)
        end
        new{IT, VT}(blocks, ℓₘᵢₙ, ℓₘₐₓ)
    end
end

@inline function check_series_block(block::AbstractWignerMatrix{IT}, i, ℓₘᵢₙ::IT) where {IT}
    ℓ(block) == ℓₘᵢₙ + (i - 1) || throw(series_block_error(block, i, ℓₘᵢₙ))
    nothing
end
check_series_block(block, i, ℓₘᵢₙ) = throw(series_block_error(block, i, ℓₘᵢₙ))
@noinline function series_block_error(block, i, ℓₘᵢₙ)
    if !(block isa AbstractWignerMatrix)
        ArgumentError(
            "The blocks of a `WignerSeries` are Wigner blocks such as `WignerMatrix`es; block "
            * "$i is a `$(typeof(block))`."
        )
    elseif !isa(ℓ(block), typeof(ℓₘᵢₙ))
        ArgumentError(
            "Block $i has ℓ=$(ℓ(block)), of type `$(typeof(ℓ(block)))`, but the bounds of the "
            * "series are of type `$(typeof(ℓₘᵢₙ))`, and so must the ℓ of every block be."
        )
    else
        ArgumentError(
            "Block $i is for ℓ=$(ℓ(block)), but in a series starting at ℓₘᵢₙ=$ℓₘᵢₙ block $i "
            * "must be for ℓ=$(ℓₘᵢₙ + (i - 1)), since block i is for ℓ = ℓₘᵢₙ + i - 1."
        )
    end
end

ℓₘᵢₙ(s::WignerSeries) = s.ℓₘᵢₙ
ℓₘₐₓ(s::WignerSeries) = s.ℓₘₐₓ
Base.parent(s::WignerSeries) = s.blocks
Base.length(s::WignerSeries) = length(s.blocks)
Base.axes(s::WignerSeries) = (WignerRange(s.ℓₘᵢₙ:s.ℓₘₐₓ),)
Base.axes(s::WignerSeries, d::Integer) = d ≤ 1 ? axes(s)[d] : Base.OneTo(1)
Base.ndims(::WignerSeries) = 1
Base.ndims(::Type{<:WignerSeries}) = 1
Base.size(s::WignerSeries) = (length(s),)
Base.size(s::WignerSeries, d::Integer) = d ≤ 1 ? length(s) : 1
Base.keys(s::WignerSeries) = s.ℓₘᵢₙ:s.ℓₘₐₓ
Base.firstindex(s::WignerSeries) = s.ℓₘᵢₙ
Base.lastindex(s::WignerSeries) = s.ℓₘₐₓ

# A series iterates as `ℓ => block` pairs, as a calculator and a `HarmonicValues` do, so
# that a loop over `D(R, ℓₘₐₓ)` reads exactly as one over `DCalculator(R, ℓₘₐₓ)`; `values`
# gives the blocks alone.  The iteration is over the series' own position, rather than
# handing an integer state to the storage's own `iterate`: for storage other than a
# `Vector`, such as a view, that state is not an integer, and a series on a view would stop
# after its first block while its `length` still counted them all — so that a comprehension
# over it would return uninitialized memory.  Indexing, `first` and `last` give blocks, as
# `s[ℓ]` does, and so do the forms of `first` and `last` that take a count, and `only`,
# which `Base` would otherwise derive from the iteration, as pairs.
function Base.iterate(s::WignerSeries, i::Int=1)
    i == 1 && check_blocks(s)
    i > length(s.blocks) && return nothing
    ((s.ℓₘᵢₙ + (i - 1)) => s.blocks[i], i + 1)
end
Base.eltype(::Type{<:WignerSeries{IT, VT}}) where {IT, VT} = Pair{IT, eltype(VT)}
Base.eltype(s::WignerSeries) = eltype(typeof(s))
Base.values(s::WignerSeries) = s.blocks
Base.pairs(s::WignerSeries) = s
Base.first(s::WignerSeries) = s[s.ℓₘᵢₙ]
Base.last(s::WignerSeries) = s[s.ℓₘₐₓ]
Base.first(s::WignerSeries, n::Integer) = (check_blocks(s); first(s.blocks, n))
Base.last(s::WignerSeries, n::Integer) = (check_blocks(s); last(s.blocks, n))
Base.only(s::WignerSeries) = (check_blocks(s); only(s.blocks))
ishalfinteger(::WignerSeries{IT}) where {IT<:Integer} = false
ishalfinteger(::WignerSeries{IT}) where {IT<:HalfOddInteger} = true

@propagate_inbounds function Base.getindex(s::WignerSeries{IT}, ℓ) where {IT}
    # Deliberately *not* inside `@boundscheck`: an index of the wrong kind, such as a whole
    # number asked of a half-integer series, is told what the series takes.
    let ℓ = container_index(IT, ℓ, s, "ℓ")
        @boundscheck begin
            if ℓ < s.ℓₘᵢₙ || ℓ > s.ℓₘₐₓ
                throw(BoundsError(s, ℓ))
            end
            check_blocks(s)
        end
        @inbounds s.blocks[Int(ℓ - s.ℓₘᵢₙ) + 1]
    end
end

# `values(s)` and `parent(s)` hand out the series' own vector of blocks, which can be
# resized, while indexing and iteration find the block of each ℓ at the position its label
# gives; a vector of any other length would put the wrong block, or none, at that position.
@inline function check_blocks(s::WignerSeries)
    if length(s.blocks) != Int(s.ℓₘₐₓ - s.ℓₘᵢₙ) + 1
        throw(blocks_error(s))
    end
    nothing
end
@noinline function blocks_error(s::WignerSeries)
    DimensionMismatch(
        "A series for ℓ ∈ $(s.ℓₘᵢₙ):$(s.ℓₘₐₓ) needs $(Int(s.ℓₘₐₓ - s.ℓₘᵢₙ) + 1) blocks, but "
        * "its vector of blocks has length $(length(s.blocks)); it was resized after the "
        * "series was built."
    )
end

Base.copy(s::WignerSeries) = WignerSeries(map(copy, s.blocks), s.ℓₘᵢₙ, s.ℓₘₐₓ)
Base.similar(s::WignerSeries) = WignerSeries(map(similar, s.blocks), s.ℓₘᵢₙ, s.ℓₘₐₓ)
# Block by block, after the range of ℓ.  The blocks' own comparisons count their labels, and
# `hash` agrees with `isequal`, since it is built from the range and the blocks' hashes.
same_range(s1::WignerSeries, s2::WignerSeries) =
    ℓₘᵢₙ(s1) == ℓₘᵢₙ(s2) && ℓₘₐₓ(s1) == ℓₘₐₓ(s2) && length(s1.blocks) == length(s2.blocks)
Base.:(==)(s1::WignerSeries, s2::WignerSeries) =
    same_range(s1, s2) && all(b1 == b2 for (b1, b2) ∈ zip(s1.blocks, s2.blocks))
Base.isequal(s1::WignerSeries, s2::WignerSeries) =
    same_range(s1, s2) && all(isequal(b1, b2) for (b1, b2) ∈ zip(s1.blocks, s2.blocks))
Base.isapprox(s1::WignerSeries, s2::WignerSeries; kwargs...) = same_range(s1, s2) &&
    all(isapprox(b1, b2; kwargs...) for (b1, b2) ∈ zip(s1.blocks, s2.blocks))
function Base.hash(s::WignerSeries, h::UInt)
    h = hash(:WignerSeries, hash(s.ℓₘᵢₙ, hash(s.ℓₘₐₓ, h)))
    foldl((h, b) -> hash(b, h), s.blocks; init=h)
end

function Base.show(io::IO, s::WignerSeries{IT}) where {IT}
    print(io, "WignerSeries{$IT} for ℓ ∈ $(s.ℓₘᵢₙ):$(s.ℓₘₐₓ)")
end
function Base.show(io::IO, ::MIME"text/plain", s::WignerSeries)
    show(io, s)
    println(io, ":")
    show_blocks(io, s, s.ℓₘᵢₙ, s.ℓₘₐₓ)
end

# The blocks of a series, or of a `HarmonicValues`, one after another with their ℓ.  Where
# the output is limited, as at the REPL, only the first two and the last two are printed
# when there are more than four, as `Base` elides the middle of a long array.
function show_blocks(io::IO, s, ℓₘᵢₙ, ℓₘₐₓ)
    n = Int(ℓₘₐₓ - ℓₘᵢₙ) + 1
    elide = get(io, :limit, false)::Bool && n > 4
    for (i, ℓ) ∈ enumerate(ℓₘᵢₙ:ℓₘₐₓ)
        if elide && 2 < i ≤ n - 2
            i == 3 && println(io, " ⋮")
            continue
        end
        println(io, " ℓ = ", ℓ, ":")
        show(io, MIME("text/plain"), s[ℓ])
        println(io)
    end
end


"""
    WignerDMatrix{IT, RT, ST}
    WignerDMatrix(parent, ℓ; m′ₘₐₓ=ℓ, m′ₘᵢₙ=-m′ₘₐₓ, mₘₐₓ=ℓ, mₘᵢₙ=-mₘₐₓ)
    WignerDMatrix(Complex{RT}, ℓ, m′ₘₐₓ=ℓ; mₘₐₓ=ℓ, mₘᵢₙ=-mₘₐₓ)

One block of Wigner's ``𝔇^{(ℓ)}_{m′,m}`` matrix: a [`WignerMatrix`](@ref) whose number type
is complex.  This is the type of each element of the [`WignerSeries`](@ref) that [`D`](@ref)
returns.

The first constructor wraps an existing `AbstractMatrix{Complex{RT}}`, which must be at
least `(m′ₘₐₓ-m′ₘᵢₙ+1) × (mₘₐₓ-mₘᵢₙ+1)`, exactly as [`WignerMatrix`](@ref) does.  The second
allocates uninitialized storage for the block with `m′ ∈ -m′ₘₐₓ:m′ₘₐₓ` and `m ∈ mₘᵢₙ:mₘₐₓ`.
Its range of `m′` is symmetric, as that of an [`HWedge`](@ref) is, so it takes the one limit
`m′ₘₐₓ`, positionally, where the first form takes both limits of `m′` as keywords.  `ℓ` and
the limits must all be integers of type `Int`, or all half-odd-integers, each a
[`HalfOddInteger`](@ref) or a `Rational{Int}` with denominator 2, and the keywords may also
be spelled `mp_max`, `mp_min`, `m_max` and `m_min` (the last two alone in the second form).

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
    WignerdMatrix(parent, ℓ; m′ₘₐₓ=ℓ, m′ₘᵢₙ=-m′ₘₐₓ, mₘₐₓ=ℓ, mₘᵢₙ=-mₘₐₓ)
    WignerdMatrix(RT, ℓ, m′ₘₐₓ=ℓ; mₘₐₓ=ℓ, mₘᵢₙ=-mₘₐₓ)

One block of Wigner's real ``d^{(ℓ)}_{m′,m}(β)`` matrix: a [`WignerMatrix`](@ref) whose
number type is real.  This is the type of each element of the [`WignerSeries`](@ref) that
[`d`](@ref) returns.  It is the real analogue of [`WignerDMatrix`](@ref) in every other
respect, including the constructors and being an alias rather than a distinct type.

# Example
```julia
w = WignerdMatrix(Float64, 3//2)   # uninitialized 4×4 block, m′, m ∈ -3//2:3//2
w[1//2, -1//2] = 0.25
```
"""
const WignerdMatrix{IT, RT, ST} = WignerMatrix{IT, RT, ST} where {IT, RT<:Real, ST<:AbstractMatrix{RT}}

# Constructors for WignerDMatrix (complex).  The storage-wrapping forms pass their limits on
# to `WignerMatrix`, which checks them; the storage of the allocating forms is sized by the
# limits after they have been checked, and `2m′ₘₐₓ + 1` and `mₘₐₓ - mₘᵢₙ + 1` are `Int`s for
# either kind of index.
@index_methods function WignerDMatrix(
    parent::AbstractMatrix{Complex{RT}}, ℓ::IndexType;
    mp_max::IndexType=ℓ, m′ₘₐₓ::IndexType=mp_max,
    mp_min::IndexType=-m′ₘₐₓ, m′ₘᵢₙ::IndexType=mp_min,
    m_max::IndexType=ℓ, mₘₐₓ::IndexType=m_max,
    m_min::IndexType=-mₘₐₓ, mₘᵢₙ::IndexType=m_min
) where {RT<:Real}
    WignerMatrix(parent, ℓ; m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
end
function WignerDMatrix(parent::AbstractMatrix{RT}, ℓ; kwargs...) where {RT<:Real}
    throw(ArgumentError(
        "WignerDMatrix only supports complex types; the input type is $RT.\n"
        * "Perhaps you meant to use WignerdMatrix?"
    ))
end
@index_methods function WignerDMatrix(
    ::Type{Complex{RT}}, ℓ::IT, m′ₘₐₓ::IT=ℓ;
    m_max::IndexType=ℓ, mₘₐₓ::IndexType=m_max,
    m_min::IndexType=-mₘₐₓ, mₘᵢₙ::IndexType=m_min
) where {RT<:Real, IT<:IndexType}
    validate_index_ranges(ℓ, m′ₘₐₓ, -m′ₘₐₓ, mₘₐₓ, mₘᵢₙ)
    parent = Matrix{Complex{RT}}(undef, 2m′ₘₐₓ + 1, (mₘₐₓ - mₘᵢₙ) + 1)
    WignerMatrix(parent, ℓ; m′ₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ, mₘₐₓ, mₘᵢₙ)
end

# Constructors for WignerdMatrix (real)
@index_methods function WignerdMatrix(
    parent::AbstractMatrix{RT}, ℓ::IndexType;
    mp_max::IndexType=ℓ, m′ₘₐₓ::IndexType=mp_max,
    mp_min::IndexType=-m′ₘₐₓ, m′ₘᵢₙ::IndexType=mp_min,
    m_max::IndexType=ℓ, mₘₐₓ::IndexType=m_max,
    m_min::IndexType=-mₘₐₓ, mₘᵢₙ::IndexType=m_min
) where {RT<:Real}
    WignerMatrix(parent, ℓ; m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
end
function WignerdMatrix(parent::AbstractMatrix{Complex{RT}}, ℓ; kwargs...) where {RT<:Real}
    throw(ArgumentError(
        "WignerdMatrix only supports real types; the input type is Complex{$RT}.\n"
        * "Perhaps you meant to use WignerDMatrix?"
    ))
end
@index_methods function WignerdMatrix(
    ::Type{RT}, ℓ::IT, m′ₘₐₓ::IT=ℓ;
    m_max::IndexType=ℓ, mₘₐₓ::IndexType=m_max,
    m_min::IndexType=-mₘₐₓ, mₘᵢₙ::IndexType=m_min
) where {RT<:Real, IT<:IndexType}
    validate_index_ranges(ℓ, m′ₘₐₓ, -m′ₘₐₓ, mₘₐₓ, mₘᵢₙ)
    parent = Matrix{RT}(undef, 2m′ₘₐₓ + 1, (mₘₐₓ - mₘᵢₙ) + 1)
    WignerMatrix(parent, ℓ; m′ₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ, mₘₐₓ, mₘᵢₙ)
end


@testitem "WignerMatrix" begin
    import SphericalFunctions: WignerDMatrix, WignerdMatrix,
        parent, ell, mp_max, mp_min, m_max, m_min, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, ℓₘᵢₙ,
        HalfOddInteger

    # Check that mixed-up types throw an error
    @test_throws ArgumentError WignerDMatrix(rand(Float64, 3, 3), 1)
    @test_throws ArgumentError WignerdMatrix(rand(ComplexF64, 3, 3), 1)
    @test_throws "WignerDMatrix only supports complex types" WignerDMatrix(rand(Float64, 3, 3), 1)
    @test_throws "WignerdMatrix only supports real types" WignerdMatrix(rand(ComplexF64, 3, 3), 1)
    @test_throws "WignerDMatrix only supports complex types" WignerDMatrix(rand(Float64, 2, 2), 1//2)
    @test_throws "WignerdMatrix only supports real types" WignerdMatrix(rand(ComplexF64, 2, 2), 1//2)

    # Check that a negative ℓ value throws an error
    @test_throws ArgumentError WignerDMatrix(rand(ComplexF64, 3, 3), -1)
    @test_throws "ℓₘₐₓ=-1 must be non-negative." WignerDMatrix(rand(ComplexF64, 3, 3), -1)
    @test_throws "ℓₘₐₓ=-1 must be non-negative." WignerdMatrix(rand(Float64, 3, 3), -1)
    @test_throws "ℓₘₐₓ=-1//2 must be non-negative." WignerDMatrix(rand(ComplexF64, 2, 2), -1//2)
    @test_throws "ℓₘₐₓ=-1//2 must be non-negative." WignerdMatrix(rand(Float64, 2, 2), -1//2)

    # A `Rational` ℓ that is not a half-odd-integer is refused before anything is built:
    # `1//3` has the wrong denominator, and `2//2`, which is `1//1`, is a whole number,
    # which is written as the integer 1.
    for W ∈ (WignerDMatrix, WignerdMatrix), n ∈ (2, 3)
        T = W === WignerDMatrix ? ComplexF64 : Float64
        @test_throws ArgumentError W(rand(T, n, n), 1//3)
        @test_throws "1//3 is neither an integer nor a half-odd-integer" W(rand(T, n, n), 1//3)
        @test_throws "1//1 is a whole number; write it as the integer 1" W(rand(T, n, n), 2//2)
    end

    ℓₘₐₓ = 2
    # Encode on twice-indices, so that the arithmetic is `Int` for both index types (a
    # `HalfOddInteger` may only be multiplied by an even integer).  `2x + 6` is in `1:11`
    # for every index used here, so base 25 keeps the encoding injective.
    code(x) = 2x + 6
    encode(ℓ, m′, m) = code(ℓ) + code(m′)*25 + code(m)*625
    for ℓ ∈ Any[collect(0:ℓₘₐₓ); HalfOddInteger.(collect(1//2:(ℓₘₐₓ+1//2)))]
        # The input must be at least as big as the block in each dimension, and may be bigger
        @test_throws "The extent of the first dimension" WignerDMatrix(Array{ComplexF64}(undef, 2ℓ, 2ℓ + 1), ℓ)
        @test_throws "The extent of the first dimension" WignerdMatrix(Array{Float64}(undef, 2ℓ, 2ℓ + 1), ℓ)
        @test_throws "The extent of the second dimension" WignerDMatrix(Array{ComplexF64}(undef, 2ℓ + 1, 2ℓ), ℓ)
        @test_throws "The extent of the second dimension" WignerdMatrix(Array{Float64}(undef, 2ℓ + 1, 2ℓ), ℓ)
        @test WignerDMatrix(Array{ComplexF64}(undef, 2ℓ + 2, 2ℓ + 3), ℓ) isa WignerDMatrix

        # Check that a data array with a dimension of 0 extent throws an error.
        @test_throws r"The extent of the second dimension.*; it is 0." WignerDMatrix(Array{ComplexF64}(undef, 2ℓ + 1, 0), ℓ)
        @test_throws r"The extent of the first dimension.*; it is 0." WignerDMatrix(Array{ComplexF64}(undef, 0, 2ℓ + 1), ℓ)
        @test_throws r"The extent of the second dimension.*; it is 0." WignerdMatrix(Array{Float64}(undef, 2ℓ + 1, 0), ℓ)
        @test_throws r"The extent of the first dimension.*; it is 0." WignerdMatrix(Array{Float64}(undef, 0, 2ℓ + 1), ℓ)

        # Every symmetric restriction of either axis, with the limits given in both
        # spellings
        for m′ₘ ∈ ℓₘᵢₙ(ℓ):ℓ, mₘ ∈ ℓₘᵢₙ(ℓ):ℓ
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
                @test mp_max(w) == m′ₘ
                @test m_max(w) == mₘ
                @test mp_min(w) == -mp_max(w)
                @test m_min(w) == -m_max(w)
                for m ∈ -mₘ:mₘ
                    for m′ ∈ -m′ₘ:m′ₘ
                        @test w[m′, m] == encode(ℓ, m′, m)
                    end
                end
                # The ASCII spellings of the keywords, and the default lower limits, which
                # are minus the upper ones, give the same block
                @test WignerMatrixType(NT.(data), ℓ; mp_max=m′ₘ, mp_min=-m′ₘ, m_max=mₘ, m_min=-mₘ) == w
                @test WignerMatrixType(NT.(data), ℓ; m′ₘₐₓ=m′ₘ, mₘₐₓ=mₘ) == w
                @test WignerMatrixType(NT.(data), ℓ; mp_max=m′ₘ, m_max=mₘ) == w
            end
        end

        for m′ₘ ∈ ℓₘᵢₙ(ℓ):ℓ
            for WignerMatrixType ∈ (WignerDMatrix, WignerdMatrix)
                data = rand(
                    WignerMatrixType<:WignerDMatrix ? ComplexF64 : Float64,
                    2m′ₘ + 1, 2ℓ + 1
                )
                w = WignerMatrixType(data, ℓ; m′ₘₐₓ=m′ₘ, m′ₘᵢₙ=-m′ₘ)

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
# indexed, and `natural_axes` gives the matching ranges.
#
# `getindex` and `setindex!` are deliberately *not* written this way, and neither are the
# accessors `ℓ`, `mₘₐₓ`, `sₘᵢₙ` and the rest.  Those are the hot path: the bounds check
# there reads the stored limits directly rather than building a range to test membership in,
# which is the difference the note above `inrange` measures at 145 ns against 0.9 ns per
# element.
#
# `HWedge` and `HAxis` are outside this union on purpose.  They are workspaces for the
# recursion rather than blocks handed to a caller, their storage is triangular, and `size`
# of one is the size of that storage rather than of any block; they keep the generic
# `AbstractWignerMatrix` methods above.

const BlockContainer = Union{
    WignerMatrix, WignerMatrixBatch,
    DegreeBlock, DegreeBlockBatch,
    SpinMatrix, SpinMatrixBatch,
}

# The roles are a property of the type, so these fold away at compile time.  A batched
# container's leading axis is an ordinary 1-based rotor position, not a natural index, which
# is why it is the one role whose range is a plain `UnitRange` rather than a `WignerRange`.
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

# The extent is taken as the difference of the endpoints rather than as `length` of the
# range, because that is an `Int` by construction for a half-odd-integer axis, where
# `Int(ℓ)` throws.
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
# The throwing form, as for an array.
function Base.checkbounds(w::BlockContainer, I...)
    checkbounds(Bool, w, I...) || throw(BoundsError(w, I))
    nothing
end

isbatched(::Union{WignerMatrixBatch, DegreeBlockBatch, SpinMatrixBatch}) = true
isbatched(::Union{WignerMatrix, DegreeBlock, SpinMatrix}) = false

# The blocks with a single m axis have no m′, and are told so rather than being shown an
# error about a missing field.
const OneAxisBlock = Union{DegreeBlock, DegreeBlockBatch, SpinMatrix, SpinMatrixBatch}
m′ₘₐₓ(w::OneAxisBlock) = throw(no_m′_axis(w))
m′ₘᵢₙ(w::OneAxisBlock) = throw(no_m′_axis(w))
@noinline no_m′_axis(w) = ArgumentError(
    "A `$(nameof(typeof(w)))` has no m′ axis; its axes are "
    * join(("`$r`" for r ∈ axis_roles(w)), ", ") * ", and the limits of the natural ones "
    * "are read with `mₘₐₓ` and `mₘᵢₙ`, and for a spin axis `sₘₐₓ` and `sₘᵢₙ`."
)

Base.Array(w::BlockContainer) = Array(stored_elements(w))
Base.collect(w::BlockContainer) = Array(w)

# The labels of a block are its kind — the roles of its axes — its ℓ, and the limits of its
# axes, including the number of rotors of a batch.  The same numbers under different labels
# are not the same block, so `==`, `isequal`, `≈` and `hash` count the labels as well as the
# numbers, which they read from the leading block of the storage without copying it.  Blocks
# of the same labels but different number types compare their numbers as arrays do, and hash
# as arrays do, so that `hash` agrees with `isequal` between them too.
block_labels(w::BlockContainer) =
    (axis_roles(w), ℓ(w), map(a -> (first(a), last(a)), natural_axes(w)))
Base.:(==)(a::BlockContainer, b::BlockContainer) =
    block_labels(a) == block_labels(b) && stored_elements(a) == stored_elements(b)
Base.isequal(a::BlockContainer, b::BlockContainer) =
    isequal(block_labels(a), block_labels(b)) && isequal(stored_elements(a), stored_elements(b))
Base.isapprox(a::BlockContainer, b::BlockContainer; kwargs...) =
    block_labels(a) == block_labels(b) &&
    isapprox(stored_elements(a), stored_elements(b); kwargs...)
Base.hash(w::BlockContainer, h::UInt) =
    hash(stored_elements(w), hash(block_labels(w), hash(:BlockContainer, h)))
