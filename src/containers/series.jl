### Series of blocks, indexed by ℓ.
#
# A `WignerSeries` holds the blocks of a Wigner matrix for a range of ℓ, and a
# `HarmonicValues` the values of the spin-weighted harmonics, whose blocks are views of one
# flat array.  Both are indexed, iterated and counted by ℓ, as the calculators are, and that
# interface is written once, on their supertype, in terms of the `ℓₘᵢₙ`, `ℓₘₐₓ` and
# `getindex(s, ℓ)` that each supplies.

"""
    AbstractDegreeSeries{IT}

Supertype of the series of blocks indexed by the degree ``ℓ``: [`WignerSeries`](@ref), the
blocks of a Wigner matrix, and [`HarmonicValues`](@ref), the values of the spin-weighted
spherical harmonics.
- `IT` is the index type, `Int` or [`HalfOddInteger`](@ref).

A series `s` holds one block for each ``ℓ`` from `ℓₘᵢₙ(s)` to `ℓₘₐₓ(s)`, and `s[ℓ]` is the
block of degree ``ℓ``; for a half-integer series `ℓ` may be written as a
[`HalfOddInteger`](@ref) or as a `Rational{Int}` with denominator 2, and for an integer
series it is an `Int`.  Like a calculator, a series iterates as `ℓ => block` pairs, so
`eltype(s)` is that pair type; `keys(s)` is the range of ``ℓ``, and `length(s)` counts the
blocks.  `first(s)` and `last(s)` are the first and last blocks, as indexing gives them, and
so are `first(s, n)` and `last(s, n)`, which give vectors of the first or last `n` blocks,
and `only(s)`, the one block of a series that has only one.  [`ishalfinteger`](@ref) says
which kind of index the series takes.

A series is not an `AbstractVector`: ``ℓ`` may be a half-odd-integer, and the `axes` of an
`AbstractVector` must be integer ranges.
"""
abstract type AbstractDegreeSeries{IT<:IntegerHalf} end

Base.keys(s::AbstractDegreeSeries) = ℓₘᵢₙ(s):ℓₘₐₓ(s)
Base.length(s::AbstractDegreeSeries) = Int(ℓₘₐₓ(s) - ℓₘᵢₙ(s)) + 1
Base.firstindex(s::AbstractDegreeSeries) = ℓₘᵢₙ(s)
Base.lastindex(s::AbstractDegreeSeries) = ℓₘₐₓ(s)
ishalfinteger(::AbstractDegreeSeries{IT}) where {IT<:Integer} = false
ishalfinteger(::AbstractDegreeSeries{IT}) where {IT<:HalfOddInteger} = true

# A series iterates as `ℓ => block` pairs, as a calculator does, so that a loop over
# `D(R, ℓₘₐₓ)` reads exactly as one over `DCalculator(R, ℓₘₐₓ)`; the element type is
# therefore that of the pairs, rather than the number type of the blocks, and a disagreement
# between the two makes `collect` throw.  Indexing, `first` and `last` give blocks, as
# `s[ℓ]` does, and so do the forms of `first` and `last` that take a count, which `Base`
# would otherwise derive from the iteration, as pairs.
@inline function Base.iterate(s::AbstractDegreeSeries{IT}, ℓ::IT=ℓₘᵢₙ(s)) where {IT}
    ℓ > ℓₘₐₓ(s) && return nothing
    (ℓ => s[ℓ], ℓ + 1)
end
Base.IteratorSize(::Type{<:AbstractDegreeSeries}) = Base.HasLength()
Base.IteratorEltype(::Type{<:AbstractDegreeSeries}) = Base.HasEltype()
Base.pairs(s::AbstractDegreeSeries) = s
Base.eltype(::Type{S}) where {IT, S<:AbstractDegreeSeries{IT}} =
    Pair{IT, Base.promote_op(getindex, S, IT)}
Base.eltype(s::AbstractDegreeSeries) = eltype(typeof(s))
Base.first(s::AbstractDegreeSeries) = s[ℓₘᵢₙ(s)]
Base.last(s::AbstractDegreeSeries) = s[ℓₘₐₓ(s)]
function Base.first(s::AbstractDegreeSeries, n::Integer)
    n < 0 && throw(ArgumentError("Number of elements must be non-negative"))
    [s[ℓₘᵢₙ(s) + (i - 1)] for i ∈ 1:min(n, length(s))]
end
function Base.last(s::AbstractDegreeSeries, n::Integer)
    n < 0 && throw(ArgumentError("Number of elements must be non-negative"))
    k = min(n, length(s))
    [s[ℓₘₐₓ(s) - (k - i)] for i ∈ 1:k]
end

# The degree `ℓ` asked of a series, converted to the series' own index type as the index
# methods convert, so that an index of the wrong kind or type is told what the series takes,
# and then checked to be one of the degrees the series holds.  Indexing a `HarmonicValues`
# begins with this.
@inline function series_ℓ(s::AbstractDegreeSeries{IT}, ℓ) where {IT}
    ℓ′ = checked_index(IT, ℓ, s, "ℓ")
    if ℓ′ < ℓₘᵢₙ(s) || ℓ′ > ℓₘₐₓ(s)
        throw(BoundsError(s, ℓ))
    end
    ℓ′
end

function Base.show(io::IO, ::MIME"text/plain", s::AbstractDegreeSeries)
    show(io, s)
    println(io, ":")
    show_blocks(io, s)
end

# The blocks of a series, one after another with their ℓ.  Where the output is limited, as
# at the REPL, only the first two and the last two are printed when there are more than
# four, as `Base` elides the middle of a long array.
function show_blocks(io::IO, s::AbstractDegreeSeries)
    n = Int(ℓₘₐₓ(s) - ℓₘᵢₙ(s)) + 1
    elide = get(io, :limit, false)::Bool && n > 4
    for (i, ℓ) ∈ enumerate(ℓₘᵢₙ(s):ℓₘₐₓ(s))
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
    WignerSeries{IT, VT}
    WignerSeries(blocks, ℓₘᵢₙ, ℓₘₐₓ)

The blocks of a Wigner matrix for every ``ℓ`` from `ℓₘᵢₙ` to `ℓₘₐₓ`, indexed by ``ℓ``:
`s[ℓ]` is the block of degree `ℓ`, and `s[ℓ][m′, m]` an element of it.  For a half-integer
series `ℓ` may be written as a [`HalfOddInteger`](@ref) or as a `Rational{Int}` with
denominator 2, and for an integer series it is an `Int`.
- `IT` is the index type, `Int` or [`HalfOddInteger`](@ref).
- `VT` is the type of the vector of blocks.

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
struct WignerSeries{IT<:IntegerHalf, VT<:AbstractVector} <: AbstractDegreeSeries{IT}
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

@inline function check_series_block(block::AbstractBlock{IT}, i, ℓₘᵢₙ::IT) where {IT}
    ℓ(block) == ℓₘᵢₙ + (i - 1) || throw(series_block_error(block, i, ℓₘᵢₙ))
    nothing
end
check_series_block(block, i, ℓₘᵢₙ) = throw(series_block_error(block, i, ℓₘᵢₙ))
@noinline function series_block_error(block, i, ℓₘᵢₙ)
    if !(block isa AbstractBlock)
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
Base.axes(s::WignerSeries) = (WignerRange(s.ℓₘᵢₙ:s.ℓₘₐₓ),)
Base.axes(s::WignerSeries, d::Integer) = d ≤ 1 ? axes(s)[d] : Base.OneTo(1)
Base.ndims(::WignerSeries) = 1
Base.ndims(::Type{<:WignerSeries}) = 1
Base.size(s::WignerSeries) = (length(s),)
Base.size(s::WignerSeries, d::Integer) = d ≤ 1 ? length(s) : 1
Base.values(s::WignerSeries) = s.blocks

# The blocks of a `WignerSeries` are held in a vector, and these methods read that vector
# directly, where the methods of every series index it ℓ by ℓ.  They give the same blocks
# while the vector has the length that the labels give, which `check_blocks` below confirms
# before any block is read; `length` counts the blocks in the vector, `first(s, n)` and
# `last(s, n)` return vectors of its element type, and `only` refuses a series of other than
# one block with the message `Base` gives for a vector.  The iteration is over the series'
# own position, rather than handing an integer state to the storage's own `iterate`: for
# storage other than a `Vector`, such as a view, that state is not an integer, and a series
# on a view would stop after its first block while its `length` still counted them all — so
# that a comprehension over it would return uninitialized memory.
Base.length(s::WignerSeries) = length(s.blocks)
function Base.iterate(s::WignerSeries, i::Int=1)
    i == 1 && check_blocks(s)
    i > length(s.blocks) && return nothing
    ((s.ℓₘᵢₙ + (i - 1)) => s.blocks[i], i + 1)
end
Base.first(s::WignerSeries, n::Integer) = (check_blocks(s); first(s.blocks, n))
Base.last(s::WignerSeries, n::Integer) = (check_blocks(s); last(s.blocks, n))
Base.only(s::WignerSeries) = (check_blocks(s); only(s.blocks))

@propagate_inbounds function Base.getindex(s::WignerSeries{IT}, ℓ) where {IT}
    # Deliberately *not* inside `@boundscheck`: an index of the wrong kind, such as a whole
    # number asked of a half-integer series, is told what the series takes.
    let ℓ = checked_index(IT, ℓ, s, "ℓ")
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


### Harmonic values, laid out in the canonical mode ordering.
#
# Two things in this package are stored as `[x(ℓ, m) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ]` — the
# weights of a spin-weighted function, and the values of the spin-weighted harmonics
# themselves — and they share the flat storage, the labels ℓₘᵢₙ and ℓₘₐₓ, and the block
# accessor `x[ℓ, :]`.  They differ in meaning, and so in the rest of their interfaces: a
# `ModeWeights` is the vector of the weights, indexed, iterated and counted by mode, while a
# `HarmonicValues` is a series, indexed, iterated and counted by ℓ, as the calculators are.
#
# The mode axis is always the *last* axis of the storage, so that a single ``ℓ`` is a view over
# a contiguous run of it with every leading axis taken whole.  That is what keeps the flat form
# usable for the products these containers exist to feed: a synthesis matrix times a vector of
# mode weights.

"""
    HarmonicValues

The values of the spin-weighted spherical harmonics ``{}_sY_{ℓ,m}`` at one or more rotors,
indexed first by ``ℓ`` and then naturally within the block:

| built for | `Y[ℓ]` is indexed |
|---|---|
| one rotor, one spin weight | `[m]` |
| many rotors, one spin weight | `[iᵣ, m]` |
| one rotor, a range of spin weights | `[s, m]` |
| many rotors, a range of spin weights | `[iᵣ, s, m]` |

The parameters of `HarmonicValues{T, IT, S, A}` are as follows:
- `T` is the number type.
- `IT` is the index type, `Int` or [`HalfOddInteger`](@ref).
- `S` is the type of the spin weights: `IT` for one, or a range of `IT` for several.
- `A` is the type of the storage, an array of `T`.

This is what [`sYlm`](@ref) returns.  The blocks are [`DegreeBlock`](@ref),
[`DegreeBlockBatch`](@ref), [`SpinMatrix`](@ref) and [`SpinMatrixBatch`](@ref) respectively,
and are views into the storage rather than copies, so writing through one writes into `Y`.

[`array_view`](@ref) gives the flat storage: a `Vector` of modes, or an array whose *last* axis is
the modes in the canonical ordering (see [`Yindex`](@ref)) and whose leading axes are the
rotors and spin weights.  That is the form a product with mode weights takes, and
[`sYlm_matrix`](@ref) is the direct name for it.

`spins(Y)` is the range of spin weights served, `spin(Y)` the single value when there is only
one, `Nᵣ(Y)` the number of rotors, [`isbatched`](@ref)`(Y)` whether there is a rotor axis, and
`ℓₘᵢₙ(Y)`/`ℓₘₐₓ(Y)` the range of ``ℓ``.

`Y[ℓ]` and `Y[ℓ, :]` both give the block of degree ``ℓ``, the second as `w[ℓ, :]` gives the
block of a [`ModeWeights`](@ref) `w`; for a half-integer `Y`, `ℓ` may be written as a
[`HalfOddInteger`](@ref) or as a `Rational{Int}` with denominator 2, and for an integer `Y` it
is an `Int`.  `first(Y)` and `last(Y)` are the first and last blocks, as indexing gives
them, and so are `first(Y, n)`, `last(Y, n)` and `only(Y)`, as for a
[`WignerSeries`](@ref).  Iterating gives `ℓ => block` pairs, as a calculator does, so
`eltype(Y)` is that pair type; the number type is `eltype(array_view(Y))`.  `length(Y)`
counts the blocks, and `keys(Y)` is the range of ``ℓ``.

`HarmonicValues(data, s, ℓₘᵢₙ, ℓₘₐₓ, Nᵣ)` wraps existing storage of one of those four
shapes without copying it.  Its labels must describe the storage: the rank must be that of
one of the shapes, the spin axis as long as the range `s`, the leading axis `Nᵣ` long
(`Nᵣ = 1` without it), and the mode axis `Ysize(ℓₘᵢₙ, ℓₘₐₓ)` long.  The indices must all be
integers of type `Int`, or all half-odd-integers, each a [`HalfOddInteger`](@ref) or a
`Rational{Int}` with denominator 2.

`copy`, `similar` and [`relabel`](@ref) keep the labels.  As for a `ModeWeights`, two
`HarmonicValues` are equal, or approximately equal, only when their labels agree as well as
their numbers, while `==`, `isequal` and `≈` against a plain array compare the numbers of
`array_view(Y)`; `hash` is therefore that of the numbers alone.  Broadcasting reads the
numbers of `array_view(Y)` and gives a plain array, and `Y .= x` writes into the storage; a
broadcast that combines two `HarmonicValues`, or writes one into another, requires their labels
to agree, since it pairs the numbers by position.

See also [`ModeWeights`](@ref), which shares this layout but holds the weights of a function
rather than the values of the harmonics.
"""
struct HarmonicValues{T, IT<:IntegerHalf, S, A<:AbstractArray{T}} <: AbstractDegreeSeries{IT}
    data::A
    s::S          # one `IT`, or an ascending range of them
    ℓₘᵢₙ::IT
    ℓₘₐₓ::IT
    Nᵣ::Int       # 1 when the container was built for a single rotor

    # The labels are checked against the storage here, once, because everything else reads
    # them in place of the storage's own shape: the rank of the storage says whether there is
    # a rotor axis and a spin axis, which must agree with `Nᵣ` and with the spin weights, and
    # the index kind of the spin weights must be that of the ℓ range, which the index
    # methods ensure.  Every block and product relies on these.
    @index_methods function HarmonicValues(
        data::A, s::IndexOrRange, ℓₘᵢₙ::IT, ℓₘₐₓ::IT, Nᵣ::Int
    ) where {T, IT<:IndexType, A<:AbstractArray{T}}
        Base.require_one_based_indexing(data)
        check_harmonic_storage(data, s, Nᵣ)
        if size(data)[end] != Ysize(ℓₘᵢₙ, ℓₘₐₓ)
            throw(ArgumentError(
                "The mode axis has length $(size(data)[end]), but "
                * "Ysize(ℓₘᵢₙ=$ℓₘᵢₙ, ℓₘₐₓ=$ℓₘₐₓ) = $(Ysize(ℓₘᵢₙ, ℓₘₐₓ))."
            ))
        end
        new{T, IT, typeof(s), A}(data, s, ℓₘᵢₙ, ℓₘₐₓ, Nᵣ)
    end
end

# The storage is `[modes]` or `[iᵣ, modes]` for one spin weight, and `[s, modes]` or
# `[iᵣ, s, modes]` for a range of them; the leading rotor axis holds `Nᵣ` rotors, and is
# present exactly when the values were computed for a vector of rotors, even one of length 1.
function check_harmonic_storage(data::AbstractArray, s, Nᵣ::Int)
    spin_axes = s isa AbstractUnitRange ? 1 : 0
    N = ndims(data)
    if !(N == spin_axes + 1 || N == spin_axes + 2)
        throw(DimensionMismatch(
            "The storage has $N dimensions, but harmonic values for "
            * (s isa AbstractUnitRange ? "the range of spin weights $s" : "one spin weight")
            * " are stored with $(spin_axes + 1), or $(spin_axes + 2) with a leading rotor axis."
        ))
    end
    batched = N == spin_axes + 2
    if s isa AbstractUnitRange
        if isempty(s)
            throw(ArgumentError(
                "The range of spin weights $s is empty; it runs from its lower limit to its "
                * "upper one."
            ))
        end
        if size(data, N - 1) != length(s)
            throw(DimensionMismatch(
                "The spin axis of the storage has length $(size(data, N - 1)), but the range "
                * "of spin weights $s has $(length(s))."
            ))
        end
    end
    if batched ? Nᵣ != size(data, 1) : Nᵣ != 1
        throw(DimensionMismatch(
            batched ?
            "The storage has a rotor axis of length $(size(data, 1)), but Nᵣ=$Nᵣ." :
            "The storage has no rotor axis, so it holds the values at one rotor, but Nᵣ=$Nᵣ."
        ))
    end
    nothing
end

ℓₘᵢₙ(Y::HarmonicValues) = Y.ℓₘᵢₙ
ℓₘₐₓ(Y::HarmonicValues) = Y.ℓₘₐₓ
Base.parent(Y::HarmonicValues) = Y.data
Nᵣ(Y::HarmonicValues) = Y.Nᵣ
# Batched when the storage has a rotor axis — one more dimension than the modes (and the spin
# weights, if there are several) need — however many rotors it holds, so that this agrees with
# the blocks, whose type is decided by the same dimensions.
isbatched(Y::HarmonicValues) = ndims(Y.data) == (Y.s isa AbstractUnitRange ? 3 : 2)
spins(Y::HarmonicValues{T, IT, S}) where {T, IT, S<:IntegerHalf} = Y.s:Y.s
spins(Y::HarmonicValues{T, IT, S}) where {T, IT, S<:AbstractUnitRange} = Y.s
spin(Y::HarmonicValues{T, IT, S}) where {T, IT, S<:IntegerHalf} = Y.s

# `length(Y)` counts the blocks, as it does for every series; the number of modes is
# `length(array_view(Y))` for the unbatched single-spin case, and `Ysize` in general.
# `only` gives the one block, as indexing does.
function Base.only(Y::HarmonicValues)
    length(Y) == 1 || throw(ArgumentError(
        "These harmonic values hold $(length(Y)) blocks, for ℓ ∈ $(ℓₘᵢₙ(Y)):$(ℓₘₐₓ(Y)), "
        * "rather than exactly one."
    ))
    Y[ℓₘᵢₙ(Y)]
end

# The positions in the flat storage that one ℓ occupies.  `Yindex` counts from `ℓₘᵢₙ`, so
# this is the same arithmetic for either kind of index.  Every caller has converted `ℓ` to
# the container's own index type, `Int` or `HalfOddInteger`, and for either of those `2ℓ` is
# an `Int`; the signature insists on it, so that an index of another type is a
# `MethodError` rather than a range of another type.  `ModeWeights` has the same method, in
# `mode_weights.jl`.
@inline function mode_range(Y::HarmonicValues{T, IT}, ℓ::IT) where {T, IT}
    i₀ = Yindex(ℓ, -ℓ, ℓₘᵢₙ(Y))
    i₀:(i₀ + 2ℓ)
end

# The four shapes.  Which one applies is fixed by the rank of the storage and by whether `S` is
# a single spin weight or a range, so each of these has a single concrete return type.
@propagate_inbounds function Base.getindex(
    Y::HarmonicValues{T, IT, S, <:AbstractVector}, ℓ
) where {T, IT, S<:IntegerHalf}
    let ℓ = series_ℓ(Y, ℓ)
        DegreeBlock(view(Y.data, mode_range(Y, ℓ)), ℓ)
    end
end
@propagate_inbounds function Base.getindex(
    Y::HarmonicValues{T, IT, S, <:AbstractMatrix}, ℓ
) where {T, IT, S<:IntegerHalf}
    let ℓ = series_ℓ(Y, ℓ)
        DegreeBlockBatch(view(Y.data, :, mode_range(Y, ℓ)), ℓ)
    end
end
@propagate_inbounds function Base.getindex(
    Y::HarmonicValues{T, IT, S, <:AbstractMatrix}, ℓ
) where {T, IT, S<:AbstractUnitRange}
    let ℓ = series_ℓ(Y, ℓ), sr = Y.s
        SpinMatrix(
            view(Y.data, :, mode_range(Y, ℓ)), ℓ;
            sₘₐₓ=last(sr), sₘᵢₙ=first(sr), mₘₐₓ=ℓ, mₘᵢₙ=-ℓ
        )
    end
end
@propagate_inbounds function Base.getindex(
    Y::HarmonicValues{T, IT, S, <:AbstractArray{T, 3}}, ℓ
) where {T, IT, S<:AbstractUnitRange}
    let ℓ = series_ℓ(Y, ℓ), sr = Y.s
        SpinMatrixBatch(
            view(Y.data, :, :, mode_range(Y, ℓ)), ℓ;
            sₘₐₓ=last(sr), sₘᵢₙ=first(sr), mₘₐₓ=ℓ, mₘᵢₙ=-ℓ
        )
    end
end

# `Y[ℓ, :]` is the same block, so that the block of one ℓ is written the same way for both mode
# containers, as `w[ℓ, :]` is for a `ModeWeights`.
@propagate_inbounds Base.getindex(Y::HarmonicValues, ℓ, ::Colon) = Y[ℓ]

Base.copy(Y::HarmonicValues) = HarmonicValues(copy(Y.data), Y.s, Y.ℓₘᵢₙ, Y.ℓₘₐₓ, Y.Nᵣ)
Base.similar(Y::HarmonicValues) = HarmonicValues(similar(Y.data), Y.s, Y.ℓₘᵢₙ, Y.ℓₘₐₓ, Y.Nᵣ)
Base.similar(Y::HarmonicValues, ::Type{S}) where {S} =
    HarmonicValues(similar(Y.data, S), Y.s, Y.ℓₘᵢₙ, Y.ℓₘₐₓ, Y.Nᵣ)

# The labels say what the numbers are the values of — which spin weights, which ℓ, how many
# rotors — so, as for `ModeWeights`, the same numbers under different labels are not equal,
# and not approximately equal either.  Against a plain array only the numbers can be compared,
# and they are compared with `array_view(Y)`, whatever its shape.
function same_labels(a::HarmonicValues, b::HarmonicValues)
    a.s == b.s && a.ℓₘᵢₙ == b.ℓₘᵢₙ && a.ℓₘₐₓ == b.ℓₘₐₓ && a.Nᵣ == b.Nᵣ
end
Base.:(==)(a::HarmonicValues, b::HarmonicValues) = same_labels(a, b) && a.data == b.data
Base.:(==)(Y::HarmonicValues, A::AbstractArray) = Y.data == A
Base.:(==)(A::AbstractArray, Y::HarmonicValues) = A == Y.data
Base.isequal(a::HarmonicValues, b::HarmonicValues) = same_labels(a, b) && isequal(a.data, b.data)
Base.isequal(Y::HarmonicValues, A::AbstractArray) = isequal(Y.data, A)
Base.isequal(A::AbstractArray, Y::HarmonicValues) = isequal(A, Y.data)
Base.isapprox(a::HarmonicValues, b::HarmonicValues; kwargs...) =
    same_labels(a, b) && isapprox(a.data, b.data; kwargs...)
Base.isapprox(Y::HarmonicValues, A::AbstractArray; kwargs...) = isapprox(Y.data, A; kwargs...)
Base.isapprox(A::AbstractArray, Y::HarmonicValues; kwargs...) = isapprox(A, Y.data; kwargs...)
# `isequal(Y, array_view(Y))` holds, so the hash must be that of the numbers alone; values
# under different labels then share a hash, which is allowed.
Base.hash(Y::HarmonicValues, h::UInt) = hash(Y.data, h)

# A batch of one rotor is still a batch, whose blocks have a rotor axis, so the count is shown
# for every batch.
function Base.show(io::IO, Y::HarmonicValues{T, IT, S}) where {T, IT, S}
    spin_text = S <: AbstractUnitRange ? "s ∈ $(Y.s)" : "s=$(Y.s)"
    rotor_text = isbatched(Y) ? ", $(Y.Nᵣ) rotor" * (Y.Nᵣ == 1 ? "" : "s") : ""
    print(io, "HarmonicValues{$T} for ℓ ∈ $(Y.ℓₘᵢₙ):$(Y.ℓₘₐₓ), $spin_text$rotor_text")
end
