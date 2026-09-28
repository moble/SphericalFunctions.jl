### Containers laid out in the canonical mode ordering.
#
# Two things in this package are stored as `[x(ℓ, m) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ]` — the
# weights of a spin-weighted function, and the values of the spin-weighted harmonics
# themselves — and they share the flat storage, the labels ℓₘᵢₙ and ℓₘₐₓ, and the block
# accessor `x[ℓ, :]`.  They differ in meaning, and so in the rest of their interfaces: a
# `ModeWeights` is the vector of the weights, indexed, iterated and counted by mode, while a
# `HarmonicValues` is indexed, iterated and counted by ℓ, as the calculators are.  The shared
# supertype is what lets the machinery that depends only on the layout be written once.
#
# The mode axis is always the *last* axis of the storage, so that a single ``ℓ`` is a view over
# a contiguous run of it with every leading axis taken whole.  That is what keeps the flat form
# usable for the products these containers exist to feed: a synthesis matrix times a vector of
# mode weights.

"""
    AbstractModeContainer{T, IT}

Supertype of the containers stored in the canonical mode ordering — [`ModeWeights`](@ref) and
[`HarmonicValues`](@ref).  `T` is the number type and `IT` the index type (`Int` or
[`HalfOddInteger`](@ref)).

The supertype promises only what the layout determines: the labels `ℓₘᵢₙ(c)` and `ℓₘₐₓ(c)`,
[`ishalfinteger`](@ref), the block of one ``ℓ`` as `c[ℓ, :]`, and the flat storage as
[`array_view`](@ref)`(c)`.  The rest of the interface differs between the two, because they
mean different things.  A `ModeWeights` is the vector of the weights of a function: `w[i]` is
the `i`-th weight in storage order, `length(w)` counts the modes, `keys(w)` are the linear
positions, and iteration yields the weights, so that `sum`, `maximum` and the other reductions
see numbers.  A `HarmonicValues` is indexed by ``ℓ``, as a calculator is: `Y[ℓ]` is the block
of degree ``ℓ``, `length(Y)` counts the blocks, `keys(Y)` is the range of ``ℓ``, and iteration
yields `ℓ => block` pairs.  `eltype` is the number type of a `ModeWeights` and that pair type
for a `HarmonicValues`.

These are *not* `AbstractArray`s.  A container indexed by ``ℓ`` cannot be one, because ``ℓ``
may be a half-odd-integer and `axes` must be integer ranges; and for the ones that could be,
being an array would let `*` and `mul!` accept them, and those return silently wrong answers
for an array with non-trivial offsets, such as an `OffsetArray`.  [`array_view`](@ref) is the
explicit route to the flat 1-based storage,
and is what the transforms and the operator matrices take.
"""
abstract type AbstractModeContainer{T, IT<:IntegerHalf} end

Base.eltype(::AbstractModeContainer{T}) where {T} = T
Base.eltype(::Type{<:AbstractModeContainer{T}}) where {T} = T
ℓₘᵢₙ(c::AbstractModeContainer) = c.ℓₘᵢₙ
ℓₘₐₓ(c::AbstractModeContainer) = c.ℓₘₐₓ
ishalfinteger(::AbstractModeContainer{T, IT}) where {T, IT<:Integer} = false
ishalfinteger(::AbstractModeContainer{T, IT}) where {T, IT<:HalfOddInteger} = true

# The positions in the flat storage that one ℓ occupies.  `Yindex` counts from `ℓₘᵢₙ`, so this
# is the same arithmetic for either kind of index.  Every caller has converted `ℓ` to the
# container's own index type, `Int` or `HalfOddInteger`, and for either of those `2ℓ` is an
# `Int`; the signature insists on it, so that an index of another type is a `MethodError`
# rather than a range of another type.
@inline function mode_range(c::AbstractModeContainer{T, IT}, ℓ::IT) where {T, IT}
    i₀ = Yindex(ℓ, -ℓ, ℓₘᵢₙ(c))
    i₀:(i₀ + 2ℓ)
end

# Shared by `getindex(c, ℓ)` on every such container: `ℓ` is converted to the container's own
# index type as the index methods convert, so that an index of the wrong kind or type is told
# what the container takes, and must then be one of the values the container holds.
@inline function check_ℓ(c::AbstractModeContainer{T, IT}, ℓ) where {T, IT}
    ℓ′ = container_index(IT, ℓ, c, "ℓ")
    if ℓ′ < ℓₘᵢₙ(c) || ℓ′ > ℓₘₐₓ(c)
        throw(BoundsError(c, ℓ))
    end
    ℓ′
end


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
struct HarmonicValues{T, IT<:IntegerHalf, S, A<:AbstractArray{T}} <: AbstractModeContainer{T, IT}
    data::A
    s::S          # one `IT`, or an ascending range of them
    ℓₘᵢₙ::IT
    ℓₘₐₓ::IT
    Nᵣ::Int       # 1 when the container was built for a single rotor

    # The labels are checked against the storage here, once, because everything else reads
    # them in place of the storage's own shape: the rank of the storage says whether there is
    # a rotor axis and a spin axis, which must agree with `Nᵣ` and with the spin weights, and
    # the index kind of the spin weights must be that of the ℓ range, which the index
    # methods ensure.  Every block, product and refill relies on these.
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

Base.parent(Y::HarmonicValues) = Y.data
Nᵣ(Y::HarmonicValues) = Y.Nᵣ
# Batched when the storage has a rotor axis — one more dimension than the modes (and the spin
# weights, if there are several) need — however many rotors it holds, so that this agrees with
# the blocks, whose type is decided by the same dimensions.
isbatched(Y::HarmonicValues) = ndims(Y.data) == (Y.s isa AbstractUnitRange ? 3 : 2)
spins(Y::HarmonicValues{T, IT, S}) where {T, IT, S<:IntegerHalf} = Y.s:Y.s
spins(Y::HarmonicValues{T, IT, S}) where {T, IT, S<:AbstractUnitRange} = Y.s
spin(Y::HarmonicValues{T, IT, S}) where {T, IT, S<:IntegerHalf} = Y.s

# `length` counts the blocks, as it does for a `WignerSeries`; the number of modes is
# `length(array_view(Y))` for the unbatched single-spin case, and `Ysize` in general.
Base.length(Y::HarmonicValues) = Int(ℓₘₐₓ(Y) - ℓₘᵢₙ(Y)) + 1
Base.keys(Y::HarmonicValues) = ℓₘᵢₙ(Y):ℓₘₐₓ(Y)
Base.firstindex(Y::HarmonicValues) = ℓₘᵢₙ(Y)
Base.lastindex(Y::HarmonicValues) = ℓₘₐₓ(Y)
# `first` and `last`, with or without a count, and `only` give blocks, as indexing does and as
# they do for a `WignerSeries`, rather than the pairs of the iteration.
Base.first(Y::HarmonicValues) = Y[ℓₘᵢₙ(Y)]
Base.last(Y::HarmonicValues) = Y[ℓₘₐₓ(Y)]
function Base.first(Y::HarmonicValues, n::Integer)
    n < 0 && throw(ArgumentError("Number of elements must be non-negative"))
    [Y[ℓₘᵢₙ(Y) + (i - 1)] for i ∈ 1:min(n, length(Y))]
end
function Base.last(Y::HarmonicValues, n::Integer)
    n < 0 && throw(ArgumentError("Number of elements must be non-negative"))
    k = min(n, length(Y))
    [Y[ℓₘₐₓ(Y) - (k - i)] for i ∈ 1:k]
end
function Base.only(Y::HarmonicValues)
    length(Y) == 1 || throw(ArgumentError(
        "These harmonic values hold $(length(Y)) blocks, for ℓ ∈ $(ℓₘᵢₙ(Y)):$(ℓₘₐₓ(Y)), "
        * "rather than exactly one."
    ))
    Y[ℓₘᵢₙ(Y)]
end

# The four shapes.  Which one applies is fixed by the rank of the storage and by whether `S` is
# a single spin weight or a range, so each of these has a single concrete return type.
@propagate_inbounds function Base.getindex(
    Y::HarmonicValues{T, IT, S, <:AbstractVector}, ℓ
) where {T, IT, S<:IntegerHalf}
    let ℓ = check_ℓ(Y, ℓ)
        DegreeBlock(view(Y.data, mode_range(Y, ℓ)), ℓ)
    end
end
@propagate_inbounds function Base.getindex(
    Y::HarmonicValues{T, IT, S, <:AbstractMatrix}, ℓ
) where {T, IT, S<:IntegerHalf}
    let ℓ = check_ℓ(Y, ℓ)
        DegreeBlockBatch(view(Y.data, :, mode_range(Y, ℓ)), ℓ)
    end
end
@propagate_inbounds function Base.getindex(
    Y::HarmonicValues{T, IT, S, <:AbstractMatrix}, ℓ
) where {T, IT, S<:AbstractUnitRange}
    let ℓ = check_ℓ(Y, ℓ), sr = Y.s
        SpinMatrix(
            view(Y.data, :, mode_range(Y, ℓ)), ℓ;
            sₘₐₓ=last(sr), sₘᵢₙ=first(sr), mₘₐₓ=ℓ, mₘᵢₙ=-ℓ
        )
    end
end
@propagate_inbounds function Base.getindex(
    Y::HarmonicValues{T, IT, S, <:AbstractArray{T, 3}}, ℓ
) where {T, IT, S<:AbstractUnitRange}
    let ℓ = check_ℓ(Y, ℓ), sr = Y.s
        SpinMatrixBatch(
            view(Y.data, :, :, mode_range(Y, ℓ)), ℓ;
            sₘₐₓ=last(sr), sₘᵢₙ=first(sr), mₘₐₓ=ℓ, mₘᵢₙ=-ℓ
        )
    end
end

# `Y[ℓ, :]` is the same block, so that the block of one ℓ is written the same way for both mode
# containers, as `w[ℓ, :]` is for a `ModeWeights`.
@propagate_inbounds Base.getindex(Y::HarmonicValues, ℓ, ::Colon) = Y[ℓ]

# Iteration yields `ℓ => block`, matching the calculators, so that a loop written against one
# reads the same against the other.
@inline function Base.iterate(Y::HarmonicValues{T, IT}, ℓ::IT=ℓₘᵢₙ(Y)) where {T, IT}
    ℓ > ℓₘₐₓ(Y) && return nothing
    (ℓ => Y[ℓ], ℓ + 1)
end
Base.IteratorSize(::Type{<:HarmonicValues}) = Base.HasLength()
Base.pairs(Y::HarmonicValues) = Y
# So the element type is that of the iteration, as it is for a calculator and a `WignerSeries`,
# rather than the number type that `AbstractModeContainer` reports for a `ModeWeights`, which
# iterates over its numbers; a disagreement makes `collect` throw.  The number type is
# `eltype(array_view(Y))`.
Base.eltype(::Type{H}) where {T, IT, H<:HarmonicValues{T, IT}} =
    Pair{IT, Base.promote_op(getindex, H, IT)}
Base.eltype(Y::HarmonicValues) = eltype(typeof(Y))
Base.IteratorEltype(::Type{<:HarmonicValues}) = Base.HasEltype()

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
function Base.show(io::IO, ::MIME"text/plain", Y::HarmonicValues)
    show(io, Y)
    println(io, ":")
    show_blocks(io, Y, ℓₘᵢₙ(Y), ℓₘₐₓ(Y))
end
