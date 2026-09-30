### Plain-array access to the labelled containers.
#
# The containers in this package are deliberately not `AbstractArray`s: their natural indices
# may be half-odd-integers, which cannot satisfy that interface.  That keeps linear algebra
# from being applied to them by accident — which matters, because an `OffsetArray` with
# non-trivial offsets, the obvious alternative, accepts `*` and `mul!` and returns silently
# wrong answers.
#
# What is offered instead is an explicit, named route to the underlying numbers.
# `array_view` hands back a 1-based array aliasing the storage — a `StridedArray` for the
# containers the package builds, on which BLAS and LAPACK work at full speed — and `relabel`
# puts the natural indices back onto the result.  `Array` and `collect`, and `Matrix` for a
# two-axis block, remain the copying forms, for results that must outlive the storage.

"""
    array_view(w)

The contents of `w` as an ordinary 1-based array, **aliasing** its storage (returning a
view) rather than copying it.

The containers in this package are deliberately not `AbstractArray`s, because their natural
indices may be half-odd-integers; this is the explicit route from one of them to the plain
numbers, with the ``ℓ``, ``m′``, ``m`` and ``s`` labels dropped and every axis starting at
1.  For a [`ModeWeights`](@ref) or the result of [`sYlm`](@ref) it is the flat array in the
canonical mode ordering of [`Yindex`](@ref), which is what the transforms take:

```julia
Y = sYlm(R, ℓₘₐₓ, s)
array_view(Y)[Yindex(ℓ, m, ℓₘᵢₙ)]      # one mode of the flat array
```

An ordinary 1-based array is returned unchanged, so a function can accept either a labelled
container or a bare array without asking which it was given.  An array with other axes, such
as an `OffsetArray`, is refused with an `ArgumentError`, because the result is indexed from
`1`.

For the containers the package builds, over ordinary storage such as a `Matrix`, the result
is a `StridedArray`, which is the second reason to want it: BLAS needs a unit stride down
the first axis and a constant stride between columns — not contiguity — so a block can go
straight to `mul!`, `lu!`, `norm` and the rest with no copy at all:

```julia
𝔇₃ = array_view(𝔇₁[ℓ]) * array_view(𝔇₂[ℓ])
```

A single rotor's slice of a batched block has a leading stride of `Nᵣ`, so BLAS cannot take
it; `mul!` then falls back to the generic implementation, which is slower but correct.  The
answer is never wrong, only sometimes slow.

!!! warning
    The result aliases the storage of `w`, so writing through it writes into `w`.  When `w`
    came from a calculator, that storage is overwritten by the next call to
    [`recurrence!`](@ref) — use `Array(w)` or `collect(w)` for a copy that survives, or
    `Matrix(w)` for a two-axis block.

See also [`relabel`](@ref), which is the return leg.
"""
function array_view end

"""
    relabel(w, A)

A container with the natural indices of `w` and the data of the plain array `A`, which must
have the shape of `w`.

This is the return leg of a trip through linear algebra, so that a result comes back
labelled rather than as a bare array:

```julia
𝔇₃ = relabel(𝔇₁[ℓ], array_view(𝔇₁[ℓ]) * array_view(𝔇₂[ℓ]))
```

`A` is used as storage, not copied, exactly as the container constructors do.

The labels are those of `w`, whatever `A` holds.  For a [`ModeWeights`](@ref) that includes
the spin weight, so the result of a spin-changing operator matrix applied to
`array_view(w)`, such as `ð(s, ℓₘᵢₙ, ℓₘₐₓ) * array_view(w)`, is labelled correctly only by
`ModeWeights(A, s + 1, ℓₘᵢₙ, ℓₘₐₓ)`; the operator itself, `ð * w`, gives the labelled result
directly.

See also [`array_view`](@ref).
"""
function relabel end


### `array_view`

# `size(w)` is the extent of the *block*, which may be smaller than the storage it sits in,
# and a tuple whose length is fixed by the type of the block, so the index tuple is known to
# the compiler and the view is free.  The storage of a `DegreeBlock` is compared with the
# block first, since it may have been resized (see `check_storage`).
@inline function array_view(w::AbstractBlock)
    check_storage(w)
    view(parent(w), map(Base.OneTo, size(w))...)
end

# The mode containers store their data flat and 1-based already, so their storage *is* the
# flat 1-based form already; there is nothing to view.  For a `ModeWeights` this is what the
# transforms in `ssht/` reach for before every `mul!`, `ldiv!` and `reshape`; for a
# `HarmonicValues` it is the synthesis array that a product with mode weights takes.
@inline array_view(w::ModeWeights) = w.data
@inline array_view(Y::HarmonicValues) = Y.data

# An ordinary array is already in that form, and is returned unchanged.  This is what lets
# the transforms accept either a labelled container or a plain array without asking which
# they were given, and it is why an array with other axes is refused here: every caller
# indexes the result from 1, much of it under `@inbounds`.  The check is resolved at compile
# time for the array types whose axes are 1-based by their type.
@inline function array_view(x::AbstractArray)
    Base.require_one_based_indexing(x)
    x
end


### `relabel`

# The block constructors accept storage larger than the block, because a block may sit in
# storage sized for a larger one, as a calculator's blocks do.  The array given to `relabel`
# is to hold exactly the block, so its shape is compared with the block's own; otherwise a
# larger array would be accepted and labelled over its leading corner.  `size(w)` is the
# extent of the block, which is also the shape of `array_view(w)`.
function check_relabel_shape(w, A)
    if size(A) != size(w)
        throw(DimensionMismatch(
            "The array has size $(size(A)), but this block has size $(size(w))."
        ))
    end
end

# The array must have the rank of the block's storage, which is that of the block.
function relabel(
    w::AbstractBlock{IT, NT, <:AbstractArray{<:Any, N}}, A::AbstractArray{<:Any, N}
) where {IT, NT, N}
    check_relabel_shape(w, A)
    rewrap(w, A)
end
function relabel(w::ModeWeights, A::AbstractVector)
    if length(A) != length(w)
        throw(DimensionMismatch(
            "The vector has length $(length(A)), but these mode weights have length "
            * "$(length(w))."
        ))
    end
    ModeWeights(A, spin(w), ℓₘᵢₙ(w), ℓₘₐₓ(w))
end
# The constructor checks only the mode axis; the leading axes, which say how many rotors and
# spin weights there are, must match too, or the blocks would be of another shape.
function relabel(Y::HarmonicValues, A::AbstractArray)
    if size(A) != size(Y.data)
        throw(DimensionMismatch(
            "The array has size $(size(A)), but these harmonic values have size "
            * "$(size(Y.data))."
        ))
    end
    HarmonicValues(A, Y.s, Y.ℓₘᵢₙ, Y.ℓₘₐₓ, Y.Nᵣ)
end


### Broadcasting over a container.
#
# A block or a `HarmonicValues` takes part in a broadcast as the numbers of `array_view`,
# wrapped together with the container itself, so that the broadcast works on a plain 1-based
# array while the labels remain at hand.  Broadcasting pairs elements by position, which means
# the same thing for two containers only when they have the same labels, so a broadcast that
# combines containers, or writes one into another with `.=`, requires their labels to agree,
# as a broadcast over `ModeWeights` does.  A plain array combined with a container has no
# labels to compare, and is paired by position as it is.  The result is a plain 1-based
# array.  A `ModeWeights` has a style of its own, in `mode_weights.jl`, because the result of
# a broadcast over mode weights is labelled; a block combined with mode weights is a plain
# vector of factors there.

# The numbers of a container, as broadcasting sees them, with the container for its labels.
# `T` is the number type, `N` the number of dimensions, `A` the type of the array of
# numbers, and `C` that of the container.
struct LabelledArray{T, N, A<:AbstractArray{T, N}, C} <: AbstractArray{T, N}
    data::A
    container::C
end
Base.size(x::LabelledArray) = size(x.data)
Base.axes(x::LabelledArray) = axes(x.data)
Base.IndexStyle(::Type{<:LabelledArray{T, N, A}}) where {T, N, A} = IndexStyle(A)
@propagate_inbounds Base.getindex(x::LabelledArray, i::Int) = x.data[i]
@propagate_inbounds Base.getindex(x::LabelledArray{T, N}, I::Vararg{Int, N}) where {T, N} =
    x.data[I...]
# A broadcast that writes into the storage it reads from makes a copy of what it reads, and
# needs to know where that storage is.
Base.dataids(x::LabelledArray) = Base.dataids(x.data)
Base.unaliascopy(x::LabelledArray) = LabelledArray(Base.unaliascopy(x.data), x.container)

const LabelledContainer = Union{AbstractBlock, HarmonicValues}
Base.Broadcast.broadcastable(c::LabelledContainer) = LabelledArray(array_view(c), c)

# The labels that two containers combined in a broadcast must share.  Those of a block are its
# kind, its ℓ and its axes (see `block_labels`), and those of harmonic values say which spin
# weights, which ℓ and how many rotors the values are for.
container_labels(w::AbstractBlock) = block_labels(w)
container_labels(Y::HarmonicValues) = (:HarmonicValues, Y.s, Y.ℓₘᵢₙ, Y.ℓₘₐₓ, Y.Nᵣ)
container_description(w::AbstractBlock) = sprint(summary, w)
container_description(Y::HarmonicValues) = sprint(show, Y)

# The containers among the operands of a broadcast, however deeply nested, as a tuple.
labelled_operands(x) = ()
labelled_operands(x::LabelledArray) = (x.container,)
labelled_operands(bc::Broadcast.Broadcasted) = concatenated(map(labelled_operands, bc.args))
concatenated(::Tuple{}) = ()
concatenated(t::Tuple) = (first(t)..., concatenated(Base.tail(t))...)

# The first of the containers whose labels differ from `labels`, or `nothing`.
first_mislabelled(labels, ::Tuple{}) = nothing
first_mislabelled(labels, t::Tuple) =
    container_labels(first(t)) == labels ? first_mislabelled(labels, Base.tail(t)) : first(t)

function check_operand_labels(c, operands::Tuple)
    other = first_mislabelled(container_labels(c), operands)
    other === nothing || throw(broadcast_labels_error(c, other))
    nothing
end
@noinline function broadcast_labels_error(a, b)
    ArgumentError(
        "Broadcasting pairs the elements of labelled containers by position, so the "
        * "containers it combines must have the same labels; got $(container_description(a)) "
        * "and $(container_description(b)).  Use `array_view` on each to combine the numbers "
        * "by position regardless of their labels."
    )
end

# The broadcast style of a `LabelledArray`; `N` is the number of dimensions.
struct LabelledStyle{N} <: Broadcast.AbstractArrayStyle{N} end
LabelledStyle{N}(::Val{M}) where {N, M} = LabelledStyle{M}()
Base.BroadcastStyle(::Type{<:LabelledArray{T, N}}) where {T, N} = LabelledStyle{N}()
Base.BroadcastStyle(::ModeWeightsStyle, ::LabelledStyle{1}) = ModeWeightsStyle()
Base.BroadcastStyle(::ModeWeightsStyle, ::LabelledStyle{N}) where {N} = LabelledStyle{N}()
# The result is what it would be for plain arrays.
Base.similar(::Broadcast.Broadcasted{LabelledStyle{N}}, ::Type{T}, dims) where {N, T} =
    similar(Array{T}, dims)
Base.similar(::Broadcast.Broadcasted{LabelledStyle{N}}, ::Type{Bool}, dims) where {N} =
    similar(BitArray, dims)

# `instantiate` is the one step that every broadcast of this style passes through, with or
# without a destination, and is documented as the place for a style to check its operands.
function Base.Broadcast.instantiate(bc::Broadcast.Broadcasted{<:LabelledStyle})
    operands = labelled_operands(bc)
    isempty(operands) || check_operand_labels(first(operands), Base.tail(operands))
    invoke(Base.Broadcast.instantiate, Tuple{Broadcast.Broadcasted}, bc)
end
# No container is zero-dimensional; this settles the ambiguity with `Base`'s method for the
# zero-dimensional styles.
Base.Broadcast.instantiate(bc::Broadcast.Broadcasted{LabelledStyle{0}}) = bc


### Broadcast assignment into a container.
#
# Writing *into* a container with `.=` needs these: `v .= x` lowers to `materialize!(v,
# broadcasted(identity, x))`, and without a method here that reaches `copyto!` on a type that
# has no `copyto!`.  The containers on the right-hand side must have the labels of the
# destination, for the reason given above.  A `ModeWeights` destination has methods of its
# own, in `mode_weights.jl`, which check its labels in its own terms.

# Two signatures, because `Base` defines both `materialize!(dest, bc)` and the more specific
# `materialize!(dest, bc::Broadcasted{<:Any})`; without the second of these, a `.=` whose
# right-hand side is already a `Broadcasted` is ambiguous rather than dispatched.
@inline function Base.Broadcast.materialize!(dest::LabelledContainer, bc)
    Base.Broadcast.materialize!(array_view(dest), bc)
    dest
end
@inline function Base.Broadcast.materialize!(
    dest::LabelledContainer, bc::Base.Broadcast.Broadcasted{<:Any}
)
    check_operand_labels(dest, labelled_operands(bc))
    Base.Broadcast.materialize!(array_view(dest), bc)
    dest
end
