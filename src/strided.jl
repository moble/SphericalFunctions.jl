### Plain-array access to the labelled containers.
#
# The containers in this package are deliberately not `AbstractArray`s: their natural indices
# may be half-odd-integers, which cannot satisfy that interface.  That keeps linear algebra
# from being applied to them by accident — which matters, because the `OffsetArray`s used
# before version 3 accepted `*` and `mul!` and returned silently wrong answers.
#
# What is offered instead is an explicit, named route to the underlying numbers.  `strided`
# hands back a 1-based `StridedArray` aliasing the storage, on which BLAS and LAPACK work at
# full speed; `relabel` puts the natural indices back onto the result.  `Matrix`, `Array` and
# `collect` remain the copying forms, for results that must outlive the storage.

"""
    strided(w)

The block `w` as an ordinary 1-based `StridedArray`, **aliasing** the storage of `w` rather
than copying it.

This is the short route from a labelled container to linear algebra.  BLAS needs a unit
stride down the first axis and a constant stride between columns — not contiguity — so the
result of an unbatched calculator can go straight to `mul!`, `lu!`, `norm` and the rest with
no copy at all:

```julia
𝔇₃ = strided(𝔇₁[ℓ]) * strided(𝔇₂[ℓ])
```

A single rotor's slice of a batched block has a leading stride of `Nᵣ`, so BLAS cannot take
it; `mul!` then falls back to the generic implementation, which is slower but correct.  The
answer is never wrong, only sometimes slow.

!!! warning
    The result aliases the storage of `w`, so writing through it writes into `w`.  When `w`
    came from a calculator, that storage is overwritten by the next call to
    [`recurrence!`](@ref) — use `Matrix(w)`, `Array(w)` or `collect(w)` for a copy that
    survives.

See also [`relabel`](@ref), which is the return leg.
"""
function strided end

"""
    relabel(w, A)

A container with the natural indices of `w` and the data of the plain array `A`, which must
have the shape of `w`.

This is the return leg of a trip through linear algebra, so that a result comes back labelled
rather than as a bare array:

```julia
𝔇₃ = relabel(𝔇₁[ℓ], strided(𝔇₁[ℓ]) * strided(𝔇₂[ℓ]))
```

`A` is used as storage, not copied, exactly as the container constructors do.

See also [`strided`](@ref).
"""
function relabel end


### `strided`

# The index tuple is built with a `Val` so that its length is known to the compiler and the
# view is free; `size(w)` is the extent of the *block*, which may be smaller than the storage
# it sits in.
@inline block_axes(w, ::Val{N}) where {N} = ntuple(d -> Base.OneTo(size(w, d)), Val(N))

@inline strided(w::WignerMatrix) = view(parent(w), block_axes(w, Val(2))...)
@inline strided(w::WignerMatrixBatch) = view(parent(w), block_axes(w, Val(3))...)
@inline strided(w::DegreeBlock) = view(parent(w), block_axes(w, Val(1))...)
@inline strided(w::DegreeBlockBatch) = view(parent(w), block_axes(w, Val(2))...)
@inline strided(w::SpinMatrix) = view(parent(w), block_axes(w, Val(2))...)
@inline strided(w::SpinMatrixBatch) = view(parent(w), block_axes(w, Val(3))...)

# The mode containers store their data flat and 1-based already, so their storage *is* the
# strided form; there is nothing to view.  For a `ModeWeights` this is what the transforms in
# `ssht/` reach for before every `mul!`, `ldiv!` and `reshape`; for a `HarmonicValues` it is the
# synthesis array that a product with mode weights takes.
@inline strided(w::ModeWeights) = w.data
@inline strided(Y::HarmonicValues) = Y.data

# An ordinary array is its own strided form.  This is what lets the transforms accept either
# a labelled container or a plain array without asking which they were given.
@inline strided(x::AbstractArray) = x


### `relabel`

function relabel(w::WignerMatrix, A::AbstractMatrix)
    WignerMatrix(A, ℓ(w); m′ₘₐₓ=m′ₘₐₓ(w), m′ₘᵢₙ=m′ₘᵢₙ(w), mₘₐₓ=mₘₐₓ(w), mₘᵢₙ=mₘᵢₙ(w))
end
function relabel(w::WignerMatrixBatch, A::AbstractArray{<:Any, 3})
    WignerMatrixBatch(A, ℓ(w); m′ₘₐₓ=m′ₘₐₓ(w), m′ₘᵢₙ=m′ₘᵢₙ(w), mₘₐₓ=mₘₐₓ(w), mₘᵢₙ=mₘᵢₙ(w))
end
relabel(w::DegreeBlock, A::AbstractVector) =
    DegreeBlock(A, ℓ(w); mₘₐₓ=mₘₐₓ(w), mₘᵢₙ=mₘᵢₙ(w))
relabel(w::DegreeBlockBatch, A::AbstractMatrix) =
    DegreeBlockBatch(A, ℓ(w); mₘₐₓ=mₘₐₓ(w), mₘᵢₙ=mₘᵢₙ(w))
function relabel(w::SpinMatrix, A::AbstractMatrix)
    SpinMatrix(A, ℓ(w); sₘₐₓ=sₘₐₓ(w), sₘᵢₙ=sₘᵢₙ(w), mₘₐₓ=mₘₐₓ(w), mₘᵢₙ=mₘᵢₙ(w))
end
function relabel(w::SpinMatrixBatch, A::AbstractArray{<:Any, 3})
    SpinMatrixBatch(A, ℓ(w); sₘₐₓ=sₘₐₓ(w), sₘᵢₙ=sₘᵢₙ(w), mₘₐₓ=mₘₐₓ(w), mₘᵢₙ=mₘᵢₙ(w))
end
relabel(w::ModeWeights, A::AbstractVector) = ModeWeights(A, spin(w), ℓₘᵢₙ(w), ℓₘₐₓ(w))


### Broadcast assignment into a container.
#
# Reading *out* of a container by broadcasting already works: `Base.Broadcast.broadcastable`
# falls back to `collect`, which is why `v .+ 1` gives an ordinary 1-based `Array`.  Writing
# *into* one is what needs this: `v .= x` lowers to `materialize!(v, broadcasted(identity, x))`,
# and without a method here that reaches `copyto!` on a type that has no `copyto!`.
const NaturalContainer = Union{
    WignerMatrix, WignerMatrixBatch, DegreeBlock, DegreeBlockBatch, SpinMatrix, SpinMatrixBatch,
    ModeWeights,
}

# Two signatures, because `Base` defines both `materialize!(dest, bc)` and the more specific
# `materialize!(dest, bc::Broadcasted{<:Any})`; without the second of these, a `.=` whose
# right-hand side is already a `Broadcasted` is ambiguous rather than dispatched.
@inline function Base.Broadcast.materialize!(dest::NaturalContainer, bc)
    Base.Broadcast.materialize!(strided(dest), bc)
    dest
end
@inline function Base.Broadcast.materialize!(
    dest::NaturalContainer, bc::Base.Broadcast.Broadcasted{<:Any}
)
    Base.Broadcast.materialize!(strided(dest), bc)
    dest
end
