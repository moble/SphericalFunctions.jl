### Iteration over ℓ
#
# The calculators hold one ``ℓ`` at a time, which is the whole point of them: the recurrence
# is what makes large ℓₘₐₓ reachable, and holding every matrix at once is what makes it
# expensive.  This file gives that shape its natural expression, so that a whole sweep is a
# `for` loop rather than a hand-written one over `recurrence!`.
#
# A calculator iterates as key–value pairs, ℓ => block.  That is the `AbstractDict` shape
# rather than the array shape, and the names follow it: `keys` is the range of ℓ, `length`
# counts ℓ values, `eltype` is a `Pair`, and `pairs` is the calculator itself.  Note that
# those last two names mean something else one level down: on a block, `length` counts matrix
# elements and `eltype` is the number type.
#
# A calculator is *not* indexed.  `recurrence!(calc, ℓ)` computes one ℓ and returns its block,
# which is the whole of the random-access story and is also what each step of the loop below
# calls; a sweep over part of the range is a `for` loop over that part.
#
# The iteration state is the *next* ℓ, in the calculator's own index type, and never the
# calculator's internal position.  That is what lets a second loop over the same calculator
# restart cleanly from ℓₘᵢₙ — including after a `break` — instead of continuing from wherever
# the previous one stopped.

# The index type of a calculator, which is also its key type.
indextype(::Type{<:WignerCalculator{IT}}) where {IT} = IT
indextype(::Type{<:HarmonicCalculator{IT}}) where {IT} = IT
indextype(::Type{<:HCalculator{IT}}) where {IT} = IT
indextype(c) = indextype(typeof(c))

# The calculators that yield one block per ℓ.  An `sYlmCalculator` belongs here because it is
# built for the spin weights it serves: its block is the whole of what it holds, whether that
# is one spin weight or several.
const IterableCalculator{IT} = Union{WignerCalculator{IT}, HarmonicCalculator{IT}}

# The concrete type of the block that `recurrence!` returns.  `promote_op` gets this exactly
# right for every combination of index type, batchedness and spin-weight count — the return
# type is settled at compile time, because the block's shape is a type parameter — and doing
# it this way means `eltype` cannot drift away from what iteration actually yields.
blocktype(::Type{C}) where {C} = Base.promote_op(recurrence!, C, indextype(C))


### Iteration proper.  Each step is one `recurrence!`, so the loop and the hand-written form
### cannot disagree about what a block is.
#
# These are `@inline`d because the extra call layer is not always seen through, and when it is
# not, the block is heap-allocated once per ℓ instead of being elided.

@inline function Base.iterate(c::WignerCalculator{IT}, ℓ::IT=ℓₘᵢₙ(c)) where {IT}
    ℓ > ℓₘₐₓ(c) && return nothing
    (ℓ => recurrence!(c, ℓ), ℓ + 1)
end
@inline function Base.iterate(c::HarmonicCalculator{IT}, ℓ::IT=ℓₘᵢₙ(c)) where {IT}
    ℓ > ℓₘₐₓ(c) && return nothing
    (ℓ => recurrence!(c, ℓ), ℓ + 1)
end

# The one calculator that is not iterable, refused by name rather than by `MethodError`,
# because the reason is worth stating and the remedy is not obvious.  `IteratorSize` is
# declared too, so that `collect` and a comprehension reach `iterate` — and that message —
# rather than failing first on the missing `length`.
Base.IteratorSize(::Type{<:HCalculator}) = Base.SizeUnknown()
function Base.iterate(w::HCalculator, args...)
    error(
        "An HCalculator is not iterable: its wedge is one mutable object handed back by "
        * "identity, rather than a view that `copy` can preserve.  Step it with "
        * "`recurrence!(calc, ℓ)`, which returns that wedge, or use a DCalculator or "
        * "dCalculator, which are iterable."
    )
end


### The container interface.

# `keys` is a plain `UnitRange`, as it is for `WignerSeries`; the half-integer `WignerRange`
# is for array axes and cannot be printed or collected.  `length` subtracts the two limits
# rather than measuring `keys`, because `Int(ℓ)` throws for a half-integer ℓ while the
# difference of two of them is an `Int` by construction.
Base.keys(c::WignerCalculator) = ℓₘᵢₙ(c):ℓₘₐₓ(c)
Base.keys(c::HarmonicCalculator) = ℓₘᵢₙ(c):ℓₘₐₓ(c)
Base.length(c::WignerCalculator) = Int(ℓₘₐₓ(c) - ℓₘᵢₙ(c)) + 1
Base.length(c::HarmonicCalculator) = Int(ℓₘₐₓ(c) - ℓₘᵢₙ(c)) + 1

Base.eltype(::Type{C}) where {C<:WignerCalculator} = Pair{indextype(C), blocktype(C)}
Base.eltype(::Type{C}) where {C<:HarmonicCalculator} = Pair{indextype(C), blocktype(C)}
Base.eltype(c::IterableCalculator) = eltype(typeof(c))

# Both are already Base's defaults for a type it knows nothing else about; they are stated
# because they are part of the promise, and to pin them against a future `size` or `ndims`
# method silently changing them.
Base.IteratorSize(::Type{<:IterableCalculator}) = Base.HasLength()
Base.IteratorEltype(::Type{<:IterableCalculator}) = Base.HasEltype()

# Iteration already yields pairs, so the calculator is its own `pairs`; this says so rather
# than leaving the generic wrapper to pair the pairs with positions.
Base.pairs(c::IterableCalculator) = c

"""
    collect(calc)

Every ``ℓ`` of a calculator, as an ordinary 1-based `Vector` of `ℓ => block` pairs, with each
block copied.

The copies are the point: the blocks a calculator yields are views into storage that the next
step overwrites, so a `collect` that did not copy would hand back `length(calc)` aliases of
one buffer, all holding the last ``ℓ``.  Each copy keeps its natural ``(m′, m)`` indices, and
generic functions that store what they are handed — `map`, `first(calc, n)`,
`Iterators.take` — do *not* copy, so use this or `copy` the block inside the loop.

Note that the result is indexed from 1, while [`D`](@ref) and [`d`](@ref) return containers
indexed by ``ℓ``; `collect(calc)[i]` is the pair whose first element is the ``ℓ``.
"""
Base.collect(c::IterableCalculator) = [ℓ => copy(block) for (ℓ, block) ∈ c]
