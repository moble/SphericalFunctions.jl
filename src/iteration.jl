### Iteration over ℓ
#
# The calculators hold one ``ℓ`` at a time, which is the whole point of them: the recurrence
# is what makes large ℓₘₐₓ reachable, and holding every matrix at once is what makes it
# expensive.  This file gives that shape its natural spelling, so that stepping through ℓ is a
# `for` loop rather than a hand-written `recurrence!`/`calc[ℓ]` pair whose two rules — call
# `recurrence!` first, and do not keep the block — nothing was enforcing.
#
# A calculator is indexed by ℓ as a *key*, not as a position, and iterates as key–value
# pairs.  That is the `AbstractDict` shape rather than the array shape, and the names follow
# it: `keys` is the range of ℓ, `length` counts ℓ values, `eltype` is a `Pair`, and `pairs`
# is the calculator itself.  Note that those last two names mean something else one level
# down: on a block, `length` counts matrix elements and `eltype` is the number type.
#
# The iteration state is the *next* ℓ, in the calculator's own index type, and never the
# calculator's internal position.  That is what lets a second loop over the same calculator
# restart cleanly from ℓₘᵢₙ — including after a `break` — instead of continuing from wherever
# the previous one stopped.

# The index type of a calculator, which is also its key type.
indextype(::Type{<:WignerCalculator{IT}}) where {IT} = IT
indextype(::Type{<:HarmonicCalculator{IT}}) where {IT} = IT
indextype(::Type{<:WignerHCalculator{IT}}) where {IT} = IT
indextype(c) = indextype(typeof(c))

# The calculators that yield one block per ℓ, and so need to be told nothing beyond the ℓ
# range.  An `sYlmCalculator` belongs here because it is built for the spin weights it serves:
# its block is the whole of what it holds, whether that is one spin weight or several.
const IterableCalculator{IT} = Union{WignerCalculator{IT}, HarmonicCalculator{IT}}

# The concrete type of the block that `calc[ℓ]` returns.  `promote_op` gets this exactly
# right for every combination of index type, batchedness and spin-weight count — the return
# type is settled at compile time, because `calc[ℓ]`'s shape is a type parameter — and doing
# it this way means `eltype` cannot drift away from what iteration actually yields.
blocktype(::Type{C}) where {C} = Base.promote_op(getindex, C, indextype(C))
blocktype(::Type{C}, ::Type{S}) where {C, S} = Base.promote_op(getindex, C, indextype(C), S)


"""
    eachℓ(calc)
    eachℓ(calc, s)
    eachell(calc, ...)

Iterator over the values of ``ℓ`` of a calculator, yielding `ℓ => block` pairs, where `block`
is what `calc[ℓ]` (or, in the second form, `calc[ℓ, s]`) returns.

The calculators iterate directly — `for (ℓ, 𝔇ˡ) ∈ calc` — so this function is needed in two
cases.  The keyword arguments `ℓₘᵢₙ` and `ℓₘₐₓ` restrict the range, which is otherwise the
calculator's own `ℓₘᵢₙ(calc):ℓₘₐₓ(calc)`; and an [`sYlmCalculator`](@ref) built for a range of
spin weights yields all of them at once, so a loop that wants one of them by itself has to say
which:

```julia
for (ℓ, 𝔇ˡ) ∈ eachℓ(calc; ℓₘᵢₙ=2)   # skip ℓ = 0, 1
    # ...
end
for (ℓ, ₛYₗ) ∈ eachℓ(calc, -2)       # one spin weight of a multi-spin sYlmCalculator
    # ...
end
```

Starting at an `ℓₘᵢₙ` above the calculator's own costs nothing in accuracy: the recurrence
runs through the intermediate values either way, and the result is bit-for-bit what a
sequential pass would give.

Each yielded block is a view into the calculator's storage and is overwritten by the next
step, so `copy` it if it must survive; `collect` on this iterator copies every block for you.
`eachell` is an ASCII alias.
"""
function eachℓ end

# `S` is `Nothing` for the calculators that yield one block per ℓ, and the spin weight's own
# type — necessarily `<:HalfInteger` — when one has been singled out of an `sYlmCalculator`.
# The methods below dispatch on `Nothing` against `<:HalfInteger` rather than on `S === IT`,
# because the latter pair would be ambiguous for the (unreachable) `IT === Nothing`, which
# Aqua's ambiguity check rightly flags.
struct EachEll{C, IT, S}
    calc::C
    ℓₘᵢₙ::IT
    ℓₘₐₓ::IT
    s::S
end

const eachell = eachℓ

function eachℓ(c::WignerCalculator{IT}; ℓₘᵢₙ=nothing, ℓₘₐₓ=nothing) where {IT}
    each_ℓ_helper(c, IT, ℓₘᵢₙ, ℓₘₐₓ, nothing)
end
function eachℓ(c::HarmonicCalculator{IT}; ℓₘᵢₙ=nothing, ℓₘₐₓ=nothing) where {IT}
    each_ℓ_helper(c, IT, ℓₘᵢₙ, ℓₘₐₓ, nothing)
end
function eachℓ(c::HarmonicCalculator{IT}, s; ℓₘᵢₙ=nothing, ℓₘₐₓ=nothing) where {IT}
    let s = convert(IT, s)
        check_spin(c, s)
        each_ℓ_helper(c, IT, ℓₘᵢₙ, ℓₘₐₓ, s)
    end
end
function eachℓ(w::WignerHCalculator, args...; kwargs...)
    error(
        "A WignerHCalculator is not iterable: its wedge is one mutable object handed back by "
        * "identity, rather than a view that `copy` can preserve.  Step it with "
        * "`recurrence!(calc, ℓ)` and read `calc.Hˡ`, or use a WignerDCalculator or "
        * "WignerdCalculator, which are iterable."
    )
end

# The limits arrive here positionally rather than as keywords, so that the accessors `ℓₘᵢₙ`
# and `ℓₘₐₓ` can be called by their own names for the defaults.
function each_ℓ_helper(c, ::Type{IT}, lo, hi, s) where {IT}
    lo = lo === nothing ? ℓₘᵢₙ(c) : convert(IT, lo)
    hi = hi === nothing ? ℓₘₐₓ(c) : convert(IT, hi)
    check_ℓ_range(c, lo, hi)
    EachEll(c, lo, hi, s)
end

function check_ℓ_range(c, lo, hi)
    if lo < ℓₘᵢₙ(c) || lo > ℓₘₐₓ(c) || hi < ℓₘᵢₙ(c) || hi > ℓₘₐₓ(c)
        error(
            "The requested range ℓ ∈ $lo:$hi is not within [$(ℓₘᵢₙ(c)), $(ℓₘₐₓ(c))] "
            * "for this calculator."
        )
    end
    nothing
end


### Iteration proper.  Bare iteration of a calculator routes through `eachℓ`, so that the two
### spellings cannot disagree.

@inline function Base.iterate(e::EachEll{C, IT, Nothing}, ℓ::IT=e.ℓₘᵢₙ) where {C, IT}
    ℓ > e.ℓₘₐₓ && return nothing
    recurrence!(e.calc, ℓ)
    (ℓ => e.calc[ℓ], ℓ + 1)
end
@inline function Base.iterate(e::EachEll{C, IT, S}, ℓ::IT=e.ℓₘᵢₙ) where {C, IT, S<:HalfInteger}
    ℓ > e.ℓₘₐₓ && return nothing
    recurrence!(e.calc, ℓ)
    (ℓ => e.calc[ℓ, e.s], ℓ + 1)
end

# `each_ℓ_all` is what `eachℓ` builds when neither limit is restricted; it is spelled out here
# because `eachℓ` takes keyword arguments, and a keyword call in the per-step path does not
# inline away — which cost 112 bytes per ℓ when this routed through `eachℓ` itself.
@inline each_ℓ_all(c::IterableCalculator) = EachEll(c, ℓₘᵢₙ(c), ℓₘₐₓ(c), nothing)
# These are `@inline`d for the same reason `each_ℓ_all` is spelled out above: the extra call
# layer is not always seen through, and when it is not, the `EachEll` and the block it yields
# are heap-allocated once per ℓ instead of being elided.  Iterating `eachℓ(calc)` was already
# allocation-free; without these, iterating the calculator itself was not.
@inline Base.iterate(c::WignerCalculator) = iterate(each_ℓ_all(c))
@inline Base.iterate(c::WignerCalculator{IT}, ℓ::IT) where {IT} = iterate(each_ℓ_all(c), ℓ)
@inline Base.iterate(c::HarmonicCalculator) = iterate(each_ℓ_all(c))
@inline Base.iterate(c::HarmonicCalculator{IT}, ℓ::IT) where {IT} = iterate(each_ℓ_all(c), ℓ)


### The container interface, for the calculators and for `EachEll` alike.

const EllIterable{IT} = Union{IterableCalculator{IT}, EachEll{<:Any, IT}}

# `keys` is a plain `UnitRange`, as it is for `WignerSeries`; the half-integer `WignerRange`
# is for array axes and cannot be printed or collected.  `length` subtracts the two limits
# rather than measuring `keys`, because `Int(ℓ)` throws for a half-integer ℓ while the
# difference of two of them is an `Int` by construction.
Base.keys(c::WignerCalculator) = ℓₘᵢₙ(c):ℓₘₐₓ(c)
Base.keys(c::HarmonicCalculator) = ℓₘᵢₙ(c):ℓₘₐₓ(c)
Base.keys(e::EachEll) = e.ℓₘᵢₙ:e.ℓₘₐₓ
Base.length(c::WignerCalculator) = Int(ℓₘₐₓ(c) - ℓₘᵢₙ(c)) + 1
Base.length(c::HarmonicCalculator) = Int(ℓₘₐₓ(c) - ℓₘᵢₙ(c)) + 1
Base.length(e::EachEll) = max(0, Int(e.ℓₘₐₓ - e.ℓₘᵢₙ) + 1)

Base.eltype(::Type{C}) where {C<:WignerCalculator} = Pair{indextype(C), blocktype(C)}
Base.eltype(::Type{C}) where {C<:HarmonicCalculator} = Pair{indextype(C), blocktype(C)}
Base.eltype(::Type{EachEll{C, IT, Nothing}}) where {C, IT} = Pair{IT, blocktype(C)}
Base.eltype(::Type{EachEll{C, IT, S}}) where {C, IT, S<:HalfInteger} = Pair{IT, blocktype(C, S)}
Base.eltype(c::EllIterable) = eltype(typeof(c))

# Both are already Base's defaults for a type it knows nothing else about; they are stated
# because they are part of the promise, and to pin them against a future `size` or `ndims`
# method silently changing them.
Base.IteratorSize(::Type{<:EllIterable}) = Base.HasLength()
Base.IteratorEltype(::Type{<:EllIterable}) = Base.HasEltype()

# Defining `keys` without this would give the generic `pairs`, which re-indexes `calc[ℓ]` for
# each key and so fails on the second step; this is what `AbstractDict` does.
Base.pairs(c::EllIterable) = c

"""
    collect(calc)
    collect(eachℓ(calc, ...))

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
Base.collect(c::EllIterable) = [ℓ => copy(block) for (ℓ, block) ∈ c]

function Base.show(io::IO, e::EachEll)
    print(io, "eachℓ(")
    show(io, e.calc)
    e.s === nothing || print(io, ", ", e.s)
    print(io, "; ℓₘᵢₙ=", e.ℓₘᵢₙ, ", ℓₘₐₓ=", e.ℓₘₐₓ, ")")
end
