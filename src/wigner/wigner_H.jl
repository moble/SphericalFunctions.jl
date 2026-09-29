### Index arithmetic
#
# Indices are `Integer`s or `HalfOddInteger`s — together, `IntegerHalf`.  Because a sum or
# difference of two `HalfOddInteger`s is an `Int`, and adding an `Int` to one gives back a
# `HalfOddInteger`, every coefficient, loop bound and storage offset below is integer
# arithmetic for *both* index types, while being written exactly as the references write it.
# For an integer index type the expressions are literally the ones a purely integer
# implementation would use, so that path is unchanged down to the last bit.

"""
    δ²(ℓ, m)

``(ℓ - m)(ℓ + m + 1) = (2 d^m_ℓ)^2``, four times the square of Gumerov and Duraiswami's
coefficient ``d^m_ℓ = \\mathrm{sgn}(m) \\sqrt{(ℓ-m)(ℓ+m+1)} / 2``; the factor cancels in
each recurrence, and the sign is applied separately (see [`sgn`](@ref)).  Both factors are
`Integer`s for integer and half-integer indices alike, so this is exact, and the value
handed to `sqrt` is the same one a purely integer implementation computes.
"""
@inline δ²(ℓ, m) = (ℓ - m) * (ℓ + m + 1)

"""
    sgn(m)

Eq. (44) in Gumerov and Duraiswami (2015).  Note that they define `sgn` differently from the
usual definition — including from Julia's `sign` — at 0, where this is ``+1``.
"""
@inline sgn(m) = ifelse(m ≥ 0, 1, -1)

"""
    ϵ(m)

Eq. (7) in Gumerov and Duraiswami (2015): ``ε(m) = (-1)^m`` for ``m > 0``, and ``1``
otherwise.  The half-integer extension ``ε(m) = (-1)^{⌊m⌋}`` for ``m > 0`` is the unique one
that keeps the ``H``-form of the ``m′`` ladder in Gumerov and Duraiswami's shape, as step 4
of the notes on the [``H`` recursion](@ref "Algorithm for computing ``H``") explains.  The
single expression below serves both index types.
"""
@inline ϵ(m) = ifelse(m > 0 && isodd(floor_int(m)), -1, 1)

"""
    HWedge{IT, RT, ST} <: AbstractWignerMatrix{IT, RT, ST}
    HWedge([RT=Float64,] Nᵣ, ℓₘₐₓ, m′ₘₐₓ=ℓₘₐₓ)

A compact, real-valued workspace holding the ``Hˡ`` matrix of one ``ℓ`` for `Nᵣ` rotors at
once.
- `IT` is the index type, `Int` or [`HalfOddInteger`](@ref).
- `RT` is the real type of the elements.
- `ST` is the type of the flat storage, a vector of `RT`.

The ``Hˡ`` matrix is critical to efficient and stable computation of the Wigner ``D`` and
``d`` matrices — in fact, it essentially *is* the ``d`` matrix with signs adjusted to avoid
numerical problems with alternating signs.  This gives it additional symmetries that reduce
the amount of data that needs to be stored to about 1/4 of the total ``d`` size.

The purpose of an `HWedge` is to provide a workspace for the Wigner recurrences that is
efficient, both in terms of the size of memory used, and the implications for vectorization.
Specifically, the data is stored as strictly `Real` values, in contiguous storage.  Indexing
is performed efficiently via precomputed row offsets.  Once the full recurrence is done, the
data can be used directly — computing phases and symmetry on the fly — or copied into a full
explicit matrix with the appropriate phases.

The recurrences require ``m`` in the full range from 0 (or 1/2) to ``ℓ``, but ``m'`` only
needs the rows ``|m'| ≤ m'_{\\mathrm{max}}``, which must include the axis ``m'=0`` for
integer indices and both rows ``m' = ±1/2`` for half-integer ones.  The range of ``m'`` is
always symmetric, so `m′ₘᵢₙ(H)` is `-m′ₘₐₓ(H)`, and only enough storage for those rows is
allocated.  Specifically, an `HWedge` stores elements in a vector as if they were components
of the `Hˡ` matrix:

    [
        Hˡ[m′, m]
        for m′ ∈ max(-ℓ, -m′ₘₐₓ):min(ℓ, m′ₘₐₓ)
        for m ∈ abs(m′):ℓ
    ]

However, for further efficiency when vectorizing over multiple rotors, the data is stored as
a 1-dimensional vector, though it can be indexed as if it were a three-dimensional array,
with the first dimension indexing `Nᵣ` different rotors, and the second dimension indexing
`m′`, and the third dimension indexing `m`.  Thus, this object can be indexed as `Hˡ[iᵣ, m′,
m]` to get the `Hˡ` value for rotor index `iᵣ`, and matrix element `(m′, m)`.  Accordingly
`ndims(Hˡ)` is 3 and `axes(Hˡ)` and `size(Hˡ, d)` describe those three indices, while
`length(Hˡ)` and `size(Hˡ)` describe the flat storage, which is what linear indexing `Hˡ[i]`
runs over.  Only the elements with ``m ≥ |m'|`` are stored; the others are read through
[`wedge_value`](@ref), which applies the symmetries of ``H``.

Because of this complicated layout, the constructor is fairly restrictive, but will do all
the allocation needed.  `ℓₘₐₓ` and `m′ₘₐₓ` may be integers, or half-integers given as
`Rational`s with denominator 2 or as [`HalfOddInteger`](@ref)s, as for the other [index
arguments](@ref "Functions that take indices").  To avoid multiple allocations, it is
advisable to first construct an instance with the maximum `ℓ` value that will be needed, and
then change the `ℓ` field as needed to compute different orders.  That is, if `H isa
HWedge`, then `H.ℓ = new_ell` can be used to change the current order being computed; this
lays out the storage for the new order without computing anything.  The constructor starts
out with the smallest `ℓ` value possible (0 or 1/2), which is the natural choice for
recurrence.  The wedge an [`HCalculator`](@ref) returns is that calculator's own workspace,
so its `ℓ` should not be reassigned; `copy(Hˡ)` gives an independent wedge holding the same
numbers, which survives the calculator's next step.

!!! warning "A single-owner workspace"

    An `HWedge` is a mutable workspace, and none of its operations may run on two tasks at
    once: in particular, changing the `ℓ` field rewrites the row offsets that every read
    uses.  Give each task its own wedge, or its own calculator.

"""
mutable struct HWedge{IT, RT<:Real, ST} <: AbstractWignerMatrix{IT, RT, ST}
    const parent::ST
    const row_index::FixedSizeVectorDefault{Int}
    const Nᵣ::Int
    const maxℓ::IT
    const maxm′ₘₐₓ::IT
    ℓ::IT
    m′ₘₐₓ::IT
    @index_methods function HWedge(
        ::Type{RT}, Nᵣ::Int, ℓₘₐₓ::IT, m′ₘₐₓ::IT=ℓₘₐₓ
    ) where {IT<:IndexType, RT<:Real}
        if Nᵣ < 1
            throw(ArgumentError("Number of rotors Nᵣ=$Nᵣ must be at least 1."))
        end
        validate_index_ranges(ℓₘₐₓ, m′ₘₐₓ, -m′ₘₐₓ)

        # Set up storage for the biggest these values will ever be
        parent = FixedSizeVector{RT}(undef, Nᵣ * HWedge_size(ℓₘₐₓ, m′ₘₐₓ, -m′ₘₐₓ))
        row_index = FixedSizeVector{Int}(undef, Int(m′ₘₐₓ - (-m′ₘₐₓ)) + 1)

        # But start out assuming ℓ is the smallest it can be
        maxℓ = ℓₘₐₓ
        maxm′ₘₐₓ = m′ₘₐₓ
        ℓ = ℓₘᵢₙ(ℓₘₐₓ)
        m′ₘₐₓ = min(ℓ, m′ₘₐₓ)
        HWedge_row_index!(row_index, Nᵣ, ℓ, m′ₘₐₓ, -m′ₘₐₓ)

        new{IT, RT, typeof(parent)}(parent, row_index, Nᵣ, maxℓ, maxm′ₘₐₓ, ℓ, m′ₘₐₓ)
    end
    @index_methods HWedge(Nᵣ::Int, ℓₘₐₓ::IT, m′ₘₐₓ::IT=ℓₘₐₓ) where {IT<:IndexType} =
        HWedge(Float64, Nᵣ, ℓₘₐₓ, m′ₘₐₓ)
end

# The range of m′ is symmetric, so its lower limits are not stored.
m′ₘᵢₙ(w::HWedge{IT}) where {IT} = -w.m′ₘₐₓ
mₘₐₓ(w::HWedge{IT}) where {IT} = ℓ(w)
mₘᵢₙ(w::HWedge{IT}) where {IT} = ℓₘᵢₙ(w)

row_index(w::HWedge{IT}) where {IT} = w.row_index
row_index(w::HWedge{IT}, m′::IT) where {IT} = row_index(w)[Int(m′ - m′ₘᵢₙ(w)) + 1]
Nᵣ(w::HWedge{IT}) where {IT} = w.Nᵣ
maxℓ(w::HWedge{IT}) where {IT} = w.maxℓ
maxm′ₘₐₓ(w::HWedge{IT}) where {IT} = w.maxm′ₘₐₓ
minm′ₘᵢₙ(w::HWedge{IT}) where {IT} = -w.maxm′ₘₐₓ

function Base.setproperty!(H::HWedge{IT}, s::Symbol, ℓ::IIT) where {IT, IIT}
    if s === :ℓ
        if IIT !== IT
            throw(ArgumentError(
                "Cannot change ℓ from type $IT to type $IIT; they must be the same."
            ))
        end
        if ℓ < ℓₘᵢₙ(IT)
            throw(ArgumentError("Cannot set ℓ=$ℓ less than ℓₘᵢₙ=$(ℓₘᵢₙ(IT))."))
        end
        if ℓ > maxℓ(H)
            throw(ArgumentError("Cannot set ℓ=$ℓ greater than maxℓ=$(maxℓ(H))."))
        end
        m′ₘₐₓ = min(ℓ, maxm′ₘₐₓ(H))
        HWedge_row_index!(row_index(H), Nᵣ(H), ℓ, m′ₘₐₓ, -m′ₘₐₓ)
        Base.setfield!(H, :ℓ, ℓ)
        Base.setfield!(H, :m′ₘₐₓ, m′ₘₐₓ)
        ℓ
    else
        throw(ArgumentError(
            "Cannot set property `$s` on HWedge; only `ℓ` is allowed to be changed."
        ))
    end
end
# A half-integer order may be written as a `Rational`, as it may at every other entry point.
function Base.setproperty!(H::HWedge{HalfOddInteger}, s::Symbol, ℓ::Rational)
    s === :ℓ || throw(ArgumentError(
        "Cannot set property `$s` on HWedge; only `ℓ` is allowed to be changed."
    ))
    setproperty!(H, :ℓ, calculator_index(HalfOddInteger, ℓ, "ℓ", "wedge"))
end

# Recompute the row offsets; called once per `H.ℓ = ℓ` assignment, hence once per
# `recurrence!`.  Row `m′` holds the `ℓ - |m′| + 1` elements `m ∈ |m′|:ℓ`, each `Nᵣ` wide.
function HWedge_row_index!(row_index, Nᵣ::Int, ℓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT) where {IT}
    index = 1
    i = 1
    for m′ ∈ m′ₘᵢₙ:m′ₘₐₓ
        @inbounds row_index[i] = index
        index += Nᵣ * ((ℓ - abs(m′)) + 1)
        i += 1
    end
    row_index
end

function HWedge_size(ℓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT) where {IT}
    let ℓₘᵢₙ = ℓₘᵢₙ(IT)
        (
            (ℓₘᵢₙ - m′ₘᵢₙ) * ((2ℓ + 1) + (m′ₘᵢₙ + ℓₘᵢₙ))
            - (ℓₘᵢₙ - m′ₘₐₓ - 1) * ((2ℓ + 2) - (m′ₘₐₓ + ℓₘᵢₙ))
        ) ÷ 2
    end
end

# A wedge is indexed by three indices, `[iᵣ, m′, m]`, which is what `axes`, `ndims` and
# `size(w, d)` describe; `length` and `size(w)` describe the flat storage that linear
# indexing runs over, since not every `(m′, m)` of the axes is stored.
function Base.axes(w::HWedge{IT}) where {IT}
    (1:Nᵣ(w), WignerRange(m′ₘᵢₙ(w):m′ₘₐₓ(w)), WignerRange(mₘᵢₙ(w):mₘₐₓ(w)))
end
Base.ndims(::HWedge) = 3
Base.ndims(::Type{<:HWedge}) = 3
Base.size(w::HWedge, d::Integer) = length(axes(w, d))

function Base.checkbounds(::Type{Bool}, w::HWedge, i::Int)
    i ≥ 1 && i ≤ length(w)
end
function Base.checkbounds(::Type{Bool}, w::HWedge{IT}, iᵣ::Int, m′::IT, m::IT) where {IT}
    iᵣ > 0 && iᵣ ≤ Nᵣ(w) && abs(m′) ≤ m ≤ ℓ(w) && m′ₘᵢₙ(w) ≤ m′ ≤ m′ₘₐₓ(w)
end

@propagate_inbounds function Base.getindex(w::HWedge, i::Int)
    @boundscheck if !checkbounds(Bool, w, i)
        throw(BoundsError(w, i))
    end
    @inbounds Base.parent(w)[i]
end
# As a block container of `wigner_matrix.jl` may be, a half-integer wedge may be indexed by
# `Rational`s, which are converted to `HalfOddInteger`s.
@propagate_inbounds Base.getindex(w::HWedge{IT}, iᵣ::Int, m′::Rational, m::Rational) where
    {IT<:HalfOddInteger} = w[iᵣ, HalfOddInteger(m′), HalfOddInteger(m)]
@propagate_inbounds Base.setindex!(w::HWedge{IT}, v, iᵣ::Int, m′::Rational, m::Rational) where
    {IT<:HalfOddInteger} = (w[iᵣ, HalfOddInteger(m′), HalfOddInteger(m)] = v)

@propagate_inbounds function Base.getindex(w::HWedge{IT}, iᵣ::Int, m′::IT, m::IT) where {IT}
    @boundscheck if !checkbounds(Bool, w, iᵣ, m′, m)
        throw(BoundsError(w, (iᵣ, m′, m)))
    end
    i = @inbounds (iᵣ - 1) + Nᵣ(w) * Int(m - abs(m′)) + row_index(w)[Int(m′ - m′ₘᵢₙ(w)) + 1]
    @inbounds Base.parent(w)[i]
end

@propagate_inbounds function Base.setindex!(w::HWedge, v, i::Int)
    @boundscheck if !checkbounds(Bool, w, i)
        throw(BoundsError(w, i))
    end
    @inbounds Base.parent(w)[i] = v
end
@propagate_inbounds function Base.setindex!(w::HWedge{IT}, v, iᵣ::Int, m′::IT, m::IT) where {IT}
    @boundscheck if !checkbounds(Bool, w, iᵣ, m′, m)
        throw(BoundsError(w, (iᵣ, m′, m)))
    end
    i = @inbounds (iᵣ - 1) + Nᵣ(w) * Int(m - abs(m′)) + row_index(w)[Int(m′ - m′ₘᵢₙ(w)) + 1]
    @inbounds Base.parent(w)[i] = v
end

# A wedge is not a matrix, so the generic `==`, which compares `Matrix` views, does not
# apply: two wedges are equal when they describe the same ℓ, range of m′ and rotors, and
# agree on every element they hold.
function Base.:(==)(w1::HWedge{IT}, w2::HWedge{IT}) where {IT}
    ℓ(w1) == ℓ(w2) && Nᵣ(w1) == Nᵣ(w2) &&
        m′ₘᵢₙ(w1) == m′ₘᵢₙ(w2) && m′ₘₐₓ(w1) == m′ₘₐₓ(w2) &&
        all(
            w1[iᵣ, m′, m] == w2[iᵣ, m′, m]
            for m′ ∈ m′ₘᵢₙ(w1):m′ₘₐₓ(w1) for m ∈ abs(m′):ℓ(w1) for iᵣ ∈ 1:Nᵣ(w1)
        )
end

# A copy is a wedge of its own, holding the same numbers and laid out for the same ℓ.  It is
# a deep copy because the storage may hold `#undef` entries (of `BigFloat`, say) beyond
# those the current ℓ uses, which `copyto!` would refuse to read.
Base.copy(w::HWedge) = deepcopy(w)

function Base.summary(io::IO, H::HWedge{IT, RT}) where {IT, RT}
    print(
        io,
        "HWedge{$IT, $RT} for ℓ=$(ℓ(H)) with m′=$(m′ₘᵢₙ(H)):$(m′ₘₐₓ(H)), ",
        "m=abs(m′):$(ℓ(H)), and iᵣ=1:$(Nᵣ(H))"
    )
end
Base.show(io::IO, H::HWedge) = summary(io, H)
function Base.show(io::IO, ::MIME"text/plain", H::HWedge{IT, RT, ST}) where {IT, RT, ST}
    summary(io, H)
    print(io, " stored in\n", summary(parent(H)), ", currently using\n")
    let ℓ = ℓ(H), m′ₘᵢₙ = m′ₘᵢₙ(H), m′ₘₐₓ = m′ₘₐₓ(H), Nᵣ = Nᵣ(H)
        i = row_index(H)[Int(m′ₘₐₓ - m′ₘᵢₙ) + 1] + Nᵣ * (Int(ℓ - abs(m′ₘₐₓ)) + 1) - 1
        # A view, not a copy, so that uninitialized `BigFloat` storage prints as `#undef`
        show(io, MIME("text/plain"), view(parent(H), firstindex(parent(H)):i))
    end
end

"""
    transpose_sign(m′, m)

Sign ``σ`` relating the transposed element of the ``H`` matrix to the original: ``H_{m,m′} =
σ H_{m′,m}`` and ``H_{-m′,-m} = σ H_{m′,m}``.  For integer indices ``σ ≡ 1``; for
half-integer indices ``σ = sgn(m) sgn(m′)``, with ``sgn(0) = 1``; see the notes on the
[``H`` recursion](@ref "Algorithm for computing ``H``").

Which rule applies is settled by the index *type*, so each specialization compiles to a
constant or to two cheap comparisons.
"""
@inline transpose_sign(m′::Integer, m::Integer) = 1
@inline transpose_sign(m′::HalfOddInteger, m::HalfOddInteger) = sgn(m) * sgn(m′)

"""
    wedge_source(m′, m, m′ₘₐₓ)

Return `(a, b, σ)` such that ``H_{m′,m} = σ H_{a,b}``, where `(a, b)` lies in the stored
wedge ``b ≥ |a|``, ``|a| ≤ m′ₘₐₓ``.

This is the definition, for the batched engine, of how the symmetries ``H_{m′,m} =
H_{-m,-m′} = σ H_{m,m′} = σ H_{-m′,-m}`` of the ``H`` matrix supply an element outside the
stored wedge, and a read of such an element should go through it or through
[`wedge_value`](@ref).  The functions that assemble the blocks of the calculators apply the
same cases, in the same order, a whole run of elements at a time, and are tested against
`wedge_value` element by element.  (The dense reference functions apply the integer
symmetries themselves, in [`recurrence_step6!`](@ref).)

An `ArgumentError` is thrown if no stored element can supply the requested one, which
happens only when both ``|m′|`` and ``|m|`` exceed `m′ₘₐₓ`.
"""
@inline function wedge_source(m′::IT, m::IT, m′ₘₐₓ::IT) where {IT}
    if abs(m′) ≤ m′ₘₐₓ
        if m ≥ abs(m′)
            return (m′, m, 1)
        elseif -m ≥ abs(m′)
            return (-m′, -m, transpose_sign(m′, m))
        end
    end
    if abs(m) ≤ m′ₘₐₓ
        if m′ ≥ abs(m)
            return (m, m′, transpose_sign(m′, m))
        elseif -m′ ≥ abs(m)
            return (-m, -m′, 1)
        end
    end
    wedge_source_error(m′, m, m′ₘₐₓ)
end

# Off the hot path, and deliberately `@noinline` so that the formatting of the indices costs
# nothing in the branch that never throws.
@noinline function wedge_source_error(m′, m, m′ₘₐₓ)
    throw(ArgumentError(
        "H[$m′, $m] cannot be obtained from a wedge with "
        * "m′ₘₐₓ=$m′ₘₐₓ; both |m′| and |m| exceed m′ₘₐₓ."
    ))
end

"""
    wedge_offset(H::HWedge, a, b)
    wedge_offset(H::HWedge, a, b, m′ₘᵢₙ)

Zero-based linear offset of the first rotor's element ``H_{a,b}`` in `parent(H)`, for a
stored wedge element (``b ≥ |a|``).  Element `iᵣ` is at `parent(H)[wedge_offset(H, a, b) +
iᵣ]`.

The four-argument form takes `m′ₘᵢₙ(H)` from the caller, which hoists it out of a loop.
"""
function wedge_offset end

@inline function wedge_offset(H::HWedge, a, b, m′ₘᵢₙ)
    @inbounds row_index(H)[(a - m′ₘᵢₙ) + 1] + Nᵣ(H) * (b - abs(a)) - 1
end
@inline wedge_offset(H::HWedge{IT}, a::IT, b::IT) where {IT} =
    wedge_offset(H, a, b, m′ₘᵢₙ(H))

"""
    wedge_value(H::HWedge, iᵣ, m′, m)

Value of ``H_{m′,m}`` for rotor `iᵣ`, for *any* ``|m′|, |m| ≤ ℓ`` (at least one of them
``≤ m′ₘₐₓ``), read from the stored wedge through [`wedge_source`](@ref), which applies the
symmetry sign ``σ`` of half-integer indices.  This is the way to read an element of ``H``
that the wedge does not store, and the only correct one for half-integer indices, where
transposing the stored elements by hand gets the sign of half of them wrong.

`m′` and `m` must be of the wedge's own kind: integers for an integer wedge, and
half-odd-integers — [`HalfOddInteger`](@ref)s, or `Rational`s with denominator 2 such as
`1//2` — for a half-integer one.  An index of the other kind is refused with an
`ArgumentError` rather than floored onto a neighboring element, and an element that the
wedge cannot supply is a `BoundsError` or, when both ``|m′|`` and ``|m|`` exceed `m′ₘₐₓ(H)`,
an `ArgumentError`.
"""
@inline function wedge_value(H::HWedge{IT}, iᵣ::Int, m′::IT, m::IT) where {IT<:IntegerHalf}
    @boundscheck if !(iᵣ > 0 && iᵣ ≤ Nᵣ(H))
        throw(BoundsError(H, (iᵣ, m′, m)))
    end
    a, b, σ = wedge_source(m′, m, m′ₘₐₓ(H))
    @boundscheck if !(abs(a) ≤ b ≤ ℓ(H) && m′ₘᵢₙ(H) ≤ a ≤ m′ₘₐₓ(H))
        throw(BoundsError(H, (iᵣ, m′, m)))
    end
    @inbounds σ * parent(H)[wedge_offset(H, a, b, m′ₘᵢₙ(H)) + iᵣ]
end
# Indices spelled otherwise, such as `1//2` or an `Int8`, are checked against the wedge's
# own kind of index and converted to it; see `check_index_kind` in `wigner_H_calculator.jl`.
@propagate_inbounds function wedge_value(
    H::HWedge{IT}, iᵣ::Integer, m′::IndexType, m::IndexType
) where {IT}
    wedge_value(
        H, Int(iᵣ), calculator_index(IT, m′, "m′", "wedge"), calculator_index(IT, m, "m", "wedge")
    )
end


# Explicit HWedge index formula, assuming no iᵣ:
# (
#     Int(ℓₘᵢₙ - m′ₘᵢₙ) * Int(2ℓ + m′ₘᵢₙ + ℓₘᵢₙ + 1)
#     -
#     Int(ℓₘᵢₙ - m′) * Int(2ℓ - abs(m′ + ℓₘᵢₙ - 1) + 2)
# ) ÷ 2 + 1


"""
    HAxis{IT, RT} <: AbstractWignerMatrix{IT, RT, FixedSizeVectorDefault{RT}}

The `HAxis` type represents the ``m'=0``, ``m≥0`` axis of the `Hˡ` matrix used in
calculation of the Wigner `D` and `d` matrices.
- `IT` is the index type, `Int` or [`HalfOddInteger`](@ref).  The two axes inside an
  [`HCalculator`](@ref) are always `Int`, even for half-integer ``ℓ``, because they hold the
  axis of the integer order ``j = ℓ - 1/2``.
- `RT` is the real type of the elements.

As with [`HWedge`](@ref), the data is stored as a 1-dimensional vector, though it can be
indexed as if it were a two-dimensional array, with the first dimension indexing `Nᵣ`
different rotors, and the second dimension indexing `m` — or alternatively as if it were a
three-dimensional array with the second dimension indexing `m′` and the third indexing `m`.
Thus, this object can be indexed as `Hˡ₀[iᵣ, m]` or `Hˡ[iᵣ, m′, m]` to get the `Hˡ` value
for rotor index `iᵣ`, and matrix element `(m′, m)`.  Its `axes` are those of the first form,
`(1:Nᵣ, m-range)`, while `length` counts the flat storage.

This is scratch space of the recurrence inside an [`HCalculator`](@ref), which holds two of
them, and is not part of the public interface.
"""
mutable struct HAxis{IT, RT} <: AbstractWignerMatrix{IT, RT, FixedSizeVectorDefault{RT}}
    const parent::FixedSizeVectorDefault{RT}
    const Nᵣ::Int
    const maxℓ::IT
    ℓ::IT
    function HAxis(::Type{RT}, Nᵣ::Int, ℓₘₐₓ::IT) where {IT, RT<:Real}
        # An axis starts at ℓₘᵢₙ, and the natural-index accessors check an index only
        # against the current ℓ before reading the storage under `@inbounds`, so the storage
        # must hold at least that one order, for at least one rotor.
        if Nᵣ < 1
            throw(ArgumentError("Number of rotors Nᵣ=$Nᵣ must be at least 1."))
        end
        if ℓₘₐₓ < ℓₘᵢₙ(IT)
            throw(ArgumentError("ℓₘₐₓ=$ℓₘₐₓ must be at least ℓₘᵢₙ=$(ℓₘᵢₙ(IT))."))
        end
        H = FixedSizeVector{RT}(undef, Nᵣ * (Int(ℓₘₐₓ - ℓₘᵢₙ(IT)) + 1))
        new{IT, RT}(H, Nᵣ, ℓₘₐₓ, ℓₘᵢₙ(IT))
    end
end

m′ₘₐₓ(w::HAxis{IT}) where {IT} = ℓₘᵢₙ(IT)
m′ₘᵢₙ(w::HAxis{IT}) where {IT} = ℓₘᵢₙ(IT)
mₘₐₓ(w::HAxis{IT}) where {IT} = ℓ(w)
mₘᵢₙ(w::HAxis{IT}) where {IT} = ℓₘᵢₙ(IT)

Nᵣ(w::HAxis{IT}) where {IT} = w.Nᵣ
maxℓ(w::HAxis{IT}) where {IT} = w.maxℓ

function Base.setproperty!(H::HAxis{IT}, s::Symbol, ℓ::IIT) where {IT, IIT}
    if s === :ℓ
        if IIT !== IT
            throw(ArgumentError(
                "Cannot change ℓ from type $IT to type $IIT; they must be the same."
            ))
        end
        if ℓ < ℓₘᵢₙ(IT)
            throw(ArgumentError("Cannot set ℓ=$ℓ less than ℓₘᵢₙ=$(ℓₘᵢₙ(IT))."))
        end
        if ℓ > maxℓ(H)
            throw(ArgumentError("Cannot set ℓ=$ℓ greater than maxℓ=$(maxℓ(H))."))
        end
        Base.setfield!(H, :ℓ, ℓ)
        ℓ
    else
        throw(ArgumentError(
            "Cannot set property `$s` on HAxis; only `ℓ` is allowed to be changed."
        ))
    end
end

# An axis is indexed `[iᵣ, m]`; as for `HWedge`, `length` and `size(w)` describe the flat
# storage.
Base.axes(w::HAxis{IT}) where {IT} = (1:Nᵣ(w), WignerRange(ℓₘᵢₙ(w):ℓ(w)))
Base.size(w::HAxis, d::Integer) = length(axes(w, d))

function Base.checkbounds(::Type{Bool}, w::HAxis, i::Int)
    i ≥ 1 && i ≤ length(w)
end
function Base.checkbounds(::Type{Bool}, w::HAxis{IT}, iᵣ::Int, m::IT) where {IT}
    iᵣ > 0 && iᵣ ≤ Nᵣ(w) && ℓₘᵢₙ(w) ≤ m ≤ ℓ(w)
end
function Base.checkbounds(::Type{Bool}, w::HAxis{IT}, iᵣ::Int, m′::IT, m::IT) where {IT}
    iᵣ > 0 && iᵣ ≤ Nᵣ(w) && m′ == ℓₘᵢₙ(w) && ℓₘᵢₙ(w) ≤ m ≤ ℓ(w)
end

@propagate_inbounds function Base.getindex(w::HAxis, i::Int)
    @boundscheck if !checkbounds(Bool, w, i)
        throw(BoundsError(w, i))
    end
    @inbounds Base.parent(w)[i]
end
@propagate_inbounds function Base.getindex(w::HAxis{IT}, iᵣ::Int, m::IT) where {IT}
    @boundscheck if !checkbounds(Bool, w, iᵣ, m)
        throw(BoundsError(w, (iᵣ, m)))
    end
    i = iᵣ + Nᵣ(w) * Int(m - ℓₘᵢₙ(w))
    @inbounds Base.parent(w)[i]
end
@propagate_inbounds function Base.getindex(w::HAxis{IT}, iᵣ::Int, m′::IT, m::IT) where {IT}
    @boundscheck if !checkbounds(Bool, w, iᵣ, m′, m)
        throw(BoundsError(w, (iᵣ, m′, m)))
    end
    i = iᵣ + Nᵣ(w) * Int(m - ℓₘᵢₙ(w))
    @inbounds Base.parent(w)[i]
end

@propagate_inbounds function Base.setindex!(w::HAxis, v, i::Int)
    @boundscheck if !checkbounds(Bool, w, i)
        throw(BoundsError(w, i))
    end
    @inbounds Base.parent(w)[i] = v
end
@propagate_inbounds function Base.setindex!(w::HAxis{IT}, v, iᵣ::Int, m::IT) where {IT}
    @boundscheck if !checkbounds(Bool, w, iᵣ, m)
        throw(BoundsError(w, (iᵣ, m)))
    end
    i = iᵣ + Nᵣ(w) * Int(m - ℓₘᵢₙ(w))
    @inbounds Base.parent(w)[i] = v
end
@propagate_inbounds function Base.setindex!(w::HAxis{IT}, v, iᵣ::Int, m′::IT, m::IT) where {IT}
    @boundscheck if !checkbounds(Bool, w, iᵣ, m′, m)
        throw(BoundsError(w, (iᵣ, m′, m)))
    end
    i = iᵣ + Nᵣ(w) * Int(m - ℓₘᵢₙ(w))
    @inbounds Base.parent(w)[i] = v
end

# As for `HWedge`: equal when they describe the same ℓ and rotors, and agree on every
# element
function Base.:(==)(w1::HAxis{IT}, w2::HAxis{IT}) where {IT}
    ℓ(w1) == ℓ(w2) && Nᵣ(w1) == Nᵣ(w2) &&
        all(w1[iᵣ, m] == w2[iᵣ, m] for m ∈ ℓₘᵢₙ(w1):ℓ(w1) for iᵣ ∈ 1:Nᵣ(w1))
end

function Base.summary(io::IO, H::HAxis{IT, RT}) where {IT, RT}
    print(io, "HAxis{$IT, $RT} for ℓ=$(ℓ(H)) with m=$(ℓₘᵢₙ(H)):$(ℓ(H)), and iᵣ=1:$(Nᵣ(H))")
end
Base.show(io::IO, H::HAxis) = summary(io, H)
function Base.show(io::IO, ::MIME"text/plain", H::HAxis{IT, RT}) where {IT, RT}
    summary(io, H)
    print(io, "\nStored in ")
    show(io, MIME("text/plain"), parent(H))
end
