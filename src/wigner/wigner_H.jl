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

``(ℓ - m)(ℓ + m + 1)``, the square of Gumerov and Duraiswami's ``d^ℓ_m`` coefficient (up to
its sign; see [`sgn`](@ref)).  Both factors are `Integer`s for integer and half-integer
indices alike, so this is exact, and the value handed to `sqrt` is the same one a purely
integer implementation computes.
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
that keeps the ``H``-form of the ``m′`` ladder in Gumerov and Duraiswami's shape (see the v3
design memo, §5.2).  The single expression below serves both index types.
"""
@inline ϵ(m) = ifelse(m > 0 && isodd(floor(Int, m)), -1, 1)

"""
    HWedge{IT, RT, ST} <: AbstractWignerMatrix{IT, RT, ST}

The ``Hˡ`` matrix is critical to efficient and stable computation of the Wigner ``D`` and
``d`` matrices — in fact, it essentially *is* the ``d`` matrix with signs adjusted to avoid
numerical problems with alternating signs.  This gives it additional symmetries that reduce
the amount of data that needs to be stored to about 1/4 of the total ``d`` size.

The purpose of an `HWedge` is to provide a workspace for the Wigner recurrences that is
efficient, both in terms of the size of memory used, and the implications for vectorization
and threading.  Specifically, the data is stored as strictly `Real` values, in contiguous
storage.  Indexing is performed efficiently via precomputed row offsets.  Once the full
recurrence is done, the data can be used directly — computing phases and symmetry on the fly
— or copied into a full explicit matrix with the appropriate phases.

The recurrences require ``m`` in the full range from 0 (or 1/2) to ``ℓ``, but ``m'`` only
needs to include the axis ``m'=0`` or 1/2.  Thus, we store `m′ₘᵢₙ`and `m′ₘₐₓ` as fields, and
only require enough storage for those ranges.  Specifically, an `HWedge` will store elements
in a vector as if they were components of the `Hˡ` matrix:

    [
        Hˡ[m′, m]
        for m′ ∈ max(-ℓ, m′ₘᵢₙ):min(ℓ, m′ₘₐₓ)
        for m ∈ abs(m′):ℓ
    ]

However, for further efficiency when vectorizing and threading over multiple rotors, the
data is stored as a 1-dimensional vector, though it can be indexed as if it were a
three-dimensional array, with the first dimension indexing `Nᵣ` different rotors, and the
second dimension indexing `m′`, and the third dimension indexing `m`.  Thus, this object can
be indexed as `Hˡ[iᵣ, m′, m]` to get the `Hˡ` value for rotor index `iᵣ`, and matrix element
`(m′, m)`.

Because of this complicated layout, the constructor is fairly restrictive, but will do all
the allocation needed.  To avoid multiple allocations, it is advisable to first construct an
instance with the maximum `ℓ` value that will be needed, and then change the `ℓ` field as
needed to compute different orders.  That is, if `H isa HWedge`, then `H.ℓ = new_ell` can be
used to change the current order being computed.  The constructor starts out with the
smallest `ℓ` value possible (0 or 1/2), which is the natural choice for recurrence.

!!! warning "Thread safety"

    The `HWedge` object is not thread safe.  Its internal storage is intended to be changed
    by different threads, but the code must be designed carefully to avoid accessing the
    same memory locations from different threads.  In particular, note that changing the `ℓ`
    field changes internal storage, and is *not* thread-safe.  It may be better to allocate
    separate `HWedge` objects for each thread.

"""
mutable struct HWedge{IT, RT<:Real, ST} <: AbstractWignerMatrix{IT, RT, ST}
    const parent::ST
    const row_index::FixedSizeVectorDefault{Int}
    const Nᵣ::Int
    const maxℓ::IT
    const maxm′ₘₐₓ::IT
    const minm′ₘᵢₙ::IT
    ℓ::IT
    m′ₘₐₓ::IT
    m′ₘᵢₙ::IT
    function HWedge(Nᵣ::Int, ℓₘₐₓ::IT, m′ₘₐₓ::IT=ℓₘₐₓ, m′ₘᵢₙ::IT=-ℓₘₐₓ) where {IT}
        HWedge(Float64, Nᵣ, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ)
    end
    function HWedge(::Type{RT}, Nᵣ::Int, ℓₘₐₓ::IT, m′ₘₐₓ::IT=ℓₘₐₓ, m′ₘᵢₙ::IT=-ℓₘₐₓ) where {IT, RT<:Real}
        if Nᵣ < 1
            error("Number of rotors Nᵣ=$Nᵣ must be at least 1.")
        end
        validate_index_ranges(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ)

        # Set up storage for the biggest these values will ever be
        parent = FixedSizeVector{RT}(undef, Nᵣ * HWedge_size(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ))
        row_index = FixedSizeVector{Int}(undef, Int(m′ₘₐₓ - m′ₘᵢₙ) + 1)

        # But start out assuming ℓ is the smallest it can be
        maxℓ = ℓₘₐₓ
        maxm′ₘₐₓ = m′ₘₐₓ
        minm′ₘᵢₙ = m′ₘᵢₙ
        ℓ = ℓₘᵢₙ(ℓₘₐₓ)
        m′ₘₐₓ = min(ℓ, m′ₘₐₓ)
        m′ₘᵢₙ = max(-ℓ, m′ₘᵢₙ)
        HWedge_row_index!(row_index, Nᵣ, ℓ, m′ₘₐₓ, m′ₘᵢₙ)

        new{IT, RT, typeof(parent)}(
            parent, row_index, Nᵣ, maxℓ, maxm′ₘₐₓ, minm′ₘᵢₙ, ℓ, m′ₘₐₓ, m′ₘᵢₙ
        )
    end
end

mₘₐₓ(w::HWedge{IT}) where {IT} = ℓ(w)
mₘᵢₙ(w::HWedge{IT}) where {IT} = ℓₘᵢₙ(w)

row_index(w::HWedge{IT}) where {IT} = w.row_index
row_index(w::HWedge{IT}, m′::IT) where {IT} = row_index(w)[Int(m′ - m′ₘᵢₙ(w)) + 1]
Nᵣ(w::HWedge{IT}) where {IT} = w.Nᵣ
maxℓ(w::HWedge{IT}) where {IT} = w.maxℓ
maxm′ₘₐₓ(w::HWedge{IT}) where {IT} = w.maxm′ₘₐₓ
minm′ₘᵢₙ(w::HWedge{IT}) where {IT} = w.minm′ₘᵢₙ

function Base.setproperty!(H::HWedge{IT}, s::Symbol, ℓ::IIT) where {IT, IIT}
    if s === :ℓ
        if IIT !== IT
            error("Cannot change ℓ from type $IT to type $IIT; they must be the same.")
        end
        if ℓ < ℓₘᵢₙ(IT)
            error("Cannot set ℓ=$ℓ less than ℓₘᵢₙ=$(ℓₘᵢₙ(IT)).")
        end
        if ℓ > maxℓ(H)
            error("Cannot set ℓ=$ℓ greater than maxℓ=$(maxℓ(H)).")
        end
        m′ₘₐₓ = min(ℓ, maxm′ₘₐₓ(H))
        m′ₘᵢₙ = max(-ℓ, minm′ₘᵢₙ(H))
        HWedge_row_index!(row_index(H), Nᵣ(H), ℓ, m′ₘₐₓ, m′ₘᵢₙ)
        Base.setfield!(H, :ℓ, ℓ)
        Base.setfield!(H, :m′ₘₐₓ, m′ₘₐₓ)
        Base.setfield!(H, :m′ₘᵢₙ, m′ₘᵢₙ)
        ℓ
    else
        error("Cannot set property `$s` on HWedge; only `ℓ` is allowed to be changed.")
    end
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

function Base.axes(w::HWedge{IT}) where {IT}
    (1:Nᵣ(w), WignerRange(m′ₘᵢₙ(w):m′ₘₐₓ(w)), WignerRange(mₘᵢₙ(w):mₘₐₓ(w)))
end

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
# See the note on `Rational` indexing in `wigner_matrix.jl`.
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
        show(io, MIME("text/plain"), parent(H)[begin:i])
    end
end

"""
    transpose_sign(m′, m)

Sign ``σ`` relating the transposed element of the ``H`` matrix to the original:
``H_{m,m′} = σ H_{m′,m}`` and ``H_{-m′,-m} = σ H_{m′,m}``.  For integer indices
``σ ≡ 1``; for half-integer indices ``σ = sgn(m) sgn(m′)``, with ``sgn(0) = 1``; see the
notes on the [``H`` recursion](@ref "Algorithm for computing ``H`` (redesigned)").

Which rule applies is settled by the index *type*, so each specialization compiles to a
constant or to two cheap comparisons.
"""
@inline transpose_sign(m′::Integer, m::Integer) = 1
@inline transpose_sign(m′::HalfOddInteger, m::HalfOddInteger) = sgn(m) * sgn(m′)

"""
    wedge_source(m′, m, m′ₘₐₓ)

Return `(a, b, σ)` such that ``H_{m′,m} = σ H_{a,b}``, where `(a, b)` lies in the stored
wedge ``b ≥ |a|``, ``|a| ≤ m′ₘₐₓ``.

This is the *only* place in the package that encodes the symmetries
``H_{m′,m} = H_{-m,-m′} = σ H_{m,m′} = σ H_{-m′,-m}`` of the ``H`` matrix; every read of an
element outside the stored wedge must go through it.

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
stored wedge element (``b ≥ |a|``).  Element `iᵣ` is at
`parent(H)[wedge_offset(H, a, b) + iᵣ]`.

The four-argument form takes `m′ₘᵢₙ(H)` from the caller, which hoists it out of a loop.
"""
@inline function wedge_offset(H::HWedge, a, b, m′ₘᵢₙ)
    @inbounds row_index(H)[(a - m′ₘᵢₙ) + 1] + Nᵣ(H) * (b - abs(a)) - 1
end
@inline wedge_offset(H::HWedge{IT}, a::IT, b::IT) where {IT} =
    wedge_offset(H, a, b, m′ₘᵢₙ(H))

"""
    wedge_value(H::HWedge, iᵣ, m′, m)

Value of ``H_{m′,m}`` for rotor `iᵣ`, for *any* ``|m′|, |m| ≤ ℓ`` (at least one of them
``≤ m′ₘₐₓ``), read from the stored wedge through [`wedge_source`](@ref).

`m′` and `m` are of the wedge's own index type, so a wrong-parity index — a whole number for
a half-integer wedge, say — cannot be expressed, let alone silently floored onto a
neighbouring element.
"""
@inline function wedge_value(H::HWedge{IT}, iᵣ::Int, m′::IT, m::IT) where {IT}
    @boundscheck if !(iᵣ > 0 && iᵣ ≤ Nᵣ(H))
        throw(BoundsError(H, (iᵣ, m′, m)))
    end
    a, b, σ = wedge_source(m′, m, m′ₘₐₓ(H))
    @boundscheck if !(abs(a) ≤ b ≤ ℓ(H) && m′ₘᵢₙ(H) ≤ a ≤ m′ₘₐₓ(H))
        throw(BoundsError(H, (iᵣ, m′, m)))
    end
    @inbounds σ * parent(H)[wedge_offset(H, a, b, m′ₘᵢₙ(H)) + iᵣ]
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

As with [`HWedge`](@ref), the data is stored as a 1-dimensional vector, though it can be
indexed as if it were a two-dimensional array, with the first dimension indexing `Nᵣ`
different rotors, and the second dimension indexing `m` — or alternatively as if it were a
three-dimensional array with the second dimension indexing `m′` and the third indexing `m`.
Thus, this object can be indexed as `Hˡ₀[iᵣ, m]` or `Hˡ[iᵣ, m′, m]` to get the `Hˡ` value
for rotor index `iᵣ`, and matrix element `(m′, m)`.

"""
mutable struct HAxis{IT, RT} <: AbstractWignerMatrix{IT, RT, FixedSizeVectorDefault{RT}}
    const parent::FixedSizeVectorDefault{RT}
    const Nᵣ::Int
    const maxℓ::IT
    ℓ::IT
    function HAxis(::Type{RT}, Nᵣ::Int, ℓₘₐₓ::IT) where {IT, RT<:Real}
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
            error("Cannot change ℓ from type $IT to type $IIT; they must be the same.")
        end
        if ℓ < ℓₘᵢₙ(IT)
            error("Cannot set ℓ=$ℓ less than ℓₘᵢₙ=$(ℓₘᵢₙ(IT)).")
        end
        if ℓ > maxℓ(H)
            error("Cannot set ℓ=$ℓ greater than maxℓ=$(maxℓ(H)).")
        end
        Base.setfield!(H, :ℓ, ℓ)
        ℓ
    else
        error("Cannot set property `$s` on HAxis; only `ℓ` is allowed to be changed.")
    end
end

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

function Base.summary(io::IO, H::HAxis{IT, RT}) where {IT, RT}
    print(io, "HAxis{$IT, $RT} for ℓ=$(ℓ(H)) with m=$(ℓₘᵢₙ(H)):$(ℓ(H)), and iᵣ=1:$(Nᵣ(H))")
end
Base.show(io::IO, H::HAxis) = summary(io, H)
function Base.show(io::IO, ::MIME"text/plain", H::HAxis{IT, RT}) where {IT, RT}
    summary(io, H)
    print(io, "\nStored in ")
    show(io, MIME("text/plain"), parent(H))
end
