

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

However, for further efficiency when vectorizing and threading over multiple rotations, the
data is stored as a 1-dimensional vector, though it can be indexed as if it were a
three-dimensional array, with the first dimension indexing `Nᵣ` different rotations, and the
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

function HWedge_row_index!(row_index, Nᵣ::Int, ℓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT) where {IT}
    index = 1
    for (i, m′) ∈ enumerate(m′ₘᵢₙ:m′ₘₐₓ)
        @inbounds row_index[i] = index
        index += Nᵣ * (Int(ℓ - abs(m′)) + 1)
    end
    row_index
end

function HWedge_size(ℓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT) where {IT}
    let ℓₘᵢₙ = ℓₘᵢₙ(IT)
        Int(
            (ℓₘᵢₙ - m′ₘᵢₙ) * (2ℓ + m′ₘᵢₙ + ℓₘᵢₙ + 1)
            - (ℓₘᵢₙ - m′ₘₐₓ - 1) * (2ℓ - m′ₘₐₓ - ℓₘᵢₙ + 2)
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
    iᵣ > 0 && iᵣ ≤ Nᵣ(w) && m ≥ abs(m′) && m′ ≥ m′ₘᵢₙ(w) && m′ ≤ m′ₘₐₓ(w)
end

@propagate_inbounds function Base.getindex(w::HWedge, i::Int)
    @boundscheck if !checkbounds(Bool, w, i)
        throw(BoundsError(w, i))
    end
    @inbounds Base.parent(w)[i]
end
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

function Base.show(io::IO, ::MIME"text/plain", H::HWedge{IT, RT, ST}) where {IT, RT, ST}
    let ℓ = ℓ(H), m′ₘᵢₙ = m′ₘᵢₙ(H), m′ₘₐₓ = m′ₘₐₓ(H), Nᵣ = Nᵣ(H)
        print(
            io,
            "SphericalFunctions.HWedge{$IT, $RT} for ℓ=$(ℓ) with m′=$(m′ₘᵢₙ:m′ₘₐₓ), ",
            "m=abs(m′):$(ℓ), and iᵣ=1:$(Nᵣ)\n",
            "Stored in ",
        )
        show(io, MIME("text/plain"), parent(H))
    end
end


# """

# An Hˡ wedge will store elements in a vector as if it were the following matrix:

#     [
#         H[ℓ, m′, m]
#         for m′ ∈ max(-ℓ, m′ₘᵢₙ):min(ℓ, m′ₘₐₓ)
#         for m ∈ abs(m′):ℓ
#     ]

# Here, m′ₘᵢₙ is a negative number and m′ₘₐₓ is a positive number.  Note that for HWedge, we
# currently impose mₘₐₓ = ℓ and mₘᵢₙ = ℓₘᵢₙ(IT), because these are all needed for the
# recurrence relations.

# This function returns the linear index into that vector that belongs to the first element
# with the given `m′` value (and therefore `m=abs(m′)`).  The formula for that index involves
# an `if` statement to account for the varying number of `m` values for each `m′` value.
# Nonetheless, it can be computed in closed form (i.e., without an explicit sum or loop).


# """



# function row_index(w::HWedge{IT}, m′::IT) where {IT}
#     let ℓ = ℓ(w), m′ₘᵢₙ = m′ₘᵢₙ(w), ℓₘᵢₙ = ℓₘᵢₙ(IT)
#         (
#             Int(ℓₘᵢₙ - m′ₘᵢₙ) * Int(2ℓ + m′ₘᵢₙ + ℓₘᵢₙ + 1)
#             -
#             Int(ℓₘᵢₙ - m′) * Int(2ℓ - abs(m′ + ℓₘᵢₙ - 1) + 2)
#         ) ÷ 2 + 1

#         # i = if m′<1
#         #     Int(m′ - m′ₘᵢₙ) * Int(2ℓ + m′ + m′ₘᵢₙ + 1) ÷ 2  # size of wedge to the left of m'
#         # else
#         #     (
#         #         # size of entire left half of wedge
#         #         Int(ℓₘᵢₙ - m′ₘᵢₙ) * Int(2ℓ + ℓₘᵢₙ + m′ₘᵢₙ + 1)
#         #         +
#         #         # size of right half of wedge to the left of m'
#         #         Int(m′ - ℓₘᵢₙ) * Int(2ℓ - ℓₘᵢₙ - m′ + 3)
#         #     ) ÷ 2
#         # end
#         # i + 1
#     end
# end


# function row_index(ℓ::IT, m′::IT) where {IT}
#     let ℓₘᵢₙ = ℓₘᵢₙ(IT)
#         i = if m′<ℓₘᵢₙ
#             # size of wedge above m′
#             Int(m′ - m′ₘᵢₙ) * Int(2ℓ + m′ + m′ₘᵢₙ + 1) ÷ 2
#         else
#             (
#                 # size of entire upper half of wedge excluding m′=ℓₘᵢₙ
#                 Int(ℓₘᵢₙ - m′ₘᵢₙ) * Int(2ℓ + ℓₘᵢₙ + m′ₘᵢₙ + 1)
#                 +
#                 # size of wedge at or below m′=ℓₘᵢₙ but above m′
#                 Int(m′ - ℓₘᵢₙ) * Int(2ℓ - ℓₘᵢₙ - m′ + 3)
#             ) ÷ 2
#         end
#         i + 1
#     end
# end

# function row_index(ℓ::IT, m′::IT, m′ₘᵢₙ::IT) where {IT}
#     let ℓₘᵢₙ = ℓₘᵢₙ(IT)
#         # size of entire upper half of wedge excluding m′=ℓₘᵢₙ
#         zero_index = Int(ℓₘᵢₙ - m′ₘᵢₙ) * Int(2ℓ + ℓₘᵢₙ + m′ₘᵢₙ + 1)

#         i = if m′<ℓₘᵢₙ
#             (
#                 zero_index
#                 +
#                 # size of wedge at or below m′ but above m′=ℓₘᵢₙ
#                 Int(m′ - ℓₘᵢₙ) * Int(2ℓ - abs(m′ + ℓₘᵢₙ - 1) + 2)
#             ) ÷ 2
#         else
#             (
#                 zero_index
#                 +
#                 # size of wedge at or below m′=ℓₘᵢₙ but above m′
#                 Int(m′ - ℓₘᵢₙ) * Int(2ℓ - abs(m′ + ℓₘᵢₙ - 1) + 2)
#             ) ÷ 2
#         end
#         i + 1
#     end
# end

# function row_index(ℓ::IT, m′::IT) where {IT}
#     let ℓₘᵢₙ = ℓₘᵢₙ(IT)#, m′ₘᵢₙ = -ℓ
#         # Size of upper half (m′ₘᵢₙ to ℓₘᵢₙ-1)
#         zero_index = Int(ℓₘᵢₙ - m′ₘᵢₙ) * Int(2ℓ + ℓₘᵢₙ + m′ₘᵢₙ + 1)

#         # Correction term (works for both m′ < ℓₘᵢₙ and m′ ≥ ℓₘᵢₙ)
#         correction = Int(m′ - ℓₘᵢₙ) * Int(2ℓ - ℓₘᵢₙ - m′ + 3)
        
#         # For m′ < ℓₘᵢₙ: correction is negative → subtract unwanted rows
#         # For m′ ≥ ℓₘᵢₙ: correction is positive → add needed rows
#         i = (zero_index + correction) ÷ 2
#         i + 1
#     end
# end
