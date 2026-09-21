### Containers laid out in the canonical mode ordering.
#
# Two things in this package are stored as `[x(ℓ, m) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ]` — the
# weights of a spin-weighted function, and the values of the spin-weighted harmonics
# themselves — and they want the same indexing, the same accessors and the same flat storage.
# They differ only in meaning, which is why they keep separate names rather than one type
# serving both.  The shared supertype is what lets the machinery be written once.
#
# The mode axis is always the *last* axis of the storage, so that a single ``ℓ`` is a view over
# a contiguous run of it with every leading axis taken whole.  That is what keeps the flat form
# usable for the products these containers exist to feed: a synthesis matrix times a vector of
# mode weights.

"""
    AbstractModeContainer{T, IT}

Supertype of the containers stored in the canonical mode ordering — [`ModeWeights`](@ref) and
[`HarmonicValues`](@ref).  `T` is the number type and `IT` the index type (an `Integer` or a
[`HalfOddInteger`](@ref)).

These are *not* `AbstractArray`s.  A container indexed by ``ℓ`` cannot be one, because ``ℓ``
may be a half-odd-integer and `axes` must be integer ranges; and for the ones that could be,
being an array is what let the `OffsetArray`s of earlier versions accept `*` and return
silently wrong answers.  [`strided`](@ref) is the explicit route to the flat 1-based storage,
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
# is the same arithmetic for either kind of index.
@inline function mode_range(c::AbstractModeContainer, ℓ)
    i₀ = Yindex(ℓ, -ℓ, ℓₘᵢₙ(c))
    i₀:(i₀ + Int(2ℓ))
end

# Shared by `getindex(c, ℓ)` on every such container: `ℓ` must be one of the values the
# container holds, and a whole number asked of a half-integer container should be told what it
# holds rather than shown a bare `InexactError` from `convert`.
@inline function check_ℓ(c::AbstractModeContainer{T, IT}, ℓ) where {T, IT}
    if !isindex(IT, ℓ)
        throw(ArgumentError(
            "ℓ=$ℓ is not one of the ℓ values of this container, which runs over "
            * "ℓ ∈ $(ℓₘᵢₙ(c)):$(ℓₘₐₓ(c)) in steps of 1."
        ))
    end
    ℓ′ = convert(IT, ℓ)
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

[`strided`](@ref) gives the flat storage: a `Vector` of modes, or an array whose *last* axis is
the modes in the canonical ordering (see [`Yindex`](@ref)) and whose leading axes are the
rotors and spin weights.  That is the form a product with mode weights takes, and
[`sYlm_matrix`](@ref) is the direct spelling of it.

`spins(Y)` is the range of spin weights served, `spin(Y)` the single value when there is only
one, `Nᵣ(Y)` the number of rotors, and `ℓₘᵢₙ(Y)`/`ℓₘₐₓ(Y)` the range of ``ℓ``.  Iterating gives
`ℓ => block` pairs, as a calculator does.

See also [`ModeWeights`](@ref), which shares this layout but holds the weights of a function
rather than the values of the harmonics.
"""
struct HarmonicValues{T, IT<:IntegerHalf, S, A<:AbstractArray{T}} <: AbstractModeContainer{T, IT}
    data::A
    s::S          # one `IT`, or an ascending range of them
    ℓₘᵢₙ::IT
    ℓₘₐₓ::IT
    Nᵣ::Int       # 1 when the container was built for a single rotor

    function HarmonicValues(
        data::A, s::S, ℓₘᵢₙ::IT, ℓₘₐₓ::IT, Nᵣ::Int
    ) where {T, IT<:IntegerHalf, S, A<:AbstractArray{T}}
        Base.require_one_based_indexing(data)
        if size(data)[end] != Ysize(ℓₘᵢₙ, ℓₘₐₓ)
            throw(ArgumentError(
                "The mode axis has length $(size(data)[end]), but "
                * "Ysize(ℓₘᵢₙ=$ℓₘᵢₙ, ℓₘₐₓ=$ℓₘₐₓ) = $(Ysize(ℓₘᵢₙ, ℓₘₐₓ))."
            ))
        end
        new{T, IT, S, A}(data, s, ℓₘᵢₙ, ℓₘₐₓ, Nᵣ)
    end
end

Base.parent(Y::HarmonicValues) = Y.data
Nᵣ(Y::HarmonicValues) = Y.Nᵣ
isbatched(Y::HarmonicValues) = Y.Nᵣ > 1
spins(Y::HarmonicValues{T, IT, S}) where {T, IT, S<:IntegerHalf} = Y.s:Y.s
spins(Y::HarmonicValues{T, IT, S}) where {T, IT, S<:AbstractUnitRange} = Y.s
spin(Y::HarmonicValues{T, IT, S}) where {T, IT, S<:IntegerHalf} = Y.s

# `length` counts the blocks, as it does for a `WignerSeries`; the number of modes is
# `length(strided(Y))` for the unbatched single-spin case, and `Ysize` in general.
Base.length(Y::HarmonicValues) = Int(ℓₘₐₓ(Y) - ℓₘᵢₙ(Y)) + 1
Base.keys(Y::HarmonicValues) = ℓₘᵢₙ(Y):ℓₘₐₓ(Y)
Base.firstindex(Y::HarmonicValues) = ℓₘᵢₙ(Y)
Base.lastindex(Y::HarmonicValues) = ℓₘₐₓ(Y)

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

# Iteration yields `ℓ => block`, matching the calculators, so that a loop written against one
# reads the same against the other.
@inline function Base.iterate(Y::HarmonicValues{T, IT}, ℓ::IT=ℓₘᵢₙ(Y)) where {T, IT}
    ℓ > ℓₘₐₓ(Y) && return nothing
    (ℓ => Y[ℓ], ℓ + 1)
end
Base.IteratorSize(::Type{<:HarmonicValues}) = Base.HasLength()
Base.pairs(Y::HarmonicValues) = Y

Base.copy(Y::HarmonicValues) = HarmonicValues(copy(Y.data), Y.s, Y.ℓₘᵢₙ, Y.ℓₘₐₓ, Y.Nᵣ)
function Base.:(==)(a::HarmonicValues, b::HarmonicValues)
    a.s == b.s && a.ℓₘᵢₙ == b.ℓₘᵢₙ && a.ℓₘₐₓ == b.ℓₘₐₓ && a.Nᵣ == b.Nᵣ && a.data == b.data
end

function Base.show(io::IO, Y::HarmonicValues{T, IT, S}) where {T, IT, S}
    spin_text = S <: AbstractUnitRange ? "s ∈ $(Y.s)" : "s = $(Y.s)"
    rotor_text = Y.Nᵣ == 1 ? "" : ", $(Y.Nᵣ) rotors"
    print(io, "HarmonicValues{$T} for ℓ ∈ $(Y.ℓₘᵢₙ):$(Y.ℓₘₐₓ), $spin_text$rotor_text")
end
function Base.show(io::IO, ::MIME"text/plain", Y::HarmonicValues)
    show(io, Y)
    println(io, ":")
    for ℓ ∈ ℓₘᵢₙ(Y):ℓₘₐₓ(Y)
        println(io, " ℓ = ", ℓ, ":")
        show(io, MIME("text/plain"), Y[ℓ])
        println(io)
    end
end
