import Base: @propagate_inbounds

"""
    AbstractWignerMatrix{IT, NT, ST}

Abstract base type for Wigner rotation‐matrix objects of a specific ``ℓ`` value.
- `IT` is the index type (an `Integer` or half-integer `Rational`), governing the allowed
  ranges of `m′` and `m`.
- `NT` is the number type (e.g., `ComplexF64` for D-matrices or `Float64` for d-matrices).
- `ST` is the storage type (typically `Matrix{NT}`, but other `AbstractMatrix{NT}` storage
  can be used).

The basic concrete subtypes (`WignerDMatrix`, `WignerdMatrix`) default to storing their data
in a `Matrix{NT}` and implement the usual `size`, `getindex` and `setindex!` so that one can
use `w[m′,m]`.  Specifically, these indices can be negative or positive, and must obey
`abs(m′) ≤ m′ₘₐₓ` and `abs(m) ≤ ℓ`.

# Methods

Methods defined for `AbstractWignerMatrix` objects include:
- `parent(w)`: the underlying data array.
- `ℓ(w)` or `ell(w)`: the value of ``ℓ``.
- `m′ₘₐₓ(w)` or `mpmax(w)`: the maximum value of ``m′``.
- `mₘₐₓ(w)` or `mmax(w)`: the maximum value of ``m``.
- `ℓₘᵢₙ(w)` or `ellmin(w)`: the minimum value of ``ℓ``, which is either 0 or 1//2.
- `isrational(w)`: whether the indices are rational (i.e., half‐integer).
- `size(w)`: the size of the underlying data array.
- `length(w)`: the length of the underlying data array.
- `getindex(w, i)`: get the value at index `i` in the underlying data array.
- `getindex(w, m′, m)`: get the value at index `(m′, m)`.
- `setindex!(w, v, i)`: set the value at index `i` in the underlying data array to `v`.
- `setindex!(w, v, m′, m)`: set the value at index `(m′, m)`.
- `axes(w)`: the axes of the matrix, which are 2-tuples of ranges for the `m′` and `m`
  indices.

# Implementation

Any new subtypes of `AbstractWignerMatrix` should inherit from this type and re-implement any of the
methods mentioned above that are not appropriate for the new type.  Specifically, the
default implementations assume that subtypes store the fields
- `parent::ST`: the underlying storage type.
- `ℓ::IT`: the value of ``ℓ``.
- `m′ₘₐₓ::IT`: the maximum value of ``m′``.

For example, if the parent Matrix is not stored as the `parent` field, then the `parent(w)`
method should be re-implemented to return the correct parent object.  The `getindex` and
`setindex!`
"""
abstract type AbstractWignerMatrix{IT<:Union{Integer,Rational}, NT, ST<:AbstractArray{NT}} <: AbstractMatrix{NT} end

### General methods for all AbstractWignerMatrix types

Base.parent(w::AbstractWignerMatrix) = w.parent

ℓ(w::AbstractWignerMatrix{IT}) where {IT} = w.ℓ
ℓₘᵢₙ(::IT) where {IT} = ℓₘᵢₙ(IT)
ℓₘᵢₙ(::Type{IT}) where {IT} = error("No method defined for `ℓₘᵢₙ(::Type{$IT})`.")
ℓₘᵢₙ(::Type{IT}) where {IT<:Integer} = zero(IT)
ℓₘᵢₙ(::Type{IT}) where {IT<:Rational} = IT(1//2)
ℓₘᵢₙ(::AbstractWignerMatrix{IT}) where {IT} = ℓₘᵢₙ(IT)

m′ₘₐₓ(w::AbstractWignerMatrix{IT}) where {IT} = w.m′ₘₐₓ
m′ₘᵢₙ(w::AbstractWignerMatrix{IT}) where {IT} = w.m′ₘᵢₙ
mₘₐₓ(w::AbstractWignerMatrix{IT}) where {IT} = w.mₘₐₓ
mₘᵢₙ(w::AbstractWignerMatrix{IT}) where {IT} = w.mₘᵢₙ

const ell = ℓ
const ellmin = ℓₘᵢₙ
const mpmax = m′ₘₐₓ
const mpmin = m′ₘᵢₙ
const mmax = mₘₐₓ
const mmin = mₘᵢₙ

isrational(::AbstractWignerMatrix{IT}) where {IT<:Integer} = false
isrational(::AbstractWignerMatrix{IT}) where {IT<:Rational} = true

Base.eltype(::AbstractWignerMatrix{IT, NT, ST}) where {IT, NT, ST} = NT
Base.size(w::AbstractWignerMatrix{IT, NT, ST}) where {IT, NT, ST} = size(parent(w))
Base.length(w::AbstractWignerMatrix{IT, NT, ST}) where {IT, NT, ST} = length(parent(w))
# function Base.axes(w::AbstractWignerMatrix{IT}) where {IT}
#     ((m′ₘᵢₙ(w):m′ₘₐₓ(w)), (mₘᵢₙ(w):mₘₐₓ(w)))
# end

struct WignerRange{T<:Union{Integer,Rational}} <: AbstractUnitRange{T}
    start::T
    stop::T

    WignerRange(r::UnitRange{T}) where {T} = new{T}(r.start, r.stop)
end
@inline Base.axes(r::WignerRange) = (axes1(r),)
@inline axes1(r::WignerRange) = WignerRange(r.start:r.stop)
if VERSION < v"1.8.2"
    Base.axes1(r::WignerRange) = axes1(r)
end
Base.inds2string(inds::NTuple{2, WignerRange}) =
    string(
        "(", inds[1].start, ":", inds[1].stop, ")",
        "×",
        "(", inds[2].start, ":", inds[2].stop, ")"
    )
Base.firstindex(r::WignerRange) = 1
Base.lastindex(r::WignerRange) = length(r)
function Base.getindex(v::WignerRange, i::Bool)
    throw(ArgumentError("invalid index: $i of type Bool"))
end
@propagate_inbounds function Base.getindex(v::WignerRange{T}, i::Integer) where {T}
    val = convert(T, v.start + (i - oneunit(i)))
    @boundscheck (i>0 && val <= v.stop && val >= v.start) || throw(BoundsError(v, i))
    val
end

function Base.axes(w::AbstractWignerMatrix{IT}) where {IT}
    (WignerRange(m′ₘᵢₙ(w):m′ₘₐₓ(w)), WignerRange(mₘᵢₙ(w):mₘₐₓ(w)))
end


# We don't have to override Base.show; most of its machinery works just fine, except that
# printing the data itself gets screwed up when the indices are Rational.  So we override
# this core part of the printing machinery to just print the parent matrix as usual.  The
# only other thing show really does is add a "summary" line, for which the `axes` and thence
# `inds2string` methods above are used.
Base.print_array(io::IO, w::AbstractWignerMatrix{<:Rational}) = Base.print_array(io, parent(w))

@propagate_inbounds function Base.getindex(w::AbstractWignerMatrix, i::Int)
    @boundscheck if i<1 || i>length(w)
        throw(BoundsError(w, i))
    end
    @inbounds Base.parent(w)[i]
end
@propagate_inbounds function Base.getindex(w::AbstractWignerMatrix{IT}, m′::IT, m::IT) where {IT}
    @boundscheck if m′ ∉ axes(w, 1) || m ∉ axes(w, 2)
        throw(BoundsError(w, (m′, m)))
    end
    @inbounds Base.parent(w)[Int(m′-m′ₘᵢₙ(w))+1, Int(m-mₘᵢₙ(w))+1]
end

@propagate_inbounds function Base.setindex!(w::AbstractWignerMatrix, v, i::Int)
    @boundscheck if i<1 || i>length(w)
        throw(BoundsError(w, i))
    end
    @inbounds Base.parent(w)[i] = v
end
@propagate_inbounds function Base.setindex!(w::AbstractWignerMatrix{IT}, v, m′::IT, m::IT) where {IT}
    @boundscheck if m′ ∉ axes(w, 1) || m ∉ axes(w, 2)
        throw(BoundsError(w, (m′, m)))
    end
    @inbounds Base.parent(w)[Int(m′-m′ₘᵢₙ(w))+1, Int(m-mₘᵢₙ(w))+1] = v
end


function validate_index_ranges(ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT) where
    {IT<:Union{Signed, Rational}}
    if IT <: Rational
        if (
            denominator(ℓₘₐₓ) ≠ 2 ||
            denominator(m′ₘᵢₙ) ≠ 2 || denominator(m′ₘₐₓ) ≠ 2 ||
            denominator(mₘᵢₙ) ≠ 2 || denominator(mₘₐₓ) ≠ 2
        )
            error(
                "For IT=$IT <: Rational, indices must have denominator 2:\n"
                * "\tℓₘₐₓ=$ℓₘₐₓ, m′ₘᵢₙ=$m′ₘᵢₙ, m′ₘₐₓ=$m′ₘₐₓ, mₘᵢₙ=$mₘᵢₙ, mₘₐₓ=$mₘₐₓ.\n"
                * "If you want an integer index type, use IT=<:Integer instead."
            )
        end
    end

    # ℓₘₐₓ must be at least as big as ℓₘᵢₙ(ℓₘₐₓ)
    if ℓₘₐₓ < ℓₘᵢₙ(ℓₘₐₓ)
        error("ℓₘₐₓ=$ℓₘₐₓ must be non-negative.")
    end

    # The m′ and m ranges must be ordered correctly
    if m′ₘₐₓ < m′ₘᵢₙ
        error("m′ₘₐₓ=$m′ₘₐₓ is less than m′ₘᵢₙ=$m′ₘᵢₙ.")
    end
    if mₘₐₓ < mₘᵢₙ
        error("mₘₐₓ=$mₘₐₓ is less than mₘᵢₙ=$mₘᵢₙ.")
    end

    # The m′ and m values must bracket ℓₘᵢₙ
    if m′ₘᵢₙ > ℓₘᵢₙ(ℓₘₐₓ)
        error("m′ₘᵢₙ=$m′ₘᵢₙ is too large for this index type, $IT.")
    end
    if m′ₘₐₓ < ℓₘᵢₙ(ℓₘₐₓ)
        error("m′ₘₐₓ=$m′ₘₐₓ is too small for this index type, $IT.")
    end
    if mₘᵢₙ > ℓₘᵢₙ(ℓₘₐₓ)
        error("mₘᵢₙ=$mₘᵢₙ is too large for this index type, $IT.")
    end
    if mₘₐₓ < ℓₘᵢₙ(ℓₘₐₓ)
        error("mₘₐₓ=$mₘₐₓ is too small for this index type, $IT.")
    end

    # The m′ and m values must be in range for ℓₘₐₓ
    if abs(m′ₘᵢₙ) > ℓₘₐₓ
        error("|m′ₘᵢₙ|=|$m′ₘᵢₙ| is too large for ℓₘₐₓ=$ℓₘₐₓ.")
    end
    if abs(m′ₘₐₓ) > ℓₘₐₓ
        error("|m′ₘₐₓ|=|$m′ₘₐₓ| is too large for ℓₘₐₓ=$ℓₘₐₓ.")
    end
    if abs(mₘᵢₙ) > ℓₘₐₓ
        error("|mₘᵢₙ|=|$mₘᵢₙ| is too large for ℓₘₐₓ=$ℓₘₐₓ.")
    end
    if abs(mₘₐₓ) > ℓₘₐₓ
        error("|mₘₐₓ|=|$mₘₐₓ| is too large for ℓₘₐₓ=$ℓₘₐₓ.")
    end

end

function validate_index_ranges(ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT) where
    {IT<:Union{Signed, Rational}}
    if IT <: Rational
        if (
            denominator(ℓₘₐₓ) ≠ 2 ||
            denominator(m′ₘᵢₙ) ≠ 2 || denominator(m′ₘₐₓ) ≠ 2
        )
            error(
                "For IT=$IT <: Rational, indices must have denominator 2:\n"
                * "\tℓₘₐₓ=$ℓₘₐₓ, m′ₘᵢₙ=$m′ₘᵢₙ, m′ₘₐₓ=$m′ₘₐₓ.\n"
                * "If you want an integer index type, use IT=<:Integer instead."
            )
        end
    end

    # ℓₘₐₓ must be at least as big as ℓₘᵢₙ(ℓₘₐₓ)
    if ℓₘₐₓ < ℓₘᵢₙ(ℓₘₐₓ)
        error("ℓₘₐₓ=$ℓₘₐₓ must be non-negative.")
    end

    # The m′ range must be ordered correctly
    if m′ₘₐₓ < m′ₘᵢₙ
        error("m′ₘₐₓ=$m′ₘₐₓ is less than m′ₘᵢₙ=$m′ₘᵢₙ.")
    end

    # The m′ values must bracket ℓₘᵢₙ
    if m′ₘᵢₙ > ℓₘᵢₙ(ℓₘₐₓ)
        error("m′ₘᵢₙ=$m′ₘᵢₙ is too large for this index type, $IT.")
    end
    if m′ₘₐₓ < ℓₘᵢₙ(ℓₘₐₓ)
        error("m′ₘₐₓ=$m′ₘₐₓ is too small for this index type, $IT.")
    end

    # The m′ values must be in range for ℓₘₐₓ
    if abs(m′ₘᵢₙ) > ℓₘₐₓ
        error("|m′ₘᵢₙ|=|$m′ₘᵢₙ| is too large for ℓₘₐₓ=$ℓₘₐₓ.")
    end
    if abs(m′ₘₐₓ) > ℓₘₐₓ
        error("|m′ₘₐₓ|=|$m′ₘₐₓ| is too large for ℓₘₐₓ=$ℓₘₐₓ.")
    end

end


"""
    WignerMatrix{IT, NT, ST} <: AbstractWignerMatrix{IT, NT, ST}

General concrete subtype of [`AbstractWignerMatrix`](@ref) for Wigner rotation matrices,
which can include D-matrices (when `NT` is complex) or d-matrices (when `NT` is real).

In general, the storage type `ST` can be any `AbstractMatrix{NT}`, but should be 1-based.
That is, the storage should generally be either a `Matrix` or a view.  That matrix will
represent a rectangular array of values representing some or all of the Wigner matrix for a
specific ``ℓ`` value.  The first dimension corresponds to the `m′` index, and the second
dimension corresponds to the `m` index.  The allowed ranges of `m′` and `m` are governed by
the fields `m′ₘₐₓ`, `m′ₘᵢₙ`, `mₘₐₓ`, and `mₘᵢₙ`, which must satisfy
```math
\begin{aligned}
-ℓₘₐₓ &≤ m′ₘᵢₙ ≤ ℓₘᵢₙ ≤ m′ₘₐₓ ≤ ℓₘₐₓ, \\
-ℓₘₐₓ &≤ mₘᵢₙ ≤ ℓₘᵢₙ ≤ mₘₐₓ ≤ ℓₘₐₓ,
\end{aligned}
```
where `ℓₘᵢₙ` is either 0 or 1//2 depending on whether `IT` is an integer or rational type.

"""
struct WignerMatrix{IT, NT, ST} <: AbstractWignerMatrix{IT, NT, ST}
    parent::ST
    ℓ::IT
    m′ₘₐₓ::IT
    m′ₘᵢₙ::IT
    mₘₐₓ::IT
    mₘᵢₙ::IT
end

function WignerMatrix(
    parent::ST, ℓ::IT;
    mp_max::IT=ℓ, mp_min::IT=-ℓ, m_max::IT=ℓ, m_min::IT=-ℓ,
    m′ₘₐₓ::IT=mp_max, m′ₘᵢₙ::IT=mp_min, mₘₐₓ::IT=m_max, mₘᵢₙ::IT=m_min
) where {IT, NT, ST<:AbstractMatrix{NT}}
    validate_index_ranges(ℓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    s₁, s₂ = size(parent)
    if s₁ < Int(m′ₘₐₓ - m′ₘᵢₙ + 1)
        error(
            "The extent of the first dimension in the input data must be at least "
            * "m′ₘₐₓ-m′ₘᵢₙ+1=$m′ₘₐₓ-$m′ₘᵢₙ+1=$(Int(m′ₘₐₓ - m′ₘᵢₙ + 1)); it is $s₁."
        )
    end
    if s₂ < Int(mₘₐₓ - mₘᵢₙ + 1)
        error(
            "The extent of the second dimension in the input data must be at least "
            * "mₘₐₓ-mₘᵢₙ+1=$mₘₐₓ-$mₘᵢₙ+1=$(Int(mₘₐₓ - mₘᵢₙ + 1)); it is $s₂."
        )
    end
    WignerMatrix{IT, NT, ST}(parent, ℓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
end


"""
    WignerDMatrix{IT, RT, ST}

Type alias for [`WignerMatrix`](@ref) with complex number type `Complex{RT}`.
Represents Wigner D-matrices (complex rotation matrices).

# Example
```julia
D = WignerDMatrix(ComplexF64, 2)  # Creates WignerMatrix with NT=ComplexF64
"""
const WignerDMatrix{IT, RT, ST} = WignerMatrix{IT, Complex{RT}, ST} where {IT, RT<:Real, ST<:AbstractMatrix{Complex{RT}}}

const WignerdMatrix{IT, RT, ST} = WignerMatrix{IT, RT, ST} where {IT, RT<:Real, ST<:AbstractMatrix{RT}}

# Constructors for WignerDMatrix (complex)
function WignerDMatrix(parent::ST, ℓ::IT; kwargs...) where {IT, RT<:Real, ST<:AbstractMatrix{Complex{RT}}}
    WignerMatrix(parent, ℓ; kwargs...)
end
function WignerDMatrix(parent::ST, ℓ::IT; kwargs...) where {IT, RT<:Real, ST<:AbstractMatrix{RT}}
    error(
        "WignerDMatrix only supports complex types; the input type is $RT.\n"
        * "Perhaps you meant to use WignerdMatrix?\n"
    )
end
function WignerDMatrix(::Type{Complex{RT}}, ℓ::IT, m′ₘₐₓ::IT=ℓ; kwargs...) where {RT<:Real, IT}
    parent = Matrix{Complex{RT}}(undef, Int(m′ₘₐₓ - (-m′ₘₐₓ) + 1), Int(2ℓ + 1))
    WignerMatrix(parent, ℓ; m′ₘₐₓ=m′ₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ, kwargs...)
end

# Constructors for WignerdMatrix (real)
function WignerdMatrix(parent::ST, ℓ::IT; kwargs...) where {IT, RT<:Real, ST<:AbstractMatrix{RT}}
    WignerMatrix(parent, ℓ; kwargs...)
end
function WignerdMatrix(parent::ST, ℓ::IT; kwargs...) where {IT, RT<:Real, ST<:AbstractMatrix{Complex{RT}}}
    error(
        "WignerdMatrix only supports real types; the input type is Complex{$RT}.\n"
        * "Perhaps you meant to use WignerDMatrix?"
    )
end
function WignerdMatrix(::Type{RT}, ℓ::IT, m′ₘₐₓ::IT=ℓ; kwargs...) where {RT<:Real, IT}
    parent = Matrix{RT}(undef, Int(m′ₘₐₓ - (-m′ₘₐₓ) + 1), Int(2ℓ + 1))
    WignerMatrix(parent, ℓ; m′ₘₐₓ=m′ₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ, kwargs...)
end


@testitem "WignerMatrix" begin
    import SphericalFunctions.Redesign: WignerDMatrix, WignerdMatrix,
        parent, ell, mpmax, mpmin, mmax, mmin, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, ℓₘᵢₙ

    # Check that mixed-up types throw an error
    @test_throws "WignerDMatrix only supports complex types" WignerDMatrix(rand(Float64, 3, 3), 1)
    @test_throws "WignerdMatrix only supports real types" WignerdMatrix(rand(ComplexF64, 3, 3), 1)
    @test_throws "WignerDMatrix only supports complex types" WignerDMatrix(rand(Float64, 2, 2), 1//2)
    @test_throws "WignerdMatrix only supports real types" WignerdMatrix(rand(ComplexF64, 2, 2), 1//2)

    # Check that a negative ℓ value throws an error
    @test_throws "ℓₘₐₓ=-1 must be non-negative." WignerDMatrix(rand(ComplexF64, 3, 3), -1)
    @test_throws "ℓₘₐₓ=-1 must be non-negative." WignerdMatrix(rand(Float64, 3, 3), -1)
    @test_throws "ℓₘₐₓ=-1//2 must be non-negative." WignerDMatrix(rand(ComplexF64, 2, 2), -1//2)
    @test_throws "ℓₘₐₓ=-1//2 must be non-negative." WignerdMatrix(rand(Float64, 2, 2), -1//2)

    # Check that a non-half-integer ℓ value throws an error
    @test_throws "For IT=Rational{Int64} <: Rational, indices must have denominator 2:" WignerDMatrix(rand(ComplexF64, 3, 3), 1//3)
    @test_throws "For IT=Rational{Int64} <: Rational, indices must have denominator 2:" WignerdMatrix(rand(Float64, 3, 3), 1//3)
    @test_throws "For IT=Rational{Int64} <: Rational, indices must have denominator 2:" WignerDMatrix(rand(ComplexF64, 2, 2), 1//3)
    @test_throws "For IT=Rational{Int64} <: Rational, indices must have denominator 2:" WignerdMatrix(rand(Float64, 2, 2), 1//3)
    @test_throws "For IT=Rational{Int64} <: Rational, indices must have denominator 2:" WignerDMatrix(rand(ComplexF64, 3, 3), 2//2)
    @test_throws "For IT=Rational{Int64} <: Rational, indices must have denominator 2:" WignerdMatrix(rand(Float64, 3, 3), 2//2)
    @test_throws "For IT=Rational{Int64} <: Rational, indices must have denominator 2:" WignerDMatrix(rand(ComplexF64, 2, 2), 2//2)
    @test_throws "For IT=Rational{Int64} <: Rational, indices must have denominator 2:" WignerdMatrix(rand(Float64, 2, 2), 2//2)

    #for ℓ ∈ Any[collect(0:8); collect(1//2:15//2)]
    ℓₘₐₓ = 2
    encode(ℓ, m′, m) = (ℓ+ℓₘₐₓ) + (m′+ℓₘₐₓ)*(4ℓₘₐₓ+1) + (m+ℓₘₐₓ)*(4ℓₘₐₓ+1)^2
    for ℓ ∈ Any[collect(0:ℓₘₐₓ); collect(1//2:(ℓₘₐₓ+1//2))]
        mₘ = ℓ

        # These tests are old; the input array can be larger than necessary now.
        # # Check that ℓ < m′ₘₐₓ and ℓ ≠ mₘₐₓ throw errors
        # @test_throws "greater than 0 and less than or equal to 2ℓ+1=" WignerDMatrix(Array{ComplexF64}(undef, Int(2ℓ)+2, Int(2ℓ)+1), ℓ)
        # @test_throws "greater than 0 and less than or equal to 2ℓ+1=" WignerdMatrix(Array{Float64}(undef, Int(2ℓ)+2, Int(2ℓ)+1), ℓ)
        # @test_throws "in the input data must be 2ℓ+1=" WignerDMatrix(Array{ComplexF64}(undef, Int(2ℓ)+1, Int(2ℓ)+2), ℓ)
        # @test_throws "in the input data must be 2ℓ+1=" WignerdMatrix(Array{Float64}(undef, Int(2ℓ)+1, Int(2ℓ)+2), ℓ)

        # # Check that the input is at least as big as needed for the given ℓ
        @test_throws "The extent of the first dimension" WignerDMatrix(Array{ComplexF64}(undef, Int(2ℓ)+0, Int(2ℓ)+1), ℓ)
        @test_throws "The extent of the first dimension" WignerdMatrix(Array{Float64}(undef, Int(2ℓ)+0, Int(2ℓ)+1), ℓ)
        @test_throws "The extent of the second dimension" WignerDMatrix(Array{ComplexF64}(undef, Int(2ℓ)+1, Int(2ℓ)+0), ℓ)
        @test_throws "The extent of the second dimension" WignerdMatrix(Array{Float64}(undef, Int(2ℓ)+1, Int(2ℓ)+0), ℓ)

        # Check that a mismatch between integer/half-integer throws an error
        if ℓ>0 && ℓ isa Int
            @test_throws "The extent of the first dimension" WignerDMatrix(rand(ComplexF64, 2ℓ, 2ℓ+1), ℓ)
            @test_throws "The extent of the first dimension" WignerdMatrix(rand(Float64, 2ℓ, 2ℓ+1), ℓ)
        elseif ℓ isa Rational
            @test_throws "The extent of the first dimension" WignerDMatrix(rand(ComplexF64, Int(2ℓ), Int(2ℓ+1)), ℓ)
            @test_throws "The extent of the first dimension" WignerdMatrix(rand(Float64, Int(2ℓ), Int(2ℓ+1)), ℓ)
        end
        @test_throws "The extent of the second dimension" WignerDMatrix(rand(ComplexF64, Int(2ℓ+1), Int(2ℓ)), ℓ)
        @test_throws "The extent of the second dimension" WignerdMatrix(rand(Float64, Int(2ℓ+1), Int(2ℓ)), ℓ)

        # Check that a data array with a dimension of 0 extent throws an error.
        @test_throws r"The extent of the second dimension.*; it is 0." WignerDMatrix(Array{ComplexF64}(undef, Int(2ℓ)+1, 0), ℓ)
        @test_throws r"The extent of the first dimension.*; it is 0." WignerDMatrix(Array{ComplexF64}(undef, 0, Int(2ℓ)+1), ℓ)
        @test_throws r"The extent of the second dimension.*; it is 0." WignerdMatrix(Array{Float64}(undef, Int(2ℓ)+1, 0), ℓ)
        @test_throws r"The extent of the first dimension.*; it is 0." WignerdMatrix(Array{Float64}(undef, 0, Int(2ℓ)+1), ℓ)

        for m′ₘ ∈ ℓₘᵢₙ(ℓ):ℓ
            # Make a big, dumb array full of the explicit indices.
            data = [
                encode(ℓ, m′, m)
                for m′ ∈ -m′ₘ:m′ₘ, m ∈ -mₘ:mₘ
            ]
            # Check that indexing works as expected.
            for (WignerMatrixType, NT) ∈ ((WignerDMatrix, ComplexF64), (WignerdMatrix, Float64))
                w = WignerMatrixType(NT.(data), ℓ; m′ₘₐₓ=m′ₘ, m′ₘᵢₙ=-m′ₘ, mₘₐₓ=mₘ, mₘᵢₙ=-mₘ)
                @test Base.parent(w) == data
                @test ell(w) == ℓ
                @test mpmax(w) == m′ₘ
                @test mmax(w) == ℓ
                @test mpmin(w) == -mpmax(w)
                @test mmin(w) == -mmax(w)
                for m ∈ -mₘ:mₘ
                    for m′ ∈ -m′ₘ:m′ₘ
                        @test w[m′, m] == encode(ℓ, m′, m)
                    end
                end
            end
        end

        for m′ₘ ∈ ℓₘᵢₙ(ℓ):ℓ
            for WignerMatrixType ∈ (WignerDMatrix, WignerdMatrix)
                data = rand(
                    WignerMatrixType<:WignerDMatrix ? ComplexF64 : Float64,
                    Int(2m′ₘ)+1, Int(2mₘ)+1
                )
                w = WignerMatrixType(data, ℓ; m′ₘₐₓ=m′ₘ, m′ₘᵢₙ=-m′ₘ, mₘₐₓ=mₘ, mₘᵢₙ=-mₘ)

                # Check that the data array is stored correctly.
                @test Base.parent(w) == data
                @test ell(w) == ℓ
                @test m′ₘₐₓ(w) == m′ₘ
                @test mₘₐₓ(w) == ℓ
                @test m′ₘᵢₙ(w) == -m′ₘₐₓ(w)
                @test mₘᵢₙ(w) == -mₘₐₓ(w)

                # The Julia docs say that the `axes` function should
                # > Return a tuple of `AbstractUnitRange{<:Integer}` of valid indices.
                # > The axes should be their own axes, that is `axes.(axes(A),1) ==
                # > axes(A)` should be satisfied.
                # https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array
                @test typeof(axes(w)) <: NTuple{2, AbstractUnitRange}
                @test axes.(axes(w),1) == axes(w)
            end
        end
    end
end


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
data is stored as a 2-dimensional array, with the first dimension indexing different
rotations, and the second dimension storing the vectorized `Hˡ` data as described above.

"""
struct HWedge{IT, RT<:Real, ST<:DenseArray{RT}} <: AbstractWignerMatrix{IT, RT, ST}
    parent::ST
    row_index::FixedSizeVectorDefault{Int}
    ℓ::IT
    m′ₘₐₓ::IT
    m′ₘᵢₙ::IT
end

function HWedge(
    parent::ST, ℓ::IT;
    mp_max::IT=ℓ, mp_min::IT=-ℓ,
    m′ₘₐₓ::IT=mp_max, m′ₘᵢₙ::IT=mp_min
) where {IT<:Union{Integer,Rational}, RT<:Real, ST<:DenseArray{RT}}
    validate_index_ranges(ℓ, m′ₘₐₓ, m′ₘᵢₙ)
    expected_size = Hwedge_size(ℓ, m′ₘₐₓ, m′ₘᵢₙ)
    if ndims(parent) ≠ 2
        error(
            "Input must be a 2-dimensional array; it has $(ndims(parent)) dimensions.\n"
            * "The first dimension is used to vectorize/thread over rotations, and can "
            * "just have extent 1.\n"
            * "The second must have length at least $(expected_size) for the input values "
            * "ℓ=$ℓ, m′ₘₐₓ=$m′ₘₐₓ, and m′ₘᵢₙ=$m′ₘᵢₙ."
        )
    end
    s = size(parent, 2)
    if s < expected_size
        error(
            "The length of the input data must be at least "
            * "(m′ₘₐₓ-m′ₘᵢₙ+1)*(ℓₘₐₓ-ℓₘᵢₙ+1)="
            * "($m′ₘₐₓ-$m′ₘᵢₙ+1)*($ℓₘₐₓ-$(ℓₘᵢₙ(ℓₘₐₓ))+1)="
            * "$(expected_size); it is $s."
        )
    end
    row_index = HWedge_row_index_array(ℓ, m′ₘₐₓ, m′ₘᵢₙ)
    HWedge{IT, RT, ST}(parent, row_index, ℓ, m′ₘₐₓ, m′ₘᵢₙ)
end

function HWedge(
    ::Type{RT}, Nᵣ::Int, ℓ::IT;
    mp_max::IT=ℓ, mp_min::IT=-ℓ,
    m′ₘₐₓ::IT=mp_max, m′ₘᵢₙ::IT=mp_min
) where {IT<:Union{Integer,Rational}, RT<:Real}
    validate_index_ranges(ℓ, m′ₘₐₓ, m′ₘᵢₙ)
    if Nᵣ < 1
        error("Number of rotors Nᵣ=$Nᵣ must be at least 1.")
    end
    parent = FixedSizeArrayDefault{RT}(undef, Nᵣ, HWedge_size(ℓ, m′ₘₐₓ, m′ₘᵢₙ))
    row_index = HWedge_row_index_array(ℓ, m′ₘₐₓ, m′ₘᵢₙ)
    HWedge{IT, RT, ST}(parent, ℓ, m′ₘₐₓ, m′ₘᵢₙ, row_index)
end

mₘₐₓ(w::HWedge{IT}) where {IT} = ℓ(w)
mₘᵢₙ(w::HWedge{IT}) where {IT} = ℓₘᵢₙ(w)

function HWedge_row_index_array(ℓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT) where {IT}
    m′range = m′ₘᵢₙ:m′ₘₐₓ
    row_index = FixedSizeVectorDefault{Int}(undef, length(m′range))
    index = 1
    for (i, m′) ∈ enumerate(m′range)
        row_index[i] = index
        index += Int(ℓ - abs(m′)) + 1
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

function HWedge_storage(::Type{RT}, Nᵣ::Int, ℓₘₐₓ::IT, m′ₘₐₓ::IT=ℓₘₐₓ, m′ₘᵢₙ::IT=ℓₘₐₓ) where {RT<:Real, IT}
    FixedSizeArrayDefault{RT}(
        undef,
        Nᵣ,
        HWedge_size(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ)
    )
end

"""

An Hˡ wedge will store elements in a vector as if it were the following matrix:

    [
        H[ℓ, m′, m]
        for m′ ∈ max(-ℓ, m′ₘᵢₙ):min(ℓ, m′ₘₐₓ)
        for m ∈ abs(m′):ℓ
    ]

Here, m′ₘᵢₙ is a negative number and m′ₘₐₓ is a positive number.  Note that for HWedge, we
currently impose mₘₐₓ = ℓ and mₘᵢₙ = ℓₘᵢₙ(IT), because these are all needed for the
recurrence relations.

This function returns the linear index into that vector that belongs to the first element
with the given `m′` value (and therefore `m=abs(m′)`).  The formula for that index involves
an `if` statement to account for the varying number of `m` values for each `m′` value.
Nonetheless, it can be computed in closed form (i.e., without an explicit sum or loop).


"""



function row_index(w::HWedge{IT}, m′::IT) where {IT}
    let ℓ = ℓ(w), m′ₘᵢₙ = m′ₘᵢₙ(w), ℓₘᵢₙ = ℓₘᵢₙ(IT)
        (
            Int(ℓₘᵢₙ - m′ₘᵢₙ) * Int(2ℓ + m′ₘᵢₙ + ℓₘᵢₙ + 1)
            -
            Int(ℓₘᵢₙ - m′) * Int(2ℓ - abs(m′ + ℓₘᵢₙ - 1) + 2)
        ) ÷ 2 + 1

        # i = if m′<1
        #     Int(m′ - m′ₘᵢₙ) * Int(2ℓ + m′ + m′ₘᵢₙ + 1) ÷ 2  # size of wedge to the left of m'
        # else
        #     (
        #         # size of entire left half of wedge
        #         Int(ℓₘᵢₙ - m′ₘᵢₙ) * Int(2ℓ + ℓₘᵢₙ + m′ₘᵢₙ + 1)
        #         +
        #         # size of right half of wedge to the left of m'
        #         Int(m′ - ℓₘᵢₙ) * Int(2ℓ - ℓₘᵢₙ - m′ + 3)
        #     ) ÷ 2
        # end
        # i + 1
    end
end


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
