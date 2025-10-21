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
function Base.axes(w::AbstractWignerMatrix{IT}) where {IT}
    ((m′ₘᵢₙ(w):m′ₘₐₓ(w)), (mₘᵢₙ(w):mₘₐₓ(w)))
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
                * "\tℓₘₐₓ=$ℓₘₐₓ, m′ₘᵢₙ=$m′ₘᵢₙ, m′ₘₐₓ=$m′ₘₐₓ, mₘᵢₙ=$mₘᵢₙ, mₘₐₓ=$mₘₐₓ."
            )
        end
    end

    if ℓₘₐₓ < ℓₘᵢₙ(ℓₘₐₓ)
        error("ℓₘₐₓ=$ℓₘₐₓ must be non-negative.")
    end

    # The m′ and m values must bracket ℓₘᵢₙ
    if m′ₘᵢₙ > ℓₘᵢₙ(ℓₘₐₓ)
        error("m′ₘᵢₙ=$m′ₘᵢₙ is too large for this index type.")
    end
    if m′ₘₐₓ < ℓₘᵢₙ(ℓₘₐₓ)
        error("m′ₘₐₓ=$m′ₘₐₓ is too small for this index type.")
    end
    if mₘᵢₙ > ℓₘᵢₙ(ℓₘₐₓ)
        error("mₘᵢₙ=$mₘᵢₙ is too large for this index type.")
    end
    if mₘₐₓ < ℓₘᵢₙ(ℓₘₐₓ)
        error("mₘₐₓ=$mₘₐₓ is too small for this index type.")
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
            * "m′ₘₐₓ-m′ₘᵢₙ+1=$(Int(m′ₘₐₓ - m′ₘᵢₙ + 1)); it is $s₁."
        )
    end
    if s₂ < Int(mₘₐₓ - mₘᵢₙ + 1)
        error(
            "The extent of the second dimension in the input data must be at least "
            * "mₘₐₓ-mₘᵢₙ+1=$(Int(mₘₐₓ - mₘᵢₙ + 1)); it is $s₂."
        )
    end
    WignerMatrix{IT, NT, ST}(parent, ℓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
end

end
function WignerDMatrix(::Type{NT}, ℓ::IT, m′::IT=ℓ) where {NT, IT}
    if complex(NT) ≢ NT
        error(
            "`WignerDMatrix` only supports complex types; the input type is $NT.\n"
            * "Perhaps you meant to use `WignerdMatrix`?"
        )
    end
    WignerDMatrix{IT, NT, Matrix{NT}}(Matrix{NT}(undef, Int(2m′)+1, Int(2ℓ)+1), ℓ)
end



"""
    WignerdMatrix{IT, NT, ST} <: AbstractWignerMatrix{IT, NT, ST}

Specialized subtype of [`AbstractWignerMatrix`](@ref) for d-matrices, which are real matrices.
"""
struct WignerdMatrix{IT, NT, ST} <: AbstractWignerMatrix{IT, NT, ST}
    parent::ST
    ℓ::IT
    m′ₘₐₓ::IT
    function WignerdMatrix{IT, NT, ST}(parent::ST, ℓ::IT) where {IT, NT, ST}
        # We want to secretly allow NTuple{3, IT} for testing purposes, so we can't just use
        # a restriction on NT in the type declaration.
        if !(NT <: NTuple{3, IT}) && real(NT) ≢ NT
            error(
                "WignerdMatrix only supports real types; the input type is $NT.\n"
                * "Perhaps you meant to use WignerDMatrix?"
            )
        end
        if ℓ < 0 || (IT <: Rational && denominator(ℓ) ≠ 2)
            error(
                "ℓ=$ℓ should be non-negative integer or half-integer.  In particular,\n"
                * "if ℓ is an integer its type must be <:Integer, not <:Rational."
            )
        end
        s₁, s₂ = size(parent)
        if s₂ ≠ Int(2ℓ + 1)
            error(
                "The extent of the second dimension in the input data must be "
                * "2ℓ+1=$(Int(2ℓ+1)); it is $s₂."
            )
        end
        if s₁ == 0 || s₁ > s₂
            error(
                "The extent of the first dimension in the input data must be greater than 0"
                * " and less than or equal to 2ℓ+1=$(Int(2ℓ+1)); it is $s₁."
            )
        end
        if IT <: Rational
            if isodd(s₁)
                error(
                    "ℓ=$ℓ is a half-integer, but the extent of the first dimension in the "
                    * "input data ($s₁) corresponds to whole-integer values of m′."
                )
            end
        else
            if iseven(s₁)
                error(
                    "ℓ=$ℓ is an integer, but the extent of the first dimension in the "
                    * "input data ($s₁) corresponds to half-integer values of m′."
                )
            end
        end
        m′ₘₐₓ = IT((s₁ - 1) // 2)
        new(parent, ℓ, m′ₘₐₓ)
    end
end

"""
    WignerdMatrix(parent, ℓ)

Construct a `WignerdMatrix` object from the given parent matrix and ``ℓ`` value.  Note that
the type of `ℓ` *must* be either `Integer` or `Rational`.  If it is `Rational`, the
denominator *must* be 2; if it is 1, you must convert to an `Int` first.  Also, the parent
matrix must have the correct size: the first dimension must be greater than 0 and less than
or equal to `2ℓ+1`, and the second dimension must be equal to `2ℓ+1`.
"""
function WignerdMatrix(parent::ST, ℓ::IT) where {IT, ST}
    WignerdMatrix{IT, eltype(ST), ST}(parent, ℓ)
end
function WignerdMatrix(::Type{NT}, ℓ::IT, m′::IT=ℓ) where {NT, IT}
    if real(NT) ≢ NT
        error(
            "`WignerdMatrix` only supports real types; the input type is $NT.\n"
            * "Perhaps you meant to use `WignerDMatrix`?"
        )
    end
    WignerdMatrix{IT, NT, Matrix{NT}}(Matrix{NT}(undef, Int(2m′)+1, Int(2ℓ)+1), ℓ)
end

"""
    Hˡrow{IT, NT, ST}

Specialized subtype of [`AbstractWignerMatrix`](@ref) intended to store one row of the ``H`` matrix
— usually the ``H^{\ell-1}_{0,m}`` or ``H^{\ell+1}_{0,m}`` components needed during the
recurrence relations.
"""
struct Hˡrow{IT, NT, ST} <: AbstractWignerMatrix{IT, NT, ST}
    parent::ST
    ℓ::IT
    m′ₘₐₓ::IT
end
function Hˡrow(parent::ST, ℓ::IT, m′::IT) where {IT, ST}
    length_m′ = 1
    length_m = Int(ℓ - ℓₘᵢₙ(ℓ)) + 1
    if size(parent,1) < length_m′ || size(parent,2) < length_m
        error(
            "The input `parent` matrix for ℓ=$ℓ must have size at least "
            * "($length_m′,$length_m); it has size $(size(parent))."
        )
    end
    Hˡrow{IT, eltype(ST), ST}(parent, ℓ, m′)
end
function Hˡrow(::Type{NT}, ℓ::IT, m′::IT) where {NT, IT}
    if real(NT) ≢ NT
        error("`Hˡrow` only supports real types; the input type is $NT.")
    end
    length_m′ = 1
    length_m = Int(ℓ - ℓₘᵢₙ(ℓ)) + 1
    Hˡrow{IT, NT, Matrix{NT}}(Matrix{NT}(undef, length_m′, length_m), ℓ, m′)
end

m′ₘᵢₙ(w::Hˡrow) = m′ₘₐₓ(w)
mₘₐₓ(w::Hˡrow) = ℓ(w)
mₘᵢₙ(w::Hˡrow) = ℓₘᵢₙ(w)


@testitem "WignerMatrix" begin
    import SphericalFunctions.Redesign: WignerDMatrix, WignerdMatrix,
        parent, ell, mpmax, mpmin, mmax, mmin, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, ℓₘᵢₙ

    # Check that mixed-up types throw an error
    @test_throws "WignerDMatrix only supports complex types" WignerDMatrix(rand(Float64, 3, 3), 1)
    @test_throws "WignerdMatrix only supports real types" WignerdMatrix(rand(ComplexF64, 3, 3), 1)
    @test_throws "WignerDMatrix only supports complex types" WignerDMatrix(rand(Float64, 2, 2), 1//2)
    @test_throws "WignerdMatrix only supports real types" WignerdMatrix(rand(ComplexF64, 2, 2), 1//2)

    # Check that a negative ℓ value throws an error
    @test_throws "should be non-negative integer or half-integer." WignerDMatrix(rand(ComplexF64, 3, 3), -1)
    @test_throws "should be non-negative integer or half-integer." WignerdMatrix(rand(Float64, 3, 3), -1)
    @test_throws "should be non-negative integer or half-integer." WignerDMatrix(rand(ComplexF64, 2, 2), -1//2)
    @test_throws "should be non-negative integer or half-integer." WignerdMatrix(rand(Float64, 2, 2), -1//2)

    # Check that a non-half-integer ℓ value throws an error
    @test_throws "should be non-negative integer or half-integer." WignerDMatrix(rand(ComplexF64, 3, 3), 1//3)
    @test_throws "should be non-negative integer or half-integer." WignerdMatrix(rand(Float64, 3, 3), 1//3)
    @test_throws "should be non-negative integer or half-integer." WignerDMatrix(rand(ComplexF64, 2, 2), 1//3)
    @test_throws "should be non-negative integer or half-integer." WignerdMatrix(rand(Float64, 2, 2), 1//3)
    @test_throws "should be non-negative integer or half-integer." WignerDMatrix(rand(ComplexF64, 3, 3), 2//2)
    @test_throws "should be non-negative integer or half-integer." WignerdMatrix(rand(Float64, 3, 3), 2//2)
    @test_throws "should be non-negative integer or half-integer." WignerDMatrix(rand(ComplexF64, 2, 2), 2//2)
    @test_throws "should be non-negative integer or half-integer." WignerdMatrix(rand(Float64, 2, 2), 2//2)

    #for ℓ ∈ Any[collect(0:8); collect(1//2:15//2)]
    for ℓ ∈ Any[collect(0:2); collect(1//2:3//2)]
        mₘ = ℓ

        # Check that ℓ < m′ₘₐₓ and ℓ ≠ mₘₐₓ throw errors
        @test_throws "greater than 0 and less than or equal to 2ℓ+1=" WignerDMatrix(Array{ComplexF64}(undef, Int(2ℓ)+2, Int(2ℓ)+1), ℓ)
        @test_throws "greater than 0 and less than or equal to 2ℓ+1=" WignerdMatrix(Array{Float64}(undef, Int(2ℓ)+2, Int(2ℓ)+1), ℓ)
        @test_throws "in the input data must be 2ℓ+1=" WignerDMatrix(Array{ComplexF64}(undef, Int(2ℓ)+1, Int(2ℓ)+2), ℓ)
        @test_throws "in the input data must be 2ℓ+1=" WignerdMatrix(Array{Float64}(undef, Int(2ℓ)+1, Int(2ℓ)+2), ℓ)
        @test_throws "in the input data must be 2ℓ+1=" WignerDMatrix(Array{ComplexF64}(undef, Int(2ℓ)+1, Int(2ℓ)+0), ℓ)
        @test_throws "in the input data must be 2ℓ+1=" WignerdMatrix(Array{Float64}(undef, Int(2ℓ)+1, Int(2ℓ)+0), ℓ)

        # Check that a mismatch between integer/half-integer throws an error
        if ℓ>0 && ℓ isa Int
            @test_throws "is an integer, but the extent of the first dimension" WignerDMatrix(rand(ComplexF64, 2ℓ, 2ℓ+1), ℓ)
            @test_throws "is an integer, but the extent of the first dimension" WignerdMatrix(rand(Float64, 2ℓ, 2ℓ+1), ℓ)
        elseif ℓ isa Rational
            @test_throws "is a half-integer, but the extent of the first dimension" WignerDMatrix(rand(ComplexF64, Int(2ℓ), Int(2ℓ+1)), ℓ)
            @test_throws "is a half-integer, but the extent of the first dimension" WignerdMatrix(rand(Float64, Int(2ℓ), Int(2ℓ+1)), ℓ)
        end
        @test_throws "in the input data must be 2ℓ+1=" WignerDMatrix(rand(ComplexF64, Int(2ℓ+1), Int(2ℓ)), ℓ)
        @test_throws "in the input data must be 2ℓ+1=" WignerdMatrix(rand(Float64, Int(2ℓ+1), Int(2ℓ)), ℓ)

        # Check that a data array with a dimension of 0 extent throws an error.
        @test_throws "in the input data must be 2ℓ+1=" WignerDMatrix(Array{ComplexF64}(undef, Int(2ℓ)+1, 0), ℓ)
        @test_throws "greater than 0 and less than or equal to 2ℓ+1=" WignerDMatrix(Array{ComplexF64}(undef, 0, Int(2ℓ)+1), ℓ)
        @test_throws "in the input data must be 2ℓ+1=" WignerdMatrix(Array{Float64}(undef, Int(2ℓ)+1, 0), ℓ)
        @test_throws "greater than 0 and less than or equal to 2ℓ+1=" WignerdMatrix(Array{Float64}(undef, 0, Int(2ℓ)+1), ℓ)

        for m′ₘ ∈ ℓₘᵢₙ(ℓ):ℓ
            # Make a big, dumb array full of the explicit indices.
            data = [
                (ℓ, m′, m)
                for m′ ∈ -m′ₘ:m′ₘ, m ∈ -mₘ:mₘ
            ]
            # Check that indexing works as expected.
            for WignerMatrixType ∈ (WignerDMatrix, WignerdMatrix)
                w = WignerMatrixType(data, ℓ)
                @test Base.parent(w) == data
                @test ell(w) == ℓ
                @test mpmax(w) == m′ₘ
                @test mmax(w) == ℓ
                @test mpmin(w) == -mpmax(w)
                @test mmin(w) == -mmax(w)
                for m ∈ -mₘ:mₘ
                    for m′ ∈ -m′ₘ:m′ₘ
                        @test w[m′, m] == (ℓ, m′, m)
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
                w = WignerMatrixType(data, ℓ)

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
