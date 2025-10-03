import Base: @propagate_inbounds

"""
    WignerMatrix{NT, IT}

Abstract base type for Wigner rotation‐matrix objects of a specific ``ℓ`` value.
- `NT` is the number type (e.g., `ComplexF64` for D-matrices or `Float64` for d-matrices).
- `IT` is the index type (an `Integer` or half‐integer `Rational`), governing the allowed
  ranges of `m′` and `m`.

The basic concrete subtypes (`WignerDMatrix`, `WignerdMatrix`) store their data in a
`Matrix{NT}` and implement the usual `size`, `getindex` and `setindex!` so that one can use
`w[m′,m]`.  Specifically, these indices can be negative or positive, and must obey `abs(m′)
≤ m′ₘₐₓ` and `abs(m) ≤ ℓ`.

# Methods

Methods defined for `WignerMatrix` objects include:
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

Any new subtypes of `WignerMatrix` should inherit from this type and re-implement any of the
methods mentioned above that are not appropriate for the new type.  Specifically, the
default implementations assume that subtypes store the fields
- `parent::Matrix{NT}`: the underlying data array.
- `ℓ::IT`: the value of ``ℓ``.
- `m′ₘₐₓ::IT`: the maximum value of ``m′``.

For example, if the parent Matrix is not stored as the `parent` field, then the `parent(w)`
method should be re-implemented to return the correct parent object.  The `getindex` and
`setindex!`
"""
abstract type WignerMatrix{NT, IT} <: AbstractMatrix{NT} end

### General methods for all WignerMatrix types

Base.parent(w::WignerMatrix{NT, IT}) where {NT, IT} = w.parent
ℓ(w::WignerMatrix{NT, IT}) where {NT, IT} = w.ℓ
m′ₘₐₓ(w::WignerMatrix{NT, IT}) where {NT, IT} = w.m′ₘₐₓ
mₘₐₓ(w::WignerMatrix{NT, IT}) where {NT, IT} = ℓ(w)

ℓₘᵢₙ(::IT) where {IT} = ℓₘᵢₙ(IT)
ℓₘᵢₙ(::Type{IT}) where {IT<:Integer} = zero(IT)
ℓₘᵢₙ(::Type{IT}) where {IT<:Rational} = IT(1//2)
ℓₘᵢₙ(::WignerMatrix{NT, IT}) where {NT, IT} = ℓₘᵢₙ(IT)

const ell = ℓ
const mpmax = m′ₘₐₓ
const mmax = mₘₐₓ
const ellmin = ℓₘᵢₙ

isrational(::WignerMatrix{NT, IT}) where {NT, IT<:Integer} = false
isrational(::WignerMatrix{NT, IT}) where {NT, IT<:Rational} = true

Base.eltype(::WignerMatrix{NT, IT}) where {NT, IT} = NT
Base.size(w::WignerMatrix{NT, IT}) where {NT, IT} = size(parent(w))
Base.length(w::WignerMatrix{NT, IT}) where {NT, IT} = length(parent(w))

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
    string(inds[1].start, ":", inds[1].stop, "×", inds[2].start, ":", inds[2].stop)

function Base.axes(w::WignerMatrix{NT, IT}) where {NT, IT}
    (WignerRange(-m′ₘₐₓ(w):m′ₘₐₓ(w)), WignerRange(-mₘₐₓ(w):mₘₐₓ(w)))
end

# We don't have to override Base.show; most of its machinery works just fine, except that
# printing the data itself gets screwed up when the indices are Rational.  So we override
# this core part of the printing machinery to just print the parent matrix as usual.  The
# only other thing show really does is add a "summary" line, for which the only
Base.print_array(io::IO, w::WignerMatrix{NT, IT}) where {NT, IT<:Rational} =
    Base.print_array(io, parent(w))

@propagate_inbounds function Base.getindex(w::WignerMatrix{NT, IT}, i::Int) where {NT, IT}
    @boundscheck begin
        if i<1 || i>length(w)
            throw(BoundsError(
                "i=$i out of bounds for WignerMatrix with length=$(length(w))."
            ))
        end
    end
    Base.parent(w)[i]
end
@propagate_inbounds function Base.getindex(w::WignerMatrix{NT, IT}, m′::IT, m::IT) where {NT, IT}
    @boundscheck begin
        if abs(m′) > m′ₘₐₓ(w)
            throw(BoundsError(
                "m′=$m′ out of bounds for WignerMatrix with m′ₘₐₓ=$(m′ₘₐₓ(w))."
            ))
        end
        if abs(m) > ℓ(w)
            throw(BoundsError(
                "m=$m out of bounds for WignerMatrix with ℓ=$(ℓ(w))."
            ))
        end
    end
    @inbounds Base.parent(w)[Int(m′+m′ₘₐₓ(w))+1, Int(m+mₘₐₓ(w))+1]
end

@propagate_inbounds function Base.setindex!(w::WignerMatrix{NT, IT}, v::NT, i::Int) where {NT, IT}
    @boundscheck begin
        if i<1 || i>length(w)
            throw(BoundsError(
                "i=$i out of bounds for WignerMatrix with length=$(length(w))."
            ))
        end
    end
    Base.parent(w)[i] = v
end
@propagate_inbounds function Base.setindex!(w::WignerMatrix{NT, IT}, v::NT, m′::IT, m::IT) where {NT, IT}
    @boundscheck begin
        if abs(m′) > m′ₘₐₓ(w)
            throw(BoundsError(
                "m′=$m′ out of bounds for WignerMatrix with m′ₘₐₓ=$(m′ₘₐₓ(w))."
            ))
        end
        if abs(m) > ℓ(w)
            throw(BoundsError(
                "m=$m out of bounds for WignerMatrix with ℓ=$(ℓ(w))."
            ))
        end
    end
    Base.parent(w)[Int(m′+m′ₘₐₓ(w))+1, Int(m+mₘₐₓ(w))+1] = v
end


### Specialize to D and d matrices

"""
    WignerDMatrix{NT, IT}

Specialized subtype of [`WignerMatrix`](@ref) for D-matrices, which are complex matrices.
"""
struct WignerDMatrix{NT, IT} <: WignerMatrix{NT, IT}
    parent::Matrix{NT}
    ℓ::IT
    m′ₘₐₓ::IT
    function WignerDMatrix{NT, IT}(parent::Matrix{NT}, ℓ::IT) where {NT, IT<:Union{Integer, Rational}}
        # We want to secretly allow NTuple{3, IT} for testing purposes, so we can't just use
        # a restriction on NT in the type declaration.
        if !(NT <: NTuple{3, IT}) && complex(NT) ≢ NT
            throw(ErrorException(
                "WignerDMatrix only supports complex types; the input type is $NT.\n"
                * "Perhaps you meant to use WignerdMatrix?"
            ))
        end
        if ℓ < 0 || (IT <: Rational && denominator(ℓ) ≠ 2)
            throw(ErrorException(
                "ℓ=$ℓ should be non-negative integer or half-integer.  In particular,\n"
                * "if ℓ is an integer its type must be <:Integer, not <:Rational."
            ))
        end
        s₁, s₂ = size(parent)
        if s₂ ≠ Int(2ℓ + 1)
            throw(ErrorException(
                "The extent of the second dimension in the input data must be "
                * "2ℓ+1=$(Int(2ℓ+1)); it is $s₂."
            ))
        end
        if s₁ == 0 || s₁ > s₂
            throw(ErrorException(
                "The extent of the first dimension in the input data must be greater than 0"
                * " and less than or equal to 2ℓ+1=$(Int(2ℓ+1)); it is $s₁."
            ))
        end
        if IT <: Rational
            if isodd(s₁)
                throw(ErrorException(
                    "ℓ=$ℓ is a half-integer, but the extent of the first dimension in the "
                    * "input data ($s₁) corresponds to whole-integer values of m′."
                ))
            end
        else
            if iseven(s₁)
                throw(ErrorException(
                    "ℓ=$ℓ is an integer, but the extent of the first dimension in the "
                    * "input data ($s₁) corresponds to half-integer values of m′."
                ))
            end
        end
        m′ₘₐₓ = IT((s₁ - 1) // 2)
        new(parent, ℓ, m′ₘₐₓ)
    end
end

"""
    WignerDMatrix(parent, ℓ)

Construct a `WignerDMatrix` object from the given parent matrix and ``ℓ`` value.  Note that
the type of `ℓ` *must* be either `Integer` or `Rational`.  If it is `Rational`, the
denominator *must* be 2; if it is 1, you must convert to an `Int` first.  Also, the parent
matrix must have the correct size: the first dimension must be greater than 0 and less than
or equal to `2ℓ+1`, and the second dimension must be equal to `2ℓ+1`.
"""
function WignerDMatrix(parent::Matrix{NT}, ℓ::IT) where {NT, IT}
    WignerDMatrix{NT, IT}(parent, ℓ)
end
function WignerDMatrix(::Type{NT}, ℓ::IT, m′::IT=ℓ) where {NT, IT}
    if complex(NT) ≢ NT
        throw(ErrorException(
            "WignerDMatrix only supports complex types; the input type is $NT.\n"
            * "Perhaps you meant to use WignerdMatrix?"
        ))
    end
    WignerDMatrix{NT, IT}(Matrix{NT}(undef, Int(2m′)+1, Int(2ℓ)+1), ℓ)
end



"""
    WignerdMatrix{NT, IT}

Specialized subtype of [`WignerMatrix`](@ref) for d-matrices, which are real matrices.
"""
struct WignerdMatrix{NT, IT} <: WignerMatrix{NT, IT}
    parent::Matrix{NT}
    ℓ::IT
    m′ₘₐₓ::IT
    function WignerdMatrix{NT, IT}(parent::Matrix{NT}, ℓ::IT) where {NT, IT<:Union{Integer, Rational}}
        # We want to secretly allow NTuple{3, IT} for testing purposes, so we can't just use
        # a restriction on NT in the type declaration.
        if !(NT <: NTuple{3, IT}) && real(NT) ≢ NT
            throw(ErrorException(
                "WignerdMatrix only supports real types; the input type is $NT.\n"
                * "Perhaps you meant to use WignerDMatrix?"
            ))
        end
        if ℓ < 0 || (IT <: Rational && denominator(ℓ) ≠ 2)
            throw(ErrorException(
                "ℓ=$ℓ should be non-negative integer or half-integer.  In particular,\n"
                * "if ℓ is an integer its type must be <:Integer, not <:Rational."
            ))
        end
        s₁, s₂ = size(parent)
        if s₂ ≠ Int(2ℓ + 1)
            throw(ErrorException(
                "The extent of the second dimension in the input data must be "
                * "2ℓ+1=$(Int(2ℓ+1)); it is $s₂."
            ))
        end
        if s₁ == 0 || s₁ > s₂
            throw(ErrorException(
                "The extent of the first dimension in the input data must be greater than 0"
                * " and less than or equal to 2ℓ+1=$(Int(2ℓ+1)); it is $s₁."
            ))
        end
        if IT <: Rational
            if isodd(s₁)
                throw(ErrorException(
                    "ℓ=$ℓ is a half-integer, but the extent of the first dimension in the "
                    * "input data ($s₁) corresponds to whole-integer values of m′."
                ))
            end
        else
            if iseven(s₁)
                throw(ErrorException(
                    "ℓ=$ℓ is an integer, but the extent of the first dimension in the "
                    * "input data ($s₁) corresponds to half-integer values of m′."
                ))
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
function WignerdMatrix(parent::Matrix{NT}, ℓ::IT) where {NT, IT}
    WignerdMatrix{NT, IT}(parent, ℓ)
end
function WignerdMatrix(::Type{NT}, ℓ::IT, m′::IT=ℓ) where {NT, IT}
    if real(NT) ≢ NT
        throw(ErrorException(
            "WignerdMatrix only supports real types; the input type is $NT.\n"
            * "Perhaps you meant to use WignerDMatrix?"
        ))
    end
    WignerdMatrix{NT, IT}(Matrix{NT}(undef, Int(2m′)+1, Int(2ℓ)+1), ℓ)
end


@testitem "WignerMatrix" begin
    import SphericalFunctions.Redesign: WignerDMatrix, WignerdMatrix,
        parent, ell, mpmax, mmax, m′ₘₐₓ, mₘₐₓ

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
