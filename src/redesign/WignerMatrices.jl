"""
    AbstractWignerMatrices{NT, IT, MT}

A container for a series of Wigner matrices (
- `NT` is the number type (e.g., `ComplexF64` for D-matrices or `Float64` for d-matrices).
- `IT` is the index type (an `Integer` or half‐integer `Rational`), governing the allowed
  ranges of `m′` and `m` in each matrix.
- `MT` is the type of the matrices.

"""
abstract type AbstractWignerMatrices{NT, IT, MT} <: AbstractVector{MT} end



struct WignerDMatrices{NT, IT, MT} <: AbstractWignerMatrices{NT, IT, MT}
    D::Vector{MT}
    ℓₘₐₓ::IT
    m′ₘₐₓ::IT
    function WignerDMatrices(
        D::Vector{WignerDMatrix{NT, IT}}
    ) where {NT, IT}
        new{NT, IT, MT}(D, ℓₘₐₓ, m′ₘₐₓ)
    end
end


# abstract type AbstractDMatrices end


# """
#     WignerDMatrices{NT, IT}

# A data structure to hold the Wigner D-matrices for a range of `ℓ` values (stored in a
# `Vector{NT}`) up to and including some `ℓₘₐₓ`, `m′ₘₐₓ`, and `mₘₐₓ` (which all have type
# `IT`).

# Indexing this object with an integer `ℓ` returns an `OffsetArray` of a view of the relevant
# part of the data vector corresponding to the `ℓ` matrix.
# """
# struct WignerDMatrices{NT, IT} <: AbstractDMatrices
#     data::Vector{NT}
#     ℓₘₐₓ::IT
#     m′ₘₐₓ::IT
# end

# data(D::WignerDMatrices) = D.data
# ℓₘᵢₙ(D::WignerDMatrices{NT, IT}) where {NT, IT<:Integer} = zero(IT)
# ℓₘᵢₙ(D::WignerDMatrices{NT, IT}) where {NT, IT<:Rational} = IT(1//2)
# ℓₘₐₓ(D::WignerDMatrices) = D.ℓₘₐₓ
# m′ₘₐₓ(D::WignerDMatrices) = D.m′ₘₐₓ
# mₘₐₓ(D::WignerDMatrices) = D.mₘₐₓ
# m′ₘₐₓ(D::WignerDMatrices{NT, IT}, ℓ::IT) where {NT, IT} = min(m′ₘₐₓ(D), ℓ)
# mₘₐₓ(D::WignerDMatrices{NT, IT}, ℓ::IT) where {NT, IT} = min(mₘₐₓ(D), ℓ)

# Base.eltype(D::WignerDMatrices) = eltype(data(D))

# isrational(D::WignerDMatrices{NT, IT}) where {NT, IT<:Integer} = false
# isrational(D::WignerDMatrices{NT, IT}) where {NT, IT<:Rational} = true


# """
#     WignerDsize(ℓₘₐₓ, m′ₘₐₓ, mₘₐₓ)

# Return the total size of the data stored in a `WignerDMatrices` object with the given sizes,
# ranging over all matrices for all ℓ values.
# """
# function WignerDsize(ℓₘₐₓ, m′ₘₐₓ, mₘₐₓ)::Int
#     m₁, m₂ = m′ₘₐₓ, mₘₐₓ
#     if m₁ > m₂
#         m₁, m₂ = m₂, m₁
#     end

#     if ℓₘₐₓ ≤ m₁
#         (2ℓₘₐₓ + 1)*(2ℓₘₐₓ + 2)*(2ℓₘₐₓ + 3) ÷ 6
#     elseif ℓₘₐₓ ≤ m₂
#         (
#             (2m₁ + 1)*(2m₁ + 2)*(2m₁ + 3) ÷ 6
#             + (ℓₘₐₓ - m₁)*(2m₁ + 1)*(ℓₘₐₓ + m₁ + 2)
#         )
#     else
#         (
#             (2m₁ + 1)*(2m₁ + 2)*(2m₁ + 3) ÷ 6
#             + (m₂ - m₁)*(2m₁ + 1)*(m₂ + m₁ + 2)
#             + (2m₁ + 1)*(2m₂ + 1)*(ℓₘₐₓ - m₂)
#         )
#     end
# end


# @testsnippet WignerDUtilities begin
#     function indices_vector(ℓₘₐₓ, m′ₘₐₓ, mₘₐₓ)
#         data = Vector{Tuple{Int64, Int64, Int64}}(undef, sum((2ℓ+1)^2 for ℓ ∈ 0:ℓₘₐₓ))
#         i=1
#         for ℓ ∈ 0:ℓₘₐₓ
#             for m ∈ -min(ℓ, mₘₐₓ):min(ℓ, mₘₐₓ)
#                 for m′ ∈ -min(ℓ, m′ₘₐₓ):min(ℓ, m′ₘₐₓ)
#                     data[i] = (ℓ, m′, m)
#                     i += 1
#                 end
#             end
#         end
#         data
#     end
# end


# @testitem "Test WignerDsize" setup=[WignerDUtilities] begin
#     import SphericalFunctions.Redesign: WignerDsize

#     for ℓₘₐₓ ∈ 0:8
#         @test WignerDsize(ℓₘₐₓ, 0, 0) == ℓₘₐₓ + 1
#         @test WignerDsize(ℓₘₐₓ, 1, 0) == 3ℓₘₐₓ + 1
#         @test WignerDsize(ℓₘₐₓ, 0, 1) == 3ℓₘₐₓ + 1
#         @test WignerDsize(ℓₘₐₓ, 1, 1) == (3^2)ℓₘₐₓ + 1
#         @test WignerDsize(ℓₘₐₓ, 2, 0) == max(1, 5ℓₘₐₓ - 1)
#         @test WignerDsize(ℓₘₐₓ, 0, 2) == max(1, 5ℓₘₐₓ - 1)
#         @test WignerDsize(ℓₘₐₓ, 2, 1) == max(1, 15ℓₘₐₓ - 5)
#         @test WignerDsize(ℓₘₐₓ, 1, 2) == max(1, 15ℓₘₐₓ - 5)
#         @test WignerDsize(ℓₘₐₓ, 2, 2) == max(1, (5^2)ℓₘₐₓ - 15)
#         @test WignerDsize(ℓₘₐₓ, ℓₘₐₓ, ℓₘₐₓ) == sum((2ℓ+1)^2 for ℓ ∈ 0:ℓₘₐₓ)

#         for mₘₐₓ ∈ 0:ℓₘₐₓ
#             for m′ₘₐₓ ∈ 0:ℓₘₐₓ
#                 @test WignerDsize(ℓₘₐₓ, m′ₘₐₓ, mₘₐₓ) == WignerDsize(ℓₘₐₓ, mₘₐₓ, m′ₘₐₓ)

#                 (m₁, m₂) = extrema((m′ₘₐₓ, mₘₐₓ))

#                 @test WignerDsize(ℓₘₐₓ, m′ₘₐₓ, mₘₐₓ) == (
#                     sum(((2ℓ+1)^2 for ℓ ∈ 0:m₁), init=0)
#                     + sum(((2m₁+1)*(2ℓ+1) for ℓ ∈ m₁+1:m₂); init=0)
#                     + sum(((2m₁+1)*(2m₂+1) for ℓ ∈ m₂+1:ℓₘₐₓ); init=0)
#                 )

#                 data = indices_vector(ℓₘₐₓ, m′ₘₐₓ, mₘₐₓ)
#                 for ℓ ∈ 0:ℓₘₐₓ-1
#                     @test data[WignerDsize(ℓ, m′ₘₐₓ, mₘₐₓ)] == (ℓ, min(m′ₘₐₓ, ℓ), min(mₘₐₓ, ℓ))
#                 end

#             end
#         end
#     end
# end


# """
#     WignerDMatrices(NT, ℓₘₐₓ; m′ₘₐₓ=ℓₘₐₓ, mₘₐₓ=ℓₘₐₓ)

# Create a `WignerDMatrices` object with the given parameters.  The data is initialized to
# zero.
# """
# function WignerDMatrices(::Type{NT}, ℓₘₐₓ::IT; m′ₘₐₓ::IT=ℓₘₐₓ, mₘₐₓ::IT=ℓₘₐₓ) where {NT, IT}
#     # Massage the inputs
#     mₘₐₓ = abs(mₘₐₓ)
#     m′ₘₐₓ = abs(m′ₘₐₓ)

#     # Check that the parameters are valid
#     if complex(NT) != NT
#         throw(ErrorException("NT=$NT must be a complex type"))
#     end
#     if ℓₘₐₓ < (limit = (IT<:Rational ? 1//2 : 0))
#         throw(ErrorException("ℓₘₐₓ < $limit"))
#     end
#     if m′ₘₐₓ > ℓₘₐₓ
#         throw(ErrorException("m′ₘₐₓ > ℓₘₐₓ"))
#     end
#     if mₘₐₓ > ℓₘₐₓ
#         throw(ErrorException("mₘₐₓ > ℓₘₐₓ"))
#     end

#     # Create the data array
#     data = zeros(NT, WignerDsize(ℓₘₐₓ, m′ₘₐₓ, mₘₐₓ))

#     return WignerDMatrices{NT, IT}(data, ℓₘₐₓ, m′ₘₐₓ, mₘₐₓ)
# end


# """
#     index(D, ℓ)

# Find the index in `data(D)` of the first element of the `WignerDMatrix` for the given ℓ
# value.
# """
# function index(D, ℓ)
#     if ℓ < ℓₘᵢₙ(D) || ℓ > ℓₘₐₓ(D)
#         throw(ErrorException("ℓ=$ℓ is out of range for D=$D"))
#     end

#     if ℓ == ℓₘᵢₙ(D)
#         1
#     else
#         WignerDsize(ℓ-1, m′ₘₐₓ(D), mₘₐₓ(D)) + 1
#     end
# end


# @testitem "Test WignerDMatrices index" setup=[WignerDUtilities] begin
#     import SphericalFunctions.Redesign: WignerDMatrices, index

#     for ℓₘₐₓ ∈ 0:8
#         for mₘₐₓ ∈ 0:ℓₘₐₓ
#             for m′ₘₐₓ ∈ 0:ℓₘₐₓ
#                 data = indices_vector(ℓₘₐₓ, m′ₘₐₓ, mₘₐₓ)
#                 D = WignerDMatrices(ComplexF64, ℓₘₐₓ; m′ₘₐₓ, mₘₐₓ)
#                 for ℓ ∈ 0:ℓₘₐₓ
#                     @test data[index(D, ℓ)] == (ℓ, -min(m′ₘₐₓ, ℓ), -min(mₘₐₓ, ℓ))
#                 end

#             end
#         end
#     end
# end


# """
#     size(D)

# Return the total size of the data stored in this WignerDMatrices object, ranging over all
# matrices for all ℓ values.  For the size of a particular matrix, use `size(D, ℓ)`.
# """
# Base.size(D::WignerDMatrices) = WignerDsize(ℓₘₐₓ(D), m′ₘₐₓ(D), mₘₐₓ(D))


# """
#     size(D, ℓ)

# Return the size of the data stored in this WignerDMatrices object for a particular ℓ value.
# For the size of all matrices combined, use `size(D)`.
# """
# function Base.size(D::WignerDMatrices{NT, IT}, ℓ::IT) where {NT, IT}
#     if ℓ < ℓₘᵢₙ(D) || ℓ > ℓₘₐₓ(D)
#         0
#     else
#         return (Int(2m′ₘₐₓ(D, ℓ)) + 1) * (Int(2mₘₐₓ(D, ℓ)) + 1)
#     end
# end

# function Base.getindex(D::WignerDMatrices{NT, IT}, ℓ::IT) where {NT, IT<:Rational}
#     throw(ErrorException("Don't yet know how to deal with Rational indices"))
# end

# function Base.getindex(D::WignerDMatrices{NT, IT}, ℓ::IT) where {NT, IT<:Integer}
#     i₁ = index(D, ℓ)
#     i₂ = i₁ + size(D, ℓ) - 1
#     m′ = m′ₘₐₓ(D, ℓ)
#     m = mₘₐₓ(D, ℓ)
#     OffsetArrays.Origin(-m′, -m)(reshape((@view data(D)[i₁:i₂]), 2m′+1, 2m+1))
# end


# @testitem "Test WignerDMatrices indices" setup=[WignerDUtilities] begin
#     import SphericalFunctions.Redesign: WignerDMatrices, index

#     for ℓₘₐₓ ∈ 0:8
#         for mₘₐₓ ∈ 0:ℓₘₐₓ
#             for m′ₘₐₓ ∈ 0:ℓₘₐₓ
#                 data = indices_vector(ℓₘₐₓ, m′ₘₐₓ, mₘₐₓ)
#                 D = WignerDMatrices{eltype(data), Int}(
#                     data, ℓₘₐₓ, m′ₘₐₓ, mₘₐₓ
#                 )

#                 for ℓ ∈ 0:ℓₘₐₓ
#                     Dˡ = D[ℓ]
#                     @test size(Dˡ) == (2min(m′ₘₐₓ, ℓ)+1, 2min(mₘₐₓ, ℓ)+1)

#                     for m ∈ -min(mₘₐₓ, ℓ):min(mₘₐₓ, ℓ)
#                         for m′ ∈ -min(m′ₘₐₓ, ℓ):min(m′ₘₐₓ, ℓ)
#                             @test Dˡ[m′, m] == (ℓ, m′, m)
#                         end
#                     end
#                 end
#             end
#         end
#     end
# end
