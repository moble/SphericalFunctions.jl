module AssociatedLegendreFunction

### The code in this module is based on a paper by Xing et al.:
### https://doi.org/10.1007/s00190-019-01331-0.  All references to equation numbers are for that
### paper.

### Note: my experiments show that
###  1. caching is about 3x faster than not caching
###  2. using a square matrix is about 23% faster than a flattened vector
###  3. using a square matrix with [m, n] ordering is about 5x faster than [n, m] ordering
###     (which can be done either in the code or by transposing the input matrix).


using OffsetArrays


# These functions implement Eqs. (13), absorbing a factor of 1/2 into b, and converting to the
# relevant float type where needed.
@inline a(n, ::Type{T}) where {T<:Real} = √((2n+1)/T(2n-1))
@inline b(n, ::Type{T}) where {T<:Real} = √((n-1)*(2n+1)/T(2n*(2n-1)))
@inline c(n, m, ::Type{T}) where {T<:Real} = √((n+m)*(n-m)*(2n+1) / T(2n-1)) / n
@inline d(n, m, ::Type{T}) where {T<:Real} = √((n-m)*(n-m-1)*(2n+1) / T(2n-1)) / (2n)
@inline e(n, m, ::Type{T}) where {T<:Real} = √((n+m)*(n+m-1)*(m!=1 ? 1 : 2)*(2n+1) / T(2n-1)) / (2n)

# These are more efficient to compute than the above, because they involve no division and fewer T-multiplications
@inline b̂(n, ::Type{T}) where {T<:Real} = √T(2n*(n-1))  # b = b̂ * a / 2n
@inline ĉ(n, m, ::Type{T}) where {T<:Real} = √T(4*(n+m)*(n-m))  # c = ĉ * a / 2n
@inline d̂(n, m, ::Type{T}) where {T<:Real} = √T((n-m)*(n-m-1))  # d = d̂ * a / 2n
@inline ê(n, m, ::Type{T}) where {T<:Real} = √T((n+m)*(n+m-1)*(m!=1 ? 1 : 2))  # e = ê * a / 2n



struct ALFRecursionCoefficients{T<:Real}
    nmax::Int
    a::Vector{T}
    b::Vector{T}
    cde::Matrix{T}
end
ALFRecursionCoefficients(nmax, T=Float64) = ALFRecursionCoefficients{T}(
    nmax,
    [a(n, T) for n in 1:nmax+1],
    [b(n, T) for n in 1:nmax+1],
    reshape([f(n, m, T) for n in 1:nmax+1 for m in 1:n for f in [c, d, e]], 3, :)
)

function Base.show(io::IO, r::ALFRecursionCoefficients)
    print(io, "$(typeof(r))($(r.nmax))")
end



function ALFrecurse!(P̄ₙ::AbstractVector{T}, P̄ₙ₋₁::AbstractVector{T}, n::Int, t::T, u::T, aₙ::T, bₙ::T, cdeₙ::AbstractMatrix{T}) where {T<:Real}
    @inbounds begin
        let m=0
            # Eq. (12a); note that we absorbed a factor of 1/2 into the definition of b
            P̄ₙ[m] = t * aₙ * P̄ₙ₋₁[m] - u * bₙ * P̄ₙ₋₁[m+1]
        end
        for m in 1:n-2
            # Eq. (12b)
            cₙₘ, dₙₘ, eₙₘ = cdeₙ[:, m]
            P̄ₙ[m] = t * cₙₘ * P̄ₙ₋₁[m] - u * (dₙₘ * P̄ₙ₋₁[m+1] - eₙₘ * P̄ₙ₋₁[m-1])
        end
        let m=n-1
            # Eq. (12b)
            cₙₘ = cdeₙ[1, m]
            eₙₘ = cdeₙ[3, m]
            P̄ₙ[m] = t * cₙₘ * P̄ₙ₋₁[m] + u * eₙₘ * P̄ₙ₋₁[m-1]
        end
        let m=n
            # Eq. (12b)
            eₙₘ = cdeₙ[3, m]
            P̄ₙ[m] = u * eₙₘ * P̄ₙ₋₁[m-1]
        end
    end
end


function ALFrecurse!(P̄ₙ::AbstractVector{T}, P̄ₙ₋₁::AbstractVector{T}, n::Int, t::T, u::T) where {T<:Real}
    @inbounds begin
        aₙ = a(n, T)
        aₙ_over_2n = aₙ / T(2n)
        bₙ = b̂(n, T) * aₙ_over_2n
        let m=0
            # Eq. (12a); note that we absorbed a factor of 1/2 into the definition of b
            P̄ₙ[m] = t * aₙ * P̄ₙ₋₁[m] - u * bₙ * P̄ₙ₋₁[m+1]
        end
        for m in 1:n-2
            # Eq. (12b)
            cₙₘ = ĉ(n, m, T) * aₙ_over_2n
            dₙₘ = d̂(n, m, T) * aₙ_over_2n
            eₙₘ = ê(n, m, T) * aₙ_over_2n
            P̄ₙ[m] = t * cₙₘ * P̄ₙ₋₁[m] - u * (dₙₘ * P̄ₙ₋₁[m+1] - eₙₘ * P̄ₙ₋₁[m-1])
        end
        let m=n-1
            # Eq. (12b)
            cₙₘ = ĉ(n, m, T) * aₙ_over_2n
            eₙₘ = ê(n, m, T) * aₙ_over_2n
            P̄ₙ[m] = t * cₙₘ * P̄ₙ₋₁[m] + u * eₙₘ * P̄ₙ₋₁[m-1]
        end
        let m=n
            # Eq. (12b)
            eₙₘ = ê(n, m, T) * aₙ_over_2n
            P̄ₙ[m] = u * eₙₘ * P̄ₙ₋₁[m-1]
        end
    end
end


"""
    ALFcompute!(P̄, expiβ, nmax, recursion_coefficients)

Compute the "fully normalized" Associated Legendre Functions.  This takes a vector P̄, to store the
data, stored in order of increasing `m` most rapidly varying and then increasing `n`.

"""
function ALFcompute!(P̄::Vector{T}, expiβ::Complex{T}, nmax::Int, recursion_coefficients::ALFRecursionCoefficients{T}) where {T<:Real}
    min_length = (nmax * (nmax + 3)) ÷ 2 + 1
    if length(P̄) < min_length
        throw(ArgumentError("Length of P̄ ($(length(P̄))) must be at least $min_length for nmax=$nmax"))
    end
    offset = -1

    @inbounds begin
        t = expiβ.re  # = cosβ
        u = expiβ.im  # = sinβ

        # Initialize with Eq. (14)
        if nmax ≥ 0
            P̄₀ = OffsetVector(P̄, offset)
            P̄₀[0] = one(T)
        end
        if nmax ≥ 1
            offset -= 1
            P̄₁ = OffsetVector(P̄, offset)
            sqrt3 = √T(3)
            P̄₁[0] = sqrt3 * t
            P̄₁[1] = sqrt3 * u
        end
        for n in 2:min(nmax, recursion_coefficients.nmax)
            aₙ = recursion_coefficients.a[n]
            bₙ = recursion_coefficients.b[n]
            cdeₙ = OffsetArray(recursion_coefficients.cde, 0, offset+1)  # -(n*(n-1))÷2
            P̄ₙ₋₁ = OffsetVector(P̄, offset)  # -(n*(n-1))÷2-1
            offset -= n
            P̄ₙ = OffsetVector(P̄, offset)  # -(n*(n+1))÷2-1
            ALFrecurse!(P̄ₙ, P̄ₙ₋₁, n, t, u, aₙ, bₙ, cdeₙ)
        end
        if nmax > recursion_coefficients.nmax
            for n in recursion_coefficients.nmax+1:nmax
                P̄ₙ₋₁ = OffsetVector(P̄, offset)  # -(n*(n-1))÷2-1
                offset -= n
                P̄ₙ = OffsetVector(P̄, offset)  # -(n*(n+1))÷2-1
                ALFrecurse!(P̄ₙ, P̄ₙ₋₁, n, t, u)
            end
        end
    end
end







# # Define some convenient abstract types
# abstract type AbstractALFCalculator{nmax, T<:Real} end
# abstract type AbstractCachedALFCalculator{nmax, T<:Real} <: AbstractALFCalculator{nmax, T} end
# abstract type AbstractUncachedALFCalculator{nmax, T<:Real} <: AbstractALFCalculator{nmax, T} end

# struct CachedALFCalculator{nmax, T<:Real} <: AbstractCachedALFCalculator{nmax, T}
#     a::Vector{T}
#     b::Vector{T}
#     c::Vector{T}
#     d::Vector{T}
#     e::Vector{T}
# end
# CachedALFCalculator(nmax, T=Float64) = CachedALFCalculator{nmax, T}(
#     [a(n, T) for n in 1:nmax+1],
#     [b(n, T) for n in 1:nmax+1],
#     [c(n, m, T) for n in 1:nmax+1 for m in 1:n],
#     [d(n, m, T) for n in 1:nmax+1 for m in 1:n],
#     [e(n, m, T) for n in 1:nmax+1 for m in 1:n]
# )


# struct UncachedALFCalculator{nmax, T<:Real} <: AbstractUncachedALFCalculator{nmax, T} end
# UncachedALFCalculator(nmax, T=Float64) = UncachedALFCalculator{nmax, T}();

# @generated T(a::AbstractALFCalculator{nmax, TA}) where {nmax, TA<:Real} = TA
# @generated nmax(a::AbstractALFCalculator{nmaxA, T}) where {nmaxA, T<:Real} = nmaxA

# a_over_2n(a, n, ::AbstractUncachedALFCalculator{nmax, T}) where {nmax, T<:Real} = a / T(2n)
# a_over_2n(a, n, ::AbstractCachedALFCalculator{nmax, T}) where {nmax, T<:Real} = 1

# @inline a(n, ::AbstractUncachedALFCalculator{nmax, T}) where {nmax, T<:Real} = a(n, T)
# @inline b(n, aₙ_over_2n, ::AbstractUncachedALFCalculator{nmax, T}) where {nmax, T<:Real} = b̂(n, T) * aₙ_over_2n
# @inline c(n, m, aₙ_over_2n, ::AbstractUncachedALFCalculator{nmax, T}) where {nmax, T<:Real} = ĉ(n, m, T) * aₙ_over_2n
# @inline d(n, m, aₙ_over_2n, ::AbstractUncachedALFCalculator{nmax, T}) where {nmax, T<:Real} = d̂(n, m, T) * aₙ_over_2n
# @inline e(n, m, aₙ_over_2n, ::AbstractUncachedALFCalculator{nmax, T}) where {nmax, T<:Real} = ê(n, m, T) * aₙ_over_2n

# @inline a(n, alf::AbstractCachedALFCalculator{nmax, T}) where {nmax, T<:Real} = alf.a[n]
# @inline b(n, aₙ_over_2n, alf::AbstractCachedALFCalculator{nmax, T}) where {nmax, T<:Real} = alf.b[n]
# @inline c(n, m, aₙ_over_2n, alf::AbstractCachedALFCalculator{nmax, T}) where {nmax, T<:Real} = alf.c[m + (n * (n - 1)) ÷ 2]
# @inline d(n, m, aₙ_over_2n, alf::AbstractCachedALFCalculator{nmax, T}) where {nmax, T<:Real} = alf.d[m + (n * (n - 1)) ÷ 2]
# @inline e(n, m, aₙ_over_2n, alf::AbstractCachedALFCalculator{nmax, T}) where {nmax, T<:Real} = alf.e[m + (n * (n - 1)) ÷ 2]


# """
#     compute!(P̄, expiβ, alf)

# This takes a full-size matrix, which does not take account of the restriction m≤n.  Note that the
# matrix is transposed relative to what you might expect: P̄ₙₘ = P̄[m, n], which is critical to speed.

# """
# @inbounds function compute!(P̄::Matrix{T}, expiβ::Complex{T}, alf::AbstractALFCalculator{nmax, T}, ::Val{true}) where {nmax, T<:Real}
#     Base.require_one_based_indexing(P̄)
#     if size(P̄)[1] < nmax || size(P̄)[2] < nmax
#         throw(ArgumentError("Shape of P̄ $(size(P̄)) must be at least (nmax, nmax)=($nmax, $nmax)"))
#     end

#     sqrt3 = √T(3)
#     t = expiβ.re  # = cosβ
#     u = expiβ.im  # = sinβ

#     # Initialize with Eq. (14)
#     P̄[0+1, 0+1] = 1
#     P̄[0+1, 1+1] = sqrt3 * t
#     P̄[1+1, 1+1] = sqrt3 * u

#     for n in 2:nmax
#         aₙ = a(n, alf)
#         aₙ_over_2n = a_over_2n(aₙ, n, alf)
#         bₙ = b(n, aₙ_over_2n, alf)
#         let m=0
#             # Eq. (12a); note that we absorbed a factor of 1/2 into the definition of b
#             P̄[0+1, n+1] = t * aₙ * P̄[0+1, n-1+1] - u * bₙ * P̄[1+1, n-1+1]
#         end
# #         for m in 1:n
# #             cₙₘ = c(n, m, aₙ_over_2n, alf)
# #             dₙₘ = d(n, m, aₙ_over_2n, alf)
# #             eₙₘ = e(n, m, aₙ_over_2n, alf)
# #             # Eq. (12b)
# #             P̄[n+1, m+1] = t * cₙₘ * P̄[n-1+1, m+1] - u * (dₙₘ * P̄[n-1+1, m+1+1] - eₙₘ * P̄[n-1+1, m-1+1])
# #         end
#         for m in 1:n-2
#             cₙₘ = c(n, m, aₙ_over_2n, alf)
#             dₙₘ = d(n, m, aₙ_over_2n, alf)
#             eₙₘ = e(n, m, aₙ_over_2n, alf)
#             # Eq. (12b)
#             P̄[m+1, n+1] = t * cₙₘ * P̄[m+1, n-1+1] - u * (dₙₘ * P̄[m+1+1, n-1+1] - eₙₘ * P̄[m-1+1, n-1+1])
#         end
#         let m=n-1
#             cₙₘ = c(n, m, aₙ_over_2n, alf)
# #             dₙₘ = d(n, m, aₙ_over_2n, alf)
#             eₙₘ = e(n, m, aₙ_over_2n, alf)
#             # Eq. (12b)
#             P̄[m+1, n+1] = t * cₙₘ * P̄[m+1, n-1+1] + u * eₙₘ * P̄[m-1+1, n-1+1]
#         end
#         let m=n
# #             cₙₘ = c(n, m, aₙ_over_2n, alf)
# #             dₙₘ = d(n, m, aₙ_over_2n, alf)
#             eₙₘ = e(n, m, aₙ_over_2n, alf)
#             # Eq. (12b)
#             P̄[m+1, n+1] = u * eₙₘ * P̄[m-1+1, n-1+1]
#         end
#     end
# end


# """
#     compute!(P̄, expiβ, alf)

# This takes a vector to store the data, stored in order of increasing `m` most rapidly varying and
# then increasing `n`.

# """
# @inbounds function compute!(P̄::Vector{T}, expiβ::Complex{T}, alf::AbstractALFCalculator{nmax, T}) where {nmax, T<:Real}
#     Base.require_one_based_indexing(P̄)
#     if length(P̄) < nmax
#         throw(ArgumentError("Length of P̄ ($(length(P̄))) must be at least nmax=$nmax"))
#     end

#     sqrt3 = √T(3)
#     t = expiβ.re  # = cosβ
#     u = expiβ.im  # = sinβ

#     # Initialize with Eq. (14)
#     P̄[1] = 1
#     P̄[2] = sqrt3 * t
#     P̄[3] = sqrt3 * u

#     let i=3, i_up=1  # Current (n, m) index, and (n-1, m) index
#         for n in 2:nmax
#             aₙ = a(n, alf)
#             aₙ_over_2n = a_over_2n(aₙ, n, alf)
#             bₙ = b(n, aₙ_over_2n, alf)
#             let m=0
#                 i += 1
#                 i_up += 1
#                 # Eq. (12a); note that we absorbed a factor of 1/2 into the definition of b
#                 P̄[i] = t * aₙ * P̄[i_up] - u * bₙ * P̄[i_up+1]
#             end
#             for m in 1:n-2
#                 cₙₘ = c(n, m, aₙ_over_2n, alf)
#                 dₙₘ = d(n, m, aₙ_over_2n, alf)
#                 eₙₘ = e(n, m, aₙ_over_2n, alf)
#                 # Eq. (12b)
#                 i += 1
#                 i_up += 1
#                 P̄[i] = t * cₙₘ * P̄[i_up] - u * (dₙₘ * P̄[i_up+1] - eₙₘ * P̄[i_up-1])
#             end
#             let m=n-1
#                 cₙₘ = c(n, m, aₙ_over_2n, alf)
#                 #dₙₘ = d(n, m, aₙ_over_2n, alf)
#                 eₙₘ = e(n, m, aₙ_over_2n, alf)
#                 # Eq. (12b)
#                 i += 1
#                 i_up += 1
#                 P̄[i] = t * cₙₘ * P̄[i_up] + u * eₙₘ * P̄[i_up-1]
#             end
#             let m=n
#                 #cₙₘ = c(n, m, aₙ_over_2n, alf)
#                 #dₙₘ = d(n, m, aₙ_over_2n, alf)
#                 eₙₘ = e(n, m, aₙ_over_2n, alf)
#                 # Eq. (12b)
#                 i += 1
#                 i_up += 1  # Actually now pointing to (n, 0), but it's nice for symmetry
#                 P̄[i] = u * eₙₘ * P̄[i_up-1]
#             end
#             i_up -= 1
#         end
#     end
# end


end  # module ALF
