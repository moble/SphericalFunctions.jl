module AssociatedLegendreFunction

import ..OffsetVec, ..OffsetMat

### The code in this module is based on a paper by Xing et al.:
### https://doi.org/10.1007/s00190-019-01331-0.  All references to equation numbers are for that
### paper.


# These functions implement Eqs. (13), absorbing a factor of 1/2 into b, and converting to the
# relevant float type where needed.
@inline a(n, ::Type{T}) where {T<:Real} = √((2n + 1) / T(2n - 1))
@inline b(n, ::Type{T}) where {T<:Real} = √((n - 1) * (2n + 1) / T(2n * (2n - 1)))
@inline c(n, m, ::Type{T}) where {T<:Real} = √((n + m) * (n - m) * (2n + 1) / T(2n - 1)) / n
@inline d(n, m, ::Type{T}) where {T<:Real} = √((n - m) * (n - m - 1) * (2n + 1) / T(2n - 1)) / (2n)
@inline e(n, m, ::Type{T}) where {T<:Real} = √((n + m) * (n + m - 1) * (m != 1 ? 1 : 2) * (2n + 1) / T(2n - 1)) / (2n)

# These are more efficient to compute than the above, because they involve no division and fewer T-multiplications
@inline b̂(n, ::Type{T}) where {T<:Real} = √T(2n * (n - 1))  # b = b̂ * a / 2n
@inline ĉ(n, m, ::Type{T}) where {T<:Real} = √T(4 * (n + m) * (n - m))  # c = ĉ * a / 2n
@inline d̂(n, m, ::Type{T}) where {T<:Real} = √T((n - m) * (n - m - 1))  # d = d̂ * a / 2n
@inline ê(n, m, ::Type{T}) where {T<:Real} = √T((n + m) * (n + m - 1) * (m != 1 ? 1 : 2))  # e = ê * a / 2n


struct ALFRecursionCoefficients{T<:Real}
    nmax::Int
    a::Vector{T}
    b::Vector{T}
    cde::Matrix{T}
end
function ALFRecursionCoefficients(nmax, ::Type{T}=Float64) where {T<:Real}
    avalues = Vector{T}(undef, nmax + 1)
    bvalues = Vector{T}(undef, nmax + 1)
    cdevalues = Matrix{T}(undef, 3, (nmax * (nmax + 3)) ÷ 2 + 1)
    i = 0
    @inbounds for n in 1:(nmax+1)
        avalues[n] = a(n, T)
        bvalues[n] = b(n, T)
        for m in 1:n
            i += 1
            cdevalues[1, i] = c(n, m, T)
            cdevalues[2, i] = d(n, m, T)
            cdevalues[3, i] = e(n, m, T)
        end
    end
    ALFRecursionCoefficients{T}(nmax, avalues, bvalues, cdevalues)
end
Base.show(io::IO, r::ALFRecursionCoefficients) = print(io, "$(typeof(r))($(r.nmax))")


function ALFrecurse!(P̄ₙ::AbstractVector{T}, P̄ₙ₋₁::AbstractVector{T}, n::Int, t::T, u::T, aₙ::T, bₙ::T, cdeₙ::AbstractMatrix{T}) where {T<:Real}
    @fastmath @inbounds begin
        let m = 0
            # Eq. (12a); note that we absorbed a factor of 1/2 into the definition of b
            P̄ₙ[m] = t * aₙ * P̄ₙ₋₁[m] - u * bₙ * P̄ₙ₋₁[m+1]
        end
        for m in 1:n-2
            # Eq. (12b)
            cₙₘ = cdeₙ[1, m]
            dₙₘ = cdeₙ[2, m]
            eₙₘ = cdeₙ[3, m]
            P̄ₙ[m] = t * cₙₘ * P̄ₙ₋₁[m] - u * (dₙₘ * P̄ₙ₋₁[m+1] - eₙₘ * P̄ₙ₋₁[m-1])
        end
        let m = n - 1
            # Eq. (12b)
            cₙₘ = cdeₙ[1, m]
            eₙₘ = cdeₙ[3, m]
            P̄ₙ[m] = t * cₙₘ * P̄ₙ₋₁[m] + u * eₙₘ * P̄ₙ₋₁[m-1]
        end
        let m = n
            # Eq. (12b)
            eₙₘ = cdeₙ[3, m]
            P̄ₙ[m] = u * eₙₘ * P̄ₙ₋₁[m-1]
        end
    end
end


function ALFrecurse!(P̄ₙ::AbstractVector{T}, P̄ₙ₋₁::AbstractVector{T}, n::Int, t::T, u::T) where {T<:Real}
    @fastmath @inbounds begin
        aₙ = a(n, T)
        aₙ_over_2n = aₙ / T(2n)
        bₙ = b̂(n, T) * aₙ_over_2n
        let m = 0
            # Eq. (12a); note that we absorbed a factor of 1/2 into the definition of b
            P̄ₙ[m] = t * aₙ * P̄ₙ₋₁[m] - u * bₙ * P̄ₙ₋₁[m+1]
        end
        for m in 1:n-2
            # Eq. (12b)
            ĉₙₘ = ĉ(n, m, T)
            d̂ₙₘ = d̂(n, m, T)
            êₙₘ = ê(n, m, T)
            P̄ₙ[m] = (t * ĉₙₘ * P̄ₙ₋₁[m] - u * (d̂ₙₘ * P̄ₙ₋₁[m+1] - êₙₘ * P̄ₙ₋₁[m-1])) * aₙ_over_2n
        end
        let m = n - 1
            # Eq. (12b)
            ĉₙₘ = ĉ(n, m, T)
            êₙₘ = ê(n, m, T)
            P̄ₙ[m] = (t * ĉₙₘ * P̄ₙ₋₁[m] + u * êₙₘ * P̄ₙ₋₁[m-1]) * aₙ_over_2n
        end
        let m = n
            # Eq. (12b)
            êₙₘ = ê(n, m, T)
            P̄ₙ[m] = u * êₙₘ * P̄ₙ₋₁[m-1] * aₙ_over_2n
        end
    end
end


"""
    ALFcompute(expiβ, nmax)
    ALFcompute!(P̄, expiβ, nmax)
    ALFcompute(expiβ, nmax, recursion_coefficients)
    ALFcompute!(P̄, expiβ, nmax, recursion_coefficients)

Compute the "fully normalized" Associated Legendre Functions up to some maximum `n` value `nmax`.

These functions can take a vector P̄, to store the data, stored in order of increasing `m` most
rapidly varying and then increasing `n`.  If not supplied, P̄ will be constructed for you and
returned.

The optional `recursion_coefficients` argument must be an `ALFRecursionCoefficients`, which stores
various constant coefficients used in the recursion.  This object requires more than 3x the memory
and more than 20x the time to compute a single P̄ vector without this argument, but passing it will
typically speed up the calculation of each P̄ by a factor of 8x or so.  Thus, if you expect to
compute P̄ more than a few times, it will take less time to pre-compute those factors, and pass them
to this function.

Note that the base real type will be inferred from the (complex) type of `expiβ`.  If present, the
base types of `P̄` and `recursion_coefficients` must agree.

"""
function ALFcompute!(P̄::Vector{T}, expiβ::Complex{T}, nmax::Int, recursion_coefficients::ALFRecursionCoefficients{T}) where {T<:Real}
    min_length = (nmax * (nmax + 3)) ÷ 2 + 1
    if length(P̄) < min_length
        throw(ArgumentError("Length of P̄ ($(length(P̄))) must be at least $min_length for nmax=$nmax"))
    end

    P̄ₙ₋₁ = OffsetVec(P̄, -1)  # -(n*(n-1))÷2-1
    P̄ₙ = OffsetVec(P̄, -1)  # -(n*(n+1))÷2-1
    cdeₙ = OffsetMat(recursion_coefficients.cde, 0, 0)  # -(n*(n-1))÷2

    @inbounds begin
        t = expiβ.re  # = cosβ
        u = expiβ.im  # = sinβ

        # Initialize with Eq. (14)
        if nmax ≥ 0
            P̄ₙ[0] = one(T)  # n=0
        end
        if nmax ≥ 1
            P̄ₙ.offset -= 1
            sqrt3 = √T(3)
            P̄ₙ[0] = sqrt3 * t
            P̄ₙ[1] = sqrt3 * u
        end
        for n in 2:min(nmax, recursion_coefficients.nmax)
            aₙ = recursion_coefficients.a[n]
            bₙ = recursion_coefficients.b[n]
            cdeₙ.offset2 -= (n - 1)
            P̄ₙ₋₁.offset = P̄ₙ.offset
            P̄ₙ.offset -= n
            ALFrecurse!(P̄ₙ, P̄ₙ₋₁, n, t, u, aₙ, bₙ, cdeₙ)
        end
        if nmax > recursion_coefficients.nmax
            for n in max(2, recursion_coefficients.nmax + 1):nmax
                P̄ₙ₋₁.offset = P̄ₙ.offset
                P̄ₙ.offset -= n
                ALFrecurse!(P̄ₙ, P̄ₙ₋₁, n, t, u)
            end
        end
    end
end
function ALFcompute!(P̄::Vector{T}, expiβ::Complex{T}, nmax::Int) where {T<:Real}
    recursion_coefficients = ALFRecursionCoefficients(0, T)
    ALFcompute!(P̄, expiβ, nmax, recursion_coefficients)
end


function ALFcompute(expiβ::Complex{T}, nmax::Int, recursion_coefficients::ALFRecursionCoefficients{T}) where {T<:Real}
    min_length = (nmax * (nmax + 3)) ÷ 2 + 1
    P̄ = Vector{T}(undef, min_length)
    ALFcompute!(P̄, expiβ, nmax, recursion_coefficients)
    P̄
end
function ALFcompute(expiβ::Complex{T}, nmax::Int) where {T<:Real}
    recursion_coefficients = ALFRecursionCoefficients(0, T)
    ALFcompute(expiβ, nmax, recursion_coefficients)
end


end  # module ALF

using .AssociatedLegendreFunction: ALFRecursionCoefficients, ALFrecurse!, ALFcompute!
