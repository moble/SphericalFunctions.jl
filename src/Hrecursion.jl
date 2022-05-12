# # Return flat index into arrray of (n, m) pairs.
# # Assumes array is ordered as
# #   [
# #     (n, m)
# #     for n in 0:n_max
# #     for m in -n:n
# # ]
# @inline nm_index(n, m) = m + n * (n + 1) + 1


# function H!(
#     Hwedge::AbstractVector{T}, Hextra::AbstractVector{T}, Hv::AbstractVector{T},
#     ℓₘᵢₙ, ℓₘₐₓ, expiβ::Complex{T}
# ) where {T<:Real}
#     H!(Hwedge, Hextra, Hv, ℓₘᵢₙ, ℓₘₐₓ, ℓₘₐₓ, expiβ)
# end

# function H!(
#     Hwedge::AbstractVector{T}, Hextra::AbstractVector{T}, Hv::AbstractVector{T},
#     ℓₘₐₓ, expiβ::Complex{T}
# ) where {T<:Real}
#     H!(Hwedge, Hextra, Hv, 0, ℓₘₐₓ, ℓₘₐₓ, expiβ)
# end

# function H(ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, expiβ::Complex{T}) where {T<:Real}
#     Hwedge = zeros(T, WignerHsize(m′ₘₐₓ, ℓₘₐₓ))
#     Hv = zeros(T, (ℓₘₐₓ + 1)^2)
#     Hextra = zeros(T, ℓₘₐₓ + 2)
#     H!(Hwedge, Hextra, Hv, ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, expiβ)
# end

# function H(ℓₘᵢₙ, ℓₘₐₓ, expiβ::Complex{T}) where {T<:Real}
#     H(ℓₘᵢₙ, ℓₘₐₓ, ℓₘₐₓ, expiβ)
# end

# function H(ℓₘₐₓ, expiβ::Complex{T}) where {T<:Real}
#     H(0, ℓₘₐₓ, ℓₘₐₓ, expiβ)
# end


# function H!(
#     H::AbstractVector{TU}, expiβ::Complex{T}, ℓₘₐₓ::Integer, m′ₘₐₓ::Integer,
#     (a,b,d), Hindex=WignerHindex
# ) where {TU<:Union{T, Complex{T}}, T<:Real}

"""
    H!(H, expiβ, ℓₘₐₓ, m′ₘₐₓ, H_rec_coeffs)
    H!(H, expiβ, ℓₘₐₓ, m′ₘₐₓ, H_rec_coeffs, Hindex)

Compute the ``H`` matrix defined by Gumerov and Duraiswami.

This computation forms the basis for computing Wigner's ``d`` and ``𝔇``
matrices via [`d!`](@ref) and [`D!`](@ref), the spin-weighted spherical
harmonics via [`Y!`](@ref), and for transforming from values of spin-weighted
spherical functions evaluated on a grid to the corresponding mode weights via
[`map2salm`](@ref).

Due to symmetries, we only need to compute ~1/4 of the elements of this matrix,
so only those elements with ``m ≥ |m′|`` are computed.  The relevant indices of
the `H` vector are computed based on the `Hindex` function — which defaults to
`WignerHindex`, but could reasonably be `WignerDindex` if the input `H` vector
contains all valid indices.  However, it is assumed that the storage scheme
used for `H` is such that the successive ``m`` values are located in successive
elements.

If ``m′ₘₐₓ < ℓₘₐₓ``, we don't even need 1/4 of the elements, and only values
with ``|m′| ≤ m′ₘₐₓ`` will be computed.  This is particularly useful for
computing spin-weighted spherical harmonics.

Note that the recursion coefficients `H_rec_coeffs` should be the quantity
returned by [`H_recursion_coefficients`](@ref).

"""
function H2!(HC::HCalculator{T}, expiβ::Complex{T}, ℓ) where {T<:Real}
    @assert ℓ ≤ HC.ℓₘₐₓ
    ℓₘₐₓ = HC.ℓₘₐₓ
    m′ₘₐₓ = HC.m′ₘₐₓ
    Hₙ₊₁⁰ = HC.Hₙ₊₁⁰
    Hₙ = HC.Hₙ
    Hₙ₋₁⁰ = HC.Hₙ₋₁⁰
    (aₙᵐ,bₙᵐ,dₙᵐ) = HC.H_rec_coeffs

    sqrt3 = √T(3)
    invsqrt2 = inv(√T(2))
    cosβ = expiβ.re
    sinβ = expiβ.im
    cosβ₊ = (1+cosβ)/2  # = cos²(β/2)
    cosβ₋ = (1-cosβ)/2  # = sin²(β/2)

    @inbounds begin
        if ℓ == 0

            # Step 1
            Hₙ[1] = 1
            # Step 2 is moot for Hₙ
            # Step 2 for Hₙ₊₁⁰
            Hₙ₊₁⁰[1] = cosβ
            Hₙ₊₁⁰[2] = invsqrt2 * sinβ
            # Steps 3—5 are moot

        elseif ℓ == 1

            # Step 1 is skipped
            # Step 2 is finished for Hₙ by copying
            # Step 2 for Hₙ₊₁⁰
            let i_100 = ifelse(m′ₘₐₓ == 0, 1, 2)
                b̄₂ = invsqrt2
                c̄₂₁ = sqrt3 / 2
                ē₂₁ = c̄₂₁ * invsqrt2
                ē₂₂ = c̄₂₁
                Hₙ₊₁⁰[1] = cosβ * Hₙ[i_100] - b̄₂ * sinβ * Hₙ[i_100+1]
                Hₙ₊₁⁰[2] = c̄₂₁ * cosβ * Hₙ[i_100+1] + sinβ * ē₂₁ * Hₙ[i_100]
                Hₙ₊₁⁰[3] = sinβ * ē₂₂ * Hₙ[i_100+1]
            end
            if m′ₘₐₓ > 0
                # Step 3
                # b^{0}_{2} H^{1, 1}_{1}
                #   = b^{−2}_{2} cosβ₋ H^{0, 2}_{2}
                #   − b^{0}_{2} cosβ₊ H^{0, 0}_{2}
                #   − a^{1}_{1} \sin β H^{0, 1}_{2}
                # Hₙ[4] = (b(2,-2) cosβ₋ Hₙ₊₁⁰[3] − b(2,0) cosβ₊ Hₙ₊₁⁰[1] − a(1,1) sinβ Hₙ₊₁⁰[2]) / b(2,0)
                #       = (b(2,-2) / b(2,0) cosβ₋ Hₙ₊₁⁰[3] − a(1,1) / b(2,0) sinβ Hₙ₊₁⁰[2]) − cosβ₊ Hₙ₊₁⁰[1]
                # a(1,1) = inv(√T(5))
                # b(2,0) = inv(√T(5)) / sqrt3 / invsqrt2
                # b(2,-2) = -2inv(√T(5))
                # b(2,-2) / b(2,0) = -2 sqrt3 * invsqrt2
                # a(1,1) / b(2,0) = sqrt3 * invsqrt2
                Hₙ[4] = (
                    invsqrt2 * sqrt3 * (-2cosβ₋ * Hₙ₊₁⁰[3] − sinβ * Hₙ₊₁⁰[2])
                    − cosβ₊ * Hₙ₊₁⁰[1]
                )
                # Step 4 is finished by step 3
                # Step 5
                Hₙ[1] = -(Hₙ[4] + Hₙ[2])
                # d(n,m) = (m<0 ? -1 : 1) * (√T((n-m)*(n+m+1))) / 2
                # d(1,0) = invsqrt2
                # d(1,-1) = -invsqrt2
                # H[1, -1, 1] = d[1, 0] (H[1, 1, 1] + H[1, 0, 0]) / d[1, -1]
            end

        else

            @info "Skipping H calculation for ℓ=$ℓ"

        end  # if ℓₘₐₓ > 0
    end  # @inbounds
    Hₙ
end
