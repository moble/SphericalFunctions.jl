dfunc(n, m) = ifelse(m<0, -1, 1) * √((n-m)*(n+m+1))

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
function H2!(HC::HCalculator{T}, expiβ::Complex{T}, n) where {T<:Real}
    @assert n ≤ HC.ℓₘₐₓ
    m′ₘₐₓ = HC.m′ₘₐₓ
    Hₙ₊₁⁰ = HC.Hₙ₊₁⁰
    Hₙ = HC.Hₙ

    sqrt3 = √T(3)
    invsqrt2 = inv(√T(2))
    inv2n = inv(T(2n))
    inv2np1 = inv(T(2n+2))
    cosβ = expiβ.re
    sinβ = expiβ.im
    cosβ₊ = (1+cosβ)/2  # = cos²(β/2)
    cosβ₋ = (1-cosβ)/2  # = sin²(β/2)

    @debug "Remove indices and replace @inbounds"
    indices = [
        (n, m′, m)
        for m′ in -min(n, m′ₘₐₓ):min(n, m′ₘₐₓ)
        for m in abs(m′):n
    ]
    #@inbounds begin
    begin
        if n == 0

            # Step 1
            Hₙ[1] = 1
            # Step 2 is moot for Hₙ
            # Step 2 for Hₙ₊₁⁰
            Hₙ₊₁⁰[1] = cosβ
            Hₙ₊₁⁰[2] = invsqrt2 * sinβ
            # Steps 3—5 are moot

        elseif n == 1

            # Step 1 is irrelevant
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
                Hₙ[4] = (
                    invsqrt2 * sqrt3 * (-2cosβ₋ * Hₙ₊₁⁰[3] − sinβ * Hₙ₊₁⁰[2])
                    − cosβ₊ * Hₙ₊₁⁰[1]
                )
                # Step 4 is finished by step 3
                # Step 5
                Hₙ[1] = -(Hₙ[4] + Hₙ[2])
            end

        else  # n > 1
            @warn "H calculation probably incomplete for n=$n"

            iₙ = WignerHindex(n, 0, 0, m′ₘₐₓ) - WignerHindex(n-1, min(m′ₘₐₓ, n-1), n-1, m′ₘₐₓ)
            iₙ0 = iₙ

            # Step 1 is irrelevant

            # Step 2 is finished for Hₙ by copying

            # Step 2 for Hₙ₊₁⁰
            b̄ₙ = √(T(n)/(n+1))
            # H^{0,0}_{n} = cosβ H^{0,0}_{n-1} - b̄ₙ sinβ H^{0,1}_{n-1}
            @show (n, iₙ)
            Hₙ₊₁⁰[1] = cosβ * Hₙ[iₙ] - b̄ₙ  * sinβ * Hₙ[iₙ+1]
            # H^{0,m}_{n} = c̄ₙₘ cosβ H^{0,m}_{n-1} - sinβ [d̄ₙₘ H^{0,m+1}_{n-1} - ēₙₘ H^{0,m-1}_{n-1}]
            for m in 1:n-1
                iₙ += 1
                c̄ₙₘ = 2inv2np1 * √T((n+1+m)*(n+1-m))
                d̄ₙₘ = inv2np1 * √T((n+1-m)*(n-m))
                ēₙₘ = inv2np1 * √T((n+1+m)*(n+m))
                Hₙ₊₁⁰[1+m] = (
                    c̄ₙₘ * cosβ * Hₙ[iₙ]
                    - sinβ * (
                        d̄ₙₘ * Hₙ[iₙ+1]
                        - ēₙₘ * Hₙ[iₙ-1]
                    )
                )
            end
            let m = n
                iₙ += 1
                c̄ₙₘ = 2inv2np1 * √T(2n+1)
                #d̄ₙₘ = 0
                ēₙₘ = inv2np1 * √T((2n+1)*(2n))
                Hₙ₊₁⁰[1+m] = (
                    c̄ₙₘ * cosβ * Hₙ[iₙ]
                    + sinβ * ēₙₘ * Hₙ[iₙ-1]
                )
            end
            let m = n+1
                iₙ += 1
                #c̄ₙₘ = 0
                #d̄ₙₘ = 0
                ēₙₘ = inv2np1 * √T((2n+2)*(2n+1))
                Hₙ₊₁⁰[1+m] = sinβ * ēₙₘ * Hₙ[iₙ-1]
            end

            if m′ₘₐₓ > 0
                # Step 3: Compute H^{1,m}_{n}(β) for m=1,...,n
                # iₙ is now pointing at H^{1, 1}_{n}
                invsqrtnnp1 = inv(√T(n*(n+1)))
                for m in 1:n
                    a = √T((n+m+1)*(n-m+1))
                    b1 = √T((n+m+1)*(n+m+2))
                    b2 = √T((n-m+1)*(n-m+2))
                    # @show (n, 1, m, iₙ)
                    @show (1, (n, 1, m), indices[iₙ])
                    Hₙ[iₙ] = -invsqrtnnp1 * (
                        b1 * cosβ₋ * Hₙ₊₁⁰[m+2]
                        + b2 * cosβ₊ * Hₙ₊₁⁰[m]
                        + a * sinβ * Hₙ₊₁⁰[m+1]
                    )
                    iₙ += 1
                end

                # Step 4: Compute H^{m′+1, m}_{n}(β) for m′=1,...,n−1, m=m′,...,n
                d1 = √T(n*(n+1))
                for m′ in 1:min(n, m′ₘₐₓ)-1
                    iₙ′ = iₙ - (n-m′)
                    # iₙ points at H^{m′+1, m}_{n}
                    # iₙ′ points at H^{m′, m}_{n}
                    # d^{m′}_{n} H^{m′+1, m}_{n}
                    #     = d^{m′−1}_{n} H^{m′−1, m}_{n}
                    #     − d^{m−1}_{n} H^{m′, m−1}_{n}
                    #     + d^{m}_{n} H^{m′, m+1}_{n}
                    # (where the last term drops out for m=n). The constants are defined by
                    # d^{m}_{n} = \frac{\mathrm{sgn}(m)}{2} \sqrt{(n-m)(n+m+1)}.
                    # d1 ≔ d^{m′}_{n} = √T((n-m′)*(n+m′+1))
                    # d2 ≔ d^{m′-1}_{n} = √T((n-m′+1)*(n+m′))
                    # d3 ≔ d^{m-1}_{n} = √T((n-m+1)*(n+m))
                    # d4 ≔ d^{m}_{n} = √T((n-m)*(n+m+1))
                    d2 = d1
                    d1 = √T((n-m′)*(n+m′+1))
                    d4 = d1
                    invd1 = inv(d1)
                    for m in m′+1:n-1
                        d3 = d4
                        d4 = √T((n-m)*(n+m+1))
                        Hₙ[iₙ] = invd1 * (
                            d2 * Hₙ[iₙ′ - (n-m′+1)]
                            - d3 * Hₙ[iₙ′-1]
                            + d4 * Hₙ[iₙ′+1]
                        )
                        iₙ += 1
                        iₙ′ += 1
                    end
                    let m = n
                        d3 = √T(2n)
                        Hₙ[iₙ] = invd1 * (
                            d2 * Hₙ[iₙ′ - (n-m′+1)]
                            - d3 * Hₙ[iₙ′-1]
                        )
                        iₙ += 1
                        iₙ′ += 1
                    end
                end  # Step 4

                # Step 5: Compute H^{m′−1, m}_{n}(β) for m′=0,...,−n+1, m=−m′,...,n
                iₙ = iₙ0 - 1
                d2 = √T(n*(n+1))
                offset = 0  # Account for 2 missing elements when m′==0
                # println()
                for m′ in 0:-1:-min(n, m′ₘₐₓ)+1
                    iₙ′ = iₙ + (n+m′) + 1
                    # @show (iₙ, (n,m′-1,n), indices[iₙ])
                    # @show (iₙ′, (n,m′,n), indices[iₙ′])
                    # println()
                    # iₙ points at H^{m′-1, m}_{n}
                    # iₙ′ points at H^{m′, m}_{n}
                    # d^{m′−1}_{n} H^{m′−1, m}_{n}
                    #   = d^{m′}_{n} H^{m′+1, m}_{n}
                    #     + d^{m−1}_{n} H^{m′, m−1}_{n}
                    #     − d^{m}_{n} H^{m′, m+1}_{n}
                    # (where the last term drops out for m=n). The constants are defined by
                    # d^{m}_{n} = \mathrm{sgn}(m) \sqrt{(n-m)(n+m+1)}
                    # d1 ≔ d^{m′-1}_{n} = -√T((n-m′+1)*(n+m′))
                    # d2 ≔ d^{m′}_{n} = √T((n-m′)*(n+m′+1))
                    # d3 ≔ d^{m-1}_{n} = √T((n-m+1)*(n+m))
                    # d4 ≔ d^{m}_{n} = √T((n-m)*(n+m+1))

                    d1 = -√T((n-m′+1)*(n+m′))
                    invd1 = inv(d1)
                    d3 = √T(2n)
                    let m = n
                        # println("A")
                        # @show (n, m′, m)
                        # @show (iₙ, (n, m′−1, m), indices[iₙ])
                        # @show (iₙ′ + (n+m′+offset), (n, m′+1, m), indices[iₙ′ + (n+m′+offset)])
                        # @show (iₙ′ + (n+m′+0), (n, m′+1, m), indices[iₙ′ + (n+m′+0)])
                        # @show (iₙ′-1, (n, m′, m−1), indices[iₙ′-1])
                        # #@show (iₙ′+1, (n, m′, m+1), indices[iₙ′+1])
                        # println()
                        # @show (d1, d2, d3)
                        # @show (dfunc(n, m′−1), dfunc(n, m′), dfunc(n, m−1))
                        # @show (invd1, d2, Hₙ[iₙ′ + (n+m′)], d3, Hₙ[iₙ′-1])
                        # println()
                        Hₙ[iₙ] = invd1 * (
                            d2 * Hₙ[iₙ′ + (n+m′+offset)]
                            + d3 * Hₙ[iₙ′-1]
                        )
                        iₙ -= 1
                        iₙ′ -= 1
                    end
                    for m in n-1:-1:1-m′
                        d4 = d3
                        d3 = √T((n-m+1)*(n+m))
                        # println("B")
                        # @show (n, m′, m)
                        # @show (iₙ, (n, m′−1, m), indices[iₙ])
                        # @show (iₙ′ + (n+m′+offset), (n, m′+1, m), indices[iₙ′ + (n+m′+offset)])
                        # @show (iₙ′-1, (n, m′, m−1), indices[iₙ′-1])
                        # @show (iₙ′+1, (n, m′, m+1), indices[iₙ′+1])
                        # println()
                        # @show (d1, d2, d3, d4)
                        # @show (dfunc(n, m′−1), dfunc(n, m′), dfunc(n, m−1), dfunc(n, m))
                        # println()
                        Hₙ[iₙ] = invd1 * (
                            d2 * Hₙ[iₙ′ + (n+m′+offset)]
                            + d3 * Hₙ[iₙ′-1]
                            - d4 * Hₙ[iₙ′+1]
                        )
                        iₙ -= 1
                        iₙ′ -= 1
                    end
                    offset = 2
                    d2 = d1
                end  # Step 5

            end


        end  # if nₘₐₓ > 0
    end  # @inbounds
    Hₙ
end
