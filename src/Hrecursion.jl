
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

"""
function H2!(WC::WCT, expiβ::Complex{T}, n) where {WCT<:WignerCalculator, T<:Real}
    @assert n ≤ WC.ℓₘₐₓ
    m′ₘₐₓ = WC.m′ₘₐₓ
    Hₙ₊₁⁰ = WC.Hₙ₊₁⁰
    Hₙ = WC.Hₙ
    sqrt3 = WC.sqrt3
    invsqrt2 = WC.invsqrt2
    sqrt2n = √T(2n)
    sqrtnnp1 = √T(n*(n+1))
    inv2np1 = inv(T(2n+2))

    cosβ = expiβ.re
    sinβ = expiβ.im
    cosβ₊ = (1+cosβ)/2  # = cos²(β/2)
    cosβ₋ = (1-cosβ)/2  # = sin²(β/2)

    # @warn "Replace inbounds; remove indices and @show statements"
    # indices = [
    #     (n, m′, m)
    #     for m′ in -min(n, m′ₘₐₓ):min(n, m′ₘₐₓ)
    #     for m in abs(m′):n
    # ]
    # @show n
    # begin
    @inbounds begin
        if n == 0

            # Step 1
            Hₙ[1] = 1
            # Step 2 is moot for Hₙ
            # Step 2 for Hₙ₊₁⁰
            Hₙ₊₁⁰[1] = cosβ
            Hₙ₊₁⁰[2] = invsqrt2 * sinβ
            # Steps 3—5 are moot

        elseif n == 1
            iₙ = ifelse(m′ₘₐₓ == 0, 1, 2)

            # Step 1 is irrelevant for n>0

            # Step 2 for Hₙ is just copying from Hₙ₊₁⁰
            @views Hₙ[iₙ:iₙ+n] .= Hₙ₊₁⁰[1:n+1]

            # Step 2 for Hₙ₊₁⁰
            b̄₂ = invsqrt2
            c̄₂₁ = sqrt3 / 2
            ē₂₁ = c̄₂₁ * invsqrt2
            ē₂₂ = c̄₂₁
            Hₙ₊₁⁰[1] = cosβ * Hₙ[iₙ] - b̄₂ * sinβ * Hₙ[iₙ+1]
            Hₙ₊₁⁰[2] = c̄₂₁ * cosβ * Hₙ[iₙ+1] + sinβ * ē₂₁ * Hₙ[iₙ]
            Hₙ₊₁⁰[3] = sinβ * ē₂₂ * Hₙ[iₙ+1]
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
            #iₙ = WignerHindex(n, 0, 0, m′ₘₐₓ) - WignerHindex(n-1, min(m′ₘₐₓ, n-1), n-1, m′ₘₐₓ)
            iₙ = offset(WC, n, 0, 0)
            iₙ0 = iₙ

            # Step 1 is irrelevant for n>0

            # Step 2 for Hₙ is just copying from Hₙ₊₁⁰
            @views Hₙ[iₙ:iₙ+n] .= Hₙ₊₁⁰[1:n+1]

            # Step 2 for Hₙ₊₁⁰
            b̄ₙ = √(T(n)/(n+1))
            # H^{0,0}_{n} = cosβ H^{0,0}_{n-1} - b̄ₙ sinβ H^{0,1}_{n-1}
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
                #c̄ₙₘ = 0
                #d̄ₙₘ = 0
                ēₙₘ = inv2np1 * √T((2n+2)*(2n+1))
                Hₙ₊₁⁰[1+m] = sinβ * ēₙₘ * Hₙ[iₙ]
            end
            # At the end of this step, iₙ is pointing at H^{0, n}_{n}

            if m′ₘₐₓ > 0
                # Step 3: Compute H^{1,m}_{n}(β) for m=1,...,n
                iₙ = iₙ + m′offset₊(WC, n, 0, n) - n + 1
                # iₙ is now pointing at H^{1, 1}_{n}
                invsqrtnnp1 = inv(sqrtnnp1)
                # println("pre")
                # println((iₙ))
                # println((n, 1, 1))
                # println(indices[iₙ])
                # println()
                for m in 1:n
                    a = √T((n+m+1)*(n-m+1))
                    b1 = √T((n+m+1)*(n+m+2))
                    b2 = √T((n-m+1)*(n-m+2))
                    Hₙ[iₙ] = -invsqrtnnp1 * (
                        b1 * cosβ₋ * Hₙ₊₁⁰[m+2]
                        + b2 * cosβ₊ * Hₙ₊₁⁰[m]
                        + a * sinβ * Hₙ₊₁⁰[m+1]
                    )
                    iₙ += 1
                    # println("0")
                    # println((iₙ))
                    # println((n, 1, m+1))
                    # if iₙ > size(indices, 1)
                    #     println("After ", indices[iₙ-1])
                    # else
                    #     println(indices[iₙ])
                    # end
                    # println()
                end
                iₙ -= 1
                # Now, iₙ is pointing at H^{1, n}_{n}

                # Step 4: Compute H^{m′+1, m}_{n}(β) for m′=1,...,n−1, m=m′,...,n
                d1 = sqrtnnp1
                # println("preA")
                # println((iₙ))
                # println((n, 1, n))
                # println(indices[iₙ])
                # println()
                iₙ′′ = iₙ - m′offset₋(WC, n, 1, n) - n + 2  # (n, 0, 2)
                iₙ′ = iₙ - n + 1  # (n, 1, 1)
                iₙ = iₙ + m′offset₊(WC, n, 1, n) - n + 2 # (n, 2, 2)
                for m′ in 1:min(n, m′ₘₐₓ)-1
                    #iₙ′ = iₙ - (n-m′)

                    # iₙ = offset(WC, n, m′+1, m′+1)
                    # iₙ′ = offset(WC, n, m′, m′)
                    # iₙ′′ = offset(WC, n, m′-1, m′+1)

                    # iₙ points at H^{m′+1, m}_{n}
                    # iₙ′ points at H^{m′, m-1}_{n}
                    # iₙ′′ points at H^{m′-1, m}_{n}

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
                        # println("A")
                        # println((n, m′, m))
                        # println((iₙ, iₙ′, iₙ′+2, iₙ′′))
                        # println("L: ", (n, m′+1, m), " ", indices[iₙ])
                        # println("R1: ", (n, m′, m-1), " ", indices[iₙ′])
                        # println("R2: ", (n, m′, m+1), " ", indices[iₙ′+2])
                        # println("R3: ", (n, m′-1, m), " ", indices[iₙ′′])
                        # println()
                        Hₙ[iₙ] = invd1 * (
                            d2 * Hₙ[iₙ′′]
                            - d3 * Hₙ[iₙ′]
                            + d4 * Hₙ[iₙ′+2]
                        )
                        iₙ += 1
                        iₙ′ += 1
                        iₙ′′ += 1
                    end
                    let m = n
                        d3 = sqrt2n
                        # println("B")
                        # println((n, m′, m))
                        # println((iₙ, iₙ′, iₙ′′))
                        # println("L: ", (n, m′+1, m), " ", indices[iₙ])
                        # println("R1: ", (n, m′, m-1), " ", indices[iₙ′])
                        # println("R3: ", (n, m′-1, m), " ", indices[iₙ′′])
                        # println()
                        Hₙ[iₙ] = invd1 * (
                            d2 * Hₙ[iₙ′′]
                            - d3 * Hₙ[iₙ′]
                        )
                        # iₙ += 1
                        # iₙ′ += 1
                    end
                    iₙ += m′offset₊(WC, n, m′+1, n) - n + m′ + 2 # → (n, m′+2, m′+2)
                    iₙ′ += m′offset₊(WC, n, m′, n-1) - (n-1) + m′ + 1 # → (n, m′+1, m′+1)
                    iₙ′′ += m′offset₊(WC, n, m′-1, n) - n + m′ + 2 # → (n, m′, m′+2)
                end  # Step 4

                # Step 5: Compute H^{m′−1, m}_{n}(β) for m′=0,...,−n+1, m=−m′,...,n
                iₙ′′ = iₙ0 + m′offset₊(WC, n, 0, 0) + n  # (n, 1, n)
                iₙ′ = iₙ0 + n - 1  # (n, 0, n-1)
                iₙ = iₙ0 - m′offset₋(WC, n, 0, 0) + n  # (n, -1, n)
                d2 = sqrtnnp1
                #offset2 = 0  # Account for 2 missing elements when m′==0
                for m′ in 0:-1:-min(n, m′ₘₐₓ)+1
                    # iₙ′ = iₙ + (n+m′) + 1
                    # # iₙ points at H^{m′-1, m}_{n}
                    # # iₙ′ points at H^{m′, m}_{n}

                    # iₙ = offset(WC, n, m′-1, n)
                    # iₙ′ = offset(WC, n, m′, n-1)
                    # iₙ′′ = offset(WC, n, m′+1, n)

                    # iₙ points at H^{m′-1, m}_{n}
                    # iₙ′ points at H^{m′, m-1}_{n}
                    # iₙ′′ points at H^{m′+1, m}_{n}

                    # d1 ≔ d^{m′-1}_{n} = -√T((n-m′+1)*(n+m′))
                    # d2 ≔ d^{m′}_{n} = √T((n-m′)*(n+m′+1))
                    # d3 ≔ d^{m-1}_{n} = √T((n-m+1)*(n+m))
                    # d4 ≔ d^{m}_{n} = √T((n-m)*(n+m+1))

                    d1 = -√T((n-m′+1)*(n+m′))
                    invd1 = inv(d1)
                    d3 = sqrt2n
                    let m = n
                        # println("A")
                        # println((n, m′, m))
                        # println([iₙ, iₙ′, iₙ′′])
                        # println("L: ", (n, m′-1, m), " ", indices[iₙ])
                        # println("R1: ", (n, m′, m-1), " ", indices[iₙ′])
                        # println("R3: ", (n, m′+1, m), " ", indices[iₙ′′])
                        # println()
                        Hₙ[iₙ] = invd1 * (
                            d2 * Hₙ[iₙ′′]
                            + d3 * Hₙ[iₙ′]
                        )
                        iₙ -= 1
                        iₙ′ -= 1
                        iₙ′′ -= 1
                    end
                    for m in n-1:-1:1-m′
                        d4 = d3
                        d3 = √T((n-m+1)*(n+m))
                        # println("B")
                        # println((n, m′, m))
                        # println([iₙ, iₙ′, iₙ′′])
                        # println("L: ", (n, m′-1, m), " ", indices[iₙ])
                        # println("R1: ", (n, m′, m-1), " ", indices[iₙ′])
                        # println("R2: ", (n, m′, m+1), " ", indices[iₙ′+2])
                        # println("R3: ", (n, m′+1, m), " ", indices[iₙ′′])
                        # println()
                        Hₙ[iₙ] = invd1 * (
                            d2 * Hₙ[iₙ′′]
                            + d3 * Hₙ[iₙ′]
                            - d4 * Hₙ[iₙ′+2]
                        )
                        iₙ -= 1
                        iₙ′ -= 1
                        iₙ′′ -= 1
                    end
                    d2 = d1
                    # iₙ points at H^{m′-1, -m′}_{n}
                    # iₙ′ points at H^{m′, -m′-1}_{n}
                    # iₙ′′ points at H^{m′+1, -m′}_{n}

                    # iₙ points at H^{m′-1, n}_{n}
                    # iₙ′ points at H^{m′, n-1}_{n}
                    # iₙ′′ points at H^{m′+1, n}_{n}

                    iₙ += -m′offset₋(WC, n, m′-1, 1-m′) + n + m′ # → (n, m′-2, n)
                    iₙ′ += -m′offset₋(WC, n, m′, -m′) + (n-1) + m′ + 1 # → (n, m′-1, n-1)
                    iₙ′′ += -m′offset₋(WC, n, m′+1, 1-m′) + n + m′ # → (n, m′, n)

                end  # Step 5

            end


        end  # if nₘₐₓ > 0
    end  # @inbounds
    Hₙ
end
