function abd(ℓₘₐₓ, T)
    a = [√((n+1+m) * (n+1-m) / T((2n+1)*(2n+3))) for n in 0:ℓₘₐₓ+1 for m in 0:n]
    b = [(m<0 ? -1 : 1) * √((n-m-1) * (n-m) / T((2n-1)*(2n+1))) for n in 0:ℓₘₐₓ+1 for m in -n:n]
    d = [(m<0 ? -1 : 1) * (√T((n-m) * (n+m+1))) / 2 for n in 0:ℓₘₐₓ+1 for m in -n:n]
    (a, b, d)
end


function H2!(
    Hwedge::AbstractVector{T}, Hextra::AbstractVector{T}, expiβ::Complex{T},
    ℓₘₐₓ, m′ₘₐₓ, (a,b,d), Hindex=WignerHindex
) where {T<:Real}
    @assert size(Hwedge) == (Hindex(ℓₘₐₓ, ℓₘₐₓ, ℓₘₐₓ, m′ₘₐₓ),)

    sqrt3 = √T(3)
    cosβ = expiβ.re
    sinβ = expiβ.im
    cosβ₊ = (1+cosβ)/2
    cosβ₋ = (1-cosβ)/2

    @inbounds begin
        # Step 1: If n=0 set H_{0}^{0,0}=1
        Hwedge[1] = 1

        if ℓₘₐₓ > 0

            if m′ₘₐₓ ≥ 0
                # Step 2: Compute H^{0,m}_{n}(β) for m=0,...,n and H^{0,m}_{n+1}(β) for m=0,...,n+1
                nₘₐₓstep2 = m′ₘₐₓ>0 ? ℓₘₐₓ+1 : ℓₘₐₓ
                # n = 1
                n0n_index = Hindex(1, 0, 1, m′ₘₐₓ)
                Hwedge[n0n_index] = sqrt3 * sinβ
                Hwedge[n0n_index-1] = sqrt3 * cosβ
                # n = 2, ..., ℓₘₐₓ, ℓₘₐₓ+1?
                for n in 2:nₘₐₓstep2
                    if n <= ℓₘₐₓ
                        n0n_index = Hindex(n, 0, n, m′ₘₐₓ)
                        H = Hwedge
                    else  # n == ℓₘₐₓ+1
                        n0n_index = n + 1
                        H = Hextra
                    end
                    nm10nm1_index = Hindex(n-1, 0, n-1, m′ₘₐₓ)
                    inv2n = inv(T(2n))
                    aₙ = √(T(2n+1)/(2n-1))
                    bₙ = aₙ * √(T(2n-2)/n)

                    # m = n
                    eₙₘ = inv2n * √T(2n*(2n+1))
                    H[n0n_index] = sinβ * eₙₘ * Hwedge[nm10nm1_index]

                    # m = n-1
                    eₙₘ = inv2n * √T((2n-2)*(2n+1))
                    cₙₘ = 2inv2n * √T(2n+1)
                    H[n0n_index-1] = (
                        cosβ * cₙₘ * Hwedge[nm10nm1_index]
                        + sinβ * eₙₘ * Hwedge[nm10nm1_index-1]
                    )

                    # m = n-2, ..., 2
                    for i in 2:n-2
                        # m = n-i
                        cₙₘ = 2inv2n * aₙ * √T((2n-i)*i)
                        dₙₘ = inv2n * aₙ * √T(i*(i-1))
                        eₙₘ = inv2n * aₙ * √T((2n-i)*(2n-i-1))
                        H[n0n_index-i] = (
                            cosβ * cₙₘ * Hwedge[nm10nm1_index-i+1]
                            - sinβ * (
                                dₙₘ * Hwedge[nm10nm1_index-i+2]
                                - eₙₘ * Hwedge[nm10nm1_index-i]
                            )
                        )
                    end

                    # m = 1
                    cₙₘ = 2inv2n * aₙ * √T((n+1)*(n-1))
                    dₙₘ = inv2n * aₙ * √T((n-1)*(n-2))
                    eₙₘ = inv2n * aₙ * √T(2n*(n+1))
                    H[n0n_index-n+1] = (
                        cosβ * cₙₘ * Hwedge[nm10nm1_index-n+2]
                        - sinβ * (
                            dₙₘ * Hwedge[nm10nm1_index-n+3]
                            - eₙₘ * Hwedge[nm10nm1_index-n+1]
                        )
                    )

                    # m = 0
                    cₙₘ = aₙ
                    dₙₘ = inv2n * aₙ * √T(n*(n-1))
                    eₙₘ = dₙₘ
                    H[n0n_index-n] = (
                        aₙ * cosβ * Hwedge[nm10nm1_index-n+1]
                        - bₙ * sinβ * Hwedge[nm10nm1_index-n+2] / 2
                    )

                end  # Step 2

                # Normalize, changing P̄ to H values
                for n in 1:nₘₐₓstep2
                    if n <= ℓₘₐₓ
                        n00_index = Hindex(n, 0, 0, m′ₘₐₓ)
                        H = Hwedge
                    else  # n == ℓₘₐₓ+1
                        n00_index = 1
                        H = Hextra
                    end
                    const0 = inv(√T(2n+1))
                    const1 = inv(√T(4n+2))
                    H[n00_index] *= const0
                    for m in 1:n
                        H[n00_index+m] *= const1
                    end
                end
            end  # if m′ₘₐₓ ≥ 0

            if m′ₘₐₓ > 0

                # Step 3: Compute H^{1,m}_{n}(β) for m=1,...,n
                for n in 1:ℓₘₐₓ
                    # m = 1, ..., n
                    i1 = Hindex(n, 1, 1, m′ₘₐₓ)
                    if n+1 <= ℓₘₐₓ
                        i2 = Hindex(n+1, 0, 0, m′ₘₐₓ)
                        H2 = Hwedge
                    else  # n+1 == ℓₘₐₓ+1
                        i2 = 1
                        H2 = Hextra
                    end
                    i3 = nm_index(n+1, 0)
                    i4 = nabsm_index(n, 1)
                    inverse_b5 = inv(b[i3])
                    for i in 0:n-1
                        b6 = b[-i+i3-2]
                        b7 = b[i+i3]
                        a8 = a[i+i4]
                        Hwedge[i+i1] = inverse_b5 * (
                            b6 * cosβ₋ * H2[i+i2+2]
                            - b7 * cosβ₊ * H2[i+i2]
                            - a8 * sinβ * H2[i+i2+1]
                        )
                    end
                end  # Step 3

                # Step 4: Compute H^{m'+1, m}_{n}(β) for m'=1,...,n−1, m=m',...,n
                for n in 2:ℓₘₐₓ
                    for mp in 1:min(n, m′ₘₐₓ)-1
                        # m = m', ..., n-1
                        i1 = Hindex(n, mp+1, mp+1, m′ₘₐₓ) - 1
                        i2 = Hindex(n, mp-1, mp, m′ₘₐₓ)
                        i3 = Hindex(n, mp, mp, m′ₘₐₓ) - 1
                        i4 = Hindex(n, mp, mp+1, m′ₘₐₓ)
                        i5 = nm_index(n, mp)
                        i6 = nm_index(n, mp-1)
                        inverse_d5 = inv(d[i5])
                        d6 = d[i6]
                        d7 = d6
                        d8 = d[i5]
                        for i in 1:n-mp-1
                            d7 = d8
                            d8 = d[i+i5]
                            Hwedge[i+i1] = inverse_d5 * (
                                d6 * Hwedge[i+i2]
                                - d7 * Hwedge[i+i3]
                                + d8 * Hwedge[i+i4]
                            )
                        end
                        # m = n
                        let i=n-mp
                            d7 = d8
                            Hwedge[i+i1] = inverse_d5 * (
                                d6 * Hwedge[i+i2]
                                - d7 * Hwedge[i+i3]
                            )
                        end
                    end
                end  # Step 4

                # Step 5: Compute H^{m'−1, m}_{n}(β) for m'=−1,...,−n+1, m=−m',...,n
                for n in 0:ℓₘₐₓ
                    for mp in 0:-1:1-min(n, m′ₘₐₓ)
                        # m = -m', ..., n-1
                        i1 = Hindex(n, mp-1, -mp+1, m′ₘₐₓ) - 1
                        i2 = Hindex(n, mp+1, -mp+1, m′ₘₐₓ) - 1
                        i3 = Hindex(n, mp, -mp, m′ₘₐₓ) - 1
                        i4 = Hindex(n, mp, -mp+1, m′ₘₐₓ)
                        i5 = nm_index(n, mp-1)
                        i6 = nm_index(n, mp)
                        i7 = nm_index(n, -mp-1)
                        i8 = nm_index(n, -mp)
                        inverse_d5 = inv(d[i5])
                        d6 = d[i6]
                        d7 = d[i7]
                        d8 = d[i8]
                        for i in 1:n+mp-1
                            d7 = d8
                            d8 = d[i+i8]
                            Hwedge[i+i1] = inverse_d5 * (
                                d6 * Hwedge[i+i2]
                                + d7 * Hwedge[i+i3]
                                - d8 * Hwedge[i+i4]
                            )
                        end
                        # m = n
                        let i=n+mp
                            d7 = d8
                            Hwedge[i+i1] = inverse_d5 * (
                                d6 * Hwedge[i+i2]
                                + d7 * Hwedge[i+i3]
                            )
                        end
                    end
                end  # Step 5

            end  # if m′ₘₐₓ > 0
        end  # if ℓₘₐₓ > 0
    end  # @inbounds
    Hwedge
end
