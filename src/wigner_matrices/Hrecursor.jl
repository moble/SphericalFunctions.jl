function abd(ℓₘₐₓ, T)
    a = [√((n+1+m)*(n+1-m)/T((2n+1)*(2n+3))) for n in 0:ℓₘₐₓ+1 for m in 0:n]
    b = [(m<0 ? -1 : 1) * √((n-m-1)*(n-m)/T((2n-1)*(2n+1))) for n in 0:ℓₘₐₓ+1 for m in -n:n]
    d = [(m<0 ? -1 : 1) * (√T((n-m)*(n+m+1))) / 2 for n in 0:ℓₘₐₓ+1 for m in -n:n]
    (a, b, d)
end


function H2!(
    Hwedge::AbstractVector{T}, expiβ::Complex{T}, ℓₘₐₓ::Integer, m′ₘₐₓ::Integer,
    (a,b,d), Hindex=WignerHindex
) where {T<:Real}
    m′ₘₐₓ = abs(m′ₘₐₓ)
    @assert m′ₘₐₓ ≤ ℓₘₐₓ
    @assert size(Hwedge) == (Hindex(ℓₘₐₓ, m′ₘₐₓ, ℓₘₐₓ, m′ₘₐₓ),)

    sqrt3 = √T(3)
    cosβ = expiβ.re
    sinβ = expiβ.im
    cosβ₊ = (1+cosβ)/2
    cosβ₋ = (1-cosβ)/2

    # Note: In step 3, we set the H^{1, m}_{ℓₘₐₓ} terms from H^{0, m+1}_{ℓₘₐₓ+1} data,
    # which of course should not appear in the final result, because it involves ℓ>ℓₘₐₓ.
    # This requires storing those extra terms somewhere temporarily during the execution
    # of this function.  Fortunately, this is only required if m′ₘₐₓ>0, which also implies
    # that there are at least ℓₘₐₓ slots in the input Hwedge array that are not needed until
    # step 5.  We use those slots for most of the temporary data.  However, there are still
    # two slots needed, which we allocate as individual variables, representing the last
    # and second-to-last elements H^{0, ℓₘₐₓ+1}_{ℓₘₐₓ+1} and H^{0, ℓₘₐₓ}_{ℓₘₐₓ+1}:
    HΩ, HΨ = zero(T), zero(T)

    Hrange = SphericalFunctions.WignerHrange(m′ₘₐₓ, ℓₘₐₓ)

    @inbounds begin
        # Step 1: If n=0 set H_{0}^{0,0}=1
        Hwedge[1] = 1

        if ℓₘₐₓ > 0

            begin
                # Step 2: Compute H^{0,m}_{n}(β) for m=0,...,n and n=1,...,ℓₘₐₓ,ℓₘₐₓ+1?
                nₘₐₓstep2 = m′ₘₐₓ>0 ? ℓₘₐₓ+1 : ℓₘₐₓ
                # n = 1
                n0n_index = Hindex(1, 0, 1, m′ₘₐₓ)
                Hwedge[n0n_index] = sqrt3 * sinβ
                Hwedge[n0n_index-1] = sqrt3 * cosβ
                # n = 2, ..., ℓₘₐₓ, ℓₘₐₓ+1?
                for n in 2:nₘₐₓstep2
                    superjacent = n > ℓₘₐₓ
                    if superjacent
                        n0n_index = Hindex(n-1, -1, n-1, m′ₘₐₓ)+2
                    else
                        n0n_index = Hindex(n, 0, n, m′ₘₐₓ)
                    end
                    nm10nm1_index = Hindex(n-1, 0, n-1, m′ₘₐₓ)
                    inv2n = inv(T(2n))
                    aₙ = √(T(2n+1)/(2n-1))
                    bₙ = aₙ * √(T(2n-2)/n)

                    # m = n
                    eₙₘ = inv2n * √T(2n*(2n+1))
                    if superjacent
                        HΩ = sinβ * eₙₘ * Hwedge[nm10nm1_index]
                    else
                        Hwedge[n0n_index] = sinβ * eₙₘ * Hwedge[nm10nm1_index]
                    end

                    # m = n-1
                    eₙₘ = inv2n * √T((2n-2)*(2n+1))
                    cₙₘ = 2inv2n * √T(2n+1)
                    if superjacent
                        HΨ = (
                            cosβ * cₙₘ * Hwedge[nm10nm1_index]
                            + sinβ * eₙₘ * Hwedge[nm10nm1_index-1]
                        )
                    else
                        Hwedge[n0n_index-1] = (
                            cosβ * cₙₘ * Hwedge[nm10nm1_index]
                            + sinβ * eₙₘ * Hwedge[nm10nm1_index-1]
                        )
                    end

                    # m = n-2, ..., 2
                    for i in 2:n-2
                        # m = n-i
                        cₙₘ = 2inv2n * aₙ * √T((2n-i)*i)
                        dₙₘ = inv2n * aₙ * √T(i*(i-1))
                        eₙₘ = inv2n * aₙ * √T((2n-i)*(2n-i-1))
                        Hwedge[n0n_index-i] = (
                            cosβ * cₙₘ * Hwedge[nm10nm1_index-i+1]
                            - sinβ * (
                                dₙₘ * Hwedge[nm10nm1_index-i+2]
                                - eₙₘ * Hwedge[nm10nm1_index-i]
                            )
                        )
                    end

                    cₙₘ = 2inv2n * aₙ * √T((n+1)*(n-1))
                    dₙₘ = inv2n * aₙ * √T((n-1)*(n-2))
                    eₙₘ = inv2n * aₙ * √T(2n*(n+1))
                    if superjacent && ℓₘₐₓ == 1
                        # This rare case repeats the calculation of HΨ, but with a different
                        # value of eₙₘ.  It only occurs when ℓₘₐₓ == 1, so there's no point
                        # trying to be too clever about avoiding the previous calculation
                        HΨ = (
                            cosβ * cₙₘ * Hwedge[nm10nm1_index-n+2]
                            + sinβ * eₙₘ * Hwedge[nm10nm1_index-n+1]
                        )
                    elseif n == 2
                        # This avoids accessing un-initialized data at
                        # `Hwedge[nm10nm1_index-n+3]`, even though dₙₘ==0 in this case.
                        # We could remove this if case, and just use the case below
                        # without any negative consequences, if we knew that element
                        # did not contain NaN or Inf.
                        Hwedge[n0n_index-n+1] = (
                            cosβ * cₙₘ * Hwedge[nm10nm1_index-n+2]
                            + sinβ * eₙₘ * Hwedge[nm10nm1_index-n+1]
                        )
                    else
                        Hwedge[n0n_index-n+1] = (
                            cosβ * cₙₘ * Hwedge[nm10nm1_index-n+2]
                            - sinβ * (
                                dₙₘ * Hwedge[nm10nm1_index-n+3]
                                - eₙₘ * Hwedge[nm10nm1_index-n+1]
                            )
                        )
                    end

                    # m = 0
                    cₙₘ = aₙ
                    dₙₘ = inv2n * aₙ * √T(n*(n-1))
                    eₙₘ = dₙₘ
                    Hwedge[n0n_index-n] = (
                        aₙ * cosβ * Hwedge[nm10nm1_index-n+1]
                        - bₙ * sinβ * Hwedge[nm10nm1_index-n+2] / 2
                    )

                end  # Step 2 recurrence

                # Normalize, changing P̄ to H values
                for n in 1:nₘₐₓstep2
                    const0 = inv(√T(2n+1))
                    const1 = inv(√T(4n+2))
                    if n <= ℓₘₐₓ
                        n00_index = Hindex(n, 0, 0, m′ₘₐₓ)
                        Hwedge[n00_index] *= const0
                        for m in 1:n
                            Hwedge[n00_index+m] *= const1
                        end
                    else  # n == ℓₘₐₓ+1  &&  m′ₘₐₓ > 0
                        n00_index = Hindex(n-1, -1, 1, m′ₘₐₓ)
                        Hwedge[n00_index] *= const0
                        for m in 1:n-2
                            Hwedge[n00_index+m] *= const1
                        end
                        HΨ *= const1
                        HΩ *= const1
                    end
                end  # Step 2 normalization

            end  # if m′ₘₐₓ ≥ 0

            if m′ₘₐₓ > 0

                # Step 3: Compute H^{1,m}_{n}(β) for m=1,...,n
                for n in 1:ℓₘₐₓ
                    # m = 1, ..., n
                    i1 = Hindex(n, 1, 1, m′ₘₐₓ)
                    if n+1 <= ℓₘₐₓ
                        i2 = Hindex(n+1, 0, 0, m′ₘₐₓ)
                        nₘₐₓstep3 = n-1
                    else  # n+1 == ℓₘₐₓ+1  &&  m′ₘₐₓ > 0
                        i2 = Hindex(n, -1, 1, m′ₘₐₓ)
                        nₘₐₓstep3 = n-3
                    end
                    i3 = nm_index(n+1, 0)
                    i4 = nabsm_index(n, 1)
                    inverse_b5 = inv(b[i3])
                    for i in 0:nₘₐₓstep3
                        b6 = b[-i+i3-2]
                        b7 = b[i+i3]
                        a8 = a[i+i4]
                        Hwedge[i+i1] = inverse_b5 * (
                            b6 * cosβ₋ * Hwedge[i+i2+2]
                            - b7 * cosβ₊ * Hwedge[i+i2]
                            - a8 * sinβ * Hwedge[i+i2+1]
                        )
                    end
                    if n == ℓₘₐₓ  # &&  m′ₘₐₓ > 0
                        if n>1
                            let i = n-2
                                b6 = b[-i+i3-2]
                                b7 = b[i+i3]
                                a8 = a[i+i4]
                                Hwedge[i+i1] = inverse_b5 * (
                                    b6 * cosβ₋ * HΨ
                                    - b7 * cosβ₊ * Hwedge[i+i2]
                                    - a8 * sinβ * Hwedge[i+i2+1]
                                )
                            end
                        end
                        let i = n-1
                            b6 = b[-i+i3-2]
                            b7 = b[i+i3]
                            a8 = a[i+i4]
                            Hwedge[i+i1] = inverse_b5 * (
                                b6 * cosβ₋ * HΩ
                                - b7 * cosβ₊ * Hwedge[i+i2]
                                - a8 * sinβ * HΨ
                            )
                        end
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
