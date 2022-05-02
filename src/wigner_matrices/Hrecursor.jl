"""Return flat index into arrray of [n, m] pairs.

Assumes array is ordered as

    [
        [n, m]
        for n in range(n_max+1)
        for m in range(-n, n+1)
    ]

"""
nm_index(n, m) = m + n * (n + 1) + 1


"""Return flat index into arrray of [n, abs(m)] pairs

Assumes array is ordered as

    [
        [n, m]
        for n in range(n_max+1)
        for m in range(n+1)
    ]

"""
nabsm_index(n, absm) = absm + (n * (n + 1)) ÷ 2 + 1


"""Return flat index into arrray of [n, mp, m]

Assumes array is ordered as

    [
        [n, mp, m]
        for n in range(n_max+1)
        for mp in range(-n, n+1)
        for m in range(-n, n+1)
    ]

"""
nmpm_index(n, mp, m) = (((4n + 6) * n + 6mp + 5) * n + 3(m + mp)) ÷ 3 + 1


"""
    abd(ℓₘₐₓ, T)

Pre-compute constants used in Wigner H recurrence.

"""
function abd(ℓₘₐₓ, T)
    a = T[√((n+1+m)*(n+1-m)/T((2n+1)*(2n+3))) for n in 0:ℓₘₐₓ+1 for m in 0:n]
    b = T[(m<0 ? -1 : 1) * √((n-m-1)*(n-m)/T((2n-1)*(2n+1))) for n in 0:ℓₘₐₓ+1 for m in -n:n]
    d = T[(m<0 ? -1 : 1) * (√T((n-m)*(n+m+1))) / 2 for n in 0:ℓₘₐₓ+1 for m in -n:n]
    (a, b, d)
end


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


function H!(
    H::AbstractVector, expiβ::Complex{T}, ℓₘₐₓ, m′ₘₐₓ, (a,b,d), Hindex=WignerHindex
) where {T<:Real}
    m′ₘₐₓ = abs(m′ₘₐₓ)
    @assert m′ₘₐₓ ≤ ℓₘₐₓ
    @assert size(H) == (Hindex(ℓₘₐₓ, m′ₘₐₓ, ℓₘₐₓ; m′ₘₐₓ=m′ₘₐₓ),)

    invsqrt2 = inv(√T(2))
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
    HΩ, HΨ = zero(eltype(H)), zero(eltype(H))

    @inbounds begin
        # Step 1: If n=0 set H_{0}^{0,0}=1
        H[1] = 1

        if ℓₘₐₓ > 0

            begin
                # Step 2: Compute H^{0,m}_{n}(β) for m=0,...,n and n=1,...,ℓₘₐₓ,ℓₘₐₓ+1?
                nₘₐₓstep2 = m′ₘₐₓ>0 ? ℓₘₐₓ+1 : ℓₘₐₓ
                # n = 1
                n0n_index = Hindex(1, 0, 0; m′ₘₐₓ=m′ₘₐₓ)
                H[n0n_index] = cosβ
                H[n0n_index+1] = invsqrt2 * sinβ
                # n = 2, ..., ℓₘₐₓ, ℓₘₐₓ+1?
                for n in 2:nₘₐₓstep2
                    superjacent = n > ℓₘₐₓ
                    if superjacent
                        n00_index = Hindex(n-1, -1, 1; m′ₘₐₓ=m′ₘₐₓ)
                    else
                        n00_index = Hindex(n, 0, 0; m′ₘₐₓ=m′ₘₐₓ)
                    end
                    nm100_index = Hindex(n-1, 0, 0; m′ₘₐₓ=m′ₘₐₓ)
                    inv2n = inv(T(2n))
                    b̄ₙ = √(T(n-1)/n)

                    # H^{0,0}_{n} = cosβ H^{0,0}_{n-1} - b̄ₙ sinβ H^{0,1}_{n-1}
                    H[n00_index] = cosβ * H[nm100_index] - b̄ₙ  * sinβ * H[nm100_index+1]

                    # H^{0,m}_{n} = c̄ₙₘ cosβ H^{0,m}_{n-1} - sinβ [ d̄ₙₘ H^{0,m+1}_{n-1} - ēₙₘ H^{0,m-1}_{n-1} ]
                    for m in 1:n-2
                        c̄ₙₘ = 2inv2n * √T((n+m)*(n-m))
                        d̄ₙₘ = inv2n * √T((n-m)*(n-m-1))
                        ēₙₘ = inv2n * √T((n+m)*(n+m-1))
                        H[n00_index+m] = (
                            c̄ₙₘ * cosβ * H[nm100_index+m]
                            - sinβ * (
                                d̄ₙₘ * H[nm100_index+m+1]
                                - ēₙₘ * H[nm100_index+m-1]
                            )
                        )
                    end
                    let m = n-1
                        c̄ₙₘ = 2inv2n * √T(2n-1)
                        d̄ₙₘ = 0
                        ēₙₘ = inv2n * √T((2n-1)*(2n-2))
                        if superjacent
                            HΨ = (
                                c̄ₙₘ * cosβ * H[nm100_index+m]
                                + sinβ * ēₙₘ * H[nm100_index+m-1]
                            )
                        else
                            H[n00_index+m] = (
                                c̄ₙₘ * cosβ * H[nm100_index+m]
                                + sinβ * ēₙₘ * H[nm100_index+m-1]
                            )
                        end
                    end
                    let m = n
                        c̄ₙₘ = 0
                        d̄ₙₘ = 0
                        ēₙₘ = inv2n * √T(2n*(2n-1))
                        if superjacent
                            HΩ = sinβ * ēₙₘ * H[nm100_index+m-1]
                        else
                            H[n00_index+m] = sinβ * ēₙₘ * H[nm100_index+m-1]
                        end
                    end

                end
            end

            if m′ₘₐₓ > 0

                # Step 3: Compute H^{1,m}_{n}(β) for m=1,...,n
                for n in 1:ℓₘₐₓ
                    # m = 1, ..., n
                    i1 = Hindex(n, 1, 1; m′ₘₐₓ=m′ₘₐₓ)
                    if n+1 <= ℓₘₐₓ
                        i2 = Hindex(n+1, 0, 0; m′ₘₐₓ=m′ₘₐₓ)
                        nₘₐₓstep3 = n-1
                    else  # n+1 == ℓₘₐₓ+1  &&  m′ₘₐₓ > 0
                        i2 = Hindex(n, -1, 1; m′ₘₐₓ=m′ₘₐₓ)
                        nₘₐₓstep3 = n-3
                    end
                    i3 = nm_index(n+1, 0)
                    i4 = nabsm_index(n, 1)
                    inverse_b5 = inv(b[i3])
                    for i in 0:nₘₐₓstep3
                        b6 = b[-i+i3-2]
                        b7 = b[i+i3]
                        a8 = a[i+i4]
                        H[i+i1] = inverse_b5 * (
                            b6 * cosβ₋ * H[i+i2+2]
                            - b7 * cosβ₊ * H[i+i2]
                            - a8 * sinβ * H[i+i2+1]
                        )
                    end
                    if n == ℓₘₐₓ  # &&  m′ₘₐₓ > 0
                        if n>1
                            let i = n-2
                                b6 = b[-i+i3-2]
                                b7 = b[i+i3]
                                a8 = a[i+i4]
                                H[i+i1] = inverse_b5 * (
                                    b6 * cosβ₋ * HΨ
                                    - b7 * cosβ₊ * H[i+i2]
                                    - a8 * sinβ * H[i+i2+1]
                                )
                            end
                        end
                        let i = n-1
                            b6 = b[-i+i3-2]
                            b7 = b[i+i3]
                            a8 = a[i+i4]
                            H[i+i1] = inverse_b5 * (
                                b6 * cosβ₋ * HΩ
                                - b7 * cosβ₊ * H[i+i2]
                                - a8 * sinβ * HΨ
                            )
                        end
                    end
                end  # Step 3

                # Step 4: Compute H^{m'+1, m}_{n}(β) for m'=1,...,n−1, m=m',...,n
                for n in 2:ℓₘₐₓ
                    for mp in 1:min(n, m′ₘₐₓ)-1
                        # m = m', ..., n-1
                        i1 = Hindex(n, mp+1, mp+1; m′ₘₐₓ=m′ₘₐₓ) - 1
                        i2 = Hindex(n, mp-1, mp; m′ₘₐₓ=m′ₘₐₓ)
                        i3 = Hindex(n, mp, mp; m′ₘₐₓ=m′ₘₐₓ) - 1
                        i4 = Hindex(n, mp, mp+1; m′ₘₐₓ=m′ₘₐₓ)
                        i5 = nm_index(n, mp)
                        i6 = nm_index(n, mp-1)
                        inverse_d5 = inv(d[i5])
                        d6 = d[i6]

                        d7 = d6
                        d8 = d[i5]
                        for i in 1:n-mp-1
                            d7 = d8
                            d8 = d[i+i5]
                            H[i+i1] = inverse_d5 * (
                                d6 * H[i+i2]
                                - d7 * H[i+i3]
                                + d8 * H[i+i4]
                            )
                        end

                        # m = n
                        let i=n-mp
                            d7 = d8
                            H[i+i1] = inverse_d5 * (
                                d6 * H[i+i2]
                                - d7 * H[i+i3]
                            )
                        end
                    end
                end  # Step 4

                # Step 5: Compute H^{m'−1, m}_{n}(β) for m'=−1,...,−n+1, m=−m',...,n
                for n in 0:ℓₘₐₓ
                    for mp in 0:-1:1-min(n, m′ₘₐₓ)
                        # m = -m', ..., n-1
                        i1 = Hindex(n, mp-1, -mp+1; m′ₘₐₓ=m′ₘₐₓ) - 1
                        i2 = Hindex(n, mp+1, -mp+1; m′ₘₐₓ=m′ₘₐₓ) - 1
                        i3 = Hindex(n, mp, -mp; m′ₘₐₓ=m′ₘₐₓ) - 1
                        i4 = Hindex(n, mp, -mp+1; m′ₘₐₓ=m′ₘₐₓ)
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
                            H[i+i1] = inverse_d5 * (
                                d6 * H[i+i2]
                                + d7 * H[i+i3]
                                - d8 * H[i+i4]
                            )
                        end
                        # m = n
                        let i=n+mp
                            d7 = d8
                            H[i+i1] = inverse_d5 * (
                                d6 * H[i+i2]
                                + d7 * H[i+i3]
                            )
                        end
                    end
                end  # Step 5

            end  # if m′ₘₐₓ > 0
        end  # if ℓₘₐₓ > 0
    end  # @inbounds
    H
end
