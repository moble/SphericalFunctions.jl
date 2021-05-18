"""Algorithm for computing H, as given by arxiv:1403.7698

H is related to Wigner's (small) d via

    dₗⁿᵐ = ϵₙ ϵ₋ₘ Hₗⁿᵐ,

where

         ⎧ 1 for k≤0
    ϵₖ = ⎨
         ⎩ (-1)ᵏ for k>0

H has various advantages over d, including the fact that it can be efficiently
and robustly valculated via recurrence relations, and the following symmetry
relations:

    H^{m', m}_n(β) = H^{m, m'}_n(β)
    H^{m', m}_n(β) = H^{-m', -m}_n(β)
    H^{m', m}_n(β) = (-1)^{n+m+m'} H^{-m', m}_n(π - β)
    H^{m', m}_n(β) = (-1)^{m+m'} H^{m', m}_n(-β)

Because of these symmetries, we only need to evaluate at most 1/4 of all the
elements.

"""

const sqrt3 = √3
const inverse_sqrt2 = inv(√2)


ϵ(m) = (m <= 0 ? 1 : (isodd(m) ? -1 : 1))


"""Return sign of input, with sign(0)=1"""
sign(m) = (m < 0 ? -1 : 1)


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


@inbounds function _step_1!(w::WignerWorkspace)
    """If n=0 set H_{0}^{0,0}=1."""
    w.Hwedge[1] = 1
end


@inbounds function _step_2!(w::WignerWorkspace, expiβ::Complex)
    """Compute values H^{0,m}_{n}(β)for m=0,...,n and H^{0,m}_{n+1}(β) for m=0,...,n+1 using Eq. (32):

        H^{0,m}_{n}(β) = (-1)^m √((n-|m|)! / (n+|m|)!) P^{|m|}_{n}(cos β)
                       = (-1)^m (sin β)^m P̂^{|m|}_{n}(cos β) / √(k (2n+1))

    This function computes the associated Legendre functions directly by recursion
    as explained by Holmes and Featherstone (2002), doi:10.1007/s00190-002-0216-2.
    Note that I had to adjust certain steps for consistency with the notation
    assumed by arxiv:1403.7698 -- mostly involving factors of (-1)**m.

    NOTE: Though not specified in arxiv:1403.7698, there is not enough information
    for step 4 unless we also use symmetry to set H^{1,0}_{n} here.  Similarly,
    step 5 needs additional information, which depends on setting H^{0, -1}_{n}
    from its symmetric equivalent H^{0, 1}_{n} in this step.

    """
    gvalues, hvalues, n_max, mp_max = g(w.W), h(w.W), ℓₘₐₓ(w.W), m′ₘₐₓ(w.W)
    Hwedge, Hextra, Hv = w.Hwedge, w.Hextra, w.Hv
    cosβ = expiβ.re
    sinβ = expiβ.im
    TW = T(w.W)
    sqrt3 = √TW(3)
    inverse_sqrt2 = inv(√TW(2))

    if n_max > 0
        # n = 1
        n0n_index = WignerHindex(1, 0, 1, mp_max)
        nn_index = nm_index(1, 1)
        Hwedge[n0n_index] = sqrt3  # Un-normalized
        Hwedge[n0n_index-1] = (gvalues[nn_index-1] * cosβ) * inverse_sqrt2  # Normalized
        # n = 2, ..., n_max+1
        for n in 2:n_max+1
            if n <= n_max
                n0n_index = WignerHindex(n, 0, n, mp_max)
                H = Hwedge
            else
                n0n_index = n + 1
                H = Hextra
            end
            nm10nm1_index = WignerHindex(n-1, 0, n-1, mp_max)
            nn_index = nm_index(n, n)
            constant = √(1 + 1/TW(2n))
            gi = gvalues[nn_index-1]
            # m = n
            H[n0n_index] = constant * Hwedge[nm10nm1_index]
            # m = n-1
            H[n0n_index-1] = gi * cosβ * H[n0n_index]
            # m = n-2, ..., 1
            for i in 2:n-1
                gi = gvalues[nn_index-i]
                hi = hvalues[nn_index-i]
                H[n0n_index-i] = gi * cosβ * H[n0n_index-i+1] - hi * sinβ^2 * H[n0n_index-i+2]
            end
            # m = 0, with normalization
            constant = 1 / √TW(4n+2)
            gi = gvalues[nn_index-n]
            hi = hvalues[nn_index-n]
            H[n0n_index-n] = (gi * cosβ * H[n0n_index-n+1] - hi * sinβ^2 * H[n0n_index-n+2]) * constant
            # Now, loop back through, correcting the normalization for this row, except for n=n element
            prefactor = constant
            for i in 1:n-1
                prefactor *= sinβ
                H[n0n_index-n+i] *= prefactor
            end
            # Supply extra edge cases as noted in docstring
            if n <= n_max
                Hv[nm_index(n, 1)] = Hwedge[WignerHindex(n, 0, 1, mp_max)]
                Hv[nm_index(n, 0)] = Hwedge[WignerHindex(n, 0, 1, mp_max)]
            end
        end
        # Correct normalization of m=n elements
        prefactor = one(TW)
        for n in 1:n_max
            prefactor *= sinβ
            Hwedge[WignerHindex(n, 0, n, mp_max)] *= prefactor / √TW(4n+2)
        end
        for n in [n_max+1]
            prefactor *= sinβ
            Hextra[n+1] *= prefactor / √TW(4n+2)
        end
        # Supply extra edge cases as noted in docstring
        Hv[nm_index(1, 1)] = Hwedge[WignerHindex(1, 0, 1, mp_max)]
        Hv[nm_index(1, 0)] = Hwedge[WignerHindex(1, 0, 1, mp_max)]
    end
end


@inbounds function _step_3!(w::WignerWorkspace, expiβ::Complex)
    """Use relation (41) to compute H^{1,m}_{n}(β) for m=1,...,n.  Using symmetry and shift
    of the indices this relation can be written as

        b^{0}_{n+1} H^{1, m}_{n} =   (b^{−m−1}_{n+1} (1−cosβ))/2 H^{0, m+1}_{n+1}
                                   − (b^{ m−1}_{n+1} (1+cosβ))/2 H^{0, m−1}_{n+1}
                                   − a^{m}_{n} sinβ H^{0, m}_{n+1}

    """
    avalues = a(w.W)
    bvalues = b(w.W)
    n_max = ℓₘₐₓ(w.W)
    mp_max = m′ₘₐₓ(w.W)
    Hwedge = w.Hwedge
    Hextra = w.Hextra
    cosβ = expiβ.re
    sinβ = expiβ.im

    if n_max > 0 && mp_max > 0
        for n in 1:n_max
            # m = 1, ..., n
            i1 = WignerHindex(n, 1, 1, mp_max)
            if n+1 <= n_max
                i2 = WignerHindex(n+1, 0, 0, mp_max)
                H2 = Hwedge
            else
                i2 = 1
                H2 = Hextra
            end
            i3 = nm_index(n+1, 0)
            i4 = nabsm_index(n, 1)
            inverse_b5 = inv(bvalues[i3])
            for i in 0:n-1
                b6 = bvalues[-i+i3-2]
                b7 = bvalues[i+i3]
                a8 = avalues[i+i4]
                Hwedge[i+i1] = inverse_b5 * (
                    (
                          b6 * (1-cosβ) * H2[i+i2+2]
                        - b7 * (1+cosβ) * H2[i+i2]
                    ) / 2
                    - a8 * sinβ * H2[i+i2+1]
                )
            end
        end
    end
end


@inbounds function _step_4!(w::WignerWorkspace)
    """Recursively compute H^{m'+1, m}_{n}(β) for m'=1,...,n−1, m=m',...,n using relation (50) resolved
    with respect to H^{m'+1, m}_{n}:

      d^{m'}_{n} H^{m'+1, m}_{n} =   d^{m'−1}_{n} H^{m'−1, m}_{n}
                                   − d^{m−1}_{n} H^{m', m−1}_{n}
                                   + d^{m}_{n} H^{m', m+1}_{n}

    (where the last term drops out for m=n).

    """
    dvalues, n_max, mp_max = d(w.W), ℓₘₐₓ(w.W), m′ₘₐₓ(w.W)
    Hwedge, Hv = w.Hwedge, w.Hv

    if n_max > 0 && mp_max > 0
        for n in 2:n_max
            for mp in 1:min(n, mp_max)-1
                # m = m', ..., n-1
                # i1 = WignerHindex(n, mp+1, mp, mp_max)
                i1 = WignerHindex(n, mp+1, mp+1, mp_max) - 1
                i2 = WignerHindex(n, mp-1, mp, mp_max)
                # i3 = WignerHindex(n, mp, mp-1, mp_max)
                i3 = WignerHindex(n, mp, mp, mp_max) - 1
                i4 = WignerHindex(n, mp, mp+1, mp_max)
                i5 = nm_index(n, mp)
                i6 = nm_index(n, mp-1)
                inverse_d5 = inv(dvalues[i5])
                d6 = dvalues[i6]
                for i in [0]
                    d7 = dvalues[i+i6]
                    d8 = dvalues[i+i5]
                    Hv[i+nm_index(n, mp+1)] = inverse_d5 * (
                          d6 * Hwedge[i+i2]
                        - d7 * Hv[i+nm_index(n, mp)]
                        + d8 * Hwedge[i+i4]
                    )
                end
                for i in 1:n-mp-1
                    d7 = dvalues[i+i6]
                    d8 = dvalues[i+i5]
                    Hwedge[i+i1] = inverse_d5 * (
                          d6 * Hwedge[i+i2]
                        - d7 * Hwedge[i+i3]
                        + d8 * Hwedge[i+i4]
                    )
                end
                # m = n
                for i in [n-mp]
                    Hwedge[i+i1] = inverse_d5 * (
                          d6 * Hwedge[i+i2]
                        - dvalues[i+i6] * Hwedge[i+i3]
                    )
                end
            end
        end
    end
end


@inbounds function _step_5!(w::WignerWorkspace)
    """Recursively compute H^{m'−1, m}_{n}(β) for m'=−1,...,−n+1, m=−m',...,n using relation (50)
    resolved with respect to H^{m'−1, m}_{n}:

      d^{m'−1}_{n} H^{m'−1, m}_{n} = d^{m'}_{n} H^{m'+1, m}_{n}
                                     + d^{m−1}_{n} H^{m', m−1}_{n}
                                     − d^{m}_{n} H^{m', m+1}_{n}

    (where the last term drops out for m=n).

    NOTE: Although arxiv:1403.7698 specifies the loop over mp to start at -1, I
    find it necessary to start at 0, or there will be missing information.  This
    also requires setting the (m',m)=(0,-1) components before beginning this loop.

    """
    dvalues, n_max, mp_max = d(w.W), ℓₘₐₓ(w.W), m′ₘₐₓ(w.W)
    Hwedge, Hv = w.Hwedge, w.Hv

    if n_max > 0 && mp_max > 0
        for n in 0:n_max
            for mp in 0:-1:1-min(n, mp_max)
                # m = -m', ..., n-1
                # i1 = WignerHindex(n, mp-1, -mp, mp_max)
                i1 = WignerHindex(n, mp-1, -mp+1, mp_max) - 1
                # i2 = WignerHindex(n, mp+1, -mp, mp_max)
                i2 = WignerHindex(n, mp+1, -mp+1, mp_max) - 1
                # i3 = WignerHindex(n, mp, -mp-1, mp_max)
                i3 = WignerHindex(n, mp, -mp, mp_max) - 1
                i4 = WignerHindex(n, mp, -mp+1, mp_max)
                i5 = nm_index(n, mp-1)
                i6 = nm_index(n, mp)
                i7 = nm_index(n, -mp-1)
                i8 = nm_index(n, -mp)
                inverse_d5 = inv(dvalues[i5])
                d6 = dvalues[i6]
                for i in [0]
                    d7 = dvalues[i+i7]
                    d8 = dvalues[i+i8]
                    if mp == 0
                        Hv[i+nm_index(n, mp-1)] = inverse_d5 * (
                              d6 * Hv[i+nm_index(n, mp+1)]
                            + d7 * Hv[i+nm_index(n, mp)]
                            - d8 * Hwedge[i+i4]
                        )
                    else
                        Hv[i+nm_index(n, mp-1)] = inverse_d5 * (
                              d6 * Hwedge[i+i2]
                            + d7 * Hv[i+nm_index(n, mp)]
                            - d8 * Hwedge[i+i4]
                        )
                    end
                end
                for i in 1:n+mp-1
                    d7 = dvalues[i+i7]
                    d8 = dvalues[i+i8]
                    Hwedge[i+i1] = inverse_d5 * (
                          d6 * Hwedge[i+i2]
                        + d7 * Hwedge[i+i3]
                        - d8 * Hwedge[i+i4]
                    )
                end
                # m = n
                for i in [n+mp]
                    Hwedge[i+i1] = inverse_d5 * (
                          d6 * Hwedge[i+i2]
                        + dvalues[i+i7] * Hwedge[i+i3]
                    )
                end
            end
        end
    end
end


function H!(w::WignerWorkspace, expiβ::Complex)
    _step_1!(w)
    _step_2!(w, expiβ)
    _step_3!(w, expiβ)
    _step_4!(w)
    _step_5!(w)
    w.Hwedge, w.Hextra
end
