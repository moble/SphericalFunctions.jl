
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


@inbounds function _step_1!(w::WignerMatrixCalculator)
    """If n=0 set H_{0}^{0,0}=1."""
    w.Hwedge[1] = 1
end


@inbounds function _step_2!(w::WignerMatrixCalculator, expiβ::Complex)
    """Compute values H^{0,m}_{n}(β)for m=0,...,n and H^{0,m}_{n+1}(β) for m=0,...,n+1

     Uses Eq. (32) of Gumerov-Duraiswami (2014) [arxiv:1403.7698]:

        H^{0,m}_{n}(β) = (-1)^m √((n-|m|)! / (n+|m|)!) P^{|m|}_{n}(cos β)
                       = (-1)^m P̄^{|m|}_{n}(cos β) / √(k (2n+1))

    Here, k=1 for m=0, and k=2 for m>0, and

        P̄ = √{k(2n+1)(n-m)!/(n+m)!} P

    We use the "fully normalized" associated Legendre functions (fnALF) P̄ because,
    as explained by Xing et al. (2020) [https://doi.org/10.1007/s00190-019-01331-0],
    it is possible to compute these values very efficiently and accurately, while
    also delaying the onset of overflow and underflow.

    NOTE: Though not specified in arxiv:1403.7698, there is not enough information
    for step 4 unless we also use symmetry to set H^{1,0}_{n} here.  Similarly,
    step 5 needs additional information, which depends on setting H^{0, -1}_{n}
    from its symmetric equivalent H^{0, 1}_{n} in this step.

    """
    n_max, mp_max, TW = ℓₘₐₓ(w), m′ₘₐₓ(w), T(w)
    Hwedge, Hextra, Hv = w.Hwedge, w.Hextra, w.Hv
    cosβ = expiβ.re
    sinβ = expiβ.im
    sqrt3 = √TW(3)

    # The general expressions for the constants are Eq. (13) of Xing et al.:
    #
    #   aₙ = √((2n+1)/TW(2n-1))
    #   bₙ = √(2*(n-1)*(2n+1)/TW(n*(2n-1)))
    #   cₙₘ = √(((n+m)*(n-m)*(2n+1)) / TW(2n-1)) / n
    #   dₙₘ = √(((n-m)*(n-m-1)*(2n+1)) / TW(2n-1)) / (2n)
    #   eₙₘ = √(((n+m)*(n+m-1)*(2n+1)) / TW(2n-1)) / (2n)
    #
    # Below, I factor aₙ out of each of these expressions, along with 1/2n where
    # relevant, to avoid divisions.
    #
    # We initialize with Eq. (14), then step through with Eq. (12), to compute all
    # values of P̄.  Finally, we normalize everything again to compute the H
    # values.


    if n_max > 0
        # n = 1
        n0n_index = WignerHindex(1, 0, 1, mp_max)
        Hwedge[n0n_index] = sqrt3 * sinβ
        Hwedge[n0n_index-1] = sqrt3 * cosβ
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
            inv2n = inv(TW(2n))
            aₙ = √((2n+1)/TW(2n-1))
            bₙ = aₙ * √((2*(n-1))/TW(n))

            # m = n
            eₙₘ = inv2n * √TW((2n)*(2n+1))
            H[n0n_index] = sinβ * eₙₘ * Hwedge[nm10nm1_index]

            # m = n-1
            eₙₘ = inv2n * √TW((2n-2)*(2n+1))
            cₙₘ = 2inv2n * √TW(2n+1)
            H[n0n_index-1] = cosβ * cₙₘ * Hwedge[nm10nm1_index] + sinβ * eₙₘ * Hwedge[nm10nm1_index-1]

            # m = n-2, ..., 2
            for i in 2:n-2
                # m = n-i
                cₙₘ = 2inv2n * aₙ * √TW((2n-i)*i)
                dₙₘ = inv2n * aₙ * √TW(i*(i-1))
                eₙₘ = inv2n * aₙ * √TW((2n-i)*(2n-i-1))
                H[n0n_index-i] = (
                    cosβ * cₙₘ * Hwedge[nm10nm1_index-i+1]
                    - sinβ * (
                        dₙₘ * Hwedge[nm10nm1_index-i+2]
                        - eₙₘ * Hwedge[nm10nm1_index-i]
                    )
                )
            end

            # m = 1
            cₙₘ = 2inv2n * aₙ * √TW((n+1)*(n-1))
            dₙₘ = inv2n * aₙ * √TW((n-1)*(n-2))
            eₙₘ = inv2n * aₙ * √TW(2n*(n+1))
            H[n0n_index-n+1] = (
                cosβ * cₙₘ * Hwedge[nm10nm1_index-n+2]
                - sinβ * (
                    dₙₘ * Hwedge[nm10nm1_index-n+3]
                    - eₙₘ * Hwedge[nm10nm1_index-n+1]
                )
            )

            # m = 0
            cₙₘ = aₙ
            dₙₘ = inv2n * aₙ * √TW(n*(n-1))
            eₙₘ = dₙₘ
            H[n0n_index-n] = (
                aₙ * cosβ * Hwedge[nm10nm1_index-n+1]
                - bₙ * sinβ * Hwedge[nm10nm1_index-n+2] / 2
            )

            # Supply extra edge cases as noted in docstring
            if n <= n_max
                Hv[nm_index(n, 1)] = Hwedge[WignerHindex(n, 0, 1, mp_max)]
                Hv[nm_index(n, 0)] = Hwedge[WignerHindex(n, 0, 1, mp_max)]
            end
        end

        # Supply extra edge cases as noted in docstring
        Hv[nm_index(1, 1)] = Hwedge[WignerHindex(1, 0, 1, mp_max)]
        Hv[nm_index(1, 0)] = Hwedge[WignerHindex(1, 0, 1, mp_max)]

        # Normalize, changing P̄ to H values
        for n in 1:n_max+1
            if n <= n_max
                n00_index = WignerHindex(n, 0, 0, mp_max)
                H = Hwedge
            else
                n00_index = 1
                H = Hextra
            end
            const0 = inv(√TW(2n+1))
            const1 = inv(√TW(4n+2))
            H[n00_index] *= const0
            for m in 1:n
                H[n00_index+m] *= const1
            end
            if n <= n_max
                Hv[nm_index(n, 1)] *= -const1
                Hv[nm_index(n, 0)] *= -const1
            end
        end
    end
end


@inbounds function _step_3!(w::WignerMatrixCalculator, expiβ::Complex)
    """Use relation (41) to compute H^{1,m}_{n}(β) for m=1,...,n.  Using symmetry and shift
    of the indices this relation can be written as

        b^{0}_{n+1} H^{1, m}_{n} =   (b^{−m−1}_{n+1} (1−cosβ))/2 H^{0, m+1}_{n+1}
                                   − (b^{ m−1}_{n+1} (1+cosβ))/2 H^{0, m−1}_{n+1}
                                   − a^{m}_{n} sinβ H^{0, m}_{n+1}

    """
    n_max = ℓₘₐₓ(w)
    mp_max = m′ₘₐₓ(w)
    Hwedge = w.Hwedge
    Hextra = w.Hextra
    cosβ = expiβ.re
    sinβ = expiβ.im
    cosβ₊ = (1+cosβ)/2
    cosβ₋ = (1-cosβ)/2
    Tw = T(w)

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
            inverse_b5 = inv(b(n+1, 0, Tw))
            for i in 0:n-1
                # b6 = bvalues[-i+i3-2]
                # b7 = bvalues[i+i3]
                # a8 = avalues[i+i4]
                b6 = b(n+1, -i-2, Tw)
                b7 = b(n+1, i, Tw)
                a8 = a(n, 1 + i, Tw)
                Hwedge[i+i1] = inverse_b5 * (
                    (
                          b6 * cosβ₋ * H2[i+i2+2]
                        - b7 * cosβ₊ * H2[i+i2]
                    )
                    - a8 * sinβ * H2[i+i2+1]
                )
            end
        end
    end
end


@inbounds function _step_4!(w::WignerMatrixCalculator)
    """Recursively compute H^{m'+1, m}_{n}(β) for m'=1,...,n−1, m=m',...,n using relation (50) resolved
    with respect to H^{m'+1, m}_{n}:

      d^{m'}_{n} H^{m'+1, m}_{n} =   d^{m'−1}_{n} H^{m'−1, m}_{n}
                                   − d^{m−1}_{n} H^{m', m−1}_{n}
                                   + d^{m}_{n} H^{m', m+1}_{n}

    (where the last term drops out for m=n).

    """
    n_max, mp_max = ℓₘₐₓ(w), m′ₘₐₓ(w)
    Hwedge, Hv = w.Hwedge, w.Hv
    Tw = T(w)

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
                # inverse_d5 = inv(dvalues[i5])
                # d6 = dvalues[i6]
                inverse_d5 = inv(d(n, mp, Tw))
                d6 = d(n, mp-1, Tw)
                let i=0
                    # d7 = dvalues[i+i6]
                    # d8 = dvalues[i+i5]
                    d7 = d(n, mp-1+i, Tw)
                    d8 = d(n, mp+i, Tw)
                    Hv[i+nm_index(n, mp+1)] = inverse_d5 * (
                          d6 * Hwedge[i+i2]
                        - d7 * Hv[i+nm_index(n, mp)]
                        + d8 * Hwedge[i+i4]
                    )
                end
                for i in 1:n-mp-1
                    # d7 = dvalues[i+i6]
                    # d8 = dvalues[i+i5]
                    d7 = d(n, mp-1+i, Tw)
                    d8 = d(n, mp+i, Tw)
                    Hwedge[i+i1] = inverse_d5 * (
                          d6 * Hwedge[i+i2]
                        - d7 * Hwedge[i+i3]
                        + d8 * Hwedge[i+i4]
                    )
                end
                # m = n
                let i=n-mp
                    Hwedge[i+i1] = inverse_d5 * (
                          d6 * Hwedge[i+i2]
                        - d(n, mp-1+i, Tw) * Hwedge[i+i3]
                    )
                end
            end
        end
    end
end


@inbounds function _step_5!(w::WignerMatrixCalculator)
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
    n_max, mp_max = ℓₘₐₓ(w), m′ₘₐₓ(w)
    Hwedge, Hv = w.Hwedge, w.Hv
    Tw = T(w)

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
                # inverse_d5 = inv(dvalues[i5])
                # d6 = dvalues[i6]
                inverse_d5 = inv(d(n, mp-1, Tw))
                d6 = d(n, mp, Tw)
                let i=0
                    # d7 = dvalues[i+i7]
                    # d8 = dvalues[i+i8]
                    d7 = d(n, -mp-1+i, Tw)
                    d8 = d(n, i-mp, Tw)
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
                    # d7 = dvalues[i+i7]
                    # d8 = dvalues[i+i8]
                    d7 = d(n, -mp-1+i, Tw)
                    d8 = d(n, i-mp, Tw)
                    Hwedge[i+i1] = inverse_d5 * (
                          d6 * Hwedge[i+i2]
                        + d7 * Hwedge[i+i3]
                        - d8 * Hwedge[i+i4]
                    )
                end
                # m = n
                let i=n+mp
                    Hwedge[i+i1] = inverse_d5 * (
                          d6 * Hwedge[i+i2]
                        + d(n, -mp-1+i, Tw) * Hwedge[i+i3]
                    )
                end
            end
        end
    end
end


"""
    H!(w, expiβ)

Compute (a quarter of) the H matrix

WARNING: The returned array will be a view into the `workspace` variable (see
below for an explanation of that).  If you need to call this function again
using the same workspace before extracting all information from the first call,
you should use `numpy.copy` to make a separate copy of the result.

Parameters
----------
expiβ : array_like
    Value of exp(i*β) on which to evaluate the H matrix.

Returns
-------
Hwedge : array
    This is a 1-dimensional array of floats; see below.
workspace : array_like, optional
    A working array like the one returned by Wigner.new_workspace().  If not
    present, this object's default workspace will be used.  Note that it is not
    safe to use the same workspace on multiple threads.  Also see the WARNING
    above.

See Also
--------
d : Compute the full Wigner d matrix
D : Compute the full Wigner 𝔇 matrix
rotate : Avoid computing the full 𝔇 matrix and rotate modes directly
evaluate : Avoid computing the full 𝔇 matrix and evaluate modes directly

Notes
-----
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
function H!(w::WignerMatrixCalculator, expiβ::Complex)
    _step_1!(w)
    _step_2!(w, expiβ)
    _step_3!(w, expiβ)
    _step_4!(w)
    _step_5!(w)
    w.Hwedge
end
