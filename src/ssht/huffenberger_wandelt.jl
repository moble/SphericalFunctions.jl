module HuffenbergerWandelt

raw"""
The Wigner ``d`` matrix corresponds to a rotation about the ``y`` axis by an angle ``β``:
``\exp\left[ \beta 𝐣 / 2 \right]``.  But computing the ``d`` matrix for a general angle is
a little awkward.  However, computing the ``d`` matrix for a rotation by ``π/2`` is somewhat
simpler, and computing the full ``D`` matrix for a rotation about the ``z`` axis is trivial
— it's just a phase factor.  Thus, we can re-express ``d`` for a general angle ``β`` in
terms of ``d`` for the angle ``π/2`` and a phase factor.

We start with the fact that the ``y`` axis equals the ``z`` axis rotated by ``π/2`` about
the ``x`` axis:
```math
𝐣 = e^{\frac{\pi}{2} 𝐢/ 2}\, 𝐤\, e^{-\frac{\pi}{2} 𝐢/ 2}.
```
So we have
```math
e^{\beta 𝐣 / 2}
=
e^{\frac{\pi}{2} 𝐢/ 2}\, e^{\beta 𝐤 / 2}\, e^{-\frac{\pi}{2} 𝐢/ 2}.
```
Unfortunately, this expression involves the ``x`` axis, which we don't want.  But we can
similarly express the ``x`` axis in terms of the ``y`` axis rotated about the ``z`` axis:
```math
𝐢 = e^{-\frac{\pi}{2} 𝐤/ 2}\, 𝐣\, e^{\frac{\pi}{2} 𝐤/ 2}.
```
And now we can use this in our first expression to find
```math
e^{\beta 𝐣 / 2}
=
e^{-\frac{\pi}{2} 𝐤/ 2}\, e^{\frac{\pi}{2} 𝐣/ 2}\, e^{\frac{\pi}{2} 𝐤/ 2}\,
e^{\beta 𝐤 / 2}\,
e^{\frac{\pi}{2} 𝐤/ 2}\, e^{-\frac{\pi}{2} 𝐣/ 2}\, e^{-\frac{\pi}{2} 𝐤/ 2}.
```
Now, we can use this expansion to find an expression for the ``d`` matrix value:
```math
\begin{align}
d^{\ell}_{m', m}(\beta)
&=
\mathfrak{D}^{\ell}_{m', m}\left(e^{\beta 𝐣 / 2}\right) \\
&=
\mathfrak{D}^{\ell}_{m', m_1}\left(e^{-\frac{\pi}{2} 𝐤/ 2}\right)\,
\mathfrak{D}^{\ell}_{m_1, m_2}\left(e^{\frac{\pi}{2} 𝐣/ 2}\right)\,
\mathfrak{D}^{\ell}_{m_2, m_3}\left(e^{\frac{\pi}{2} 𝐤/ 2}\right)\,
\mathfrak{D}^{\ell}_{m_3, m_4}\left(e^{\beta 𝐤 / 2}\right)\, \\
&\quad \times
\mathfrak{D}^{\ell}_{m_4, m_5}\left(e^{\frac{\pi}{2} 𝐤/ 2}\right)\,
\mathfrak{D}^{\ell}_{m_5, m_6}\left(e^{-\frac{\pi}{2} 𝐣/ 2}\right)\,
\mathfrak{D}^{\ell}_{m_6, m}\left(e^{-\frac{\pi}{2} 𝐤/ 2}\right) \\
&=
\delta_{m', m_1} e^{im'\frac{\pi}{2}}\,
\mathfrak{D}^{\ell}_{m_1, m_2}\left(e^{\frac{\pi}{2} 𝐣/ 2}\right)\,
\delta_{m_2, m_3} e^{-im_2\frac{\pi}{2}}\,
\mathfrak{D}^{\ell}_{m_3, m_4}\left(e^{\beta 𝐤 / 2}\right)\, \\
&\quad \times
\delta_{m_4, m_5} e^{-im_4\frac{\pi}{2}}\,
\mathfrak{D}^{\ell}_{m_5, m_6}\left(e^{-\frac{\pi}{2} 𝐣/ 2}\right)\,
\delta_{m_6, m} e^{im\frac{\pi}{2}} \\
&=
e^{im'\frac{\pi}{2}}\, e^{-im''\frac{\pi}{2}}\,
e^{-im'''\frac{\pi}{2}}\, e^{im\frac{\pi}{2}}\,
&\quad \times
\mathfrak{D}^{\ell}_{m', m''}\left(e^{\frac{\pi}{2} 𝐣/ 2}\right)\,
\mathfrak{D}^{\ell}_{m'', m'''}\left(e^{\beta 𝐤 / 2}\right)\, \\
\mathfrak{D}^{\ell}_{m''', m}\left(e^{-\frac{\pi}{2} 𝐣/ 2}\right) \\
&=
e^{im'\frac{\pi}{2}}\, e^{-im''\frac{\pi}{2}}\,
e^{-im'''\frac{\pi}{2}}\, e^{im\frac{\pi}{2}}\,
&\quad \times
\mathfrak{D}^{\ell}_{m', m''}\left(e^{\frac{\pi}{2} 𝐣/ 2}\right)\,
e^{-im''\beta}\,
\mathfrak{D}^{\ell}_{m'', m}\left(e^{-\frac{\pi}{2} 𝐣/ 2}\right) \\
&=
i^{m'+m-2m''}\,
d^{\ell}_{m', m''}\left(e^{\frac{\pi}{2} 𝐣/ 2}\right)\,
e^{-im''\beta}\,
d^{\ell}_{m'', m}\left(e^{-\frac{\pi}{2} 𝐣/ 2}\right) \\
&=
i^{m'+m}(-1)^{m''}\,
d^{\ell}_{m', m''}\left(\frac{\pi}{2}\right)\,
e^{-im''\beta}\,
d^{\ell}_{m'', m}\left(-\frac{\pi}{2}\right) \\
&=
i^{m'+m}(-1)^{m}\,
d^{\ell}_{m', m''}\left(\frac{\pi}{2}\right)\,
e^{-im''\beta}\,
d^{\ell}_{m'', m}\left(\frac{\pi}{2}\right) \\
&=
i^{m'-m}\,
d^{\ell}_{m', m''}\left(\frac{\pi}{2}\right)\,
e^{-im''\beta}\,
d^{\ell}_{m'', m}\left(\frac{\pi}{2}\right)
\end{align}
```

"""



function deducelmax(salm)
    N = length(salm)
    lmax = floor(Int, sqrt(N) - 1)
    if (lmax + 1)^2 != N
        error("Cannot deduce lmax from salm length $N")
    end
    return lmax
end


"""
    map2salm(𝒯, f)

Map function values `f` to spin-weighted spherical-harmonic mode weights `f̃` using the
`SSHT` object `𝒯`.
"""
function map2salm end
function map2salm! end

"""
    salm2map(salm, s, [lmax=deducelmax(salm), Ntheta=2lmax+1, Nphi=2lmax+1])

Compatibility function for `spinsfast` package.  Converts mode weights of spin-weighted
function, stored in the array `salm`, to values on an equiangular grid with `Ntheta` points
in θ and `Nphi` points in ϕ.  The spin weight is `s` and the maximum ℓ value in the input
data is `lmax`.

The optional arguments `lmax`, `Ntheta`, and `Nphi` are available for backwards
compatibility, but if they are given values that are inconsistent with the defaults, an
error will be thrown.

See also [`map2salm`](@ref) for the (rough) inverse operation.


# Notes

The input `salm` data should be given in increasing order of ℓ value, always starting with
`(ℓ, m) = (0, 0)` even if `s` is nonzero, proceeding to `(1, -1)`, `(1, 0)`, `(1, 1)`, etc.
Explicitly, the ordering should match this:

    [f̃[ℓ, m] for ℓ ∈ 0:lmax for m ∈ -ℓ:ℓ]

The output data are presented on this grid of spherical coordinates:

    [f(θ, ϕ) for θ ∈ LinRange(0, π, Ntheta), ϕ ∈ LinRange(0, 2π, Nphi+1)[begin:end-1]]

Note that `map2salm` and `salm2map` are not true inverses of each other for several reasons.
First, modes with `ell < |s|` should always be zero; they are simply assumed to be zero on
input to `salm2map`.  It is also possible to define a `map` function that violates this
assumption -- for example, having a nonzero average value over the sphere; if the function
has nonzero spin `s`, this is impossible.  Also, it is possible to define a map of a
function with so much angular dependence that it cannot be captured with the given `lmax`
value.  For example, a discontinuous function will never be perfectly resolved.

"""
function salm2map(salm, s, lmax=deducelmax(salm), Ntheta=2lmax+1, Nphi=2lmax+1)
    # Not implemented
    error("salm2map is not yet implemented")
end

using FFTW
using LinearAlgebra

"""
    lm_index(ℓ, m, lmax) -> Int

Map (ℓ, m) to the flattened index matching spinsfast ordering.
Inspired by Python wrapper ordering in `python/__init__.py` lines 31–49 and C helper `ind_lm` in `cextension.c` lines ~20–53.
"""
lm_index(ℓ, m, lmax) = ℓ^2 + (m + ℓ) + 1  # 1-based for Julia

"""
    ind_lm(idx, lmax) -> (ℓ, m)

Inverse of `lm_index`.
"""
function ind_lm(idx, lmax)
    i0 = idx - 1
    ℓ = floor(Int, sqrt(i0))
    while (ℓ + 1)^2 <= i0
        ℓ += 1
    end
    m = i0 - ℓ^2 - ℓ
    return ℓ, m
end

"""
    salm2map(salm, s, lmax; Nθ = 2lmax+1, Nφ = 2lmax+1)

Spin-weighted (s) spherical-harmonic synthesis to grid samples.
Follows structure of `spinsfast_salm2map` in `code/spinsfast_backward_transform.c` lines 131–200 and wrapper in `python/cextension.c` lines 56–90.
"""
function salm2map(salm::AbstractVector{<:Complex}, s::Integer, lmax::Integer;
                  Nθ::Integer = 2lmax + 1, Nφ::Integer = 2lmax + 1)

    Nm = 2lmax + 1

    # Placeholder Wigner d(π/2); replace with stable helper like `wdhp_TN_helper` (see `spinsfast_backward_transform.c` lines 151–160).
    wig_d_halfpi(ℓ, m, mp) = wigner_d_halfpi(ℓ, m, mp)

    # G matrix (mode-coupled intermediate), cf. `Gmm` in `spinsfast_backward_transform.c` lines 151–186.
    G = zeros(ComplexF64, Nm, Nm)

    # Build G via Δ products, cf. `spinsfast_backward_Gmm` call in `spinsfast_salm2map` lines 151–160.
    for ℓ in 0:lmax
        for m in -ℓ:ℓ
            alm = salm[lm_index(ℓ, m, lmax)]
            for mp in -ℓ:ℓ
                Δ1 = wig_d_halfpi(ℓ, m, mp)
                Δ2 = wig_d_halfpi(ℓ, mp, s)
                G[m + lmax + 1, mp + lmax + 1] += alm * Δ1 * Δ2
            end
        end
    end

    # FFT work array F, matching the quadrant placement in `spinsfast_backward_transform.c` lines 151–186.
    F = zeros(ComplexF64, Nθ, Nφ)
    limit = lmax
    for mp in 0:limit
        for m in 0:limit
            F[mp + 1, m + 1] = G[mp + 1, m + 1]                  # ++
            if m > 0
                F[mp + 1, Nφ - m + 1] = G[mp + 1, Nm - m + 1]    # +-
            end
            if mp > 0
                F[Nθ - mp + 1, m + 1] = G[Nm - mp + 1, m + 1]    # -+
            end
            if mp > 0 && m > 0
                F[Nθ - mp + 1, Nφ - m + 1] = G[Nm - mp + 1, Nm - m + 1]  # --
            end
        end
    end

    # Inverse FFT (2D), mirroring `fftw_execute` usage in `spinsfast_backward_transform.c` lines 188–198.
    f = ifft(ifft(F, 1), 2) .* Nm  # scale to align with C conventions
    return f
end

"""
    map2salm(f, s, lmax)

Spin-weighted analysis: grid samples -> spherical-harmonic coefficients.
Inspired by forward path (`spinsfast_map2salm`) inverse of the above; quadrant extraction mirrors `spinsfast_backward_transform.c` lines 151–186.
"""
function map2salm(f::AbstractMatrix{<:Complex}, s::Integer, lmax::Integer)
    Nθ, Nφ = size(f)
    Nm = 2lmax + 1

    # Forward FFT (2D), inverse of salm2map’s iFFT; normalization mirrors scaling above.
    F = fft(fft(f, 1), 2) ./ Nm

    # Reconstruct G from quadrants, inverse of the packing in `salm2map`.
    G = zeros(ComplexF64, Nm, Nm)
    for mp in 0:lmax
        for m in 0:lmax
            G[mp + 1, m + 1] = F[mp + 1, m + 1]                     # ++
            if m > 0
                G[mp + 1, Nm - m + 1] = F[mp + 1, Nφ - m + 1]       # +-
            end
            if mp > 0
                G[Nm - mp + 1, m + 1] = F[Nθ - mp + 1, m + 1]       # -+
            end
            if mp > 0 && m > 0
                G[Nm - mp + 1, Nm - m + 1] = F[Nθ - mp + 1, Nφ - m + 1]  # --
            end
        end
    end

    # Placeholder Wigner d(π/2); replace with stable helper (see `spinsfast_backward_transform.c` lines 151–160).
    wig_d_halfpi(ℓ, m, mp) = wigner_d_halfpi(ℓ, m, mp)

    # Project G back to alm, paralleling the inverse of `spinsfast_backward_Gmm`.
    salm = zeros(ComplexF64, (lmax + 1)^2)
    for ℓ in 0:lmax
        for m in -ℓ:ℓ
            acc = 0.0 + 0im
            for mp in -ℓ:ℓ
                Δ1 = wig_d_halfpi(ℓ, m, mp)
                Δ2 = wig_d_halfpi(ℓ, mp, s)
                acc += conj(Δ1) * conj(Δ2) * G[m + lmax + 1, mp + lmax + 1]
            end
            salm[lm_index(ℓ, m, lmax)] = acc
        end
    end
    return salm
end


end # module HuffenbergerWandelt
