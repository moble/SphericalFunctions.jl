raw"""
Formulas and conventions from [NIST's Digital Library of Mathematical Functions](@cite
NIST_DLMF).

The DLMF distinguishes between Ferrer's function (of the first kind) ``\mathup{P}_\nu^\mu``
and the associated Legendre function (of the first kind) ``P_\nu^\mu``.  Here, ``\nu`` is
called the "degree" and ``\mu`` is called the "order".  We can see from their definitions in
Eqs. [14.3.1](http://dlmf.nist.gov/14.3#E1) and [14.3.6](http://dlmf.nist.gov/14.3#E6),
respectively, that they differ only by a factor of ``(-1)^{\mu/2}``.

For integer degree and order, we have [Eq. 14.7.10](http://dlmf.nist.gov/14.7#E10)
```math
    \mathsf{P}^{m}_{n}\left(x\right)
    =
    (-1)^{m+n}
    \frac{\left(1-x^{2}\right)^{m/2}}{2^{n}n!}
    \frac{{\mathrm{d}}^{m+n}}{{\mathrm{d}x}^{m+n}}
    \left(1-x^{2}\right)^{n}
```
or [Eq. 14.7.14](http://dlmf.nist.gov/14.7#E14)
```math
    P^{m}_{n}\left(x\right)
    =
    \frac{\left(x^{2}-1\right)^{m/2}}{2^{n}n!}
    \frac{{\mathrm{d}}^{m+n}}{{\mathrm{d}x}^{m+n}}
    \left(x^{2}-1\right)^{n}.
```
And for the spherical harmonics, [Eq. 14.30.1](http://dlmf.nist.gov/14.30#E1) gives
```math
    Y_{\ell, m}\left(\theta,\phi\right)
    =
    \left(\frac{(\ell-m)!(2\ell+1)}{4\pi(\ell+m)!}\right)^{1/2}
    \mathsf{e}^{im\phi}
    \mathsf{P}_{\ell}^{m}\left(\cos\theta\right).
```

"""

@testmodule NIST_DLMF begin

import FastDifferentiation

const 𝒾 = im

include("../utilities/naive_factorial.jl")
import .NaiveFactorials: ❗

function P(x::T, n, m) where {T}
    if m > n
        zero(x)
    else
        (-1)^(m+n) * (1-x^2)^(m/2) / T(2^n * (n)❗) *
        FastDifferentiation.derivative(x -> (1-x^2)^n, m+n)
    end
end

function 𝘗(x::T, n, m) where {T}
    if m > n
        zero(x)
    else
        (x^2-1)^(m/2) / T(2^n * (n)❗) *
        FastDifferentiation.derivative(x -> (x^2-1)^n, m+n)
    end
end

function Y(ℓ, m, θ::T, φ::T) where {T}
    if m > ℓ
        zero(Complex{T})
    else
        let π = T(π)
            √T((ℓ-m)❗ * (2ℓ+1) / (4π * (ℓ+m)❗)) *
            exp(𝒾*m*φ) * P(cos(θ), ℓ, m)
        end
    end
end

end  # @testmodule NIST_DLMF


@testitem "NIST_DLMF conventions" setup=[Utilities, NIST_DLMF] begin
    using Random
    using Quaternionic: from_spherical_coordinates

    Random.seed!(1234)
    const T = Float64
    const ℓₘₐₓ = 3
    ϵₐ = 8eps(T)
    ϵᵣ = 20eps(T)

    # TODO: Add tests

end  # @testitem NIST_DLMF
