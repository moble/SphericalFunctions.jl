md"""
# NIST DLMF (2026)

!!! info "Summary"
    The NIST Digital Library of Mathematical Functions defines the spherical harmonics in
    terms of Ferrers functions in a way that agrees with the definition used in the
    `SphericalFunctions` package.

The [NIST Digital Library of Mathematical Functions](@cite NIST_DLMF) is the successor to
Abramowitz and Stegun, and is probably the most careful modern reference for the special
functions.  The statements below are quoted from Release 1.2.7 of 2026-06-15, accessed
2026-09-08.

The DLMF distinguishes between Ferrers' function (of the first kind) ``\mathsf{P}_\nu^\mu``
and the associated Legendre function (of the first kind) ``P_\nu^\mu``.  Here, ``\nu`` is
called the "degree" and ``\mu`` is called the "order".  We can see from their definitions in
Eqs. [14.3.1](https://dlmf.nist.gov/14.3#E1) and [14.3.6](https://dlmf.nist.gov/14.3#E6),
respectively, that they differ only by a phase ``e^{\pm iπ\mu/2}`` (the sign depending on
the side from which the cut along ``(-1, 1)`` is approached; see [Eq.
14.23.1](https://dlmf.nist.gov/14.23#E1)).

For integer degree and order, we have [Eq. 14.7.10](https://dlmf.nist.gov/14.7#E10)
```math
    \mathsf{P}^{m}_{n}\left(x\right)
    =
    (-1)^{m+n}
    \frac{\left(1-x^{2}\right)^{m/2}}{2^{n}n!}
    \frac{{\mathrm{d}}^{m+n}}{{\mathrm{d}x}^{m+n}}
    \left(1-x^{2}\right)^{n},
```
or [Eq. 14.7.14](https://dlmf.nist.gov/14.7#E14)
```math
    P^{m}_{n}\left(x\right)
    =
    \frac{\left(x^{2}-1\right)^{m/2}}{2^{n}n!}
    \frac{{\mathrm{d}}^{m+n}}{{\mathrm{d}x}^{m+n}}
    \left(x^{2}-1\right)^{n}.
```
Note that Eq. 14.7.10 includes the Condon–Shortley phase ``(-1)^m`` in the Ferrers function
itself.  The spherical harmonics are then defined in [Eq.
14.30.1](https://dlmf.nist.gov/14.30#E1), for ``0 \leq θ \leq π`` and ``0 \leq ϕ \leq 2π``,
as
```math
    Y_{l, m}\left(θ,ϕ\right)
    =
    \left(\frac{(l-m)!(2l+1)}{4π(l+m)!}\right)^{1/2}
    e^{imϕ}
    \mathsf{P}_{l}^{m}\left(\cos θ\right),
```
using the Ferrers function.  [Eq. 14.30.6](https://dlmf.nist.gov/14.30#E6) records the
conjugation relation
```math
    Y_{l,-m}(θ, ϕ) = (-1)^m\, \overline{Y_{l,m}(θ, ϕ)}.
```
Since Eq. 14.7.10 is only stated for ``m \geq 0``, we use Eq. 14.30.6 for negative ``m``.
This is the standard Condon–Shortley convention, so we expect agreement with [our spherical
harmonics](@ref summary_spherical_harmonics).  The DLMF does not define Wigner's ``D``
matrices.

## Implementing formulas

We begin by writing code that implements the formulas from the DLMF.  We encapsulate the
formulas in a module so that we can test them against the `SphericalFunctions` package.
"""

using TestItems: @testitem  #hide
@testitem "NIST DLMF conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide

module NIST_DLMF
#+

# We'll use some predefined utilities to make the code look more like the equations,
# including `∂ⁿ`, which computes derivatives symbolically so that we can transcribe the
# formulas literally.
import ..ConventionsUtilities: 𝒾, ❗, ∂ⁿ
#+

# Ferrers' function of the first kind, Eq. 14.7.10:
function 𝖯(n, m, x)
    (-1)^(m+n) * (1 - x^2)^(m/2) / (2^n * (n)❗) * ∂ⁿ(x -> (1 - x^2)^n, m+n)(x)
end
#+

# The associated Legendre function of the first kind, Eq. 14.7.14.  For ``|x| < 1`` the
# factor ``(x^2-1)^{m/2}`` is complex; with the principal branch of the power (which
# corresponds to approaching the cut from above), it differs from the Ferrers function by
# ``e^{-iπm/2}``.  We include it only to check that relation.
function P(n, m, x)
    (complex(x^2 - 1))^(m/2) / (2^n * (n)❗) * ∂ⁿ(x -> (x^2 - 1)^n, m+n)(x)
end
#+

# The spherical harmonics, Eq. 14.30.1 for ``m \geq 0`` and Eq. 14.30.6 for ``m < 0``.  We
# capture the floating-point type `T` to ensure that we don't lose precision when converting
# π and the factorials to floating-point numbers.
function Y(l, m, θ::T, ϕ::T) where {T<:Real}
    if m < 0
        return (-1)^m * conj(Y(l, -m, θ, ϕ))
    end
    √T((l-m)❗ * (2l+1) / (4big(π) * (l+m)❗)) * exp(𝒾 * m * ϕ) * T(𝖯(l, m, cos(θ)))
end
#+

end  # module NIST_DLMF
#+

# ## Tests
#
# We can now test the functions against the equivalent functions from the
# `SphericalFunctions` package.  We will need to test approximate floating-point equality,
# so we set absolute and relative tolerances (respectively) in terms of the machine epsilon:
ϵₐ = 100eps()
ϵᵣ = 1000eps()
#+

# We only test up to
ℓₘₐₓ = 4
#+
# because the formulas are slow, and this will be sufficient to sort out any sign or
# normalization differences, which are the most likely source of error.

# We will also need the imaginary unit from the utilities module.
import .ConventionsUtilities: 𝒾
#+

# First, the relation between the Ferrers function and the associated Legendre function,
# ``P_n^m(x) = e^{-iπm/2}\, \mathsf{P}_n^m(x)`` for ``-1 < x < 1`` with the principal
# branch of the power:
for x ∈ range(-0.9, 0.9, length=7)
    for n ∈ 0:ℓₘₐₓ, m ∈ 0:n
        @test NIST_DLMF.P(n, m, x) ≈ exp(-𝒾 * π * m / 2) * NIST_DLMF.𝖯(n, m, x) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Now, the spherical harmonics agree with ours:
for (θ, ϕ) ∈ θϕrange(Float64, 7)
    for (ℓ, m) ∈ ℓmrange(ℓₘₐₓ)
        @test NIST_DLMF.Y(ℓ, m, θ, ϕ) ≈ ConventionsUtilities.Y(ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# This successful test shows that the spherical harmonics defined by the DLMF agree with the
# spherical harmonics defined by the `SphericalFunctions` package.

end  #hide
