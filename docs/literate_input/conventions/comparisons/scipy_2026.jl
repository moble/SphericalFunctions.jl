md"""
# SciPy (2026)

!!! info "Summary"
    SciPy's `sph_harm_y` agrees with the spherical harmonics used in the
    `SphericalFunctions` package.  Note that its arguments are ordered `(n, m, theta, phi)`
    with `theta` the *polar* angle, unlike the legacy `sph_harm(m, n, theta, phi)`, in which
    `theta` was the azimuthal angle.

[SciPy](@cite SciPy_2026) is the standard scientific library for Python.  Its spherical
harmonics are provided by
[`scipy.special.sph_harm_y`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_harm_y.html),
whose documentation (SciPy 1.18.1, accessed 2026-09-08) gives the signature
```python
scipy.special.sph_harm_y(n, m, theta, phi, *, diff_n=0)
```
and the definition
```math
Y_n^m(θ, φ)
=
\sqrt{\frac{2n+1}{4π} \frac{(n-m)!}{(n+m)!}}\; P_n^m(\cos θ)\, e^{imφ},
```
"where ``P_n^m`` denotes the associated Legendre polynomials (unnormalized)".  The
documentation notes that

> In SciPy `theta` is the polar angle and `phi` is the azimuthal angle.  It is common to see
> the opposite convention.

and that

> SciPy's spherical harmonics include the Condon–Shortley phase because it is part of
> `sph_legendre_p`.

That is, the associated Legendre polynomial carries the Condon–Shortley phase, ``P_n^m(x) =
(-1)^m (1-x^2)^{m/2} \frac{d^m}{dx^m} P_n(x)`` with ``P_n`` the Legendre polynomial, so we
expect agreement with [our spherical harmonics](@ref summary_spherical_harmonics).  For
negative ``m`` we use the standard relation [DLMF 14.9.3](https://dlmf.nist.gov/14.9#E3),
``P_n^{-m}(x) = (-1)^m \frac{(n-m)!}{(n+m)!} P_n^m(x)``, which is what `sph_legendre_p`
evaluates to.

The legacy function
[`scipy.special.sph_harm`](https://docs.scipy.org/doc/scipy-1.15.0/reference/generated/scipy.special.sph_harm.html)
had the signature `sph_harm(m, n, theta, phi)`, with `theta` the *azimuthal* angle and `phi`
the *polar* angle — the definition there reads ``Y_n^m(θ, φ) = \sqrt{\ldots}\, e^{imθ}
P_n^m(\cos φ)`` — and was deprecated in SciPy 1.15.0 ("This function is deprecated and will
be removed in SciPy 1.17.0.  Please use `scipy.special.sph_harm_y` instead.").  Apart from
the argument order, the two functions agree.  SciPy does not provide Wigner's ``D``
matrices.

## Implementing formulas

We begin by writing code that implements the formulas from the SciPy documentation.  We
encapsulate the formulas in a module so that we can test them against the
`SphericalFunctions` package.  A second test — which is skipped by default because it needs
a Python installation, and is run by the scheduled CI workflow — calls the actual SciPy
function through `PythonCall`.
"""

using TestItems: @testitem  #hide
@testitem "SciPy conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide

module SciPy
#+

# We'll use some predefined utilities to make the code look more like the equations,
# including `∂ⁿ`, which computes derivatives symbolically so that we can transcribe the
# Legendre formulas literally.
import ..ConventionsUtilities: 𝒾, ❗, ∂ⁿ
#+

# The Legendre polynomial (Rodrigues' formula) and the associated Legendre polynomial with
# the Condon–Shortley phase, with DLMF 14.9.3 for negative ``m``:
function P(n, x)
    1 / (2^n * (n)❗) * ∂ⁿ(x -> (x^2 - 1)^n, n)(x)
end
function P(n, m, x)
    if m < 0
        return (-1)^m * (n+m)❗ / (n-m)❗ * P(n, -m, x)
    end
    (-1)^m * (1 - x^2)^(m/2) * ∂ⁿ(x -> P(n, x), m)(x)
end
#+

# `sph_harm_y(n, m, theta, phi)`.  We capture the floating-point type `T` to ensure that we
# don't lose precision when converting π and the factorials to floating-point numbers.
function sph_harm_y(n, m, θ::T, φ::T) where {T<:Real}
    √T((2n+1) * (n-m)❗ / (4big(π) * (n+m)❗)) * T(P(n, m, cos(θ))) * exp(𝒾 * m * φ)
end
#+

# The legacy `sph_harm(m, n, theta, phi)`, with the azimuthal angle first:
function sph_harm(m, n, θ::T, φ::T) where {T<:Real}
    √T((2n+1) * (n-m)❗ / (4big(π) * (n+m)❗)) * exp(𝒾 * m * θ) * T(P(n, m, cos(φ)))
end
#+

end  # module SciPy
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

# `sph_harm_y` agrees with ours, and the legacy `sph_harm` agrees with `sph_harm_y` once the
# arguments are reordered:
for (θ, ϕ) ∈ θϕrange(Float64, 7)
    for (ℓ, m) ∈ ℓmrange(ℓₘₐₓ)
        @test SciPy.sph_harm_y(ℓ, m, θ, ϕ) ≈ ConventionsUtilities.Y(ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
        @test SciPy.sph_harm(m, ℓ, ϕ, θ) ≈ SciPy.sph_harm_y(ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# This successful test shows that the spherical harmonics as documented by SciPy agree with
# the spherical harmonics defined by the `SphericalFunctions` package.

end  #hide

md"""
## Cross-check against the actual SciPy library

The test above implements the *documented* formula.  To make sure the documentation matches
the code, the following test calls `scipy.special.sph_harm_y` itself, via `PythonCall`.  It is
tagged `:python` (and `:skipci`) so that it runs only when explicitly requested — e.g., with
`julia --project=. scripts/test.jl :python` — or in the scheduled CI workflow, which installs
SciPy through `CondaPkg`.
"""

@testitem "SciPy cross-check" tags=[:python, :skipci] setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide
using PythonCall
special = pyimport("scipy.special")
scipy_version = pyconvert(String, pyimport("scipy").__version__)
@info "Cross-checking against SciPy version $scipy_version"
ϵₐ = 100eps()
ϵᵣ = 1000eps()
for (θ, ϕ) ∈ θϕrange(Float64, 7)
    for (ℓ, m) ∈ ℓmrange(6)
        Y_scipy = pyconvert(ComplexF64, pybuiltins.complex(special.sph_harm_y(ℓ, m, θ, ϕ)))
        @test Y_scipy ≈ ConventionsUtilities.Y(ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
end  #hide
