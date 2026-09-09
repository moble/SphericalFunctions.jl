md"""
# Thorne (1980)

!!! info "Summary"
    Thorne's definition of the scalar spherical harmonics agrees with the definition used
    in the `SphericalFunctions` package.

[Thorne_1980](@citet) is the standard reference for multipole expansions of gravitational
radiation, and its conventions for the spherical harmonics are inherited by much of the
gravitational-wave literature.  Thorne uses the standard (physicists') spherical
coordinates ``(θ, ϕ)``, and writes the indices of the spherical harmonics as superscripts:
``Y^{ℓ m}``.

Thorne's Eq. (2.7) gives the scalar spherical harmonics for ``m \geq 0`` as
```math
Y^{ℓ m}(θ, ϕ)
=
C^{ℓ m}\, \left(e^{iϕ} \sin θ\right)^m
\sum_{j=0}^{\lfloor (ℓ-m)/2 \rfloor} a^{ℓ m j}\, (\cos θ)^{ℓ-m-2j},
```
with the coefficients given in Eq. (2.8),
```math
C^{ℓ m}
=
(-1)^m \left[ \frac{(2ℓ+1)\,(ℓ-m)!}{4π\,(ℓ+m)!} \right]^{1/2},
\qquad
a^{ℓ m j}
=
\frac{(-1)^j}{2^ℓ\, j!\, (ℓ-j)!}\, \frac{(2ℓ-2j)!}{(ℓ-m-2j)!},
```
and the harmonics with negative ``m`` are defined through the conjugation relation of Eq.
(2.9b),
```math
Y^{ℓ, -m} = (-1)^m\, \overline{Y^{ℓ m}}.
```
The factor ``(-1)^m`` in ``C^{ℓ m}`` is the Condon–Shortley phase, so we expect agreement
with [our spherical harmonics](@ref summary_spherical_harmonics).  Thorne's paper also
defines the "pure-spin" vector and tensor harmonics built from these scalar harmonics, which
are outside the scope of this package.

## Implementing formulas

We begin by writing code that implements the formulas from Thorne.  We encapsulate the
formulas in a module so that we can test them against the `SphericalFunctions` package.
"""

using TestItems: @testitem  #hide
@testitem "Thorne conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide

module Thorne
#+

# We'll use some predefined utilities to make the code look more like the equations.
import ..ConventionsUtilities: 𝒾, ❗
#+

# The coefficients of Eq. (2.8).  We capture the floating-point type `T` to ensure that we
# don't lose precision when converting π and the factorials to floating-point numbers.
function C(ℓ, m, ::Type{T}) where {T}
    (-1)^m * √T(((2ℓ+1) * (ℓ-m)❗) / (4big(π) * (ℓ+m)❗))
end
function a(ℓ, m, j, ::Type{T}) where {T}
    T((-1)^j / (2^ℓ * (j)❗ * (ℓ-j)❗) * ((2ℓ-2j)❗ / (ℓ-m-2j)❗))
end
#+

# Equation (2.7) then gives the spherical harmonics for ``m \geq 0``, and Eq. (2.9b) extends
# them to negative ``m``:
function Y(ℓ, m, θ::T, ϕ::T) where {T<:Real}
    if m < 0
        return (-1)^m * conj(Y(ℓ, -m, θ, ϕ))
    end
    sinθ, cosθ = sincos(θ)
    C(ℓ, m, T) * (exp(𝒾 * ϕ) * sinθ)^m * sum(
        a(ℓ, m, j, T) * cosθ^(ℓ-m-2j)
        for j ∈ 0:(ℓ-m)÷2;
        init=zero(T)
    )
end
#+

end  # module Thorne
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
ℓₘₐₓ = 5
#+
# because the formulas are slow, and this will be sufficient to sort out any sign or
# normalization differences, which are the most likely source of error.

# First, we check that Thorne's conjugation relation (2.9b) is consistent with the explicit
# formula applied directly to negative ``m`` — that is, that the ``(-1)^m`` in ``C^{ℓ m}``
# and the ``e^{imϕ}`` factor really do reproduce the relation.  (The formula for negative
# ``m`` is only meaningful when ``(ℓ - |m| - 2j)!`` is defined, so we compare with the
# ``m>0`` form via the relation itself.)
for (θ, ϕ) ∈ θϕrange()
    for ℓ ∈ 0:ℓₘₐₓ
        for m ∈ 1:ℓ
            @test Thorne.Y(ℓ, -m, θ, ϕ) ≈ (-1)^m * conj(Thorne.Y(ℓ, m, θ, ϕ)) atol=ϵₐ rtol=ϵᵣ
        end
    end
end
#+

# Now we compare to the `SphericalFunctions` package:
for (θ, ϕ) ∈ θϕrange()
    for (ℓ, m) ∈ ℓmrange(ℓₘₐₓ)
        @test Thorne.Y(ℓ, m, θ, ϕ) ≈ ConventionsUtilities.Y(ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# This successful test shows that the spherical harmonics defined by Thorne agree with the
# spherical harmonics defined by the `SphericalFunctions` package.

end  #hide
