md"""
# Torres del Castillo (2003)

!!! info "Summary"
    Torres del Castillo's definitions of the spin-weighted spherical harmonics and of
    Wigner's ``d`` and ``D`` matrices agree with the definitions used in the
    `SphericalFunctions` package.

[TorresDelCastillo_2003](@citet) is a monograph on 3-D spinors and spin-weighted functions,
and one of the few book-length treatments of the spin-weighted spherical harmonics.

## Rotations and Wigner's ``D`` matrix

He starts by defining a rotation ``ℛ`` as transforming a point ``x_i`` into another point
with coordinates ``x_i' = a_{ij}x_j``.  Under that rotation, any scalar function ``f``
transforms into another function ``f' = ℛ f`` defined by [Eq. (2.43)]
```math
f'\big(x_i\big) = f\big( a^{-1}_{ij} x_j \big).
```
In particular, ``f'(x'_i) = f(x_i)``.  This is exactly [our rotation operator
``U(𝐑)``](@ref summary_wigner_D), which rotates the field and hence the argument of the
function inversely.  He then defines Wigner's D-matrix to satisfy [Eq. (2.45)]
```math
ℛ Y_{l,m} = \sum_{m'} D^l_{m',m}(ℛ) Y_{l,m'}.
```
Including the arguments to the spherical harmonics, this becomes
```math
Y_{l,m}\big(ℛ^{-1} R_{θ, ϕ}\big)
=
\sum_{m'} D^l_{m',m}(ℛ) Y_{l,m'}\big(R_{θ, ϕ}\big),
```
which is [our rotation law](@ref summary_spherical_harmonics).  In this form, we have [Eq.
(2.46)]
```math
D^l_{m'',m}(ℛ_1 ℛ_2)
=
\sum_{m'} D^l_{m'',m'}(ℛ_1) D^l_{m',m}(ℛ_2).
```
He computes [Eq. (2.53)]
```math
D^l_{m',m}(ϕ, θ, \chi)
=
e^{-i m' ϕ} d^l_{m',m}(θ) e^{-i m \chi},
```
where the ``d`` matrix is given in the second equation below Eq. (2.53) by
```math
d^l_{m',m}(θ)
=
\sqrt{(l+m)!(l-m)!(l+m')!(l-m')!}
\sum_{k} \frac{
  (-1)^k
  (\sin \tfrac{1}{2} θ)^{m-m'+2k}
  (\cos \tfrac{1}{2} θ)^{2l-m+m'-2k}
} {
  k!(l+m'-k)!(l-m-k)!(m-m'+k)!
},
```
and the spin-weighted spherical harmonic is related to ``D`` in the equation following Eq.
(2.53) by
```math
{}_{s}Y_{j,m}(θ, ϕ)
=
(-1)^m
\sqrt{\frac{2j+1}{4π}}
d^j_{-m,s}(θ)
e^{i m ϕ}.
```
The phases in ``D`` are those of [our definition](@ref summary_wigner_D), and the last
expression is the integer-index form of [our definition of
``{}_sY_{ℓ,m}``](@ref summary_swsh), so we expect complete agreement.

## Implementing formulas

We begin by writing code that implements the formulas from Torres del Castillo.  We
encapsulate the formulas in a module so that we can test them against the
`SphericalFunctions` package.
"""

using TestItems: @testitem  #hide
@testitem "Torres del Castillo conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide

module TorresDelCastillo
#+

# We'll use some predefined utilities to make the code look more like the equations.
import ..ConventionsUtilities: 𝒾, ❗
#+

# The ``d`` matrix.  The sum runs over all ``k`` for which the factorials have non-negative
# arguments, ``\max(0, m'-m) \leq k \leq \min(l+m', l-m)``.
function d(l, m′, m, θ::T) where {T<:Real}
    √T((l+m)❗ * (l-m)❗ * (l+m′)❗ * (l-m′)❗) *
    sum(
        (-1)^k * sin(θ/2)^(m-m′+2k) * cos(θ/2)^(2l-m+m′-2k)
        / T((k)❗ * (l+m′-k)❗ * (l-m-k)❗ * (m-m′+k)❗)
        for k ∈ max(0, m′-m):min(l+m′, l-m);
        init=zero(T)
    )
end
#+

# Equation (2.53):
function D(l, m′, m, ϕ, θ, χ)
    exp(-𝒾 * m′ * ϕ) * d(l, m′, m, θ) * exp(-𝒾 * m * χ)
end
#+

# The spin-weighted spherical harmonics, from the equation following Eq. (2.53):
function Y(s, j, m, θ::T, ϕ::T) where {T<:Real}
    (-1)^m * √((2j+1) / (4T(π))) * d(j, -m, s, θ) * exp(𝒾 * m * ϕ)
end
#+

end  # module TorresDelCastillo
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
# and
sₘₐₓ = 2
#+
# because the formulas are slow, and this will be sufficient to sort out any sign or
# normalization differences, which are the most likely source of error.  For the same reason
# we use a modest grid of Euler angles.
αβγs = αβγrange(Float64, 5)
#+

# The ``d`` and ``D`` matrices agree with ours:
for β ∈ βrange()
    for (l, m′, m) ∈ ℓm′mrange(ℓₘₐₓ)
        @test TorresDelCastillo.d(l, m′, m, β) ≈ ConventionsUtilities.d(l, m′, m, β) atol=ϵₐ rtol=ϵᵣ
    end
end
for (α, β, γ) ∈ αβγs
    for (l, m′, m) ∈ ℓm′mrange(ℓₘₐₓ)
        @test TorresDelCastillo.D(l, m′, m, α, β, γ) ≈ ConventionsUtilities.D(l, m′, m, α, β, γ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# And the spin-weighted spherical harmonics agree with ours:
for (θ, ϕ) ∈ θϕrange()
    for (s, ℓ, m) ∈ sℓmrange(ℓₘₐₓ, sₘₐₓ)
        @test TorresDelCastillo.Y(s, ℓ, m, θ, ϕ) ≈ ConventionsUtilities.Y(s, ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# These successful tests show that Torres del Castillo's ``d`` and ``D`` matrices and
# spin-weighted spherical harmonics agree with the corresponding functions defined by the
# `SphericalFunctions` package.

end  #hide
