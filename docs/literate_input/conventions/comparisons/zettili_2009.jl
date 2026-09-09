md"""
# Zettili (2009)

!!! info "Summary"
    Zettili's definitions of the spherical harmonics and of Wigner's ``d`` and ``D``
    matrices agree with the definitions used in the `SphericalFunctions` package, and his
    rotation law for the spherical harmonics is ours.

[Zettili_2009](@citet) is a relatively recent textbook that seems to be gaining popularity.
(Note that there is a 3rd edition from 2022, but I do not have access to it; all the
references here are to the 2nd edition from 2009.)

## Coordinates and operators

In Appendix B.1, we find that the spherical coordinates are related to Cartesian coordinates
in the usual (physicist's) way.  Equation (5.132) gives the angular-momentum operator
```math
\hat{L}_z = -i \hbar \frac{\partial}{\partial φ},
```
which agrees with [our expression](@ref "``L`` operators in spherical coordinates").  This is
followed by equation (5.134):
```math
\hat{L}_{\pm}
=
\hat{L}_x \pm i \hat{L}_y
=
\pm \hbar e^{\pm iφ} \left(
  \frac{\partial}{\partial θ}
  \pm i \cot θ \frac{\partial}{\partial φ}
\right),
```
which also agrees with [our results.](@ref "``L_{\pm}`` operators in spherical coordinates")

## Spherical harmonics

Equation (5.180) gives the spherical harmonics as
```math
Y_{l, m}(θ, φ)
=
\frac{(-1)^l}{2^l l!}
\sqrt{\frac{2l+1}{4π} \frac{(l+m)!}{(l-m)!}}
e^{imφ}
\frac{1}{\sin^m θ}
\frac{d^{l-m}}{d(\cos θ)^{l-m}}
(\sin θ)^{2l},
```
which is the Condon–Shortley form, valid for all ``m``.

## Rotations and Wigner's ``D`` matrix

Section 7.2.1 denotes by ``\hat{R}_z(δ ϕ)`` the

> rotation of the coordinates of a *spinless* particle over an *infinitesimal* angle ``δ
> ϕ`` about the ``z``-axis

and shows its action [Eq. (7.16)]
```math
\hat{R}_z (δ ϕ) ψ(r, θ, ϕ)
=
ψ(r, θ, ϕ - δ ϕ).
```

> We may generalize this relation to a rotation of angle ``δ ϕ`` about an arbitrary axis
> whose direction is given by the unit vector ``\vec{n}``:

```math
\hat{R}(δ ϕ)
=
1 - \frac{i}{\hbar} δ ϕ \vec{n} \cdot \hat{\vec{L}}.
```
Both of these are precisely [our ``U(𝐑)``](@ref summary_wigner_D): the field is rotated,
so its argument is rotated inversely.  This extends to finite rotation by defining the
operator [Eq. (7.48)]
```math
\hat{R}(α, β, γ)
=
e^{-iα J_z / \hbar} e^{-iβ J_y / \hbar} e^{-iγ J_z / \hbar}.
```
Equation (7.52) then defines
```math
D^{(j)}_{m', m}(α, β, γ)
=
\langle j, m' | \hat{R}(α, β, γ) | j, m \rangle,
```
so that [Eq. (7.54)]
```math
D^{(j)}_{m', m}(α, β, γ)
=
e^{-i (m' α + m γ)} d^{(j)}_{m', m}(β),
```
where [Eq. (7.55)]
```math
d^{(j)}_{m', m}(β)
=
\langle j, m' | e^{-iβ J_y / \hbar} | j, m \rangle.
```
The explicit expression for ``d`` is [Eq. (7.56)]
```math
d^{(j)}_{m', m}(β)
=
\sum_k (-1)^{k+m'-m}
\frac{\sqrt{(j+m)!(j-m)!(j+m')!(j-m')!}}
{(j-m'-k)!(j+m-k)!(k+m'-m)!k!}
\left(\cos\frac{β}{2}\right)^{2j+m-m'-2k}
\left(\sin\frac{β}{2}\right)^{m'-m+2k}.
```
In Sec. 7.2.6, we find that if the operator ``\hat{R}(α, β, γ)`` rotates a vector pointing
in the ``(θ, ϕ)`` direction to a vector pointing in the ``(θ', ϕ')`` direction, then the
spherical harmonics transform as [Eq. (7.70)]
```math
Y_{ℓ, m}^\ast (θ', ϕ')
=
\sum_{m'} D^{(ℓ)}_{m, m'}(α, β, γ) Y_{ℓ, m'}^\ast (θ, ϕ).
```
Taking the complex conjugate, this is the second form of [our rotation law](@ref
summary_spherical_harmonics), ``Y_{ℓ,m}(𝐑\,𝐐) = \sum_{m'} \overline{𝔇_{m,m'}(𝐑)}\,
Y_{ℓ,m'}(𝐐)``, evaluated at a rotated *point*.

## Implementing formulas

We begin by writing code that implements the formulas from Zettili.  We encapsulate the
formulas in a module so that we can test them against the `SphericalFunctions` package.
"""

using TestItems: @testitem  #hide
@testitem "Zettili conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide

module Zettili
#+

# We'll use some predefined utilities to make the code look more like the equations.
import ..ConventionsUtilities: 𝒾, ❗, dʲsin²ᵏθdcosθʲ
#+

# Equation (5.180).  We capture the floating-point type `T` to ensure that we don't lose
# precision when converting π and the factorials to floating-point numbers.
function Y(l, m, θ::T, φ::T) where {T<:Real}
    (-1)^l / T(2^l * (l)❗) * √T((2l+1) * (l+m)❗ / (4big(π) * (l-m)❗)) *
    exp(𝒾 * m * φ) * (1 / sin(θ)^m) * dʲsin²ᵏθdcosθʲ(j=l-m, k=l, θ=θ)
end
#+

# Equation (7.56).  The sum runs over all ``k`` for which the factorials have non-negative
# arguments, ``\max(0, m-m') \leq k \leq \min(j-m', j+m)``.
function d(j, m′, m, β::T) where {T<:Real}
    sum(
        (-1)^(k+m′-m) * T(
            √((j+m)❗ * (j-m)❗ * (j+m′)❗ * (j-m′)❗)
            / ((j-m′-k)❗ * (j+m-k)❗ * (k+m′-m)❗ * (k)❗)
        )
        * cos(β/2)^(2j+m-m′-2k) * sin(β/2)^(m′-m+2k)
        for k ∈ max(0, m-m′):min(j-m′, j+m);
        init=zero(T)
    )
end
#+

# Equation (7.54):
function D(j, m′, m, α, β, γ)
    exp(-𝒾 * (m′ * α + m * γ)) * d(j, m′, m, β)
end
#+

end  # module Zettili
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
# normalization differences, which are the most likely source of error.  For the same reason
# we use a modest grid of Euler angles.
αβγs = αβγrange(Float64, 5)
#+

# First, the spherical harmonics agree with ours.  The general formula has a factor of
# ``1/\sin^m θ``, so we avoid the poles.
for (θ, ϕ) ∈ θϕrange(; avoid_poles=ϵₐ/40)
    for (ℓ, m) ∈ ℓmrange(ℓₘₐₓ)
        @test Zettili.Y(ℓ, m, θ, ϕ) ≈ ConventionsUtilities.Y(ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# The ``d`` and ``D`` matrices agree with ours:
for β ∈ βrange()
    for (j, m′, m) ∈ ℓm′mrange(ℓₘₐₓ)
        @test Zettili.d(j, m′, m, β) ≈ ConventionsUtilities.d(j, m′, m, β) atol=ϵₐ rtol=ϵᵣ
    end
end
for (α, β, γ) ∈ αβγs
    for (j, m′, m) ∈ ℓm′mrange(ℓₘₐₓ)
        @test Zettili.D(j, m′, m, α, β, γ) ≈ ConventionsUtilities.D(j, m′, m, α, β, γ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Finally, the rotation law of Eq. (7.70).  We rotate the unit vector ``𝐧(θ, ϕ)`` by the
# rotor ``𝐑_{α,β,γ}`` to obtain ``𝐧(θ', ϕ')``, and check that the spherical harmonics at
# the rotated point are given by the stated combination of the harmonics at the original
# point.  We use our own spherical harmonics on both sides, since we have just shown that
# they agree with Zettili's, and avoid the poles so that ``ϕ'`` is well defined.
import Quaternionic: from_euler_angles, from_spherical_coordinates, imz, components
for (α, β, γ) ∈ αβγrange(Float64, 3)
    R = from_euler_angles(α, β, γ)
    for (θ, ϕ) ∈ θϕrange(Float64, 3; avoid_poles=1e-3)
        𝐧 = from_spherical_coordinates(θ, ϕ) * imz * conj(from_spherical_coordinates(θ, ϕ))
        𝐧′ = R * 𝐧 * conj(R)
        _, x, y, z = components(𝐧′)
        θ′, ϕ′ = atan(hypot(x, y), z), atan(y, x)  # well conditioned near the poles
        for (ℓ, m) ∈ ℓmrange(ℓₘₐₓ)
            @test conj(ConventionsUtilities.Y(ℓ, m, θ′, ϕ′)) ≈ sum(
                Zettili.D(ℓ, m, m′, α, β, γ) * conj(ConventionsUtilities.Y(ℓ, m′, θ, ϕ))
                for m′ ∈ -ℓ:ℓ
            ) atol=ϵₐ rtol=ϵᵣ
        end
    end
end
#+

# These successful tests show that Zettili's spherical harmonics, ``d`` and ``D`` matrices,
# and rotation law agree with the corresponding definitions in the `SphericalFunctions`
# package.

end  #hide
