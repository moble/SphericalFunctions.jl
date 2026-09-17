md"""
# Sakurai (1994)

!!! info "Summary"
    Sakurai's definitions of the spherical harmonics and of Wigner's ``d`` and ``𝒟``
    matrices agree with the definitions used in the `SphericalFunctions` package.

[Sakurai_1994](@citet) is the standard graduate textbook on quantum mechanics, and its
Chapter 3 is one of the clearest treatments of rotations.  (Page and equation numbers below
refer to the Revised Edition of 1994.)  Sakurai's conventions are the ones adopted on our
[conventions pages](@ref summary_wigner_D), so this page mostly serves to make that
agreement explicit and to test it.

## Rotation operators

Sakurai is careful to state that a rotation "affects the physical system itself, [...] while
the coordinate axes remain *unchanged*" (p. 154), and defines the rotation operator by
```math
|α\rangle_R = 𝒟(R) |α\rangle,
```
"where ``|α\rangle_R`` and ``|α\rangle`` stand for the kets of the rotated and original
system, respectively" (p. 156).  For an infinitesimal rotation, Eq. (3.1.15) gives
```math
𝒟(\hat{𝐧}, dϕ) = 1 - i \left( \frac{𝐉 \cdot \hat{𝐧}}{\hbar} \right) dϕ,
```
which is exactly [our ``U(e^{ϵ𝐮/2}) = 1 - iϵ L_𝐮``](@ref summary_wigner_D).  His Euler
angles (Sec. 3.3, p. 173) are defined as ours, and the matrix elements of the rotation
operator define the ``𝒟`` matrix [Eq. (3.5.42), p. 192]
```math
𝒟^{(j)}_{m',m}(R) = \langle j, m' | \exp\left( \frac{-i 𝐉 \cdot \hat{𝐧} ϕ}{\hbar} \right) | j, m \rangle,
```
which obeys the representation property [Eq. (3.5.46)]
```math
𝒟^{(j)}_{m'',m}(R_1 R_2) = \sum_{m'} 𝒟^{(j)}_{m'',m'}(R_1)\, 𝒟^{(j)}_{m',m}(R_2).
```
In terms of Euler angles [Eqs. (3.5.50)–(3.5.51), p. 194],
```math
𝒟^{(j)}_{m',m}(α, β, γ)
=
\langle j, m' | \exp\left(\frac{-iJ_z α}{\hbar}\right)
  \exp\left(\frac{-iJ_y β}{\hbar}\right)
  \exp\left(\frac{-iJ_z γ}{\hbar}\right) | j, m \rangle
=
e^{-i(m'α + mγ)}\, d^{(j)}_{m',m}(β),
\qquad
d^{(j)}_{m',m}(β) = \langle j, m' | \exp\left(\frac{-iJ_y β}{\hbar}\right) | j, m \rangle,
```
and Wigner's formula for ``d`` is given as [Eq. (3.8.33), p. 223]
```math
d^{(j)}_{m',m}(β)
=
\sum_k (-1)^{k-m+m'}
\frac{\sqrt{(j+m)!\,(j-m)!\,(j+m')!\,(j-m')!}}
     {(j+m-k)!\,k!\,(j-k-m')!\,(k-m+m')!}
\left(\cos\frac{β}{2}\right)^{2j-2k+m-m'}
\left(\sin\frac{β}{2}\right)^{2k-m+m'}.
```

## Spherical harmonics

Sakurai writes the spherical harmonics as ``Y_ℓ^m(θ, ϕ)`` (note the upper index ``m``) and
relates them to the ``𝒟`` matrix by [Eq. (3.6.51), p. 203]
```math
𝒟^{(ℓ)}_{m,0}(α, β, γ=0)
=
\sqrt{\frac{4π}{2ℓ+1}}\; \left. Y_ℓ^{m\ast}(θ, ϕ) \right|_{θ=β,\, ϕ=α},
```
which is [our relation](@ref summary_spherical_harmonics).  Appendix A gives the explicit
forms [Eq. (A.5.7), p. 451]
```math
\begin{aligned}
Y_0^0 &= \frac{1}{\sqrt{4π}}, &
Y_1^{\pm 1} &= \mp\sqrt{\frac{3}{8π}} \sin θ\, e^{\pm iϕ}, &
Y_1^0 &= \sqrt{\frac{3}{4π}} \cos θ, \\
Y_2^{\pm 2} &= \sqrt{\frac{15}{32π}} \sin^2 θ\, e^{\pm 2iϕ}, &
Y_2^{\pm 1} &= \mp\sqrt{\frac{15}{8π}} \sin θ \cos θ\, e^{\pm iϕ}, &
Y_2^0 &= \sqrt{\frac{5}{16π}} \left(3\cos^2 θ - 1\right),
\end{aligned}
```
together with the conjugation relation [Eq. (A.5.7b)]
```math
Y_ℓ^{-m}(θ, ϕ) = (-1)^m\, \left[Y_ℓ^m(θ, ϕ)\right]^\ast.
```

## Implementing formulas

We begin by writing code that implements the formulas from Sakurai.  We encapsulate the
formulas in a module so that we can test them against the `SphericalFunctions` package.
"""

using TestItems: @testitem  #hide
@testitem "Sakurai conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide

module Sakurai
#+

# We'll use some predefined utilities to make the code look more like the equations.
import ..ConventionsUtilities: 𝒾, ❗
#+

# Wigner's formula for ``d``, Eq. (3.8.33).  The sum runs over all ``k`` for which the
# factorials have non-negative arguments, ``\max(0, m-m') \leq k \leq \min(j+m, j-m')``.
function d(j, m′, m, β::T) where {T<:Real}
    sum(
        (-1)^(k-m+m′) * T(
            √((j+m)❗ * (j-m)❗ * (j+m′)❗ * (j-m′)❗)
            / ((j+m-k)❗ * (k)❗ * (j-k-m′)❗ * (k-m+m′)❗)
        )
        * cos(β/2)^(2j-2k+m-m′) * sin(β/2)^(2k-m+m′)
        for k ∈ max(0, m-m′):min(j+m, j-m′);
        init=zero(T)
    )
end
#+

# The ``𝒟`` matrix in terms of Euler angles, Eq. (3.5.50):
function 𝒟(j, m′, m, α, β, γ)
    exp(-𝒾 * (m′ * α + m * γ)) * d(j, m′, m, β)
end
#+

# The spherical harmonics, obtained by solving Eq. (3.6.51) for ``Y_ℓ^m``:
function Y(ℓ, m, θ::T, ϕ::T) where {T<:Real}
    conj(√((2ℓ+1) / (4T(π))) * 𝒟(ℓ, m, 0, ϕ, θ, zero(T)))
end
#+

# The explicit formulas of Eq. (A.5.7):
Y₀⁰(θ, ϕ) = 1 / √(4π)
Y₁⁻¹(θ, ϕ) = +√(3/(8π)) * sin(θ) * exp(-𝒾*ϕ)
Y₁⁰(θ, ϕ) = √(3/(4π)) * cos(θ)
Y₁⁺¹(θ, ϕ) = -√(3/(8π)) * sin(θ) * exp(+𝒾*ϕ)
Y₂⁻²(θ, ϕ) = √(15/(32π)) * sin(θ)^2 * exp(-2𝒾*ϕ)
Y₂⁻¹(θ, ϕ) = +√(15/(8π)) * sin(θ) * cos(θ) * exp(-𝒾*ϕ)
Y₂⁰(θ, ϕ) = √(5/(16π)) * (3cos(θ)^2 - 1)
Y₂⁺¹(θ, ϕ) = -√(15/(8π)) * sin(θ) * cos(θ) * exp(+𝒾*ϕ)
Y₂⁺²(θ, ϕ) = √(15/(32π)) * sin(θ)^2 * exp(+2𝒾*ϕ)
#+

end  # module Sakurai
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
# normalization differences, which are the most likely source of error.  For the same reason
# we use a modest grid of Euler angles.
αβγs = αβγrange(Float64, 5)
#+

# First, Sakurai's own explicit spherical harmonics and conjugation relation are consistent
# with his general expression via ``𝒟``:
for (θ, ϕ) ∈ θϕrange()
    @test Sakurai.Y₀⁰(θ, ϕ) ≈ Sakurai.Y(0, 0, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Sakurai.Y₁⁻¹(θ, ϕ) ≈ Sakurai.Y(1, -1, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Sakurai.Y₁⁰(θ, ϕ) ≈ Sakurai.Y(1, 0, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Sakurai.Y₁⁺¹(θ, ϕ) ≈ Sakurai.Y(1, 1, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Sakurai.Y₂⁻²(θ, ϕ) ≈ Sakurai.Y(2, -2, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Sakurai.Y₂⁻¹(θ, ϕ) ≈ Sakurai.Y(2, -1, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Sakurai.Y₂⁰(θ, ϕ) ≈ Sakurai.Y(2, 0, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Sakurai.Y₂⁺¹(θ, ϕ) ≈ Sakurai.Y(2, 1, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Sakurai.Y₂⁺²(θ, ϕ) ≈ Sakurai.Y(2, 2, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    for ℓ ∈ 0:ℓₘₐₓ, m ∈ -ℓ:-1
        @test Sakurai.Y(ℓ, m, θ, ϕ) ≈ (-1)^m * conj(Sakurai.Y(ℓ, -m, θ, ϕ)) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Now the spherical harmonics agree with ours:
for (θ, ϕ) ∈ θϕrange()
    for (ℓ, m) ∈ ℓmrange(ℓₘₐₓ)
        @test Sakurai.Y(ℓ, m, θ, ϕ) ≈ ConventionsUtilities.Y(ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# The ``d`` matrix agrees with ours:
for β ∈ βrange()
    for (j, m′, m) ∈ ℓm′mrange(ℓₘₐₓ)
        @test Sakurai.d(j, m′, m, β) ≈ ConventionsUtilities.d(j, m′, m, β) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# And the ``𝒟`` matrix agrees with ours — with no complex conjugation, transposition, or
# reordering of the Euler angles:
for (α, β, γ) ∈ αβγs
    for (j, m′, m) ∈ ℓm′mrange(ℓₘₐₓ)
        @test Sakurai.𝒟(j, m′, m, α, β, γ) ≈ ConventionsUtilities.D(j, m′, m, α, β, γ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# These successful tests show that Sakurai's spherical harmonics and Wigner ``d`` and ``𝒟``
# matrices agree with the corresponding functions defined by the `SphericalFunctions`
# package.

end  #hide
