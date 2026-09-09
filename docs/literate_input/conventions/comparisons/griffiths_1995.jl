md"""
# Griffiths (1995)

!!! info "Summary"
    Griffiths' definition of the spherical harmonics agrees with the definition used in the
    `SphericalFunctions` package.

Griffiths' ["Introduction to Quantum Mechanics"](@cite Griffiths_1995) is probably the most
common introductory text used in undergraduate physics programs, so it would be useful to
compare.  (Equation and table numbers refer to the first edition, 1995.)  Griffiths uses the
standard (physicists') spherical coordinates, and writes the spherical harmonics as
``Y_ℓ^m``.

## Spherical harmonics

Equation (4.27) gives the associated Legendre function as
```math
P_{ℓ}^{m}(x)
=
(1-x^2)^{|m|/2} \left(\frac{d}{dx}\right)^{|m|} P_{ℓ}(x),
```
and (4.28) gives the Legendre polynomial as
```math
P_{ℓ}(x)
=
\frac{1}{2^ℓ ℓ!} \left(\frac{d}{dx}\right)^ℓ (x^2-1)^ℓ.
```
Then, (4.32) gives the spherical harmonics as
```math
Y_{ℓ}^{m}(θ, ϕ)
=
ϵ
\sqrt{\frac{2ℓ+1}{4π} \frac{(ℓ-|m|)!}{(ℓ+|m|)!}}
e^{imϕ} P_{ℓ}^{m}(\cos θ),
```
where ``ϵ = (-1)^m`` for ``m\geq 0`` and ``ϵ = 1`` for ``m\leq 0``.  Since ``P_ℓ^m``
depends only on ``|m|``, this places the Condon–Shortley phase entirely in ``ϵ``; the result
is the standard Condon–Shortley convention, so we expect agreement with [ours](@ref
summary_spherical_harmonics).  In Table 4.2, he explicitly lists the first few spherical
harmonics:
```math
\begin{aligned}
  Y_{0}^{0} &= \left(\frac{1}{4π}\right)^{1/2},\\
  Y_{1}^{0} &= \left(\frac{3}{4π}\right)^{1/2} \cos θ,\\
  Y_{1}^{\pm 1} &= \mp \left(\frac{3}{8π}\right)^{1/2} \sin θ e^{\pm iϕ},\\
  Y_{2}^{0} &= \left(\frac{5}{16π}\right)^{1/2} \left(3\cos^2θ - 1\right),\\
  Y_{2}^{\pm 1} &= \mp \left(\frac{15}{8π}\right)^{1/2} \sin θ \cos θ e^{\pm iϕ},\\
  Y_{2}^{\pm 2} &= \left(\frac{15}{32π}\right)^{1/2} \sin^2θ e^{\pm 2iϕ},\\
  Y_{3}^{0} &= \left(\frac{7}{16π}\right)^{1/2} \left(5\cos^3θ - 3\cos θ\right),\\
  Y_{3}^{\pm 1} &= \mp \left(\frac{21}{64π}\right)^{1/2} \sin θ \left(5\cos^2θ - 1\right) e^{\pm iϕ},\\
  Y_{3}^{\pm 2} &= \left(\frac{105}{32π}\right)^{1/2} \sin^2θ \cos θ e^{\pm 2iϕ},\\
  Y_{3}^{\pm 3} &= \mp \left(\frac{35}{64π}\right)^{1/2} \sin^3θ e^{\pm 3iϕ}.
\end{aligned}
```

## Angular-momentum operators

In Eqs. (4.127)—(4.129), he gives the angular-momentum operators in terms of spherical
coordinates:
```math
\begin{aligned}
L_x &= \frac{\hbar}{i} \left(
    -\sin ϕ \frac{\partial} {\partial θ}
    - \cos ϕ \cot θ \frac{\partial} {\partial ϕ}
\right), \\
L_y &= \frac{\hbar}{i} \left(
    \cos ϕ \frac{\partial} {\partial θ}
    - \sin ϕ \cot θ \frac{\partial} {\partial ϕ}
\right), \\
L_z &= -i \hbar \frac{\partial} {\partial ϕ},
\end{aligned}
```
which agree (up to the factor of ``\hbar``) with [ours](@ref summary_L_R_euler), since
``\hbar/i = -i\hbar``.  Griffiths does not discuss Wigner's ``D`` matrices.

## Implementing formulas

We begin by writing code that implements the formulas from Griffiths.  We encapsulate the
formulas in a module so that we can test them against the `SphericalFunctions` package.
"""

using TestItems: @testitem  #hide
@testitem "Griffiths conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide

module Griffiths
#+

# We'll use some predefined utilities to make the code look more like the equations,
# including `∂ⁿ`, which computes derivatives symbolically so that we can transcribe Eqs.
# (4.27) and (4.28) literally.
import ..ConventionsUtilities: 𝒾, ❗, ∂ⁿ
#+

# Equation (4.28):
function P(ℓ, x)
    1 / (2^ℓ * (ℓ)❗) * ∂ⁿ(x -> (x^2 - 1)^ℓ, ℓ)(x)
end
#+

# Equation (4.27):
function P(ℓ, m, x)
    (1 - x^2)^(abs(m)/2) * ∂ⁿ(x -> P(ℓ, x), abs(m))(x)
end
#+

# Equation (4.32).  We capture the floating-point type `T` to ensure that we don't lose
# precision when converting π and the factorials to floating-point numbers.
function Y(ℓ, m, θ::T, ϕ::T) where {T<:Real}
    ϵ = m ≥ 0 ? (-1)^m : 1
    ϵ * √T((2ℓ+1) * (ℓ-abs(m))❗ / (4big(π) * (ℓ+abs(m))❗)) * exp(𝒾 * m * ϕ) * T(P(ℓ, m, cos(θ)))
end
#+

# Table 4.2:
Y₀⁰(θ, ϕ) = √(1/(4π))
Y₁⁰(θ, ϕ) = √(3/(4π)) * cos(θ)
Y₁⁺¹(θ, ϕ) = -√(3/(8π)) * sin(θ) * exp(+𝒾*ϕ)
Y₁⁻¹(θ, ϕ) = +√(3/(8π)) * sin(θ) * exp(-𝒾*ϕ)
Y₂⁰(θ, ϕ) = √(5/(16π)) * (3cos(θ)^2 - 1)
Y₂⁺¹(θ, ϕ) = -√(15/(8π)) * sin(θ) * cos(θ) * exp(+𝒾*ϕ)
Y₂⁻¹(θ, ϕ) = +√(15/(8π)) * sin(θ) * cos(θ) * exp(-𝒾*ϕ)
Y₂⁺²(θ, ϕ) = √(15/(32π)) * sin(θ)^2 * exp(+2𝒾*ϕ)
Y₂⁻²(θ, ϕ) = √(15/(32π)) * sin(θ)^2 * exp(-2𝒾*ϕ)
Y₃⁰(θ, ϕ) = √(7/(16π)) * (5cos(θ)^3 - 3cos(θ))
Y₃⁺¹(θ, ϕ) = -√(21/(64π)) * sin(θ) * (5cos(θ)^2 - 1) * exp(+𝒾*ϕ)
Y₃⁻¹(θ, ϕ) = +√(21/(64π)) * sin(θ) * (5cos(θ)^2 - 1) * exp(-𝒾*ϕ)
Y₃⁺²(θ, ϕ) = √(105/(32π)) * sin(θ)^2 * cos(θ) * exp(+2𝒾*ϕ)
Y₃⁻²(θ, ϕ) = √(105/(32π)) * sin(θ)^2 * cos(θ) * exp(-2𝒾*ϕ)
Y₃⁺³(θ, ϕ) = -√(35/(64π)) * sin(θ)^3 * exp(+3𝒾*ϕ)
Y₃⁻³(θ, ϕ) = +√(35/(64π)) * sin(θ)^3 * exp(-3𝒾*ϕ)
#+

end  # module Griffiths
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

# First, Griffiths' Table 4.2 agrees with his general formula (4.32):
for (θ, ϕ) ∈ θϕrange(Float64, 7)
    @test Griffiths.Y₀⁰(θ, ϕ) ≈ Griffiths.Y(0, 0, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Griffiths.Y₁⁰(θ, ϕ) ≈ Griffiths.Y(1, 0, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Griffiths.Y₁⁺¹(θ, ϕ) ≈ Griffiths.Y(1, 1, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Griffiths.Y₁⁻¹(θ, ϕ) ≈ Griffiths.Y(1, -1, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Griffiths.Y₂⁰(θ, ϕ) ≈ Griffiths.Y(2, 0, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Griffiths.Y₂⁺¹(θ, ϕ) ≈ Griffiths.Y(2, 1, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Griffiths.Y₂⁻¹(θ, ϕ) ≈ Griffiths.Y(2, -1, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Griffiths.Y₂⁺²(θ, ϕ) ≈ Griffiths.Y(2, 2, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Griffiths.Y₂⁻²(θ, ϕ) ≈ Griffiths.Y(2, -2, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Griffiths.Y₃⁰(θ, ϕ) ≈ Griffiths.Y(3, 0, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Griffiths.Y₃⁺¹(θ, ϕ) ≈ Griffiths.Y(3, 1, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Griffiths.Y₃⁻¹(θ, ϕ) ≈ Griffiths.Y(3, -1, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Griffiths.Y₃⁺²(θ, ϕ) ≈ Griffiths.Y(3, 2, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Griffiths.Y₃⁻²(θ, ϕ) ≈ Griffiths.Y(3, -2, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Griffiths.Y₃⁺³(θ, ϕ) ≈ Griffiths.Y(3, 3, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Griffiths.Y₃⁻³(θ, ϕ) ≈ Griffiths.Y(3, -3, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
end
#+

# Now, the general formula agrees with ours:
for (θ, ϕ) ∈ θϕrange(Float64, 7)
    for (ℓ, m) ∈ ℓmrange(ℓₘₐₓ)
        @test Griffiths.Y(ℓ, m, θ, ϕ) ≈ ConventionsUtilities.Y(ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# This successful test shows that the spherical harmonics defined by Griffiths agree with
# the spherical harmonics defined by the `SphericalFunctions` package.

end  #hide
