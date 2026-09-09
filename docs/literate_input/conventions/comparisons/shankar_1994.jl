md"""
# Shankar (1994)

!!! info "Summary"
    Shankar's definition of the spherical harmonics, and the matrix elements of his rotation
    operator ``U[R(α, β, γ)]``, agree with the spherical harmonics and Wigner ``𝔇`` matrices
    used in the `SphericalFunctions` package.

[Shankar_1994](@citet) is a widely used graduate textbook on quantum mechanics.  Shankar
uses the standard (physicists') spherical coordinates, and writes the spherical harmonics as
``Y_ℓ^m``.

## Angular-momentum operators

Just below Eq. (12.5.27), Shankar gives the angular-momentum operators in spherical
coordinates as
```math
\begin{aligned}
L_x &= i \hbar \left(
    \sin ϕ \frac{\partial} {\partial θ}
    + \cos ϕ \cot θ \frac{\partial} {\partial ϕ}
\right),
\\
L_y &= i \hbar \left(
    -\cos ϕ \frac{\partial} {\partial θ}
    + \sin ϕ \cot θ \frac{\partial} {\partial ϕ}
\right),
\\
L_z &= -i \hbar \frac{\partial} {\partial ϕ},
\end{aligned}
```
which agree (up to the factor of ``\hbar``) with [our expressions](@ref summary_L_R_euler).

## Spherical harmonics

Equation (12.5.35) gives the spherical harmonics for ``m \geq 0`` as
```math
Y_{ℓ}^{m}(θ, ϕ)
=
(-1)^ℓ
\left[ \frac{(2ℓ+1)!}{4π} \right]^{1/2}
\frac{1}{2^ℓ ℓ!}
\left[ \frac{(ℓ+m)!}{(2ℓ)!(ℓ-m)!} \right]^{1/2}
e^{i m ϕ}
(\sin θ)^{-m}
\frac{d^{ℓ-m}}{d(\cos θ)^{ℓ-m}}
(\sin θ)^{2ℓ},
```
with Eq. (12.5.40) extending this to negative ``m``,
```math
Y_{ℓ}^{-m}(θ, ϕ)
=
(-1)^m \left( Y_{ℓ}^{m}(θ, ϕ) \right)^\ast.
```
Equation (12.5.39) lists the first few explicitly:
```math
\begin{aligned}
Y_0^0 &= \frac{1}{\sqrt{4π}}, &
Y_1^{\pm 1} &= \mp\sqrt{\frac{3}{8π}} \sin θ\, e^{\pm iϕ}, &
Y_1^0 &= \sqrt{\frac{3}{4π}} \cos θ, \\
Y_2^{\pm 2} &= \sqrt{\frac{15}{32π}} \sin^2 θ\, e^{\pm 2iϕ}, &
Y_2^{\pm 1} &= \mp\sqrt{\frac{15}{8π}} \sin θ \cos θ\, e^{\pm iϕ}, &
Y_2^0 &= \sqrt{\frac{5}{16π}} \left(3\cos^2 θ - 1\right).
\end{aligned}
```
This is the Condon–Shortley form, so we expect agreement with
[ours](@ref summary_spherical_harmonics).

## Rotation operator

In Exercise 12.5.7, the rotation operator is defined by
```math
U\left[ R(α, β, γ) \right]
=
e^{-i α J_z/\hbar}
e^{-i β J_y/\hbar}
e^{-i γ J_z/\hbar}.
```
Shankar never actually uses notation like ``D^{(j)}_{m', m}``, but he does talk about
``\langle j, m' | D^{(j)}\left[ R(α, β, γ) \right] | j, m \rangle``, the matrix elements of
``U`` restricted to the states ``|j, m\rangle`` for a given ``j``.  This is precisely [our
definition of ``𝔇``](@ref summary_wigner_D), so we expect agreement.  To test it without
relying on any further formula from Shankar, we construct the matrices of ``J_z`` and
``J_y = (J_+ - J_-)/2i`` in the ``|j, m\rangle`` basis from the standard ladder relations
(which Shankar derives in Sec. 12.5), exponentiate them numerically, and take the matrix
elements of ``U``.

## Implementing formulas

We begin by writing code that implements the formulas from Shankar.  We encapsulate the
formulas in a module so that we can test them against the `SphericalFunctions` package.
"""

using TestItems: @testitem  #hide
@testitem "Shankar conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide

module Shankar
#+

# We'll use some predefined utilities to make the code look more like the equations, and
# the matrix exponential from `LinearAlgebra` for the rotation operator.
import ..ConventionsUtilities: 𝒾, ❗, dʲsin²ᵏθdcosθʲ
import LinearAlgebra: exp, Diagonal
#+

# Equation (12.5.35) for ``m \geq 0``, and Eq. (12.5.40) for ``m < 0``.  We capture the
# floating-point type `T` to ensure that we don't lose precision when converting π and the
# factorials to floating-point numbers.
function Y(ℓ, m, θ::T, ϕ::T) where {T<:Real}
    if m < 0
        return (-1)^m * conj(Y(ℓ, -m, θ, ϕ))
    end
    (-1)^ℓ * √T((2ℓ+1)❗ / (4big(π))) * (1 / T(2^ℓ * (ℓ)❗)) * √T((ℓ+m)❗ / ((2ℓ)❗ * (ℓ-m)❗)) *
    exp(𝒾 * m * ϕ) * sin(θ)^(-m) * dʲsin²ᵏθdcosθʲ(j=ℓ-m, k=ℓ, θ=θ)
end
#+

# The explicit formulas of Eq. (12.5.39):
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

# The matrices of ``J_z``, ``J_\pm``, and ``J_y`` in the ``|j, m\rangle`` basis (with
# ``\hbar = 1``), ordered so that row and column `k` correspond to ``m = -j + k - 1``:
function Jz(j)
    Diagonal([m for m ∈ -j:j])
end
function J₊(j)
    M = zeros(2j+1, 2j+1)
    for (k, m) ∈ enumerate(-j:j-1)
        M[k+1, k] = √((j-m) * (j+m+1))  # ⟨j, m+1| J₊ |j, m⟩
    end
    M
end
J₋(j) = J₊(j)'
Jy(j) = (J₊(j) - J₋(j)) / 2𝒾
#+

# The rotation operator of Exercise 12.5.7, and its matrix element ``\langle j, m' | U | j,
# m \rangle``:
function U(j, α, β, γ)
    exp(-𝒾 * α * Matrix(Jz(j))) * exp(-𝒾 * β * Matrix(Jy(j))) * exp(-𝒾 * γ * Matrix(Jz(j)))
end
function D(j, m′, m, α, β, γ)
    U(j, α, β, γ)[m′ + j + 1, m + j + 1]
end
#+

end  # module Shankar
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

# First, the explicit formulas of Eq. (12.5.39) agree with the general formula (12.5.35).
# The general formula has a factor of ``1/\sin^m θ``, so we avoid the poles.
for (θ, ϕ) ∈ θϕrange(; avoid_poles=ϵₐ/40)
    @test Shankar.Y₀⁰(θ, ϕ) ≈ Shankar.Y(0, 0, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Shankar.Y₁⁻¹(θ, ϕ) ≈ Shankar.Y(1, -1, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Shankar.Y₁⁰(θ, ϕ) ≈ Shankar.Y(1, 0, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Shankar.Y₁⁺¹(θ, ϕ) ≈ Shankar.Y(1, 1, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Shankar.Y₂⁻²(θ, ϕ) ≈ Shankar.Y(2, -2, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Shankar.Y₂⁻¹(θ, ϕ) ≈ Shankar.Y(2, -1, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Shankar.Y₂⁰(θ, ϕ) ≈ Shankar.Y(2, 0, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Shankar.Y₂⁺¹(θ, ϕ) ≈ Shankar.Y(2, 1, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test Shankar.Y₂⁺²(θ, ϕ) ≈ Shankar.Y(2, 2, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
end
#+

# Next, the spherical harmonics agree with ours:
for (θ, ϕ) ∈ θϕrange(; avoid_poles=ϵₐ/40)
    for (ℓ, m) ∈ ℓmrange(ℓₘₐₓ)
        @test Shankar.Y(ℓ, m, θ, ϕ) ≈ ConventionsUtilities.Y(ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Finally, the matrix elements of Shankar's rotation operator agree with our ``𝔇``:
for (α, β, γ) ∈ αβγs
    for (j, m′, m) ∈ ℓm′mrange(ℓₘₐₓ)
        @test Shankar.D(j, m′, m, α, β, γ) ≈ ConventionsUtilities.D(j, m′, m, α, β, γ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# These successful tests show that Shankar's spherical harmonics and rotation operator
# agree with the corresponding functions defined by the `SphericalFunctions` package.

end  #hide
