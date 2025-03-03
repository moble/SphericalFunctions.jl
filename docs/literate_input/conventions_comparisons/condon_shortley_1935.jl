md"""
# Condon-Shortley (1935)

!!! info "Summary"
    Condon and Shortley's definition of the spherical harmonics agrees with the definition
    used in the `SphericalFunctions` package.

[Condon and Shortley's "The Theory Of Atomic Spectra"](@cite CondonShortley_1935) is the
standard reference for the "Condon-Shortley phase convention".  Though some references are
not very clear about precisely what they mean by that phrase, it seems clear that the
original meaning revolved around the idea that the angular-momentum raising and lowering
operators have eigenvalues that are *real and positive* when acting on the spherical
harmonics.  Specifically, they discuss the phase ambiguity of the eigenfunction ``\psi`` —
which includes spherical harmonics indexed by ``j`` and ``m`` for the angular part — in
section 3³ (page 48).  This culminates in Eq. (3) of that section, which is as explicit as
they get:
```math
\left( J_x \pm i J_y \right) \psi(\gamma j m)
=
\hbar \sqrt{(j \mp m)(j \pm m + 1)} \psi(\gamma j m \pm 1).
```
This eliminates any *relative* phase ambiguity between modes with neighboring ``m`` values,
and specifically determines what factors of ``(-1)^m`` should be included in the definition
of the spherical harmonics.

To avoid re-introducing ambiguity, we can just look at the actual spherical harmonics they
define.  The method we use here is as direct and explicit as possible.  In particular,
Condon and Shortley provide a formula for the φ=0 part in terms of iterated derivatives of a
power of sin(θ).  Rather than expressing these derivatives in terms of the Legendre
polynomials — which would subject us to another round of ambiguity — the functions in this
module use automatic differentiation to compute the derivatives explicitly.

Condon and Shortley are not very explicit about the meaning of the spherical coordinates,
but they do describe them as "spherical polar coordinates ``r, \theta, \varphi``".
Immediately before equation (1) of section 4³ (page 50), they define the angular-momentum
operator
```math
L_z = -i \hbar \frac{\partial}{\partial \varphi},
```
which agrees with [our expression](@ref "``L`` operators in spherical coordinates").  This
is followed by equation (8):
```math
\begin{aligned}
L_x + i L_y &= \hbar e^{i\varphi} \left(
  \frac{\partial}{\partial \theta}
  + i \cot\theta \frac{\partial}{\partial \varphi}
\right) \\
L_x - i L_y &= \hbar e^{-i\varphi} \left(
  -\frac{\partial}{\partial \theta}
  + i \cot\theta \frac{\partial}{\partial \varphi}
\right),
\end{aligned}
```
which also agrees with [our results.](@ref "``L_{\pm}`` operators in spherical coordinates")
We can infer that the definitions of the spherical coordinates are consistent with ours.

Condon and Shortley do not give an expression for the Wigner D-matrices.


## Implementing formulas

We begin by writing code that implements the formulas from Condon-Shortley.  We encapsulate
the formulas in a module so that we can test them against the `SphericalFunctions` package.
"""

using TestItems: @testitem  #hide
@testitem "Condon-Shortley conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide

module CondonShortley
#+

# We'll also use some predefined utilities to make the code look more like the equations.
import ..ConventionsUtilities: 𝒾, ❗, dʲsin²ᵏθdcosθʲ
#+

# Equation (12) of section 4³ (page 51) writes the solution to the three-dimensional Laplace
# equation in spherical coordinates as
# ```math
# \psi(\gamma, \ell, m_\ell)
# =
# B(\gamma, \ell) \Theta(\ell, m_\ell) \Phi(m_\ell),
# ```
# where ``B`` is independent of ``\theta`` and ``\varphi``, and ``\gamma`` represents any
# number of eigenvalues required to specify the state.  More explicitly, below Eq. (5) of
# section 5⁵ (page 127), they specifically define the spherical harmonics as
# ```math
# \phi(\ell, m_\ell) = \Theta(\ell, m_\ell) \Phi(m_\ell).
# ```
# One quirk of their notation is that the dependence on ``\theta`` and ``\varphi`` is
# implicit in their functions; we make it explicit, as Julia requires:
function 𝜙(ℓ, mₗ, 𝜃, φ)
    Θ(ℓ, mₗ, 𝜃) * Φ(mₗ, φ)
end
#+

# The ``\varphi`` part is given by equation (5) of section 4³ (page 50):
# ```math
# \Phi(m_\ell)
# =
# \frac{1}{\sqrt{2\pi}} e^{i m_\ell \varphi}.
# ```
# Again, we make the dependence on ``\varphi`` explicit, and we capture its type to ensure
# that we don't lose precision when converting π to a floating-point number.
function Φ(mₗ, φ::T) where {T}
    1 / √(2T(π)) * exp(𝒾 * mₗ * φ)
end
#+

# Equation (15) of section 4³ (page 52) gives the ``\theta`` dependence as
# ```math
# \Theta(\ell, m)
# =
# (-1)^\ell
# \sqrt{\frac{(2\ell+1)}{2} \frac{(\ell+m)!}{(\ell-m)!}}
# \frac{1}{2^\ell \ell!}
# \frac{1}{\sin^m \theta}
# \frac{d^{\ell-m}}{d(\cos\theta)^{\ell-m}} \sin^{2\ell}\theta.
# ```
# Again, we make the dependence on ``\theta`` explicit, and we capture its type to ensure
# that we don't lose precision when converting the factorials to a floating-point number.
function Θ(ℓ, m, 𝜃::T) where {T}
    (-1)^ℓ * T(√(((2ℓ+1) * (ℓ+m)❗) / (2 * (ℓ - m)❗)) * (1 / (2^ℓ * (ℓ)❗))) *
    (1 / sin(𝜃)^T(m)) * dʲsin²ᵏθdcosθʲ(j=ℓ-m, k=ℓ, θ=𝜃)
end
#+

# It may be helpful to check some values against explicit formulas for the first few
# spherical harmonics as given by Condon-Shortley in the footnote to Eq. (15) of Sec. 4³
# (page 52).  Note the subtle difference between the character `Θ` defining the function above
# and the character `ϴ` defining the function below.
ϴ(ℓ, m, 𝜃) = ϴ(Val(ℓ), Val(m), 𝜃)
ϴ(::Val{0}, ::Val{0}, 𝜃) = √(1/2)
ϴ(::Val{1}, ::Val{0}, 𝜃) = √(3/2) * cos(𝜃)
ϴ(::Val{2}, ::Val{0}, 𝜃) = √(5/8) * (2cos(𝜃)^2 - sin(𝜃)^2)
ϴ(::Val{3}, ::Val{0}, 𝜃) = √(7/8) * (2cos(𝜃)^3 - 3cos(𝜃)sin(𝜃)^2)
ϴ(::Val{1}, ::Val{+1}, 𝜃) = -√(3/4) * sin(𝜃)
ϴ(::Val{1}, ::Val{-1}, 𝜃) = +√(3/4) * sin(𝜃)
ϴ(::Val{2}, ::Val{+1}, 𝜃) = -√(15/4) * cos(𝜃) * sin(𝜃)
ϴ(::Val{2}, ::Val{-1}, 𝜃) = +√(15/4) * cos(𝜃) * sin(𝜃)
ϴ(::Val{3}, ::Val{+1}, 𝜃) = -√(21/32) * (4cos(𝜃)^2*sin(𝜃) - sin(𝜃)^3)
ϴ(::Val{3}, ::Val{-1}, 𝜃) = +√(21/32) * (4cos(𝜃)^2*sin(𝜃) - sin(𝜃)^3)
ϴ(::Val{2}, ::Val{+2}, 𝜃) = √(15/16) * sin(𝜃)^2
ϴ(::Val{2}, ::Val{-2}, 𝜃) = √(15/16) * sin(𝜃)^2
ϴ(::Val{3}, ::Val{+2}, 𝜃) = √(105/16) * cos(𝜃) * sin(𝜃)^2
ϴ(::Val{3}, ::Val{-2}, 𝜃) = √(105/16) * cos(𝜃) * sin(𝜃)^2
ϴ(::Val{3}, ::Val{+3}, 𝜃) = -√(35/32) * sin(𝜃)^3
ϴ(::Val{3}, ::Val{-3}, 𝜃) = +√(35/32) * sin(𝜃)^3
#+

# Condon and Shortley do not give an expression for the Wigner D-matrices, but the
# convention for spherical harmonics is what they are known for, so this will suffice.

end  # module CondonShortley
#+

# ## Tests
#
# We can now test the functions against the equivalent functions from the
# `SphericalFunctions` package.  We will need to test approximate floating-point equality,
# so we set absolute and relative tolerances (respectively) in terms of the machine epsilon:
ϵₐ = 100eps()
ϵᵣ = 1000eps()
#+

# The explicit formulas will be a good preliminary test.  In this case, the formulas are
# only given up to
ℓₘₐₓ = 3
#+
# so we test up to that point, and just compare the general form to the explicit formulas —
# again, noting the subtle difference between the characters `Θ` and `ϴ`.  Note that the
# ``1/\sin\theta`` factor in the general form will cause problems at the poles, so we avoid
# the poles by using `βrange` with a small offset:
for θ ∈ θrange(; avoid_poles=ϵₐ/10)
    for (ℓ, m) ∈ ℓmrange(ℓₘₐₓ)
        @test CondonShortley.ϴ(ℓ, m, θ) ≈ CondonShortley.Θ(ℓ, m, θ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Finally, we can test Condon-Shortley's full expressions for spherical harmonics against
# the `SphericalFunctions` package.  We will only test up to
ℓₘₐₓ = 4
#+
# because the formulas are very slow, and this will be sufficient to sort out any sign or
# normalization differences, which are the most likely source of error.
for (θ, ϕ) ∈ θϕrange(; avoid_poles=ϵₐ/40)
    for (ℓ, m) ∈ ℓmrange(ℓₘₐₓ)
        @test CondonShortley.𝜙(ℓ, m, θ, ϕ) ≈ SphericalFunctions.Y(ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# This successful test shows that the function ``\phi`` defined by Condon and Shortley
# agrees with the spherical harmonics defined by the `SphericalFunctions` package.

end  #hide
