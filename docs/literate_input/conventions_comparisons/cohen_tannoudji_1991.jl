md"""
# Cohen-Tannoudji (1991)

!!! info "Summary"
    Cohen-Tannoudji's definition of the spherical harmonics agrees with the definition used
    in the `SphericalFunctions` package.

    TODO: Compare angular-momentum operators and rotation operator.

[CohenTannoudji_1991](@citet), by a Nobel-prize winner and collaborators, is an extensive
two-volume set on quantum mechanics that is widely used in graduate courses.

They define spherical coordinates in the usual (physicist's) way in Chapter VI.  They then
compute the angular-momentum operators as [Eqs.  (D-5)]
```math
\begin{aligned}
L_x &= i \hbar \left(
    \sin\phi \frac{\partial} {\partial \theta}
    + \frac{\cos\phi}{\tan\theta} \frac{\partial} {\partial \phi}
\right),
\\
L_y &= i \hbar \left(
    -\cos\phi \frac{\partial} {\partial \theta}
    + \frac{\sin\phi}{\tan\theta} \frac{\partial} {\partial \phi}
\right),
\\
L_z &= \frac{\hbar}{i} \frac{\partial} {\partial \phi}.
\end{aligned}
```
In Complement ``\mathrm{B}_{\mathrm{VI}}`` they define a rotation operator ``R`` as acting
on a state such that [Eq. (21)]
```math
\langle \mathbf{r} | R | \psi \rangle
=
\langle \mathscr{R}^{-1} \mathbf{r} | \psi \rangle.
```
For an infinitesimal rotation through angle ``d\alpha`` about the axis ``\mathbf{u}``, he
shows [Eq. (49)]
```math
R_{\mathbf{u}}(d\alpha) = 1 - \frac{i}{\hbar} d\alpha \mathbf{L}.\mathbf{u}.
```


## Implementing formulas

We begin by writing code that implements the formulas from Cohen-Tannoudji.  We encapsulate
the formulas in a module so that we can test them against the `SphericalFunctions` package.
"""

using TestItems: @testitem  #hide
@testitem "Cohen-Tannoudji conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide

module CohenTannoudji
#+

# We'll also use some predefined utilities to make the code look more like the equations.
import ..ConventionsUtilities: 𝒾, ❗, dʲsin²ᵏθdcosθʲ
#+

# They derive the spherical harmonics in two ways and get two different, but equivalent,
# expressions in Complement ``\mathrm{A}_{\mathrm{VI}}``.  The first is Eq. (26)
# ```math
# Y_{l}^{m}(\theta, \phi)
# =
# \frac{(-1)^l}{2^l l!} \sqrt{\frac{(2l+1)}{4\pi} \frac{(l+m)!}{(l-m)!}}
# e^{i m \phi} (\sin \theta)^{-m}
# \frac{d^{l-m}}{d(\cos \theta)^{l-m}} (\sin \theta)^{2l},
# ```
function Y₁(l, m, θ::T, ϕ::T) where {T<:Real}
    (
        (-1)^l / (2^l * (l)❗)
        * √((2l + 1) / (4T(π)) * (l + m)❗ / (l - m)❗)
        * exp(𝒾 * m * ϕ) * sin(θ)^(-m)
        * dʲsin²ᵏθdcosθʲ(j=l-m, k=l, θ=θ)
    )
end
#+

# while the second is Eq. (30)
# ```math
# Y_{l}^{m}(\theta, \phi)
# =
# \frac{(-1)^{l+m}}{2^l l!} \sqrt{\frac{(2l+1)}{4\pi} \frac{(l-m)!}{(l+m)!}}
# e^{i m \phi} (\sin \theta)^m
# \frac{d^{l+m}}{d(\cos \theta)^{l+m}} (\sin \theta)^{2l}.
# ```
function Y₂(l, m, θ::T, ϕ::T) where {T<:Real}
    (
        (-1)^(l+m) / (2^l * (l)❗)
        * √((2l + 1) / (4T(π)) * (l - m)❗ / (l + m)❗)
        * exp(𝒾 * m * ϕ) * sin(θ)^m
        * dʲsin²ᵏθdcosθʲ(j=l+m, k=l, θ=θ)
    )
end
#+

# Cohen-Tannoudji do not give an expression for the Wigner D-matrices, but the comparisons
# of the definitions of the angular-momentum operators and the rotation operator are also
# useful for comparison, and comparing the spherical harmonics is also important.

end  # module CohenTannoudji
#+

# ## Tests
#
# We can now test the functions against the equivalent functions from the
# `SphericalFunctions` package.  We will need to test approximate floating-point equality,
# so we set absolute and relative tolerances (respectively) in terms of the machine epsilon:
ϵₐ = 100eps()
ϵᵣ = 1000eps()
#+

# We will only test up to
ℓₘₐₓ = 2
#+
#
# because the formulas are very slow, and this will be sufficient to sort out any sign or
# normalization differences, which are the most likely source of error.  Also, the formulas
# are singular at the poles, so we avoid evaluating there.
for (θ, ϕ) ∈ θϕrange(; avoid_poles=ϵₐ/40)
    for (ℓ, m) ∈ ℓmrange(ℓₘₐₓ)
        @test CohenTannoudji.Y₁(ℓ, m, θ, ϕ) ≈ SphericalFunctions.Y(ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
        @test CohenTannoudji.Y₂(ℓ, m, θ, ϕ) ≈ SphericalFunctions.Y(ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# This successful test shows that both versions of the spherical harmonics given by
# Cohen-Tannoudji agree with the spherical harmonics defined by the `SphericalFunctions`
# package.

end  #hide
