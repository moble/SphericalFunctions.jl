md"""
# Cohen-Tannoudji (1991)

!!! info "Summary"
    Cohen-Tannoudji's angular-momentum operators and definition of the spherical harmonics
    agree with the definition used in the `SphericalFunctions` package.

[CohenTannoudji_1991](@citet), by a Nobel-prize winner and collaborators, is an extensive
two-volume set on quantum mechanics that is widely used in graduate courses.

They define spherical coordinates in the usual (physicist's) way in Chapter VI.  They then
compute the angular-momentum operators as [Eqs.  (D-5)]
```math
\begin{aligned}
L_x &= i \hbar \left(
    \sin ϕ \frac{\partial} {\partial θ}
    + \frac{\cos ϕ}{\tan θ} \frac{\partial} {\partial ϕ}
\right),
\\
L_y &= i \hbar \left(
    -\cos ϕ \frac{\partial} {\partial θ}
    + \frac{\sin ϕ}{\tan θ} \frac{\partial} {\partial ϕ}
\right),
\\
L_z &= \frac{\hbar}{i} \frac{\partial} {\partial ϕ},
\end{aligned}
```
which agree precisely with [the results in our conventions](@ref
L-operators-in-spherical-coordinates).

In Complement ``\mathrm{B}_{\mathrm{VI}}`` they define a rotation operator ``R`` as acting
on a state such that [Eq. (21)]
```math
\langle 𝐫 | R | ψ \rangle
=
\langle \mathscr{R}^{-1} 𝐫 | ψ \rangle.
```
For an infinitesimal rotation through angle ``dα`` about the axis ``𝐮``, he
shows [Eq. (49)]
```math
R_{𝐮}(dα) = 1 - \frac{i}{\hbar} dα 𝐋.𝐮.
```


## Implementing formulas

We begin by writing code that implements the formulas from Cohen-Tannoudji.  We encapsulate
the formulas in a module so that we can test them against the `SphericalFunctions` package.
"""

using TestItems: @testitem  #hide
@testitem "Cohen-Tannoudji conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide
import .Utilities: ℓmrange, θϕrange  #hide

module CohenTannoudji
#+

# We'll also use some predefined utilities to make the code look more like the equations,
# and `ForwardDiff` to evaluate the derivatives in the angular-momentum operators.
import ..ConventionsUtilities: 𝒾, ❗, dʲsin²ᵏθdcosθʲ
import ForwardDiff
#+

# Cohen-Tannoudji include ``\hbar``, so we will include it in the expressions, but we will
# set it to 1 to match the conventions of the `SphericalFunctions` package.
const ħ = 1
#+

# They derive the spherical harmonics in two ways and get two different, but equivalent,
# expressions in Complement ``\mathrm{A}_{\mathrm{VI}}``.  The first is Eq. (26)
# ```math
# Y_{l}^{m}(θ, ϕ)
# =
# \frac{(-1)^l}{2^l l!} \sqrt{\frac{(2l+1)}{4π} \frac{(l+m)!}{(l-m)!}}
# e^{i m ϕ} (\sin θ)^{-m}
# \frac{d^{l-m}}{d(\cos θ)^{l-m}} (\sin θ)^{2l},
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
# Y_{l}^{m}(θ, ϕ)
# =
# \frac{(-1)^{l+m}}{2^l l!} \sqrt{\frac{(2l+1)}{4π} \frac{(l-m)!}{(l+m)!}}
# e^{i m ϕ} (\sin θ)^m
# \frac{d^{l+m}}{d(\cos θ)^{l+m}} (\sin θ)^{2l}.
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

# Eqs. (D-5) and (D-6) give the angular-momentum operators.  Each operator takes a function
# `f(θ, ϕ)` and returns a new function, with the derivatives evaluated by forward-mode
# automatic differentiation.
∂θ(f) = (θ, ϕ) -> ForwardDiff.derivative(θ′ -> f(θ′, ϕ), θ)
∂ϕ(f) = (θ, ϕ) -> ForwardDiff.derivative(ϕ′ -> f(θ, ϕ′), ϕ)
L_z(f) = (θ, ϕ) -> ħ / 𝒾 * ∂ϕ(f)(θ, ϕ)
L₊(f) = (θ, ϕ) -> ħ * exp(𝒾 * ϕ) * (∂θ(f)(θ, ϕ) + 𝒾 * cot(θ) * ∂ϕ(f)(θ, ϕ))
L₋(f) = (θ, ϕ) -> ħ * exp(-𝒾 * ϕ) * (-∂θ(f)(θ, ϕ) + 𝒾 * cot(θ) * ∂ϕ(f)(θ, ϕ))
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
ℓₘₐₓ = 6
#+
#
# because the symbolic derivatives in the formulas become expensive to compute at higher
# orders, and this will be sufficient to sort out any sign or normalization differences,
# which are the most likely source of error.  Also, the formulas are singular at the poles,
# so we avoid evaluating there.
for (θ, ϕ) ∈ θϕrange(rng; avoid_poles=ϵₐ/40)
    for (ℓ, Yˡ) ∈ SphericalFunctions.YlmCalculator(θ, ϕ, ℓₘₐₓ), m ∈ -ℓ:ℓ
        @test CohenTannoudji.Y₁(ℓ, m, θ, ϕ) ≈ Yˡ[m] atol=ϵₐ rtol=ϵᵣ
        @test CohenTannoudji.Y₂(ℓ, m, θ, ϕ) ≈ Yˡ[m] atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Their operators (D-5) agree with ours if they act on the spherical harmonics as [our
# summary page](@ref summary_swsh) says ours do (with ``s = 0``): ``L_z Y_{ℓ,m} = m
# Y_{ℓ,m}`` and ``(L_x \pm i L_y) Y_{ℓ,m} = \sqrt{(ℓ \mp m)(ℓ \pm m + 1)}\, Y_{ℓ,m \pm 1}``,
# where the raised or lowered function vanishes at the edges ``m = ±ℓ``.  We apply them to
# the spherical harmonics of Eq. (26), which we have just shown to be ours.  The operators
# involve ``\cot θ``, so we avoid the poles, and we evaluate in `BigFloat` arithmetic to
# keep the rounding errors amplified near the poles below `Float64` precision.
for (θ, ϕ) ∈ ((big(θ), big(ϕ)) for (θ, ϕ) ∈ θϕrange(rng, Float64, 5; avoid_poles=1e-3))
    for (ℓ, m) ∈ ℓmrange(4)
        Y = (θ, ϕ) -> CohenTannoudji.Y₁(ℓ, m, promote(θ, ϕ)...)
        L₊Y = CohenTannoudji.L₊(Y)(θ, ϕ)
        L₋Y = CohenTannoudji.L₋(Y)(θ, ϕ)
        @test CohenTannoudji.L_z(Y)(θ, ϕ) ≈ m * Y(θ, ϕ) atol=ϵₐ rtol=ϵᵣ
        if m < ℓ
            @test L₊Y ≈ √((ℓ-m) * (ℓ+m+1)) * CohenTannoudji.Y₁(ℓ, m+1, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
        else
            @test L₊Y ≈ 0 atol=ϵₐ
        end
        if m > -ℓ
            @test L₋Y ≈ √((ℓ+m) * (ℓ-m+1)) * CohenTannoudji.Y₁(ℓ, m-1, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
        else
            @test L₋Y ≈ 0 atol=ϵₐ
        end
    end
end
#+

# These successful tests show that both versions of the spherical harmonics given by
# Cohen-Tannoudji agree with the spherical harmonics defined by the `SphericalFunctions`
# package, and that their angular-momentum operators agree with ours.

end  #hide
