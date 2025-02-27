md"""
# Condon-Shortley (1935)

[Condon and Shortley's "The Theory Of Atomic Spectra"](@cite CondonShortley_1935) is the
standard reference for the "Condon-Shortley phase convention".  Though some references are
not very clear about precisely what they mean by that phrase, it seems clear that the
original meaning included the idea that the angular-momentum raising and lowering operators
have eigenvalues that are *real and positive* when acting on the spherical harmonics.  To
avoid ambiguity, we can just look at the actual spherical harmonics they define.

The method we use here is as direct and explicit as possible.  In particular, Condon and
Shortley provide a formula for the φ=0 part in terms of iterated derivatives of a power of
sin(θ).  Rather than expressing these derivatives in terms of the Legendre polynomials —
which would subject us to another round of ambiguity — the functions in this module use
automatic differentiation to compute the derivatives explicitly.

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

The result is that the original Condon-Shortley spherical harmonics agree perfectly with the
ones computed by this package.

(Condon and Shortley do not give an expression for the Wigner D-matrices.)

"""

using TestItems: @testmodule, @testitem  #hide

# ## Function definitions
#
# We begin with some basic code

@testmodule CondonShortley1935 begin  #hide

import FastDifferentiation
const 𝒾 = im
struct Factorial end
Base.:*(n::Integer, ::Factorial) = factorial(big(n))
const ❗ = Factorial()

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
function ϕ(ℓ, mₗ, θ, φ)
    Θ(ℓ, mₗ, θ) * Φ(mₗ, φ)
end

#+
# The ``\varphi`` part is given by equation (5) of section 4³ (page 50):
# ```julia
# 1 / √(2T(π)) * exp(𝒾 * mₗ * φ)
# ```
# ```math
# \Phi(m_\ell)
# =
# \frac{1}{\sqrt{2\pi}} e^{i m_\ell \varphi}.
# ```
# The dependence on ``\varphi`` is implicit, but we make it explicit here:
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
# Again, the dependence on ``\theta`` is implicit, but we make it explicit here:
function Θ(ℓ, m, θ::T) where {T}
    (-1)^ℓ * T(√(((2ℓ+1) * (ℓ+m)❗) / (2 * (ℓ - m)❗)) * (1 / (2^ℓ * (ℓ)❗))) *
    (1 / sin(θ)^T(m)) * dʲsin²ᵏθdcosθʲ(ℓ-m, ℓ, θ)
end

#+
# We can use `FastDifferentiation` to compute the derivative term:
function dʲsin²ᵏθdcosθʲ(j, k, θ)
    if j < 0
        throw(ArgumentError("j=$j must be non-negative"))
    end
    if j == 0
        return sin(θ)^(2k)
    end
    x = FastDifferentiation.make_variables(:x)[1]
    ∂ₓʲfᵏ = FastDifferentiation.derivative((1 - x^2)^k, (x for _ ∈ 1:j)...)
    return FastDifferentiation.make_function([∂ₓʲfᵏ,], [x,])(cos(θ))[1]
end

#+

# It may be helpful to check some values against explicit formulas for the first few
# spherical harmonics as given by Condon-Shortley in the footnote to Eq. (15) of Sec. 4³
# (page 52):
ϴ(ℓ, m, θ) = ϴ(Val(ℓ), Val(m), θ)
ϴ(::Val{0}, ::Val{0}, θ) = √(1/2)
ϴ(::Val{1}, ::Val{0}, θ) = √(3/2) * cos(θ)
ϴ(::Val{2}, ::Val{0}, θ) = √(5/8) * (2cos(θ)^2 - sin(θ)^2)
ϴ(::Val{3}, ::Val{0}, θ) = √(7/8) * (2cos(θ)^3 - 3cos(θ)sin(θ)^2)
ϴ(::Val{1}, ::Val{+1}, θ) = -√(3/4) * sin(θ)
ϴ(::Val{1}, ::Val{-1}, θ) = +√(3/4) * sin(θ)
ϴ(::Val{2}, ::Val{+1}, θ) = -√(15/4) * cos(θ) * sin(θ)
ϴ(::Val{2}, ::Val{-1}, θ) = +√(15/4) * cos(θ) * sin(θ)
ϴ(::Val{3}, ::Val{+1}, θ) = -√(21/32) * (4cos(θ)^2*sin(θ) - sin(θ)^3)
ϴ(::Val{3}, ::Val{-1}, θ) = +√(21/32) * (4cos(θ)^2*sin(θ) - sin(θ)^3)
ϴ(::Val{2}, ::Val{+2}, θ) = √(15/16) * sin(θ)^2
ϴ(::Val{2}, ::Val{-2}, θ) = √(15/16) * sin(θ)^2
ϴ(::Val{3}, ::Val{+2}, θ) = √(105/16) * cos(θ) * sin(θ)^2
ϴ(::Val{3}, ::Val{-2}, θ) = √(105/16) * cos(θ) * sin(θ)^2
ϴ(::Val{3}, ::Val{+3}, θ) = -√(35/32) * sin(θ)^3
ϴ(::Val{3}, ::Val{-3}, θ) = +√(35/32) * sin(θ)^3

#+
# Condon and Shortley do not give an expression for the Wigner D-matrices, but the
# convention for spherical harmonics is what they are known for, so this will suffice.

end  #hide


# ## Tests

@testitem "Condon-Shortley conventions" setup=[Utilities, CondonShortley] begin  #hide

using Random
using Quaternionic: from_spherical_coordinates
#const check = NaNChecker.NaNCheck

Random.seed!(1234)
const T = Float64
const ℓₘₐₓ = 4
ϵₐ = 4eps(T)
ϵᵣ = 1000eps(T)

## Tests for Y(ℓ, m, θ, ϕ)
let Y=CondonShortley.ϕ, Θ=CondonShortley.Θ, ϴ=CondonShortley.ϴ, ϕ=zero(T)
    for θ ∈ βrange(T)
        if abs(sin(θ)) < ϵₐ
            continue
        end

        ## # Find where NaNs are coming from
        ## for ℓ ∈ 0:ℓₘₐₓ
        ##     for m ∈ -ℓ:ℓ
        ##         Θ(ℓ,  m, check(θ))
        ##     end
        ## end

        ## Test footnote to Eq. (15) of Sec. 4³ of Condon-Shortley
        let Y = ₛ𝐘(0, 3, T, [from_spherical_coordinates(θ, ϕ)])[1,:]
            for ℓ ∈ 0:3
                for m ∈ -ℓ:ℓ
                    @test ϴ(ℓ, m, θ) / √(2π) ≈ Y[Yindex(ℓ, m)] atol=ϵₐ rtol=ϵᵣ
                end
            end
        end

        ## Test Eq. (18) of Sec. 4³ of Condon-Shortley
        for ℓ ∈ 0:ℓₘₐₓ
            for m ∈ -ℓ:ℓ
                @test Θ(ℓ, m, θ) ≈ (-1)^(m) * Θ(ℓ, -m, θ) atol=ϵₐ rtol=ϵᵣ
            end
        end

        ## Compare to SphericalHarmonics Y
        let s = 0
            Y₁ = ₛ𝐘(s, ℓₘₐₓ, T, [from_spherical_coordinates(θ, ϕ)])[1,:]
            Y₂ = [Y(ℓ, m, θ, ϕ) for ℓ ∈ abs(s):ℓₘₐₓ for m ∈ -ℓ:ℓ]
            @test Y₁ ≈ Y₂ atol=ϵₐ rtol=ϵᵣ
        end
    end
end

end  #hide
