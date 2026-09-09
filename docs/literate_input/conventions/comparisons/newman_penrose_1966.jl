md"""
# Newman-Penrose (1966)

!!! info "Summary"
    Newman and Penrose's definition of spin weight, and their spin-raising and -lowering
    operators ``\eth`` and ``\bar{\eth}``, agree with the definitions used in the
    `SphericalFunctions` package: ``\eth = R_+`` and ``\bar{\eth} = -R_-``, with
    ``\eth\, {}_sY_{ℓ,m} = \sqrt{(ℓ-s)(ℓ+s+1)}\, {}_{s+1}Y_{ℓ,m}`` and
    ``\bar{\eth}\, {}_sY_{ℓ,m} = -\sqrt{(ℓ+s)(ℓ-s+1)}\, {}_{s-1}Y_{ℓ,m}``.  Their
    (unnormalized) expression for the spin-weighted spherical harmonics in terms of the
    stereographic coordinate ``ζ`` is proportional to ours, with the constant of
    proportionality ``(-1)^ℓ e^{-isϕ} \sqrt{4π / [(2ℓ+1)(ℓ+m)!(ℓ-m)!]}``; the phase
    ``e^{-isϕ}`` reflects the fact that their expression is referred to the stereographic
    tangent basis rather than the ``(\boldsymbol{θ}, \boldsymbol{ϕ})`` basis.

In their 1966 paper, [Newman_1966](@citet), Newman and Penrose first introduced the
spin-weighted spherical harmonics, ``{}_sY_{ℓ m}``.  They use the standard (physicists')
convention for spherical coordinates and introduce the stereographic coordinate
```math
ζ = e^{iϕ} \cot\frac{θ}{2}.
```

## Spin weight

They are a little ambiguous about the relationship of the complex basis vector ``m^\mu``
to the coordinates:

> The vectors ``\Re(m^\mu)`` and ``\Im(m^\mu)`` may be regarded as orthogonal tangent
> vectors (of length ``2^{-1/2}``) at each point of the surface. [...] If spherical polar
> coordinates are used, a natural choice for ``m^\mu`` is to make ``\Re(m^\mu)`` and
> ``\Im(m^\mu)`` tangential, respectively, to the curves ``ϕ = \mathrm{const}`` and ``θ =
> \mathrm{const}.``

The ambiguity is in the sign implied by "tangential", but the natural choice is to assume
they mean that the components are *positive* multiples of ``\partial_θ`` and
``\partial_ϕ`` respectively, in which case we have
```math
m^\mu = \frac{1}{\sqrt{2}} \left[ \partial_θ + i \csc θ \partial_ϕ \right],
```
which is precisely [our ``𝐦``](@ref summary_euler_angles).  They define the spin weight in
terms of behavior of a quantity under rotation of ``m^\mu`` in its own plane as
```math
(m^\mu)' = e^{iψ} m^\mu
=
\frac{1}{\sqrt{2}}
\left[
  \left(\cos ψ\partial_θ - \sin ψ\csc θ \partial_ϕ\right)
  + i \left(\cos ψ\csc θ \partial_ϕ + \sin ψ\partial_θ\right)
\right],
```
and a quantity ``\eta`` has spin weight ``s`` if it transforms as
```math
\eta' = e^{i s ψ} \eta.
```
Raising the spherical coordinates ``(θ, ϕ)`` to Euler angles ``(ϕ, θ, -ψ)``, we see that
the rotor ``𝐑_{ϕ, θ, -ψ}`` rotates the ``𝐳`` basis vector to the point ``(θ, ϕ)``, and it
rotates ``(𝐱 + i 𝐲) / \sqrt{2}`` onto ``(m^\mu)'``.  Supposing that these quantities are
functions of Euler angles, we can write
```math
\eta(ϕ, θ, -ψ) = e^{i s ψ} \eta(ϕ, θ, 0),
\qquad \text{or} \qquad
\eta(ϕ, θ, γ) = e^{-i s γ} \eta(ϕ, θ, 0).
```
Thus, the operator with eigenvalue ``s`` is ``i \partial_γ``, which is exactly [our
``R_z``](@ref summary_spin_weight).  This is the definition of spin weight adopted by this
package, so there is nothing to test numerically here; the rest of this page tests the
operators and functions built on it.

## The operators ``\eth`` and ``\bar{\eth}``

Newman and Penrose define the spin-raising operator ``\eth`` acting on a function of spin
weight ``s`` in their Eq. (3.8) as
```math
\eth \eta
=
-\left(\sin θ\right)^s
\left\{
    \frac{\partial}{\partial θ}
    + \frac{i}{\sin θ} \frac{\partial}{\partial ϕ}
\right\} \left\{\left(\sin θ\right)^{-s} \eta\right\},
```
along with its complex conjugate,
```math
\bar{\eth} \eta
=
-\left(\sin θ\right)^{-s}
\left\{
    \frac{\partial}{\partial θ}
    - \frac{i}{\sin θ} \frac{\partial}{\partial ϕ}
\right\} \left\{\left(\sin θ\right)^{s} \eta\right\},
```
which lowers the spin weight.  On the [Summary](@ref summary_spin_weight) page — and in
detail on the [calculation page](@ref euler_R_S2) — we show that these are precisely
``\eth = R_+ = R_x + i R_y`` and ``\bar{\eth} = -R_- = -(R_x - i R_y)`` in terms of our
right angular-momentum operators, so that
```math
\eth\, {}_sY_{ℓ,m} = \sqrt{(ℓ-s)(ℓ+s+1)}\, {}_{s+1}Y_{ℓ,m},
\qquad
\bar{\eth}\, {}_sY_{ℓ,m} = -\sqrt{(ℓ+s)(ℓ-s+1)}\, {}_{s-1}Y_{ℓ,m}.
```
The minus sign in the second relation is real, and is the reason ``\bar{\eth} \neq R_-``.
Below, we implement Newman and Penrose's differential operators literally (using automatic
differentiation) and verify these relations on the spin-weighted spherical harmonics of this
package.

## Spin-weighted spherical harmonics

Newman and Penrose then compute the spin-weighted spherical harmonics, up to normalization,
as
```math
{}_sY_{ℓ, m}
\propto
\frac{1}{\left[(ℓ-s)! (ℓ+s)!\right]^{1/2}}
\left(1 + ζ \bar{ζ}\right)^{-ℓ}
\sum_p ζ^p (-\bar{ζ})^{p+s-m}
\binom{ℓ-s}{p} \binom{ℓ+s}{p+s-m},
```
where the sum is over all integers ``p`` such that the binomial coefficients are nonzero.
Note that, because ``ζ`` carries a factor ``e^{iϕ}``, this expression is proportional to
``e^{i(m-s)ϕ}`` rather than ``e^{imϕ}``: it is referred to the tangent basis defined by the
stereographic coordinate, which is rotated relative to the ``(\boldsymbol{θ},
\boldsymbol{ϕ})`` basis by the angle ``ϕ``.  Since Newman and Penrose only give this
expression up to proportionality, the test below establishes the exact relationship to our
``{}_sY_{ℓ,m}``.

## Implementing formulas

We begin by writing code that implements the formulas from Newman and Penrose.  We
encapsulate the formulas in a module so that we can test them against the
`SphericalFunctions` package.
"""

# TODO: Confirm Newman-Penrose equation numbers for the spin-weight definition and the ζ-sum.  #src

using TestItems: @testitem  #hide
@testitem "Newman-Penrose conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide

module NewmanPenrose
#+

# We'll use some predefined utilities to make the code look more like the equations, and
# `ForwardDiff` to evaluate the derivatives in the differential operators exactly.
import ..ConventionsUtilities: 𝒾, ❗
import ForwardDiff
#+

# The stereographic coordinate is
ζ(θ, ϕ) = exp(𝒾 * ϕ) * cot(θ / 2)
#+

# The spin-raising operator ``\eth`` of Eq. (3.8) acts on a function `η(θ, ϕ)` of spin
# weight `s`, returning a new function of ``(θ, ϕ)``.  We evaluate the partial derivatives
# with forward-mode automatic differentiation.
function ð(η, s)
    function ðη(θ, ϕ)
        η̃(θ, ϕ) = sin(θ)^(-s) * η(θ, ϕ)
        ∂θ = ForwardDiff.derivative(θ′ -> η̃(θ′, ϕ), θ)
        ∂ϕ = ForwardDiff.derivative(ϕ′ -> η̃(θ, ϕ′), ϕ)
        -sin(θ)^s * (∂θ + 𝒾 / sin(θ) * ∂ϕ)
    end
end
#+

# The spin-lowering operator ``\bar{\eth}`` is the complex conjugate of ``\eth``, meaning
# that it is obtained by conjugating every explicit ``i`` — including the one hidden in
# ``s``, whose sign flips under conjugation because the conjugate of a spin-weight-``s``
# quantity has spin weight ``-s``:
function ð̄(η, s)
    function ð̄η(θ, ϕ)
        η̃(θ, ϕ) = sin(θ)^(s) * η(θ, ϕ)
        ∂θ = ForwardDiff.derivative(θ′ -> η̃(θ′, ϕ), θ)
        ∂ϕ = ForwardDiff.derivative(ϕ′ -> η̃(θ, ϕ′), ϕ)
        -sin(θ)^(-s) * (∂θ - 𝒾 / sin(θ) * ∂ϕ)
    end
end
#+

# Newman and Penrose's unnormalized expression for the spin-weighted spherical harmonics in
# terms of ``ζ`` is a sum over all ``p`` for which both binomial coefficients are nonzero,
# which means ``\max(0, m-s) \leq p \leq \min(ℓ-s, ℓ+m)``:
function ₛYₗₘ_ζ(s, ℓ, m, θ, ϕ)
    z = ζ(θ, ϕ)
    (1 / √((ℓ-s)❗ * (ℓ+s)❗)) * (1 + z * conj(z))^(-ℓ) *
    sum(
        z^p * (-conj(z))^(p+s-m) * binomial(big(ℓ-s), p) * binomial(big(ℓ+s), p+s-m)
        for p ∈ max(0, m-s):min(ℓ-s, ℓ+m);
        init=zero(Complex{typeof(θ)})
    )
end
#+

# Finally, to apply the differential operators we need a version of *our* spin-weighted
# spherical harmonics that can be differentiated.  The functions provided by
# `SphericalFunctions` are evaluated with recursions in floating-point arithmetic, which is
# not suitable for automatic differentiation, so we transcribe the explicit formula from the
# [Summary](@ref summary_swsh) page:
# ```math
# {}_{s}Y_{ℓ,m}(θ, ϕ)
# =
# (-1)^s\sqrt{\frac{2ℓ+1}{4π}}\, e^{imϕ}
# \sum_{k = k_1}^{k_2}
# \frac{(-1)^k[(ℓ+m)!(ℓ-m)!(ℓ-s)!(ℓ+s)!]^{1/2}}
# {(ℓ+m-k)!\,(ℓ+s-k)!\,k!\,(k-s-m)!}
# \left(\cos\frac{θ}{2}\right)^{2ℓ+m+s-2k}
# \left(\sin\frac{θ}{2}\right)^{2k-s-m},
# ```
# with ``k_1 = \max(0, m+s)`` and ``k_2 = \min(ℓ+m, ℓ+s)``.  We check below that this
# agrees with the package, before using it in the differential tests.
#
# We write the phase ``e^{imϕ}`` as `cis(m * ϕ)` rather than `exp(𝒾 * m * ϕ)` because
# `Base.exp(::Complex)` short-circuits to `Complex(exp(real(z)), imag(z))` whenever
# `iszero(imag(z))`.  That branch is correct for numbers, but it is only accurate to *first*
# order in the imaginary part, so under automatic differentiation at ``ϕ = 0`` it silently
# discards every derivative beyond the first — which would break the second-order
# ``\bar{\eth}\eth`` test below.  `cis` computes `Complex(cos, sin)` unconditionally, and
# so differentiates correctly to all orders.
function ₛYₗₘ(s, ℓ, m, θ, ϕ)
    T = float(typeof(θ))
    (-1)^s * √((2ℓ+1) / (4T(π))) * cis(m * ϕ) *
    sum(
        (-1)^k * T(√((ℓ+m)❗ * (ℓ-m)❗ * (ℓ-s)❗ * (ℓ+s)❗) / ((ℓ+m-k)❗ * (ℓ+s-k)❗ * (k)❗ * (k-s-m)❗)) *
        cos(θ/2)^(2ℓ+m+s-2k) * sin(θ/2)^(2k-s-m)
        for k ∈ max(0, m+s):min(ℓ+m, ℓ+s);
        init=zero(Complex{T})
    )
end
#+

end  # module NewmanPenrose
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
# normalization differences, which are the most likely source of error.  For the same
# reason we use a modest grid of points:
θϕs = θϕrange(Float64, 7)
#+
# The differential operators involve ``1/\sin θ``, so for them we avoid the poles.  Even so,
# the factor ``\sin^{-s} θ`` inside ``\eth`` amplifies rounding errors near the poles by a
# factor of order ``1/\sin^{|s|} θ``, so we evaluate the differential tests in `BigFloat`
# arithmetic (the automatic differentiation is exact, so this leaves only the analytic
# comparison), while still requiring agreement to `Float64` precision:
θϕs_big = [(big(θ), big(ϕ)) for (θ, ϕ) ∈ θϕrange(Float64, 5; avoid_poles=1e-3)]
#+

# We will also need the imaginary unit from the utilities module.
import .ConventionsUtilities: 𝒾
#+

# First, we check that our transcription of the explicit formula agrees with the package:
for (θ, ϕ) ∈ θϕs
    for (s, ℓ, m) ∈ sℓmrange(ℓₘₐₓ, sₘₐₓ)
        @test NewmanPenrose.ₛYₗₘ(s, ℓ, m, θ, ϕ) ≈ ConventionsUtilities.Y(s, ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Next, we apply Newman and Penrose's ``\eth`` and ``\bar{\eth}`` to our spin-weighted
# spherical harmonics, and verify the raising and lowering relations — including the
# minus sign in the lowering relation.  When the raised or lowered spin weight would exceed
# ``ℓ`` in magnitude, the result must vanish.
for (θ, ϕ) ∈ θϕs_big
    for (s, ℓ, m) ∈ sℓmrange(ℓₘₐₓ, sₘₐₓ)
        Y(θ, ϕ) = NewmanPenrose.ₛYₗₘ(s, ℓ, m, θ, ϕ)
        ðY = NewmanPenrose.ð(Y, s)(θ, ϕ)
        ð̄Y = NewmanPenrose.ð̄(Y, s)(θ, ϕ)
        if abs(s+1) ≤ ℓ
            @test ðY ≈ √((ℓ-s) * (ℓ+s+1)) * NewmanPenrose.ₛYₗₘ(s+1, ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
        else
            @test ðY ≈ 0 atol=ϵₐ
        end
        if abs(s-1) ≤ ℓ
            @test ð̄Y ≈ -√((ℓ+s) * (ℓ-s+1)) * NewmanPenrose.ₛYₗₘ(s-1, ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
        else
            @test ð̄Y ≈ 0 atol=ϵₐ
        end
    end
end
#+

# A consequence of these two relations is that ``\bar{\eth}\eth`` acts on ``{}_sY_{ℓ,m}``
# as multiplication by ``-(ℓ-s)(ℓ+s+1)``, which is the eigenvalue relation Newman and Penrose
# use to characterize the spin-weighted spherical harmonics.  We check it directly:
for (θ, ϕ) ∈ θϕs_big
    for (s, ℓ, m) ∈ sℓmrange(ℓₘₐₓ, sₘₐₓ)
        Y(θ, ϕ) = NewmanPenrose.ₛYₗₘ(s, ℓ, m, θ, ϕ)
        ð̄ðY = NewmanPenrose.ð̄(NewmanPenrose.ð(Y, s), s+1)(θ, ϕ)
        @test ð̄ðY ≈ -(ℓ-s) * (ℓ+s+1) * Y(θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Finally, we compare Newman and Penrose's ``ζ``-expression to the package.  Since they
# only give it up to proportionality, we determine the constant of proportionality — which
# turns out to be ``(-1)^ℓ e^{-isϕ} \sqrt{4π / [(2ℓ+1)(ℓ+m)!(ℓ-m)!]}`` — and assert it
# exactly.  The ``e^{-isϕ}`` is the rotation from the stereographic tangent basis to the
# ``(\boldsymbol{θ}, \boldsymbol{ϕ})`` basis, as discussed above.  Because ``ζ`` itself
# diverges at the north pole, we again use the pole-avoiding `BigFloat` grid.
for (θ, ϕ) ∈ θϕs_big
    for (s, ℓ, m) ∈ sℓmrange(ℓₘₐₓ, sₘₐₓ)
        c = (-1)^ℓ * exp(-𝒾 * s * ϕ) * √(4big(π) / ((2ℓ+1) * factorial(ℓ+m) * factorial(ℓ-m)))
        @test NewmanPenrose.ₛYₗₘ_ζ(s, ℓ, m, θ, ϕ) ≈ c * ConventionsUtilities.Y(s, ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# These successful tests show that Newman and Penrose's ``\eth`` and ``\bar{\eth}`` act on
# this package's spin-weighted spherical harmonics exactly as the conventions pages claim,
# and that their stereographic expression for ``{}_sY_{ℓ,m}`` is proportional to ours with
# the constant given above.

end  #hide
