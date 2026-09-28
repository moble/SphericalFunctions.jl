md"""
# Edmonds (1960)

!!! info "Summary"
    Edmonds' Euler angles, spherical coordinates, angular-momentum operators, and spherical
    harmonics agree with those used in the `SphericalFunctions` package.  His rotation
    matrices are those of the *inverse* rotation, so that
    ```math
    𝒟^{(j)}_{m'm}(α, β, γ)\big|_{\text{Edmonds}}
    = 𝔇^{(j)}_{m',m}(-γ, -β, -α)
    = \overline{𝔇^{(j)}_{m,m'}(α, β, γ)},
    ```
    and correspondingly his ``d^{(j)}_{m'm}(β)`` is the transpose of ours.

[Edmonds_2016](@citet) is a standard reference for the theory of angular momentum.  (The
equation numbers below refer to the 2016 Princeton reprint of the 1960 second edition, which
is the edition cited in the bibliography.)

## Euler angles and spherical coordinates

In Sec. 1.3 he actually does a fair job of defining the Euler angles.  The upshot is that
his definition agrees with ours, though he uses the "active" definition style.  That is, the
rotations are to be performed successively in order:

> 1. A rotation ``α(0 \leq α < 2π)`` about the ``z``-axis, bringing the frame of axes from
>    the initial position ``S`` into the position ``S'``.  The axis of this rotation is
>    commonly called the *vertical*.
>
> 2. A rotation ``β(0 \leq β < π)`` about the ``y``-axis of the frame ``S'``, called the
>    *line of nodes*.  Note that its position is in general different from the initial
>    position of the ``y``-axis of the frame ``S``. The resulting position of the frame of
>    axes is symbolized by ``S''``.
>
> 3. A rotation ``γ(0 \leq γ < 2π)`` about the ``z``-axis of the frame of axes ``S''``,
>    called the *figure axis*; the position of this axis depends on the previous rotations
>    ``α`` and ``β``.  The final position of the frame is symbolized by ``S'''``.

I would simply write the "``y``-axis of the frame ``S'``" as ``y'``, and so on.  In
quaternionic language, I would write these rotations as ``\exp[γ 𝐤''/2]\, \exp[β 𝐣'/2]\,
\exp[α 𝐤/2]``.  But we also have
```math
\exp[β 𝐣'/2] = \exp[α 𝐤/2]\, \exp[β 𝐣/2]\, \exp[-α 𝐤/2]
```
so we can just swap the ``α`` rotation with the ``β`` rotation while dropping the prime from
``𝐣'``.  We can do a similar trick swapping the ``α`` and ``β`` rotations with the ``γ``
rotation while dropping the double prime from ``𝐤''``.  That is, an easy calculation shows
that
```math
\exp[γ 𝐤''/2]\, \exp[β 𝐣'/2]\, \exp[α 𝐤/2]
=
\exp[α 𝐤/2]\, \exp[β 𝐣/2]\, \exp[γ 𝐤/2],
```
which is precisely [our definition](@ref summary_euler_angles).

The spherical coordinates are implicitly defined by this statement:

> It should be noted that the polar coordinates ``φ, θ`` with respect to the original frame
> ``S`` of the ``z``-axis in its final position are identical with the Euler angles ``α, β``
> respectively.

Again, this agrees with our definition.

## Angular-momentum operators

His expression for the angular-momentum operator in Euler angles — Eq. (2.2.2) — agrees with
[ours](@ref summary_L_R_euler):
```math
\begin{aligned}
L_x &= -i \hbar \left\{
    -\frac{\cos α}{\tan β} \frac{\partial} {\partial α}
    - \sin α \frac{\partial} {\partial β}
    + \frac{\cos α}{\sin β} \frac{\partial} {\partial γ}
\right\},
\\
L_y &= -i \hbar \left\{
    -\frac{\sin α}{\tan β} \frac{\partial} {\partial α}
    + \cos α \frac{\partial} {\partial β}
    +\frac{\sin α}{\sin β} \frac{\partial} {\partial γ}
\right\},
\\
L_z &= -i \hbar \frac{\partial} {\partial α}.
\end{aligned}
```
(The corresponding restriction to spherical coordinates also precisely agrees with our
results, with the extra factor of ``\hbar``.)

## Spherical harmonics

Equation (2.5.5) gives the spherical harmonics as
```math
Y_{ℓm}(θ, φ)
=
\frac{(-1)^{ℓ+m}}{2^ℓ ℓ!}
\sqrt{\frac{(2ℓ+1)\,(ℓ-m)!}{4π\,(ℓ+m)!}}\;
\sin^m θ\,
\frac{d^{ℓ+m}}{d(\cos θ)^{ℓ+m}} \sin^{2ℓ} θ\;
e^{imφ},
```
which is the Condon–Shortley form written with ``ℓ+m`` derivatives rather than ``ℓ-m``, and
we expect it to agree with [ours](@ref summary_spherical_harmonics).

## Rotation matrices

Here is where the disagreement lies.  Edmonds defines the rotation operator in Eq. (4.1.9)
as
```math
𝒟(α β γ)
=
\exp\left( \frac{iγ}{\hbar} J_z\right)
\exp\left( \frac{iβ}{\hbar} J_y\right)
\exp\left( \frac{iα}{\hbar} J_z\right),
```
and its matrix elements in Eq. (4.1.12) as
```math
𝒟^{(j)}_{m'm}(α β γ)
=
e^{im'γ}\, d^{(j)}_{m'm}(β)\, e^{imα},
```
with the explicit formula of Eq. (4.1.15),
```math
d^{(j)}_{m'm}(β)
=
\left[\frac{(j+m')!\,(j-m')!}{(j+m)!\,(j-m)!}\right]^{1/2}
\sum_σ \binom{j+m}{j-m'-σ} \binom{j-m}{σ}
(-1)^{j-m'-σ}
\left(\cos\frac{β}{2}\right)^{2σ+m'+m}
\left(\sin\frac{β}{2}\right)^{2j-2σ-m'-m}.
```
Compared with [our ``U(𝐑_{α,β,γ}) = e^{-iαL_z} e^{-iβL_y} e^{-iγL_z}``](@ref
summary_wigner_D), Edmonds' operator is exactly ``U(𝐑_{α,β,γ}^{-1}) = U(𝐑_{-γ,-β,-α})``:
the signs in the exponents are reversed *and* the order of the angles is reversed.  This is
the natural consequence of his "active" description, in which the Euler angles move the
*frame of axes* to its new position — that is, the rotation of the coordinates rather than
of the field.  Since ``𝔇(𝐑^{-1}) = 𝔇(𝐑)^\dagger``, we expect
```math
𝒟^{(j)}_{m'm}(α, β, γ)\big|_{\text{Edmonds}}
= 𝔇^{(j)}_{m',m}(-γ, -β, -α)
= \overline{𝔇^{(j)}_{m,m'}(α, β, γ)},
\qquad
d^{(j)}_{m'm}(β)\big|_{\text{Edmonds}} = d^{(j)}_{m,m'}(β),
```
which the tests below confirm.  Note that this is *not* simply the complex conjugate of our
``𝔇``: the indices are also transposed.  (Several later sources — e.g., [SymPy](@ref "SymPy
(2026)") — cite Edmonds for their conventions but transcribe him inexactly.)

## Implementing formulas

We begin by writing code that implements the formulas from Edmonds.  We encapsulate the
formulas in a module so that we can test them against the `SphericalFunctions` package.
"""

using TestItems: @testitem  #hide
@testitem "Edmonds conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide
import .Utilities: ℓm′mrange, βrange, αβγrange, θϕrange  #hide

module Edmonds
#+

# We'll use some predefined utilities to make the code look more like the equations, and
# `ForwardDiff` to evaluate the derivatives in the angular-momentum operators.
import ..ConventionsUtilities: 𝒾, ❗, dʲsin²ᵏθdcosθʲ
import ForwardDiff
#+

# Edmonds includes ``\hbar``, so we will include it in the expressions, but we will
# set it to 1 to match the conventions of the `SphericalFunctions` package.
const ħ = 1
#+

# Equation (2.5.5).  We capture the floating-point type `T` to ensure that we don't lose
# precision when converting π and the factorials to floating-point numbers.
function Y(ℓ, m, θ::T, φ::T) where {T<:Real}
    (-1)^(ℓ+m) / T(2^ℓ * (ℓ)❗) * √T((2ℓ+1) * (ℓ-m)❗ / (4big(π) * (ℓ+m)❗)) *
    sin(θ)^m * dʲsin²ᵏθdcosθʲ(j=ℓ+m, k=ℓ, θ=θ) * exp(𝒾 * m * φ)
end
#+

# Equation (4.1.15).  The sum runs over all ``σ`` for which the binomial coefficients are
# nonzero, ``\max(0, -m'-m) \leq σ \leq \min(j-m', j-m)``.
function d(j, m′, m, β::T) where {T<:Real}
    √T((j+m′)❗ * (j-m′)❗ / ((j+m)❗ * (j-m)❗)) *
    sum(
        T(binomial(big(j+m), j-m′-σ) * binomial(big(j-m), σ)) * (-1)^(j-m′-σ)
        * cos(β/2)^(2σ+m′+m) * sin(β/2)^(2j-2σ-m′-m)
        for σ ∈ max(0, -m′-m):min(j-m′, j-m);
        init=zero(T)
    )
end
#+

# Equation (4.1.12):
function 𝒟(j, m′, m, α, β, γ)
    exp(𝒾 * m′ * γ) * d(j, m′, m, β) * exp(𝒾 * m * α)
end
#+

# Equation (2.2.2), the angular-momentum operators in Euler angles.  Note that we defined
# ``\hbar = 1`` above.  Each operator takes a function `f(α, β, γ)` and returns a new
# function, with the derivatives evaluated by forward-mode automatic differentiation.
∂α(f) = (α, β, γ) -> ForwardDiff.derivative(α′ -> f(α′, β, γ), α)
∂β(f) = (α, β, γ) -> ForwardDiff.derivative(β′ -> f(α, β′, γ), β)
∂γ(f) = (α, β, γ) -> ForwardDiff.derivative(γ′ -> f(α, β, γ′), γ)
L_x(f) = (α, β, γ) -> -𝒾 * ħ * (-cos(α)*cot(β) * ∂α(f)(α, β, γ) - sin(α) * ∂β(f)(α, β, γ) + cos(α)/sin(β) * ∂γ(f)(α, β, γ))
L_y(f) = (α, β, γ) -> -𝒾 * ħ * (-sin(α)*cot(β) * ∂α(f)(α, β, γ) + cos(α) * ∂β(f)(α, β, γ) + sin(α)/sin(β) * ∂γ(f)(α, β, γ))
L_z(f) = (α, β, γ) -> -𝒾 * ħ * ∂α(f)(α, β, γ)
#+

end  # module Edmonds
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
# we use modest grids of points.
αβγs = αβγrange(rng, Float64, 5)
#+

# First, the spherical harmonics agree with ours.  The formula involves ``\sin^m θ`` with
# ``m`` possibly negative, so we avoid the poles.  Close to them, the formula divides a nearly
# vanishing derivative by a power of the small quantity ``\sin θ``, which amplifies its
# rounding errors enormously — in `Float64` arithmetic at ``θ ≈ 10^{-8}``, to an absolute
# error of about ``2 × 10^{-8}`` — so we evaluate it in `BigFloat` arithmetic.
for (θ, ϕ) ∈ θϕrange(rng, Float64, 7; avoid_poles=ϵₐ/40)
    for (ℓ, Yˡ) ∈ SphericalFunctions.YlmCalculator(θ, ϕ, ℓₘₐₓ), m ∈ -ℓ:ℓ
        @test Edmonds.Y(ℓ, m, big(θ), big(ϕ)) ≈ Yˡ[m] atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Edmonds' ``d`` is the transpose of ours:
for β ∈ βrange(rng)
    for (j, dʲ) ∈ SphericalFunctions.dCalculator(β, ℓₘₐₓ), m′ ∈ -j:j, m ∈ -j:j
        @test Edmonds.d(j, m′, m, β) ≈ dʲ[m, m′] atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# And Edmonds' ``𝒟`` is our ``𝔇`` of the inverse rotation — equivalently, the Hermitian
# conjugate of ours:
for (α, β, γ) ∈ αβγs
    𝔇 = SphericalFunctions.DCalculator(α, β, γ, ℓₘₐₓ)
    𝔇⁻¹ = SphericalFunctions.DCalculator(-γ, -β, -α, ℓₘₐₓ)  # the inverse rotation
    for ((j, 𝔇ʲ), (_, 𝔇ʲ⁻¹)) ∈ zip(𝔇, 𝔇⁻¹), m′ ∈ -j:j, m ∈ -j:j
        @test Edmonds.𝒟(j, m′, m, α, β, γ) ≈ 𝔇ʲ⁻¹[m′, m] atol=ϵₐ rtol=ϵᵣ
        @test Edmonds.𝒟(j, m′, m, α, β, γ) ≈ conj(𝔇ʲ[m, m′]) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Finally, his operators (2.2.2) agree with ours if they act on our ``𝔇`` as [our summary
# page](@ref summary_wigner_D) says ours do: ``L_z 𝔇^{(j)}_{m',m} = -m' 𝔇^{(j)}_{m',m}`` and
# ``(L_x \pm i L_y) 𝔇^{(j)}_{m',m} = -\sqrt{(j \pm m')(j \mp m' + 1)}\, 𝔇^{(j)}_{m' \mp
# 1,m}``, where the lowered or raised function vanishes at the edges.  We apply them to our
# ``𝔇`` in the form ``\overline{𝒟^{(j)}_{m m'}}`` just established, which automatic
# differentiation can evaluate.  The operators involve ``1/\sin β``, so we use a few
# generic Euler-angle triples away from the poles.
𝔇(j, m′, m, α, β, γ) = conj(Edmonds.𝒟(j, m, m′, α, β, γ))
for (α, β, γ) ∈ [(0.7, 1.1, 2.3), (2.9, 0.4, 5.1), (4.0, 2.2, 0.3), (1.3, 2.9, 4.7)]
    for (j, m′, m) ∈ ℓm′mrange(3)
        f = (α, β, γ) -> 𝔇(j, m′, m, α, β, γ)
        L₊f = Edmonds.L_x(f)(α, β, γ) + im * Edmonds.L_y(f)(α, β, γ)
        L₋f = Edmonds.L_x(f)(α, β, γ) - im * Edmonds.L_y(f)(α, β, γ)
        @test Edmonds.L_z(f)(α, β, γ) ≈ -m′ * f(α, β, γ) atol=ϵₐ rtol=ϵᵣ
        @test L₊f ≈ (m′ > -j ? -√((j+m′) * (j-m′+1)) * 𝔇(j, m′-1, m, α, β, γ) : 0) atol=ϵₐ rtol=ϵᵣ
        @test L₋f ≈ (m′ < j ? -√((j-m′) * (j+m′+1)) * 𝔇(j, m′+1, m, α, β, γ) : 0) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# These successful tests show that Edmonds' spherical harmonics and angular-momentum
# operators agree with ours, while his rotation matrices are those of the inverse rotation:
# the Hermitian conjugates of the ``𝔇`` matrices defined by the `SphericalFunctions`
# package.

end  #hide
