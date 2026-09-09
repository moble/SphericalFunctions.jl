md"""
# Wikipedia (2026)

!!! info "Summary"
    Wikipedia's definitions of the Wigner ``D`` and ``d`` matrices, the spherical harmonics,
    and the spin-weighted spherical harmonics all agree with the definitions used in the
    `SphericalFunctions` package.  Wikipedia's "body-fixed" angular-momentum operators
    ``\mathcal{P}`` are the negatives of our ``R`` operators.

Wikipedia is not a primary source, but it is probably the first place most people look, so
it is worth recording exactly how its conventions relate to ours.  The relevant articles
are [Wigner D-matrix](@cite Wikipedia_WignerD), [Spin-weighted spherical
harmonics](@cite Wikipedia_SWSH), [Spherical harmonics](@cite Wikipedia_SphericalHarmonics),
and [Euler angles](@cite Wikipedia_EulerAngles).  All quotations below were taken from the
versions of those pages current on 2026-09-08 (the Wigner D-matrix article was last edited
2026-05-04 and the spin-weighted spherical harmonics article 2026-02-09 at that time).
Wikipedia articles change, so a later reader should expect to re-check them.

## Euler angles

The [Euler angles](@cite Wikipedia_EulerAngles) article describes our convention as the
"proper Euler angles" in the ``z``-``y``-``z`` (or ``Z``-``Y'``-``Z''``) form, and notes that
an *extrinsic* sequence of rotations about fixed axes is equivalent to an *intrinsic*
sequence about the moving axes taken in the reverse order — which is exactly the
[observation](@ref summary_euler_angles) that lets us write
``𝐑_{α,β,γ} = e^{α𝐤/2} e^{β𝐣/2} e^{γ𝐤/2}`` either way.  There is nothing to test here.

## Wigner ``D``-matrix

The [Wigner D-matrix](@cite Wikipedia_WignerD) article defines the rotation operator
```math
\mathcal{R}(α, β, γ) = e^{-iα J_z} e^{-iβ J_y} e^{-iγ J_z},
```
and then
```math
D^j_{m'm}(α, β, γ)
\equiv
\langle j m' | \mathcal{R}(α, β, γ) | j m \rangle
=
e^{-im'α}\, d^j_{m'm}(β)\, e^{-imγ},
```
which is precisely [our definition](@ref summary_wigner_D).  The article gives Wigner's
formula for the ``d``-matrix as
```math
d^j_{m'm}(β)
=
\left[(j+m')!\,(j-m')!\,(j+m)!\,(j-m)!\right]^{1/2}
\sum_{s=s_{\min}}^{s_{\max}}
\left[
  \frac{(-1)^{m'-m+s}
  \left(\cos\frac{β}{2}\right)^{2j+m-m'-2s}
  \left(\sin\frac{β}{2}\right)^{m'-m+2s}}
  {(j+m-s)!\, s!\, (m'-m+s)!\, (j-m'-s)!}
\right],
```
with ``s_{\min} = \max(0, m-m')`` and ``s_{\max} = \min(j+m, j-m')``, and lists the
explicit elements for ``j = 1/2, 1, 2``, of which we transcribe those for ``j = 1`` and
``j = 2``.  It records the symmetries
```math
D^j_{m'm}(α, β, γ) = (-1)^{m'-m}\, \overline{D^j_{-m',-m}(α, β, γ)},
\qquad
d^j_{m',m}(-β) = d^j_{m,m'}(β) = (-1)^{m'-m}\, d^j_{m',m}(β),
```
and the relation to the spherical harmonics
```math
D^{ℓ}_{m0}(α, β, γ) = \sqrt{\frac{4π}{2ℓ+1}}\; \overline{Y_{ℓ}^{m}(β, α)},
```
all of which agree with [ours](@ref summary_wigner_D).

The section "Properties of the Wigner D-matrix" introduces "space-fixed" operators
``\hat{\mathcal{J}}`` and "body-fixed" operators ``\hat{\mathcal{P}}`` in terms of
Euler-angle derivatives, with
```math
\hat{\mathcal{J}}_3 = -i \frac{\partial}{\partial α},
\qquad
\hat{\mathcal{P}}_3 = -i \frac{\partial}{\partial γ},
```
and states that it is the *complex conjugate* ``D^{j\ast}_{m'm}`` that is the
simultaneous eigenfunction, with ``\mathcal{J}_3 D^{j\ast}_{m'm} = m' D^{j\ast}_{m'm}`` and
``\mathcal{P}_3 D^{j\ast}_{m'm} = m D^{j\ast}_{m'm}``.  Comparing with [our operators](@ref
summary_L_R_euler), ``\hat{\mathcal{J}}_3 = L_z`` while ``\hat{\mathcal{P}}_3 = -R_z``; more
generally ``\mathcal{P} = -R`` (see [the calculation page](@ref euler_R_S3)).  Conjugating
Wikipedia's eigenvalue statements gives ``L_z D_{m'm} = -m' D_{m'm}`` and ``R_z D_{m'm} = m
D_{m'm}``, which is exactly what [we find](@ref summary_wigner_D).

## Spherical harmonics

The [Spherical harmonics](@cite Wikipedia_SphericalHarmonics) article uses standard
physicists' coordinates (``θ`` the colatitude, ``ϕ`` the azimuth), and gives the
quantum-mechanical definition, including the Condon–Shortley phase, as
```math
Y_ℓ^m(θ, ϕ)
=
(-1)^m \sqrt{\frac{2ℓ+1}{4π} \frac{(ℓ-m)!}{(ℓ+m)!}}\;
\tilde{P}_ℓ^m(\cos θ)\, e^{imϕ},
```
where ``\tilde{P}_ℓ^m`` is the associated Legendre function *without* the Condon–Shortley
phase,
```math
\tilde{P}_ℓ^m(x) = (1-x^2)^{m/2} \frac{d^m}{dx^m} P_ℓ(x),
\qquad
P_ℓ(x) = \frac{1}{2^ℓ ℓ!} \frac{d^ℓ}{dx^ℓ} (x^2-1)^ℓ,
```
so as to "avoid counting the phase twice"; and
```math
\overline{Y_ℓ^m(θ, ϕ)} = (-1)^m\, Y_ℓ^{-m}(θ, ϕ).
```
We implement this with the derivatives evaluated by automatic differentiation, exactly as
on the [Condon–Shortley page](@ref "Condon-Shortley (1935)").

## Spin-weighted spherical harmonics

The [Spin-weighted spherical harmonics](@cite Wikipedia_SWSH) article follows Newman and
Penrose: a quantity of spin weight ``s`` picks up ``e^{isθ}`` when the tangent basis is
rotated by ``θ``, and the operators
```math
\eth \eta = -(\sin θ)^s \left\{ \frac{\partial}{\partial θ} + \frac{i}{\sin θ}
\frac{\partial}{\partial ϕ} \right\} \left[ (\sin θ)^{-s} \eta \right],
\qquad
\bar{\eth} \eta = -(\sin θ)^{-s} \left\{ \frac{\partial}{\partial θ} - \frac{i}{\sin θ}
\frac{\partial}{\partial ϕ} \right\} \left[ (\sin θ)^{s} \eta \right],
```
act on the harmonics as
```math
\eth\, {}_sY_{ℓm} = +\sqrt{(ℓ-s)(ℓ+s+1)}\, {}_{s+1}Y_{ℓm},
\qquad
\bar{\eth}\, {}_sY_{ℓm} = -\sqrt{(ℓ+s)(ℓ-s+1)}\, {}_{s-1}Y_{ℓm},
```
which are the relations tested on the [Newman–Penrose page](@ref "Newman-Penrose (1966)").
The article gives the explicit formula
```math
{}_sY_{ℓm}(θ, ϕ)
=
(-1)^{ℓ+m-s}
\sqrt{\frac{(ℓ+m)!\,(ℓ-m)!\,(2ℓ+1)}{4π\,(ℓ+s)!\,(ℓ-s)!}}\;
\sin^{2ℓ}\left(\frac{θ}{2}\right) e^{imϕ}
\sum_r (-1)^r \binom{ℓ-s}{r} \binom{ℓ+s}{r+s-m} \cot^{2r+s-m}\left(\frac{θ}{2}\right),
```
the relation to the Wigner ``D``-matrix
```math
D^ℓ_{-m\,s}(ϕ, θ, -ψ) = (-1)^m \sqrt{\frac{4π}{2ℓ+1}}\; {}_sY_{ℓm}(θ, ϕ)\, e^{isψ},
```
and the conjugation relation ``\overline{{}_sY_{ℓm}} = (-1)^{s+m}\, {}_{-s}Y_{ℓ,-m}``.
Note that the explicit formula differs from [Goldberg et al.'s Eq. (3.1)](@ref "Goldberg et
al. (1967)") by exactly the factor ``(-1)^m``, so we expect it to agree with
[ours](@ref summary_swsh); the ``D``-matrix relation is the integer-index form of our
definition, ``{}_sY_{ℓ,m} = (-1)^m \sqrt{(2ℓ+1)/4π}\; 𝔇^{(ℓ)}_{-m,s}``, extended to
``γ = -ψ``.

## Implementing formulas

We begin by writing code that implements the formulas from Wikipedia.  We encapsulate the
formulas in a module so that we can test them against the `SphericalFunctions` package.
"""

using TestItems: @testitem  #hide
@testitem "Wikipedia conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide

module Wikipedia
#+

# We'll use some predefined utilities to make the code look more like the equations.
import ..ConventionsUtilities: 𝒾, ❗, ∂ⁿ
#+

# Wigner's formula for the ``d``-matrix, as quoted above:
function d(j, m′, m, β::T) where {T<:Real}
    √T((j+m′)❗ * (j-m′)❗ * (j+m)❗ * (j-m)❗) *
    sum(
        (-1)^(m′-m+s) * cos(β/2)^(2j+m-m′-2s) * sin(β/2)^(m′-m+2s)
        / T((j+m-s)❗ * (s)❗ * (m′-m+s)❗ * (j-m′-s)❗)
        for s ∈ max(0, m-m′):min(j+m, j-m′);
        init=zero(T)
    )
end
#+

# The ``D``-matrix, as the matrix element of ``\mathcal{R}(α, β, γ)``:
function D(j, m′, m, α, β, γ)
    exp(-𝒾 * m′ * α) * d(j, m′, m, β) * exp(-𝒾 * m * γ)
end
#+

# The explicit elements listed in the article for ``j=1`` and ``j=2`` (we transcribe the
# rows with ``m' \geq |m|``; the others follow from the symmetries also given in the
# article, which we apply here):
function d_explicit(j, m′, m, β)
    if abs(m′) < abs(m)
        return (-1)^(m-m′) * d_explicit(j, m, m′, β)
    end
    if m′ < 0
        return (-1)^(m-m′) * d_explicit(j, -m′, -m, β)
    end
    if (j, m′, m) == (1, 1, 1)
        (1 + cos(β)) / 2
    elseif (j, m′, m) == (1, 1, 0)
        -sin(β) / √2
    elseif (j, m′, m) == (1, 1, -1)
        (1 - cos(β)) / 2
    elseif (j, m′, m) == (1, 0, 0)
        cos(β)
    elseif (j, m′, m) == (2, 2, 2)
        (1 + cos(β))^2 / 4
    elseif (j, m′, m) == (2, 2, 1)
        -sin(β) * (1 + cos(β)) / 2
    elseif (j, m′, m) == (2, 2, 0)
        √(3/8) * sin(β)^2
    elseif (j, m′, m) == (2, 2, -1)
        -sin(β) * (1 - cos(β)) / 2
    elseif (j, m′, m) == (2, 2, -2)
        (1 - cos(β))^2 / 4
    elseif (j, m′, m) == (2, 1, 1)
        (2cos(β)^2 + cos(β) - 1) / 2
    elseif (j, m′, m) == (2, 1, 0)
        -√(3/8) * sin(2β)
    elseif (j, m′, m) == (2, 1, -1)
        (-2cos(β)^2 + cos(β) + 1) / 2
    elseif (j, m′, m) == (2, 0, 0)
        (3cos(β)^2 - 1) / 2
    else
        error("No explicit formula transcribed for (j, m′, m) = ($j, $m′, $m)")
    end
end
#+

# The Legendre polynomial and the associated Legendre function (without Condon–Shortley
# phase), transcribed literally from Rodrigues' formula using the `∂ⁿ` utility, which
# computes derivatives symbolically:
function P(ℓ, x)
    1 / (2^ℓ * (ℓ)❗) * ∂ⁿ(x -> (x^2 - 1)^ℓ, ℓ)(x)
end
function P̃(ℓ, m, x)
    (1 - x^2)^(m/2) * ∂ⁿ(x -> P(ℓ, x), m)(x)
end
#+

# The spherical harmonics, for ``m \geq 0`` from the definition, and for ``m < 0`` from the
# conjugation relation:
function Y(ℓ, m, θ::T, ϕ::T) where {T<:Real}
    if m < 0
        return (-1)^m * conj(Y(ℓ, -m, θ, ϕ))
    end
    (-1)^m * √T((2ℓ+1) * (ℓ-m)❗ / (4big(π) * (ℓ+m)❗)) * T(P̃(ℓ, m, cos(θ))) * exp(𝒾 * m * ϕ)
end
#+

# The explicit formula for the spin-weighted spherical harmonics.  The sum runs over all
# ``r`` for which the binomial coefficients are nonzero.  The cotangent is singular at the
# north pole, so we will avoid the poles when testing this expression.
function ₛYₗₘ(s, ℓ, m, θ::T, ϕ::T) where {T<:Real}
    (-1)^(ℓ+m-s) * √T((ℓ+m)❗ * (ℓ-m)❗ * (2ℓ+1) / (4big(π) * (ℓ+s)❗ * (ℓ-s)❗)) *
    sin(θ/2)^(2ℓ) * exp(𝒾 * m * ϕ) *
    sum(
        (-1)^r * T(binomial(big(ℓ-s), r) * binomial(big(ℓ+s), r+s-m)) * cot(θ/2)^(2r+s-m)
        for r ∈ max(0, m-s):min(ℓ-s, ℓ+m);
        init=zero(Complex{T})
    )
end
#+

end  # module Wikipedia
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
# normalization differences, which are the most likely source of error.  For the same reason
# we use modest grids of points.  The explicit formula for ``{}_sY_{ℓm}`` involves
# ``\cot(θ/2)``, so we avoid the poles by a small amount.
θϕs = θϕrange(Float64, 7; avoid_poles=1e-3)
αβγs = αβγrange(Float64, 5)
#+

# We will also need the imaginary unit from the utilities module.
import .ConventionsUtilities: 𝒾
#+

# First, Wikipedia's explicit ``d`` elements agree with Wikipedia's general formula, and
# both agree with ours:
for β ∈ βrange(Float64, 15)
    for (j, m′, m) ∈ ℓm′mrange(2)
        j == 0 && continue  # the article lists no j=0 element (it is just 1)
        @test Wikipedia.d_explicit(j, m′, m, β) ≈ Wikipedia.d(j, m′, m, β) atol=ϵₐ rtol=ϵᵣ
    end
    for (j, m′, m) ∈ ℓm′mrange(ℓₘₐₓ)
        @test Wikipedia.d(j, m′, m, β) ≈ ConventionsUtilities.d(j, m′, m, β) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# The ``D``-matrix agrees with ours, and satisfies the conjugation symmetry quoted above:
for (α, β, γ) ∈ αβγs
    for (j, m′, m) ∈ ℓm′mrange(ℓₘₐₓ)
        @test Wikipedia.D(j, m′, m, α, β, γ) ≈ ConventionsUtilities.D(j, m′, m, α, β, γ) atol=ϵₐ rtol=ϵᵣ
        @test Wikipedia.D(j, m′, m, α, β, γ) ≈ (-1)^(m′-m) * conj(Wikipedia.D(j, -m′, -m, α, β, γ)) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# The spherical harmonics agree with ours, and satisfy Wikipedia's relation to ``D``:
for (θ, ϕ) ∈ θϕrange(Float64, 7)
    for (ℓ, m) ∈ ℓmrange(ℓₘₐₓ)
        @test Wikipedia.Y(ℓ, m, θ, ϕ) ≈ ConventionsUtilities.Y(ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
        @test Wikipedia.D(ℓ, m, 0, ϕ, θ, zero(θ)) ≈ √(4π/(2ℓ+1)) * conj(Wikipedia.Y(ℓ, m, θ, ϕ)) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Finally, the spin-weighted spherical harmonics agree with ours, and satisfy Wikipedia's
# relation to the ``D``-matrix and its conjugation symmetry:
for (θ, ϕ) ∈ θϕs
    for (s, ℓ, m) ∈ sℓmrange(ℓₘₐₓ, sₘₐₓ)
        @test Wikipedia.ₛYₗₘ(s, ℓ, m, θ, ϕ) ≈ ConventionsUtilities.Y(s, ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
        @test conj(Wikipedia.ₛYₗₘ(s, ℓ, m, θ, ϕ)) ≈ (-1)^(s+m) * Wikipedia.ₛYₗₘ(-s, ℓ, -m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
        for ψ ∈ (0.0, 0.7, 2.9)
            @test Wikipedia.D(ℓ, -m, s, ϕ, θ, -ψ) ≈
                (-1)^m * √(4π/(2ℓ+1)) * Wikipedia.ₛYₗₘ(s, ℓ, m, θ, ϕ) * exp(𝒾 * s * ψ) atol=ϵₐ rtol=ϵᵣ
        end
    end
end
#+

# These successful tests show that Wikipedia's ``d``, ``D``, ``Y_ℓ^m``, and
# ``{}_sY_{ℓm}`` all agree with the corresponding functions defined by the
# `SphericalFunctions` package.

end  #hide
