md"""
# Goldberg et al. (1967)

!!! info "Summary"
    Goldberg et al. use the same conventions as we do for spherical
    coordinates and Euler angles.  Their Wigner ``D`` matrix is
    defined as a function of the *inverse* rotation, and is related to
    ours by
    ```math
    D^{j}_{m',m}(α, β, γ)\big|_{\text{Goldberg}} =
    \overline{𝔇^{(j)}_{m',m}(γ, β, α)} = (-1)^{m+m'}\,
    \overline{𝔇^{(j)}_{m,m'}(α, β, γ)}.
    ```
    Their spin-raising and -lowering operators ``\eth`` and
    ``\bar{\eth}`` agree with ours, including the minus sign in
    ``\bar{\eth}\, {}_sY_{ℓ,m} = -\sqrt{(ℓ+s)(ℓ-s+1)}\,
    {}_{s-1}Y_{ℓ,m}``.

[GoldbergEtAl_1967](@citet) presented the first paper specifically
about spin-weighted spherical harmonics (after [Newman_1966](@citet)
introduced them; see [the previous page](@ref "Newman-Penrose
(1966)")), and the first to relate them to the Wigner D-matrices.
They also used functional symmetry to introduce what we recognize
geometrically as the right-Lie derivative operator.

## Coordinates and derivative operators

Goldberg et al. start out by saying that they work in
"three-dimensional Euclidean space with polar  coordinates ``r, θ,
ϕ``".  They write down Newman-Penrose's operator acting on a field of
spin weight ``s`` as
```math
\eth \eta = - (\sin\theta)^s \left[ \frac{\partial}{\partial\theta}
+i \csc\theta \frac{\partial}{\partial\phi} \right] (\sin\theta)^{-s}
\eta.
```

They then "define a rotation ``R(α, β, γ)`` of Euler angles ``α, β,
γ`` as being composed of ``γ`` about ``OZ`` followed by ``β`` about
``OY`` and then ``α`` about ``OZ``", noting that "this procedure is
clearly equivalent to the more usual one of a  rotation ``γ`` around
``OZ``, followed by ``β`` around ``OY'`` and finally ``α`` around
``OZ''``."  This agrees precisely with our definition.

They then go on to define operators "familiar from the theory of the
symmetric top".  Unfortunately, this section contains a small mess of
sign errors.  Their *corrected* Eqs. (3.16) read
```math
L_z = -i \frac{\partial}{\partial\alpha},
\qquad
L_{\pm} = \pm e^{\pm i \alpha} \left(
    \frac{\partial}{\partial\beta}
    \pm i \cot \beta \frac{\partial}{\partial\alpha}
    \mp i \csc \beta \frac{\partial}{\partial\gamma}
\right),
```
where the correction is to flip the sign of that last term in
parentheses.  With that correction, we have perfect agreement between
[this package's ``L`` operators in Euler coordinates](@ref euler_L_S3)
and Goldberg's.

They then construct a second set of operators by swapping ``α`` and
``-γ``, though even that statement has a typo in the original, and the
swapping is done incorrectly.  Their *corrected* Eqs. (3.18) read
```math
K_z = i \frac{\partial}{\partial\gamma},
\qquad
K_{\pm} = \pm e^{\mp i \gamma} \left(
    \frac{\partial}{\partial\beta}
    \mp i \cot \beta \frac{\partial}{\partial\gamma}
    \pm i \csc \beta \frac{\partial}{\partial\alpha}
\right),
```
where the corrections are to flip the sign in the exponent and in the
middle term in parentheses.  We those corrections, we still have some
strange disagreement.  With ``R`` as [our right-Lie derivative](@ref
euler_R_S3), we have
```math
K_z = R_z, \qquad K_{\pm} = -R_{\pm}.
```

It's also worth noting that there is an incorrect factor of ``i`` in
the upper line of Eq. (3.20), but the second line is correct.  They
then find the relationship between ``K_+`` and ``\eth`` — but have yet
another sign error.  Their *corrected* [Eq. (3.21)] reads
```math
[K_+ D^{l}_{-sm}]_{\alpha=\phi, \beta=\theta, \gamma=0}
=
-\eth D^{l}_{-sm} (\phi\theta 0),
```
where the correction is the sign on the right-hand side.  Including
this corrected sign, and remembering that ``K_+ = -R_+``, we see that
this agrees with [this package's relationship ``R_+ = \eth``](@ref
summary_spin_weight).

## Spin-weighted spherical harmonics

They give the explicit expression [Eq. (3.1)]
```math
\begin{aligned}
{}_sY_{ℓ,m}(θ, ϕ)
&=
\left[ \frac{(ℓ+m)!\,(ℓ-m)!\,(2ℓ+1)}{(ℓ+s)!\,(ℓ-s)!\,4π} \right]^{1/2}
\left(\sin\tfrac{θ}{2}\right)^{2ℓ} \\
&\times
\sum_r \binom{ℓ-s}{r} \binom{ℓ+s}{r+s-m}
(-1)^{ℓ-r-s}\, e^{imϕ} \left(\cot\tfrac{θ}{2}\right)^{2r+s-m},
\end{aligned}
```
where the sum runs over all ``r`` for which the binomial coefficients are nonzero.  They
also record the conjugation symmetry [Eq. (2.6), though with a missing minus sign on the
last index]
```math
\overline{{}_sY_{ℓ,m}} = (-1)^{m+s}\, {}_{-s}Y_{ℓ,-m},
```
and the action of the Newman–Penrose operators, with the minus sign in the lowering
relation [Eqs. (2.7a) and (2.7b)]
```math
\begin{gather}
\eth\, {}_sY_{ℓ,m} = \sqrt{(ℓ-s)(ℓ+s+1)}\, {}_{s+1}Y_{ℓ,m},
\\
\bar{\eth}\, {}_sY_{ℓ,m} = -\sqrt{(ℓ+s)(ℓ-s+1)}\, {}_{s-1}Y_{ℓ,m},
\end{gather}
```
so that ``\bar{\eth}\eth\, {}_sY_{ℓ,m} = -(ℓ-s)(ℓ+s+1)\, {}_sY_{ℓ,m}``.  These are exactly
the relations on [our summary page](@ref summary_swsh) — but the explicit expression (3.1)
differs from [ours](@ref summary_swsh) by a factor of ``(-1)^m``, which the tests below
confirm.  Since ``(-1)^m`` is independent of ``s``, it does not affect the ladder relations.

## Wigner's ``D`` matrix

If we relate two vectors by a rotation matrix as ``x'^k = R^{kl} x^l``, then Goldberg et
al. define ``D`` by its action on spherical harmonics [Eq. (3.3)]
```math
Y_{ℓ,m}(x') = \sum_{m'} Y_{ℓ,m'}(x) D^{ℓ}_{m',m}\left( R^{-1} \right).
```
They then define the Euler angles as we do, and write [Eq. (3.4)]
```math
D^{ℓ}_{m', m}(α, β, γ)
\equiv
D^{ℓ}_{m', m}\left( R(α β γ)^{-1} \right)
=
e^{i m' γ} d^{ℓ}_{m', m}(β) e^{i m α}.
```
Note the two differences from [our definition](@ref summary_wigner_D): the argument is the
*inverse* rotation, and the phases are ``e^{+im'γ}`` and ``e^{+imα}`` rather than
``e^{-im'α}`` and ``e^{-imγ}``.  Finally, they derive [Eq. (3.9)]
```math
D^{j}_{m', m}(α, β, γ)
=
\left[\frac{(j+m)!(j-m)!}{(j+m')!(j-m')!}\right]^{1/2}
(\sin \tfrac{1}{2}β)^{2j}
\sum_r \binom{j+m'}{r} \binom{j-m'}{r-m-m'}
(-1)^{j+m'-r}
e^{imα}
(\cot \tfrac{1}{2}β)^{2r-m-m'}
e^{im'γ}.
```
Comparing to our ``𝔇``, we find (and test below) that
```math
D^{j}_{m',m}(α, β, γ)\big|_{\text{Goldberg}}
= \overline{𝔇^{(j)}_{m',m}(γ, β, α)}
= (-1)^{m+m'}\, \overline{𝔇^{(j)}_{m,m'}(α, β, γ)}
= 𝔇^{(j)}_{m,m'}(-γ, -β, -α),
```
where the second and third forms follow from [the symmetries of
``𝔇``](@ref summary_wigner_D).  That is, their ``D`` is the complex conjugate of ours with
the roles of the first and last Euler angles swapped — as expected from a definition in
terms of the inverse rotation ``R(αβγ)^{-1} = R(-γ, -β, -α)``.

Finally, they relate the two quantities [Eq. (3.11)]:
```math
  {}_sY_{ℓ, m}(θ, ϕ)
  =
  \left[ \left(2ℓ+1\right) / 4π \right]^{1/2}
  D^{ℓ}_{-s,m}(ϕ, θ, 0).
```
This naturally extends to
```math
  {}_sY_{ℓ, m}(θ, ϕ, γ)
  =
  \left[ \left(2ℓ+1\right) / 4π \right]^{1/2}
  D^{ℓ}_{-s,m}(ϕ, θ, γ),
```
where Eq. (3.4) shows that ``D^{ℓ}_{m', m}(α, β, γ) = D^{ℓ}_{m', m}(α, β, 0) e^{i m' γ}``,
so we have
```math
  {}_sY_{ℓ, m}(θ, ϕ, γ)
  =
  {}_sY_{ℓ, m}(θ, ϕ)\, e^{-i s γ}.
```
This is the most natural extension of the standard spin-weighted spherical harmonics to
``\mathrm{Spin}(3)``, and it is precisely [our definition of spin weight](@ref
summary_spin_weight): the spin-weight operator is ``i \partial_γ = R_z``.

## Implementing formulas

We begin by writing code that implements the formulas from Goldberg et al.  We encapsulate
the formulas in a module so that we can test them against the `SphericalFunctions` package.
"""

using TestItems: @testitem  #hide
@testitem "Goldberg et al. conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide

module GoldbergEtAl
#+

# We'll use some predefined utilities to make the code look more like the equations, and
# `ForwardDiff` to evaluate the derivatives in ``\eth`` and ``\bar{\eth}``.
import ..ConventionsUtilities: 𝒾, ❗
import ForwardDiff
#+

# Equation (3.1) gives the spin-weighted spherical harmonics.  The sum runs over all ``r``
# for which the binomial coefficients are nonzero, which means ``\max(0, m-s) \leq r \leq
# \min(ℓ-s, ℓ+m)``.  Note that the cotangent is singular at the north pole, so we will
# avoid the poles when testing this expression.
function ₛYₗₘ(s, ℓ, m, θ, ϕ)
    T = float(promote_type(typeof(θ), typeof(ϕ)))
    √T(((ℓ+m)❗ * (ℓ-m)❗ * (2ℓ+1)) / ((ℓ+s)❗ * (ℓ-s)❗ * 4big(π))) * sin(θ/2)^(2ℓ) *
    sum(
        T(binomial(big(ℓ-s), r) * binomial(big(ℓ+s), r+s-m)) * (-1)^(ℓ-r-s)
        * exp(𝒾 * m * ϕ) * cot(θ/2)^(2r+s-m)
        for r ∈ max(0, m-s):min(ℓ-s, ℓ+m);
        init=zero(Complex{T})
    )
end
#+

# Equation (3.9) gives the ``D`` matrix.  Again, the sum runs over all ``r`` for which the
# binomial coefficients are nonzero, ``\max(0, m+m') \leq r \leq \min(j+m', j+m)``, and
# again the cotangent is singular at ``β = 0``, so we will avoid the poles when testing.
function D(j, m′, m, α, β, γ)
    T = float(promote_type(typeof(α), typeof(β), typeof(γ)))
    √T(((j+m)❗ * (j-m)❗) / ((j+m′)❗ * (j-m′)❗)) * sin(β/2)^(2j) *
    sum(
        T(binomial(big(j+m′), r) * binomial(big(j-m′), r-m-m′)) * (-1)^(j+m′-r)
        * exp(𝒾 * m * α) * cot(β/2)^(2r-m-m′) * exp(𝒾 * m′ * γ)
        for r ∈ max(0, m+m′):min(j+m′, j+m);
        init=zero(Complex{T})
    )
end
#+

# The operators ``\eth`` and ``\bar{\eth}`` are those of Newman and Penrose; we implement
# them exactly as on [the previous page](@ref "Newman-Penrose (1966)").
function ð(η, s)
    function ðη(θ, ϕ)
        η̃(θ, ϕ) = sin(θ)^(-s) * η(θ, ϕ)
        ∂θ = ForwardDiff.derivative(θ′ -> η̃(θ′, ϕ), θ)
        ∂ϕ = ForwardDiff.derivative(ϕ′ -> η̃(θ, ϕ′), ϕ)
        -sin(θ)^s * (∂θ + 𝒾 / sin(θ) * ∂ϕ)
    end
end
function ð̄(η, s)
    function ð̄η(θ, ϕ)
        η̃(θ, ϕ) = sin(θ)^(s) * η(θ, ϕ)
        ∂θ = ForwardDiff.derivative(θ′ -> η̃(θ′, ϕ), θ)
        ∂ϕ = ForwardDiff.derivative(ϕ′ -> η̃(θ, ϕ′), ϕ)
        -sin(θ)^(-s) * (∂θ - 𝒾 / sin(θ) * ∂ϕ)
    end
end
#+

end  # module GoldbergEtAl
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
# reason we use modest grids of points.  Because Goldberg et al.'s expressions involve
# ``\cot(θ/2)`` and ``\cot(β/2)``, they are singular at ``θ = 0`` and ``β = 0`` (where the
# limits are finite), so we avoid the poles by a small amount:
θϕs = θϕrange(Float64, 7; avoid_poles=1e-3)
αβγs = αβγrange(Float64, 5; avoid_poles=1e-3)
#+

# First, the internal consistency of Goldberg et al.'s own expressions: the conjugation
# symmetry of Eq. (2.6), and the relation (3.11) between ``{}_sY_{ℓ,m}`` and ``D``:
for (θ, ϕ) ∈ θϕs
    for (s, ℓ, m) ∈ sℓmrange(ℓₘₐₓ, sₘₐₓ)
        @test conj(GoldbergEtAl.ₛYₗₘ(s, ℓ, m, θ, ϕ)) ≈
            (-1)^(m+s) * GoldbergEtAl.ₛYₗₘ(-s, ℓ, -m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
        @test GoldbergEtAl.ₛYₗₘ(s, ℓ, m, θ, ϕ) ≈
            √((2ℓ+1) / (4π)) * GoldbergEtAl.D(ℓ, -s, m, ϕ, θ, zero(θ)) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Next, we compare their spin-weighted spherical harmonics to ours, and find the factor of
# ``(-1)^m``:
for (θ, ϕ) ∈ θϕs
    for (s, ℓ, m) ∈ sℓmrange(ℓₘₐₓ, sₘₐₓ)
        @test GoldbergEtAl.ₛYₗₘ(s, ℓ, m, θ, ϕ) ≈
            (-1)^m * ConventionsUtilities.Y(s, ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Now the ``D`` matrix.  We test both forms of the relation given above: the conjugate of
# ours with ``α`` and ``γ`` swapped, and ``(-1)^{m+m'}`` times the conjugate transpose.
for (α, β, γ) ∈ αβγs
    for (j, m′, m) ∈ ℓm′mrange(ℓₘₐₓ)
        @test GoldbergEtAl.D(j, m′, m, α, β, γ) ≈
            conj(ConventionsUtilities.D(j, m′, m, γ, β, α)) atol=ϵₐ rtol=ϵᵣ
        @test GoldbergEtAl.D(j, m′, m, α, β, γ) ≈
            (-1)^(m+m′) * conj(ConventionsUtilities.D(j, m, m′, α, β, γ)) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# The extension to nonzero ``γ`` just multiplies by ``e^{-isγ}``, the defining property of
# a function of spin weight ``s``:
for (α, β, γ) ∈ αβγs
    for (s, ℓ, m) ∈ sℓmrange(ℓₘₐₓ, sₘₐₓ)
        @test GoldbergEtAl.D(ℓ, -s, m, α, β, γ) ≈
            GoldbergEtAl.D(ℓ, -s, m, α, β, zero(γ)) * exp(-im * s * γ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Finally, the ladder relations of Eqs. (2.7), applied to their own ``{}_sY_{ℓ,m}``.  The
# operators involve ``1/\sin θ``, so we avoid the poles, and — as on the Newman–Penrose
# page — evaluate in `BigFloat` arithmetic to keep the rounding errors amplified by
# ``\sin^{-s} θ`` below `Float64` precision.  When the raised or lowered spin weight would
# exceed ``ℓ`` in magnitude, the result must vanish.
for (θ, ϕ) ∈ ((big(θ), big(ϕ)) for (θ, ϕ) ∈ θϕrange(Float64, 5; avoid_poles=1e-3))
    for (s, ℓ, m) ∈ sℓmrange(ℓₘₐₓ, sₘₐₓ)
        Y(θ, ϕ) = GoldbergEtAl.ₛYₗₘ(s, ℓ, m, θ, ϕ)
        ðY = GoldbergEtAl.ð(Y, s)(θ, ϕ)
        ð̄Y = GoldbergEtAl.ð̄(Y, s)(θ, ϕ)
        if abs(s+1) ≤ ℓ
            @test ðY ≈ √((ℓ-s) * (ℓ+s+1)) * GoldbergEtAl.ₛYₗₘ(s+1, ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
        else
            @test ðY ≈ 0 atol=ϵₐ
        end
        if abs(s-1) ≤ ℓ
            @test ð̄Y ≈ -√((ℓ+s) * (ℓ-s+1)) * GoldbergEtAl.ₛYₗₘ(s-1, ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
        else
            @test ð̄Y ≈ 0 atol=ϵₐ
        end
    end
end
#+

# These successful tests show that Goldberg et al.'s spin-weighted spherical harmonics are
# ``(-1)^m`` times ours, that their ``D`` matrix is the complex conjugate of ours with the
# first and last Euler angles interchanged, and that their ``\eth`` and ``\bar{\eth}``
# relations hold exactly as stated on our conventions pages.

end  #hide
