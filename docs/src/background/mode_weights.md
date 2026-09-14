# Mode weights

On the [previous page](@ref sYlm_and_Dlmpm), we introduced the
eigenfunctions of [the differential operators](@ref
background_differential_operators) defined on ``\mathrm{Spin}(3)``.
These eigenfunctions are the spin-weighted spherical harmonics (SWSHs)
``{}_{s}Y_{ℓ,m}(R)``, or equivalently Wigner's 𝔇 matrices.

Now that we have introduced the spin-weighted spherical harmonics
(SWSHs) as eigenfunctions of the relevant differential operators, we
can define mode weights of a general spin-weighted function in terms
of these harmonics.

These eigenfunctions form a complete *orthogonal* basis for the space
of square-integrable functions defined on ``\mathrm{Spin}(3)``.  They
are orthogonal rather than orthonormal, because they are normalized
(but not orthogonal) on ``𝕊²`` instead, which leaves ``\left\|
{}_{s}Y_{ℓ,m} \right\|²_{\mathrm{Spin}(3)} = π/2`` (see the [previous
page](@ref Integration-and-normalization)).  Thus, *any*
square-integrable function ``f(R)`` can be expressed as a linear
combination of these harmonics:
```math
f(R) = \sum_{ℓ=0}^{∞} \sum_{s=-ℓ}^{ℓ} \sum_{m=-ℓ}^{ℓ}
{}_{s}f_{ℓ,m}\, {}_{s}Y_{ℓ,m}(R),
```
where the coefficients ``{}_{s}f_{ℓ,m}`` are called the *mode weights*
of the function ``f``.  These mode weights can be computed from the
function using the orthogonality of the SWSHs:
```math
{}_{s}f_{ℓ,m}
= \frac{2}{π} \int_{\mathrm{Spin}(3)} f(R)\, {}_{s}\bar{Y}_{ℓ,m}(R)\, dR,
```
where the factor in front of the integral is the reciprocal of
``\left\| {}_{s}Y_{ℓ,m} \right\|²_{\mathrm{Spin}(3)} = π/2``, which is
explained [here](@ref Integration-and-normalization).  Note that we
have not restricted the spin weight ``s`` of the function ``f``; a
general function on ``\mathrm{Spin}(3)`` can have contributions from
SWSHs of any spin weight, so ``s`` was included in the sum above.  In
fact, those spin weights may have half-integral values as well, in
which case the sum over ``ℓ`` must include all positive half-integral
values as well as the integral values.

However, if we restrict to functions with a *specific* spin weight
``s``, only the SWSHs with *that* spin weight will contribute to the
expansion, and we can simplify the expressions above to
```math
f(R) = \sum_{ℓ=|s|}^{∞} \sum_{m=-ℓ}^{ℓ}
{}_{s}f_{ℓ,m}\, {}_{s}Y_{ℓ,m}(R)
```
and
```math
{}_{s}f_{ℓ,m}
= \frac{2}{π} \int_{\mathrm{Spin}(3)} f(R)\, {}_{s}\bar{Y}_{ℓ,m}(R)\, dR.
```

In practice the integral is almost always written over the sphere
``𝕊²`` instead, but the relationship is simple.  The integrand ``f\,
{}_{s}\bar{Y}_{ℓ,m}`` is the product of a function of spin weight
``s`` with the conjugate of another, so *that product* has spin weight
``0.``  A function of spin weight 0 does not depend on the third Euler
angle at all, so it [pushes forward](@ref "Pushing forward to
``𝕊²``") to an ordinary function on ``𝕊²``, and — by exactly the
argument that relates the two norms [above](@ref
Integration-and-normalization) — its integral over
``\mathrm{Spin}(3)`` is ``π/2`` times its integral over ``𝕊²``.  The
factor ``2/π`` therefore cancels, leaving the familiar expression
```math
{}_{s}f_{ℓ,m}
= \int_{𝕊²} f(θ, ϕ)\, {}_{s}\bar{Y}_{ℓ,m}(θ, ϕ)\, dΩ,
\qquad dΩ = \sin θ\, dθ\, dϕ,
```
in which ``f(θ, ϕ)`` means ``f`` evaluated at the rotor that
[spherical coordinates](@ref "Pushing forward to ``𝕊²``") assign to
the point ``(θ, ϕ)``.  Note that this step relies on ``f`` having a
definite spin weight; for a function with contributions from several
spin weights — where the third Euler angle genuinely matters — the
integral over ``\mathrm{Spin}(3)`` in the previous expression is the
one to use.  Functions whose domain is really ``I × 𝕊¹`` rather than
``𝕊²`` are discussed [here](@ref "Pulling back to ``I×𝕊¹``").

Half-integer spin weight is *not* an exception here: the cancellation
above only needs ``f`` and ``{}_{s}Y_{ℓ,m}`` to have the *same* spin
weight, so the product is ``γ``-independent — and in particular
unchanged under ``𝐑 → -𝐑`` — whatever the parity of ``s`` may be.
What is peculiar to half-integer ``s`` is that each factor separately
is double-valued on ``𝕊²``, changing sign under ``ϕ → ϕ + 2π``; the
integral must therefore be taken with one fixed rotor assignment ``(θ,
ϕ) ↦ 𝐐``, under which the integrand is single-valued and the identity
holds unchanged.


## Differential operators

Mode weights and functions transform in related but different ways
under the action of the differential operators.

One important point to note is that mode weights transform
"contravariantly" (very loosely speaking) relative to the
spin-weighted spherical functions under some operators.  For example,
take the action of the ``L_+`` operator, which acts on a SWSH as
```math
L_+ \left\{{}_{s}Y_{ℓ,m}\right\} (R)
= \sqrt{(ℓ-m)(ℓ+m+1)}\ {}_{s}Y_{ℓ,m+1}(R).
```
We can use this to derive mode weights of a general spin-weighted
function ``f`` under the action of this operator:[^1]
```math
\begin{aligned}
\left\{L_+ f\right\}_{ℓ,m}
&=
\int \left\{L_+ f(R)\right\}\, {}_{s}\bar{Y}_{ℓ,m}(R)\, dR \\
&=
\int \left\{L_+ \sum_{ℓ',m'}f_{ℓ',m'}\, {}_{s}Y_{ℓ',m'}(R)\right\}\, {}_{s}\bar{Y}_{ℓ,m}(R)\, dR \\
&=
\int \sum_{ℓ',m'} f_{ℓ',m'}\, \left\{L_+ {}_{s}Y_{ℓ',m'}(R)\right\}\, {}_{s}\bar{Y}_{ℓ,m}(R)\, dR \\
&=
\sum_{ℓ',m'} f_{ℓ',m'}\, \int \left\{L_+ {}_{s}Y_{ℓ',m'}(R)\right\}\, {}_{s}\bar{Y}_{ℓ,m}(R)\, dR \\
&=
\sum_{ℓ',m'} f_{ℓ',m'}\, \int \left\{\sqrt{(ℓ'-m')(ℓ'+m'+1)} {}_{s}Y_{ℓ',m'+1}(R)\right\}\, {}_{s}\bar{Y}_{ℓ,m}(R)\, dR \\
&=
\sum_{ℓ',m'} f_{ℓ',m'}\, \sqrt{(ℓ'-m')(ℓ'+m'+1)} \int {}_{s}Y_{ℓ',m'+1}(R)\, {}_{s}\bar{Y}_{ℓ,m}(R)\, dR \\
&=
\sum_{ℓ',m'} f_{ℓ',m'}\, \sqrt{(ℓ'-m')(ℓ'+m'+1)} δ_{ℓ,ℓ'} δ_{m,m'+1} \\
&=
f_{ℓ,m-1}\, \sqrt{(ℓ-m+1)(ℓ+m)}
\end{aligned}
```
Note that this expression (and in particular its signs) more resembles
the expression for ``L_- \left\{{}_{s}Y_{ℓ,m}\right\}`` than for
``L_+ \left\{{}_{s}Y_{ℓ,m}\right\}``.  Similar relations hold for
the action of ``L_-``.

[^1]:
    A technical note about the integrals above: the integrals should
    be taken over the appropriate space and with the appropriate
    weight such that the SWSHs are orthonormal.  In general, this
    integral should be over ``\mathrm{Spin}(3)`` and weighted by
    ``2/π`` so that the result will be either ``0`` or ``1``; in
    general the SWSHs are not truly orthonormal when integrated over
    an ``𝕊²`` subspace (nor even is the integral invariant).
    However, if we know that the spins are the same in both cases, it
    *is* possible to integrate over an ``𝕊²`` subspace.

However, it is important to note that the same "contravariance" is not
present for the spin-raising and -lowering operators:
```math
\begin{aligned}
\left\{\eth f\right\}_{s+1,ℓ,m}
&=
\int \left\{\eth f(R)\right\}\, {}_{s+1}\bar{Y}_{ℓ,m}(R)\, dR \\
&=
\int \left\{\eth \sum_{ℓ',m'}f_{ℓ',m'}\, {}_{s}Y_{ℓ',m'}(R)\right\}\, {}_{s+1}\bar{Y}_{ℓ,m}(R)\, dR \\
&=
\sum_{ℓ',m'} f_{ℓ',m'}\, \int \left\{\eth {}_{s}Y_{ℓ',m'}(R)\right\}\, {}_{s+1}\bar{Y}_{ℓ,m}(R)\, dR \\
&=
\sum_{ℓ',m'} f_{ℓ',m'}\, \sqrt{(ℓ'-s)(ℓ'+s+1)} \int {}_{s+1}Y_{ℓ',m'}(R)\, {}_{s+1}\bar{Y}_{ℓ,m}(R)\, dR \\
&=
\sum_{ℓ',m'} f_{ℓ',m'}\, \sqrt{(ℓ'-s)(ℓ'+s+1)} δ_{ℓ,ℓ'} δ_{m,m'} \\
&=
\left\{f\right\}_{s,ℓ,m}\, \sqrt{(ℓ-s)(ℓ+s+1)}
\end{aligned}
```
The operators ``R_\pm`` obey this same, more "covariant" form of
transformation, and so does ``\bar{\eth}`` — except that the latter
includes Newman and Penrose's minus sign, ``\bar{\eth} = -R_-``, so
that
```math
\left\{\bar{\eth} f\right\}_{s-1,ℓ,m}
= -\sqrt{(ℓ+s)(ℓ-s+1)}\, \left\{f\right\}_{s,ℓ,m}.
```
See the [conventions summary](@ref summary_spin_weight) for the
identification ``\eth = R_+`` and ``\bar{\eth} = -R_-``.
