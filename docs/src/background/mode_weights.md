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

These eigenfunctions form a complete orthonormal basis for the space
of square-integrable functions defined on ``\mathrm{Spin}(3)``.  Thus,
*any* such function ``f(R)`` can be expressed as a linear combination
of these harmonics:
```math
f(R) = \sum_{ℓ=0}^{∞} \sum_{s=-ℓ}^{ℓ} \sum_{m=-ℓ}^{ℓ}
{}_{s}f_{ℓ,m}\, {}_{s}Y_{ℓ,m}(R),
```
where the coefficients ``{}_{s}f_{ℓ,m}`` are called the *mode weights*
of the function ``f``.  These mode weights can be computed from the
function using the orthogonality of the SWSHs:
```math
{}_{s}f_{ℓ,m}
= \frac{π}{2} \int_{\mathrm{Spin}(3)} f(R)\, {}_{s}\bar{Y}_{ℓ,m}(R)\, dR,
```
where the factor in front of the integral is explained [here](@ref
Integration-and-normalization).  Note that we have not restricted the
spin weight ``s`` of the function ``f``; a general function on
``\mathrm{Spin}(3)`` can have contributions from SWSHs of any spin
weight, so ``s`` was included in the sum above.  In fact, those spin
weights may have half-integral values as well, in which case the sum
over ``ℓ`` must include all positive half-integral values as well as
the integral values.

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
= \frac{π}{2} \int_{\mathrm{Spin}(3)} f(R)\, {}_{s}\bar{Y}_{ℓ,m}(R)\, dR.
```

!!! danger "#TODO"

    Finish this section, citing the section on [functions 
    of spherical coordinates](@ref Pulling-back-to-I𝕊) and
    reducing to the integral over ``𝕊²``.  Refine the 
    following section.


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
function ``f`` under the action of this operator:[^2]
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

[^2]:
    A technical note about the integrals above: the integrals should be taken
    over the appropriate space and with the appropriate weight such that the
    SWSHs are orthonormal.  In general, this integral should be over
    ``\mathrm{Spin}(3)`` and weighted by ``1/2π`` so that the result will be
    either ``0`` or ``1``; in general the SWSHs are not truly orthonormal when
    integrated over an ``𝕊²`` subspace (nor even is the integral invariant).
    However, if we know that the spins are the same in both cases, it *is*
    possible to integrate over an ``𝕊²`` subspace.

However, it is important to note that the same "contravariance" is not
present for the spin-raising and -lowering operators:
```math
\begin{aligned}
\left\{\eth f\right\}_{ℓ,m}
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
f_{ℓ,m}\, \sqrt{(ℓ-s)(ℓ+s+1)}
\end{aligned}
```
Similarly ``\bar{\eth}`` — and ``R_\pm`` of course — obey this more
"covariant" form of transformation.

