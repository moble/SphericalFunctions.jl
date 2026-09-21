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

The same result holds for half-integer spins: the cancellation above
only needs ``f`` and ``{}_{s}Y_{ℓ,m}`` to have the *same* spin weight,
so the product is ``γ``-independent — and in particular unchanged
under ``𝐑 → -𝐑`` — whatever the parity of ``s`` may be.  What is
peculiar to half-integer ``s`` is that each factor separately is
double-valued on ``𝕊²``, changing sign under ``ϕ → ϕ + 2π``; the
integral must therefore be taken with one fixed rotor assignment ``(θ,
ϕ) ↦ 𝐐``, under which the integrand is single-valued and the identity
holds unchanged.


## Differential operators

Mode weights and functions transform in related but different ways
under the action of the differential operators.

One important point to note is that mode weights transform
"oppositely" (loosely speaking) relative to the spin-weighted
spherical functions under some operators.  For example, take the
action of the ``L_+`` operator, which acts on a SWSH as
```math
L_+ \left\{{}_{s}Y_{ℓ,m}\right\} (R)
= \sqrt{(ℓ-m)(ℓ+m+1)}\ {}_{s}Y_{ℓ,m+1}(R).
```
We can use this to derive mode weights of a general spin-weighted
function ``f`` under the action of this operator:
```math
\begin{aligned}
\left\{L_+ f\right\}_{ℓ,m}
&=
\frac{2}{\pi} \int \left\{L_+ f(R)\right\}\, {}_{s}\bar{Y}_{ℓ,m}(R)\, dR \\
&=
\frac{2}{\pi} \int \left\{L_+ \sum_{ℓ',m'}f_{ℓ',m'}\, {}_{s}Y_{ℓ',m'}(R)\right\}\, {}_{s}\bar{Y}_{ℓ,m}(R)\, dR \\
&=
\frac{2}{\pi} \int \sum_{ℓ',m'} f_{ℓ',m'}\, \left\{L_+ {}_{s}Y_{ℓ',m'}(R)\right\}\, {}_{s}\bar{Y}_{ℓ,m}(R)\, dR \\
&=
\sum_{ℓ',m'} f_{ℓ',m'}\, \frac{2}{\pi} \int \left\{L_+ {}_{s}Y_{ℓ',m'}(R)\right\}\, {}_{s}\bar{Y}_{ℓ,m}(R)\, dR \\
&=
\sum_{ℓ',m'} f_{ℓ',m'}\, \frac{2}{\pi} \int \left\{\sqrt{(ℓ'-m')(ℓ'+m'+1)} {}_{s}Y_{ℓ',m'+1}(R)\right\}\, {}_{s}\bar{Y}_{ℓ,m}(R)\, dR \\
&=
\sum_{ℓ',m'} f_{ℓ',m'}\, \sqrt{(ℓ'-m')(ℓ'+m'+1)} \frac{2}{\pi} \int {}_{s}Y_{ℓ',m'+1}(R)\, {}_{s}\bar{Y}_{ℓ,m}(R)\, dR \\
&=
\sum_{ℓ',m'} f_{ℓ',m'}\, \sqrt{(ℓ'-m')(ℓ'+m'+1)} δ_{ℓ,ℓ'} δ_{m,m'+1} \\
&=
f_{ℓ,m-1}\, \sqrt{(ℓ-m+1)(ℓ+m)}
\end{aligned}
```
Note that this expression (and in particular its signs) more resembles
the expression for ``L_- \left\{{}_{s}Y_{ℓ,m}\right\}`` than for
``L_+ \left\{{}_{s}Y_{ℓ,m}\right\}``.  Similar relations hold for
the action of ``L_-``:
```math
\begin{gathered}
L_- \left\{{}_{s}Y_{ℓ,m}\right\} (R)
= \sqrt{(ℓ+m)(ℓ-m+1)}\ {}_{s}Y_{ℓ,m-1}(R),
\\
\left\{L_- f\right\}_{ℓ,m}
= f_{ℓ,m+1}\, \sqrt{(ℓ+m+1)(ℓ-m)}.
\end{gathered}
```
Seeing this, it is obvious that this "opposite" transformation is just
down to the fact that ``L_+`` and ``L_-`` are dual operators, and we
have
```math
\left\langle {}_{s}Y_{ℓ,m} \middle| L_+ f \right\rangle
= \left\langle L_- {}_{s}Y_{ℓ,m} \middle| f \right\rangle.
``` 

However, it is important to note that the same duality is not
*apparent* for the spin-raising and -lowering operator laws:
```math
\begin{aligned}
\left\{\eth f\right\}_{s+1,ℓ,m}
&=
\frac{2}{\pi} \int \left\{\eth f(R)\right\}\, {}_{s+1}\bar{Y}_{ℓ,m}(R)\, dR \\
&=
\frac{2}{\pi} \int \left\{\eth \sum_{ℓ',m'}f_{ℓ',m'}\, {}_{s}Y_{ℓ',m'}(R)\right\}\, {}_{s+1}\bar{Y}_{ℓ,m}(R)\, dR \\
&=
\sum_{ℓ',m'} f_{ℓ',m'}\, \frac{2}{\pi} \int \left\{\eth {}_{s}Y_{ℓ',m'}(R)\right\}\, {}_{s+1}\bar{Y}_{ℓ,m}(R)\, dR \\
&=
\sum_{ℓ',m'} f_{ℓ',m'}\, \sqrt{(ℓ'-s)(ℓ'+s+1)} \frac{2}{\pi} \int {}_{s+1}Y_{ℓ',m'}(R)\, {}_{s+1}\bar{Y}_{ℓ,m}(R)\, dR \\
&=
\sum_{ℓ',m'} f_{ℓ',m'}\, \sqrt{(ℓ'-s)(ℓ'+s+1)} δ_{ℓ,ℓ'} δ_{m,m'} \\
&=
\left\{f\right\}_{s,ℓ,m}\, \sqrt{(ℓ-s)(ℓ+s+1)}.
\end{aligned}
```
The reason for this apparent asymmetry between ``L_\pm`` and ``\eth``
is actually just an asymmetry in the way we treat modes as they vary
over ``m`` and over ``s``, belied by a sleight of hand we played in
the notation above.  We generally assemble mode weights for a *single*
spin weight ``s``, while varying the index ``m`` over its full range.
So the expression for ``\left\{L_+ f\right\}_{ℓ,m}`` has ``m`` on the
left-hand side and ``m-1`` on the right-hand side, whereas the
expression for ``\left\{\eth f\right\}_{s+1,ℓ,m}`` has ``s+1`` on the
left-hand side and ``s`` on the right-hand side.  That is, we assume
that ``f`` has a fixed spin weight ``s``, so we don't care about the
``s`` component of ``\eth f`` — it is automatically zero.  On the
other hand, we assume that ``f`` has a range of ``m`` components, so
we *do* care about the ``m`` component of ``L_+ f``.  Perhaps a
simpler way to see this is to write the modes as vectors, and look at
the operators in matrix form, as we do in the next section.

# Vector/matrix forms

We conventionally assemble the mode weights of a function ``f`` with a
fixed spin weight ``s`` into a single vector of data, with ``ℓ``
increasing, but ``m`` running from ``-ℓ`` to ``ℓ`` between each
increment of ``ℓ``.  Schematically, this looks like
```math
f_{ℓ,m} \leftrightarrow
[f_{0,0}, f_{1,-1}, f_{1,0}, f_{1,1}, f_{2,-2}, f_{2,-1}, f_{2,0}, f_{2,1}, f_{2,2}, \ldots]^T.
```
This is the form used throughout this package — and in particular the
[`ModeWeights`](@ref) type.  (Optionally, the ``ℓ < |s|`` modes can be
omitted since they are automatically zero, but on this page we will
assume they are present, for simplicity.)

The differential operators can then be represented as matrices acting
on these vectors.  For example, ``L_z`` is diagonal in this basis,
with the ``m`` values along the diagonal:
```math
L_z = \begin{pmatrix}
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \cdots \\
0 & -1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \cdots \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \cdots \\
0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & \cdots \\
0 & 0 & 0 & 0 & -2 & 0 & 0 & 0 & 0 & \cdots \\
0 & 0 & 0 & 0 & 0 & -1 & 0 & 0 & 0 & \cdots \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \cdots \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & \cdots \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 2 & \cdots \\
\vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \ddots
\end{pmatrix}.
```
The ``L_+`` operator is represented by a matrix that is *not*
diagonal, but rather has nonzero entries just below the diagonal, with
the square-root factors from the previous section:
```math
L_+ = \begin{pmatrix}
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \cdots \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \cdots \\
0 & \sqrt{2} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \cdots \\
0 & 0 & \sqrt{2} & 0 & 0 & 0 & 0 & 0 & 0 & \cdots \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \cdots \\
0 & 0 & 0 & 0 & 2 & 0 & 0 & 0 & 0 & \cdots \\
0 & 0 & 0 & 0 & 0 & \sqrt{6} & 0 & 0 & 0 & \cdots \\
0 & 0 & 0 & 0 & 0 & 0 & \sqrt{6} & 0 & 0 & \cdots \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 2 & 0 & \cdots \\
\vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \ddots
\end{pmatrix}.
```
Similarly, ``L_-`` has nonzero entries just above the diagonal.  Note
that these matrices are independent of ``s`` (though entries with ``ℓ
< |s|`` will never be used).

The ``R`` matrices are different: they are *all* diagonal, and they
all depend on ``s``.  ``R_z`` is simply ``s`` times the identity (up
to irrelevant factors in the ``ℓ < |s|`` entries):
```math
R_z = \begin{pmatrix}
s & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \cdots \\
0 & s & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \cdots \\
0 & 0 & s & 0 & 0 & 0 & 0 & 0 & 0 & \cdots \\
0 & 0 & 0 & s & 0 & 0 & 0 & 0 & 0 & \cdots \\
0 & 0 & 0 & 0 & s & 0 & 0 & 0 & 0 & \cdots \\
0 & 0 & 0 & 0 & 0 & s & 0 & 0 & 0 & \cdots \\
0 & 0 & 0 & 0 & 0 & 0 & s & 0 & 0 & \cdots \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & s & 0 & \cdots \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & s & \cdots \\
\vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \ddots
\end{pmatrix}.
```
But now, ``R_+`` is also diagonal.  For ``s=1``, for example:
```math
R_+ = \begin{pmatrix}
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \cdots \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \cdots \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \cdots \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \cdots \\
0 & 0 & 0 & 0 & 2 & 0 & 0 & 0 & 0 & \cdots \\
0 & 0 & 0 & 0 & 0 & 2 & 0 & 0 & 0 & \cdots \\
0 & 0 & 0 & 0 & 0 & 0 & 2 & 0 & 0 & \cdots \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 2 & 0 & \cdots \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 2 & \cdots \\
\vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \ddots
\end{pmatrix}.
```
The contrast with the sub-diagonal ``L_+`` is striking.  The
difference is that the matrix-vector product ``L_+ f`` represents
modes of a function with the *same* spin weight as ``f``, while ``R_+
f`` implicitly represents modes of a function with *different* spin
weight — a fact that is not readily apparent in matrix-vector
notation.

Because of this subtlety, we *cannot* represent ``R_x = (R_+ +
R_-)/2`` and ``R_y = (R_+ - R_-)/(2i)`` as matrices in the same way.
These operators result in functions with *indefinite* spin weight, so
they cannot be represented as a single vector of mode weights.  This
package defines functions for [`L₊`](@ref), [`L₋`](@ref),
[`R₊`](@ref), and [`R₋`](@ref), but only [`Lx`](@ref) and [`Ly`](@ref)
— not `Rx` or `Ry`.
