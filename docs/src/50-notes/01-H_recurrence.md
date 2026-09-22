# Algorithm for computing ``H``

The ``H`` array, as given by [Gumerov_2015](@citet), is related to Wigner's (small) ``d`` matrices —
which is itself related to the (big) ``𝔇`` matrices and the various spin-weighted
spherical harmonics ``{}_{s}Y_{ℓ,m}`` — via

```math
d^{(ℓ)}_{m',m} = ϵ_{m'} ϵ_{-m} H^{ℓ}_{m',m},
```

where

```math
ϵ_k =
  \begin{cases}
    1 & k\leq 0, \\
    (-1)^{\lfloor k \rfloor} & k > 0.
  \end{cases}
```

The floor matters only for half-integer indices (see below); for integer ``k`` this is the
familiar ``(-1)^k``.  Writing ``(-1)^k`` for half-integer ``k`` would be wrong, because
under this package's principal branch ``(-1)^k ≡ e^{iπk}`` is ``\pm i`` there, whereas
``ϵ`` is always a real sign.

(Note that I have swapped superscripts and subscripts on ``H``
compared to Gumerov and Duraiswami's paper, to be consistent with the
rest of this documentation.)

``H`` has various advantages over ``d`` and ``𝔇``,
including the fact that it can be efficiently and robustly calculated
via recurrence relations, and the following symmetry relations:

```math
\begin{aligned}
  H_{m', m}^ℓ(β) &= σ\, H_{m, m'}^ℓ(β) \\
  H_{m', m}^ℓ(β) &= σ\, H_{-m', -m}^ℓ(β) \\
  H_{m', m}^ℓ(β) &= (-1)^{ℓ-m}\, ϵ_{|m'|}\, H_{-m', m}^ℓ(π - β) \\
  H_{m', m}^ℓ(β) &= (-1)^{m'-m} H_{m', m}^ℓ(-β)
\end{aligned}
```

where

```math
σ =
  \begin{cases}
    1 & \text{integer } ℓ, m', m, \\
    \mathrm{sgn}(m)\, \mathrm{sgn}(m') & \text{half-integer } ℓ, m', m,
  \end{cases}
\qquad \mathrm{sgn}(0) = +1.
```

Both cases are the single expression ``σ = (-1)^{m'-m}\, ϵ_{|m'|}\, ϵ_{|m|}``, which is what
the familiar relation ``d^{(ℓ)}_{m',m} = (-1)^{m'-m} d^{(ℓ)}_{m,m'}``
becomes when it is written in terms of ``H``.  For integer indices
``ϵ_{|m'|} = (-1)^{m'}`` and ``(-1)^{ℓ-m} = (-1)^{ℓ+m}``, so the first
two relations reduce to plain symmetry and the last two to the forms
usually quoted, ``(-1)^{ℓ+m+m'}`` and ``(-1)^{m+m'}``.  For
half-integer indices ``σ`` is genuinely ``-1`` whenever
``\mathrm{sgn}(m) ≠ \mathrm{sgn}(m')`` — for example
``H^{1/2}_{1/2,-1/2} = -\sin(β/2)`` while ``H^{1/2}_{-1/2,1/2} =
+\sin(β/2)`` — and ``ℓ+m+m'`` is not even an integer, so the exponents
above have to be written in the forms ``ℓ-m`` and ``m'-m``, which are.
No other choice of the signs ``ϵ`` could remove ``σ``: on the
anti-diagonal ``m = -m'`` the ``ϵ`` factors appear squared, so that
``H`` and ``d`` coincide there, and ``d^{(1/2)}_{1/2,-1/2} =
-d^{(1/2)}_{-1/2,1/2}``.  Composing the first two relations shows that
``H_{m', m}^ℓ = H_{-m, -m'}^ℓ`` holds without any sign, for both kinds
of index.

!!! warning "Do not hand-roll these symmetries"
    In the code, ``σ`` is applied for you by [`wedge_value`](@ref
    SphericalFunctions.wedge_value) and [`wedge_source`](@ref
    SphericalFunctions.wedge_source), which are the only places the
    symmetries are encoded.  Reading the stored wedge of a
    [`HCalculator`](@ref) directly and transposing it by hand gives
    the wrong sign for every half-integer element with
    ``\mathrm{sgn}(m) ≠ \mathrm{sgn}(m')``.

Because of these symmetries — specifically the first two — we only
need to evaluate about 1/4 of all the elements for a given value of
``β``.


## Steps to compute ``H``

The following describes various details that are not entirely spelled
out by [Gumerov_2015](@citet).  All equation numbers refer to that
paper unless otherwise noted.

Because of the symmetries noted above, we only compute ``H_{m', m}^ℓ``
with ``m ≥ |m'|`` — roughly one quarter of all possible values.
Furthermore, for spin-weighted spherical harmonics of weight ``s``, we
only need to compute values with ``|m'| ≤ |s|``, which constitutes a
dramatic savings when ``|s| ≪ ℓₘₐₓ``.  The limit on ``|m'|`` is the
`m′ₘₐₓ` argument of [`HCalculator`](@ref).

The results for a given ``ℓ`` are held in two kinds of storage.  The
wedge ``m ≥ |m'|``, ``|m'| ≤ m'_{\mathrm{max}}`` is stored in an
[`HWedge`](@ref SphericalFunctions.HWedge),
row by row, and is overwritten when the next value of ``ℓ`` is
computed.  The ``m'=0`` axis ``H^{n}_{0,m}`` for ``m = 0, \ldots, n``
is stored separately, in two [`HAxis`](@ref SphericalFunctions.HAxis)
buffers that hold successive orders ``n`` and ``n+1``.  Each advance
of step 2 below overwrites the older of the two with the next order,
and the two then exchange roles.  Two buffers are needed because step
3 reads the axis of order ``ℓ+1`` while building the wedge of order
``ℓ``.

Gumerov and Duraiswami consider only integer ``ℓ``, but the same
algorithm works for half-integer ``ℓ, m', m`` with one substantive
change.  As explained [below](@ref "Why the recurrences hold for
half-integer indices"), the relation behind steps 4 and 5 makes no
assumption that the indices are integers; only the seed from which
those steps start does.  For half-integer ``ℓ``, the ``m'=0`` axis of
steps 1 and 2 is run at the *integer* order ``j = ℓ - 1/2``, and step
3 — which climbs from the ``m'=0`` row, which does not exist for
half-integers, to the ``m'=1`` row — is replaced by a seed that
produces both rows ``m' = \pm 1/2`` at once from that axis.  Steps 4
and 5 then proceed from those rows.  It is convenient to write
``ℓ_\mathrm{min}`` for the smallest allowed value of ``ℓ``, which is
``0`` for integer indices and ``1/2`` for half-integer ones; the loop
bounds below are written in terms of it, and specialize to the integer
bounds when ``ℓ_\mathrm{min} = 0``.


### Why the recurrences hold for half-integer indices

Gumerov and Duraiswami construct their recurrences from the regular
solutions ``j_n(kr)\, Y_n^m(θ, φ)`` of the Helmholtz equation with
wavenumber ``k``.  A rotated solution of degree ``n`` is a combination
of unrotated solutions of the same degree, with coefficients given by
the rotation matrix.  Applying the operator ``\tfrac{1}{k}∇`` to both
sides of that relation, and using the fact that the gradient of a
rotated function is the rotated gradient of the function, gives
relations between rotation coefficients of *different* degrees,
because each component of ``\tfrac{1}{k}∇`` maps a solution of degree
``n`` to a combination of solutions of degrees ``n-1`` and ``n+1``.
Relation (41), used in step 3, is of this kind: it computes ``H^{ℓ}``
from ``H^{ℓ+1}``.

That construction cannot say anything about half-integer indices,
because it rests on scalar functions in three-dimensional space, and
those exist only for integer degree — there is no solution ``j_n
Y_n^m`` with half-integer ``n`` and ``m``.  Relation (50), on the
other hand, which is the only relation used in steps 4 and 5, involves
just one value of ``ℓ``, and has a much more elementary proof that
uses nothing but the representation theory of rotations.  In the basis
``|ℓ, m⟩`` of an irreducible representation, the ``d`` matrix is the
matrix exponential
```math
d^{(ℓ)}(β) = e^{-iβ J_y},
\qquad
J_y = \frac{J_+ - J_-}{2i},
```
where the ladder operators act as
```math
J_+ |ℓ, m⟩ = λ^{m}_{ℓ}\, |ℓ, m+1⟩,
\qquad
J_- |ℓ, m⟩ = λ^{m-1}_{ℓ}\, |ℓ, m-1⟩,
\qquad
λ^{m}_{ℓ} = \sqrt{(ℓ-m)(ℓ+m+1)}.
```
A matrix commutes with its own exponential, so ``J_y\, d^{(ℓ)}(β) =
d^{(ℓ)}(β)\, J_y``.  The ``(m', m)`` element of the left-hand side is
``\left[λ^{m'-1}_{ℓ} d^{(ℓ)}_{m'-1,m} - λ^{m'}_{ℓ}
d^{(ℓ)}_{m'+1,m}\right]/2i``, and that of the right-hand side is
``\left[λ^{m}_{ℓ} d^{(ℓ)}_{m',m+1} - λ^{m-1}_{ℓ}
d^{(ℓ)}_{m',m-1}\right]/2i``, so that
```math
λ^{m'}_{ℓ}\, d^{(ℓ)}_{m'+1, m}
  = λ^{m'-1}_{ℓ}\, d^{(ℓ)}_{m'-1, m}
  + λ^{m-1}_{ℓ}\, d^{(ℓ)}_{m', m-1}
  - λ^{m}_{ℓ}\, d^{(ℓ)}_{m', m+1}.
```
This is relation (50), written in terms of ``d`` rather than ``H``;
the conversion between the two forms is described in step 4 below.
Nothing in the argument depends on whether ``ℓ`` is an integer.  The
ladder operators have the same matrix elements in every irreducible
representation, integer or half-integer, and ``d^{(ℓ)}(β) = e^{-iβ
J_y}`` holds in every one of them.  (For half-integer ``ℓ`` it is this
exponential that makes ``d^{(ℓ)}(β + 2π) = -d^{(ℓ)}(β)``.)  The
argument does rely on the standard phase convention, in which the
matrix elements of ``J_\pm`` are real and positive; this package uses
that convention.  It does not rely on the sign of ``β`` in the
exponent, since replacing ``e^{-iβ J_y}`` by ``e^{+iβ J_y}`` gives a
matrix that commutes with ``J_y`` just as well.

The recurrences that change ``ℓ`` are the ones that actually depend on
the kind of index.  In the language of representation theory,
``\tfrac{1}{k}∇`` is a vector operator, and Gumerov and Duraiswami's
relations between degrees ``n`` and ``n \pm 1`` are instances of the
Clebsch–Gordan series for the product of the spin-``1`` representation
with the spin-``n`` representation.  Coupling with spin ``1`` can
never lead from an integer representation to a half-integer one, but
coupling with spin ``1/2`` can, and that is what Eqs. 4.8.2(14) and
(15) of [Varshalovich_1988](@citet) express.  The half-integer seed
described below therefore takes the place of step 3, with coupling to
spin ``1/2`` in place of the gradient.

There is also a sense in which both kinds of index are solutions of
one higher-dimensional problem.  The functions ``𝔇^{(ℓ)}_{m',m}``,
for all integer and half-integer ``ℓ`` together, are the matrix
elements of the irreducible representations of the group of unit
quaternions, which is the three-sphere.  By the Peter–Weyl theorem
they form a complete orthogonal basis for functions on that sphere,
and each ``𝔇^{(ℓ)}_{m',m}`` is the restriction to it of a harmonic
polynomial of degree ``2ℓ`` in the four components of the quaternion —
a four-dimensional hyperspherical harmonic.  The integer-``ℓ``
functions are those that are even under ``R \to -R``, and those with
``m' = 0`` or ``m = 0`` reduce to ordinary spherical harmonics on the
two-sphere.  An argument like Gumerov and Duraiswami's could
presumably be made in that four-dimensional setting, and would cover
both kinds of index at once, but the commutation argument above makes
it unnecessary.


### Step 1: Initialize ``H^{0}_{0,0}``

Set ``H^{0}_{0,0}=1``.  This is also the starting point of the integer
axis used for half-integer ``ℓ``.


### Step 2: ``H_{0,m}^{n-1} \to H_{0,m}^{n}`` for ``m \geq 0``

Compute values ``H^{n}_{0,m}(β)`` for ``m=0,\ldots,n`` from
``H^{n-1}_{0,m}(β)``, one order at a time.  Note that on this axis
``H`` and ``d`` coincide, because ``ϵ_0 ϵ_{-m} = 1`` for ``m ≥ 0``.
Using Eq. (32), we see that within Gumerov and Duraiswami's
conventions
```math
\begin{aligned}
  H^{n}_{0,m}(β) &= (-1)^m \sqrt{\frac{(n-|m|)!}{(n+|m|)!}} P^{|m|}_{n}(\cos β) \\
                 &= \frac{1}{\sqrt{k_m (2n+1)}} P̄_{n,|m|}(\cos β).
\end{aligned}
```
Here, ``k_0=1`` and ``k_m=2`` for ``m>0``, and ``P̄`` is defined as
```math
  P̄_{n,|m|} = \sqrt{\frac{k_m(2n+1)(n-m)!}{(n+m)!}} P_{n,|m|}.
```
Note that the factor of ``(-1)^m`` in the first equation above is
different from the convention used here, and is related to the
[Condon-Shortley
phase](https://en.wikipedia.org/wiki/Spherical_harmonics#Condon%E2%80%93Shortley_phase).
Note that Gumerov and Duraiswami use the notation ``P^{|m|}_{n}``,
whereas we are using the notation ``P_{n,|m|}`` — which usually differ
by a factor of ``(-1)^m``.

We use the "fully normalized" associated Legendre functions (fnALF)
``P̄`` because, as explained by [Xing_2019](@citet), it is possible to
compute these values very efficiently and accurately, while also
delaying the onset of overflow and underflow.

The algorithm Xing et al. describe as the best for computing ``P̄`` is
due to [Strakhov_1980](@citet) via [Belikov_1991](@citet), and is
given by them as
```math
\begin{aligned}
  P̄_{0,0} &= 1 \\
  P̄_{1,0} &= \sqrt{3} \cos β \\
  P̄_{1,1} &= \sqrt{3} \sin β \\
  P̄_{n,0} &= a_n \cos β P̄_{n-1,0} - b_n \frac{\sin β}{2} P̄_{n-1,1} \\
  P̄_{n,m} &=
    c_{n,m} \cos β P̄_{n-1,m}
    - \sin β \left[ d_{n,m} P̄_{n-1,m+1} - e_{n,m} P̄_{n-1,m-1} \right],
\end{aligned}
```
where the coefficients are given by
```math
\begin{aligned}
  a_n &= \sqrt{\frac{2n+1}{2n-1}} \\
  b_n &= \sqrt{\frac{2(n-1)(2n+1)}{n(2n-1)}} \\
  c_{n,m} &= \frac{1}{n} \sqrt{\frac{(n+m)(n-m)(2n+1)}{2n-1}} \\
  d_{n,m} &= \frac{1}{2n} \sqrt{\frac{(n-m)(n-m-1)(2n+1)}{2n-1}} \\
  e_{n,m} &= \frac{1}{2n} \sqrt{\frac{2}{2-δ_0^{m-1}}} \sqrt{\frac{(n+m)(n+m-1)(2n+1)}{2n-1}}.
\end{aligned}
```

Now, we can directly obtain a recurrence relation for ``H^{n}_{0,m} =
P̄_{n,|m|} / \sqrt{k_m (2n+1)} `` from those expressions:
```math
\begin{aligned}
  H^{0}_{0,0} &= 1 \\
  H^{1}_{0,0} &= \cos β \\
  H^{1}_{0,1} &= \sqrt{1/2} \sin β \\
  H^{n}_{0,0} &= \cos β H^{n-1}_{0,0} - b̄_n \sin β H^{n-1}_{0,1} \\
  H^{n}_{0,m} &=
    c̄_{n,m} \cos β H^{n-1}_{0,m}
    - \sin β \left[ d̄_{n,m} H^{n-1}_{0,m+1} - ē_{n,m} H^{n-1}_{0,m-1} \right],
\end{aligned}
```
where the coefficients are given by
```math
\begin{aligned}
  b̄_n &= \sqrt{\frac{n-1}{n}} \\
  c̄_{n,m} &= \frac{1}{n} \sqrt{(n+m)(n-m)} \\
  d̄_{n,m} &= \frac{1}{2n} \sqrt{(n-m)(n-m-1)} \\
  ē_{n,m} &= \frac{1}{2n} \sqrt{(n+m)(n+m-1)}.
\end{aligned}
```
Note that the coefficients all simplified (in fact, ``a_n``
disappeared), without any increase in the complexity of the recurrence
relations themselves.  Rewriting Belikov's algorithm explicitly in
terms of the ``H^{n}_{0,m}`` also allows us to avoid an extra
normalization step.

The terms involving ``H^{n-1}_{0,n}`` and ``H^{n-1}_{0,n+1}``, which
lie outside the axis of order ``n-1``, vanish: at ``m = n-1`` the
coefficient ``d̄_{n,m}`` is zero, and at ``m = n`` so is ``c̄_{n,m}``.
The code treats these two values of ``m`` separately, so that it never
reads storage beyond the end of the axis.

For integer ``ℓ``, the axis of order ``ℓ`` is copied into the ``m'=0``
row of the wedge, and the axis of order ``ℓ+1`` is used by step 3.
For half-integer ``ℓ``, only the axis of order ``j = ℓ - 1/2`` is
used, by the seed described below.  Either way, the axis is always
indexed by integer orders; half-integer values never appear in this
step.


### Step 3: ``H_{0,m}^{ℓ+1} \to H_{1,m}^{ℓ}`` for ``m \geq 1``

Compute ``H^{ℓ}_{1,m}(β)`` for ``m=1,\ldots,ℓ`` using relation (41).
Symmetry and shift of the indices allow this relation to be written as
```math
b^{0}_{ℓ+1} H^{ℓ}_{1, m}
  = \frac{b^{−m−1}_{ℓ+1} (1−\cos β)}{2} H^{ℓ+1}_{0, m+1}
  − \frac{b^{ m−1}_{ℓ+1} (1+\cos β)}{2} H^{ℓ+1}_{0, m−1}
  − a^{m}_{ℓ} \sin β H^{ℓ+1}_{0, m}.
```
Here the constants are defined by
```math
a^{m}_{ℓ} = \sqrt{\frac{(ℓ+m+1)(ℓ-m+1)} {(2ℓ+1)(2ℓ+3)}},
```
```math
b^{m}_{ℓ} = \mathrm{sgn}(m) \sqrt{\frac{(ℓ-m-1)(ℓ-m)} {(2ℓ-1)(2ℓ+1)}}.
```
Note that all values are assumed to be zero whenever ``|m| > ℓ``, we use
``\mathrm{sgn}(0)=1`` (unlike the common convention that
``\mathrm{sgn}(0)=0``), and we have ``a^{m}_{ℓ} = a^{-m}_{ℓ}``.  Also
note that these coefficients *only* appear in this step, and because
of how they appear (specifically, because ``b`` always appears with
argument ``ℓ+1``), we can factor out the denominators in the
definitions of the constants.  The signs can be resolved as well,
because for ``m ≥ 1`` we have ``\mathrm{sgn}(-m-1) = -1`` and
``\mathrm{sgn}(m-1) = \mathrm{sgn}(0) = +1``, so that all three terms
on the right-hand side enter with the same sign.  We obtain this
simplified formula
```math
H^{ℓ}_{1, m}
  = -\frac{1}{\sqrt{ℓ(ℓ+1)}} \left[
      \frac{\bar{b}^{−m−1}_{ℓ+1} (1−\cos β)}{2} H^{ℓ+1}_{0, m+1}
      + \frac{\bar{b}^{ m−1}_{ℓ+1} (1+\cos β)}{2} H^{ℓ+1}_{0, m−1}
      + \bar{a}^{m}_{ℓ} \sin β H^{ℓ+1}_{0, m}
    \right],
```
with
```math
\bar{a}^{m}_{ℓ} = \sqrt{(ℓ+m+1)(ℓ-m+1)},
```
```math
\bar{b}^{m}_{ℓ+1} = \sqrt{(ℓ-m)(ℓ-m+1)}.
```

This step is skipped when ``ℓ = 0`` or ``m'_{\mathrm{max}} = 0``,
since there is then no ``m'=1`` row to fill.


### Step 3 for half-integer ``ℓ``: ``H_{0,m}^{ℓ-1/2} \to H_{\pm 1/2,m}^{ℓ}``

For half-integer ``ℓ`` there is neither an ``m'=0`` row to copy the
axis into nor an ``m'=1`` row for step 3 to produce.  Instead, the two
rows ``m' = \pm 1/2`` are computed directly from the integer axis of
order ``j = ℓ - 1/2``, using Eqs. 4.8.2(14) and (15) of
[Varshalovich_1988](@citet), which relate an element with half-integer
``ℓ`` to elements with ``ℓ - 1/2``.  Applied at ``M' = ∓1/2``, where
the elements on the right-hand side fall on the ``m'=0`` axis, and
rewritten in this package's conventions and ``(m', m)`` order, they
give
```math
\begin{aligned}
  H^{ℓ}_{+1/2, m}
    &= \frac{1}{\sqrt{ℓ + 1/2}} \left[
      \sqrt{ℓ+m}\, \cos\tfrac{β}{2}\, H^{j}_{0, m-1/2}
      - \sqrt{ℓ-m}\, \sin\tfrac{β}{2}\, H^{j}_{0, m+1/2}
    \right], \\
  H^{ℓ}_{-1/2, m}
    &= \frac{1}{\sqrt{ℓ + 1/2}} \left[
      \sqrt{ℓ+m}\, \sin\tfrac{β}{2}\, H^{j}_{0, m-1/2}
      + \sqrt{ℓ-m}\, \cos\tfrac{β}{2}\, H^{j}_{0, m+1/2}
    \right],
\end{aligned}
```
for ``m = 1/2, \ldots, ℓ``.  At ``m = ℓ`` the second term vanishes
along with its coefficient; it would refer to ``H^{j}_{0, j+1}``,
which does not exist, so the code treats that value of ``m``
separately.  On these two rows ``H`` and ``d`` coincide, because
``ϵ_{\pm 1/2}\, ϵ_{-m} = 1`` for ``m > 0``.  The denominator is
``\sqrt{j+1} ≥ 1``, so nothing here is singular, and the cost is
``O(ℓ)`` per rotor — less than the integer step 3, which needs the
axis one order higher.  The [comparison with Varshalovich et al.](@ref
"Varshalovich et al. (1988)") checks Eqs. 4.8.2(14) and (15), as
transcribed from the book, against the closed form of their Eq.
4.3.1(2).

Both rows are required, even when only ``m' ≥ 0`` is wanted: the
corner element ``H^{ℓ}_{-1/2, 1/2}`` cannot be reached from the
``+1/2`` row by steps 4 and 5 without leaving the wedge.  This is why
``m'_{\mathrm{max}}`` can be no smaller than ``1/2`` for half-integer
indices.

Unlike every other step, the seed needs the half angles ``\cos(β/2)``
and ``\sin(β/2)`` rather than just ``e^{iβ}``, so the calculator
stores those for each rotor.  How they are obtained depends on the
form in which ``β`` was given.  From a rotor ``R = W + X𝐢 + Y𝐣 +
Z𝐤``, they are ``\sqrt{W^2+Z^2}/|R|`` and ``\sqrt{X^2+Y^2}/|R|``,
which are accurate near both poles; both are non-negative, so ``β ∈
[0, π]``, and the double-cover sign of the rotor enters only through
the phases of step 7.  From an angle ``β`` they are computed directly,
which respects the ``4π`` periodicity of half-integer ``d``.  A bare
phase ``e^{iβ}`` determines ``β`` only modulo ``2π``, and hence ``d``
only up to the sign ``(-1)^{2ℓ}``; the branch ``β ∈ (-π, π]`` is used
in that case.


### Step 4: ``H_{m',m-1}^{ℓ}, H_{m'-1,m}^{ℓ}, H_{m',m+1}^{ℓ} \to H_{m'+1,m}^{ℓ}`` for ``m' \geq 1`` and ``m > m'``

Recursively compute ``H^{ℓ}_{m'+1, m}(β)`` for ``m' = 1 -
ℓ_\mathrm{min}, \ldots, \min(ℓ, m'_{\mathrm{max}}) - 1`` and ``m =
m'+1, \ldots, ℓ`` using relation (50) resolved with respect to
``H^{ℓ}_{m'+1, m}``:
```math
d^{m'}_{ℓ} H^{ℓ}_{m'+1, m}
  = d^{m'−1}_{ℓ} H^{ℓ}_{m'−1, m}
  − d^{m−1}_{ℓ} H^{ℓ}_{m', m−1}
  + d^{m}_{ℓ} H^{ℓ}_{m', m+1}
```
(where the last term drops out for ``m=ℓ``).  The constants are
defined by
```math
d^{m}_{ℓ} = \frac{\mathrm{sgn}(m)}{2} \sqrt{(ℓ-m)(ℓ+m+1)}.
```
(These ``d^{m}_{ℓ}`` are Gumerov and Duraiswami's coefficients, and
are unrelated to both the Wigner ``d`` matrix and the coefficient
``d̄_{n,m}`` of step 2.)  The factor of ``1/2`` cancels throughout, so
the code drops it.  The lower limit on ``m`` is ``m'+1`` rather than
``m'``, because the element being computed, ``H^{ℓ}_{m'+1, m}``, lies
in the wedge only when ``m ≥ m'+1``.  For integer ``ℓ`` the loop
starts at ``m'=1``, reading the ``m'=0`` row copied from the axis and
the ``m'=1`` row from step 3; for half-integer ``ℓ`` it starts at
``m'=1/2``, reading the two rows of the seed.

The signs of the coefficients need some care.  On the ``m`` side they
are all ``+1``, because ``m ≥ m'+1 ≥ 3/2`` throughout.  On the ``m'``
side, for integer indices, both ``\mathrm{sgn}(m')`` and
``\mathrm{sgn}(m'-1)`` are ``+1`` as well, so this step is entirely
sign-free.  That is not true for half-integer indices: at ``m' = 1/2``
the coefficient of ``H^{ℓ}_{m'-1, m} = H^{ℓ}_{-1/2, m}`` is
``d^{-1/2}_{ℓ}``, and ``\mathrm{sgn}(-1/2) = -1``.

The form of relation (50) in terms of ``d``, derived [above](@ref "Why
the recurrences hold for half-integer indices"), involves no signs of
this kind, and holds for every ``ℓ``.  That it takes exactly Gumerov
and Duraiswami's form, with their ``m'``-side signs, when rewritten in
terms of ``H`` is a consequence of the definition of ``ϵ``:
substituting ``d^{(ℓ)}_{m',m} = ϵ_{m'} ϵ_{-m} H^{ℓ}_{m',m}``, noting
that ``λ^{m}_{ℓ} = 2\,|d^{m}_{ℓ}|``, and dividing through by ``ϵ_{m'}
ϵ_{-m}`` gives the relation above, because ``ϵ_{k+1}/ϵ_{k} =
-\mathrm{sgn}(k)`` for every ``k``.  The extension ``ϵ_k =
(-1)^{\lfloor k \rfloor}`` for half-integer ``k > 0`` is the only one,
given ``ϵ_k = 1`` for ``k ≤ 0``, for which that ratio still holds.  On
the ``m`` side the same calculation gives factors of
``-\mathrm{sgn}(-m)`` on the ``H^{ℓ}_{m', m−1}`` term and
``-\mathrm{sgn}(-m-1)`` on the ``H^{ℓ}_{m', m+1}`` term.  For integer
``m`` these equal ``\mathrm{sgn}(m-1)`` and ``\mathrm{sgn}(m)``, as
Gumerov and Duraiswami write them, but for half-integer ``m`` they
equal ``\mathrm{sgn}(m)`` and ``\mathrm{sgn}(m+1)``, which differ from
those at ``m = \pm 1/2``.  Neither step 4 nor step 5 ever visits ``m <
3/2``, so this difference has no effect here; it would matter only if
the recurrence were extended outside the wedge.


### Step 5: ``H_{m',m-1}^{ℓ}, H_{m'+1,m}^{ℓ}, H_{m',m+1}^{ℓ} \to H_{m'-1,m}^{ℓ}`` for ``m' \leq 0`` and ``m > -m'``

Recursively compute ``H^{ℓ}_{m'−1, m}(β)`` for ``m' = -ℓ_\mathrm{min},
\ldots, -\min(ℓ, m'_{\mathrm{max}}) + 1`` (in decreasing order) and
``m = 1-m', \ldots, ℓ`` using relation (50) resolved with respect to
``H^{ℓ}_{m'−1, m}``:
```math
d^{m'−1}_{ℓ} H^{ℓ}_{m'−1, m}
  = d^{m'}_{ℓ} H^{ℓ}_{m'+1, m}
  + d^{m−1}_{ℓ} H^{ℓ}_{m', m−1}
  − d^{m}_{ℓ} H^{ℓ}_{m', m+1}
```
(where the last term drops out for ``m=ℓ``).  As in step 4, the signs
on the ``m`` side are all ``+1``, because ``m ≥ 1-m' ≥ 1``.  On the
``m'`` side they are not: ``\mathrm{sgn}(m'-1) = -1`` always, and
``\mathrm{sgn}(m') = -1`` except at ``m' = 0``.

Although Gumerov and Duraiswami specify the loop over ``m'`` to start
at ``-1``, that would require the row ``m' = -1`` to exist before the
loop begins, and nothing earlier produces it.  The loop here starts
instead at ``m' = -ℓ_\mathrm{min}`` — that is, at ``m' = 0`` for
integer ``ℓ`` and at ``m' = -1/2`` for half-integer ``ℓ`` — so that
its first pass produces the row ``m' = -1`` (or ``-3/2``) from rows
that steps 2 and 3 (or the seed) have already filled.  Also, the lower
limit on ``m`` is ``1-m'``, the smallest value for which
``H^{ℓ}_{m'−1, m}`` lies in the wedge.  An earlier version of these
notes started the loop over ``m`` at ``-m'``, which computed one
element outside the wedge at each ``m'``, and thereby required the
elements ``H^{n}_{0, -1}`` to be set in advance; with the limits given
here, nothing outside the wedge is ever read or written.


### Step 6: Use symmetries to fill in the rest of ``H``

Every element of ``H^{ℓ}`` with ``|m'| ≤ m'_{\mathrm{max}}`` or ``|m|
≤ m'_{\mathrm{max}}`` is equal, up to the sign ``σ``, to an element of
the wedge computed above.  Specifically, [`wedge_source`](@ref
SphericalFunctions.wedge_source) maps a requested element ``H^{ℓ}_{m',
m}`` to a stored one as follows:
```math
H^{ℓ}_{m', m} =
  \begin{cases}
    H^{ℓ}_{m', m} & |m'| ≤ m'_{\mathrm{max}} \text{ and } m ≥ |m'|, \\
    σ\, H^{ℓ}_{-m', -m} & |m'| ≤ m'_{\mathrm{max}} \text{ and } -m ≥ |m'|, \\
    σ\, H^{ℓ}_{m, m'} & |m| ≤ m'_{\mathrm{max}} \text{ and } m' ≥ |m|, \\
    H^{ℓ}_{-m, -m'} & |m| ≤ m'_{\mathrm{max}} \text{ and } -m' ≥ |m|.
  \end{cases}
```
If both ``|m'|`` and ``|m|`` exceed ``m'_{\mathrm{max}}``, the element
cannot be obtained, and an error is thrown.

The batched engine never actually runs this step.  It leaves the wedge
alone, and applies the symmetries on the fly, as each element is read,
through [`wedge_value`](@ref SphericalFunctions.wedge_value) or
`wedge_source`.  Only the unbatched, integer-only reference
implementation, [`recurrence_step6!`](@ref
SphericalFunctions.recurrence_step6!), fills in the rest of the matrix
explicitly.


### Step 7: Include phases to obtain ``d`` or ``𝔇``

The ``d`` matrix is obtained from ``H`` by the signs of the defining
relation, ``d^{(ℓ)}_{m',m} = ϵ_{m'} ϵ_{-m} H^{ℓ}_{m',m}``, together
with ``σ`` when the element of ``H`` was read by symmetry.  The ``𝔇``
matrix then follows from the Euler angles as
```math
𝔇^{(ℓ)}_{m',m}(R) = e^{-im'α}\, d^{(ℓ)}_{m',m}(β)\, e^{-imγ}.
```
For half-integer ``m'`` and ``m`` the phases are not integer powers of
``e^{iα}`` and ``e^{iγ}``, but no square roots of those phases are
needed, because the sum and difference ``m' \pm m`` are always
integers.  With
```math
z_{+} = e^{i(α+γ)/2},
\qquad
z_{-} = e^{i(α-γ)/2},
```
which are computed directly from the components of the rotor, we have
```math
e^{-i(m'α + mγ)} = \overline{z_{+}^{\,m'+m}\, z_{-}^{\,m'-m}}
```
with integer exponents, for integer and half-integer indices alike.
For half-integer indices one of the two exponents is odd, so that
``𝔇(-R) = -𝔇(R)`` follows automatically, with no choice of branch.
At the poles one of ``z_\pm`` is undefined and is set to ``1``; the
elements of ``d`` it would multiply vanish there, so the choice has no
effect.


## Batching over rotors

Every index this recursion touches is sequential for a single
rotation: each ``ℓ`` is built from the one before, and within a given
``ℓ`` the ladders of steps 4 and 5 read the rows that the previous
iteration has just written.  The rotation itself is the one free
dimension — the same arithmetic on unrelated data — which is why the
implementation runs a whole batch of rotors through every step at
once, with the rotor index varying fastest in memory so that the
innermost loop is the one that can be vectorized.  What that means for
the interface is described in the section on [iterating over ``ℓ`` and
reusing the storage](@ref "Iterating over ``ℓ`` and reusing the
storage").


## Pre-computing constants versus computing on the fly

Each of the constants in steps 2 through 5 involves a square root, and
most involve a division, which can be very costly to compute.  It can
be advantageous to pre-compute the constants, and simply index the
pre-computed arrays rather than re-computing them on each recursion.

Measurements on the earlier, whole-array implementation of this
recursion found that, *if* we include the cost of computing all these
constants in a single call to the ``H`` recurrence, it can be much
cheaper to compute each constant as needed within the algorithm,
rather than computing them all at once at the beginning of the
algorithm — but only for very small computations, such as those
involving ``ℓ_{\mathrm{max}} ≈ 10``.  Beyond this, despite the storage
penalties for all those constants, it turned out to be better to
pre-compute them.  However, it should be noted that the fractional
cost of storing the constants is ``\sim 3/ℓ_{\mathrm{max}}`` compared
to just storing ``H`` itself, so this will never be a very significant
amount of space.  On the other hand, if we can pre-compute the
constants just once, and store them between multiple calls to the
``H`` recurrence, then it was always advantageous to do so — typically
by factors of 2 or 3 in speed.

The current implementation nonetheless computes every constant on the
fly, because batching changes the balance.  Each constant is computed
once for a given ``(ℓ, m', m)`` and then used for every rotor in the
batch, so its cost is divided among all of them.  Every transform in
the package uses the batched path, and there the constants are a small
part of the total.  A single rotor, on the other hand, pays the full
cost of every constant: per element and per rotor, a full sweep with
one rotor was measured to be between 7 and 19 times slower than the
same sweep with a batch of 512, and most of that difference is the
work of computing constants, which a cache could recover.  Such a
cache has not been added, because it would give the calculators —
which currently allocate nothing after construction and have no state
to invalidate — both ``O(ℓ_{\mathrm{max}})`` storage and invalidation
logic.  It could be added later without changing the interface, since
it would be entirely internal to [`recurrence!`](@ref).
