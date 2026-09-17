# Normalization

For any fixed values of ``ℓ`` and ``s``, the spin-weighted spherical
harmonics normalized as usual, satisfy the relation
```math
\sum_{m} |{}_{s}Y_{ℓ,m}(𝐐)|² = \frac{2ℓ+1}{4π},
```
for every ``𝐐 ∈ \mathrm{Spin}(3)``.  This result is crucial for
relating values of the spin-weighted spherical harmonics to Wigner
``𝔇`` matrices, and we will derive it here.

---

Begin by fixing ``ℓ`` and ``s``, and define the space of functions
with these eigenvalues:
```math
ℋ_{ℓ,s} = \left\{
  f : \mathrm{Spin}(3) → ℂ
  \middle|
  L² f = ℓ(ℓ+1) f
  \mathrm{\ \ and\ \ }
  R_z f = s f
\right\}.
```
Now, for a given point ``𝐏 ∈ \mathrm{Spin}(3)``, we define the
["reproducing
kernel"](https://en.wikipedia.org/wiki/Reproducing_kernel_Hilbert_space)
``K_𝐏`` as a function in ``ℋ_{ℓ,s}`` such that for *every*
square-integrable function ``f ∈ ℋ_{ℓ,s}``,
```math
f(𝐏) = \int_{\mathrm{Spin}(3)} K̄_𝐏(𝐐)\, f(𝐐)\, d𝐐.
```
Now, we want to expand this kernel in terms of the basis functions
``{}_{s}Y_{ℓ,m}`` of this space.  Recall that these are not
ortho*normal* over ``\mathrm{Spin}(3)``; they are only orthogonal.
They are normalized over ``𝕊²`` in the restricted sense discussed
[here](@ref sYlm_and_Dlmpm) so that, when integrating over
``\mathrm{Spin}(3)``, we get an extra factor of ``π/2``:
```math
\int_{\mathrm{Spin}(3)} {}_{s}Ȳ_{ℓ,m'}(𝐐)\, {}_{s}Y_{ℓ,m}(𝐐)\, d𝐐
= \frac{π}{2} δ_{m',m}.
```
Now, if we expand ``K_𝐏`` and ``f`` in terms of these basis
functions, we can calculate
```math
\begin{aligned}
f(𝐏)
&= \int_{\mathrm{Spin}(3)} K̄_𝐏(𝐐)\, f(𝐐)\, d𝐐 \\
&= \int_{\mathrm{Spin}(3)}
    \sum_{m',m} K̄_{𝐏,m'}\, {}_{s}Ȳ_{ℓ,m'}(𝐐)\,
    f_{m}\, {}_{s}Y_{ℓ,m}(𝐐)\, d𝐐 \\
&= \frac{π}{2} \sum_{m} K̄_{𝐏,m}\, f_{m},
\end{aligned}
```
the last of which implies that
```math
K_{𝐏,m} = \frac{2}{π} {}_{s}Ȳ_{ℓ,m}(𝐏)
```
for every ``m``.  That is,
```math
K_𝐏(𝐐) = \sum_m \frac{2}{π} {}_{s}Ȳ_{ℓ,m}(𝐏)\; {}_{s}Y_{ℓ,m}(𝐐).
```
Now we take the norm of the kernel function:
```math
\begin{aligned}
\|K_𝐏\|²_{\mathrm{Spin}(3)}
&= \int_{\mathrm{Spin}(3)} |K_𝐏(𝐐)|²\, d𝐐 \\
&= \int_{\mathrm{Spin}(3)}
    \sum_{m',m} \frac{4}{π²}
    {}_{s}Y_{ℓ,m'}(𝐏)\, {}_{s}Ȳ_{ℓ,m'}(𝐐)\, {}_{s}Ȳ_{ℓ,m}(𝐏)\, {}_{s}Y_{ℓ,m}(𝐐)\, d𝐐 \\
&= \frac{2}{π} \sum_{m} |{}_{s}Y_{ℓ,m}(𝐏)|².
\end{aligned}
```
Now, we are integrating with respect to a [Haar
measure](https://en.wikipedia.org/wiki/Haar_measure) on the group
``\mathrm{Spin}(3)``, which means that the integral must be invariant
under group actions.  In particular, this means that the norm of the
kernel function cannot depend on the choice of ``𝐏``, which means
that the sum in the last expression is independent of ``𝐏``:
```math
\sum_{m} |{}_{s}Y_{ℓ,m}(𝐏)|² = \sum_{m} |{}_{s}Y_{ℓ,m}(𝐏')|²
```
for *all* ``𝐏, 𝐏' ∈ \mathrm{Spin}(3)``.

The last trick is to just integrate this expression over
``\mathrm{Spin}(3)`` again, and evaluate it in two different ways.  On
one hand, we have
```math
\int_{\mathrm{Spin}(3)} \sum_{m} |{}_{s}Y_{ℓ,m}(𝐐)|²\, d𝐐
= \sum_{m} \int_{\mathrm{Spin}(3)} |{}_{s}Y_{ℓ,m}(𝐐)|²\, d𝐐
= \sum_{m} \frac{π}{2}
= (2ℓ+1) \frac{π}{2},
```
since there are ``2ℓ+1`` values of ``m`` in the sum.  On the other
hand, since the sum is independent of ``𝐐``, we can change to
integrating over a dummy variable and pull the sum out of the
integral:
```math
\int_{\mathrm{Spin}(3)} \sum_{m} |{}_{s}Y_{ℓ,m}(𝐐)|²\, d𝐐
= \sum_{m} |{}_{s}Y_{ℓ,m}(1)|²\, \int_{\mathrm{Spin}(3)} d𝐐'
= 2π² \sum_{m} |{}_{s}Y_{ℓ,m}(1)|².
```
Equating these two expressions, we find that
```math
\sum_{m} |{}_{s}Y_{ℓ,m}(𝐐)|² = \frac{2ℓ+1}{4π},
```
for arbitrary ``ℓ``, ``s``, and ``𝐐 ∈ \mathrm{Spin}(3)``.

Now, we use this to relate the spin-weighted spherical harmonics to
the Wigner ``𝔇`` matrices.  Recall that the SWSHs are [defined in
terms of their eigenvalues](@ref "``{}_{s}Y_{ℓ,m}`` as eigenfunctions"):
``{}_{s}Y_{ℓ,m}`` has ``L_z`` eigenvalue ``m`` and ``R_z`` eigenvalue
``s``.  The ``𝔇`` matrices [have eigenvalues](@ref "Defining ``𝔇^{(ℓ)}_{m', m}``")
``L_z 𝔇^{(ℓ)}_{m',m} = -m'\, 𝔇^{(ℓ)}_{m',m}`` and ``R_z 𝔇^{(ℓ)}_{m',m}
= m\, 𝔇^{(ℓ)}_{m',m}``, so it is the *conjugate*
``\overline{𝔇^{(ℓ)}_{m, -s}}`` that has ``L_z`` eigenvalue ``m`` and
``R_z`` eigenvalue ``s``.  Both functions therefore lie in the same
one-dimensional eigenspace, so
```math
{}_{s}Y_{ℓ,m}(𝐐) = κ_{ℓ,s}\, \overline{𝔇^{(ℓ)}_{m, -s}(𝐐)}
```
for some constant ``κ_{ℓ,s}`` independent of ``𝐐`` and of ``m`` — the
``L_±`` ladder operators are linear, so they commute with the map
between ``Y`` and ``\bar{𝔇}`` and force the same constant for every
``m``.

To find its magnitude, use the unitarity of ``𝔇^{(ℓ)}``, which makes
each of its columns a unit vector:
```math
\sum_{m} |{}_{s}Y_{ℓ,m}(𝐐)|²
= \frac{2ℓ+1}{4π}
= |κ_{ℓ,s}|² \sum_{m} |𝔇^{(ℓ)}_{m, -s}(𝐐)|²
= |κ_{ℓ,s}|².
```
Thus, up to a phase factor,
```math
|κ_{ℓ,s}| = \sqrt{\frac{2ℓ+1}{4π}}.
```
At ``s = 0`` that phase factor is 1, because of the near-universal
agreement that ``Y_{ℓ,0}(θ, ϕ)`` should be real-valued, which we
interpret to mean that ``{}_{0}Y_{ℓ,0}(𝟏)`` is real-valued; and
``𝔇^{(ℓ)}_{m', m}(𝟏) = δ_{m',m}`` is real, so the constant must be
real and positive.

The dependence on ``s`` is then fixed by the ``R_±`` ladder
operators.  Conjugation reverses them — a short calculation from the
definitions gives ``R_+ \bar{g} = -\overline{R_- g}`` — so
```math
R_+ \left\{\overline{𝔇^{(ℓ)}_{m, -s}}\right\}
= -\overline{R_- 𝔇^{(ℓ)}_{m, -s}}
= -\sqrt{(ℓ-s)(ℓ+s+1)}\; \overline{𝔇^{(ℓ)}_{m, -(s+1)}},
```
while the convention for the SWSHs is ``R_+ \left\{{}_{s}Y_{ℓ,m}\right\}
= \sqrt{(ℓ-s)(ℓ+s+1)}\, {}_{s+1}Y_{ℓ,m}`` with a *positive*
coefficient.  Comparing the two gives ``κ_{ℓ,s+1} = -κ_{ℓ,s}``, and
so, with the value at ``s = 0`` already fixed, this completely
determines the relationship between the spin-weighted spherical
harmonics and the Wigner ``𝔇`` matrices:
```math
{}_{s}Y_{ℓ,m}(𝐐)
= (-1)^s \sqrt{\frac{2ℓ+1}{4π}}\; \overline{𝔇^{(ℓ)}_{m, -s}(𝐐)}
```
for every ``ℓ``, ``m``, ``s``, and ``𝐐``.  This is the definition
recorded in the [conventions summary](@ref summary_swsh).  (For
half-integer ``s`` the factor ``(-1)^s`` means the principal branch
``e^{iπs}``.)
