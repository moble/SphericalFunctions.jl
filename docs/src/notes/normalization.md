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
