# Sampling theorems and transformations of spin-weighted spherical harmonics

```@meta
CurrentModule = SphericalFunctions
```

[McEwenWiaux_2011](@citet) (MW) provide a very thorough review of the
literature on sampling theorems related to spin-weighted spherical
harmonics up to 2011.  [Reinecke_2013](@citet) (RS) outlined one of
the more efficient and accurate implementations of spin-weighted
spherical harmonic transforms (``s``SHT) currently available as
`libsharp`, but their algorithm is ``∼4L²``, whereas McEwen and
Wiaux's is ``∼2L²``, while [Elahi_2018](@citet) (EKKM) have obtained
the optimal result that scales as ``∼L²``.

The downside of the EKKM algorithm is that the ``θ`` values at which
to sample have to be obtained by iteratively minimizing the condition
numbers of various matrices (which are involved in the computation
itself).  This expensive step only has to be performed once per choice
of spin ``s`` and maximum ``ℓ`` value ``L``.  Otherwise, the algorithm
is accurate at small ``L``, but its sample points become badly
conditioned as ``L`` grows: in double precision, a round trip loses
about 3.5 digits by ``L = 32``, 6 by ``L = 48`` and 8 by ``L = 64`` for
``s = 0``, and about 9 by ``L = 32`` for ``s = 2``.  This does not
compare favorably with the MW algorithm, which has slowly growing
errors through ``L = 4096``.

## EKKM analysis

The EKKM analysis looks like the following (with some notational
changes).  We begin by defining
```math
  {}_{s}\tilde{f}_{θ}(m) := \int_0^{2π} {}_sf(θ, ϕ)\, e^{-imϕ}\, dϕ.
```
We will denote the vector of these quantities for all values of
``θ`` as ``{}_{s}\tilde{𝐟}_m``.  Inserting the
``{}_sY_{ℓ,m}`` expansion for ``{}_sf(θ, ϕ)``, and
performing the integration using orthogonality of complex
exponentials, we can find that
```math
  {}_{s}\tilde{f}_{θ}(m) = (-1)^s\, 2π \sum_{ℓ=\Delta}^L \sqrt{\frac{2ℓ+1}{4π}}\, d_{m,-s}^{ℓ}(θ)\, {}_sf_{ℓ,m},
```
where ``\Delta = \max(|m|, |s|)`` is the smallest ``ℓ`` that has a
mode with this ``m`` and spin weight ``s``.
Now, denoting the vector of ``{}_sf_{ℓ,m}`` for all values of
``ℓ`` as ``{}_s𝐟_m``, we can write this as a matrix-vector
equation:
```math
  {}_{s}\tilde{𝐟}_m = (-1)^s\, 2π\, {}_s𝐝_{m}\, {}_s𝐟_m.
```
We are effectively measuring the ``{}_{s}\tilde{𝐟}_m``
values, we can easily construct the ``{}_s𝐝_{m}`` matrix, and
we are seeking the ``{}_s𝐟_m`` values, so we can just invert
this equation to solve for the latter.


## Discretizing the Fourier transform

Now, the only flaw in this analysis is that we have undersampled
everywhere except ``ℓ = L``, which means that the second equation
(re-expressing the Fourier transforms as a sum using orthogonality of
complex exponentials) isn't quite right; in general there is some
folding due to aliasing of higher-frequency modes, so we need an
additional sum over ``|m'|>|m|``.  Or perhaps more precisely, the
first equation isn't actually what we implement.  It should look more
like this:
```math
  {}_{s}\tilde{f}_{j}(m) := \sum_{k=0}^{2j} {}_sf(θ_j, ϕ_k)\, e^{-imϕ_k}\, \Delta ϕ,
```
where ``ϕ_k = \frac{2π k}{2j+1}``, and ``\Delta ϕ =
\frac{2π}{2j+1}``.  (Recall the subtle notational distinction common
in time-frequency analysis that ``\tilde{s}(t_j) = \Delta t
\tilde{s}_j``, which would suggest we use ``{}_{s}\tilde{f}_{j}(m) =
\Delta ϕ\, {}_{s}\tilde{f}_{j,m}``.)  Next, we can insert the
expansion for ``{}_sf(θ, ϕ)``:

```math
\begin{aligned}
    {}_{s}\tilde{f}_{j}(m)
    &= \sum_{k=0}^{2j} \sum_{ℓ,m'} {}_sf_{ℓ,m'}\, {}_sY_{ℓ,m'}(θ_j, ϕ_k)\, e^{-imϕ_k}\, \Delta ϕ \\
    &= \sum_{k=0}^{2j} \sum_{ℓ,m'} {}_sf_{ℓ,m'}\, (-1)^{s}\, \sqrt{\frac{2ℓ+1}{4π}}\, d_{ℓ}^{m',-s}(θ_j) e^{i m' ϕ_k}\, e^{-imϕ_k}\, \frac{2π}{2j+1} \\
    &= (-1)^{s}\, \frac{2π}{2j+1} \sum_{ℓ,m'} {}_sf_{ℓ,m'}\, \sqrt{\frac{2ℓ+1}{4π}}\, d_{ℓ}^{m',-s}(θ_j) \sum_{k=0}^{2j}e^{i (m'-m) ϕ_k}.
\end{aligned}
```
We can evaluate this last sum easily:
```math
  \sum_{k=0}^{2j}e^{i (m'-m) ϕ_k} = \begin{cases}
    2j+1 & m'-m = n(2j+1)\ \mathrm{for}\ n\in\mathbb{Z}, \\
    0 & \mathrm{otherwise}.
  \end{cases}
```
This allows us to simplify as

```math
\begin{aligned}
    {}_{s}\tilde{f}_{j}(m) = (-1)^{s}\, 2π \sum_{ℓ,m'} {}_sf_{ℓ,m'}\, \sqrt{\frac{2ℓ+1}{4π}}\, d_{ℓ}^{m',-s}(θ_j),
\end{aligned}
```
where ``m'`` ranges over ``m + n(2j+1)`` for all ``n\in \mathbb{Z}`` such that ``|m + n(2j+1)| \leq ℓ``
— that is, all ``n\in \mathbb{Z}`` such that
```math
  \left \lceil \frac{-ℓ-m}{2j+1} \right \rceil \leq n \leq \left \lfloor \frac{ℓ-m}{2j+1} \right \rfloor.
```


## Matrix representation

Usually, we would take the sum over ``ℓ`` ranging from ``\mathrm{max}(|m|,|s|)`` to ``L``, and the sum
over ``m'`` ranging over ``m + n(2j+1)`` for all ``n\in \mathbb{Z}`` such that ``|m + n(2j+1)| \leq ℓ``.
However, we can also consider these sums to range over all possible
values of ``ℓ, m'``, and just set the coefficient to zero whenever
these conditions are not satisfied.  In that case, we can again think
of this as a (much larger) vector-matrix equation reading
```math
  {}_s\tilde{𝐟} = (-1)^s\, 2π\, {}_s𝐝\, {}_s𝐟,
```
where the index on ``{}_s\tilde{𝐟}`` loops over ``j`` and
``m``, the index on ``{}_s𝐟`` loops over ``ℓ`` and ``m'``,
and the indices on ``{}_s𝐝`` loop over each of those pairs.


## De-aliasing

While it is *far* simpler to simply invert the full ``{}_s𝐝``
matrix, its size scales as ``L^4``, which means that it very quickly
becomes impractical to store and manipulate the full matrix.  In CMB
astronomy, for example, it is not uncommon to use ``L`` into the tens
of thousands, which would make the full matrix utterly impractical to
use.

However, the matrix has a fairly sparse structure, with the number of
*nonzero* elements scaling as ``L^3``.  More particularly, the
sparsity has a fairly special structure, where the full matrix is
mostly block diagonal, along with some sparse upper triangular
elements.  Of course, the goal is to solve the linear equation.  For
that, the first obvious choice is an LU decomposition.  Unfortunately,
the L and U components are *not* sparse.  A second obvious choice is
the QR decomposition, which is more tailored to the structure of this
matrix — the Q factor being essentially just the block diagonal, and
the R factor being a somewhat less sparse upper triangle.

In principle, this alone could delay the impracticality threshold —
though still not enough for CMB astronomy.  We can use the unusual
structure to solve the linear equation in a more piecewise fashion,
with fairly low memory overhead.  Essentially, we start with the
highest-``|k|`` values, and solve for the corresponding
highest-``|m|`` values.  Those harmonics will alias to other
frequencies in ``θ_j`` rings with ``j < |k|``.  But crucially, we
know *how* they alias, and can simply remove them from the Fourier
transforms of those rings.  We then repeat, solving for the
next-highest ``|k|`` values, and so on.

## Spin weights other than zero

The scheme above places a ring of ``2j+1`` points for each ``j ∈
|s|:L``, and for ``s = 0`` that works well.  For any other spin weight
it is badly conditioned, and no choice of the colatitudes cures it.
The reason is visible in the behavior of the harmonics near the poles:
``{}_{s}λ_{ℓ,m}(θ)`` is proportional to ``\sin^{|m+s|}(θ/2)\,
\cos^{|m-s|}(θ/2)``, so near the north pole a function of spin weight
``s`` is dominated by the modes with ``m`` near ``-s``, and near the
south pole by those near ``+s``.  A small ring near the north pole
measures the frequencies ``|m| ≤ j``, but it sees the ones near ``m =
+j`` only weakly, while the aliases of ``m = -(j+1), -(j+2), …`` land on
its coefficients at nearly full strength.  Each step of the de-aliasing
then amplifies the errors of the steps before it, and the error of a
round trip grows by more than an order of magnitude with each unit of
``L``: at ``s = 2`` in double precision, half the digits are gone by
``L = 10``, and all of them by ``L = 14``.

The fix is to center each polar ring's window of frequencies on the
modes that dominate near its pole.  What the analysis needs is that
every ``m`` be measured by as many rings as there are modes with that
``m``, which is ``L - \max(|m|, |s|) + 1``.  The windows ``|m| ≤ a``
and ``|m| ≤ a+2d`` together cover every ``m`` exactly as often as the
windows ``|m+d| ≤ a+d`` and ``|m-d| ≤ a+d`` do, so pairs of the original
windows can be replaced by pairs of windows of equal size, one centered
on ``-d\,\mathrm{sign}(s)`` for a ring in the northern hemisphere and
one on ``+d\,\mathrm{sign}(s)`` for a ring in the southern.  With ``d =
|s|`` wherever possible, the condition number at ``s = 2`` and ``L =
12`` drops from ``7×10^{12}`` to about 400.  The number of rings and of
points is unchanged; [`minimal_rings`](@ref) gives the details.

The price is in the de-aliasing.  A frequency ``m′`` outside a ring's
window must be removed from that ring's coefficients before the
frequency it aliases to can be solved for.  With every window centered
on 0 this ordering is simply that of decreasing ``|m|``, but with two
windows of equal size centered on ``±d`` some frequencies in each alias
into the other, and neither can be solved first.  Those frequencies
must be solved together.  The groups that must be solved together are
the strongly connected components of the graph whose edges run from
each such ``m′`` to the frequency it aliases to, and solving the groups
in topological order restores the triangular structure.  In every case
I have measured the groups hold at most ``4|s|-1`` values of ``m``,
however large ``L`` is, so the cost remains ``O(L^3)``.

Even so, the sample points become badly conditioned as ``L`` grows, for
every spin weight — the error of a round trip in double precision is
about ``10^{-12}`` at ``L = 32``, ``10^{-10}`` at ``L = 48`` and
``10^{-8}`` at ``L = 64`` for ``s = 0``, and grows faster for larger
``|s|`` — which is the limitation mentioned at the top of this page.


## Implementation

[`SSHTMinimal`](@ref) precomputes, at construction, the table
`Λ[i, r]` of ``{}_{s}λ_{ℓ,m}(θ_r)`` for every mode ``i = (ℓ, m)`` on
every ring ``r``, along with the groups of ``m`` values and the LU
decomposition of the matrix of each.  A group's matrix couples its
modes to the Fourier coefficients that measure its ``m`` values on
every ring whose window includes them; a mode enters a coefficient
whenever its ``m`` is congruent, modulo the size of the ring, to the
frequency that coefficient measures.  (Evaluating the
``{}_{s}λ_{ℓ,m}`` on the fly instead, recomputing the recursion once
per ring per ``m``, costs more than storing them.)

The following pseudo-code summarizes the analysis algorithm:
```julia
# Fourier coefficients of each ring, normalized as (1/N) Σₖ f(ϕₖ) exp(-imϕₖ)
for r ∈ rings
    F[r] = fft(f[pixels of r]) / Nϕ[r]
end

for group ∈ groups  # in topological order
    # Every other mode that reaches these coefficients has already been removed
    rhs = [F[r][mod(m, Nϕ[r]) + 1] for (r, m) ∈ coefficients(group)]
    f̃[modes(group)] = lu(group) \ rhs

    # Remove this group's modes from the coefficients of every ring they reach
    for r ∈ rings, i ∈ modes(group)
        F[r][mod(m[i], Nϕ[r]) + 1] -= f̃[i] * Λ[i, r]
    end
end
```

Synthesis needs no ordering at all, because every mode's contribution
to every ring is known:
```julia
for r ∈ rings
    F[r] .= 0
    for i ∈ modes  # aliased or not
        F[r][mod(m[i], Nϕ[r]) + 1] += f̃[i] * Λ[i, r]
    end
    f[pixels of r] = bfft(F[r])  # Σₘ Fₘ exp(imϕₖ)
end
```
