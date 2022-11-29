# Sampling theorems

[McEwen and Wiaux](https://arxiv.org/abs/1110.6298) provide a very thorough review of the literature
on sampling theorems related to spin-weighted spherical harmonics up to 2011, and presents a
``∼2L²`` algorithm, though [Elahi, Khalid, Kennedy, and McEwen](https://arxiv.org/abs/1809.01321)
have obtained the optimal result that scales as ``∼L²``.

The downside of the EKKM algorithm is that the ``θ`` values at which to sample have to be obtained
by iteratively minimizing the condition numbers of various matrices (which are involved in the
computation itself).  This expensive step only has to be performed once per choice of spin ``s`` and
maximum ``ℓ`` value ``L``.  Otherwise, the results of this algorithm seem to be relatively good — at
least for ``L`` up to 64.  This does not compare favorably with the MW algorithm, which has slowly
growing errors through ``L = 4096``.

## EKKM analysis

The EKKM analysis looks like the following (with some notational changes).  We begin by defining
```math
  {}_{s}\tilde{f}_{\theta}(m) := \int_0^{2\pi} {}_sf(\theta, \phi)\, e^{-im\phi}\, d\phi.
```
We will denote the vector of these quantities for all values of $\theta$ as
${}_{s}\tilde{\mathbf{f}}_m$.  Inserting the ${}_sY_{\ell,m}$ expansion for ${}_sf(\theta, \phi)$,
and performing the integration using orthogonality of complex exponentials, we can find that
```math
  {}_{s}\tilde{f}_{\theta}(m) = (-1)^s\, 2\pi \sum_{\ell=\Delta}^L \sqrt{\frac{2\ell+1}{4\pi}}\, d_{m,-s}^{\ell}(\theta)\, {}_sf_{\ell,m}.
```
Now, denoting the vector of ${}_sf_{\ell,m}$ for all values of $\ell$ as ${}_s\mathbf{f}_m$, we can
write this as a matrix-vector equation:
```math
  {}_{s}\tilde{\mathbf{f}}_m = (-1)^s\, 2\pi\, {}_s\mathbf{d}_{m}\, {}_s\mathbf{f}_m.
```
We are effectively measuring the ${}_{s}\tilde{\mathbf{f}}_m$ values, we can easily construct the
${}_s\mathbf{d}_{m}$ matrix, and we are seeking the ${}_s\mathbf{f}_m$ values, so we can just invert
this equation to solve for the latter.


## Discretize Fourier transform

Now, the only flaw in this analysis is that we have undersampled everywhere except $\ell = L$, which
means that the second equation (re-expressing the Fourier transforms as a sum using orthogonality of
complex exponentials) isn't quite right; in general there is some folding due to aliasing of
higher-frequency modes, so we need an additional sum over $|m'|>|m|$.  Or perhaps more precisely,
the first equation isn't actually what we implement.  It should look more like this:
```math
  {}_{s}\tilde{f}_{j}(m) := \sum_{k=0}^{2j} {}_sf(\theta_j, \phi_k)\, e^{-im\phi_k}\, \Delta \phi,
```
where $\phi_k = \frac{2\pi k}{2j+1}$, and $\Delta \phi = \frac{2\pi}{2j+1}$.  (Recall the subtle
notational distinction common in time-frequency analysis that $\tilde{s}(t_j) = \Delta t
\tilde{s}_j$, which would suggest we use ${}_{s}\tilde{f}_{j}(m) = \Delta \phi\,
{}_{s}\tilde{f}_{j,m}$.)  Next, we can insert the expansion for ${}_sf(\theta, \phi)$:

```math
\begin{aligned}
    {}_{s}\tilde{f}_{j}(m)
    &= \sum_{k=0}^{2j} \sum_{\ell,m'} {}_sf_{\ell,m'}\, {}_sY_{\ell,m'}(\theta_j, \phi_k)\, e^{-im\phi_k}\, \Delta \phi \\
    &= \sum_{k=0}^{2j} \sum_{\ell,m'} {}_sf_{\ell,m'}\, (-1)^{s}\, \sqrt{\frac{2\ell+1}{4\pi}}\, d_{\ell}^{m',-s}(\theta_j) e^{i m' \phi_k}\, e^{-im\phi_k}\, \frac{2\pi}{2j+1} \\
    &= (-1)^{s}\, \frac{2\pi}{2j+1} \sum_{\ell,m'} {}_sf_{\ell,m'}\, \sqrt{\frac{2\ell+1}{4\pi}}\, d_{\ell}^{m',-s}(\theta_j) \sum_{k=0}^{2j}e^{i (m'-m) \phi_k}.
\end{aligned}
```
We can evaluate this last sum easily:
```math
  \sum_{k=0}^{2j}e^{i (m'-m) \phi_k} = \begin{cases}
    2j+1 & m'-m = n(2j+1)\ \mathrm{for}\ n\in\mathbb{Z}, \\
    0 & \mathrm{otherwise}.
  \end{cases}
```
This allows us to simplify as

```math
\begin{aligned}
    {}_{s}\tilde{f}_{j}(m) = (-1)^{s}\, 2\pi \sum_{\ell,m'} {}_sf_{\ell,m'}\, \sqrt{\frac{2\ell+1}{4\pi}}\, d_{\ell}^{m',-s}(\theta_j),
\end{aligned}
```
where $m'$ ranges over $m + n(2j+1)$ for all $n\in \mathbb{Z}$ such that $|m + n(2j+1)| \leq \ell$
— that is, all $n\in \mathbb{Z}$ such that
```math
  \left \lceil \frac{-\ell-m}{2j+1} \right \rceil \leq n \leq \left \lfloor \frac{\ell-m}{2j+1} \right \rfloor.
```


## Matrix representation

Usually, we would take the sum over $\ell$ ranging from $\mathrm{max}(|m|,|s|)$ to $L$, and the sum
over $m'$ ranging over $m + n(2j+1)$ for all $n\in \mathbb{Z}$ such that $|m + n(2j+1)| \leq \ell$.
However, we can also consider these sums to range over all possible values of $\ell, m'$, and just
set the coefficient to zero whenever these conditions are not satisfied.  In that case, we can again
think of this as a (much larger) vector-matrix equation reading
```math
  {}_s\tilde{\mathbf{f}} = (-1)^s\, 2\pi\, {}_s\mathbf{d}\, {}_s\mathbf{f},
```
where the index on ${}_s\tilde{\mathbf{f}}$ loops over $j$ and $m$, the index on ${}_s\mathbf{f}$
loops over $\ell$ and $m'$, and the indices on ${}_s\mathbf{d}$ loop over each of those pairs.


## De-aliasing

While it is *far* simpler to simply invert the full ${}_s\mathbf{d}$ matrix, its size scales as
$L^4$, which means that it very quickly becomes impractical to store and manipulate the full matrix.
In CMB astronomy, for example, it is not uncommon to use $L$ into the tens of thousands, which would
make the full matrix utterly impractical to use.

However, the matrix has a fairly sparse structure, with the number of *nonzero* elements scaling as
$L^3$.  More particularly, the sparsity has a fairly special structure, where the full matrix is
mostly block diagonal, along with some sparse upper triangular elements.  Of course, the goal is to
solve the linear equation.  For that, the first obvious choice is an LU decomposition.
Unfortunately, the L and U components are *not* sparse.  A second obvious choice is the QR
decomposition, which is more tailored to the structure of this matrix — the Q factor being
essentially just the block diagonal, and the R factor being a somewhat less sparse upper triangle.

In principle, this alone could delay the impracticality threshold — though still not enough for CMB
astronomy.  We can use the unusual structure to solve the linear equation in a more piecewise
fashion, with fairly low memory overhead.  Essentially, we start with the highest-$|k|$ values, and
solve for the corresponding highest-$|m|$ values.  Those harmonics will alias to other frequencies
in $\theta_j$ rings with $j < |k|$.  But crucially, we know *how* they alias, and can simply remove
them from the Fourier transforms of those rings.  We then repeat, solving for the next-highest $|k|$
values, and so on.
