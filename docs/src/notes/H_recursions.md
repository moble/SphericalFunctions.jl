# Algorithm for computing ``H``

The ``H`` array, as given by [Gumerov and
Duraiswami](https://arxiv.org/abs/1403.7698), is related to Wigner's (small)
``d`` matrices — which is itself related to the (big) ``\mathfrak{D}`` matrices and
the various spin-weighted spherical harmonics ``{}_{s}Y_{\ell,m}`` — via

```math
d_{\ell}^{n,m} = \epsilon_n \epsilon_{-m} H_{\ell}^{n,m},
```

where

```math
\epsilon_k =
  \begin{cases}
  1 & k\leq 0, \\
  (-1)^k & k > 0.
\end{cases}
```

``H`` has various advantages over ``d``, including the fact that it can be efficiently
and robustly valculated via recurrence relations, and the following symmetry
relations:

```math
\begin{aligned}
  H^{m', m}_n(β) &= H^{m, m'}_n(β) \\
  H^{m', m}_n(β) &= H^{-m', -m}_n(β) \\
  H^{m', m}_n(β) &= (-1)^{n+m+m'} H^{-m', m}_n(π - β) \\
  H^{m', m}_n(β) &= (-1)^{m+m'} H^{m', m}_n(-β)
\end{aligned}
```

Because of these symmetries, we only need to evaluate at most 1/4 of all the
elements.

## Steps to compute $H$

The following describes various details that are not spelled out correctly by
[Gumerov and Duraiswami](https://arxiv.org/abs/1403.7698).  All equation
numbers refer to that paper.

Because of the symmetries noted above, we only compute ``H^{m', m}_n`` with ``m
≥ |m'|`` — roughly one quarter of all possible values.  Furthermore, for
computations of spin-weighted spherical harmonics of weight ``s``, we only need
to compute values with ``|m'| ≤ |s|``, which constitutes a dramatic savings
when ``|s| ≪ ℓₘₐₓ``.  The data are stored in the array `Hwedge`.

However, some parts of this calculation require terms corresponding to ``m =
|m'|-1``, which are retained in an additional array `Hv` (because they track
the V-shaped bottom of `Hwedge`).

Finally, some parts of this calculation require calculating terms with
``m=n+1`` — whereas such elements of ``d`` and ``\mathfrak{D}`` are considered
zero.  For this purpose, we define additional storage: `Hextra`.


### Step 1

Set ``H_{0}^{0,0}=1``.


### Step 2

Compute values ``H^{0,m}_{n}(β)`` for ``m=0,\ldots,n`` and ``H^{0,m}_{n+1}(β)``
for ``m=0,\ldots,n+1``.  Using Eq. (32), we see

```math
\begin{aligned}
  H^{0,m}_{n}(β) &= (-1)^m \sqrt{\frac{(n-|m|)!}{(n+|m|)!}} P^{|m|}_{n}(\cos β) \\
                 &= \frac{(-1)^m}{\sqrt{k (2n+1)}} P̄^{|m|}_{n}(\cos β).
\end{aligned}
```

Here, ``k=1`` for ``m=0``, and ``k=2`` for ``m>0``, and ``P̄`` is defined as

```math
  P̄ = \sqrt{\frac{k(2n+1)(n-m)!}{(n+m)!}} P.
```

We use the "fully normalized" associated Legendre functions (fnALF) ``P̄``
because, as explained by [Xing et
al. (2020)](https://doi.org/10.1007/s00190-019-01331-0), it is possible to
compute these values very efficiently and accurately, while also delaying the
onset of overflow and underflow.

NOTE: Though not specified by Gumerov and Duraiswami, there is not enough
information for step 4 unless we also use symmetry to set ``H^{1,0}_{n}`` here.
Similarly, step 5 needs additional information, which depends on setting
``H^{0, -1}_{n}`` from its symmetric equivalent ``H^{0, 1}_{n}`` in this step.


### Step 3
Compute ``H^{1,m}_{n}(β)`` for ``m=1,\ldots,n`` using relation (41).  Symmetry
and shift of the indices allow this relation to be written as
```math
b^{0}_{n+1} H^{1, m}_{n}
  = \frac{b^{−m−1}_{n+1} (1−\cos β)}{2} H^{0, m+1}_{n+1}
  − \frac{b^{ m−1}_{n+1} (1+\cos β)}{2} H^{0, m−1}_{n+1}
  − a^{m}_{n} \sin β H^{0, m}_{n+1}.
```
Here the constants are defined by
```math
a^{m}_{n} = \sqrt{\frac{(n+m+1)(n-m+1)} {(2n+1)(2n+3)}},
```
```math
b^{m}_{n} = \mathrm{sgn}(m) \sqrt{\frac{(n-m-1)(n-m)} {(2n-1)(2n+1)}}.
```
Note that all values are assumed to be zero whenever ``|m| > n``, we use
``\mathrm{sgn}(0)=1`` (unlike the common convention that
``\mathrm{sgn}(0)=0``), and we have ``a^{m}_{n} = a^{-m}_{n}``.

### Step 4
Recursively compute ``H^{m'+1, m}_{n}(β)`` for ``m'=1,\ldots,n−1``,
``m=m',...,n`` using relation (50) resolved with respect to ``H^{m'+1,
m}_{n}``:
```math
d^{m'}_{n} H^{m'+1, m}_{n}
  = d^{m'−1}_{n} H^{m'−1, m}_{n}
  − d^{m−1}_{n} H^{m', m−1}_{n}
  + d^{m}_{n} H^{m', m+1}_{n}
```
(where the last term drops out for ``m=n``).  The constants are defined by
```math
d^{m}_{n} = \frac{\mathrm{sgn}(m)}{2} \sqrt{(n-m)(n+m+1)}.
```


### Step 5
Recursively compute ``H^{m'−1, m}_{n}(β)`` for ``m'=−1,\ldots,−n+1``,
``m=−m',\ldots,n`` using relation (50) resolved with respect to ``H^{m'−1,
m}_{n}``:
```math
d^{m'−1}_{n} H^{m'−1, m}_{n}
  = d^{m'}_{n} H^{m'+1, m}_{n}
  + d^{m−1}_{n} H^{m', m−1}_{n}
  − d^{m}_{n} H^{m', m+1}_{n}
```
(where the last term drops out for ``m=n``).

NOTE: Although Gumerov and Duraiswami specify the loop over ``m'`` to start at
-1, I find it necessary to start at 0, or there will be missing information.
This also requires setting the ``H^{0, -1}_{n}`` components (for all ``n``)
before beginning this loop.


## Pre-computing constants versus computing on the fly

Each of the constants ``a^{m}_{n}``, ``b^{m}_{n}``, and ``c^{m}_{n}`` involves
divisions and square-roots, which can be very costly to compute.  It can be
advantageous to pre-compute the constants, and simply index the pre-computed
arrays rather than re-computing them on each recursion.

