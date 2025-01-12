# Algorithm for computing ``H``

The ``H`` array, as given by [Gumerov_2015](@citet), is related to Wigner's (small) ``d`` matrices —
which is itself related to the (big) ``\mathfrak{D}`` matrices and the various spin-weighted
spherical harmonics ``{}_{s}Y_{\ell,m}`` — via

```math
d_{\ell}^{m',m} = \epsilon_{m'} \epsilon_{-m} H_{\ell}^{m',m},
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

## Steps to compute ``H``

The following describes various details that are not spelled out correctly by
[Gumerov_2015](@citet).  All equation numbers refer to that paper.

Because of the symmetries noted above, we only compute ``H^{m', m}_n`` with ``m
≥ |m'|`` — roughly one quarter of all possible values.  Furthermore, for
computations of spin-weighted spherical harmonics of weight ``s``, we only need
to compute values with ``|m'| ≤ |s|``, which constitutes a dramatic savings
when ``|s| ≪ ℓₘₐₓ``.  The data are stored in the array `Hwedge`.

However, some parts of this calculation require calculating terms with
``m=n+1`` — whereas such elements of ``d`` and ``\mathfrak{D}`` are considered
zero.  For this purpose, we need additional storage.  Rather than allocating
extra space, or requiring some additional workspace to be passed in, we can
actually use parts of the input ``H`` data space for temporary storage while
these extra terms are needed, which is before those parts of the storage are
needed.  Specifically, we need this additional storage for
``H^{0, m}_{n_\mathrm{max}+1}`` with ``m \in [0, n_\mathrm{max}+1]``, and we
can use the storage destined for ``H^{-1, m}_{n_\mathrm{max}}`` with
``m \in [1, n_\mathrm{max}]``.  But this leaves two more indices, which we just
store as individual variables — `HΩ` and `HΨ` — representing the last and
second-to-last of these additional elements stored.


### Step 1

Set ``H_{0}^{0,0}=1``.


### Step 2

Compute values ``H^{0,m}_{n}(β)`` for ``m=0,\ldots,n`` and ``H^{0,m}_{n+1}(β)``
for ``m=0,\ldots,n+1``.  Using Eq. (32), we see that within Gumerov and
Duraiswami's conventions
```math
\begin{aligned}
  H^{0,m}_{n}(β) &= (-1)^m \sqrt{\frac{(n-|m|)!}{(n+|m|)!}} P^{|m|}_{n}(\cos β) \\
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
  e_{n,m} &= \frac{1}{2n} \sqrt{\frac{2}{2-\delta_0^{m-1}}} \sqrt{\frac{(n+m)(n+m-1)(2n+1)}{2n-1}}.
\end{aligned}
```

Now, we can directly obtain a recurrence relation for
``H^{0,m}_{n} = P̄_{n,|m|} / \sqrt{k_m (2n+1)} `` from those expressions:
```math
\begin{aligned}
  H^{0,0}_{0} &= 1 \\
  H^{0,0}_{1} &= \cos β \\
  H^{0,1}_{1} &= \sqrt{1/2} \sin β \\
  H^{0,0}_{n} &= \cos β H^{0,0}_{n-1} - b̄_n \sin β H^{0,1}_{n-1} \\
  H^{0,m}_{n} &=
    c̄_{n,m} \cos β H^{0,m}_{n-1}
    - \sin β \left[ d̄_{n,m} H^{0,m+1}_{n-1} - ē_{n,m} H^{0,m-1}_{n-1} \right],
\end{aligned}
```
where the coefficients are given by
```math
\begin{aligned}
  b̄_n &= \sqrt{\frac{n-1}{n}} \\
  c̄_{n,m} &= \frac{1}{n} \sqrt{(n+m)(n-m)} \\
  d̄_{n,m} &= \frac{1}{2n} \sqrt{(n-m)(n-m-1)} \\
  ē_{n,m} &= \frac{1}{2n} \sqrt{(n+m)(n+m-1)}.
\end{aligned}
```
Note that the coefficients all simplified (in fact, ``a_n`` disappeared), without any
increase in the complexity of the recurrence relations themselves.  Rewriting Belikov's
algorithm explicitly in terms of the ``H^{0,m}_{n}`` also allows us to avoid an extra
normalization step.

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
``\mathrm{sgn}(0)=0``), and we have ``a^{m}_{n} = a^{-m}_{n}``.  Also note
that these coefficients *only* appear in this step, and because of how they
appear (specifically, because ``b`` always appears with argument ``n+1``), we
can factor out the denominators in the definitions of the constants.  We
obtain this simplified formula
```math
H^{1, m}_{n}
  = -\frac{1}{\sqrt{n(n+1)}} \left[
      \frac{\bar{b}^{−m−1}_{n+1} (1−\cos β)}{2} H^{0, m+1}_{n+1}
      + \frac{\bar{b}^{ m−1}_{n+1} (1+\cos β)}{2} H^{0, m−1}_{n+1}
      + \bar{a}^{m}_{n} \sin β H^{0, m}_{n+1}
    \right],
```
with
```math
\bar{a}^{m}_{n} = \sqrt{(n+m+1)(n-m+1)},
```
```math
\bar{b}^{m}_{n+1} = \sqrt{(n-m)(n-m+1)}.
```


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
Note that we can drop the factor of ``1/2``, and *for this case only* the sign
is always +1.


### Step 5
Recursively compute ``H^{m'−1, m}_{n}(β)`` for ``m'=0,\ldots,−n+1``,
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

*If* we include the cost of computing all these constants in a single call to
the ``H`` recurrence, it can be much cheaper to compute each constant as needed
within the algorithm, rather than computing them all at once at the beginning
of the algorithm — but only for very small computations, such as those
involving ``n_{\mathrm{max}} ≈ 10``.  Beyond this, despite the storage
penalties for all those constants, it turns out to be better to pre-compute
them.  However, it should be noted that the fractional cost of storing the
constants is ``\sim 3/n_{\mathrm{max}}`` compared to just storing ``H`` itself,
so this will never be a very significant amount of space.

On the other hand, if we can pre-compute the constants just once, and store
them between multiple calls to the ``H`` recurrence, then it is always
advantageous to do so — typically by factors of 2 or 3 in speed.  The only
difficulty here is ensuring that each call to the recurrence has access to the
constants, which can be a little awkward when using multiple processes and/or
threads.  However, it should be thread safe, since we only need to read those
constants within the ``H`` recurrence.  All in all, I conclude that it is
probably not worth the effort to maintain separate versions of the recurrence
for pre-computed and on-the-fly constants.
