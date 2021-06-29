# Spherical Functions

[![Test Status](https://github.com/moble/Spherical.jl/workflows/tests/badge.svg)](https://github.com/moble/Spherical.jl/actions)
[![Test Coverage](https://codecov.io/gh/moble/Spherical.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/moble/Spherical.jl)
[![Documentation
Status](https://github.com/moble/Spherical.jl/workflows/docs/badge.svg)](https://moble.github.io/Spherical.jl/dev)


Julia package for evaluating and transforming Wigner's 𝔇 matrices, Wigner's 3-j symbols, and
spin-weighted (and scalar) spherical harmonics.  These functions are evaluated directly in terms of
quaternions, as well as in the more standard forms of spherical coordinates and Euler
angles.<sup>[1](#1-euler-angles-are-awful)</sup>

These quantities are computed using recursion relations, which makes it possible to compute to very
high ℓ values.  Unlike direct evaluation of individual elements, which will generally cause overflow
or underflow beyond ℓ≈30 when using double precision, these recursion relations should be valid for
far higher ℓ values.  More precisely, `Inf` values appear starting at ℓ=22 for `Float16`, ℓ=183 for
`Float32`, and ℓ=1474 for `Float64`.  `BigFloat` also works, and presumably will not overflow for
any ℓ value that could reasonably fit into computer memory — though it is far slower.  Also note
that [`DoubleFloats`](https://github.com/JuliaMath/DoubleFloats.jl) will work, and achieve
significantly greater accuracy (but no greater range) than `Float64`.  The results are accurate to
roughly ℓ times the precision of the input quaternion.

The conventions for this package are described in detail on [this
page](http://moble.github.io/spherical/).

## Installation

```bash
julia -e 'using Pkg; pkg"add https://github.com/moble/Spherical.jl.git"'
```

## References

The most important routine in this package is the computation of the 𝔇 matrices — or more
specifically, of terms proportional to parts of the 𝔇 matrices.  This mostly follows the treatment
of [Gumerov and Duraiswami (2014)](https://arxiv.org/abs/1403.7698).  To seed the recursions they
present, we also need to calculate the associated Legendre functions.  Currently, this is done using
the "modified forward row method" of [Holmes and Featherstone
(2002)](https://doi.org/10.1007/s00190-002-0216-2).  Note that this is apparently the source of
overflow noted above.  Two other (though more complicated) methods have appeared more recently in
the literature, which could presumably extend these limits much further.  [Fukushima
(2012)](https://doi.org/10.1007/s00190-011-0519-2) showed that using "X-numbers" (wherein the
exponent is stored as a separate integer) in the core of the recursion could increase the range to
ℓ≈2³².  [Xing et al. (2020)](https://doi.org/10.1007/s00190-019-01331-0) showed that Fukushima's
results exhibited increased error for certain angles, whereas their Eqs. (12)—(14) could be used
directly to obtain results with greater accuracy for those certain angles, and comparable accuracy
for other angles.

It may be worthwhile simply implementing X-numbers as a subtype of `AbstractFloat`, and simply
passing them to these algorithms if higher orders are needed.  (Though using `BigFloat` would
achieve a similar objective, it would probably be far slower.)  The actual recommendation of
Fukushima is more sophisticated — just using X-numbers in the core calculation — but it looks like
this probably wouldn't be *too* much slower.
