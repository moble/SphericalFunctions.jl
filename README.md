# Spherical Functions in Julia

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://moble.github.io/SphericalFunctions.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://moble.github.io/SphericalFunctions.jl/dev)
[![Build Status](https://github.com/moble/SphericalFunctions.jl/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/moble/SphericalFunctions.jl/actions/workflows/tests.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/moble/SphericalFunctions.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/moble/SphericalFunctions.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/S/SphericalFunctions.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/report.html)

[![DOI](https://zenodo.org/badge/381490836.svg)](https://zenodo.org/badge/latestdoi/381490836)


Julia package for evaluating and transforming Wigner's 𝔇 matrices, and spin-weighted spherical
harmonics (which includes the ordinary scalar spherical harmonics).  These functions are evaluated
directly in terms of quaternions, as well as in the more standard forms of spherical coordinates and
Euler angles.<sup>[1](#1-euler-angles-are-inadequate)</sup> Among other applications, those
functions permit "synthesis" (evaluation of the spin-weighted spherical functions) of spin-weighted
spherical harmonic coefficients on regular or distorted grids.  This package also includes functions
enabling efficient "analysis" (decomposition into mode coefficients) of functions evaluated on
regular grids to high order and accuracy.

These quantities are computed using recursion relations, which makes it possible to compute to very
high ℓ values.  Unlike direct evaluation of individual elements, which will generally cause overflow
or underflow beyond ℓ≈30 when using double precision, these recursion relations should be valid for
far higher ℓ values.  More precisely, when using *this* package, `Inf` values appear starting at
ℓ=22 for `Float16`, ℓ=183 for `Float32`, and ℓ=1474 for `Float64`.  `BigFloat` also works, and
presumably will not overflow for any ℓ value that could reasonably fit into computer memory — though
it is far slower.  Also note that [`DoubleFloats`](https://github.com/JuliaMath/DoubleFloats.jl)
will work, and achieve significantly greater accuracy but no greater ℓ range than `Float64`.  The
results are accurate to roughly ℓ times the precision of the input quaternion.

The conventions for this package are described in detail on [this
page](https://moble.github.io/spherical_functions/).

Note that numerous other packages cover some of these use cases, including
[`FastTransforms.jl`](https://JuliaApproximation.github.io/FastTransforms.jl/),
[`FastSphericalHarmonics.jl`](https://eschnett.github.io/FastSphericalHarmonics.jl/dev/),
[`WignerSymbols.jl`](https://github.com/Jutho/WignerSymbols.jl), and
[`WignerFamilies.jl`](https://github.com/xzackli/WignerFamilies.jl).  However, I need support for
quaternions (via [`Quaternionic.jl`](https://github.com/moble/Quaternionic.jl)) and for
higher-precision numbers — even at the cost of a very slight decrease in speed in some cases — which
are what this package provides.


## Installation

```bash
using Pkg
Pkg.add("SphericalFunctions")
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

It may be worthwhile simply passing [X-numbers](https://github.com/moble/XNumbers.jl) to these
algorithms if higher orders are needed.  (Though using `BigFloat` *may* achieve a similar objective,
it would probably be far slower.)  The actual recommendation of Fukushima is more sophisticated —
just using X-numbers in the core calculation — but it looks like the simpler approach wouldn't be
*too* much slower.

The other major functionality of this package is `map2salm`, which decomposes function values on
regular grids into mode weights (coefficients).  The approach used here is taken from [Reinecke and
Seljebotn](https://dx.doi.org/10.1051/0004-6361/201321494), with weights based on [Waldvogel's
method](https://doi.org/10.1007/s10543-006-0045-4).


<br/>

---

###### <sup>1</sup> Euler angles are inadequate

Euler angles are quite generally a very poor choice for computing with rotations.  (The only context
in which they may be preferred is when *analytically* integrating some analytically known
functions.)  Almost universally, it is best to use quaternions when computing with rotations.  All
the computations done within this package use quaternions; the user interfaces involving Euler
angles essentially convert to/from quaternions.  While the calculations needed for those conversions
would still need to be done if this package used Euler angles internally — meaning that this
approach is as efficient as any — that work can be avoided entirely if you work with quaternions
directly.
