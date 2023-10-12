# Introduction

This is a Julia package for evaluating and transforming Wigner's 𝔇 matrices,
and spin-weighted spherical harmonics ``{}_{s}Y_{\ell,m}`` (which includes the
ordinary scalar spherical harmonics).  Because [*both* 𝔇 *and* the harmonics
are most correctly considered](@cite Boyle_2016) functions on the rotation group
``𝐒𝐎(3)`` — or more generally, the spin group ``𝐒𝐩𝐢𝐧(3)`` that covers it —
these functions are evaluated directly in terms of quaternions.  Concessions are
also made for more standard forms of spherical coordinates and Euler angles.[^1]
Among other applications, those functions permit "synthesis" (evaluation of the
spin-weighted spherical functions) of spin-weighted spherical harmonic
coefficients on regular or distorted grids.  This package also includes
functions enabling efficient "analysis" (decomposition into mode coefficients)
of functions evaluated on regular grids to high order and accuracy.

These quantities are computed using recursion relations, which makes it possible
to compute to very high ℓ values.  Unlike direct evaluation of individual
elements, which would generally cause overflow or underflow beyond ℓ≈30 when
using double precision, these recursion relations should be valid for far higher
ℓ values.  More precisely, when using *this* package, `Inf` values appear
starting at ℓ=128 for `Float16`, but I have not yet found any for values up to
at least ℓ=1024 with `Float32`, and presumably far higher for `Float64`.
`BigFloat` also works, and presumably will not overflow for any ℓ value that
could reasonably fit into computer memory — though it is far slower.  Also note
that [`DoubleFloats`](https://github.com/JuliaMath/DoubleFloats.jl) will work,
and achieve significantly greater accuracy (but no greater ℓ range) than
`Float64`.  In all cases, results are typically accurate to roughly ℓ times the
precision of the input quaternion.

The conventions for this package are mostly inherited from — and are described
in detail by — its predecessors found
[here](https://moble.github.io/spherical_functions/) and
[here](https://moble.github.io/spherical/).

Note that numerous other packages cover some of these use cases, including
[`FastTransforms.jl`](https://JuliaApproximation.github.io/FastTransforms.jl/),
[`FastSphericalHarmonics.jl`](https://eschnett.github.io/FastSphericalHarmonics.jl/dev/),
[`WignerSymbols.jl`](https://github.com/Jutho/WignerSymbols.jl), and
[`WignerFamilies.jl`](https://github.com/xzackli/WignerFamilies.jl).  However, I
need support for quaternions (via
[`Quaternionic.jl`](https://github.com/moble/Quaternionic.jl)) and for
higher-precision numbers — even at the cost of a very slight decrease in speed
in some cases — which are what this package provides.


[^1]:
    Euler angles are quite generally a very poor choice for computing with
    rotations.  (The only context in which they may be preferred is when
    *analytically* integrating some analytically known functions.)  Almost
    universally, it is best to use quaternions when computing with rotations.
    All the computations done within this package use quaternions; the user
    interfaces involving Euler angles essentially convert to/from quaternions.
    While the calculations needed for those conversions would still need to be
    done if this package used Euler angles internally — meaning that this
    approach is as efficient as any — that work can be avoided entirely if you
    work with quaternions directly.
