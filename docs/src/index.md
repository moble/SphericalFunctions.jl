# Introduction

Julia package for evaluating and transforming Wigner's 𝔇 matrices, and spin-weighted spherical
harmonics (which includes the ordinary scalar spherical harmonics).  These functions are evaluated
directly in terms of quaternions, as well as in the more standard forms of spherical coordinates and
Euler angles.[^1]

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

Note that numerous other packages cover some of these use cases, including
[`FastTransforms.jl`](https://JuliaApproximation.github.io/JuliaApproximation/FastTransforms.jl/),
[`FastSphericalHarmonics.jl`](https://eschnett.github.io/FastSphericalHarmonics.jl/), and [this
package](https://github.com/Jutho/WignerSymbols.jl).  However, I need support for higher-precision
numbers — even at the cost of some speed — which is what this package provides.


## Contents

```@contents
Depth = 4
```

## Function list

The following list contains the public functions inside the `SphericalFunctions` module.

```@index
Modules = [SphericalFunctions]
```


[^1]:
    Euler angles are quite generally a very poor choice for computing with rotations.  (The only
    context in which they may be preferred is when *analytically* integrating some analytically
    known functions.)  Almost universally, it is best to use quaternions when computing with
    rotations.  All the computations done within this package use quaternions; the user interfaces
    involving Euler angles essentially convert to/from quaternions.  While the calculations needed
    for those conversions would still need to be done if this package used Euler angles internally —
    meaning that this approach is as efficient as any — that work can be avoided entirely if you
    work with quaternions directly.
