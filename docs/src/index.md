# Introduction

Julia package for evaluating and transforming Wigner's 𝔇 matrices, Wigner's 3-j symbols, and
spin-weighted (and scalar) spherical harmonics.  These functions are evaluated directly in terms of
quaternions, as well as in the more standard forms of spherical coordinates and Euler
angles.[^1]

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
numbers — even at the cost of speed — which is what this package provides.


## Contents

```@contents
Depth = 4
```

## Function list

The following list contains the public functions inside the `Spherical` module.

```@index
Modules = [Spherical]
```


[^1]:
    Euler angles are pretty much [the worst things
    ever](http://moble.github.io/spherical/#euler-angles) and it makes me feel bad even supporting
    them.  Quaternions are faster, more accurate, basically free of singularities, more intuitive,
    and generally easier to understand.  You can work entirely without Euler angles (I certainly
    do).  You absolutely never need them.  But if you're so old fashioned that you really can't give
    them up, they are fully supported.
