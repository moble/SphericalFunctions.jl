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


This is a Julia package for evaluating and transforming Wigner's 𝔇 matrices, the associated Legendre
functions, and spin-weighted spherical harmonics (which includes the ordinary scalar spherical
harmonics).  See [the documentation](https://moble.github.io/SphericalFunctions.jl/) for more
details.


## Installation

```bash
using Pkg
Pkg.add("SphericalFunctions")
```


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
