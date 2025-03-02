# Internal functions

There are various functions that are only used internally, some of which are likely
to be deprecated in the near future.  These are documented here for completeness.

## ``H`` recursion and ALFs

The fundamental algorithm is the ``H`` recursion, which is the core computation
needed for Wigner's ``d`` and ``𝔇`` matrices, and the spin-weighted spherical
harmonics ``{}_{s}Y_{\ell,m}``, as well as `map2salm` functions.

```@autodocs
Modules = [SphericalFunctions]
Pages   = ["Hrecursion.jl"]
```

Internally, the ``H`` recursion relies on calculation of the Associated Legendre
Functions (ALFs), which can also be called on their own:

```@autodocs
Modules = [SphericalFunctions]
Pages   = ["associated_legendre.jl"]
```

The function ``{}_{s}\lambda_{\ell,m}`` is defined as essentially ``{}_{s}Y_{\ell,0}``, and is important internally for computing the ALFs.  We have some important utilities for computing it:

```@docs
SphericalFunctions.λ_recursion_initialize
SphericalFunctions.λ_iterator
SphericalFunctions.AlternatingCountdown
```



## ``Y``, ``d``, and ``D``

Various `d`, `D`, and `sYlm` functions are important in the main API.  Their
names and signatures have been tweaked from older versions of this package.  The
only one with remaining documentation is [`ₛ𝐘`](@ref), which could probably be
replaced by [`sYlm_values`](@ref), except that the default pixelization is
[`golden_ratio_spiral_rotors`](@ref), which makes it very convenient for
interacting with [`SSHT`](@ref).

```@docs
ₛ𝐘
SphericalFunctions.Y
SphericalFunctions.d
SphericalFunctions.D
```


# Transformation

The newer [`SSHT`](@ref) interface is more efficient for most purposes, but this
package used to use functions named `map2salm`, which is still present, but may
be deprecated.

```@autodocs
Modules = [SphericalFunctions]
Pages   = ["map2salm.jl"]
```
