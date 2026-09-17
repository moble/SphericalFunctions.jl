# Utilities

While not usually the star of the show, the following utilities can be quite helpful for actually
using the rest of the code.


## Complex powers

One common task we find when working with spherical functions is the computation of a range of
integer powers of some complex number — so much so that it can be best to pre-compute the powers and
cache their values.  While a naive approach is generally quite accurate, and reasonably fast, we can
do a little better with a specialized routine.

```@autodocs
Modules = [SphericalFunctions]
Pages   = ["utilities/complex_powers.jl"]
Order   = [:module, :type, :constant, :function, :macro]
```


## Sizes of and indexing into ``Y`` data

By ``Y`` data, we mean anything indexed like ``Y_{ℓ, m}`` modes: a
single vector, ordered by ``ℓ`` and then by ``m``, holding one entry
per mode.  This is the ordering that [`sYlm`](@ref),
[`sYlm_matrix`](@ref), [`ModeWeights`](@ref) and the transforms all
use.  Wigner's ``𝔇`` and ``d`` matrices are *not* stored this way in
version 3: they are indexed by ``(ℓ, m', m)`` through the views a
[`WignerCalculator`](@ref) returns, so no separate index arithmetic is
needed for them.

```@autodocs
Modules = [SphericalFunctions]
Pages   = ["mode_weights/indexing.jl", "mode_weights/mode_weights.jl"]
Order   = [:module, :type, :constant, :function, :macro]
```


## Combinatorics

Spherical functions frequently involve binomial coefficients and similar terms, with arguments
proportional to ``ℓ``, which we aim to allow to be very large — of order 1,000 or more.
Unfortunately, due to combinatorical explosions, this is frequently infeasible with naive
methods.  Here, we collect any specialized methods that help us beat the limits.

```@autodocs
Modules = [SphericalFunctions]
Pages   = ["utilities/utils.jl"]
```


