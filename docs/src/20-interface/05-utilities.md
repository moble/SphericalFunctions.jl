# Utilities

While not usually the star of the show, the following utilities can be
quite helpful for actually using the rest of the code.


## Complex powers

One common task we find when working with spherical functions is the
computation of a range of integer powers of some complex number — so
much so that it can be best to pre-compute the powers and cache their
values.  While a naive approach is generally quite accurate, and
reasonably fast, we can do a little better with a specialized routine.

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
use.  It is the same ordering for half-integer indices, with ``ℓ``
starting at ``1/2`` rather than ``0``: the closed forms below hold for
both kinds of index, and the indices in one call must all be of one
kind.  Wigner's ``𝔇`` and ``d`` matrices are *not* stored this way:
they are indexed by ``(ℓ, m', m)`` through the views a
[`WignerCalculator`](@ref) returns, so no separate index arithmetic is
needed for them.

```@autodocs
Modules = [SphericalFunctions]
Pages   = ["indices/mode_ordering.jl", "mode_weights/mode_weights.jl"]
Order   = [:module, :type, :constant, :function, :macro]
```


## Types of rotor data

These functions report how the package reads the rotor data given to a
calculator: the floating-point type it computes in, and the number of
rotors.  [`floattype`](@ref SphericalFunctions.floattype), which is
documented with the accessors, and `nrotors` are public; the others
are the internal helpers that check the data and word the refusals.
The same file holds `spinor_phases` and `half_angles`, which compute
the phases that a calculator stores for each rotor; they are
documented among the [internal functions](@ref "Rotor phases").

```@autodocs
Modules = [SphericalFunctions]
Pages   = ["calculators/rotors.jl"]
Filter  = f -> f ∉ (SphericalFunctions.spinor_phases, SphericalFunctions.half_angles)
```
