# Complete function list

The following list contains the functions, types, and constants of the
`SphericalFunctions` module, each linked to its docstring.  Every
docstring in the package appears somewhere in this manual — those of
the public interface on the pages that describe it, and those of
internal helpers on the [Internal functions](@ref) page or beside the
public functions they serve — because the build refuses a docstring
that is not included.

`using SphericalFunctions` brings in the names that most code needs:
the functions `D`, `d`, `sYlm`, and `Ylm` and their calculators, the
containers they return, `ModeWeights`, the differential operators, the
transforms, and the common pixelizations and quadrature weights.  The
rest of the public interface is declared `public` but not exported,
and is reached as `SphericalFunctions.name`, or brought in with
`using SphericalFunctions: name`.  That group holds the accessors,
such as `ℓₘᵢₙ` and `spins`, whose names are too generic to export; the
types that most code never names, such as `HarmonicCalculator`,
`HWedge`, and `HalfOddInteger`; the machinery for defining functions
of indices, `IndexType`, `IndexRange`, `IndexOrRange`, and
`@index_methods`; `sλlmCalculator`, the calculator of the real
harmonics on which the ring-based transforms are built; the
specialized pixelizations, such as `driscoll_healy_pixels` and
`minimal_rings`; a few helpers, such as `wedge_value` and
`map2salm_plan`; and the ASCII aliases of the names that are not
ASCII, such as `ell_min`, `L2`, and `set_beta!`.  Anything else is
internal, as described on the [Internal functions](@ref) page.

```@index
Modules = [SphericalFunctions]
```
