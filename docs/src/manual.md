# Function calculators

Typically, when calculating special functions, we will use recursion relations along with some
coefficients โ which frequently requires significant setup processing.  That processing can be
cached, so that the calculations themselves consist primarily of memory accesses and simple
arithmetic.

Similarly, some computations require a certain amount of "workspace" โ significant amounts of memory
that are needed for the computation, but are not part of the core result (and don't even need to be
returned).  Rather than allocating this memory on each call to the computation, we can allocate once
and pass the workspace in as an argument.

These two patterns โ pre-computed coefficients and pre-allocated workspaces โ will appear as
possible call signatures in each of the following functions.  Note that we also have helper
functions like `Dprep`, `Dstorage`, etc., which construct the necessary coefficients and/or
workspaces, as well as the actual results for in-place computations, to make it easier to call the
corresponding calculator functions.

Note that these functions impose assumptions about the ordering of arrays as vectors.  For example,
the series of Wigner ``๐`` matrices will be stored as a single vector.  And for particular values of
``(โ, n, m)``, the element ``๐หกโ,โ`` will have to be accessed as some linearly indexed element of
that vector.  These formats, along with utility functions for accessing particular elements, are
described [here](/utilities/#Y-and-D-data).

The fundamental algorithm is the ``H`` recursion, which is the core computation needed for Wigner's
``d`` and ``๐`` matrices, and the spin-weighted spherical harmonics ``{}_{s}Y_{\ell,m}``, as well as
`map2salm` functions.

```@autodocs
Modules = [SphericalFunctions]
Pages   = ["wigner_matrices/Hrecursor.jl"]
```

Internally, the ``H`` recursion relies on calculation of the Associated Legendre Functions (ALFs),
which can also be called on their own:

```@autodocs
Modules = [SphericalFunctions]
Pages   = ["associated_legendre.jl"]
```

Based on those, we have the `d!`, `D!`, and `Y!` functions:

```@autodocs
Modules = [SphericalFunctions]
Pages   = ["wigner_matrices/evaluate.jl"]
```

# Transformation

The functions above โ especially `Y!` โ allow us to evaluate mode weights (often denoted `salm`) of
a spin-weighted spherical-harmonic expansion as a function of a particular point (often denoted
`map`).  This process is frequently referred to as "synthesis", and is also called `salm2map`.

We can also perform the reverse process of going from the function evaluated on some kind of grid,
back to the corresponding mode weights โ also called "analysis" or `map2salm`.

```@autodocs
Modules = [SphericalFunctions]
Pages   = ["map2salm.jl"]
```
