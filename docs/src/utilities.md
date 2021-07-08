# Utilities

While not usually the star of the show, the following utilities can be quite helpful for actually
using the rest of the code.


## Complex powers

One common task we find when working with spherical functions is the computation of a range of
integer powers of some complex number — so much so that it can be best to pre-compute the powers and
cache their values.  While a naive approach is generally quite accurate, and reasonably fast, we can
do a little better with a specialized routine.

```@autodocs
Modules = [Spherical]
Pages   = ["complex_powers.jl"]
Order   = [:module, :type, :constant, :function, :macro]
```


## Quadrature weights

```@autodocs
Modules = [Spherical]
Pages   = ["weights.jl"]
Order   = [:module, :type, :constant, :function, :macro]
```


## ``Y`` data

By ``Y`` data, we mean anything indexed like ``Y_{\ell, m}`` modes.

```@autodocs
Modules = [Spherical]
Pages   = ["indexing.jl"]
Order   = [:module, :type, :constant, :function, :macro]
```


## ``D`` data

By ``D`` data, we mean anything indexed like Wigner's ``\mathfrak{D}`` matrices, or special subsets
of them, like the ``H`` matrices.


```@autodocs
Modules = [Spherical]
Pages   = ["wigner_matrices/indexing.jl", "wigner_matrices/Hrecursions.jl"]
Order   = [:module, :type, :constant, :function, :macro]
```


