# Calculators

Typically, when calculating special functions, we will use recurrence relations along with some
coefficients — which frequently requires significant setup processing.  That processing can be
cached, so that the calculations themselves consist primarily of memory accesses and simple
arithmetic.  For this reason, we use "calculator" objects, which will be constructed with some
indication of the largest indices you expect to use.  The calculators can then be called repeatedly
for specific values of the arguments, which will compute the function values for all (or some subset
of) indices.

```@autodocs
Modules = [Spherical.AssociatedLegendreFunction]
Pages   = ["src/associated_legendre/calculator.jl"]
```

```@autodocs
Modules = [WignerMatrices]
Pages   = ["wigner_matrices/calculator.jl", "wigner_matrices/evaluate.jl"]
```


# Transformation


```@autodocs
Modules = [Spherical]
Pages   = ["map2salm.jl"]
```
