# References

The most important routine in this package is the computation of the 𝔇 matrices — or more
specifically, of terms proportional to parts of the 𝔇 matrices.  This mostly follows the treatment
of [Gumerov_2015](@citet) (with minor modifications to account for errors in their presentation, as
described [here](/SphericalFunctions/notes/H_recursions/#Steps-to-compute-H)).
To seed the recursions they present, we also need to calculate the associated Legendre functions.
This is now done using the "fully normalized column-wise recurrence formula" (fnCWF) given by
Eqs. (12)—(14) of [Xing_2019](@citet).  This improves significantly over the older implementation
using the "modified forward row method" of [Holmes_2002](@citet), for which the results would fail
to be finite starting at ℓ=22 for `Float16`, ℓ=183 for `Float32`, and ℓ=1474 for `Float64`.  Another
approach that was never precisely implemented in this package was due to [Fukushima_2011](@citet),
who showed that using "X-numbers", wherein the exponent is stored as a separate integer,
(implemented in [this package](https://github.com/moble/XNumbers.jl)) in the core of the recursion
could increase the range to ℓ≈2³².  Xing et al. showed that Fukushima's results exhibited increased
error for certain angles, whereas their Eqs. (12)—(14) could be used directly to obtain results with
greater accuracy for those certain angles, and comparable accuracy for other angles.

The other major functionality of this package is `map2salm` / `salm2map`, which decomposes function
values on regular grids into mode weights (coefficients), and vice versa.  The approach used here is
taken from [Reinecke_2013](@citet), with weights based on the method by [Waldvogel_2006](@citet).
However, this interface has been superseded by the `SSHT` object, which implements several
approaches, including the Reinecke-Seljebotn-Waldvogel method, as well as the optimal-dimensionality
method due to [Elahi_2018](@citet), as well as a new unpublished optimal-dimensionality method I
(Mike Boyle) created.

# Bibliography

```@bibliography
*
```

