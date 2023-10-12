# ``{}_{s}Y_{\ell,m}`` functions

The spin-weighted spherical harmonics are an [important set of functions defined
on](@cite Boyle_2016) the rotation group ``𝐒𝐎(3)``, or more generally, the
spin group ``𝐒𝐩𝐢𝐧(3)`` that covers it.  They are eigenfunctions of [the
left- and right-Lie derivatives](@ref "Differential operators"), and are
particularly useful in describing the angular dependence of polarized fields,
like the electromagnetic field and gravitational-wave field.  Originally
introduced by [Newman_1966](@citet), they are essentially components of Wigner's
``\frak{D}`` matrices:
```math
{}_{s}Y_{\ell,m}(\mathbf{R})
  = (-1)^s \sqrt{\frac{2\ell+1}{4\pi}} \, \frak{D}^{(\ell)}_{m, -s}(\mathbf{R}).
```
As such, they can be computed using the same [``H`` recursion](@ref "Algorithm
for computing ``H``") algorithm as the Wigner ``\frak{D}^{(\ell)}_{m, -s}``
matrices.  But because not all values of ``s \in -\ell:\ell`` are used, we can
be much more efficient in both storage and computation time.

The user interface is very similar to the one for [Wigner's ``𝔇`` and ``d``
matrices](@ref):
```julia
using Quaternionic
using SphericalFunctions

R = randn(RotorF64)
ℓₘₐₓ = 8
s = -2
Y = sYlm_values(R, ℓₘₐₓ, s)
```
Again, the results can take up a lot of memory, so for maximum efficiency when
calling this function repeatedly with different `R` values, it is best to
pre-allocate the necessary memory with the [`sYlm_prep`](@ref) function, and the
pass that in as an argument to [`sYlm_values!`](@ref):
```julia
Y_storage = sYlm_prep(ℓₘₐₓ, s)
Y = sYlm_values!(Y_storage, R, s)
```
(Beware that, as noted in the documentation for [`sYlm_values!`](@ref), the
output `Y` is just a reference to part of the `Y_storage` object, so you should
not reuse `Y_storage` until you have copied or otherwise finished using `Y`.)

The output `Y` is a single vector of `Complex` numbers with the same base type
as `R`.  The ordering of the elements is described in the documentation for
[`sYlm_values!`](@ref).  It is also possible to efficiently view slices of this
vector as a series of individual vectors using a [`sYlm_iterator`](@ref):
```julia
for (ℓ, Yˡ) in zip(0:ℓₘₐₓ, sYlm_iterator(Y, ℓₘₐₓ))
    # Do something with the matrix Yˡ[ℓ+m′+1, ℓ+m+1]
end
```

## Docstrings

```@docs
sYlm_values
sYlm_values!
sYlm_prep
sYlm_iterator
```
