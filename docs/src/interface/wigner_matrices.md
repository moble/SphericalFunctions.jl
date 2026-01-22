# Wigner's ``𝔇`` and ``d`` matrices

```@meta
CurrentModule = SphericalFunctions.Deprecated
```

Wigner's ``𝔇`` matrices — and to a lesser extent, the related ``d`` matrices —
are extremely important in the theory of rotations.  Each element is, itself, a
special function of the rotation group: in particular, an eigenfunction of [the
left- and right-Lie derivatives](@ref "Differential operators"), and
thus a spin-weighted spherical function.  Collectively, they describe how
spin-weighted spherical functions transform under rotation.  But their accurate
and efficient computation is surprisingly subtle.  This package implements the
current state-of-the-art techniques for their fast and accurate computation,
based on the [``H`` recursion](@ref "Algorithm for computing ``H``").

The actual computations can be done with the [`D_matrices`](@ref) function:
```julia
using Quaternionic
using SphericalFunctions
using SphericalFunctions.Deprecated

R = randn(RotorF64)
ℓₘₐₓ = 8
𝔇 = D_matrices(R, ℓₘₐₓ)
```
However, the matrices can take up a lot of memory.  So for maximum efficiency
when calling this function repeatedly with different `R` values, it is best to
pre-allocate the necessary memory with the [`D_prep`](@ref) function, and the
pass that in as an argument to [`D_matrices!`](@ref):
```julia
D_storage = D_prep(ℓₘₐₓ)
𝔇 = D_matrices!(D_storage, R)
```
(Beware that, as noted in the documentation for [`D_matrices!`](@ref), the
output `𝔇` is just a reference to part of the `D_storage` object, so you should
not reuse `D_storage` until you have copied or otherwise finished using `𝔇`.)

The output `𝔇` is a (linear!) vector of `Complex` numbers with the same base
type as `R`.  The ordering of the elements is described in the documentation for
[`D_matrices`](@ref).  It is also possible to efficiently view slices of this
vector as a series of individual matrices using a [`D_iterator`](@ref):
```julia
for (ℓ, Dˡ) in zip(0:ℓₘₐₓ, D_iterator(𝔇, ℓₘₐₓ))
    # Do something with the matrix Dˡ[ℓ+m′+1, ℓ+m+1]
end
```

For the ``d`` matrices, we have almost the same interface, except that instead of
the input quaternion `R` we only need the angle `β` (or its complex angle `expiβ`,
which can be computed directly in some cases), and the output is real-valued:
```julia
using Quaternionic
using SphericalFunctions
using SphericalFunctions.Deprecated

β = π * rand(Float64)
ℓₘₐₓ = 8
d = d_matrices(β, ℓₘₐₓ)
```
Again, for repeated calls, it is best to pre-allocate storage:
```julia
d_storage = d_prep(ℓₘₐₓ)
d = d_matrices!(d_storage, β)
```
The output `d` is a vector of numbers of the same type as `β`, ordered in the
same way as the output of [`D_matrices`](@ref).  And similarly, we can iterate
over the individual matrices using a [`d_iterator`](@ref).


## Docstrings

```@docs
D_matrices
D_matrices!
D_prep
D_iterator
d_matrices
d_matrices!
d_prep
d_iterator
```


## Workspaces

```@meta
CurrentModule = SphericalFunctions
```

```@docs
HWedge
HAxis
```
