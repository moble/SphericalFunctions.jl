# ``s``-SHT Transformations

Any square-integrable function on the sphere 𝕊² or 𝕊³ can be
represented as an expansion in spherical harmonics or spin-weighted
spherical harmonics, respectively.  For a particular spin-weight
``s``, we can restrict the spin-weighted spherical harmonics to 𝕊²
(because its behavior on the remaining [𝕊¹ factor of
𝕊³](https://en.wikipedia.org/wiki/Hopf_fibration) is determined by
the spin).  So, for example, if ``f`` is a function with spin weight
``s``, we will frequently write
```math
\begin{gathered}
f(θ, ϕ) = \sum_{ℓ = |s|}^{ℓₘₐₓ} \sum_{m = -ℓ}^{ℓ} f̃_{ℓ, m}
    \, {}_{s}Y_{ℓ, m}(θ, ϕ) \\
f(𝐑) = \sum_{ℓ = |s|}^{ℓₘₐₓ} \sum_{m = -ℓ}^{ℓ} f̃_{ℓ, m}
    \, {}_{s}Y_{ℓ, m}(𝐑),
\end{gathered}
```
where ``(θ, ϕ)`` are spherical coordinates on 𝕊², and we use the
rotor ``𝐑`` to describe a point on 𝕊³.  The upper limit ``ℓₘₐₓ`` is
— in principle — infinite, but because we are finite ``ℓₘₐₓ`` will
also be finite in all applications here.  Similarly, the set of points
on which we evaluate the function will also be finite.  A little
terminology will be helpful:

  1. Values ``f`` of the function evaluated on a set of points or
     "pixels" in the domain of the function — sometimes called "nodes"
     or the "map" values.
  2. Values ``f̃`` of the mode weights (coefficients) of the expansion
     — usually called "modes" or sometimes "salm" (for ``_sa_{ℓ,m}``).

Here, we are concerned with transformations between these two
representations of the function.  In the literature, the
transformation ``f ↦ f̃`` is usually called "analysis" or `map2salm`,
while the inverse transformation ``f̃ ↦ f`` is called "synthesis" or
`salm2map`.  These are both referred to as spin-spherical-harmonic
transforms, or ``s``-SHTs.

## Synthesis

The expressions written above already show us one way to transform
*from* modes ``f̃`` *to* function values ``f``.  If the pixels on
which we evaluate are indexed by ``k``, then we can write
```math
f(𝐑_k) = \sum_{ℓ,m} f̃_{ℓ, m}\, {}_{s}Y_{ℓ, m}(𝐑_k).
```
Considering the ``(ℓ,m)`` pairs as a single index, this is just a
matrix-vector multiplication, where the matrix elements are given by
```math
𝒯_{k, (ℓ,m)} = {}_{s}Y_{ℓ, m}(𝐑_k),
```
and we might more compactly write
```math
f = 𝒯 \, f̃.
```
In this expression, ``f̃`` is being treated as a column vector of mode
weights, and ``f`` is being treated as a column vector of function
values at the sequence of pixels.  For instance, we will frequently
store the values ``f̃`` as the (column) vector
```julia
f̃ = [mode_weight(ℓ, m) for ℓ ∈ abs(s):ℓₘₐₓ for m ∈ -ℓ:ℓ]
```
(Here, `mode_weight` is just for illustration.)  If this
transformation will be performed repeatedly, it can be very efficient
to pre-compute the matrix ``𝒯``, and capitalize on the impressive
efficiency of linear-algebra libraries to perform the transformation.
And indeed, this is the approach taken by the [`SSHTDirect`](@ref)
method.

However, we must consider the memory requirements of this approach.
The number of nonzero mode weights is ``M = (ℓₘₐₓ+1)^2 - s^2``.  If
there are ``N`` pixels in the pixelization, then the matrix ``𝒯`` has
size ``N × M``.  Typically, we will use roughly the same number of
pixels as there are mode weights to be able to express roughly the
same number of degrees of freedom, so that ``N ≈ M``.  Therefore, the
size of this matrix will typically be ``N × M ∼ ℓₘₐₓ^4``.  For
standard complex numbers stored in 128 bits, this will exceed 1 GiB of
memory for ``ℓₘₐₓ ≳ 90``, and grow rapidly.  The computational cost of
the matrix-vector multiplication will scale just as poorly.  This
full-storage-matrix approach is very attractive, but only for
relatively small ``ℓₘₐₓ``.

Another approach would be to compute the matrix elements as needed.
Given that the SWSH values are best computed via recurrence relations,
it would make sense to compute them row-by-row.  This would
essentially eliminate the memory requirements.  Though the
computational cost would still scale poorly, it would parallelize very
well, since each row is independent until a final summation.  This is
not currently implemented, but the pieces are all in place, and it
will be implemented in the future.

However, we can achieve significantly better performance at high
``ℓₘₐₓ`` if we select the pixelization carefully.  Recall that a
spherical harmonic of *any* spin weight still varies azimuthally as
``e^{i m ϕ}``.  Therefore, if we select a pixelization made up of a
series of rings, each of which is at a constant ``θ`` and has pixels
equally spaced in ``ϕ``, then we can use Fast Fourier Transforms
(FFTs) to perform the azimuthal integration very quickly.




# Analysis

Analytically, we use orthogonality of the spin-weighted spherical
harmonics to compute the mode weights from the function values as
```math
\begin{gathered}
f̃_{ℓ, m} = \int_{𝕊²} f(θ, ϕ)\, {}_{s}Ȳ_{ℓ, m}(θ, ϕ) \, d^2Ω, \\
f̃_{ℓ, m} = \int_{𝕊³} f(𝐑)\, {}_{s}Ȳ_{ℓ, m}(𝐑) \, d^3Ω.
\end{gathered}
```
But — again because we are finite — we will only be able to evaluate
the function at a finite number of points, and so we will need to use
discrete quadrature to evaluate the integrals.



To describe the mode weights of a spin-``s`` function up to (and
including) some maximum angular resolution ``\ell_\mathrm{max}``,
there are ``(\ell_\mathrm{max}+1)^2 - s^2`` mode weights.  We assume
throughout that the values `f̃` are stored as the (column) vector
```julia
f̃ = [mode_weight(ℓ, m) for ℓ ∈ abs(s):ℓₘₐₓ for m ∈ -ℓ:ℓ]
```
(Here, `mode_weight` is a made-up function for schematic purposes.) In
particular, the ``m`` index varies most rapidly, and the ``\ell``
index varies most slowly.  Correspondingly, there must be *at least*
``(\ell_\mathrm{max}+1)^2 - s^2`` function values `f`.  However, some
``s``-SHT algorithms require more function values — usually by a
factor of 2 or 4 — trading off between speed and memory usage.

The `SSHT` object implements these transformations, storing
pre-computed constants and pre-allocated workspace for the
transformations.  The interface is designed to be similar to that of
`AbstractFFTs.jl`, whereby an `SSHT` object `𝒯` can be used to
perform the transformation as either
```julia
f = 𝒯 * f̃
```
or
```julia
f̃ = 𝒯 \ f
```
Currently, there are three algorithms implemented, each having different
advantages and disadvantages:

  1. The "Direct" algorithm (introduced here for the first time), which is the
     default, but should only be used up to ``\ell_\mathrm{max} \lesssim 50``
     because its intermediate storage requirements scale as
     ``\ell_\mathrm{max}^4``.  This algorithm is the fastest for small
     ``\ell_\mathrm{max}``, it can be used with arbitrary (non-degenerate)
     pixelizations, and achieves optimal dimensionality.
  2. The "Minimal" algorithm due to [Elahi_2018](@citet), with some minor
     improvements.  This algorithm is fast and — as the name implies — also
     achieves optimal dimensionality, and its storage scales as
     ``\ell_\mathrm{max}^3``.  However, its pixelization is restricted, and its
     accuracy at very high ``\ell_\mathrm{max}`` is not as good as the "RS"
     algorithm.  The algorithm itself is not actually fully specified by Elahi
     et al., and leaves out some relatively simple improvements, so I have had
     to take some liberties with my interpretation.
  3. The "RS" algorithm due to [Reinecke_2013](@citet).  This forms the basis
     for the [`libsharp`](https://gitlab.mpcdf.mpg.de/mtr/libsharp) and
     [`ducc.sht`](https://gitlab.mpcdf.mpg.de/mtr/ducc#duccsht) packages.  It
     requires pixelizations on "iso-latitude rings", and does not achieve
     optimal dimensionality.  However, it is very fast, and its accuracy is
     excellent at extremely high ``\ell_\mathrm{max}``.


## `SSHT` objects

```@autodocs
Modules = [SphericalFunctions]
Pages   = ["ssht.jl", "ssht/direct.jl", "ssht/minimal.jl", "ssht/rs.jl"]
```

## Pixelizations

The algorithms implemented here require pixelizations.  While the
"Direct" algorithm can be used with arbitrary pixelizations, the
"Minimal" and "RS" algorithms require more specific choices, as noted
in their docstrings.

Typically, "pixelization" refers exclusively to a choice of points on
the sphere 𝕊² at which to compute function values.  Of course, as
mentioned [elsewhere](@cite Boyle_2016), it is not *technically
possible* to define spin-weighted functions as functions of a point on
𝕊² alone; we also need some sense of reference direction in the
tangent space.  Quite generally, we can define spin-weighted functions
on the group 𝐒𝐎(3) or 𝐒𝐩𝐢𝐧(3), so we will also refer to a choice
of a set of points in 𝐒𝐩𝐢𝐧(3) (which is essentially the group of
unit quaternions) as a "pixelization".  However, assuming spherical
coordinates, a choice of *coordinates* on the sphere almost everywhere
induces a choice of the reference direction in the tangent space, so
it is *almost* possible to define pixelizations just in terms of
points on 𝕊².  But using spherical coordinates is actually enough to
fully specify the pixelization, because the degeneracies at the poles
also allow us to define the reference direction.

In principle, we could be concerned about the choice of reference
direction in the tangent space.  That is, we might expect to care
about pixelizations over 𝕊³.  However, we are dealing with
spin-weighted functions, which are eigenfunctions of a final rotation
about the reference direction.  This means that once we choose any
reference direction at each point, we know the function values for any
other reference direction at those points.  In particular, an
important property of a pixelization is the condition number of the
transformation matrix between the function values and the mode
weights.  If we rotate the reference direction at a single point, this
is equivalent to multiplying the matrix by a diagonal matrix with
entries of 1 everywhere except the entry corresponding to that point,
where the entry is some complex phase.  This does not change the
condition number of the matrix, so we can ignore the choice of
reference direction at every point.  For other situations, where we
might care about the choice of reference direction, it might be
interesting to consider [this work by Marc
Alexa](https://github.com/moble/superfibonacci), and references
therein.

Interesting discussions of various pixelizations and metrics can be
found in [Saff and Kuijlaars (1997)](@cite SaffKuijlaars_1997) and
[Brauchart and Grabner (2015)](@cite BrauchartGrabner_2015), as well
as blog posts
[here](https://web.archive.org/web/20220303150307/https://www.maths.unsw.edu.au/about/distributing-points-sphere)
and
[here](https://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/).
Note that the "equal-area" pixelizations of
[Healpix](https://healpix.sourceforge.io/) are very restrictive—only
being available for very specific numbers of points—and do not provide
any obvious advantages over the more flexible pixelizations available
here.

The various pixelizations may be computed as follows:

```@autodocs
Modules = [SphericalFunctions]
Pages   = ["pixelizations.jl"]
```


## Quadrature weights

The "RS" algorithm requires quadrature weights corresponding to the input
pixelization.  Though there is a working default choice, it is possible to use
others.  There are several that are currently implemented, along with their
corresponding pixelizations:

```@autodocs
Modules = [SphericalFunctions]
Pages   = ["weights.jl"]
Order   = [:module, :type, :constant, :function, :macro]
```
