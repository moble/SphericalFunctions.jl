# ``s``-SHT Transformations

One important capability of this package is the transformation between the two
representations of a spin-weighted spherical function:

  1. Values `f` of the function evaluated on a set of points or "pixels" in the
     domain of the function.
  2. Values `f̃` of the mode weights (coefficients) of an expansion in the
     standard spin-weighted spherical-harmonic basis.

In the literature, the transformation `f` ↦ `f̃` is usually called "analysis" or
`map2salm`, while the inverse transformation `f` ↦ `f̃` is called "synthesis" or
`salm2map`.  These are both referred to as spin-spherical-harmonic transforms,
or ``s``-SHTs.

To describe the values of a spin-``s`` function up to some maximum angular
resolution ``\ell_\mathrm{max}``, we need ``(\ell_\mathrm{max}+1)^2 - s^2`` mode
weights.  We assume throughout that the values `f̃` are stored as
```julia
f̃ = [mode_weight(ℓ, m) for ℓ ∈ abs(s):ℓₘₐₓ for m ∈ -ℓ:ℓ]
```
(Here, `mode_weight` is a made-up function intended to provide a schematic.)  In
particular, the ``m`` index varies most rapidly, and the ``\ell`` index varies
most slowly.  Correspondingly, there must be *at least*
``(\ell_\mathrm{max}+1)^2 - s^2`` function values `f`.  However, some ``s``-SHT
algorithms require more function values — usually by a factor of 2 or 4 —
trading off between speed and memory usage.

The `SSHT` object implements these transformations, storing pre-computed
constants and pre-allocated workspace for the transformations.  The interface is
designed to be similar to that of `FFTW.jl`, whereby an `SSHT` object `𝒯` can
be used to perform the transformation as either
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

The algorithms implemented here require pixelizations.  While the "Direct"
algorithm can be used with arbitrary pixelizations, the "Minimal" and "RS"
algorithms require more specific choices, as noted in their docstrings.

Typically, "pixelization" refers exclusively to a choice of points on the sphere
𝕊² at which to compute function values.  Of course, as mentioned
[elsewhere](@cite Boyle_2016), it is not *technically possible* to define
spin-weighted functions as functions of a point on 𝕊² alone; we also need some
sense of reference direction in the tangent space.  Quite generally, we can
define spin-weighted functions on the group 𝐒𝐎(3) or 𝐒𝐩𝐢𝐧(3), so we will
also refer to a choice of a set of points in 𝐒𝐩𝐢𝐧(3) (which is essentially
the group of unit quaternions) as a "pixelization".  However, assuming spherical
coordinates, a choice of *coordinates* on the sphere almost everywhere induces a
choice of the reference direction in the tangent space, so it is *almost*
possible to define pixelizations just in terms of points on 𝕊².

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
