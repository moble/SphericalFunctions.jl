# Spherical Functions

Julia package for evaluating and transforming Wigner's 𝔇 matrices, Wigner's 3-j symbols, and
spin-weighted (and scalar) spherical harmonics.  These functions are evaluated directly in terms of
quaternions, as well as in the more standard forms of spherical coordinates and Euler
angles.<sup>[1](#1-euler-angles-are-awful)</sup>

These quantities are computed using recursion relations, which makes it possible to compute to very
high ℓ values.  Unlike direct evaluation of individual elements, which will generally cause overflow
or underflow beyond ℓ≈30 when using double precision, these recursion relations should be valid for
far higher ℓ values.  More precisely, `Inf` values appear starting at ℓ=22 for `Float16`, ℓ=183 for
`Float32`, and ℓ=1474 for `Float64`.  `BigFloat` also works, and presumably will not overflow for
any ℓ value that could reasonably fit into computer memory — though it is far slower.  Also note
that [`DoubleFloats`](https://github.com/JuliaMath/DoubleFloats.jl) will work, and achieve
significantly greater accuracy (but no greater range) than `Float64`.  The results are accurate to
roughly ℓ times the precision of the input quaternion.

The conventions for this package are described in detail on [this
page](http://moble.github.io/spherical/).

## Installation

```bash
julia -e 'using Pkg; pkg"add https://github.com/moble/Spherical.jl.git"'
```

## References

The most important routine in this package is the computation of the 𝔇 matrices — or more
specifically, of terms proportional to parts of the 𝔇 matrices.  This mostly follows the treatment
of [Gumerov and Duraiswami (2014)](https://arxiv.org/abs/1403.7698).  To seed the recursions they
present, we also need to calculate the associated Legendre functions.  Currently, this is done using
the "modified forward row method" of [Holmes and Featherstone
(2002)](https://doi.org/10.1007/s00190-002-0216-2).  Note that this is apparently the source of
overflow noted above.  Two other (though more complicated) methods have appeared more recently in
the literature, which could presumably extend these limits much further.  [Fukushima
(2012)](https://doi.org/10.1007/s00190-011-0519-2) showed that using "X-numbers" (wherein the
exponent is stored as a separate integer) in the core of the recursion could increase the range to
ℓ≈2³².  [Xing et al. (2020)](https://doi.org/10.1007/s00190-019-01331-0) showed that Fukushima's
results exhibited increased error for certain angles, whereas their Eqs. (12)—(14) could be used
directly to obtain results with greater accuracy for those certain angles, and comparable accuracy
for other angles.

It may be worthwhile simply implementing X-numbers as a subtype of `AbstractFloat`, and simply
passing them to these algorithms if higher orders are needed.  (Though using `BigFloat` would
achieve a similar objective, it would probably be far slower.)  The actual recommendation of
Fukushima is more sophisticated — just using X-numbers in the core calculation — but it looks like
this probably wouldn't be *too* much slower.


## Usage

#### Functions of angles or rotations

Currently, due to the nature of recursions, this module does not allow
calculation of individual elements, but returns ranges of results.  For
example, when computing Wigner's 𝔇 matrix, all matrices up to a given ℓ will be
returned; when evaluating a spin-weighted spherical harmonic, all harmonics up
to a given ℓ will be returned.  Fortunately, this is usually what is required
in any case.

To calculate Wigner's d or 𝔇 matrix or spin-weighted spherical harmonics, first
construct a `Wigner` object.

```python
using Spherical
ell_max = 16  # Use the largest ℓ value you expect to need
wigner = Wigner(ell_max)
```

This module takes input as quaternions.  The `Quaternionic` module has [various
ways of constructing
quaternions](https://quaternionic.readthedocs.io/en/latest/#rotations),
including direct construction or conversion from rotation matrices, axis-angle
representation, Euler angles,<sup>[1](#euler-angles-are-awful)</sup> or
spherical coordinates, among others:

```python
R = quaternionic.array([1, 2, 3, 4]).normalized
R = quaternionic.array.from_axis_angle(vec)
R = quaternionic.array.from_euler_angles(alpha, beta, gamma)
R = quaternionic.array.from_spherical_coordinates(theta, phi)
```

Mode weights can be rotated as

```python
wigner.rotate(modes, R)
```

or evaluated as

```python
wigner.evaluate(modes, R)
```

We can compute the 𝔇 matrix as

```python
D = wigner.D(R)
```

which can be indexed as

```python
D[wigner.Dindex(ell, mp, m)]
```

or we can compute the spin-weighted spherical harmonics as

```python
Y = wigner.sYlm(s, R)
```

which can be indexed as

```python
Y[wigner.Yindex(ell, m)]
```

Note that, if relevant, it is probably more efficient to use the `rotate` and
`evaluate` methods than to use `D` or `Y`.

The packages
[`FastTransforms.jl`](https://JuliaApproximation.github.io/JuliaApproximation/FastTransforms.jl/)
and
[`FastSphericalHarmonics.jl`](https://eschnett.github.io/FastSphericalHarmonics.jl/)
(which is largely a wrapper for the former) also have many of these functions —
possibly in more efficient forms — for `Float64` types.



#### Clebsch-Gordan and 3-j symbols

It is possible to compute individual values of the 3-j or Clebsch-Gordan
symbols:

```python
w3j = spherical.Wigner3j(j_1, j_2, j_3, m_1, m_2, m_3)
cg = spherical.clebsch_gordan(j_1, m_1, j_2, m_2, j_3, m_3)
```

However, when more than one element is needed (as is typically the case), it is
much more efficient to compute a range of values:

```python
calc3j = spherical.Wigner3jCalculator(j2_max, j3_max)
w3j = calc3j.calculate(j2, j3, m2, m3)
```

Note that [this package](https://github.com/Jutho/WignerSymbols.jl) seems to do
a clever job of evaluating these quantities in terms of [rational
roots](https://github.com/Jutho/RationalRoots.jl).  I haven't used or examined
the package otherwise, so I can't really vouch for it, but this seems like a
clever idea.


## Acknowledgments

I very much appreciate Barry Wardell's help in sorting out the relationships
between my conventions and those of other people and software packages
(especially Mathematica's crazy conventions).

This code is, of course, hosted on github.  Because it is an open-source
project, the hosting is free, and all the wonderful features of github are
available, including free wiki space and web page hosting, pull requests, a
nice interface to the git logs, etc.

The work of creating this code was supported in part by the Sherman Fairchild
Foundation and by NSF Grants No. PHY-1306125 and AST-1333129.


<br/>

---

###### <sup>1</sup> Euler angles are awful

Euler angles are pretty much
[the worst things ever](http://moble.github.io/spherical/#euler-angles)
and it makes me feel bad even supporting them.  Quaternions are
faster, more accurate, basically free of singularities, more
intuitive, and generally easier to understand.  You can work entirely
without Euler angles (I certainly do).  You absolutely never need
them.  But if you're so old fashioned that you really can't give them
up, they are fully supported.
