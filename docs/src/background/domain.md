# [Domain](@id background_domain)

This package deals with standard spherical harmonics, spin-weighted
spherical harmonics, and Wigner's 𝔇 matrices.  The key question is
what domains these functions are defined on — what their arguments
are.  We will discover that it's best to use quaternions for all
three.

Usually, these are written as functions of spherical coordinates for
both types of harmonics, and Euler angles for the 𝔇 matrices:
```math
\begin{gathered}
Y_{ℓ,m}(θ, ϕ), \\
{}_sY_{ℓ,m}(θ, ϕ), \\
𝔇^{(ℓ)}_{m', m}(α, β, γ).
\end{gathered}
```
While that's good enough for many purposes, it obscures the true
nature of these functions, and makes certain manipulations more
difficult or even impossible.  We can make a few observations:

  1. Standard (scalar) spherical harmonics are just a special case of
     spin-weighted spherical harmonics with spin weight ``s=0``.
  2. Spin-weighted spherical harmonics *cannot* be defined as
     functions on the sphere ``𝕊²`` alone; they depend on an
     additional choice of frame at each point on the sphere — a choice
     of basis for the tangent space at that point.  This extra
     structure is necessary to define what "spin weight" means at all.
     (We'll define it [soon](@ref background_differential_operators).)
  3. Euler angles are a bad way to deal with rotations.  Besides the
     coordinate singularities at the poles, they make composing
     rotations far more difficult.  More importantly, they make it
     hard to deal with *generators* of rotations, which are crucial
     for angular-momentum operators.

The second point especially requires a little more explanation.  In
fact, I wrote an entire paper that gets very deep into the details
about this point [Boyle_2016](@cite).  It may seem odd that
spin-weighted spherical functions cannot be defined on the sphere
``𝕊²``, when they're written in terms of spherical coordinates.  The
key subtlety here is that the spherical coordinates ``(θ, ϕ)``
implicitly define a conventional choice of tangent basis.  But spin
weight is *defined* in terms of what happens to a function when you
rotate that tangent basis.  Thus, if we start from spherical
coordinates, we are naturally led to introduce a third angle ``ψ``
describing rotation of the tangent basis about the radial direction.
And those three angles together are just Euler angles.  In fact, ``(θ,
ϕ, ψ)`` is *precisely* the set of "Eulerian" angles defined by the
very influential [book by Whittaker](@ref Whittaker-(1904))
[Whittaker_1947](@cite).  It turns out [GoldbergEtAl_1967](@cite) that
with this extra angle included explicitly, the original spin-weighted
spherical harmonics [Newman_1966](@cite) are simply proportional to
Wigner's 𝔇 matrices,
```math
{}_{s}Y_{ℓ,m}(θ, ϕ, ψ)
= (-1)^s \sqrt{\frac{2ℓ+1}{4\pi}} \, 𝔇^{(ℓ)}_{m, -s}(ϕ, θ, ψ),
```
so we might as well think of them as being parameterized in the same
way.

It might be helpful to pause for a moment and consider a more
intuitive physical picture.  Consider a telescope.  We can point it
any which way we want, which corresponds to choosing a point on the
sphere ``𝕊²``.  Now, imagine an ideal pixel at the very center of the
image formed by that telescope.  The entire sky can be described in
terms of that pixel as a function of where the telescope is pointed.
That is, it's a function of ``𝕊²``.  But now, suppose we put a
polarizer on the telescope.  The pixel's value will now depend not
only on where the telescope is pointed, but also on the orientation of
the polarizer about the optical axis.  This extra degree of freedom is
exactly the extra angle ``ψ`` we mentioned above.  But really, what
we're doing is rotating our telescope in three-dimensional space, so
the pixel's value is a function on the *rotation group*
``\mathrm{SO}(3)``, not just the sphere ``𝕊²``.  Now, with light, we
know every possible polarized value once we measure the value with one
polarization and the value with the orthogonal polarization, which are
combined into a single complex number.  (Compare [Jones
vectors](https://en.wikipedia.org/wiki/Jones_calculus).)  As we rotate
our choice of the first polarization by an angle ``ψ``, the complex
combination varies as ``e^{iψ}``.  This is exactly what it means to
have spin weight ``s=1`` — which is, not coincidentally, the spin of a
photon.  Imagining that we had a polarizing telescope measuring some
other kind of field with spin ``s``, the complex combination would
vary as ``e^{isψ}``.  And since we could have half-integer spins, we
actually need to consider not just the rotation group
``\mathrm{SO}(3)``, but its double cover ``\mathrm{Spin}(3) \cong
\mathrm{SU}(2)``, to fully capture the behavior of general fields.

So, at this point, we see that all three of the functions we care
about — standard spherical harmonics, spin-weighted spherical
harmonics, and Wigner's 𝔇 matrices — are most naturally thought of as
functions on ``\mathrm{Spin}(3)``.  The traditional way to
parameterize this would be with (extended) Euler angles.  However, as
mentioned above, Euler angles are a poor choice for many purposes.
First, we have the problem of coordinate singularities at the "poles"
(aka ["gimbal lock"](https://en.wikipedia.org/wiki/Gimbal_lock)).
This is not a problem if we use quaternions instead.  Second, we have
the problem of composing rotations; given rotations described by Euler
angles ``(α, β, γ)`` and ``(α', β', γ')``, finding another set of
angles ``(α'', β'', γ'')`` representing their composition is a nasty,
nonlinear problem.  With quaternions, composition is just quaternion
multiplication: ``R'' = R'\, R``.  Third, and most importantly, we have
the problem of generators of rotations.  With Euler angles,
infinitesimal rotations are described by complicated derivatives with
bizarre trigonometric factors.  With quaternions, infinitesimal
rotations are described by "pure-vector" quaternions that can be
exponentiated almost as simply as imaginary numbers.  That is, if
``𝐮`` is a unit quaternion with 0 scalar part, and ``θ`` is a real
number, then
```math
\exp \left[ 𝐮 \frac{θ}{2} \right]
= \cos \frac{θ}{2} + 𝐮 \sin \frac{θ}{2}
```
is a unit quaternion representing a rotation by ``θ`` about ``𝐮``.
Note the similarity to Euler's formula for complex exponentials.
Alternatively, going from a unit quaternion to its generator just uses
the logarithm, which is only slightly more complicated than the
complex logarithm.

So the conclusion is clear: we should represent the arguments of all
these functions as unit quaternions, which is what this package does:
```math
\begin{gathered}
Y_{ℓ,m}(𝐐), \\
{}_sY_{ℓ,m}(𝐐), \\
𝔇^{(ℓ)}_{m', m}(𝐐).
\end{gathered}
```
We can still provide convenience methods that convert from spherical
coordinates or Euler angles to quaternions, but internally everything
is done with quaternions, and the documentation will generally be
written in terms of quaternions.  See [Boyle_2016](@citet) for full
details.
