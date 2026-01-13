# [Domain](@id background_domain)

This package deals with standard spherical harmonics, spin-weighted
spherical harmonics, and Wigner's 𝔇 matrices.  The key question is
what domains these functions are defined on — what their arguments
are.  We will discover that it's best to use quaternions for all
three.

## Three functions on two domains

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
= (-1)^s \sqrt{\frac{2ℓ+1}{4π}} \, 𝔇^{(ℓ)}_{m, -s}(ϕ, θ, ψ),
```
so we might as well think of them as being parameterized in the same
way.

## The geometric picture

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
exactly the extra angle ``ψ`` we mentioned above.

But really, what we're doing is rotating our telescope in
three-dimensional space, so the pixel's value is a function on the
*rotation group* ``\mathrm{SO}(3)``, not just the sphere ``𝕊²``.
Now, with light, we know every possible polarized value once we
measure the value with one polarization and the value with the
orthogonal polarization, which are combined into a single complex
number.  (Compare [Jones
vectors](https://en.wikipedia.org/wiki/Jones_calculus).)  As we rotate
our choice of the first polarization by an angle ``ψ``, the complex
combination varies as ``e^{iψ}``.  This is exactly what it means to
have spin weight ``s=1`` — which is, not coincidentally, the spin of a
photon.

Imagining that we had a polarizing telescope measuring some other kind
of field with spin ``s``, the complex combination would vary as
``e^{isψ}``.  (This is important for gravitational-wave astronomy with
spin-2 fields.  In principle, we could also consider neutrino
telescopes with spin-1/2 polarization.)  And since we could have
half-integer spins, we actually need to consider not just the rotation
group ``\mathrm{SO}(3)``, but its double cover ``\mathrm{Spin}(3)
\cong \mathrm{SU}(2)``, to fully capture the behavior of general
fields.

## Unification in ``\mathrm{Spin}(3)``

So, at this point, we see that all three of the functions we care
about — standard spherical harmonics, spin-weighted spherical
harmonics, and Wigner's 𝔇 matrices — are naturally thought of as
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
multiplication: ``R'' = R'\, R``.  Third, and most importantly, we
have the problem of generators of rotations.  With Euler angles,
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

## Pushing forward to ``𝕊²``

Assuming the functions have been defined on ``\mathrm{Spin}(3)``, we
can push them forward to functions on the sphere ``𝕊²``.[^1]  Given a
choice of a special point in ``𝕊²`` — conventionally the north pole
``𝐳`` — we can map ``𝐐 ∈ \mathrm{Spin}(3)`` to ``𝕊²`` simply by
using it to rotate ``𝐳`` to ``π(𝐐) = 𝐐\, 𝐳\, 𝐐⁻¹``.  For any
particular point ``𝐧 ∈ 𝕊²``, the set of all rotors that map to that
point (its preimage) is of the form
```math
π^{-1}(𝐧) = \left\{𝐐\, e^{𝐤 θ/2} \,\middle|\, θ ∈ [0, 4π)\right\},
```
where ``𝐐`` is any particular rotor that maps ``𝐳`` to ``𝐧``, and
``𝐤`` is the generator of rotations about the ``z`` axis.  In terms
of the telescope analogy presented above, ``θ`` represents rotation of
the telescope polarizer about the optical axis.

[^1]: Okay, technically we'll use more structure than just ``𝕊²``.
    We're picking out a special point in that space, and assuming the
    action of ``\mathrm{Spin}(3)`` on that point to define the
    mapping.  Elsewhere, I criticize the usual approach because it
    ignores the fundamental importance of the choice of tangent basis
    at each point.  This pushforward to ``𝕊²`` also uses some extra
    structure — perhaps reinforcing the notion that spin-weighted
    functions really should just be thought of as functions on
    ``\mathrm{Spin}(3)``.

Now, for any function ``f(𝐐)`` on ``\mathrm{Spin}(3)``, we can define
the pushforward function on ``𝕊²`` by taking a point to the average
value of ``f`` over the preimage.  (Such "integration along the fiber"
is a standard concept in differential geometry [BottTu_1982](@cite).)
By abuse of notation, we will use the same symbol ``f`` for this new
function, with the understanding that the symbol will be reinterpreted
as needed.  Thus, for ``𝐧 ∈ 𝕊²``, we define
```math
f(𝐧) = \frac{1}{4π} \int_0^{4π} f\left(𝐐\, e^{𝐤 θ/2}\right) \, dθ,
```
where ``𝐐`` is any rotor such that ``π(𝐐) = 𝐧``.  The choice of
``𝐐`` does not matter because the integral averages over all possible
choices.  Also, given the behavior ``e^{isψ}`` described above, we
know that only fields with ``s=0`` will have nonzero values under this
operation.  Thus, only functions with spin weight 0 can be pushed
forward to nontrivial functions on ``𝕊²``.  And again, because the
standard scalar spherical harmonics are precisely the spin-weighted
spherical harmonics with ``s=0``, they can be pushed forward in this
way.

This is effectively the only continuous way to push functions forward
to *all of* ``𝕊²``.  For functions with nonzero spin-weight, we need
a reference direction to specify the tangent basis.  However, the
[hairy-ball theorem](https://en.wikipedia.org/wiki/Hairy_ball_theorem)
tells us this cannot be done continuously over the entire sphere.
Originally, the spin-weighted spherical harmonics were not defined on
``𝕊²``, but the spherical coordinates — which are topologically the
cylinder ``I×𝕊¹``.

## Pulling back to ``I×𝕊¹``

Finally, we can come back to spherical coordinates: the coordinates
``(θ, ϕ)``, with ``θ ∈ [0, π]`` and ``ϕ ∈ [0, 2π)``.  Note that this
range for ``ϕ`` is half open, only because it is parameterizing a
circle; ``0`` and ``2π`` represent the same point.  Thus, we see that
the spherical coordinates do not actually parameterize a sphere, but
rather a cylinder.  Of course, that cylinder is usually mapped onto
the sphere, shrinking the top and bottom edges to points which
represent the poles of the sphere and therefore the singularities of
the spherical coordinates.  Nonetheless, a function defined on
spherical coordinates is really a function on the cylinder ``I×𝕊¹``.

When we write *spin-weighted* functions in terms of spherical
coordinates, we are not only quietly using the "wrong" domain, but
also implicitly choosing a tangent basis at each point.  Specifically,
the ``\boldsymbol{θ}`` unit vector points directly down the cylinder,
mapping to a vector field that everywhere points toward the south pole
on the sphere.  This mapping fails to give a unique tangent vector at
both poles, of course, which is consistent with the hairy-ball
theorem.  Still, it gives us a well defined function *almost
everywhere* on the sphere.  Interestingly, the function is well
defined *everywhere* on the cylinder.  These two considerations — the
actual topology and the implicit choice of tangent basis — show that
it really is more correct to think of spin-weighted spherical
functions written in spherical coordinates as being defined on the
cylinder ``I×𝕊¹``.

Now, interestingly, once we've picked out a basis, we can actually
define the *unique* rotor
```math
𝐐(θ, ϕ) = e^{ϕ𝐤/2}\, e^{θ𝐣/2}.
```
Conflating spherical coordinates with geometry of ``𝕊²`` again for a
moment, we can think of this mapping as taking the north pole ``𝐳``
to the point on the sphere with coordinates ``(θ, ϕ)``, and the vector
``𝐱`` onto ``\boldsymbol{θ}`` at that point.

The key point is that, in this case, the mapping doesn't *start* with
``\mathrm{Spin}(3)``, but rather *ends* there.  Thus, we can define a
pullback of a function ``f(𝐐)`` to a function ``f(θ, ϕ)`` in the
usual way — simple composition:
```math
f(θ, ϕ) = f\big(𝐐(θ, ϕ)\big).
```
This is defined for all spin weights, because we've chosen a tangent
basis at each point.  However, if we then try to "wrap" the cylinder
back onto the sphere, we cannot define the function at the poles for
nonzero spin weights, because the tangent basis is not uniquely
defined there.  That is, the function value will not be independent of
the choice of ``ϕ`` at the poles, so the "wrapping" is not well
defined.

This is how we can quite reasonably write spin-weighted spherical
harmonics as functions of spherical coordinates, even though they
cannot be defined as functions on the sphere itself.
