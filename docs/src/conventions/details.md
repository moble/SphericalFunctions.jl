# Details

This page carefully works through all the conventions used in this
package, starting from first principles to motivate the choices and
ensure that each step is on firm footing.  The [previous page](@ref
"Summary") collects the results in a more concise form.

Note that we will use Euler angles and spherical coordinates here, but
*they are not used internally in this package* — though conversion
functions are available.  It is almost always a bad idea to use Euler
angles in *computing*; quaternions are clearly the preferred
representation for numerous reasons.  However, Euler angles are
important for (a) comparing to other sources, and (b) performing
*analytic* integrations.  These are the only two uses we will make of
Euler angles.


## [Three-dimensional space](@id conv_three_dimensional_space)

The space we are working in is naturally three-dimensional Euclidean
space, so we start with a
[right-handed](https://en.wikipedia.org/wiki/Right-hand_rule)
Cartesian coordinate system ``(x, y, z)``.  These also give us the
unit basis vectors ``(𝐱, 𝐲, 𝐳)``.  Note that these basis vectors
are assumed to have unit norm, but we omit the hats just to keep the
notation simple.  Any vector in this space can be written as
```math
𝐯 = v_x 𝐱 + v_y 𝐲 + v_z 𝐳,
```
in which case the Euclidean norm is given by
```math
\| 𝐯 \| = \sqrt{v_x^2 + v_y^2 + v_z^2}.
```
Equivalently, we can write the components of the Euclidean metric as
```math
g_{ij} = \left( \begin{array}{ccc}
  1 & 0 & 0 \\
  0 & 1 & 0 \\
  0 & 0 & 1
\end{array} \right)_{ij}.
```
Note that, because the points of the space are in one-to-one
correspondence with the vectors, we will frequently use a vector to
label a point in space.

### [Spherical coordinates](@id conv_spherical_coordinates)

We will be working on the sphere, so it will be very convenient to use
spherical coordinates ``(r, θ, ϕ)``.  We choose the standard
"physics" conventions for these, in which we relate to the Cartesian
coordinates by
```math
\begin{aligned}
r &= \sqrt{x^2 + y^2 + z^2} &&\in [0, \infty), \\
θ &= \arccos\left(\frac{z}{r}\right) &&\in [0, π], \\
ϕ &= \arctan\left(\frac{y}{x}\right) &&\in [0, 2π),
\end{aligned}
```
where we assume the ``\arctan`` in the expression for ``ϕ`` is
really the two-argument form that gives the correct quadrant.  The
inverse transformation is given by
```math
\begin{aligned}
x &= r \sin θ \cos ϕ, \\
y &= r \sin θ \sin ϕ, \\
z &= r \cos θ.
\end{aligned}
```
We can use this to find the components of the metric in spherical
coordinates:
```math
g_{i'j'}
= \sum_{i,j} \frac{\partial x^i}{\partial x^{i'}} \frac{\partial x^j}{\partial x^{j'}} g_{ij}
= \left( \begin{array}{ccc}
  1 & 0 & 0 \\
  0 & r^2 & 0 \\
  0 & 0 & r^2 \sin^2θ
\end{array} \right)_{i'j'}.
```
The unit coordinate vectors in spherical coordinates are then
```math
\begin{aligned}
𝐧 &= \sin θ \cos ϕ 𝐱 + \sin θ \sin ϕ 𝐲 + \cos θ 𝐳, \\
\boldsymbol{θ} &= \cos θ \cos ϕ 𝐱 + \cos θ \sin ϕ 𝐲 - \sin θ 𝐳, \\
\boldsymbol{ϕ} &= -\sin ϕ 𝐱 + \cos ϕ 𝐲,
\end{aligned}
```
where, again, we omit the hats on the unit vectors to keep the
notation simple.  Conversely, we can express the Cartesian basis
vectors in terms of the spherical basis vectors as
```math
\begin{aligned}
𝐱 &= \sin θ \cos ϕ 𝐧 + \cos θ \cos ϕ \boldsymbol{θ} - \sin ϕ \boldsymbol{ϕ}, 
\\
𝐲 &= \sin θ \sin ϕ 𝐧 + \cos θ \sin ϕ \boldsymbol{θ} + \cos ϕ \boldsymbol{ϕ},
\\
𝐳 &= \cos θ 𝐧 - \sin θ \boldsymbol{θ}.
\end{aligned}
```

One seemingly obvious — but extremely important — fact is that the
unit basis frame ``(𝐱, 𝐲, 𝐳)`` can be rotated onto
``(\boldsymbol{θ}, \boldsymbol{ϕ}, 𝐧)`` by first
rotating through the "polar" angle ``θ`` about the ``𝐲``
axis, and then through the "azimuthal" angle ``ϕ`` about the
``𝐳`` axis.  This becomes important when we consider
spin-weighted functions.

Integration in Cartesian coordinates is, of course, trivial as
```math
\int_{\mathbb{R}^3} f\, d^3𝐫 = \int_{-\infty}^{\infty} \int_{-\infty}^{\infty} \int_{-\infty}^{\infty} f\, dx\, dy\, dz.
```
In spherical coordinates, the integrand involves the square-root of
the determinant of the metric, so we have
```math
\int_{\mathbb{R}^3} f\, d^3𝐫 = \int_0^\infty \int_0^π \int_0^{2π} f\, r^2 \sin θ\, dr\, dθ\, dϕ.
```
Restricting to the unit sphere, we obtain the usual surface element
```math
\int_{𝕊²} f\, d^2\Omega = \int_0^π \int_0^{2π} f\, \sin θ\, dθ\, dϕ.
```
Note that ``\int_{𝕊²} d^2\Omega = 4π``.


## [Four-dimensional space: Quaternions and rotations](@id conv_quaternions)

### Geometric algebra

Given the basis vectors ``(𝐱, 𝐲, 𝐳)`` and the Euclidean norm, we
can define the [geometric
algebra](https://en.wikipedia.org/wiki/Geometric_algebra).  The key
feature is the geometric product, which we could define for any pair
of vectors as ``𝐯`` and ``𝐰`` as
```math
𝐯 𝐰 = 𝐯 ⋅ 𝐰 + 𝐯 ∧ 𝐰,
```
where the dot product is the usual scalar product and the wedge
product is the antisymmetric part of the tensor product — acting just
like the standard [exterior
product](https://en.wikipedia.org/wiki/Exterior_algebra) from the
algebra of [differential
forms](https://en.wikipedia.org/wiki/Differential_form).  The
geometric product is linear, associative, distributive, and has the
property that
```math
𝐯𝐯 = \| 𝐯 \|^2.
```
The most useful properties of the geometric product are that parallel
vectors commute with each other, while orthogonal vectors anticommute.
Since the geometric product is linear, the product of any two vectors
can be decomposed into parallel and orthogonal parts.

The basis for this entire space is then the set
```math
\begin{gather}
𝟏, \\
𝐱, 𝐲, 𝐳,\\
𝐱𝐲, 𝐱𝐳, 𝐲𝐳, \\
𝐱𝐲𝐳.
\end{gather}
```
The standard presentation of quaternions (including the confused
historical development) uses different symbols for these last four
basis elements:
```math
\begin{gather}
𝐢 = 𝐳𝐲 = -𝐲𝐳, \\
𝐣 = 𝐱𝐳 = -𝐳𝐱, \\
𝐤 = 𝐲𝐱 = -𝐱𝐲, \\
𝐈 = 𝐱𝐲𝐳.
\end{gather}
```
Note that each of these squares to -1.  For example, recalling that
orthogonal vectors anticommute, the product is associative, and the
product of a vector with itself is just its squared norm, we have
```math
𝐱𝐲𝐱𝐲 = -𝐱𝐲𝐲𝐱 = -𝐱(𝐲𝐲)𝐱 = -𝐱𝐱 = -1.
```
Any of these could act like the unit imaginary; ``𝐱𝐲`` is probably
the canonical choice.

``𝐈`` is sometimes called the pseudoscalar.  Its inverse is ``𝐈^{-1}
= 𝐳𝐲𝐱 = -𝐱𝐲𝐳``, which can also serve as something very much like
the [Hodge star
operator](https://en.wikipedia.org/wiki/Hodge_star_operator),[^1]
mapping elements to their "dual" elements.  In particular, we have
```math
\begin{aligned}
𝐢 &= 𝐈^{-1}𝐱, \\
𝐣 &= 𝐈^{-1}𝐲, \\
𝐤 &= 𝐈^{-1}𝐳.
\end{aligned}
```
We will see that ``𝐢`` generates right-handed rotations in the
positive sense about ``𝐱``, ``𝐣`` about ``𝐲``, and ``𝐤`` about
``𝐳``.  Moreover, this mapping between ``(𝐱, 𝐲, 𝐳)`` and ``(𝐢,
𝐣, 𝐤)`` is a vector-space isomorphism.  In fact, the reader who is
not familiar with geometric algebra but is familiar with quaternions
may be able to read an expression like ``𝐣 𝐱 𝐣⁻¹`` as if it is just
an abuse of notation, and mentally replace ``𝐱`` with ``𝐢`` to read
those symbols as a valid quaternion expression; both viewpoints are
equally correct by the isomorphism.

[^1]: Note that quaternions will only be spanned by elements made from
      an even number of the basis vectors.  It turns out that those
      with an odd number will produce reflections, rather than
      rotations, when acting on a vector — as discussed below.  This
      explains why quaternions are restricted to just those elements
      with an even number to represent rotations.  For details see any
      geometric algebra text, like [Doran and Lasenby](@cite
      DoranLasenby_2010).

### [Quaternions and Euler angles](@id Quaternions-and-Euler-angles)

Note that there are different conventions for the signs of the ``(𝐢,
𝐣, 𝐤)`` basis.  Everyone agrees that ``𝐢² = 𝐣² = 𝐤² = -1``, but
we could easily flip the sign of any basis element, and these would
still be satisfied.  The identifications we chose above are made to
ensure that ``𝐢`` generates rotations about ``𝐱``, and so on, but
even that depends on how we define quaternions as acting on vectors
(to be discussed below).  A different choice of the latter would
result in all flipping the sign of all three basis elements, which is
a convention that is commonly used — though almost exclusively in
aerospace.  The key expressions that eliminate ambiguity are the
multiplications
```math
\begin{aligned}
𝐢 𝐣 &= 𝐤, \\
𝐣 𝐤 &= 𝐢, \\
𝐤 𝐢 &= 𝐣.
\end{aligned}
```
We can also use these rules above to determine ``𝐢𝐣𝐤 = -𝟏``.  All
four of these equations have flipped signs in other conventions.  See
[Sommer et al.](@cite SommerEtAl_2018) for a discussion of the
different conventions.

We use coordinates ``(W, X, Y, Z)`` on the space of quaternions, so
that a quaternion would be written as
```math
𝐐 = W𝟏 + X𝐢 + Y𝐣 + Z𝐤,
```
though we usually omit the ``𝟏``.  The space of all quaternions is
thus four dimensional.  The norm is just the standard Euclidean norm,
so that the norm of a quaternion is
```math
\| 𝐐 \| = \sqrt{W^2 + X^2 + Y^2 + Z^2}.
```
An important operation is the conjugate, which is defined as
```math
\overline{𝐐} = W - X𝐢 - Y𝐣 - Z𝐤.
```
Note that the squared norm can be written as the quaternion times its
conjugate.  Any nonzero quaternion has an inverse, which is just the
conjugate divided by the squared norm:
```math
𝐐^{-1} = \frac{\overline{𝐐}}{𝐐\overline{𝐐}} = \frac{\overline{𝐐}}{\| 𝐐 \|^2}.
```
The other important operation is exponentiation.  Since a scalar
commutes with any quaternion, including a nonzero scalar component in
the quaternion will simply multiply the result by the exponential of
that scalar component.  Moreover, we will not have any use for such an
exponential, so we assume that the argument to the exponential
function is a "pure" quaternion — that is, one with zero scalar
component.  Moreover, we write it as a unit quaternion ``𝐮`` times
some real number ``\sigma``.  In particular, note that ``𝐮^2 = -1``,
so that it acts like the imaginary unit, which means we already know
how to exponentiate it:
```math
\exp(𝐮\, \sigma) = \cos\sigma + 𝐮\, \sin\sigma.
```
Note that the inverse of the result can be obtained simply by negating
the argument, as usual.

Much as with standard three-dimensional space, we could introduce a
generalization of spherical coordinates, though we use a slight
variant: extended Euler coordinates.  We will see below how to
interpret these as a series of rotations.  For now, we simply state
the relation:
```math
\begin{aligned}
R &= \sqrt{W^2 + X^2 + Y^2 + Z^2} &&\in [0, \infty), \\
α &= \arctan\frac{Z}{W} + \arctan\frac{-X}{Y} &&\in [0, 2π), \\
β &= 2\arccos\sqrt{\frac{W^2+Z^2}{W^2+X^2+Y^2+Z^2}} &&\in [0, π], \\
γ &= \arctan\frac{Z}{W} - \arctan\frac{-X}{Y} &&\in [0, 4π),
\end{aligned}
```
where we again assume the ``\arctan`` in the expressions for ``α`` and
``γ`` is really the two-argument form that gives the correct quadrant,
and if relevant, we use `mod` to limit the values on output.  Note
that here, ``γ`` ranges up to ``4π`` rather than just ``2π``, as
in the standard Euler angles.  This is because we are describing the
space of quaternions, rather than just the space of rotations.  If we
restrict to quaternions with magnitude ``R=1``, we have exactly the
group of unit quaternions ``\mathrm{Spin}(3)=\mathrm{SU}(2)``, which
is a double cover of the rotation group ``\mathrm{SO}(3)``.  This
extended range for ``γ`` is necessary to cover the entire space of
quaternions; if we further restrict to ``[0, 2π)``, we would only
cover the space of rotations.  This and the inclusion of ``R``
identify precisely how this coordinate system extends the standard
Euler angles.

Note that it would also be reasonable to limit ``γ`` to ``2π``, while
allowing ``β`` to range up to ``2π`` to cover the entire space of
quaternions.  This is just somewhat more delicate to compute, and is
simply not conventional.  Also, using ``γ ∈ [0,4π)`` integrates nicely
with our [framework of a telescope](@ref background_domain) with ``γ``
representing the rotation about its line of sight; a full ``4π``
rotation is required for the polarizer to explore the full range of
states of half-integer spin fields.

The inverse transformation is given by
```math
\begin{aligned}
  W &= R\, \cos\frac{β}{2} \cos\frac{α+γ}{2}, \\
  X &= -R\, \sin\frac{β}{2} \sin\frac{α-γ}{2}, \\
  Y &= R\, \sin\frac{β}{2} \cos\frac{α-γ}{2}, \\
  Z &= R\, \cos\frac{β}{2} \sin\frac{α+γ}{2}.
\end{aligned}
```
As with the spherical coordinates, we can use this to find the
components of the metric in our extended Euler coordinates:
```math
g_{i'j'}
= \sum_{i,j} \frac{\partial X^i}{\partial X^{i'}} \frac{\partial X^j}{\partial X^{j'}} g_{ij}
= \left( \begin{array}{cccc}
  1 & 0 & 0 & 0 \\
  0 & \frac{R^2}{4} & 0 & \frac{R^2 \cos β}{4} \\
  0 & 0 & \frac{R^2}{4} & 0 \\
  0 & \frac{R^2 \cos β}{4} & 0 & \frac{R^2}{4}
\end{array} \right)_{i'j'}.
```
The unit basis vectors in extended Euler coordinates in terms of the
unit basis vectors in quaternion coordinates are
```math
\begin{aligned}
𝐑 &= \frac{1}{R} \left(
  \cos \frac{β}{2} \cos \frac{α+γ}{2} 𝟏
  - \sin \frac{β}{2} \sin \frac{α-γ}{2} 𝐢
  + \sin \frac{β}{2} \cos \frac{α-γ}{2} 𝐣
  + \cos \frac{β}{2} \sin \frac{α+γ}{2} 𝐤
\right), \\
\boldsymbol{α} &= \frac{R}{2} \left(
  -\cos \frac{β}{2} \sin \frac{α+γ}{2} 𝟏
  - \sin \frac{β}{2} \cos \frac{α-γ}{2} 𝐢
  - \sin \frac{β}{2} \sin \frac{α-γ}{2} 𝐣
  + \cos \frac{β}{2} \cos \frac{α+γ}{2} 𝐤
\right), \\
\boldsymbol{β} &= \frac{R}{2} \left(
  -\sin \frac{β}{2} \cos \frac{α+γ}{2} 𝟏
  - \cos \frac{β}{2} \sin \frac{α-γ}{2} 𝐢
  + \cos \frac{β}{2} \cos \frac{α-γ}{2} 𝐣
  - \sin \frac{β}{2} \sin \frac{α+γ}{2} 𝐤
\right), \\
\boldsymbol{γ} &= \frac{R}{2} \left(
  -\cos \frac{β}{2} \sin \frac{α+γ}{2} 𝟏
  + \sin \frac{β}{2} \cos \frac{α-γ}{2} 𝐢
  - \sin \frac{β}{2} \cos \frac{α-γ}{2} 𝐣
  - \cos \frac{β}{2} \sin \frac{α+γ}{2} 𝐤
\right).
\end{aligned}
```

### [Invariant measure](@id conv_haar_measure)

Again, integration involves a square-root of the determinant of the
metric, which reduces to ``R^3 \sin β / 8``.  The integral over the
entire space of quaternions is then
```math
\int_{\mathbb{R}^4} f\, d^4𝐐
= \int_{-\infty}^\infty \int_{-\infty}^\infty \int_{-\infty}^\infty \int_{-\infty}^\infty f\, dW\, dX\, dY\, dZ
= \int_0^\infty \int_0^{2π} \int_0^{π} \int_0^{4π} f\, \frac{R^3}{8} \sin β\, dR\, dα\, dβ\, dγ.
```
Restricting to the unit sphere ``R=1``, we obtain the measure on
``\mathrm{Spin}(3)``,
```math
\int_{\mathrm{Spin}(3)} f\, d^3\Omega
= \frac{1}{8} \int_0^{2π} \int_0^{π} \int_0^{4π} f\, \sin β\, dα\, dβ\, dγ,
```
where ``\int_{\mathrm{Spin}(3)} d^3\Omega = 2π^2`` is the volume of
the unit 3-sphere.  This measure is inherited from the Euclidean
metric of ``\mathbb{R}^4``, which is invariant under multiplication by
unit quaternions on either side; it is therefore the (bi-invariant)
Haar measure of the group, and it is the measure we use whenever we
integrate over ``\mathrm{Spin}(3)`` — in particular for the
orthogonality of the Wigner ``𝔇`` matrices [below](@ref
conv_D_orthogonality).  Finally, restricting to the space of rotations
by limiting ``γ`` to ``[0, 2π)``, we can further simplify this to
```math
\int_{\mathrm{SO}(3)} f\, d^3\Omega
= \frac{1}{8} \int_0^{2π} \int_0^{π} \int_0^{2π} f\, \sin β\, dα\, dβ\, dγ,
```
where ``\int_{\mathrm{SO}(3)} d^3\Omega = π^2``.  (Many references
instead normalize the measure on ``\mathrm{SO}(3)`` to ``8π^2`` by
omitting the factor of ``1/8``, or to ``1``; only the relative
normalizations matter below.)  These volume factors are verified
symbolically on the [metrics and integration](@ref
metrics_and_integration) page.

## [Rotations](@id conv_rotations)

We restrict to a unit quaternion ``𝐑``, for which ``W^2 + X^2 + Y^2 +
Z^2 = 1``.  Given this constraint we can, without loss of generality,
write the quaternion as
```math
𝐑
= \exp\left(\frac{\rho}{2} \hat{𝔯}\right)
= \cos\frac{\rho}{2} + \sin\frac{\rho}{2}\, \hat{𝔯},
```
where ``\rho`` is an angle of rotation and ``\hat{𝔯}`` is a
unit "pure-vector" quaternion.  We can multiply a vector ``𝐯`` as
```math
𝐑\, 𝐯\, 𝐑^{-1}.
```
Splitting ``𝐯 = 𝐯_⟂ + 𝐯_∥`` into components perpendicular and
parallel to ``\hat{𝔯}``, we see that ``𝐯_∥`` commutes with
``𝐑`` and ``𝐑^{-1}``, while ``𝐯_⟂`` anticommutes with
``\hat{𝔯}``.  To find the full rotation, we expand the
product:
```math
\begin{aligned}
𝐑\, 𝐯\, 𝐑^{-1}
&= 𝐯_∥
   + \left(\cos\frac{\rho}{2} + \sin\frac{\rho}{2}\, \hat{𝔯}\right)
     𝐯_⟂
     \left(\cos\frac{\rho}{2} - \sin\frac{\rho}{2}\, \hat{𝔯}\right) \\
&= 𝐯_∥
   + \left(\cos\frac{\rho}{2}\, 𝐯_⟂ + \sin\frac{\rho}{2}\, \hat{𝔯}\, 𝐯_⟂\right)
     \left(\cos\frac{\rho}{2} - \sin\frac{\rho}{2}\, \hat{𝔯}\right) \\
&= 𝐯_∥
   + \cos^2\frac{\rho}{2}\, 𝐯_⟂ + \sin\frac{\rho}{2}\, \cos\frac{\rho}{2}\, \hat{𝔯}\, 𝐯_⟂
   - \sin\frac{\rho}{2}\, \cos\frac{\rho}{2}\, 𝐯_⟂ \, \hat{𝔯} - \sin^2\frac{\rho}{2}\, \hat{𝔯}\, 𝐯_⟂\, \hat{𝔯} \\
&= 𝐯_∥
   + \cos^2\frac{\rho}{2}\, 𝐯_⟂ + \sin\frac{\rho}{2}\, \cos\frac{\rho}{2}\, [\hat{𝔯}, 𝐯_⟂] - \sin^2\frac{\rho}{2}\, 𝐯_⟂ \\
&= 𝐯_∥
   + \cos\rho\, 𝐯_⟂ + \sin\rho\, \hat{𝔯}\times 𝐯_⟂
\end{aligned}
```
The final expression shows that this is precisely what we expect when
rotating ``𝐯`` through an angle ``\rho`` (in a positive, right-handed
sense) about the axis ``\hat{𝔯}``.

Note that the presence of two factors of ``𝐑`` in the expression for
rotating a vector explains two things.  First, it explains why the
angle of rotation is twice the angle of the quaternion: one factor of
``𝐑`` either commutes and cancels or anti-commutes and combines with
the the other factor.  Second, it explains why the quaternion group is
a double cover of the rotation group: negating ``𝐑`` results in the
same rotation.  Thus, for any rotation, there are two (precisely
opposite) quaternions that represent it.

### [Euler angles and spherical coordinates](@id conv_euler_angles)

Now that we understand how rotations work, we can provide geometric
intuition for the expressions given above for Euler angles.  The Euler
angles *in our convention* represent an initial rotation through
``γ`` about the ``𝐳`` axis, followed by a rotation through
``β`` about the ``𝐲`` axis, and finally a rotation through
``α`` about the ``𝐳`` axis.  Note that the axes are fixed, and
not subject to any preceding rotations.  More precisely, we can write
the unit quaternion as
```math
𝐑 = \exp\left(\frac{α}{2} 𝐤\right)
    \exp\left(\frac{β}{2} 𝐣\right)
    \exp\left(\frac{γ}{2} 𝐤\right).
```
One of the more important interpretations of a rotor is considering
what it does to the basis triad ``(𝐱, 𝐲, 𝐳)``.  In particular, the
vector ``𝐳`` is rotated onto the point given by spherical coordinates
``(θ, ϕ) = (β, α)``, while ``𝐱`` and ``𝐲`` are
rotated into the plane spanned by the unit basis vectors
``\boldsymbol{θ}`` and ``\boldsymbol{ϕ}`` corresponding to
that point.  If ``γ = 0`` the rotation is precise, meaning that
``𝐱`` is rotated onto ``\boldsymbol{θ}`` and ``𝐲`` onto
``\boldsymbol{ϕ}``; if ``γ ≠ 0`` then they are rotated within
that plane by the angle ``γ`` about the ``𝐧`` axis.
Thus, we identify the spherical coordinates ``(θ, ϕ)`` with
the Euler angles ``(α, β, γ) = (ϕ, θ, 0)``.


## [Rotation and angular-momentum operators](@id conv_operators)

### Complex-valued functions

Starting with Cartesian coordinates and the Euclidean norm on
``\mathbb{R}^3``, we have *constructed* the geometric algebra over
that space, as well as the spaces ``\mathrm{Spin}(3) =
\mathrm{SU}(2)`` (topologically ``𝕊³``), ``\mathrm{SO}(3)``
(topologically ``\mathbb{RP}^3``), and ``𝕊²``.  We will be defining
complex-valued functions on these spaces, and defining operators to
construct and classify them.  In particular, because we have
constructed the spaces, they are naturally supplied with coordinates
that are effectively inherited from the original Cartesian system.  We
will be using these coordinate systems to construct both the operators
and functions.  However, it is important to note that the coordinate
systems may have singularities, which means that the spaces of
coordinates may have different topologies than the spaces they
represent.  For example, Euler angles have topology ``𝕊¹ \times I
\times 𝕊¹`` instead of the ``𝕊³`` and ``\mathbb{RP}^3`` topologies
of the spaces they represent; spherical coordinates have topology
``𝕊¹ \times I`` instead of ``𝕊²``.

Defining functions on the coordinate system of a space is subtly
different from defining functions on the space itself.  For example,
spin-weighted functions are generally written as functions of
(``𝕊²``) spherical coordinates.  However, they *cannot* be defined as
functions on ``𝕊²`` itself; some notion of a reference tangent
direction is needed at each point.  The difference is that spherical
*coordinates* supply a natural choice for the reference tangent
direction: the unit vector in the ``\boldsymbol{θ}`` direction.
This supplies just enough information to define the spin-weighted
functions — though this ends up not being a useful form when more
general transformations or deeper understanding are needed.

Because of this variety of spaces, we will need to use function
composition in several ways; functions defined on one space can be
"lifted" or "lowered" to another via maps between the spaces.  In the
diagram below, the function ``F`` can be used to define the function
``f`` via the mapping ``m`` as ``f = m \circ F``.
```@raw html
<div class="composition-diagram">
<?xml version='1.0' encoding='UTF-8'?>
<!-- This file was generated by dvisvgm 3.2.2 -->
<svg version='1.1' xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink' width='135.540542pt' height ='77.088398pt' viewBox='-77.707098 -79.853543 145.540542 99.088398'>
<g id='page1'>
<text class='f1' x='8.094022' y='-98.111224' transform='matrix(1 0 0 1 -75.80112 37.3849)'>A</text>
<text class='f1' x='135.870232' y='-98.111224' transform='matrix(1 0 0 1 -76.0934 37.3849)'>B</text>
<text class='f0' x='71.835977' y='-30.150075' transform='matrix(1 0 0 1 -75.66273 37.38493)'>C</text>
<path d='M-55.742164-63.2187H54.886736' stroke-width='.39848' stroke-miterlimit='10'/>
<path d='M53.0156-65.60941C53.3906-64.17582 54.234347-63.496133 55.085909-63.218789C54.234347-62.937539 53.3906-62.261758 53.0156-60.82816' stroke-width='.39848' stroke-miterlimit='10' stroke-linecap='round' stroke-linejoin='round'/>
<text class='f1' x='.064136' y='-30.150075' transform='matrix(1 0 0 1 -4.666016 -35.414)'>m</text>
<path d='M-58.027364-56.8789L-8.207034-3.7695' stroke-width='.39848' stroke-miterlimit='10' stroke-dasharray='2.78941 1.59395'/>
<path d='M-7.742183-6.769534C-8.531245-5.511715-8.449215-4.433593-8.070309-3.621092C-8.855465-4.05078-9.925779-4.203125-11.230462-3.499999' stroke-width='.39848' stroke-miterlimit='10' stroke-linecap='round' stroke-linejoin='round'/>
<text class='f1' x='.064136' y='-30.150075' transform='matrix(1 0 0 1 -41.341736 9.2371)'>f</text>
<path d='M57.839836-56.8789L7.785156-3.7656' stroke-width='.39848' stroke-miterlimit='10'/>
<path d='M10.808601-3.492188C9.507814-4.199216 8.433598-4.050779 7.648441-3.621093C8.027348-4.433593 8.113288-5.511721 7.33204-6.769528' stroke-width='.39848' stroke-miterlimit='10' stroke-linecap='round' stroke-linejoin='round'/>
<text class='f1' x='.064136' y='-30.150075' transform='matrix(1 0 0 1 34.957364 9.1264)'>F</text>
</g>
</svg>
</div>
```
For example, ``A`` could be the space of spherical coordinates, ``B``
could be ``\mathrm{Spin}(3)``, and ``F`` could be a spin-weighted
function.  There are many maps from spherical coordinates into
``\mathrm{Spin}(3)``; we expect that all such maps will be related by
rotations from ``\mathrm{SO}(3)``, and in some sense equivalent via
some universality relation.  However, for singular maps — such as
coordinate singularities where multiple coordinate values correspond
to a single "physical" point — we find exceptions to the universality.
These compositions will be useful, in that we can define functions on
the "largest" available space, and extend them to any space that maps
into the first.

In principle, our functions should be defined on ``\mathrm{Spin}(3)``
or even the quaternions in general, though in practice we will define
them on the space of coordinates on those spaces.  In any case, we
will classify the functions by their behavior with respect to actions
of ``\mathrm{Spin}(3)`` on the argument to the function.  Therefore,
we need to consider the general behavior of functions under such
actions.

### [Finite rotations](@id conv_finite_rotations)

We work with functions ``f: A \to \mathbb{C}``, where ``A`` is either
the group of unit quaternions, or the full algebra of quaternions.
Any non-zero quaternion can be expressed as ``e^𝔤`` for
some finite quaternion ``𝔤``, which is referred to as the
"generator" of the action of ``e^𝔤``.  This can act on a
function ``f`` by multiplying the argument by ``e^𝔤``.
However, there is an ambiguity: we could multiply either on the left
or the right:[^2]
```math
f\left(𝐐\right) \mapsto f\left(e^𝔤 𝐐\right)
\qquad \text{or} \qquad
f\left(𝐐\right) \mapsto f\left(𝐐 e^𝔤\right).
```
There is an additional ambiguity, in that this action rotates the
*argument* of the function, whereas we will often prefer to think in
terms of rotating the *function* itself.  For example, our function
may describe the measurement of some field in a particular coordinate
system.  Here, the argument ``𝐐`` describes a particular
value of the coordinates, and ``e^𝔤`` changes the point
under consideration.  If, on the other hand ``e^𝔤``
describes how the field itself is rotated, then we can write the
rotated field as a function ``f'`` which is related to the original
function ``f`` by
```math
f'\left(𝐐\right) = f\left(e^{-𝔤} 𝐐\right)
\qquad \text{or} \qquad
f'\left(𝐐\right) = f\left(𝐐 e^{-𝔤}\right).
```
Note that the exponent is negated, because the action of
``e^𝔤`` on the argument is the inverse of the action of
``e^{-𝔤}`` on the function.  This is a general property of
the action of a group on a space, and is a consequence of the group
action being a homomorphism.

[^2]: In group theory, this type of transformation is often referred
      to as a "translation", even when — as in this case — we would
      usually describe these as rotations.

To validate the signs here, it may be helpful to work through a simple
example involving the sphere ``𝕊²``.  We define a function on
spherical coordinates as
```math
f(θ, ϕ) = \sin θ \sin ϕ.
```
Recall that we can map the spherical coordinates into the Euler
angles, and the Euler angles into the quaternion
```math
(θ, ϕ) \mapsto (ϕ, θ, 0) \mapsto 𝐐
=
\exp\left(\frac{ϕ}{2} 𝐤\right)
\exp\left(\frac{θ}{2} 𝐣\right).
```
It is straightforward to see that we can write ``f`` as a function of
``𝐐`` as
```math
f(𝐐) = \left\langle 𝐐\, 𝐤\, 𝐐^{-1} \right\rangle_{𝐣},
```
where the angle brackets and subscript indicate that we are taking the
``𝐣`` component.  That is, ``f`` is the ``y`` component of
the vector ``𝐳`` rotated by ``𝐐``. 

Now, we imagine rotating the field by an angle ``α`` in the
positive sense about the ``z`` axis.  Visualizing the situation, we
can see that the rotated field should be represented by
```math
f'(θ, ϕ) = \sin θ \sin(ϕ - α).
```
For example, the rotated field evaluated at the point ``(θ, ϕ)
= (π/2, 0)`` along the positive ``x`` axis should correspond to the
original field evaluated at the point ``(θ, ϕ) = (π/2,
-α)``.  This rotation is generated by ``𝔤 = α
𝐤 / 2``, which allows us to immediately calculate
```math
\begin{aligned}
f(e^𝔤 𝐐) &= \sin θ \sin(ϕ + α) &&&
f(𝐐 e^𝔤) &= \sin θ \sin ϕ \\
f(e^{-𝔤} 𝐐) &= \sin θ \sin(ϕ - α) &&&
f(𝐐 e^{-𝔤}) &= \sin θ \sin ϕ.
\end{aligned}
```
Thus, we see that left-multiplication by ``e^{-𝔤}``
corresponds to rotation of the field while leaving the coordinates
fixed; left-multiplication by ``e^𝔤`` corresponds to
rotation of the coordinates while leaving the field fixed; and
right-multiplication by either doesn't affect this function at all.

Of course, right-multiplication using other choices for
``𝔤`` could certainly have some effect on this function,
and this choice of ``𝔤`` could have an effect on other
functions.  Note that right-multiplication can also be interpreted as
left-multiplication, where the generator itself is rotated by the
argument to the function.  That is,
```math
\begin{aligned}
f(𝐐 e^𝔤)
  &= f(𝐐 e^{𝔤} 𝐐^{-1} 𝐐)
  = f(e^{𝔤'} 𝐐) \\
f(𝐐 e^{-𝔤})
  &= f(𝐐 e^{-𝔤} 𝐐^{-1} 𝐐)
  = f(e^{-𝔤'} 𝐐),
\end{aligned}
```
where ``𝔤' = 𝐐 𝔤 𝐐^{-1}``.  In
this example, ``𝔤'`` generates a rotation by an angle
``α`` about the point in question, which leaves that point fixed,
and since this is a scalar function it has no effect on the value.  Of
course, we will see below that changing by a phase proportional to
``α`` is the defining feature of a *spin-weighted* function.

### [Differential rotations](@id conv_L_R_definitions)

We now define a pair of operators that differentiate a function with
respect to infinitesimal rotations we apply to the functions
themselves:
```math
\begin{aligned}
L_{𝔤} f(𝐐) &= \lambda \left. \frac{\partial} {\partial θ} f \left( e^{-θ 𝔤 / 2} 𝐐 \right) \right|_{θ=0}, \\
R_{𝔤} f(𝐐) &= \rho \left. \frac{\partial} {\partial θ} f \left( 𝐐 e^{-θ 𝔤 / 2} \right) \right|_{θ=0}.
\end{aligned}
```
Here, we have introduced the constants ``\lambda`` and ``\rho``
because we will actually be able to derive their values — up to signs
— based on the requirement that raising and lowering operators exist
for each.  Finally, we will choose the signs based on demands that
these operators correspond as naturally as possible to the standard
canonical angular-momentum operators.

Note that when composing operators, it is critical to keep track of
the order of operations, which may look slightly unnatural:
```math
\begin{aligned}
  L_𝔤 L_𝔥 f(𝐐)
  % &= \left. \lambda \frac{\partial} {\partial γ} f'\left(e^{-γ 𝔤 / 2} 𝐐 \right) \right|_{γ=0}, \\
  &= \left. \lambda^2 \frac{\partial} {\partial γ} \frac{\partial} {\partial \eta} f\left(e^{-\eta 𝔥 / 2} e^{-γ 𝔤 / 2} 𝐐 \right) \right|_{γ=\eta=0}, \\
  R_𝔤 R_𝔥 f(𝐐)
  % &= \rho \left. \frac{\partial} {\partial γ} f' \left( 𝐐 e^{-γ 𝔤 / 2} \right) \right|_{γ=0} \\
  &= \left. \rho^2 \frac{\partial} {\partial γ} \frac{\partial} {\partial \eta} f\left( 𝐐 e^{-γ 𝔤 / 2} e^{-\eta 𝔥 / 2} \right) \right|_{γ=\eta=0}.
\end{aligned}
```
We can prove the first of these, for example, by defining
``f'(𝐐) = L_𝔥 f(𝐐)``, then applying the
definition of ``L_𝔤`` to ``f'(𝐐)``, and finally
substituting the definition of ``f'`` back in.  If we failed to use
the correct order of operations, we would get sign errors when trying
to evaluate the commutators.

These operators have some nice properties.  For any scalar ``s``, we have
```math
\begin{aligned}
L_{s 𝔤} &= s L_{𝔤}, \\
R_{s 𝔤} &= s R_{𝔤}.
\end{aligned}
```
Given any basis ``𝐞_n`` for the quaternions, we can use
the multivariable chain rule to expand the operators in terms of
components:
```math
\begin{aligned}
L_{𝔤} &= \sum_n g_n\, L_{𝐞_n}, \\
R_{𝔤} &= \sum_n g_n\, R_{𝐞_n}.
\end{aligned}
```
This implies that vector addition holds more generally:
```math
\begin{aligned}
L_{𝔤 + 𝔥} &= L_{𝔤} + L_{𝔥} \\
R_{𝔤 + 𝔥} &= R_{𝔤} + R_{𝔥}.
\end{aligned}
```
Moreover, we can show that these operators form a Lie algebra with the
commutator as the Lie bracket.  That is, we have
```math
\begin{aligned}
[L_{𝔤}, L_{𝔥}]
    &= \frac{\lambda}{2} L_{[𝔤, 𝔥]},
\\
[R_{𝔤}, R_{𝔥}]
    &= -\frac{\rho}{2} R_{[𝔤, 𝔥]},
\\
[L_{𝔤}, R_{𝔥}] &= 0.
\end{aligned}
```

Conventionally, we single out the ``𝐳`` axis — or
equivalently the generator ``𝐤 = 𝐲𝐱`` — as
a sort of fiducial axis, and ``L_z = L_𝐤`` and ``R_z =
R_𝐤`` as the fiducial operators.  Then, *by definition*,
their raising operators ``L_+`` and ``R_+`` and lowering operators
``L_-`` and ``R_-`` satisfy 
```math
\begin{aligned}
[L_z, L_\pm] &= \pm L_\pm, \\
[R_z, R_\pm] &= \pm R_\pm.
\end{aligned}
```
Assuming that the raising and lowering operators can be written as
linear combinations of the basis operators, these equations imply that
they have no component proportional ``L_𝐳``, and that both of
the remaining components must be nonzero.  This actually allows us to
deduce that ``\lambda^2 = \rho^2 = -1``.  This, in turn, allows us to
deduce the values of the raising and lowering operators up to an
overall factor.  Conventionally the factor is chosen so that
```math
\begin{aligned}
L_\pm &= L_𝐱 \pm i L_𝐲, \\
R_\pm &= R_𝐱 \pm i R_𝐲.
\end{aligned}
```

To pin down ``\lambda`` and ``\rho``, suppose the raising operator for
``L_z`` is some linear combination ``L_+ = a L_𝐱 + b L_𝐲 + c L_𝐳``.
Using the commutator relation above together with ``[𝐤, 𝐢] = 2𝐣``
and ``[𝐤, 𝐣] = -2𝐢``, we find
```math
[L_z, L_+] = \lambda \left( a L_𝐲 - b L_𝐱 \right).
```
Setting this equal to ``L_+`` requires ``c = 0``, ``a = -\lambda b``,
and ``b = \lambda a``, which together imply
```math
\lambda^2 = -1
\qquad \text{and} \qquad
L_+ \propto L_𝐱 + \lambda L_𝐲.
```
Thus the mere *existence* of ladder operators forces ``\lambda = \pm
i``, and the sign of ``\lambda`` is tied to which of the two
combinations we choose to call the raising operator.  Adopting the
universal convention
```math
L_\pm = L_𝐱 \pm i L_𝐲
```
therefore fixes ``\lambda = i``.  The same argument for the right
operators, with the extra minus sign in their commutator relation,
gives ``[R_z, R_+] = -\rho\, (a R_𝐲 - b R_𝐱)``, hence ``\rho^2 =
-1`` and ``R_+ \propto R_𝐱 - \rho R_𝐲``.  Adopting the analogous
convention
```math
R_\pm = R_𝐱 \pm i R_𝐲
```
fixes ``\rho = -i``.  We therefore have
```math
\begin{aligned}
\lambda &= i, \\
\rho &= -i,
\end{aligned}
```
which are the values used in the definitions given in the
[Summary](@ref summary_L_R_definitions).  Two consequences deserve
emphasis, both of which will be borne out below.  First, ``L`` then
coincides with the standard angular-momentum operator of quantum
mechanics (``L_z = -i\partial_\phi`` on the sphere), with the standard
commutation relations ``[L_a, L_b] = i\epsilon_{abc} L_c``.  Second,
``R`` obeys exactly the *same* commutation relations, and its fiducial
component is ``R_z = i\partial_γ``, which we will find to have
eigenvalue ``+s`` on a function of spin weight ``s``.  The opposite
sign for ``\rho`` — as used, for example, in the body-fixed operators
of
[Wikipedia](https://en.wikipedia.org/wiki/Wigner_D-matrix#Properties_of_the_Wigner_D-matrix)
and in earlier versions of this package — leads to "anomalous"
commutation relations with an extra minus sign and to an eigenvalue of
``-s``.

Once the constants are fixed, the commutator relations read
```math
[L_𝔤, L_𝔥] = \frac{i}{2} L_{[𝔤, 𝔥]},
\qquad
[R_𝔤, R_𝔥] = \frac{i}{2} R_{[𝔤, 𝔥]},
\qquad
[L_𝔤, R_𝔥] = 0,
```
and the Casimir operators ``L^2 = L_𝐱^2 + L_𝐲^2 + L_𝐳^2`` and
``R^2`` commute with everything.  In fact ``L^2 = R^2``, because both
equal the Laplacian on ``\mathrm{Spin}(3)`` up to a constant factor —
see the [section on Laplacians](@ref conv_laplacians).

### [Angular-momentum operators in Euler angles](@id conv_L_R_euler)

Having defined ``L`` and ``R`` in terms of the group structure, we can
express them in any coordinate system on ``\mathrm{Spin}(3)``.  The
most useful is the system of Euler angles, because it makes contact
with the standard expressions in the literature.  The procedure is
mechanical: write ``e^{-ϵ𝐮/2}\, 𝐑_{α, β, γ}`` (or ``𝐑_{α, β, γ}\,
e^{-ϵ𝐮/2}``) in terms of its components, extract the new Euler angles
``(α', β', γ')`` as functions of ``ϵ``, differentiate at ``ϵ=0``, and
apply the chain rule.  This is carried out symbolically on the
[``L_j`` and ``R_j`` with Euler angles](@ref euler_angular_momentum)
page, which also verifies the commutation relations claimed above.
The results are collected in the [Summary](@ref summary_L_R_euler);
the essential ones are
```math
L_𝐤 = -i \frac{\partial}{\partial α}
\qquad \text{and} \qquad
R_𝐤 = i \frac{\partial}{\partial γ}.
```
Restricting to spin-weight-0 functions via ``(α, β, γ) = (ϕ, θ, 0)``
reproduces the textbook expressions for ``L_x``, ``L_y``, ``L_z``, and
``L_\pm`` in spherical coordinates, as shown [there](@ref "``L``
operators in spherical coordinates").  The ``R`` operators have no
analog on the 2-sphere: their expressions retain derivatives with
respect to ``γ``, which is the first sign that spin-weighted functions
cannot be defined on ``𝕊²``.

## [Wigner's 𝔇 matrices](@id conv_wigner_D)

### [Rotation operators](@id conv_rotation_operator)

[Sakurai_1994](@citet) says that

> Because rotations affect physical systems, the state ket
> corresponding to a rotated system is expected to look different from
> the state ket corresponding to the original unrotated system.  Given
> a rotation operation ``R``, characterized by a ``3×3`` orthogonal
> matrix ``R``, we associate an operator ``𝒟(R)`` in the appropriate
> ket space such that
> ```math
> |α\rangle_R = 𝒟(R) |α\rangle,
> ```
> where ``|α\rangle_R`` and ``|α\rangle`` stand for the kets of
> the rotated and original system, respectively.

In our setting the "kets" are complex-valued functions on
``\mathrm{Spin}(3)``, and we saw [above](@ref conv_finite_rotations)
that rotating the *field* by ``𝐑`` while leaving the coordinates
fixed means rotating the *argument* by ``𝐑^{-1}``.  We therefore
define the rotation operator ``U(𝐑)`` by
```math
\left[U(𝐑) f\right](𝐐) = f\left(𝐑^{-1}\, 𝐐\right).
```
A short calculation shows that ``U(𝐑_1)\, U(𝐑_2) = U(𝐑_1 𝐑_2)``, so
``U`` is a representation of ``\mathrm{Spin}(3)`` on the space of
functions.  It is unitary with respect to the inner product defined by
the [invariant measure](@ref conv_haar_measure), and it commutes with
every right operator ``R_𝐮`` (left and right multiplication commute).
For an infinitesimal rotation we have
```math
\begin{aligned}
\left[U\left(e^{ϵ 𝐮/2}\right) f\right](𝐑)
&=
f\left(e^{-ϵ 𝐮/2}𝐑\right) \\
&\approx
f\left(𝐑\right) + ϵ \left. \frac{d}{dϵ} \right|_{ϵ=0}
f\left(e^{-ϵ 𝐮/2}𝐑\right) \\
&=
f\left(𝐑\right) - i ϵ L_𝐮 f\left(𝐑\right),
\end{aligned}
```
which is precisely Sakurai's Eq. (3.1.15),
```math
𝒟\left(\hat{𝐧}, dϕ \right)
=
1 - i \left( 𝐉 \cdot \hat{𝐧} \right) dϕ,
```
with ``𝐉 \to 𝐋``.  Exponentiating, ``U(e^{ϑ 𝐮/2}) = \exp(-i ϑ
L_𝐮)``, and in particular for Euler angles
```math
U(𝐑_{α, β, γ})
= \exp[-iα L_z]\, \exp[-iβ L_y]\, \exp[-iγ L_z].
```
This is the operator whose matrix elements every source in our
priority list (LALSuite, Wikipedia, Sakurai, Shankar, Zettili,
Varshalovich et al.) uses to define Wigner's ``𝔇`` matrices; we do
the same.

### [Definition](@id conv_wigner_D_definition)

Let ``|ℓ, m\rangle`` denote an orthonormal set of simultaneous
eigenfunctions of ``L^2``, ``L_z``, and ``R_z``, with eigenvalues
``ℓ(ℓ+1)``, ``m``, and some fixed ``s``.  (The value of ``s`` is
irrelevant here because ``U(𝐑)`` commutes with ``R_z``; for integer
``ℓ`` we may take ``s=0``, so that the ``|ℓ, m\rangle`` are just the
ordinary spherical harmonics defined [below](@ref
conv_spherical_harmonics).)  Because ``U(𝐑)`` commutes with ``L^2``,
it maps the ``(2ℓ+1)``-dimensional eigenspace of ``L^2`` to itself,
and we define Wigner's ``𝔇`` matrix as the matrix of ``U(𝐑)`` in this
basis — Sakurai's Eq. (3.5.42):
```math
𝔇^{(ℓ)}_{m',m}(𝐑)
=
\langle ℓ, m' | U(𝐑) | ℓ, m \rangle,
\qquad \text{equivalently} \qquad
U(𝐑)\, |ℓ, m\rangle = \sum_{m'} |ℓ, m'\rangle\, 𝔇^{(ℓ)}_{m',m}(𝐑).
```
Here ``ℓ \in \{0, \tfrac{1}{2}, 1, \tfrac{3}{2}, \ldots\}`` and ``m',
m \in \{-ℓ, -ℓ+1, \ldots, ℓ\}``, so that ``ℓ``, ``m'``, and ``m`` are
either all integers or all half-integers.  Three points about this
definition deserve emphasis, because they are exactly the points on
which the literature disagrees:

1. **Argument.** ``𝔇`` is a function of the rotation ``𝐑`` itself —
   the rotation applied to the *field* — not of its inverse.
2. **Index order.** The first index ``m'`` is attached to the bra,
   and will be found below to pair with the *first* Euler angle ``α``
   and with the *left* operators ``L``; the second index ``m`` is
   attached to the ket, and pairs with the *last* Euler angle ``γ``
   and the *right* operators ``R``.
3. **Conjugation.** With ``U`` as defined above, ``𝔇`` carries phases
   ``e^{-i m' α}`` and ``e^{-i m γ}`` — see below.  Several sources
   (including [Wigner](@cite Wigner_1959), [Edmonds](@cite
   Edmonds_2016), [Goldberg et al.](@cite GoldbergEtAl_1967), and
   [Boyle (2016)](@cite Boyle_2016), whose convention was used by
   versions of this package before 3.0) define a matrix that is the
   complex conjugate of this one, possibly with additional signs
   ``(-1)^{m'-m}`` or transposition.  The comparisons pages record the
   precise relationship for each source.

### [Basic properties](@id conv_D_symmetries)

Because ``U`` is a unitary representation, ``𝔇`` is too — Sakurai's
Eq. (3.5.46):
```math
𝔇^{(ℓ)}(𝐑_1\, 𝐑_2) = 𝔇^{(ℓ)}(𝐑_1)\, 𝔇^{(ℓ)}(𝐑_2),
\qquad
𝔇^{(ℓ)}(𝟏) = 𝟙,
\qquad
𝔇^{(ℓ)}(𝐑^{-1}) = 𝔇^{(ℓ)}(𝐑)^{-1} = 𝔇^{(ℓ)}(𝐑)^\dagger,
```
where the products are matrix products in the indices ``(m', m)``.
Since ``e^{π 𝐮} = -𝟏`` for any unit vector ``𝐮``, and ``U(e^{π𝐮/2}) =
\exp(-iπ L_𝐮)`` has eigenvalues ``e^{-iπ m}`` on the eigenspace of
``L^2``, we also have
```math
𝔇^{(ℓ)}(-𝐑) = (-1)^{2ℓ}\, 𝔇^{(ℓ)}(𝐑).
```
That is, for integer ``ℓ`` the two quaternions representing a given
rotation give the same matrix, so ``𝔇^{(ℓ)}`` is a representation of
``\mathrm{SO}(3)``; for half-integer ``ℓ`` it is a genuine
representation of ``\mathrm{Spin}(3)`` only, and changes sign under a
rotation through ``2π``.  This is the group-theoretic reason that the
Euler angle ``γ`` must range over ``[0, 4π)`` to cover
``\mathrm{Spin}(3)``.

The other basic symmetry follows from the reality of the ``d`` matrix
introduced below and the structure of the phases:
```math
\overline{𝔇^{(ℓ)}_{m',m}(𝐑)}
=
(-1)^{m'-m}\, 𝔇^{(ℓ)}_{-m',-m}(𝐑).
```
Note that ``m'-m`` is always an integer, so this holds for
half-integer ``ℓ`` as well.

### [Euler-angle form and the ``d`` matrix](@id conv_wigner_d_formula)

Using the Euler-angle factorization of ``U`` and the eigenvalue
property of ``L_z`` on the basis kets, we find — Sakurai's Eq.
(3.5.50) —
```math
\begin{aligned}
𝔇^{(ℓ)}_{m',m}(α, β, γ)
&=
\langle ℓ, m' |
    \exp[-iL_z α]\exp[-iL_y β]\exp[-iL_z γ]
| ℓ, m \rangle \\
&=
e^{-i m' α}\, d^{(ℓ)}_{m',m}(β)\, e^{-i m γ},
\end{aligned}
```
where the "small" ``d`` matrix is
```math
d^{(ℓ)}_{m',m}(β)
=
\langle ℓ, m' | \exp[-iL_y β] | ℓ, m \rangle.
```
We derive the phases very explicitly.  Recall that in general
```math
\left(\left\langle ψ | A\, B\, C | \chi \right\rangle\right)^\ast
=
\left\langle \chi | C^\dag\, B^\dag\, A^\dag | ψ \right\rangle,
\qquad
\left( e^{-i ϵ L_u} \right)^\dag
=
e^{i ϵ L_u^\dag}
=
e^{i ϵ L_u},
```
so that ``\langle ℓ, m' | e^{-iαL_z} = \left(e^{iαL_z} |ℓ, m'\rangle
\right)^\dagger = e^{-i m' α} \langle ℓ, m'|``, and similarly on the
right.  Because ``L_y = (L_+ - L_-)/2i`` has purely imaginary matrix
elements in the standard basis (whose ladder coefficients are real and
positive — the Condon–Shortley convention discussed below), ``d`` is
real.  Wigner's explicit formula for it is
```math
d^{(ℓ)}_{m',m}(β)
=
\sum_{k}
(-1)^{k - m + m'}
\frac{\sqrt{(ℓ+m)!\,(ℓ-m)!\,(ℓ+m')!\,(ℓ-m')!}}
     {(ℓ+m-k)!\,k!\,(ℓ-m'-k)!\,(k-m+m')!}
\left(\cos\frac{β}{2}\right)^{2ℓ+m-m'-2k}
\left(\sin\frac{β}{2}\right)^{2k-m+m'},
```
where the sum runs over all integers ``k`` for which the factorials
have non-negative arguments, namely ``\max(0, m-m') \leq k \leq
\min(ℓ+m, ℓ-m')``.  This is the form given by [Sakurai (Eq.
3.8.33)](@cite Sakurai_1994), [Zettili (Eq. 7.56)](@cite
Zettili_2009), and
[Wikipedia](https://en.wikipedia.org/wiki/Wigner_D-matrix#Wigner_(small)_d-matrix),
and it is the form implemented by [LALSuite](@cite LALSuite_2018).
The simplest cases are
```math
d^{(1/2)}(β) =
\begin{pmatrix}
  \cos\frac{β}{2} & -\sin\frac{β}{2} \\
  \sin\frac{β}{2} & \phantom{-}\cos\frac{β}{2}
\end{pmatrix},
\qquad
d^{(1)}(β) =
\begin{pmatrix}
  \frac{1+\cos β}{2} & -\frac{\sin β}{\sqrt{2}} & \frac{1-\cos β}{2} \\
  \frac{\sin β}{\sqrt{2}} & \cos β & -\frac{\sin β}{\sqrt{2}} \\
  \frac{1-\cos β}{2} & \frac{\sin β}{\sqrt{2}} & \frac{1+\cos β}{2}
\end{pmatrix},
```
with rows indexed by ``m'`` and columns by ``m``, both *decreasing*
from ``ℓ`` to ``-ℓ`` (the usual textbook layout).  The ``d`` matrix
satisfies
```math
\begin{aligned}
d^{(ℓ)}_{m',m}(β) &= (-1)^{m'-m}\, d^{(ℓ)}_{m,m'}(β) = d^{(ℓ)}_{-m,-m'}(β), \\
d^{(ℓ)}_{m',m}(-β) &= d^{(ℓ)}_{m,m'}(β), \\
d^{(ℓ)}_{m',m}(π - β) &= (-1)^{ℓ+m'}\, d^{(ℓ)}_{m',-m}(β), \\
d^{(ℓ)}_{m',m}(π) &= (-1)^{ℓ-m}\, δ_{m',-m}, \\
d^{(ℓ)}_{m',m}(2π + β) &= (-1)^{2ℓ}\, d^{(ℓ)}_{m',m}(β),
\end{aligned}
```
which are the standard relations found, e.g., on Wikipedia.  Note
that this package does not compute ``d`` from the formula above, which
is numerically disastrous for large ``ℓ``; the formula is the
*definition* against which the recursive algorithms (documented in
the Notes) are tested.

### [Differential relations](@id conv_D_differential)

Viewed as functions on ``\mathrm{Spin}(3)``, the matrix elements of
``𝔇`` are eigenfunctions of our differential operators.  The
representation property gives ``𝔇_{m',m}(𝐑\, e^{-ϵ𝐮/2}) = \sum_k
𝔇_{m',k}(𝐑)\, 𝔇_{k,m}(e^{-ϵ𝐮/2})``, and ``𝔇_{k,m}(e^{-ϵ𝐮/2}) =
\langle ℓ, k| e^{iϵL_𝐮} |ℓ, m\rangle``, so differentiating the
definitions of ``R_𝐮`` and (similarly) ``L_𝐮`` yields
```math
R_𝐮\, 𝔇^{(ℓ)}_{m',m} = \sum_k 𝔇^{(ℓ)}_{m',k}\, (J_𝐮)_{k,m},
\qquad
L_𝐮\, 𝔇^{(ℓ)}_{m',m} = -\sum_k (J_𝐮)_{m',k}\, 𝔇^{(ℓ)}_{k,m},
```
where ``(J_𝐮)_{k,m} = \langle ℓ, k | L_𝐮 | ℓ, m \rangle`` are the
standard ``(2ℓ+1)``-dimensional angular-momentum matrices.  Inserting
the standard matrix elements of ``J_z`` and ``J_\pm`` gives
```math
\begin{aligned}
L^2\, 𝔇^{(ℓ)}_{m',m} &= R^2\, 𝔇^{(ℓ)}_{m',m} = ℓ(ℓ+1)\, 𝔇^{(ℓ)}_{m',m}, \\
L_z\, 𝔇^{(ℓ)}_{m',m} &= -m'\, 𝔇^{(ℓ)}_{m',m},
&
R_z\, 𝔇^{(ℓ)}_{m',m} &= m\, 𝔇^{(ℓ)}_{m',m}, \\
L_\pm\, 𝔇^{(ℓ)}_{m',m} &= -\sqrt{(ℓ \pm m')(ℓ \mp m' + 1)}\, 𝔇^{(ℓ)}_{m' \mp 1, m},
&
R_\pm\, 𝔇^{(ℓ)}_{m',m} &= \sqrt{(ℓ \mp m)(ℓ \pm m + 1)}\, 𝔇^{(ℓ)}_{m', m \pm 1}.
\end{aligned}
```
These can be checked directly on the Euler-angle form: ``L_z =
-i\partial_α`` acting on ``e^{-im'α}`` gives ``-m'``, while ``R_z =
i\partial_γ`` acting on ``e^{-imγ}`` gives ``m``.  The minus sign and
the reversed direction of the ``L`` ladder are unavoidable
consequences of defining ``𝔇`` through matrix elements of ``U(𝐑)``:
as a function of its argument, ``𝔇_{m',m}`` behaves like a *row* of
expansion coefficients rather than like a basis function.  It is the
complex conjugate ``\overline{𝔇_{m',m}}`` that behaves like a
wavefunction, with ``L_z`` eigenvalue ``+m'`` and ``R_z`` eigenvalue
``-m`` — which is exactly why the spin-weighted spherical harmonics
below are proportional to ``\overline{𝔇_{m,-s}}``.  (Wikipedia states
the same fact in the language of the symmetric top: it is
``D^{j\ast}_{m'm}`` that is the eigenfunction of the space-fixed
``𝒥_z`` with eigenvalue ``m'`` and of the body-fixed ``𝒫_z = -R_z``
with eigenvalue ``m``.)

### [Orthogonality and completeness](@id conv_D_orthogonality)

By the Peter–Weyl theorem (Schur orthogonality), the matrix elements
of the irreducible unitary representations form a complete orthogonal
basis of ``L^2(\mathrm{Spin}(3))``.  With the [invariant measure](@ref
conv_haar_measure) normalized to total volume ``2π^2``, the
orthogonality relation reads
```math
\int_{\mathrm{Spin}(3)}
  \overline{𝔇^{(ℓ')}_{m'_1, m_1}(𝐑)}\,
  𝔇^{(ℓ)}_{m', m}(𝐑)\,
  d^3Ω
=
\frac{2π^2}{2ℓ+1}\, δ_{ℓ', ℓ}\, δ_{m'_1, m'}\, δ_{m_1, m}.
```
Restricting the sum over ``ℓ`` to integers gives the corresponding
statement for ``\mathrm{SO}(3)``, whose volume in this normalization
is ``π^2``.  This is the sense in which the ``𝔇`` matrices generalize
Fourier analysis from the circle to the rotation group, and it is the
basis of the spin-weighted spherical-harmonic transforms implemented
in this package.

## [Spherical harmonics and spin-weighted spherical harmonics](@id conv_harmonics)

### [Spherical harmonics](@id conv_spherical_harmonics)

Fortunately, there does not seem to be any disagreement in the physics
literature about the definition of the spherical harmonics ``Y_{ℓ,
m}(θ, ϕ)``; everyone uses the Condon–Shortley convention, and our
``Y_{ℓ,m}`` is the standard one.  The only place where care is needed
is that many sources *define* ``Y_{ℓ,m}`` in terms of associated
Legendre functions, for which there is a sign ambiguity ``(-1)^m``
that is sometimes absorbed into the Legendre function and sometimes
not.  To be unambiguous, we go back to the original.

The [Condon–Shortley](@cite CondonShortley_1935) phase convention is
a choice of phase factors in the definition of the spherical harmonics
that requires the coefficients in
```math
L_{\pm} |ℓ,m\rangle = α^{\pm}_{ℓ,m} |ℓ, m \pm 1\rangle
```
to be real and positive, together with the requirement that
``Y_{ℓ,0}`` be real and positive on the positive ``𝐳`` axis.  The
reasoning behind this choice is explained more fully in Section 2 of
[Ufford and Shortley (1932)](@cite UffordShortley_1932).  As a
practical matter, Condon and Shortley's Eq. (15) of section 4³ (page
52) gives
```math
\Theta(ℓ, m) = (-1)^ℓ \sqrt{\frac{2ℓ+1}{2} \frac{(ℓ+m)!}{(ℓ-m)!}}
\frac{1}{2^ℓ ℓ!} \frac{1}{\sin^mθ}
\frac{d^{ℓ-m}}{d(\cos θ)^{ℓ-m}} \sin^{2ℓ}θ,
```
and the spherical harmonic is ``Y_{ℓ,m}(θ, ϕ) = \Theta(ℓ, m)\,
\Phi(m)`` with ``\Phi(m) = e^{imϕ} / \sqrt{2π}`` from their Eq. (5).
This expression is tested directly against this package (see the
Condon–Shortley comparison page), so we can definitively say that
***the spherical-harmonic functions provided by this package obey the
Condon–Shortley phase convention***.  For example,
```math
Y_{0,0} = \frac{1}{\sqrt{4π}},
\qquad
Y_{1,0} = \sqrt{\frac{3}{4π}} \cos θ,
\qquad
Y_{1,\pm 1} = \mp \sqrt{\frac{3}{8π}} \sin θ\, e^{\pm iϕ}.
```

An equivalent closed form that avoids Legendre functions altogether —
and is the ``s=0`` case of the general expression for spin-weighted
spherical harmonics given below — is
```math
\begin{aligned}
  Y_{ℓ,m}(θ, ϕ)
  &=
  \sqrt{\frac{2ℓ+1}{4π}}\, e^{imϕ}
  \sum_{k = k_1}^{k_2}
  \frac{(-1)^k\, ℓ!\, [(ℓ+m)!(ℓ-m)!]^{1/2}}
  {(ℓ+m-k)!\,(ℓ-k)!\,k!\,(k-m)!}
  \left(\cos\frac{θ}{2}\right)^{2ℓ+m-2k}
  \left(\sin\frac{θ}{2}\right)^{2k-m},
\end{aligned}
```
where ``k_1 = \max(0, m)`` and ``k_2 = \min(ℓ+m, ℓ)``.

Now, ``Y_{ℓ,m}`` is a function on ``𝕊²``, but we have argued that the
natural domain for everything in this package is ``\mathrm{Spin}(3)``.
We lift it in the obvious way: for a unit quaternion ``𝐑`` that
rotates ``𝐳`` onto the direction with spherical coordinates ``(θ,
ϕ)`` — that is, for ``𝐑 = 𝐑_{θ, ϕ}\, e^{γ 𝐤 / 2}`` with any ``γ`` —
we define ``Y_{ℓ,m}(𝐑) = Y_{ℓ,m}(θ, ϕ)``.  The lifted function is
independent of ``γ``, so ``R_z Y_{ℓ,m} = 0``: the spherical harmonics
have spin weight zero.  Comparing the closed form above with Wigner's
formula for ``d``, we find
```math
Y_{ℓ,m}(𝐑)
=
\sqrt{\frac{2ℓ+1}{4π}}\, \overline{𝔇^{(ℓ)}_{m,0}(𝐑)},
```
which is the standard relation (e.g., on
[Wikipedia](https://en.wikipedia.org/wiki/Wigner_D-matrix#Relation_to_spherical_harmonics_and_Legendre_polynomials)),
and — recalling the [differential relations](@ref conv_D_differential)
— is consistent with ``L_z Y_{ℓ,m} = m Y_{ℓ,m}`` and ``R_z Y_{ℓ,m} =
0``.  The complex conjugate is essential; without it the sign of ``m``
would be wrong.

### [Rotation of spherical harmonics](@id conv_rotation_law)

The definition of ``𝔇`` as the matrix of ``U(𝐑)`` in the basis
``|ℓ, m\rangle`` is precisely the statement of how spherical
harmonics transform under rotation:
```math
\left[U(𝐑)\, Y_{ℓ,m}\right](𝐐)
=
Y_{ℓ,m}\left(𝐑^{-1}\, 𝐐\right)
=
\sum_{m'} 𝔇^{(ℓ)}_{m',m}(𝐑)\, Y_{ℓ,m'}(𝐐).
```
In words: the field ``Y_{ℓ,m}`` rotated by ``𝐑`` is a combination of
the unrotated ``Y_{ℓ,m'}`` with coefficients given by the *column*
``m`` of ``𝔇(𝐑)``.  Using unitarity, the same law can be written
with the rotation in the other position,
```math
Y_{ℓ,m}\left(𝐑\, 𝐐\right)
=
\sum_{m'} \overline{𝔇^{(ℓ)}_{m,m'}(𝐑)}\, Y_{ℓ,m'}(𝐐),
```
which is the form given by Wigner (his Eq. A.8) and is the one that
arises when evaluating the harmonics at a rotated *point*.  The first
form — with ``𝐑^{-1}`` in the argument and no conjugate — is the one
used by Sakurai, Le Bellac, Torres del Castillo, and Zettili, and we
take it as canonical.  Exactly the same law holds for the
spin-weighted spherical harmonics defined below, with ``Y`` replaced
by ``{}_sY``, because ``U(𝐑)`` commutes with ``R_z``.

### [Spin-weighted functions](@id conv_spin_weight)

[Newman_1966](@citet) define the spherical tangent basis vectors as
```math
m^\mu = \frac{1}{\sqrt{2}} \left(
    \boldsymbol{θ} + i \boldsymbol{ϕ}
\right)^\mu
```
and discuss spin weight in terms of the rotation
```math
(m^\mu)' = e^{iψ} m^\mu,
```
where the tangent basis rotates but we are "keeping the coordinates
fixed".  They then define a function ``\eta`` to have spin weight
``s`` if it transforms as
```math
\eta' = e^{isψ} \eta.
```
Such functions are generally the result of contracting a tensor field
with some number of ``m^\mu`` and some number of ``\bar{m}^\mu``
vectors (though spinor extensions resulting in half-integer spin
weights are also possible).

This definition shows that it is *impossible* to define
spin-weighted functions on the 2-sphere alone; the 2-sphere includes
no information about the directions of basis vectors in its tangent
space.  Instead, spin-weighted functions live on the "unit tangent
bundle" of the 2-sphere, which is homeomorphic to the 3-sphere — the
space of unit quaternions.  Thus, we think of spin-weighted functions
as functions on ``\mathrm{Spin}(3)``, and frequently discuss them in
terms of Euler angles.

As we saw [above](@ref conv_euler_angles), the basis ``m^\mu``
corresponds to the Euler angles ``(ϕ, θ, 0)``, while ``(m^\mu)'``
corresponds to ``(ϕ, θ, -ψ)``.  Written as a function of Euler
angles, the defining transformation becomes
```math
\eta(ϕ, θ, -ψ) = e^{isψ} \eta(ϕ, θ, 0),
\qquad \text{or} \qquad
\eta(α, β, γ) = e^{-isγ} \eta(α, β, 0),
```
and, independently of the choice of Euler angles,
```math
\eta\left(𝐐\, e^{γ 𝐤/2}\right) = e^{-isγ}\, \eta(𝐐).
```
This is the crucial definition giving us the behavior of spin-weighted
functions: applying ``R_z = i\partial_γ`` we find
```math
R_z\, \eta = s\, \eta,
```
so ***spin-weighted functions are eigenfunctions of ``R_z`` with
eigenvalue equal to the spin weight***.  It was to make this sign
come out positive that we chose ``\rho = -i`` in the definition of
``R``.

The spin-raising and -lowering operators — canonically denoted ``\eth``
and ``\bar{\eth}`` — were introduced by [Newman and Penrose](@cite
Newman_1966) in their Eq. (3.8) as
```math
\begin{aligned}
\eth \eta &= -\sin^s θ \left\{
        \frac{\partial}{\partial θ}
        + \frac{i}{\sin θ} \frac{\partial}{\partial ϕ}
    \right\} \left(\eta \sin^{-s} θ\right), \\
\bar{\eth} \eta &= -\sin^{-s} θ \left\{
        \frac{\partial}{\partial θ}
        - \frac{i}{\sin θ} \frac{\partial}{\partial ϕ}
    \right\} \left(\eta \sin^{s} θ\right).
\end{aligned}
```
These expressions look like operators on the 2-sphere, but that
appearance is misleading: they depend on ``s``, which encodes the
missing information about the tangent basis.  Writing out the
Euler-angle expressions for ``R_𝐱`` and ``R_𝐲``, replacing
``\partial_γ`` by ``-is`` (as is legitimate when acting on a function
of definite spin weight), and converting the remaining angles to
spherical coordinates — see the [calculation page](@ref euler_R_S2) —
we find that these are just the ladder operators of ``R_z``:
```math
\eth = R_+ = R_𝐱 + i R_𝐲,
\qquad
\bar{\eth} = -R_- = -\left(R_𝐱 - i R_𝐲\right).
```
The minus sign in the second relation is real and worth remembering.
It arises because ``\bar{\eth}`` is defined as the complex-conjugate
*operator* of ``\eth`` (``\bar{\eth}\eta = \overline{\eth\bar{\eta}}``),
and ``R_\pm`` contain an explicit factor of ``i``, so that
``\overline{R_+ \bar{\eta}} = -R_- \eta``.  As a result ``R_-`` has
the Condon–Shortley-positive coefficient ``R_- {}_sY_{ℓ,m} =
+\sqrt{(ℓ+s)(ℓ-s+1)}\, {}_{s-1}Y_{ℓ,m}``, whereas Newman–Penrose's
operator has ``\bar{\eth}\, {}_sY_{ℓ,m} = -\sqrt{(ℓ+s)(ℓ-s+1)}\,
{}_{s-1}Y_{ℓ,m}``, exactly as in Eq. (2.7b) of [Goldberg et
al.](@cite GoldbergEtAl_1967), and ``\bar{\eth}\eth\, {}_sY_{ℓ,m} =
-(ℓ-s)(ℓ+s+1)\, {}_sY_{ℓ,m}``.  Since the literature on spin-weighted
functions universally uses Newman–Penrose's ``\bar{\eth}``, we keep
their sign and simply note that ``\bar{\eth} \neq R_-``.  Either way,
``[R_z, \eth] = \eth`` and ``[R_z, \bar{\eth}] = -\bar{\eth}``.

### [Spin-weighted spherical harmonics](@id conv_swsh)

We now have everything needed to define the spin-weighted spherical
harmonics as functions on ``\mathrm{Spin}(3)``.  We require ``{}_sY_{ℓ,m}``
to be a simultaneous eigenfunction of ``L^2``, ``L_z``, and ``R_z``
with eigenvalues ``ℓ(ℓ+1)``, ``m``, and ``s``, normalized in the
standard way, with phases fixed by the Condon–Shortley requirement
that *both* sets of ladder operators have real, positive coefficients:
```math
\begin{aligned}
L_\pm\, {}_sY_{ℓ,m} &= \sqrt{(ℓ \mp m)(ℓ \pm m + 1)}\, {}_sY_{ℓ,m \pm 1}, \\
R_\pm\, {}_sY_{ℓ,m} &= \sqrt{(ℓ \mp s)(ℓ \pm s + 1)}\, {}_{s \pm 1}Y_{ℓ,m}.
\end{aligned}
```
The [differential relations](@ref conv_D_differential) show that
``\overline{𝔇^{(ℓ)}_{m,-s}}`` has exactly the required eigenvalues, so
``{}_sY_{ℓ,m} = c_{ℓ,s}\, \overline{𝔇^{(ℓ)}_{m,-s}}`` for some
constants ``c_{ℓ,s}``, which the ladder conditions relate to one
another.  Because ``R_\pm`` contain an explicit factor of ``i``, they
do not commute with complex conjugation; rather ``R_+ \bar{g} =
-\overline{R_- g}``, and therefore
```math
R_+ \overline{𝔇^{(ℓ)}_{m,-s}}
= -\overline{R_- 𝔇^{(ℓ)}_{m,-s}}
= -\sqrt{(ℓ-s)(ℓ+s+1)}\; \overline{𝔇^{(ℓ)}_{m,-s-1}},
```
while ``L_+ \overline{𝔇^{(ℓ)}_{m,-s}} = +\sqrt{(ℓ-m)(ℓ+m+1)}\,
\overline{𝔇^{(ℓ)}_{m+1,-s}}``.  Positivity of the ``L_+`` coefficient
is thus automatic, but positivity of the ``R_+`` coefficient requires
``c_{ℓ,s+1} = -c_{ℓ,s}``: the phase must alternate with ``s``.
Identifying the ``s=0`` case with the standard spherical harmonics
fixes ``c_{ℓ,0} = \sqrt{(2ℓ+1)/4π}``, and we obtain the definition
```math
\boxed{
{}_sY_{ℓ,m}(𝐑)
=
(-1)^s \sqrt{\frac{2ℓ+1}{4π}}\; \overline{𝔇^{(ℓ)}_{m,-s}(𝐑)}
}
```
The factor ``(-1)^s`` is therefore not an arbitrary phase: it is the
Condon–Shortley convention applied to ``R_\pm``, and it is exactly the
factor that appears in the standard expressions of the
gravitational-wave literature — [LALSuite](@cite LALSuite_2018), the
NINJA paper [AjithEtAl_2011](@cite), and (up to their differing
``𝔇`` convention) [Goldberg et al.](@cite GoldbergEtAl_1967).  In
Euler angles, ``\overline{𝔇_{m,-s}(ϕ, θ, γ)} = e^{imϕ}\, d_{m,-s}(θ)\,
e^{-isγ}``, so restricting to ``γ=0`` gives the familiar
```math
{}_sY_{ℓ,m}(θ, ϕ)
=
(-1)^s \sqrt{\frac{2ℓ+1}{4π}}\; d^{(ℓ)}_{m,-s}(θ)\, e^{imϕ},
```
and inserting Wigner's formula for ``d`` gives the explicit expression
```math
\begin{aligned}
  {}_{s}Y_{ℓ,m}(θ, ϕ)
  &=
  (-1)^s\sqrt{\frac{2ℓ+1}{4π}}\, e^{imϕ}
  \sum_{k = k_1}^{k_2}
  \frac{(-1)^k[(ℓ+m)!(ℓ-m)!(ℓ-s)!(ℓ+s)!]^{1/2}}
  {(ℓ+m-k)!\,(ℓ+s-k)!\,k!\,(k-s-m)!}
  \left(\cos\frac{θ}{2}\right)^{2ℓ+m+s-2k}
  \left(\sin\frac{θ}{2}\right)^{2k-s-m},
\end{aligned}
```
where ``k_1 = \max(0, m+s)`` and ``k_2 = \min(ℓ+m, ℓ+s)``.  This is
Eq. (II.7)–(II.8) of the NINJA paper and the expression coded in
LALSuite; both are tested against this package on the comparisons
pages.  Again, this package does not *use* this form; it is shown here
to make comparison with other sources easy.

Summarizing the properties of ``{}_sY_{ℓ,m}`` that follow directly
from those of ``𝔇``:
```math
\begin{aligned}
L^2\, {}_sY_{ℓ,m} &= ℓ(ℓ+1)\, {}_sY_{ℓ,m},
\qquad
L_z\, {}_sY_{ℓ,m} = m\, {}_sY_{ℓ,m},
\qquad
R_z\, {}_sY_{ℓ,m} = s\, {}_sY_{ℓ,m}, \\
{}_sY_{ℓ,m}\left(𝐑\, e^{γ𝐤/2}\right) &= e^{-isγ}\, {}_sY_{ℓ,m}(𝐑),
\qquad
{}_sY_{ℓ,m}\left(𝐑^{-1}\, 𝐐\right) = \sum_{m'} 𝔇^{(ℓ)}_{m',m}(𝐑)\, {}_sY_{ℓ,m'}(𝐐), \\
\overline{{}_sY_{ℓ,m}} &= (-1)^{m+s}\, {}_{-s}Y_{ℓ,-m},
\qquad
{}_sY_{ℓ,m}(-𝐑) = (-1)^{2ℓ}\, {}_sY_{ℓ,m}(𝐑), \\
\int_{\mathrm{Spin}(3)} \overline{{}_{s'}Y_{ℓ',m'}}\; {}_sY_{ℓ,m}\; d^3Ω
&= \frac{π}{2}\, δ_{ℓ',ℓ}\, δ_{m',m}\, δ_{s',s},
\qquad
\int_{𝕊^2} \overline{{}_{s}Y_{ℓ',m'}(θ,ϕ)}\; {}_sY_{ℓ,m}(θ,ϕ)\; \sin θ\, dθ\, dϕ
= δ_{ℓ',ℓ}\, δ_{m',m}.
\end{aligned}
```
The first orthogonality relation is the natural one on
``\mathrm{Spin}(3)``; the second is the conventional one, which
integrates over the 2-sphere at fixed ``γ`` and is valid only when
both functions have the *same* spin weight (functions of different
spin weight are not orthogonal on ``𝕊²``, and the integral is not even
coordinate-independent).  The conjugation relation is the one given by
Goldberg et al. (their Eq. 2.6); note that the exponent ``m+s`` is
always an integer.

### [Half-integer indices](@id conv_swsh_half_integer)

Everything above is valid when ``ℓ``, ``m``, and ``s`` are all
half-integers, with one exception: the factor ``(-1)^s`` is then ``\pm
i``, and a branch must be chosen.  It is worth being precise about
what is and is not forced here.  The ladder conditions determine all
the ``{}_sY_{ℓ,m}`` for a given ``ℓ`` up to one overall phase ``c_ℓ``,
and there is no ``s=0`` member to anchor it.  Requiring the
conjugation relation ``\overline{{}_sY_{ℓ,m}} = (-1)^{m+s}\,
{}_{-s}Y_{ℓ,-m}`` to hold as written forces ``c_ℓ`` to be real, so
``c_ℓ = \pm 1``; requiring the standard product formula (with real
Clebsch–Gordan coefficients) to hold across integer and half-integer
``ℓ`` forces the same sign for every half-integer ``ℓ``.  Nothing
forces that one remaining global sign, because ``(-1)^{2s}`` is itself
multiplicative.  This is the exact analog of the arbitrary choice
"``Y_{ℓ,0}`` is positive at the north pole" in the integer case.  We
choose the principal branch,
```math
(-1)^s \equiv e^{iπ s} = i^{2s},
```
so that the definition above is well defined for all ``s``, and the
anchor condition can be stated uniformly as
```math
{}_sY_{ℓ,-s}(𝟏) = i^{2s} \sqrt{\frac{2ℓ+1}{4π}}
```
(recall that ``𝔇^{(ℓ)}(𝟏)`` is the identity, so only ``m=-s``
survives at ``𝐑=𝟏``).  The same branch is used for the prefactor of
the explicit formula above.

A related trap: using the conjugation symmetry of ``𝔇``, the
definition can be rewritten without a complex conjugate as
```math
{}_sY_{ℓ,m}(𝐑)
=
(-1)^s (-1)^{m+s} \sqrt{\frac{2ℓ+1}{4π}}\; 𝔇^{(ℓ)}_{-m,s}(𝐑)
=
e^{iπ(m+2s)} \sqrt{\frac{2ℓ+1}{4π}}\; 𝔇^{(ℓ)}_{-m,s}(𝐑).
```
For integer indices this is simply ``(-1)^m \sqrt{(2ℓ+1)/4π}\,
𝔇^{(ℓ)}_{-m,s}(𝐑)``, which is a perfectly good alternative form.  For
half-integer ``s``, however, the naive replacement ``e^{iπ(m+2s)} \to
(-1)^m = e^{iπm}`` is off by ``(-1)^{2s} = -1``.  We therefore treat
the conjugate form as *the* definition and the ``𝔇_{-m,s}`` form as an
integer-index corollary.

## [Laplacians](@id conv_laplacians)

[Bander_1966](@citet) show that Wigner's D matrices (extended to the
full space of quaternions with arbitrary norm) are harmonic with
respect to the Laplacian of the full 4-D space.  We also know that
```math
\Delta_{𝕊^{n-1}} f(x) = \Delta_{\mathbb{R}^n} f(x/|x|),
```
and
```math
\Delta_{\mathbb{R}^n} f(x)
=
\frac{1}{r^{n-1}} \frac{\partial}{\partial r} \left( r^{n-1} \frac{\partial f}{\partial r} \right)
+
\frac{1}{r^2} \Delta_{𝕊^{n-1}} f.
```
These imply that the restriction to the space of unit quaternions is
not harmonic with respect to the Laplacian on the 3-sphere, but is an
eigenfunction with eigenvalue ``-2ℓ(2ℓ+2) = -4ℓ(ℓ+1)``, because the
matrix elements of ``𝔇^{(ℓ)}`` are homogeneous polynomials of degree
``2ℓ`` in ``(W, X, Y, Z)``.  (Recall that ``𝔇^{(1/2)}`` is linear in
the components of the quaternion.)

```math
\frac{1}{r^{n-1}} \frac{\partial}{\partial r} \left( r^{n-1} \frac{\partial f}{\partial r} \right)
=
\frac{1}{r^{n-1}} \left( r^{n-1} \frac{\partial}{\partial r} \frac{\partial f}{\partial r} \right)
+
\frac{1}{r^{n-1}} \frac{\partial}{\partial r} \left( r^{n-1} \right) \frac{\partial f}{\partial r}
=
\frac{\partial^2 f}{\partial r^2}
+
\frac{n-1}{r^{n-1}} r^{n-2} \frac{\partial f}{\partial r}
=
\frac{\partial^2 f}{\partial r^2}
+
\frac{n-1}{r} \frac{\partial f}{\partial r}
```

```math
\frac{\partial^2 f}{\partial r^2}
+
\frac{n-1}{r} \frac{\partial f}{\partial r}
=
\frac{f}{r^ℓ} \frac{\partial^2 r^ℓ}{\partial r^2}
+
\frac{f}{r^ℓ} \frac{n-1}{r} \frac{\partial r^ℓ}{\partial r}
=
ℓ(ℓ-1) \frac{f}{r^ℓ} r^{ℓ-2}
+
ℓ \frac{f}{r^ℓ} \frac{n-1}{r} r^{ℓ-1}
=
ℓ(ℓ-1) \frac{f}{r^2}
+
ℓ (n-1) \frac{f}{r^2}
=
ℓ(ℓ+n-2) \frac{f}{r^2}
\to
ℓ(ℓ+2) \frac{f}{r^2}
```

Note that [Lee_2012](@citet) points out that there is a sign ambiguity
in the Laplacian.  As I see it, the geometry community skews toward
including a negative sign (which means that all eigenvalues are
non-negative), while the physics community skews toward excluding it
(which means that all eigenvalues are non-positive).  It's also easy
to prove that on a closed and connected manifold, eigenfunctions with
distinct eigenvalues are orthogonal, since
```math
(\lambda_u - \lambda_v) \int f_u f_v
= \int (\lambda_u f_u) f_v - \int f_u (\lambda_v f_v)
= \int (\Delta f_u) f_v - \int f_u (\Delta f_v) = 0
```
(the last equality by Green's theorem).  Since the eigenvalues are
distinct, this can only be true if ``\int f_u f_v=0``.

[BoydPetschek_2014](@citet) produced an interesting discussion with
numerous little insights into the use of special functions on
different spaces.  In particular, they show why associated Legendre
functions are preferred to Chebyshev polynomials for the spherical
harmonics.  They also mention that since the Laplacian measures
curvature, and spherical harmonics of a given degree have the same
Laplacian eigenvalue, they all have the same measure of curvature.
So, for example, the ``ℓ = m`` mode varies most rapidly with
longitude but not at all with latitude, while the ``ℓ = 0`` mode
varies just as rapidly with latitude but not at all with longitude.

[Vasil_2019](@citet) use spin-weighted spherical harmonics to do
tensor calculus in the 3-ball, and have a lot formulas for
derivatives, as a result.


Carrying the general-``n`` calculation above through with degree
``2ℓ`` in ``n=4`` dimensions gives ``2ℓ(2ℓ+2) f / r^2``, so on the unit
3-sphere ``\Delta_{\mathrm{Spin}(3)}\, 𝔇^{(ℓ)}_{m',m} = -4ℓ(ℓ+1)\,
𝔇^{(ℓ)}_{m',m}``.  Comparing with the eigenvalue ``ℓ(ℓ+1)`` of the
Casimir operators, we see that
```math
\Delta_{\mathrm{Spin}(3)} = -4 L^2 = -4 R^2,
```
which is the group-theoretic reason that ``L^2 = R^2``: both are the
Laplacian of the bi-invariant metric, and the Laplacian commutes with
both left and right multiplication.

!!! note "Not conventions"
    Nothing in this section is a convention; it is included only
    because the Laplacian provides a useful cross-check on the
    normalizations above.  The relationship between the Laplacian on
    ``𝕊²`` and ``L^2`` (namely ``\Delta_{𝕊^2} = -L^2`` on
    spin-weight-0 functions) follows from the restriction of the
    result above.
