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


## Three-dimensional space

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

We will be working on the sphere, so it will be very convenient to use
spherical coordinates ``(r, θ, ϕ)``.  We choose the standard
"physics" conventions for these, in which we relate to the Cartesian
coordinates by
```math
\begin{aligned}
r &= \sqrt{x^2 + y^2 + z^2} &&\in [0, \infty), \\
θ &= \arccos\left(\frac{z}{r}\right) &&\in [0, \pi], \\
ϕ &= \arctan\left(\frac{y}{x}\right) &&\in [0, 2\pi),
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
𝐫 &= \sin θ \cos ϕ 𝐱 + \sin θ \sin ϕ 𝐲 + \cos θ 𝐳, \\
\boldsymbol{θ} &= \cos θ \cos ϕ 𝐱 + \cos θ \sin ϕ 𝐲 - \sin θ 𝐳, \\
\boldsymbol{ϕ} &= -\sin ϕ 𝐱 + \cos ϕ 𝐲,
\end{aligned}
```
where, again, we omit the hats on the unit vectors to keep the
notation simple.  Conversely, we can express the Cartesian basis
vectors in terms of the spherical basis vectors as
```math
\begin{aligned}
𝐱 &= \sin θ \cos ϕ 𝐫 + \cos θ \cos ϕ \boldsymbol{θ} - \sin ϕ \boldsymbol{ϕ}, 
\\
𝐲 &= \sin θ \sin ϕ 𝐫 + \cos θ \sin ϕ \boldsymbol{θ} + \cos ϕ \boldsymbol{ϕ},
\\
𝐳 &= \cos θ 𝐫 - \sin θ \boldsymbol{θ}.
\end{aligned}
```

One seemingly obvious — but extremely important — fact is that the
unit basis frame ``(𝐱, 𝐲, 𝐳)`` can be rotated onto
``(\boldsymbol{θ}, \boldsymbol{ϕ}, 𝐫)`` by first
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
\int_{\mathbb{R}^3} f\, d^3𝐫 = \int_0^\infty \int_0^\pi \int_0^{2\pi} f\, r^2 \sin θ\, dr\, dθ\, dϕ.
```
Restricting to the unit sphere, and normalizing so that the integral
of 1 over the sphere is 1, we can simplify this to
```math
\int_{𝕊²} f\, d^2\Omega = \frac{1}{4\pi} \int_0^\pi \int_0^{2\pi} f\, \sin θ\, dθ\, dϕ.
```


## Four-dimensional space: Quaternions and rotations

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

### Quaternions and Euler angles

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
α &= \arctan\frac{Z}{W} + \arctan\frac{-X}{Y} &&\in [0, 2\pi), \\
β &= 2\arccos\sqrt{\frac{W^2+Z^2}{W^2+X^2+Y^2+Z^2}} &&\in [0, 2\pi], \\
γ &= \arctan\frac{Z}{W} - \arctan\frac{-X}{Y} &&\in [0, 2\pi),
\end{aligned}
```
where we again assume the ``\arctan`` in the expressions for
``α`` and ``γ`` is really the two-argument form that gives
the correct quadrant.  Note that here, ``β`` ranges up to ``2\pi``
rather than just ``\pi``, as in the standard Euler angles.  This is
because we are describing the space of quaternions, rather than just
the space of rotations.  If we restrict to ``R=1``, we have exactly
the group of unit quaternions ``\mathrm{Spin}(3)=\mathrm{SU}(2)``,
which is a double cover of the rotation group ``\mathrm{SO}(3)``.
This extended range for ``β`` is necessary to cover the entire
space of quaternions; if we further restrict to ``[0, \pi)``, we would
only cover the space of rotations.  This and the inclusion of ``R``
identify precisely how this coordinate system extends the standard
Euler angles.

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

Again, integration involves a square-root of the determinant of the
metric, which reduces to ``R^3 |\sin β| / 8``.  Note that — unlike
with standard spherical coordinates — the absolute value is necessary
because ``β`` ranges over the entire interval ``[0, 2\pi]``.  The
integral over the entire space of quaternions is then
```math
\int_{\mathbb{R}^4} f\, d^4𝐐
= \int_{-\infty}^\infty \int_{-\infty}^\infty \int_{-\infty}^\infty \int_{-\infty}^\infty f\, dW\, dX\, dY\, dZ
= \int_0^\infty \int_0^{2\pi} \int_0^{2\pi} \int_0^{2\pi} f\, \frac{R^3}{8} |\sin β|\, dR\, dα\, dβ\, dγ.
```
Restricting to the unit sphere, and normalizing so that the integral
of 1 over the sphere is 1, we can simplify this to
```math
\int_{\mathrm{Spin}(3)} f\, d^3\Omega
= \frac{1}{16\pi^2} \int_0^{2\pi} \int_0^{2\pi} \int_0^{2\pi} f\, |\sin β|\, dα\, dβ\, dγ.
```
Finally, restricting to the space of rotations, we can further
simplify this to
```math
\int_{\mathrm{SO}(3)} f\, d^3\Omega
= \frac{1}{8\pi^2} \int_0^{2\pi} \int_0^{\pi} \int_0^{2\pi} f\, \sin β\, dα\, dβ\, dγ.
```

## Rotations

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

### Euler angles and spherical coordinates

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
that plane by the angle ``γ`` about the ``𝐫`` axis.
Thus, we identify the spherical coordinates ``(θ, ϕ)`` with
the Euler angles ``(α, β, γ) = (ϕ, θ, 0)``.


## Rotation and angular-momentum operators

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

### Finite rotations

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
= (\pi/2, 0)`` along the positive ``x`` axis should correspond to the
original field evaluated at the point ``(θ, ϕ) = (\pi/2,
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

### Differential rotations

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

* TODO: Impose ``R_z = s``
* TODO: Impose Condon-Shortley condition (positive, real eigenvalues of ``R_\pm``)
* TODO: Show how the following happens

Using these relations, we can actually solve for the constants
``\lambda`` and ``\rho`` up to a sign.  We find that
```math
\begin{aligned}
\lambda &= i, \\
\rho &= -i.
\end{aligned}
```


### Angular-momentum operators in Euler angles

  - Express angular momentum operators in terms of Euler angles
    - We just rewrite the ``R`` in the Lie definitions in terms of
      Euler angles, multiply by ``\exp(θ/2)``, rederive the new
      Euler angles from that result, and use the chain rule
  - Show for both the three- and two-spheres
  - Show how they act on functions on the three-sphere

The idea here is to express, e.g., $e^{θ 𝐞_i /
2}𝐑_{α, β, γ}$ in quaternion components, then
solve for the new Euler angles $𝐑_{α', β', γ'}$
in terms of the quaternion components, where these new angles all
depend on $θ$.  We then use the chain rule to express
$\partial_θ$ in terms of $\partial_{α'}$, etc., which become
$\partial_α$, etc., when $θ=0$.


```math
\begin{aligned}
  L_i f(𝐑_{α, β, γ})
  &=
  \left. -𝐳 \frac{\partial} {\partial θ} f \left( e^{θ 𝐞_i / 2} 𝐑_{α, β, γ} \right) \right|_{θ=0} \\
  &=
  \left. -𝐳 \frac{\partial} {\partial θ} f \left( 𝐑_{α', β', γ'} \right) \right|_{θ=0} \\
  &=
  \left. -𝐳 \left[ \frac{\partial α'} {\partial θ}\frac{\partial} {\partial α'} + \frac{\partial β'} {\partial θ}\frac{\partial} {\partial β'} + \frac{\partial γ'} {\partial θ}\frac{\partial} {\partial γ'} \right] f \left( 𝐑_{α', β', γ'} \right) \right|_{θ=0} \\
  &=
  -𝐳 \left[ \frac{\partial α'} {\partial θ}\frac{\partial} {\partial α} + \frac{\partial β'} {\partial θ}\frac{\partial} {\partial β} + \frac{\partial γ'} {\partial θ}\frac{\partial} {\partial γ} \right]_{θ=0} f \left( 𝐑_{α, β, γ} \right) \\
  K_i f(𝐑_{α, β, γ})
  &=
  -𝐳 \left[ \frac{\partial α''} {\partial θ}\frac{\partial} {\partial α} + \frac{\partial β''} {\partial θ}\frac{\partial} {\partial β} + \frac{\partial γ''} {\partial θ}\frac{\partial} {\partial γ} \right]_{θ=0} f \left( 𝐑_{α, β, γ} \right),
\end{aligned}
```

```math
\begin{aligned}
𝐑_{α, β, γ}
&=
  R\, \cos\frac{β}{2} \cos\frac{α+γ}{2}
  -R\, \sin\frac{β}{2} \sin\frac{α-γ}{2} 𝐢
  + R\, \sin\frac{β}{2} \cos\frac{α-γ}{2} 𝐣
  + R\, \cos\frac{β}{2} \sin\frac{α+γ}{2} 𝐤.
\\
e^{θ 𝐮 / 2} 𝐑_{α, β, γ}
&= \left(\cos\frac{θ}{2} + 𝐮 \sin\frac{θ}{2}\right) 𝐑_{α, β, γ}
\\
&=
  R\, \cos\frac{θ}{2} \cos\frac{β}{2} \cos\frac{α+γ}{2}
  -R\, \cos\frac{θ}{2} \sin\frac{β}{2} \sin\frac{α-γ}{2} 𝐢
  + R\, \cos\frac{θ}{2} \sin\frac{β}{2} \cos\frac{α-γ}{2} 𝐣
  + R\, \cos\frac{θ}{2} \cos\frac{β}{2} \sin\frac{α+γ}{2} 𝐤
\\
&\quad +
  R\, \sin\frac{θ}{2}\cos\frac{β}{2} \cos\frac{α+γ}{2} 𝐮
  -R\, \sin\frac{θ}{2}\sin\frac{β}{2} \sin\frac{α-γ}{2} 𝐮𝐢
  + R\, \sin\frac{θ}{2}\sin\frac{β}{2} \cos\frac{α-γ}{2} 𝐮𝐣
  + R\, \sin\frac{θ}{2}\cos\frac{β}{2} \sin\frac{α+γ}{2} 𝐮𝐤
\end{aligned}
```

```math
\begin{aligned}
α &= \arctan\frac{Z}{W} + \arctan\frac{-X}{Y} &&\in [0, 2\pi), \\
β &= 2\arccos\sqrt{\frac{W^2+Z^2}{W^2+X^2+Y^2+Z^2}} &&\in [0, 2\pi], \\
γ &= \arctan\frac{Z}{W} - \arctan\frac{-X}{Y} &&\in [0, 2\pi),
\end{aligned}
```


### Laplacians

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
eigenfunction with eigenvalue ``-ℓ(ℓ+2)``.

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

* TODO: Show the relationship between the spherical Laplacian and the
  angular momentum operator.
* TODO: Show how ``D`` matrices are harmonic with respect to the
  Laplacian on the 3-sphere.


## Wigner's 𝔇 matrices

[Sakurai_1994](@citet) says that

> [...] rotations affect physical systems, the state ket corresponding
> to a rotated system is expected to look different from the state ket
> corresponding to the original unrotated system. Given a rotation
> operation ``R``, characterized by a ``3\times 3`` orthogonal matrix
> ``R``, we associate an operator ``\mathscr{D}(R)`` in the
> appropriate ket space such that
> ```math
> |α\rangle_R = \mathscr{D}(R) |α\rangle,
> ```
> ``|α\rangle_R`` and ``|α\rangle`` stand for the kets of
> the rotated and original system, respectively.

If the field is represented as a function ``f(𝐑)``, then rotating the
field by ``e^{ϵ 𝐮/2}`` is equivalent to rotating the argument
of the function by ``e^{-ϵ 𝐮/2}``:
```math
\begin{aligned}
f\left(𝐑\right)
&\to
f\left(e^{-ϵ 𝐮/2}𝐑\right) \\
&\approx
f\left(𝐑\right) + ϵ \left. \frac{d}{dϵ} \right|_{ϵ=0}
f\left(e^{-ϵ 𝐮/2}𝐑\right) \\
&=
f\left(𝐑\right) - i ϵ L_𝐮 f\left(𝐑\right).
\end{aligned}
```
This final expression is precisely equivalent to Sakurai's Eq. (3.1.15):
```math
\mathscr{D}\left(\hat{𝐧}, dϕ \right)
=
1 - i \left( 𝐉 \cdot \hat{𝐧} \right) dϕ.
```

Now, we can write the eigenkets of ``L^2`` and ``L_z`` as ``|ℓ,
m\rangle``, where the eigenvalues are ``ℓ(ℓ+1)`` and ``m``,
respectively.  Finally, define the 𝔇 matrix as (Eq. 3.5.42)
```math
𝔇^{(ℓ)}_{m',m}(R)
=
\langle ℓ, m' | 𝔇(R) | ℓ, m \rangle.
```
Sakurai notes the important result that (Eq. 3.5.46)
```math
𝔇^{(ℓ)}_{m'',m}(R_1\, R_2)
=
\sum_{m'} 𝔇^{(ℓ)}_{m'',m'}(R_1) 𝔇^{(ℓ)}_{m',m}(R_2),
```
and we can readily find the essential behavior with respect to the
first and last Euler angles (Eq. 3.5.50):
```math
\begin{aligned}
𝔇^{(ℓ)}_{m',m}(α, β, γ)
&=
\langle ℓ, m' |
    \exp[-iL_z α]\exp[-iL_y β]\exp[-iL_z γ]
| ℓ, m \rangle \\
&=
\exp[-i(m' α+mγ)]
\langle ℓ, m' | \exp[-iL_y β] | ℓ, m \rangle.
\end{aligned}
```
To belabor this point, recall that in general
```math
\left(\left\langle ψ | A\, B\, C | \chi \right\rangle\right)^\ast
=
\left\langle \chi | C^\dag\, B^\dag\, A^\dag | ψ \right\rangle,
```
and
```math
\left( e^{-i ϵ L_u} \right)^\dag
=
e^{i ϵ L_u^\dag}
=
e^{i ϵ L_u}.
```
Together with the eigenvalue property for the ``L_z`` operator acting
on a ket, this allows us to derive the above result by factoring out
the first and last operators.

Now we are left with the middle operator, which we use to define
```math
\begin{aligned}
d^{(ℓ)}_{m',m}(β)
&=
\langle ℓ, m' | \exp[-iL_y β] | ℓ, m \rangle.
\end{aligned}
```
Using
```math
L_y = (L₊ − L₋) / (2i)
```
we can expand
```math
\exp[-iL_y β]
=
Σ_k (-iL_y β)^k / k!
=
Σ_k (L₋ - L₊)^k (β/2)^k / k!
```

Now, writing ``d_+(X) = [L_+, X]``, Eq. (9) of https://arxiv.org/pdf/1707.03861 says
```math
(L₋ - L₊)^k = \sum_{j=0}^k \binom{k}{j} \left((L₋ - d_+)^j 1\right) (-L₊)^{k-j}
```
The sum will automatically be zero unless ``m+k-j ≤ ℓ`` — which means ``j ≥ m+k-ℓ``
```math
(-L₊)^{k-j}|ℓ,m\rangle = (-1)^{k-j} \sqrt{\frac{(ℓ+m+k-j)!}{(ℓ+m)!}\,\frac{(ℓ-m)!}{(ℓ-m-k+j)!}} |ℓ,m+k-j\rangle
```

``[L₊, L₋] = 2 L_z``

``[L_z, L_\pm] = \pm L_\pm``

I wonder if there's a nicer approach using the symmetry transformation
Edmonds notes in Sec. 4.5 (and credits to Wigner) — or the presumably
equivalent one McEwen and Wiaux use (and credit to Risbo):
```math
\exp\left[ β 𝐣 / 2 \right]
=
\exp\left[ \pi 𝐤 / 4 \right]
\exp\left[ \pi 𝐣 / 4 \right]
\exp\left[ β 𝐤 / 2 \right]
\exp\left[ -\pi 𝐣 / 4 \right]
\exp\left[ -\pi 𝐤 / 4 \right]
```
The 𝔇 matrices corresponding to the ``𝐤`` rotations are simple
phases, which converts the problem into one of finding the 𝔇 matrices
for the ``𝐣`` rotations through angles of ``\pm\pi/2`` — which are
presumably simpler to compute.  See, e.g., Varshalovich's Eq.
4.16.(5), where they are given by purely combinatorial terms.


##  Representation theory / harmonic analysis
  - Representations show up in Fourier analysis on groups
  - Peter-Weyl theorem
    - Generalizes Fourier analysis to compact groups
    - Has three parts, [as given by Wikipedia](https://en.wikipedia.org/wiki/Peter%E2%80%93Weyl_theorem):
      1. "The matrix coefficients of irreducible representations of
         ``G`` are dense in the space ``C(G)`` of continuous
         complex-valued functions on ``G``, and thus also in the space
         ``L^2(G)`` of square-integrable functions."
      2. Unitary representations of ``G`` are completely reducible.
      3. "The regular representation of ``G`` on ``L^2(G)`` decomposes
         as the direct sum of all irreducible unitary representations.
         Moreover, the matrix coefficients of the irreducible unitary
         representations form an orthonormal basis of ``L^2(G)``."
  - Representation theory of ``\mathbf{Spin}(3)``
    - Show how the Lie algebra is represented by the angular-momentum operators
    - Show how the Lie group is represented by the Wigner D-matrices
    - Demonstrate that ``𝔇`` is a representation
    - Demonstrate its behavior under left and right rotation
    - Demonstrate orthonormality
  - Representation theory of ``\mathbf{SO}(3)``
    - There are several places in [Folland](@cite Folland_2016) (e.g.,
      above corollary 5.48) where he mentions that representations of
      a quotient group are just representations that are trivial
      (evidently meaning mapping everything to the identity matrix) on
      the factor.  I can't find anywhere that he explains this
      explicitly, but it seems easy enough to show.  He might do it
      using characters.
    - For ``\mathbf{Spin}(3)`` and ``\mathbf{SO}(3)``, the factor
      group is just ``\{1, -1\}``.  Presumably, every representation
      acting on ``1`` will give the identity matrix, so that's
      trivial.  So we just need a criterion for when a representation
      is trivial on ``-1``.  Noting that ``\exp(\pi \vec{v}) = -1``
      for any ``\vec{v}``, I think we can show that this requires
      ``m \in \mathbb{Z}``.
    - Basically, the point is that the representations of
      ``\mathbf{SO}(3)`` are just the integer representations of
      ``\mathbf{Spin}(3)``.
  - Restrict to homogeneous space (S³ -> S²)
    - The circle group is a closed (normal?) subgroup of
      ``\mathbf{Spin}(3)``, which we might implement as initial
      multiplication about a particular axis.
    - In Eq. (2.47) [Folland (2016)](@cite Folland_2016) defines a
      functional taking a function on the group to a function on the
      homogeneous space by integrating over the factor (the circle
      group).  This gives you the spherical harmonics, but *not* the
      spin-weighted spherical harmonics — because the spin-weighted
      spherical harmonics cannot be defined on the 2-sphere.
    - Spin weight comes from Fourier analysis on the subgroup.
    - Representation matrices transfer to the homogeneous space, with
      sparsity patterns

Theorem 2.16 of [Hanson-Yakovlev](@cite HansonYakovlev_2002) says that
an orthonormal basis of a product of ``L^2`` spaces is given by the
product of the orthonormal bases of the individual spaces.
Furthermore, on page 354, they point out that ``\{(1/\sqrt{2\pi})
e^{imϕ}\}`` is an orthonormal basis of ``L^2(0,2\pi)``, while the
set ``\{1/c_{n,m} P_n^m(\cos θ)`` is an orthonormal basis of
``L^2(0, \pi)`` in the ``θ`` coordinate.  Therefore, the product
of these two sets is an orthonormal basis of the product space
``L^2\left((0,2\pi) \times (0, \pi)\right)``, which forms a coordinate
space for ``𝕊²``.  I would probably modify this to point out that
``(0,2\pi)`` is really ``𝕊¹``, and then we could extend it to point
out that you can throw on another factor of ``𝕊¹`` to cover ``𝕊³``,
which happens to give us the Wigner D-matrices.

## Recursion relations

[Gumerov and Duraiswami (2001)](@cite Gumerov_2001) derive their
recursion relations by differentiating solutions of the Helmholtz
equation ``\nabla^2 ψ + k^2 ψ = 0`` as ``\tfrac{1}{k} \nabla
ψ``.  More precisely, they differentiate both sides of the equation
relating one solution to its rotated form — which naturally involves
Wigner's ``𝔇`` matrix.  Using orthogonal basis functions
for the solution, this allows them to equate terms on the two sides
proportional to a given basis function, which leaves them with
expressions involving sums of only the ``𝔇`` matrices and
some coefficients depending on the indices of the basis functions (and
hence of ``𝔇``) on both sides of the equation.  Since
``\nabla`` is a 3-vector operator, this gives them three relations.

This, of course, is happening in 3-D space, since ``ψ`` is a
function of location in the Helmholtz equation.  It seems likely to
me, however, that we could use the 4-D (quaternionic) version of the
functions.  Note that G&D use ``\partial_z`` and ``\partial_x \pm i
\partial_y`` as their operators to differentiate the functions — that
is, the derivatives are with respect to Cartesian coordinates, which
may be more similar to the right-derivative defined above.  However, I
don't know that we'll necessarily be able to achieve the same results
with just angular-momentum operators, since their operators do involve
moving off of the sphere.  Maybe we'd need to move off of the sphere
in 4-D space to get comparable results.  Or maybe just use something
like ``𝐫 ∧ L``, which should also have 3 degrees of freedom.

The SWSHs/``𝔇`` functions can be naturally promoted to
functions not just on the 3-sphere, but also in 4-D space just by
allowing the quaternions to be non-unit quaternions.
