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
\mathbf{v} = v_x \mathbf{𝐱} + v_y \mathbf{𝐲} + v_z \mathbf{𝐳},
```
in which case the Euclidean norm is given by
```math
\| \mathbf{v} \| = \sqrt{v_x^2 + v_y^2 + v_z^2}.
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
spherical coordinates ``(r, \theta, \phi)``.  We choose the standard
"physics" conventions for these, in which we relate to the Cartesian
coordinates by
```math
\begin{aligned}
r &= \sqrt{x^2 + y^2 + z^2} &&\in [0, \infty), \\
\theta &= \arccos\left(\frac{z}{r}\right) &&\in [0, \pi], \\
\phi &= \arctan\left(\frac{y}{x}\right) &&\in [0, 2\pi),
\end{aligned}
```
where we assume the ``\arctan`` in the expression for ``\phi`` is
really the two-argument form that gives the correct quadrant.  The
inverse transformation is given by
```math
\begin{aligned}
x &= r \sin\theta \cos\phi, \\
y &= r \sin\theta \sin\phi, \\
z &= r \cos\theta.
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
  0 & 0 & r^2 \sin^2\theta
\end{array} \right)_{i'j'}.
```
The unit coordinate vectors in spherical coordinates are then
```math
\begin{aligned}
\mathbf{𝐫} &= \sin\theta \cos\phi \mathbf{𝐱} + \sin\theta \sin\phi \mathbf{𝐲} + \cos\theta \mathbf{𝐳}, \\
\boldsymbol{\theta} &= \cos\theta \cos\phi \mathbf{𝐱} + \cos\theta \sin\phi \mathbf{𝐲} - \sin\theta \mathbf{𝐳}, \\
\boldsymbol{\phi} &= -\sin\phi \mathbf{𝐱} + \cos\phi \mathbf{𝐲},
\end{aligned}
```
where, again, we omit the hats on the unit vectors to keep the
notation simple.  Conversely, we can express the Cartesian basis
vectors in terms of the spherical basis vectors as
```math
\begin{aligned}
\mathbf{𝐱} &= \sin\theta \cos\phi \mathbf{𝐫} + \cos\theta \cos\phi \boldsymbol{\theta} - \sin\phi \boldsymbol{\phi}, 
\\
\mathbf{𝐲} &= \sin\theta \sin\phi \mathbf{𝐫} + \cos\theta \sin\phi \boldsymbol{\theta} + \cos\phi \boldsymbol{\phi},
\\
\mathbf{𝐳} &= \cos\theta \mathbf{𝐫} - \sin\theta \boldsymbol{\theta}.
\end{aligned}
```

One seemingly obvious — but extremely important — fact is that the
unit basis frame ``(𝐱, 𝐲, 𝐳)`` can be rotated onto
``(\boldsymbol{\theta}, \boldsymbol{\phi}, \mathbf{r})`` by first
rotating through the "polar" angle ``\theta`` about the ``\mathbf{y}``
axis, and then through the "azimuthal" angle ``\phi`` about the
``\mathbf{z}`` axis.  This becomes important when we consider
spin-weighted functions.

Integration in Cartesian coordinates is, of course, trivial as
```math
\int_{\mathbb{R}^3} f\, d^3\mathbf{r} = \int_{-\infty}^{\infty} \int_{-\infty}^{\infty} \int_{-\infty}^{\infty} f\, dx\, dy\, dz.
```
In spherical coordinates, the integrand involves the square-root of
the determinant of the metric, so we have
```math
\int_{\mathbb{R}^3} f\, d^3\mathbf{r} = \int_0^\infty \int_0^\pi \int_0^{2\pi} f\, r^2 \sin\theta\, dr\, d\theta\, d\phi.
```
Restricting to the unit sphere, and normalizing so that the integral
of 1 over the sphere is 1, we can simplify this to
```math
\int_{S^2} f\, d^2\Omega = \frac{1}{4\pi} \int_0^\pi \int_0^{2\pi} f\, \sin\theta\, d\theta\, d\phi.
```


## Four-dimensional space: Quaternions and rotations

### Geometric algebra

Given the basis vectors ``(𝐱, 𝐲, 𝐳)`` and the Euclidean norm, we
can define the [geometric
algebra](https://en.wikipedia.org/wiki/Geometric_algebra).  The key
feature is the geometric product, which is defined for any pair of
vectors as ``𝐯`` and ``𝐰`` as
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
\alpha &= \arctan\frac{Z}{W} + \arctan\frac{-X}{Y} &&\in [0, 2\pi), \\
\beta &= 2\arccos\sqrt{\frac{W^2+Z^2}{W^2+X^2+Y^2+Z^2}} &&\in [0, 2\pi], \\
\gamma &= \arctan\frac{Z}{W} - \arctan\frac{-X}{Y} &&\in [0, 2\pi),
\end{aligned}
```
where we again assume the ``\arctan`` in the expressions for
``\alpha`` and ``\gamma`` is really the two-argument form that gives
the correct quadrant.  Note that here, ``\beta`` ranges up to ``2\pi``
rather than just ``\pi``, as in the standard Euler angles.  This is
because we are describing the space of quaternions, rather than just
the space of rotations.  If we restrict to ``R=1``, we have exactly
the group of unit quaternions ``\mathrm{Spin}(3)=\mathrm{SU}(2)``,
which is a double cover of the rotation group ``\mathrm{SO}(3)``.
This extended range for ``\beta`` is necessary to cover the entire
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
  0 & \frac{R^2}{4} & 0 & \frac{R^2 \cos\beta}{4} \\
  0 & 0 & \frac{R^2}{4} & 0 \\
  0 & \frac{R^2 \cos\beta}{4} & 0 & \frac{R^2}{4}
\end{array} \right)_{i'j'}.
```
The unit basis vectors in extended Euler coordinates in terms of the
unit basis vectors in quaternion coordinates are
```math
\begin{aligned}
\mathbf{𝐑} &= \frac{1}{R} \left(
  \cos \frac{\beta}{2} \cos \frac{\alpha+\gamma}{2} 𝟏
  - \sin \frac{\beta}{2} \sin \frac{\alpha-\gamma}{2} 𝐢
  + \sin \frac{\beta}{2} \cos \frac{\alpha-\gamma}{2} 𝐣
  + \cos \frac{\beta}{2} \sin \frac{\alpha+\gamma}{2} 𝐤
\right), \\
\boldsymbol{\alpha} &= \frac{R}{2} \left(
  -\cos \frac{\beta}{2} \sin \frac{\alpha+\gamma}{2} 𝟏
  - \sin \frac{\beta}{2} \cos \frac{\alpha-\gamma}{2} 𝐢
  - \sin \frac{\beta}{2} \sin \frac{\alpha-\gamma}{2} 𝐣
  + \cos \frac{\beta}{2} \cos \frac{\alpha+\gamma}{2} 𝐤
\right), \\
\boldsymbol{\beta} &= \frac{R}{2} \left(
  -\sin \frac{\beta}{2} \cos \frac{\alpha+\gamma}{2} 𝟏
  - \cos \frac{\beta}{2} \sin \frac{\alpha-\gamma}{2} 𝐢
  + \cos \frac{\beta}{2} \cos \frac{\alpha-\gamma}{2} 𝐣
  - \sin \frac{\beta}{2} \sin \frac{\alpha+\gamma}{2} 𝐤
\right), \\
\boldsymbol{\gamma} &= \frac{R}{2} \left(
  -\cos \frac{\beta}{2} \sin \frac{\alpha+\gamma}{2} 𝟏
  + \sin \frac{\beta}{2} \cos \frac{\alpha-\gamma}{2} 𝐢
  - \sin \frac{\beta}{2} \cos \frac{\alpha-\gamma}{2} 𝐣
  - \cos \frac{\beta}{2} \sin \frac{\alpha+\gamma}{2} 𝐤
\right).
\end{aligned}
```

Again, integration involves a square-root of the determinant of the
metric, which reduces to ``R^3 |\sin\beta| / 8``.  Note that — unlike
with standard spherical coordinates — the absolute value is necessary
because ``\beta`` ranges over the entire interval ``[0, 2\pi]``.  The
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
= \exp\left(\frac{\rho}{2} \hat{\mathfrak{r}}\right)
= \cos\frac{\rho}{2} + \sin\frac{\rho}{2}\, \hat{\mathfrak{r}},
```
where ``\rho`` is an angle of rotation and ``\hat{\mathfrak{r}}`` is a
unit "pure-vector" quaternion.  We can multiply a vector ``𝐯`` as
```math
𝐑\, 𝐯\, 𝐑^{-1}.
```
Splitting ``𝐯 = 𝐯_⟂ + 𝐯_∥`` into components perpendicular and
parallel to ``\hat{\mathfrak{r}}``, we see that ``𝐯_∥`` commutes with
``𝐑`` and ``𝐑^{-1}``, while ``𝐯_⟂`` anticommutes with
``\hat{\mathfrak{r}}``.  To find the full rotation, we expand the
product:
```math
\begin{aligned}
𝐑\, 𝐯\, 𝐑^{-1}
&= 𝐯_∥
   + \left(\cos\frac{\rho}{2} + \sin\frac{\rho}{2}\, \hat{\mathfrak{r}}\right)
     𝐯_⟂
     \left(\cos\frac{\rho}{2} - \sin\frac{\rho}{2}\, \hat{\mathfrak{r}}\right) \\
&= 𝐯_∥
   + \left(\cos\frac{\rho}{2}\, 𝐯_⟂ + \sin\frac{\rho}{2}\, \hat{\mathfrak{r}}\, 𝐯_⟂\right)
     \left(\cos\frac{\rho}{2} - \sin\frac{\rho}{2}\, \hat{\mathfrak{r}}\right) \\
&= 𝐯_∥
   + \cos^2\frac{\rho}{2}\, 𝐯_⟂ + \sin\frac{\rho}{2}\, \cos\frac{\rho}{2}\, \hat{\mathfrak{r}}\, 𝐯_⟂
   - \sin\frac{\rho}{2}\, \cos\frac{\rho}{2}\, 𝐯_⟂ \, \hat{\mathfrak{r}} - \sin^2\frac{\rho}{2}\, \hat{\mathfrak{r}}\, 𝐯_⟂\, \hat{\mathfrak{r}} \\
&= 𝐯_∥
   + \cos^2\frac{\rho}{2}\, 𝐯_⟂ + \sin\frac{\rho}{2}\, \cos\frac{\rho}{2}\, [\hat{\mathfrak{r}}, 𝐯_⟂] - \sin^2\frac{\rho}{2}\, 𝐯_⟂ \\
&= 𝐯_∥
   + \cos\rho\, 𝐯_⟂ + \sin\rho\, \hat{\mathfrak{r}}\times 𝐯_⟂
\end{aligned}
```
The final expression shows that this is precisely what we expect when
rotating ``𝐯`` through an angle ``\rho`` (in a positive, right-handed
sense) about the axis ``\hat{\mathfrak{r}}``.

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
``\gamma`` about the ``𝐳`` axis, followed by a rotation through
``\beta`` about the ``𝐲`` axis, and finally a rotation through
``\alpha`` about the ``𝐳`` axis.  Note that the axes are fixed, and
not subject to any preceding rotations.  More precisely, we can write
the unit quaternion as
```math
𝐑 = \exp\left(\frac{\alpha}{2} 𝐤\right)
    \exp\left(\frac{\beta}{2} 𝐣\right)
    \exp\left(\frac{\gamma}{2} 𝐤\right).
```
One of the more important interpretations of a rotor is considering
what it does to the basis triad ``(𝐱, 𝐲, 𝐳)``.  In particular, the
vector ``𝐳`` is rotated onto the point given by spherical coordinates
``(\theta, \phi) = (\beta, \alpha)``, while ``𝐱`` and ``𝐲`` are
rotated into the plane spanned by the unit basis vectors
``\boldsymbol{\theta}`` and ``\boldsymbol{\phi}`` corresponding to
that point.  If ``\gamma = 0`` the rotation is precise, meaning that
``𝐱`` is rotated onto ``\boldsymbol{\theta}`` and ``𝐲`` onto
``\boldsymbol{\phi}``; if ``\gamma ≠ 0`` then they are rotated within
that plane by the angle ``\gamma`` about the ``\mathbf{r}`` axis.
Thus, we identify the spherical coordinates ``(\theta, \phi)`` with
the Euler angles ``(\alpha, \beta, \gamma) = (\phi, \theta, 0)``.


## Rotation and angular-momentum operators

### Complex-valued functions

Starting with Cartesian coordinates and the Euclidean norm on
``\mathbb{R}^3``, we have *constructed* the geometric algebra over
that space, as well as the spaces ``\mathrm{Spin}(3) =
\mathrm{SU}(2)`` (topologically ``S^3``), ``\mathrm{SO}(3)``
(topologically ``\mathbb{RP}^3``), and ``S^2``.  We will be defining
complex-valued functions on these spaces, and defining operators to
construct and classify them.  In particular, because we have
constructed the spaces, they are naturally supplied with coordinates
that are effectively inherited from the original Cartesian system.  We
will be using these coordinate systems to construct both the operators
and functions.  However, it is important to note that the coordinate
systems may have singularities, which means that the spaces of
coordinates may have different topologies than the spaces they
represent.  For example, Euler angles have topology ``S^1 \times I
\times S^1`` instead of the ``S^3`` and ``\mathbb{RP}^3`` topologies
of the spaces they represent; spherical coordinates have topology
``S^1 \times I`` instead of ``S^2``.

Defining functions on the coordinate system of a space is subtly
different from defining functions on the space itself.  For example,
spin-weighted functions are generally written as functions of
(``S^2``) spherical coordinates.  However, they *cannot* be defined as
functions on ``S^2`` itself; some notion of a reference tangent
direction is needed at each point.  The difference is that spherical
*coordinates* supply a natural choice for the reference tangent
direction: the unit vector in the ``\boldsymbol{\theta}`` direction.
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
Any non-zero quaternion can be expressed as ``e^\mathfrak{g}`` for
some finite quaternion ``\mathfrak{g}``, which is referred to as the
"generator" of the action of ``e^\mathfrak{g}``.  This can act on a
function ``f`` by multiplying the argument by ``e^\mathfrak{g}``.
However, there is an ambiguity: we could multiply either on the left
or the right:[^2]
```math
f\left(\mathbf{Q}\right) \mapsto f\left(e^\mathfrak{g} \mathbf{Q}\right)
\qquad \text{or} \qquad
f\left(\mathbf{Q}\right) \mapsto f\left(\mathbf{Q} e^\mathfrak{g}\right).
```
There is an additional ambiguity, in that this action rotates the
*argument* of the function, whereas we will often prefer to think in
terms of rotating the *function* itself.  For example, our function
may describe the measurement of some field in a particular coordinate
system.  Here, the argument ``\mathbf{Q}`` describes a particular
value of the coordinates, and ``e^\mathfrak{g}`` changes the point
under consideration.  If, on the other hand ``e^\mathfrak{g}``
describes how the field itself is rotated, then we can write the
rotated field as a function ``f'`` which is related to the original
function ``f`` by
```math
f'\left(\mathbf{Q}\right) = f\left(e^{-\mathfrak{g}} \mathbf{Q}\right)
\qquad \text{or} \qquad
f'\left(\mathbf{Q}\right) = f\left(\mathbf{Q} e^{-\mathfrak{g}}\right).
```
Note that the exponent is negated, because the action of
``e^\mathfrak{g}`` on the argument is the inverse of the action of
``e^{-\mathfrak{g}}`` on the function.  This is a general property of
the action of a group on a space, and is a consequence of the group
action being a homomorphism.

[^2]: In group theory, this type of transformation is often referred
      to as a "translation", even when — as in this case — we would
      usually describe these as rotations.

To validate the signs here, it may be helpful to work through a simple
example involving the sphere ``S^2``.  We define a function on
spherical coordinates as
```math
f(\theta, \phi) = \sin\theta \sin\phi.
```
Recall that we can map the spherical coordinates into the Euler
angles, and the Euler angles into the quaternion
```math
(\theta, \phi) \mapsto (\phi, \theta, 0) \mapsto \mathbf{Q}
=
\exp\left(\frac{\phi}{2} \mathbf{k}\right)
\exp\left(\frac{\theta}{2} \mathbf{j}\right).
```
It is straightforward to see that we can write ``f`` as a function of
``\mathbf{Q}`` as
```math
f(\mathbf{Q}) = \left\langle \mathbf{Q}\, \mathbf{k}\, \mathbf{Q}^{-1} \right\rangle_{\mathbf{j}},
```
where the angle brackets and subscript indicate that we are taking the
``\mathbf{j}`` component.  That is, ``f`` is the ``y`` component of
the vector ``\mathbf{z}`` rotated by ``\mathbf{Q}``. 

Now, we imagine rotating the field by an angle ``\alpha`` in the
positive sense about the ``z`` axis.  Visualizing the situation, we
can see that the rotated field should be represented by
```math
f'(\theta, \phi) = \sin\theta \sin(\phi - \alpha).
```
For example, the rotated field evaluated at the point ``(\theta, \phi)
= (\pi/2, 0)`` along the positive ``x`` axis should correspond to the
original field evaluated at the point ``(\theta, \phi) = (\pi/2,
-\alpha)``.  This rotation is generated by ``\mathfrak{g} = \alpha
\mathbf{k} / 2``, which allows us to immediately calculate
```math
\begin{aligned}
f(e^\mathfrak{g} \mathbf{Q}) &= \sin\theta \sin(\phi + \alpha) &&&
f(\mathbf{Q} e^\mathfrak{g}) &= \sin\theta \sin\phi \\
f(e^{-\mathfrak{g}} \mathbf{Q}) &= \sin\theta \sin(\phi - \alpha) &&&
f(\mathbf{Q} e^{-\mathfrak{g}}) &= \sin\theta \sin\phi.
\end{aligned}
```
Thus, we see that left-multiplication by ``e^{-\mathfrak{g}}``
corresponds to rotation of the field while leaving the coordinates
fixed; left-multiplication by ``e^\mathfrak{g}`` corresponds to
rotation of the coordinates while leaving the field fixed; and
right-multiplication by either doesn't affect this function at all.

Of course, right-multiplication using other choices for
``\mathfrak{g}`` could certainly have some effect on this function,
and this choice of ``\mathfrak{g}`` could have an effect on other
functions.  Note that right-multiplication can also be interpreted as
left-multiplication, where the generator itself is rotated by the
argument to the function.  That is,
```math
\begin{aligned}
f(\mathbf{Q} e^\mathfrak{g})
  &= f(\mathbf{Q} e^{\mathfrak{g}} \mathbf{Q}^{-1} \mathbf{Q})
  = f(e^{\mathfrak{g}'} \mathbf{Q}) \\
f(\mathbf{Q} e^{-\mathfrak{g}})
  &= f(\mathbf{Q} e^{-\mathfrak{g}} \mathbf{Q}^{-1} \mathbf{Q})
  = f(e^{-\mathfrak{g}'} \mathbf{Q}),
\end{aligned}
```
where ``\mathfrak{g}' = \mathbf{Q} \mathfrak{g} \mathbf{Q}^{-1}``.  In
this example, ``\mathfrak{g}'`` generates a rotation by an angle
``\alpha`` about the point in question, which leaves that point fixed,
and since this is a scalar function it has no effect on the value.  Of
course, we will see below that changing by a phase proportional to
``\alpha`` is the defining feature of a *spin-weighted* function.

### Differential rotations

We now define a pair of operators that differentiate a function with
respect to infinitesimal rotations we apply to the functions
themselves:
```math
\begin{aligned}
L_{\mathfrak{g}} f(\mathbf{Q}) &= \lambda \left. \frac{\partial} {\partial \theta} f \left( e^{-\theta \mathfrak{g} / 2} \mathbf{Q} \right) \right|_{\theta=0}, \\
R_{\mathfrak{g}} f(\mathbf{Q}) &= \rho \left. \frac{\partial} {\partial \theta} f \left( \mathbf{Q} e^{-\theta \mathfrak{g} / 2} \right) \right|_{\theta=0}.
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
  L_\mathfrak{g} L_\mathfrak{h} f(\mathbf{Q})
  % &= \left. \lambda \frac{\partial} {\partial \gamma} f'\left(e^{-\gamma \mathfrak{g} / 2} \mathbf{Q} \right) \right|_{\gamma=0}, \\
  &= \left. \lambda^2 \frac{\partial} {\partial \gamma} \frac{\partial} {\partial \eta} f\left(e^{-\eta \mathfrak{h} / 2} e^{-\gamma \mathfrak{g} / 2} \mathbf{Q} \right) \right|_{\gamma=\eta=0}, \\
  R_\mathfrak{g} R_\mathfrak{h} f(\mathbf{Q})
  % &= \rho \left. \frac{\partial} {\partial \gamma} f' \left( \mathbf{Q} e^{-\gamma \mathfrak{g} / 2} \right) \right|_{\gamma=0} \\
  &= \left. \rho^2 \frac{\partial} {\partial \gamma} \frac{\partial} {\partial \eta} f\left( \mathbf{Q} e^{-\gamma \mathfrak{g} / 2} e^{-\eta \mathfrak{h} / 2} \right) \right|_{\gamma=\eta=0}.
\end{aligned}
```
We can prove the first of these, for example, by defining
``f'(\mathbf{Q}) = L_\mathfrak{h} f(\mathbf{Q})``, then applying the
definition of ``L_\mathfrak{g}`` to ``f'(\mathbf{Q})``, and finally
substituting the definition of ``f'`` back in.  If we failed to use
the correct order of operations, we would get sign errors when trying
to evaluate the commutators.

These operators have some nice properties.  For any scalar ``s``, we have
```math
\begin{aligned}
L_{s \mathfrak{g}} &= s L_{\mathfrak{g}}, \\
R_{s \mathfrak{g}} &= s R_{\mathfrak{g}}.
\end{aligned}
```
Given any basis ``\mathbf{e}_n`` for the quaternions, we can use
the multivariable chain rule to expand the operators in terms of
components:
```math
\begin{aligned}
L_{\mathfrak{g}} &= \sum_n g_n\, L_{\mathbf{e}_n}, \\
R_{\mathfrak{g}} &= \sum_n g_n\, R_{\mathbf{e}_n}.
\end{aligned}
```
This implies that vector addition holds more generally:
```math
\begin{aligned}
L_{\mathfrak{g} + \mathfrak{h}} &= L_{\mathfrak{g}} + L_{\mathfrak{h}} \\
R_{\mathfrak{g} + \mathfrak{h}} &= R_{\mathfrak{g}} + R_{\mathfrak{h}}.
\end{aligned}
```
Moreover, we can show that these operators form a Lie algebra with the
commutator as the Lie bracket.  That is, we have
```math
\begin{aligned}
[L_{\mathfrak{g}}, L_{\mathfrak{h}}]
    &= \frac{\lambda}{2} L_{[\mathfrak{g}, \mathfrak{h}]},
\\
[R_{\mathfrak{g}}, R_{\mathfrak{h}}]
    &= -\frac{\rho}{2} R_{[\mathfrak{g}, \mathfrak{h}]},
\\
[L_{\mathfrak{g}}, R_{\mathfrak{h}}] &= 0.
\end{aligned}
```

Conventionally, we single out the ``\mathbf{z}`` axis — or
equivalently the generator ``\mathbf{k} = \mathbf{y}\mathbf{x}`` — as
a sort of fiducial axis, and ``L_z = L_\mathbf{k}`` and ``R_z =
R_\mathbf{k}`` as the fiducial operators.  Then, *by definition*,
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
they have no component proportional ``L_\mathbf{z}``, and that both of
the remaining components must be nonzero.  This actually allows us to
deduce that ``\lambda^2 = \rho^2 = -1``.  This, in turn, allows us to
deduce the values of the raising and lowering operators up to an
overall factor.  Conventionally the factor is chosen so that
```math
\begin{aligned}
L_\pm &= L_\mathbf{x} \pm i L_\mathbf{y}, \\
R_\pm &= R_\mathbf{x} \pm i R_\mathbf{y}.
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
      Euler angles, multiply by ``\exp(\theta/2)``, rederive the new
      Euler angles from that result, and use the chain rule
  - Show for both the three- and two-spheres
  - Show how they act on functions on the three-sphere

The idea here is to express, e.g., $e^{\theta \mathbf{e}_i /
2}\mathbf{R}_{\alpha, \beta, \gamma}$ in quaternion components, then
solve for the new Euler angles $\mathbf{R}_{\alpha', \beta', \gamma'}$
in terms of the quaternion components, where these new angles all
depend on $\theta$.  We then use the chain rule to express
$\partial_\theta$ in terms of $\partial_{\alpha'}$, etc., which become
$\partial_\alpha$, etc., when $\theta=0$.


```math
\begin{aligned}
  L_i f(\mathbf{R}_{\alpha, \beta, \gamma})
  &=
  \left. -\mathbf{z} \frac{\partial} {\partial \theta} f \left( e^{\theta \mathbf{e}_i / 2} \mathbf{R}_{\alpha, \beta, \gamma} \right) \right|_{\theta=0} \\
  &=
  \left. -\mathbf{z} \frac{\partial} {\partial \theta} f \left( \mathbf{R}_{\alpha', \beta', \gamma'} \right) \right|_{\theta=0} \\
  &=
  \left. -\mathbf{z} \left[ \frac{\partial \alpha'} {\partial \theta}\frac{\partial} {\partial \alpha'} + \frac{\partial \beta'} {\partial \theta}\frac{\partial} {\partial \beta'} + \frac{\partial \gamma'} {\partial \theta}\frac{\partial} {\partial \gamma'} \right] f \left( \mathbf{R}_{\alpha', \beta', \gamma'} \right) \right|_{\theta=0} \\
  &=
  -\mathbf{z} \left[ \frac{\partial \alpha'} {\partial \theta}\frac{\partial} {\partial \alpha} + \frac{\partial \beta'} {\partial \theta}\frac{\partial} {\partial \beta} + \frac{\partial \gamma'} {\partial \theta}\frac{\partial} {\partial \gamma} \right]_{\theta=0} f \left( \mathbf{R}_{\alpha, \beta, \gamma} \right) \\
  K_i f(\mathbf{R}_{\alpha, \beta, \gamma})
  &=
  -\mathbf{z} \left[ \frac{\partial \alpha''} {\partial \theta}\frac{\partial} {\partial \alpha} + \frac{\partial \beta''} {\partial \theta}\frac{\partial} {\partial \beta} + \frac{\partial \gamma''} {\partial \theta}\frac{\partial} {\partial \gamma} \right]_{\theta=0} f \left( \mathbf{R}_{\alpha, \beta, \gamma} \right),
\end{aligned}
```

```math
\begin{aligned}
\mathbf{R}_{\alpha, \beta, \gamma}
&=
  R\, \cos\frac{β}{2} \cos\frac{α+γ}{2}
  -R\, \sin\frac{β}{2} \sin\frac{α-γ}{2} \mathbf{i}
  + R\, \sin\frac{β}{2} \cos\frac{α-γ}{2} \mathbf{j}
  + R\, \cos\frac{β}{2} \sin\frac{α+γ}{2} \mathbf{k}.
\\
e^{\theta \mathbf{u} / 2} \mathbf{R}_{\alpha, \beta, \gamma}
&= \left(\cos\frac{\theta}{2} + \mathbf{u} \sin\frac{\theta}{2}\right) \mathbf{R}_{\alpha, \beta, \gamma}
\\
&=
  R\, \cos\frac{\theta}{2} \cos\frac{β}{2} \cos\frac{α+γ}{2}
  -R\, \cos\frac{\theta}{2} \sin\frac{β}{2} \sin\frac{α-γ}{2} \mathbf{i}
  + R\, \cos\frac{\theta}{2} \sin\frac{β}{2} \cos\frac{α-γ}{2} \mathbf{j}
  + R\, \cos\frac{\theta}{2} \cos\frac{β}{2} \sin\frac{α+γ}{2} \mathbf{k}
\\
&\quad +
  R\, \sin\frac{\theta}{2}\cos\frac{β}{2} \cos\frac{α+γ}{2} \mathbf{u}
  -R\, \sin\frac{\theta}{2}\sin\frac{β}{2} \sin\frac{α-γ}{2} \mathbf{u}\mathbf{i}
  + R\, \sin\frac{\theta}{2}\sin\frac{β}{2} \cos\frac{α-γ}{2} \mathbf{u}\mathbf{j}
  + R\, \sin\frac{\theta}{2}\cos\frac{β}{2} \sin\frac{α+γ}{2} \mathbf{u}\mathbf{k}
\end{aligned}
```

```math
\begin{aligned}
\alpha &= \arctan\frac{Z}{W} + \arctan\frac{-X}{Y} &&\in [0, 2\pi), \\
\beta &= 2\arccos\sqrt{\frac{W^2+Z^2}{W^2+X^2+Y^2+Z^2}} &&\in [0, 2\pi], \\
\gamma &= \arctan\frac{Z}{W} - \arctan\frac{-X}{Y} &&\in [0, 2\pi),
\end{aligned}
```


### Laplacians

[Bander_1966](@citet) show that Wigner's D matrices (extended to the
full space of quaternions with arbitrary norm) are harmonic with
respect to the Laplacian of the full 4-D space.  We also know that
```math
\Delta_{S^{n-1}} f(x) = \Delta_{\mathbb{R}^n} f(x/|x|),
```
and
```math
\Delta_{\mathbb{R}^n} f(x)
=
\frac{1}{r^{n-1}} \frac{\partial}{\partial r} \left( r^{n-1} \frac{\partial f}{\partial r} \right)
+
\frac{1}{r^2} \Delta_{S^{n-1}} f.
```
These imply that the restriction to the space of unit quaternions is
not harmonic with respect to the Laplacian on the 3-sphere, but is an
eigenfunction with eigenvalue ``-\ell(\ell+2)``.

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
\frac{f}{r^\ell} \frac{\partial^2 r^\ell}{\partial r^2}
+
\frac{f}{r^\ell} \frac{n-1}{r} \frac{\partial r^\ell}{\partial r}
=
\ell(\ell-1) \frac{f}{r^\ell} r^{\ell-2}
+
\ell \frac{f}{r^\ell} \frac{n-1}{r} r^{\ell-1}
=
\ell(\ell-1) \frac{f}{r^2}
+
\ell (n-1) \frac{f}{r^2}
=
\ell(\ell+n-2) \frac{f}{r^2}
\to
\ell(\ell+2) \frac{f}{r^2}
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
> |\alpha\rangle_R = \mathscr{D}(R) |\alpha\rangle,
> ```
> ``|\alpha\rangle_R`` and ``|\alpha\rangle`` stand for the kets of
> the rotated and original system, respectively.

If the field is represented as a function ``f(𝐑)``, then rotating the
field by ``e^{\epsilon 𝐮/2}`` is equivalent to rotating the argument
of the function by ``e^{-\epsilon 𝐮/2}``:
```math
\begin{aligned}
f\left(𝐑\right)
&\to
f\left(e^{-\epsilon 𝐮/2}𝐑\right) \\
&\approx
f\left(𝐑\right) + \epsilon \left. \frac{d}{d\epsilon} \right|_{\epsilon=0}
f\left(e^{-\epsilon 𝐮/2}𝐑\right) \\
&=
f\left(𝐑\right) - i \epsilon L_𝐮 f\left(𝐑\right).
\end{aligned}
```
This final expression is precisely equivalent to Sakurai's Eq. (3.1.15):
```math
\mathscr{D}\left(\hat{\mathbf{n}}, d\phi \right)
=
1 - i \left( \mathbf{J} \cdot \hat{\mathbf{n}} \right) d\phi.
```

Now, we can write the eigenkets of ``L^2`` and ``L_z`` as ``|\ell,
m\rangle``, where the eigenvalues are ``\ell(\ell+1)`` and ``m``,
respectively.  Finally, define the 𝔇 matrix as (Eq. 3.5.42)
```math
𝔇^{(\ell)}_{m',m}(R)
=
\langle \ell, m' | 𝔇(R) | \ell, m \rangle.
```
Sakurai notes the important result that (Eq. 3.5.46)
```math
𝔇^{(\ell)}_{m'',m}(R_1\, R_2)
=
\sum_{m'} 𝔇^{(\ell)}_{m'',m'}(R_1) 𝔇^{(\ell)}_{m',m}(R_2),
```
and we can readily find the essential behavior with respect to the
first and last Euler angles (Eq. 3.5.50):
```math
\begin{aligned}
𝔇^{(\ell)}_{m',m}(\alpha, \beta, \gamma)
&=
\langle \ell, m' |
    \exp[-iL_z \alpha]\exp[-iL_y \beta]\exp[-iL_z \gamma]
| \ell, m \rangle \\
&=
\exp[-i(m' \alpha+m\gamma)]
\langle \ell, m' | \exp[-iL_y \beta] | \ell, m \rangle.
\end{aligned}
```


Using
```math
L_y = (L₊ − L₋) / (2i)
```
we can expand
```math
exp[-iL_y β]
=
Σ_k (-iL_y β)^k / k!
=
Σ_k (L₋ - L₊)^k (β/2)^k / k!
```

Now, writing ``d_+(X) = [L_+, X]``, Eq. (9) of https://arxiv.org/pdf/1707.03861 says
```math
(L₋ - L₊)^k = \sum_{j=0}^k \binom{k, j} ((L₋ - d_+)^j 1) (-L₊)^{k-j}
```
The sum will automatically be zero unless ``m+k-j ≤ ℓ`` — which means ``j ≥ m+k-ℓ``
```math
(-L₊)^{k-j}|ℓ,m\rangle = (-1)^{k-j} \sqrt{\frac{(\ell+m+k-j)!}{(\ell+m)!},\frac{(\ell-m)!}{(\ell-m-k+j)!}} |ℓ,m+k-j\rangle
```

``[L₊, L₋] = 2 L_z``

``[L_z, L_\pm] = \pm L_\pm``

I wonder if there's a nicer approach using the symmetry transformation
Edmonds notes in Sec. 4.5 (and credits to Wigner) — or the presumably
equivalent one McEwan and Wieux use (and credit Risbo):
```math
\exp\left[ \beta 𝐣 / 2 \right]
=
\exp\left[ \pi 𝐤 / 4 \right]
\exp\left[ \pi 𝐣 / 4 \right]
\exp\left[ \beta 𝐤 / 2 \right]
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
    - Demonstrate that ``\mathfrak{D}`` is a representation
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
e^{im\phi}\}`` is an orthonormal basis of ``L^2(0,2\pi)``, while the
set ``\{1/c_{n,m} P_n^m(\cos\theta)`` is an orthonormal basis of
``L^2(0, \pi)`` in the ``\theta`` coordinate.  Therefore, the product
of these two sets is an orthonormal basis of the product space
``L^2\left((0,2\pi) \times (0, \pi)\right)``, which forms a coordinate
space for ``S^2``.  I would probably modify this to point out that
``(0,2\pi)`` is really ``S^1``, and then we could extend it to point
out that you can throw on another factor of ``S^1`` to cover ``S^3``,
which happens to give us the Wigner D-matrices.

## Recursion relations

[Gumerov and Duraiswami (2001)](@cite Gumerov_2001) derive their
recursion relations by differentiating solutions of the Helmholtz
equation ``\nabla^2 \psi + k^2 \psi = 0`` as ``\tfrac{1}{k} \nabla
\psi``.  More precisely, they differentiate both sides of the equation
relating one solution to its rotated form — which naturally involves
Wigner's ``\mathfrak{D}`` matrix.  Using orthogonal basis functions
for the solution, this allows them to equate terms on the two sides
proportional to a given basis function, which leaves them with
expressions involving sums of only the ``\mathfrak{D}`` matrices and
some coefficients depending on the indices of the basis functions (and
hence of ``\mathfrak{D}``) on both sides of the equation.  Since
``\nabla`` is a 3-vector operator, this gives them three relations.

This, of course, is happening in 3-D space, since ``\psi`` is a
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

The SWSHs/``\mathfrak{D}`` functions can be naturally promoted to
functions not just on the 3-sphere, but also in 4-D space just by
allowing the quaternions to be non-unit quaternions.
