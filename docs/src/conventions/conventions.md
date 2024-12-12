# Conventions

Here, we work through all the conventions used in this package,
starting from first principles to motivate the choices and ensure that
each step is on firm footing.

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
notation simple.

One seemingly obvious — but extremely important — fact is that the
unit basis frame ``(𝐱, 𝐲, 𝐳)`` can be rotated onto
``(\boldsymbol{\theta}, \boldsymbol{\phi}, \mathbf{r})`` by first
rotating through the "polar" angle ``\theta`` about the ``\mathbf{y}``
axis, and then through the "azimuthal" angle ``\phi`` about the
``\mathbf{z}`` axis.  This becomes important when we consider
spin-weighted functions.

Integration in Cartesian coordinates is, of course, trivial as
```math
\int f\, d^3\mathbf{r} = \int_{-\infty}^{\infty} \int_{-\infty}^{\infty} \int_{-\infty}^{\infty} f\, dx\, dy\, dz.
```
In spherical coordinates, the integrand involves the square-root of
the determinant of the metric, so we have
```math
\int f\, d^3\mathbf{r} = \int_0^\infty \int_0^\pi \int_0^{2\pi} f\, r^2 \sin\theta\, dr\, d\theta\, d\phi.
```
Restricting to the unit sphere, and normalizing so that the integral
of 1 over the sphere is 1, we can simplify this to
```math
\int f\, d^2\Omega = \frac{1}{4\pi} \int_0^\pi \int_0^{2\pi} f\, \sin\theta\, d\theta\, d\phi.
```


## Four-dimensional space: Quaternions and rotations

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
geometric product is associative, distributive, and has the property
that
```math
𝐯𝐯 = \| 𝐯 \|^2.
```
The basis for this entire space is then the set
```math
\begin{gather}
𝟏, \\
𝐱, 𝐲, 𝐳,\\
𝐱𝐲, 𝐱𝐳, 𝐲𝐳, \\
𝐱𝐲𝐳.
\end{gather}
```
It's useful to note that the first four of these square to 1, while
the last four square to -1 — meaning that they could serve as a unit
imaginary to generate the complex numbers.  The more standard symbols
— and the ones we will use — are
```math
\begin{gather}
𝟏, \\
𝐱, 𝐲, 𝐳,\\
𝐢, 𝐣, 𝐤, \\
𝐈.
\end{gather}
```
The interpretation of these is that ``𝟏`` represents the scalars;
``𝐱, 𝐲, 𝐳`` span the vectors; ``𝐢, 𝐣, 𝐤`` are the standard
quaternion components; and ``𝐈`` is the pseudoscalar, which can also
serve as the [Hodge
dual](https://en.wikipedia.org/wiki/Hodge_star_operator).  (Note that
quaternions will only be spanned by elements made from an even number
of the basis vectors.  It turns out that those with an odd number will
produce reflections, rather than rotations, when acting on a vector —
as discussed below.  This explains why quaternions are restricted to
just those elements with an even number to represent rotations.  For
details see any geometric algebra text, like [Doran and Lasenby](@cite
DoranLasenby_2010).)

The key expressions that help to determine the arbitrary choices we
have made thus far are the multiplications
```math
\begin{aligned}
𝐢 𝐣 &= 𝐤, \\
𝐣 𝐤 &= 𝐢, \\
𝐤 𝐢 &= 𝐣.
\end{aligned}
```
Everyone agrees that ``𝐢² = 𝐣² = 𝐤² = -1``, so we can also use the
rules above to determine ``𝐢𝐣𝐤 = -𝟏``.  Different conventions are
sometimes used (almost exclusively in aerospace) so that this last
equation and the three displayed above have a flipped sign.  See
[Sommer et al.](@cite SommerEtAl_2018) for a discussion of the
different conventions.

We use coordinates ``(W, X, Y, Z)`` on the space of quaternions, so
that such a quaternion would be written as
```math
W𝟏 + X𝐢 + Y𝐣 + Z𝐤,
```
though we usually omit the ``𝟏``.  As with standard three-dimensional
space, we could introduce spherical coordinates, though we use a
slight variant: extended Euler coordinates.  In our conventions, we
have
```math
\begin{aligned}
R &= \sqrt{W^2 + X^2 + Y^2 + Z^2} &&\in [0, \infty), \\
\alpha &= \arctan\frac{Z}{W} + \arctan\frac{-X}{Y} &&\in [0, 2\pi], \\
\beta &= 2\arccos\sqrt{\frac{W^2+Z^2}{W^2+X^2+Y^2+Z^2}} &&\in [0, 2\pi), \\
\gamma &= \arctan\frac{Z}{W} - \arctan\frac{-X}{Y} &&\in [0, 2\pi],
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
cover the space of rotations.  This and the inclusion of ``R``
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
Again, integration involves a square-root of the determinant of the
metric, which reduces to ``R^3 |\sin\beta| / 8``.
```math
\int f\, d^4𝐑
= \int_{-\infty}^\infty \int_{-\infty}^\infty \int_{-\infty}^\infty \int_{-\infty}^\infty f\, dW\, dX\, dY\, dZ
= \int_0^\infty \int_0^{2\pi} \int_0^{2\pi} \int_0^{2\pi} f\, \frac{R^3}{8} |\sin β|\, dR\, dα\, dβ\, dγ.
```
Restricting to the unit sphere, and normalizing so that the integral
of 1 over the sphere is 1, we can simplify this to
```math
\int f\, d^3\Omega = \frac{1}{16\pi^2} \int_0^{2\pi} \int_0^{2\pi} \int_0^{2\pi} f\, |\sin β|\, dα\, dβ\, dγ.
```

## Rotations


## Euler angles and spherical coordinates


