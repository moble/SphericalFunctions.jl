# Conventions

Here, we work through all the conventions used in this package,
starting from first principles to motivate the choices and ensure that
each step is on firm footing.

## Three-dimensional space

The space we are working in is naturally three-dimensional Euclidean
space, so we start with Cartesian coordinates ``(x, y, z)``.  These
also give us the unit basis vectors ``(𝐱, 𝐲, 𝐳)``.  Note that these
basis vectors are assumed to have unit norm, but we omit the hats just
to keep the notation simple.  Any vector in this space can be written
as
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
If we restrict to just the unit sphere, we can simplify this to
```math
\int f\, d^2\Omega = \int_0^\pi \int_0^{2\pi} f\, \sin\theta\, d\theta\, d\phi.
```


## Four-dimensional space: Quaternions and rotations


## Rotations


## Euler angles and spherical coordinates


