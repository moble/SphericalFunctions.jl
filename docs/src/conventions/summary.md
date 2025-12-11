# Summary

This page lists the most important conventions used in this package.
The [following page](@ref "Details") derives all of these conventions
from the very basics (i.e., starting from Cartesian coordinates of
3-dimensional space).

Note that we will use Euler angles and spherical coordinates here, but
*they are not used internally in this package* — though conversion
functions are available.  It is almost always a bad idea to use Euler
angles in *computing*; quaternions are clearly the preferred
representation for numerous reasons.  However, Euler angles are
important for (a) comparing to other sources, and (b) performing
*analytic* integrations.  These are the only two uses we will make of
Euler angles.

## Fundamental coordinates
We use standard right-handed Cartesian coordinates ``(x, y, z)`` and
unit basis vectors ``(𝐱, 𝐲, 𝐳)``.

## Spherical coordinates
We define spherical coordinates ``(r, \theta, \phi)`` and unit basis
vectors ``(𝐧, \boldsymbol{\theta}, \boldsymbol{\phi})``.  The "polar
angle" ``\theta \in [0, \pi]`` measures the angle between the
specified direction and the positive ``𝐳`` axis.  The "azimuthal
angle" ``\phi \in [0, 2\pi)`` measures the angle between the
projection of the specified direction onto the ``𝐱``-``𝐲`` plane and
the positive ``𝐱`` axis, with the positive ``𝐲`` axis corresponding
to the positive angle ``\phi = \pi/2``.

## Quaternions
A quaternion is written ``𝐐 = W + X𝐢 + Y𝐣 + Z𝐤``, where ``𝐢^2 =
𝐣^2 = 𝐤^2 = 𝐢𝐣𝐤 = -1``.  In the code, this quaternion is
represented by ``(W, X, Y, Z)``.

We will frequently depict a three-dimensional vector ``𝐯 = v_x 𝐱 +
v_y 𝐲 + v_z 𝐳`` interchangeably as a quaternion ``v_x 𝐢 + v_y 𝐣 +
v_z 𝐤``.  Even though they really belong to different spaces, there
is a (vector-space) isomorphism between them, and the subspaces are
even isomorphic under the *algebra* isomorphism given by duality in
the geometric algebra.  This translation allows us to operate on
vectors as if they were quaternions, and vice versa.

## Quaternion rotations
A rotation represented by the unit quaternion ``𝐑`` acts on a vector
``𝐯`` as ``𝐑\, 𝐯\, 𝐑^{-1}``.  Where relevant, rotations will be
assumed to be right-handed, so that a quaternion characterizing the
rotation through an angle ``\vartheta`` about a unit vector ``𝐮`` can
be expressed as ``𝐑 = \exp(\vartheta 𝐮/2)``.  Note that ``-𝐑``
would deliver the same *rotation*, which means that the group of unit
quaternions ``\mathrm{Spin}(3) = \mathrm{SU}(2)`` is a *double cover*
of the group of rotations ``\mathrm{SO}(3)``.  Nonetheless, ``𝐑`` and
``-𝐑`` are distinct quaternions, and represent distinct "spinors".

## Spherical coordinates as quaternions
A point on the unit sphere with spherical coordinates ``(\theta,
\phi)`` can be represented by the unit quaternion
```math
𝐑_{\theta, \phi}
=
\exp(\phi 𝐤/2)\, \exp(\theta 𝐣/2).
```
This not only takes the positive ``𝐳`` axis to the specified
direction, but also takes the ``𝐱`` and ``𝐲`` axes onto the unit
basis vectors of the spherical coordinate system:
```math
\begin{aligned}
𝐧 &= 𝐑_{\theta, \phi}\, 𝐳\, 𝐑_{\theta, \phi}^{-1}, \\
\boldsymbol{\theta} &= 𝐑_{\theta, \phi}\, 𝐱\, 𝐑_{\theta, \phi}^{-1}, \\
\boldsymbol{\phi} &= 𝐑_{\theta, \phi}\, 𝐲\, 𝐑_{\theta, \phi}^{-1}.
\end{aligned}
```

## Euler angles (and spherical coordinates)
Euler angles parametrize a unit quaternion as
```math
𝐑_{\alpha, \beta, \gamma}
=
\exp(\alpha 𝐤/2)\, \exp(\beta 𝐣/2)\, \exp(\gamma 𝐤/2).
```
The angles ``\alpha`` and ``\gamma`` take values in ``[0, 2\pi)``.
The angle ``\beta`` takes values in ``[0, 2\pi]`` to parametrize the
group of unit quaternions ``\mathrm{Spin}(3) = \mathrm{SU}(2)``, or in
``[0, \pi]`` to parametrize the group of rotations ``\mathrm{SO}(3)``.

By comparison, we can immediately see that spherical coordinates
``(\theta, \phi)`` can be represented as Euler angles with the
equivalence ``(\alpha, \beta, \gamma) = (\phi, \theta, 0)``.  In
particular, any function of spherical coordinates can be promoted to a
function on Euler angles using this identification.

It's worth noting that the action of Euler angles on the Cartesian
basis is similar to the action of the spherical-coordinate quaternion,
but rotates the tangent basis ($\boldsymbol{\theta},
\boldsymbol{\phi}$).  That is, we still have
```math
𝐧 = 𝐑_{\phi, \theta, \gamma}\, 𝐳\, 𝐑_{\phi, \theta, \gamma}^{-1},
```
but the action on the ``𝐱`` and ``𝐲`` axes is a little more
complicated due to the initial rotation by ``\exp(\gamma 𝐤/2)``,
which is equivalent to a *final* rotation through ``\gamma`` about
``𝐧``.  It's easier to write this down if we form the combination
```math
𝐦 = \frac{1}{\sqrt{2}} \left(
    \boldsymbol{\theta} + i \boldsymbol{\phi}
\right),
```
and find that
```math
𝐦 = e^{-i\gamma} 𝐑_{\phi, \theta, \gamma}\, \frac{1}{\sqrt{2}} \left(
    𝐱 + i 𝐲
\right)\, 𝐑_{\phi, \theta, \gamma}^{-1}.
```

## Left and right angular-momentum operators
For a complex-valued function ``f(𝐑)``, we define two operators, the
left and right angular-momentum operators:
```math
L_𝐮 f(𝐑) = \left.i \frac{d}{d\epsilon}\right|_{\epsilon=0}
f\left(e^{-\epsilon 𝐮/2}\, 𝐑\right)
\qquad \text{and} \qquad
R_𝐮 f(𝐑) = -\left.i \frac{d}{d\epsilon}\right|_{\epsilon=0}
f\left(𝐑\, e^{-\epsilon 𝐮/2}\right),
```
where ``𝐮`` can be any quaternion, though unit pure-vector
quaternions are the most common.  In particular, ``L`` represents the
standard angular-momentum operators, and we can compute the
expressions in Euler angles for the basis vectors:
```math
\begin{aligned}
L_𝐢 &= i \left\{
    \frac{\cos\alpha}{\tan\beta} \frac{\partial} {\partial \alpha}
    + \sin\alpha \frac{\partial} {\partial \beta}
    - \frac{\cos\alpha}{\sin\beta} \frac{\partial} {\partial \gamma}
\right\},
&
R_𝐢 &= i \left\{
    -\frac{\cos\gamma}{\sin\beta} \frac{\partial} {\partial \alpha}
    +\sin\gamma \frac{\partial} {\partial \beta}
    +\frac{\cos\gamma}{\tan\beta} \frac{\partial} {\partial \gamma}
\right\},
\\
L_𝐣 &= i \left\{
    \frac{\sin\alpha}{\tan\beta} \frac{\partial} {\partial \alpha}
    - \cos\alpha \frac{\partial} {\partial \beta}
    -\frac{\sin\alpha}{\sin\beta} \frac{\partial} {\partial \gamma}
\right\},
&
R_𝐣 &= i \left\{
    \frac{\sin\gamma}{\sin\beta} \frac{\partial} {\partial \alpha}
    +\cos\gamma \frac{\partial} {\partial \beta}
    -\frac{\sin\gamma}{\tan\beta} \frac{\partial} {\partial \gamma}
\right\},
\\
L_𝐤 &= -i \frac{\partial} {\partial \alpha},
&
R_𝐤 &= i \frac{\partial} {\partial \gamma}.
\end{aligned}
```
These correspond precisely to the standard expressions for the
angular-momentum operators, with ``𝐢 \leftrightarrow 𝐱``, etc.  We
also obtain a generalization of the usual commutator relations and
find that
```math
[L_𝐮, L_𝐯] = \frac{i}{2} L_{[𝐮,𝐯]}
\qquad
[R_𝐮, R_𝐯] = \frac{i}{2} R_{[𝐮,𝐯]}
\qquad
[L_𝐮, R_𝐯] = 0.
```
Restricting to just the basis vectors, indexed as ``a,b,c``, the first
of these reduces to ``[L_a, L_b] = i \epsilon_{abc} L_c``, which is
precisely the standard result.  We can also lift any function on
``𝕊²`` to a function on ``𝕊³`` — or more precisely any function on
spherical coordinates to a function on the space of Euler angles — by
the correspondence ``(\theta, \phi) \mapsto (\alpha, \beta, \gamma) =
(\phi, \theta, 0)``.  We can then express the angular-momentum
operators in their more common form, in terms of spherical
coordinates:
```math
L_x = i \left\{
    \frac{\cos\phi}{\tan\theta} \frac{\partial} {\partial \phi}
    + \sin\phi \frac{\partial} {\partial \theta}
\right\}
\qquad
L_y = i \left\{
    \frac{\sin\phi}{\tan\theta} \frac{\partial} {\partial \phi}
    - \cos\phi \frac{\partial} {\partial \theta}
\right\}
\qquad
L_z = -i \frac{\partial} {\partial \phi}
```
The ``R`` operators make less sense for a function of spherical
coordinates, because of their inherent dependence on ``\gamma.``  We
will come back to them, however, when we consider spin-weighted
functions — which are inherently ill-defined on the 2-sphere, but can
be interpreted as restrictions of functions on the 3-sphere with this
special "weight" property.

## Spherical harmonics
There is essentially no disagreement in the literature about the
definitions of the spherical harmonics, so we adopt a function that is
consistent with the standard expressions.  More specifically, this
package defines the spherical harmonics in terms of Wigner's 𝔇
matrices, by way of the spin-weighted spherical harmonics, as a
function of a quaternion:
```math
Y_{l,m}(\mathbf{Q}) = \sqrt{\frac{2\ell+1}{4\pi}} e^{im\phi}
    D^{(l)}_{m,0}(\mathbf{Q}),
```
where ``D^{(l)}_{m,0}`` is the Wigner 𝔇 matrix.  This is a

For concreteness, however, we can write the standard expression in
terms of spherical coordinates.  This is what our definition will
reduce to for spin weight ``s=0``, and transforming the spherical
coordinates into a quaternion in the way given above.  Explicitly, in
terms of spherical coordinates, that expression is
```math
\begin{align}
  Y_{l,m}
  &=
  \sqrt{\frac{2\ell+1}{4\pi}} e^{im\phi}
  \sum_{k = k_1}^{k_2}
  \frac{(-1)^k \ell! [(\ell+m)!(\ell-m)!]^{1/2}}
  {(\ell+m-k)!(\ell-k)!k!(k-m)!}
  \\ &\qquad \times
  \left(\cos\left(\frac{\theta}{2}\right)\right)^{2\ell+m-2k}
  \left(\sin\left(\frac{\theta}{2}\right)\right)^{2k-m}
\end{align}
```
where ``k_1 = \textrm{max}(0, m)`` and ``k_2=\textrm{min}(\ell+m,
\ell)``.  Again, we must emphasize that this package does not actually
use this form; it is just shown here to make it easier to compare to
other sources.

## Spin-weighted functions
[Newman_1966](@citet) define the spherical tangent basis vectors
as
```math
m^\mu = \frac{1}{\sqrt{2}} \left(
    \boldsymbol{\theta} + i \boldsymbol{\phi}
\right)^\mu
```
and discuss spin weight in terms of the rotation
```math
(m^\mu)' = e^{i\psi} m^\mu,
```
where the tangent basis rotates but we are "keeping the
coordinates fixed".  They then define a function to have
spin weight ``s`` if it transforms as
```math
\eta' = e^{is\psi} \eta.
```
Such functions are generally the result of contracting a tensor field
with some number of ``m^\mu`` and some number of ``\bar{m}^\mu``
vectors (though spinor extensions resulting in half-integer spin
weights are also possible).

Note that this definition shows that it is clearly *impossible* to
define spin-weighted functions on the 2-sphere; the 2-sphere alone
includes no information about the directions of basis vectors in its
tangent space.  Instead, we *must* think of spin-weighted functions as
defined on the "unit tangent bundle over the 2-sphere" so that this
behavior with respect to rotation of tangent basis can possibly have
any effect.  This unit tangent bundle happens to be homeomorphic to
the 3-sphere, which is also the space of unit quaternions.  Thus, we
think of spin-weighted functions as defined on the group of unit
quaternions ``\mathrm{Spin}(3)=\mathrm{SU}(2)``, and frequently
discuss them in terms of Euler angles.

As we saw [above](@ref "Euler angles (and spherical coordinates)"),
``m^\mu`` corresponds to the Euler angles ``(\phi, \theta, 0)``, while
``(m^\mu)'`` corresponds to the Euler angles ``(\phi, \theta,
-\psi)``.  The function, written in terms of Euler angles, becomes
```math
\eta(\phi, \theta, -\psi) = e^{is\psi} \eta(\phi, \theta, 0),
```
or
```math
\eta(\alpha, \beta, \gamma) = e^{-is\gamma} \eta(\alpha, \beta, 0).
```
This is the crucial definition giving us the behavior of
spin-weighted functions: they are eigenfunctions of the operator
``R_z = i \partial_\gamma`` with eigenvalue ``s``.

We can make this a little less dependent on the choice of Euler
angles by writing ``\eta`` not as a function of Euler angles, but as
a function of a quaternion.  We then have
```math
\eta(\mathbf{Q}\, e^{\gamma 𝐤/2}) = e^{-is\gamma} \eta(\mathbf{Q}),
```
which means that spin-weighted functions are eigenfunctions of the
operator ``R_𝐤`` with eigenvalue ``s``.

We can also immediately find the spin-raising and -lowering operators
— canonically denoted ``\eth`` and ``\bar{\eth}`` — from the
commutator relations for ``R``:
```math
\begin{aligned}
\eth \eta &= \left(R_x + i R_y\right)\eta
    = -\sin^s \theta \left\{
        \frac{\partial}{\partial \theta}
        + \frac{i}{\sin\theta} \frac{\partial}{\partial \phi}
    \right\} \left(\eta \sin^{-s} \theta\right), \\
\bar{\eth} \eta &= \left(R_x - i R_y\right)\eta
    = -\sin^s \theta \left\{
        \frac{\partial}{\partial \theta}
        - \frac{i}{\sin\theta} \frac{\partial}{\partial \phi}
    \right\} \left(\eta \sin^{-s} \theta\right).
\end{aligned}
```
Here, we have used the full expressions for ``R_x`` and ``R_y``
given above in terms of Euler angles, replacing the derivatives
with respect to ``\gamma`` by a factor of ``-i s``, and converting
the remaining Euler angles to spherical coordinates.  This allows
us to write them as if they were operators on the 2-sphere, even
though this is mathematically ill-defined and spin-weighted
functions really must be defined on the 3-sphere.

## Spin-weighted spherical harmonics

Given the (scalar) spherical harmonics, and the spin-raising and
-lowering operators, we can now define the spin-weighted spherical
harmonics.  These are obtained by applying the relevant operator to
the scalar spherical harmonics the specified number of times, and
normalizing.  Again, this results in a function of a quaternion, but
we can write it in terms of spherical coordinates purely for the sake
of comparison with other sources.  The expression is
```math
\begin{align}
  {}_{s}Y_{l,m}
  &=
  (-1)^s\sqrt{\frac{2\ell+1}{4\pi}} e^{im\phi}
  \sum_{k = k_1}^{k_2}
  \frac{(-1)^k[(\ell+m)!(\ell-m)!(\ell-s)!(\ell+s)!]^{1/2}}
  {(\ell+m-k)!(\ell+s-k)!k!(k-s-m)!}
  \\ &\qquad \times
  \left(\cos\left(\frac{\iota}{2}\right)\right)^{2\ell+m+s-2k}
  \left(\sin\left(\frac{\iota}{2}\right)\right)^{2k-s-m}
\end{align}
```
where ``k_1 = \textrm{max}(0, m+s)`` and ``k_2=\textrm{min}(\ell+m,
\ell+s)``.


## Wigner D-matrices


