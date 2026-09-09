# Summary

This page lists the most important conventions used in this package.
The [following page](@ref "Details") derives all of these conventions
from the very basics (i.e., starting from Cartesian coordinates of
3-dimensional space), and each section below links to the
corresponding derivation.  The [comparisons](@ref "Comparisons") pages
then relate these conventions to the literature and to other software.

Note that we will use Euler angles and spherical coordinates here, but
*they are not used internally in this package* — though conversion
functions are available.  It is almost always a bad idea to use Euler
angles in *computing*; quaternions are clearly the preferred
representation for numerous reasons.  However, Euler angles are
important for (a) comparing to other sources, and (b) performing
*analytic* integrations.  These are the only two uses we will make of
Euler angles.

Throughout, ``ℓ``, ``m'``, ``m``, and ``s`` are either all integers or
all half-integers; every statement below holds in both cases unless
noted otherwise.

## [Fundamental coordinates](@id summary_cartesian)
We use standard right-handed Cartesian coordinates ``(x, y, z)`` and
unit basis vectors ``(𝐱, 𝐲, 𝐳)``.
[Details.](@ref conv_three_dimensional_space)

## [Spherical coordinates](@id summary_spherical_coordinates)
We define spherical coordinates ``(r, θ, ϕ)`` and unit basis
vectors ``(𝐧, \boldsymbol{θ}, \boldsymbol{ϕ})`` in the standard
"physics" convention:
```math
x = r \sin θ \cos ϕ,
\qquad
y = r \sin θ \sin ϕ,
\qquad
z = r \cos θ.
```
The "polar angle" ``θ \in [0, π]`` measures the angle between the
specified direction and the positive ``𝐳`` axis.  The "azimuthal
angle" ``ϕ \in [0, 2π)`` measures the angle between the projection of
the specified direction onto the ``𝐱``-``𝐲`` plane and the positive
``𝐱`` axis, with the positive ``𝐲`` axis corresponding to the
positive angle ``ϕ = π/2``.  The surface element on the unit sphere is
``\sin θ\, dθ\, dϕ``, with total area ``4π``.  [Details.](@ref
conv_spherical_coordinates)

## [Quaternions](@id summary_quaternions)
A quaternion is written ``𝐐 = W + X𝐢 + Y𝐣 + Z𝐤``, where
```math
𝐢^2 = 𝐣^2 = 𝐤^2 = -1,
\qquad
𝐢𝐣 = 𝐤, \quad 𝐣𝐤 = 𝐢, \quad 𝐤𝐢 = 𝐣,
\qquad
𝐢𝐣𝐤 = -1.
```
In the code, this quaternion is represented by the components ``(W,
X, Y, Z)``, in that order.  The conjugate is ``\overline{𝐐} = W - X𝐢
- Y𝐣 - Z𝐤``, the norm is ``\|𝐐\| = \sqrt{W^2+X^2+Y^2+Z^2}``, and
for a unit "pure-vector" quaternion ``𝐮`` we have ``\exp(𝐮\,θ) =
\cos θ + 𝐮 \sin θ``.

We will frequently depict a three-dimensional vector ``𝐯 = v_x 𝐱 +
v_y 𝐲 + v_z 𝐳`` interchangeably as a quaternion ``v_x 𝐢 + v_y 𝐣 +
v_z 𝐤``.  Even though they really belong to different spaces, there
is a (vector-space) isomorphism between them — given by duality in the
geometric algebra, ``𝐢 = 𝐈^{-1}𝐱``, etc. — which allows us to
operate on vectors as if they were quaternions, and vice versa.
[Details.](@ref conv_quaternions)

## [Quaternion rotations](@id summary_rotations)
A rotation represented by the unit quaternion ``𝐑`` acts on a vector
``𝐯`` as ``𝐑\, 𝐯\, 𝐑^{-1}``.  Rotations are right-handed, so that
the quaternion characterizing the rotation through an angle ``ϑ``
about a unit vector ``𝐮`` is ``𝐑 = \exp(ϑ 𝐮/2)``; in particular
``𝐢``, ``𝐣``, ``𝐤`` generate positive rotations about ``𝐱``,
``𝐲``, ``𝐳``.  Note that ``-𝐑`` delivers the same *rotation*, which
means that the group of unit quaternions ``\mathrm{Spin}(3) =
\mathrm{SU}(2)`` is a *double cover* of the group of rotations
``\mathrm{SO}(3)``.  Nonetheless, ``𝐑`` and ``-𝐑`` are distinct
quaternions, and represent distinct "spinors".  [Details.](@ref
conv_rotations)

## [Spherical coordinates as quaternions](@id summary_spherical_quaternions)
A point on the unit sphere with spherical coordinates ``(θ,
ϕ)`` can be represented by the unit quaternion
```math
𝐑_{θ, ϕ}
=
\exp(ϕ 𝐤/2)\, \exp(θ 𝐣/2).
```
This not only takes the positive ``𝐳`` axis to the specified
direction, but also takes the ``𝐱`` and ``𝐲`` axes onto the unit
basis vectors of the spherical coordinate system:
```math
\begin{aligned}
𝐧 &= 𝐑_{θ, ϕ}\, 𝐳\, 𝐑_{θ, ϕ}^{-1}, \\
\boldsymbol{θ} &= 𝐑_{θ, ϕ}\, 𝐱\, 𝐑_{θ, ϕ}^{-1}, \\
\boldsymbol{ϕ} &= 𝐑_{θ, ϕ}\, 𝐲\, 𝐑_{θ, ϕ}^{-1}.
\end{aligned}
```
[Details.](@ref conv_euler_angles)

## [Euler angles](@id summary_euler_angles)
Euler angles parametrize a unit quaternion as
```math
𝐑_{α, β, γ}
=
\exp(α 𝐤/2)\, \exp(β 𝐣/2)\, \exp(γ 𝐤/2),
```
i.e., a rotation through ``γ`` about ``𝐳``, followed by ``β`` about
the *fixed* ``𝐲`` axis, followed by ``α`` about the *fixed* ``𝐳``
axis.  The ranges are
```math
α \in [0, 2π),
\qquad
β \in [0, π],
\qquad
γ \in [0, 4π)
```
to cover the group of unit quaternions ``\mathrm{Spin}(3) =
\mathrm{SU}(2)`` exactly once (up to subsets of measure zero);
restricting to ``γ \in [0, 2π)`` covers the rotation group
``\mathrm{SO}(3)`` instead.  With the components ``(W, X, Y, Z)`` of
``𝐑_{α,β,γ}`` we have
```math
\begin{aligned}
  W &= \cos\tfrac{β}{2} \cos\tfrac{α+γ}{2}, &
  X &= -\sin\tfrac{β}{2} \sin\tfrac{α-γ}{2}, \\
  Y &= \sin\tfrac{β}{2} \cos\tfrac{α-γ}{2}, &
  Z &= \cos\tfrac{β}{2} \sin\tfrac{α+γ}{2}.
\end{aligned}
```
The invariant (Haar) measure inherited from ``\mathbb{R}^4`` is
```math
\int_{\mathrm{Spin}(3)} f\, d^3\Omega
= \frac{1}{8} \int_0^{2π} \int_0^{π} \int_0^{4π} f\, \sin β\, dα\, dβ\, dγ,
\qquad
\int_{\mathrm{Spin}(3)} d^3\Omega = 2π^2,
```
and the corresponding integral over ``\mathrm{SO}(3)`` (with ``γ \in
[0, 2π)``) has total volume ``π^2``.  [Details](@ref
Quaternions-and-Euler-angles) and [measure](@ref conv_haar_measure).

By comparison, we can immediately see that spherical coordinates ``(θ,
ϕ)`` can be represented as Euler angles with the equivalence ``(α, β,
γ) = (ϕ, θ, 0)``.  In particular, any function of spherical
coordinates can be promoted to a function on Euler angles — or on
``\mathrm{Spin}(3)`` — using this identification.  More generally,
``𝐑_{ϕ, θ, γ}`` still takes ``𝐳`` to ``𝐧``, but rotates the tangent
basis: writing ``𝐦 = (\boldsymbol{θ} + i \boldsymbol{ϕ})/\sqrt{2}``,
```math
𝐦 = e^{-iγ}\, 𝐑_{ϕ, θ, γ}\, \frac{𝐱 + i 𝐲}{\sqrt{2}}\, 𝐑_{ϕ, θ, γ}^{-1}.
```
[Details.](@ref conv_euler_angles)

## [Left and right angular-momentum operators](@id summary_L_R_definitions)
For a complex-valued function ``f(𝐑)`` on ``\mathrm{Spin}(3)``, we
define two operators, the left and right angular-momentum operators:
```math
L_𝐮 f(𝐑) = \left.i \frac{d}{dϵ}\right|_{ϵ=0}
f\left(e^{-ϵ 𝐮/2}\, 𝐑\right)
\qquad \text{and} \qquad
R_𝐮 f(𝐑) = -\left.i \frac{d}{dϵ}\right|_{ϵ=0}
f\left(𝐑\, e^{-ϵ 𝐮/2}\right),
```
where ``𝐮`` can be any quaternion, though unit pure-vector
quaternions are the most common; we write ``L_x = L_𝐢``, etc.  The
signs (``i`` on the left, ``-i`` on the right) are fixed by requiring
that the raising and lowering operators take the conventional forms
```math
L_\pm = L_x \pm i L_y,
\qquad
R_\pm = R_x \pm i R_y,
\qquad
[L_z, L_\pm] = \pm L_\pm,
\qquad
[R_z, R_\pm] = \pm R_\pm.
```
Both sets of operators then obey the *standard* commutation relations,
```math
[L_𝐮, L_𝐯] = \frac{i}{2} L_{[𝐮,𝐯]},
\qquad
[R_𝐮, R_𝐯] = \frac{i}{2} R_{[𝐮,𝐯]},
\qquad
[L_𝐮, R_𝐯] = 0,
```
which for the basis vectors reduce to ``[L_a, L_b] = i ϵ_{abc} L_c``
and ``[R_a, R_b] = i ϵ_{abc} R_c``.  The Casimir operators coincide:
``L^2 = R^2``.  ``L`` is the standard angular-momentum operator of
quantum mechanics; ``R`` is (minus) the same operator expressed in the
body-fixed frame, ``R_𝐮 = -L_{𝐑 𝐮 𝐑^{-1}}``.  [Details.](@ref
conv_L_R_definitions)

### [Euler-angle and spherical-coordinate expressions](@id summary_L_R_euler)
In Euler angles the operators are
```math
\begin{aligned}
L_x &= i \left\{
    \frac{\cos α}{\tan β} \frac{\partial} {\partial α}
    + \sin α \frac{\partial} {\partial β}
    - \frac{\cos α}{\sin β} \frac{\partial} {\partial γ}
\right\},
&
R_x &= i \left\{
    -\frac{\cos γ}{\sin β} \frac{\partial} {\partial α}
    +\sin γ \frac{\partial} {\partial β}
    +\frac{\cos γ}{\tan β} \frac{\partial} {\partial γ}
\right\},
\\
L_y &= i \left\{
    \frac{\sin α}{\tan β} \frac{\partial} {\partial α}
    - \cos α \frac{\partial} {\partial β}
    -\frac{\sin α}{\sin β} \frac{\partial} {\partial γ}
\right\},
&
R_y &= i \left\{
    \frac{\sin γ}{\sin β} \frac{\partial} {\partial α}
    +\cos γ \frac{\partial} {\partial β}
    -\frac{\sin γ}{\tan β} \frac{\partial} {\partial γ}
\right\},
\\
L_z &= -i \frac{\partial} {\partial α},
&
R_z &= i \frac{\partial} {\partial γ}.
\end{aligned}
```
Lifting a function on ``𝕊²`` to ``\mathrm{Spin}(3)`` via ``(α, β, γ)
= (ϕ, θ, 0)``, the ``L`` operators reduce to their familiar
spherical-coordinate forms:
```math
L_x = i \left\{
    \frac{\cos ϕ}{\tan θ} \frac{\partial} {\partial ϕ}
    + \sin ϕ \frac{\partial} {\partial θ}
\right\},
\qquad
L_y = i \left\{
    \frac{\sin ϕ}{\tan θ} \frac{\partial} {\partial ϕ}
    - \cos ϕ \frac{\partial} {\partial θ}
\right\},
\qquad
L_z = -i \frac{\partial} {\partial ϕ},
```
```math
L_\pm = e^{\pm iϕ} \left\{
    \pm \frac{\partial}{\partial θ}
    + \frac{i}{\tan θ} \frac{\partial}{\partial ϕ}
\right\}.
```
The ``R`` operators have no counterpart on ``𝕊²``, because they
retain derivatives with respect to ``γ``; they act naturally only on
spin-weighted functions, discussed next.  All of these expressions are
derived symbolically on the [``L_j`` and ``R_j`` with Euler
angles](@ref euler_angular_momentum) page.  [Details.](@ref
conv_L_R_euler)

## [Spin-weighted functions](@id summary_spin_weight)
Following [Newman_1966](@citet), a function ``\eta`` has spin weight
``s`` if it picks up a phase ``e^{isψ}`` when the tangent basis ``𝐦 =
(\boldsymbol{θ} + i\boldsymbol{ϕ})/\sqrt{2}`` is rotated by ``e^{iψ}``
at fixed coordinates.  Since that rotation corresponds to
right-multiplication of the quaternion argument by ``e^{-ψ𝐤/2}``,
this means
```math
\eta\left(𝐐\, e^{γ 𝐤/2}\right) = e^{-isγ}\, \eta(𝐐)
\qquad \Longleftrightarrow \qquad
R_z\, \eta = s\, \eta.
```
Spin-weighted functions are therefore *eigenfunctions of ``R_z`` with
eigenvalue ``s``*; they cannot be defined on ``𝕊²`` itself, only on
``\mathrm{Spin}(3)`` (or, equivalently, on spherical coordinates
together with the reference tangent direction ``\boldsymbol{θ}``).
Newman and Penrose's spin-raising and -lowering operators are the
ladder operators of ``R_z``, up to a sign in the second case:
```math
\eth = R_+ = R_x + i R_y,
\qquad
\bar{\eth} = -R_- = -\left(R_x - i R_y\right).
```
In their conventional (but coordinate-dependent) spherical form,
```math
\eth \eta = -\sin^s θ \left\{
        \frac{\partial}{\partial θ}
        + \frac{i}{\sin θ} \frac{\partial}{\partial ϕ}
    \right\} \left(\eta \sin^{-s} θ\right),
\qquad
\bar{\eth} \eta = -\sin^{-s} θ \left\{
        \frac{\partial}{\partial θ}
        - \frac{i}{\sin θ} \frac{\partial}{\partial ϕ}
    \right\} \left(\eta \sin^{s} θ\right).
```
[Details.](@ref conv_spin_weight)

## [Wigner 𝔇 matrices](@id summary_wigner_D)
Rotating a field ``f`` by ``𝐑`` gives ``\left[U(𝐑) f\right](𝐐) =
f(𝐑^{-1}𝐐)``, with ``U(e^{ϑ𝐮/2}) = \exp(-iϑ L_𝐮)`` and hence
``U(𝐑_{α,β,γ}) = e^{-iαL_z} e^{-iβL_y} e^{-iγL_z}``.  Wigner's ``𝔇``
matrix is the matrix of ``U(𝐑)`` in an orthonormal basis ``|ℓ,
m\rangle`` of eigenfunctions of ``L^2`` and ``L_z`` (with any fixed
spin weight):
```math
𝔇^{(ℓ)}_{m',m}(𝐑) = \langle ℓ, m' | U(𝐑) | ℓ, m \rangle,
\qquad
U(𝐑)\, |ℓ, m\rangle = \sum_{m'} |ℓ, m'\rangle\, 𝔇^{(ℓ)}_{m',m}(𝐑).
```
It is a function of the rotation ``𝐑`` itself (not its inverse), the
first index ``m'`` pairs with the first Euler angle and the left
operators, and the second index ``m`` with the last Euler angle and
the right operators:
```math
𝔇^{(ℓ)}_{m',m}(α, β, γ)
= e^{-i m' α}\, d^{(ℓ)}_{m',m}(β)\, e^{-i m γ},
\qquad
d^{(ℓ)}_{m',m}(β) = \langle ℓ, m' | e^{-iβL_y} | ℓ, m \rangle,
```
where the real matrix ``d`` is given by Wigner's formula
```math
d^{(ℓ)}_{m',m}(β)
=
\sum_{k=\max(0, m-m')}^{\min(ℓ+m, ℓ-m')}
(-1)^{k - m + m'}
\frac{\sqrt{(ℓ+m)!\,(ℓ-m)!\,(ℓ+m')!\,(ℓ-m')!}}
     {(ℓ+m-k)!\,k!\,(ℓ-m'-k)!\,(k-m+m')!}
\left(\cos\frac{β}{2}\right)^{2ℓ+m-m'-2k}
\left(\sin\frac{β}{2}\right)^{2k-m+m'}.
```
This agrees with LALSuite, Wikipedia, Sakurai, Shankar, Zettili, and
Varshalovich et al.; it is the complex conjugate of the matrices used
by Wigner (who also includes ``(-1)^{m'-m}``), Edmonds, Goldberg et
al., Boyle (2016), and versions of this package before 3.0.

Basic properties:
```math
\begin{gathered}
𝔇^{(ℓ)}(𝐑_1 𝐑_2) = 𝔇^{(ℓ)}(𝐑_1)\, 𝔇^{(ℓ)}(𝐑_2),
\qquad
𝔇^{(ℓ)}(𝐑^{-1}) = 𝔇^{(ℓ)}(𝐑)^\dagger,
\qquad
𝔇^{(ℓ)}(-𝐑) = (-1)^{2ℓ}\, 𝔇^{(ℓ)}(𝐑),
\\
\overline{𝔇^{(ℓ)}_{m',m}} = (-1)^{m'-m}\, 𝔇^{(ℓ)}_{-m',-m},
\qquad
d^{(ℓ)}_{m',m}(β) = (-1)^{m'-m}\, d^{(ℓ)}_{m,m'}(β) = d^{(ℓ)}_{-m,-m'}(β) = d^{(ℓ)}_{m,m'}(-β),
\\
L^2\, 𝔇^{(ℓ)}_{m',m} = R^2\, 𝔇^{(ℓ)}_{m',m} = ℓ(ℓ+1)\, 𝔇^{(ℓ)}_{m',m},
\qquad
L_z\, 𝔇^{(ℓ)}_{m',m} = -m'\, 𝔇^{(ℓ)}_{m',m},
\qquad
R_z\, 𝔇^{(ℓ)}_{m',m} = m\, 𝔇^{(ℓ)}_{m',m},
\\
L_\pm\, 𝔇^{(ℓ)}_{m',m} = -\sqrt{(ℓ \pm m')(ℓ \mp m' + 1)}\, 𝔇^{(ℓ)}_{m' \mp 1, m},
\qquad
R_\pm\, 𝔇^{(ℓ)}_{m',m} = \sqrt{(ℓ \mp m)(ℓ \pm m + 1)}\, 𝔇^{(ℓ)}_{m', m \pm 1},
\\
\int_{\mathrm{Spin}(3)}
  \overline{𝔇^{(ℓ')}_{m'_1, m_1}}\; 𝔇^{(ℓ)}_{m', m}\; d^3Ω
= \frac{2π^2}{2ℓ+1}\, δ_{ℓ', ℓ}\, δ_{m'_1, m'}\, δ_{m_1, m}.
\end{gathered}
```
Note the signs in the ``L`` relations: as a function on
``\mathrm{Spin}(3)``, ``𝔇_{m',m}`` behaves like a row of expansion
coefficients; its complex conjugate is what behaves like a
wavefunction, with ``L_z`` eigenvalue ``+m'``.  [Details.](@ref
conv_wigner_D)

## [Spherical harmonics](@id summary_spherical_harmonics)
There is essentially no disagreement in the literature about the
spherical harmonics: everyone uses the Condon–Shortley phase
convention, in which the coefficients of ``L_\pm`` are real and
positive and ``Y_{ℓ,0}`` is positive on the ``+𝐳`` axis.  So do we.
Explicitly, with ``k_1 = \max(0, m)`` and ``k_2 = \min(ℓ+m, ℓ)``,
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
so that, e.g., ``Y_{1,\pm1} = \mp\sqrt{3/8π}\, \sin θ\, e^{\pm iϕ}``.
Lifted to ``\mathrm{Spin}(3)`` (as a function of spin weight 0), the
spherical harmonic is
```math
Y_{ℓ,m}(𝐑) = \sqrt{\frac{2ℓ+1}{4π}}\; \overline{𝔇^{(ℓ)}_{m,0}(𝐑)},
```
and it transforms under rotation of the field by ``𝐑`` as
```math
Y_{ℓ,m}\left(𝐑^{-1}\, 𝐐\right) = \sum_{m'} 𝔇^{(ℓ)}_{m',m}(𝐑)\, Y_{ℓ,m'}(𝐐),
\qquad \text{equivalently} \qquad
Y_{ℓ,m}\left(𝐑\, 𝐐\right) = \sum_{m'} \overline{𝔇^{(ℓ)}_{m,m'}(𝐑)}\, Y_{ℓ,m'}(𝐐).
```
[Details](@ref conv_spherical_harmonics) and
[rotation law](@ref conv_rotation_law).

## [Spin-weighted spherical harmonics](@id summary_swsh)
The spin-weighted spherical harmonics are the simultaneous
eigenfunctions of ``L^2``, ``L_z``, ``R_z`` with eigenvalues
``ℓ(ℓ+1)``, ``m``, ``s``, normalized as below, with phases fixed by
requiring the ladder coefficients of *both* ``L_\pm`` and ``R_\pm`` to
be real and positive and anchoring at the ordinary spherical harmonics
for ``s=0``.  This gives the definition
```math
{}_sY_{ℓ,m}(𝐑)
=
(-1)^s \sqrt{\frac{2ℓ+1}{4π}}\; \overline{𝔇^{(ℓ)}_{m,-s}(𝐑)},
\qquad
(-1)^s \equiv e^{iπs},
```
where the branch of ``(-1)^s`` matters only for half-integer ``s``
(for which it is ``\pm i``); equivalently, the phase is anchored by
``{}_sY_{ℓ,-s}(𝟏) = i^{2s}\sqrt{(2ℓ+1)/4π}``.  For integer indices
this can also be written ``{}_sY_{ℓ,m} = (-1)^m \sqrt{(2ℓ+1)/4π}\;
𝔇^{(ℓ)}_{-m,s}``; for half-integer indices that form is off by
``(-1)^{2s} = -1``, so the conjugate form is the definition.  The
``(-1)^s`` is not arbitrary: it is the Condon–Shortley condition
applied to ``R_+``, and it is exactly the factor appearing in the
LALSuite and NINJA expressions.  In spherical coordinates (``γ=0``),
with ``k_1 = \max(0, m+s)`` and ``k_2 = \min(ℓ+m, ℓ+s)``,
```math
\begin{aligned}
  {}_{s}Y_{ℓ,m}(θ, ϕ)
  &=
  (-1)^s \sqrt{\frac{2ℓ+1}{4π}}\; d^{(ℓ)}_{m,-s}(θ)\, e^{imϕ}
  \\
  &=
  (-1)^s\sqrt{\frac{2ℓ+1}{4π}}\, e^{imϕ}
  \sum_{k = k_1}^{k_2}
  \frac{(-1)^k[(ℓ+m)!(ℓ-m)!(ℓ-s)!(ℓ+s)!]^{1/2}}
  {(ℓ+m-k)!\,(ℓ+s-k)!\,k!\,(k-s-m)!}
  \left(\cos\frac{θ}{2}\right)^{2ℓ+m+s-2k}
  \left(\sin\frac{θ}{2}\right)^{2k-s-m}.
\end{aligned}
```
Again, we must emphasize that this package does not actually use this
form; it is shown here to make it easier to compare to other sources.
The essential properties are
```math
\begin{gathered}
L_\pm\, {}_sY_{ℓ,m} = \sqrt{(ℓ \mp m)(ℓ \pm m + 1)}\, {}_sY_{ℓ,m \pm 1},
\qquad
R_\pm\, {}_sY_{ℓ,m} = \sqrt{(ℓ \mp s)(ℓ \pm s + 1)}\, {}_{s \pm 1}Y_{ℓ,m},
\\
{}_sY_{ℓ,m}\left(𝐑\, e^{γ𝐤/2}\right) = e^{-isγ}\, {}_sY_{ℓ,m}(𝐑),
\qquad
{}_sY_{ℓ,m}\left(𝐑^{-1}\, 𝐐\right) = \sum_{m'} 𝔇^{(ℓ)}_{m',m}(𝐑)\, {}_sY_{ℓ,m'}(𝐐),
\\
\overline{{}_sY_{ℓ,m}} = (-1)^{m+s}\, {}_{-s}Y_{ℓ,-m},
\qquad
{}_sY_{ℓ,m}(-𝐑) = (-1)^{2ℓ}\, {}_sY_{ℓ,m}(𝐑),
\\
\int_{\mathrm{Spin}(3)} \overline{{}_{s'}Y_{ℓ',m'}}\; {}_sY_{ℓ,m}\; d^3Ω
= \frac{π}{2}\, δ_{ℓ',ℓ}\, δ_{m',m}\, δ_{s',s},
\qquad
\int_{𝕊^2} \overline{{}_{s}Y_{ℓ',m'}}\; {}_sY_{ℓ,m}\; \sin θ\, dθ\, dϕ
= δ_{ℓ',ℓ}\, δ_{m',m}.
\end{gathered}
```
[Details](@ref conv_swsh) and [half-integer indices](@ref
conv_swsh_half_integer).
