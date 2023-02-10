# Differential operators

Spin-weighted spherical functions *cannot* be defined on the sphere $S^2$, but
are well defined on the group ``\mathrm{Spin}(3) \cong \mathrm{SU}(2)`` or the
rotation group ``\mathrm{SO}(3)``.  (See [this
paper](https://arxiv.org/abs/1604.08140) for the explanation.)  However, this
also allows us to define a variety of differential operators acting on these
functions, relating to infinitesimal motions in these groups, acting either from
the left or the right on their arguments.  Right or left matters because the
groups mentioned above are all non-commutative groups.

In general, the *left* Lie derivative of a function ``f(Q)`` over the unit
quaternions with respect to a generator of rotation ``g`` is defined as
```math
L_g(f)\{Q\} := \frac{1}{2i} \left. \frac{f\left(\exp(t\,g)\, Q\right)}{dt} \right|_{t=0}.
```
Note that the exponential multiplies ``Q`` *on the left* — hence the name.

So, for example, a rotation about the ``z`` axis has the quaternion ``z`` as its
generator of rotation, and ``L_z`` defined in this way agrees with [the usual
angular-momentum
operator](https://en.wikipedia.org/wiki/Angular_momentum_operator) ``L_z``
familiar from spherical-harmonic theory, and reduces to it when the function has
spin weight 0, but also applies to functions of general spin weight.  Similarly,
we can compute ``L_x`` and ``L_y``, and take appropriate combinations to find
[the usual raising and lowering (ladder)
operators](https://en.wikipedia.org/wiki/Ladder_operator#Angular_momentum)
``L_+`` and ``L_-``.

In just the same way, we can define the *right* Lie derivative of a function
``f(Q)`` over the unit quaternions with respect to a generator of rotation ``g``
as
```math
R_g(f)\{Q\} := \frac{1}{2i} \left. \frac{f\left(Q\, \exp(t\,g)\right)}{dt} \right|_{t=0}.
```
Note that the exponential multiplies ``Q`` *on the right* — hence the name.

This operator is less common in physics, because it represents the dependence of
the function on the choice of frame (or coordinate system), which is not usually
of interest. Multiplication on the left represents a rotation of the physical
system, while rotation on the right represents a rotation of the coordinate
system.  However, this dependence on coordinate system is precisely what defines
the *spin weight* of a function, so this class of operators is relevant in
discussions of spin-weighted spherical functions.  In particular, the operators
``R_\pm`` correspond (up to a sign) to the spin-raising and -lowering operators
``\eth`` and ``\bar{\eth}`` originally introduced by [Newman and
Penrose](https://dx.doi.org/10.1063/1.1931221), as explained in greater detail
in [this paper](https://arxiv.org/abs/1604.08140).


## Commutators

In general, for generators ``a`` and ``b``, we have the commutator relations
```math
\left[ L_a, L_b \right] = -\frac{1}{2i} L_{[a,b]}
\qquad
\left[ R_a, R_b \right] = \frac{1}{2i} R_{[a,b]},
```
where ``[a,b]`` is the commutator of the two generators, which can be obtained
directly as the commutator of the corresponding quaternions.  Note the sign
difference between these two equations.  The factors of ``1/2i`` are inherited
directly from the definitions of ``L_g`` and ``R_g`` given above; the sign
difference *could* be absorbed there instead by defining those operators with
opposite signs.[^1]  The arbitrary sign choices used here are purely for
historical reasons.

[^1]:
    In fact, we can define the left and right Lie derivative operators quite
    generally, for functions on *any* Lie group and for the corresponding Lie
    algebra.  And in all cases (at least for finite-dimensional Lie algebras) we
    obtain the same commutator relations. The only potential difference is that
    it may not make sense to use the coefficient ``1/2i`` in general; it was
    chosen here for consistency with the standard angular-momentum operators.
    If that coefficient is changed in the definitions of the Lie derivatives,
    the only change to the commutator relations would the substitution of that
    coefficient.

The raising and lowering operators relative to ``L_z`` and ``R_z`` satisfy — by
definition of raising and lowering operators — the relations
```math
[L_z, L_\pm] = \pm L_\pm
\qquad
[R_z, R_\pm] = \pm R_\pm.
```
These allow us to solve — up to an overall factor — for those operators in terms
of the basic generators (again, noting the sign difference):
```math
L_\pm = L_x \pm i L_y
\qquad
R_\pm = R_x \mp i R_y.
```
In particular, this results in the commutator relations
```math
[L_+, L_-] = 2L_z
\qquad
[R_+, R_-] = 2R_z.
```
Here, the signs are *similar* because the two sign differences noted above
essentially cancel each other out.

In the functions [listed below](#Module-functions), these operators are returned
as matrices acting on vectors of mode weights.  As such, we can actually
evaluate these commutators as given to cross-validate the expressions and those
functions.


## Transformations of functions vs. mode weights

One important point to note is that mode weights transform "contravariantly"
(very loosely speaking) relative to the spin-weighted spherical functions under
some operators.  For example, take the action of the ``L_+`` operator, which
acts on a SWSH as
```math
L_+ \left\{{}_{s}Y_{\ell,m}\right\} (R) = \sqrt{(\ell-m)(\ell+m+1)} {}_{s}Y_{\ell,m+1}(R).
```
We can use this to derive the mode weights of a general spin-weighted function
``f`` under the action of this operator:[^2]
```math
\begin{aligned}
\left\{L_+ f\right\}_{\ell,m}
&=
\int \left\{L_+ f(R)\right\}\, {}_{s}\bar{Y}_{\ell,m}(R)\, dR \\
&=
\int \left\{L_+ \sum_{\ell',m'}f_{\ell',m'}\, {}_{s}Y_{\ell',m'}(R)\right\}\, {}_{s}\bar{Y}_{\ell,m}(R)\, dR \\
&=
\int \sum_{\ell',m'} f_{\ell',m'}\, \left\{L_+ {}_{s}Y_{\ell',m'}(R)\right\}\, {}_{s}\bar{Y}_{\ell,m}(R)\, dR \\
&=
\sum_{\ell',m'} f_{\ell',m'}\, \int \left\{L_+ {}_{s}Y_{\ell',m'}(R)\right\}\, {}_{s}\bar{Y}_{\ell,m}(R)\, dR \\
&=
\sum_{\ell',m'} f_{\ell',m'}\, \int \left\{\sqrt{(\ell'-m')(\ell'+m'+1)} {}_{s}Y_{\ell',m'+1}(R)\right\}\, {}_{s}\bar{Y}_{\ell,m}(R)\, dR \\
&=
\sum_{\ell',m'} f_{\ell',m'}\, \sqrt{(\ell'-m')(\ell'+m'+1)} \int {}_{s}Y_{\ell',m'+1}(R)\, {}_{s}\bar{Y}_{\ell,m}(R)\, dR \\
&=
\sum_{\ell',m'} f_{\ell',m'}\, \sqrt{(\ell'-m')(\ell'+m'+1)} \delta_{\ell,\ell'} \delta_{m,m'+1} \\
&=
f_{\ell,m-1}\, \sqrt{(\ell-m+1)(\ell+m)}
\end{aligned}
```
Note that this expression (and in particular its signs) more resembles the
expression for ``L_- \left\{{}_{s}Y_{\ell,m}\right\}`` than for ``L_+
\left\{{}_{s}Y_{\ell,m}\right\}``.  Similar relations hold for the action of
``L_-``.

[^2]:
    A technical note about the integrals above: the integrals should be taken
    over the appropriate space and with the appropriate weight such that the
    SWSHs are orthonormal.  In general, this integral should be over
    ``\mathrm{Spin}(3)`` and weighted by ``1/2\pi`` so that the result will be
    either ``0`` or ``1``; in general the SWSHs are not truly orthonormal when
    integrated over an ``S^2`` subspace (nor even is the integral invariant).
    However, if we know that the spins are the same in both cases, it *is*
    possible to integrate over an ``S^2`` subspace.

However, it is important to note that the same "contravariance" is not present
for the spin-raising and -lowering operators:
```math
\begin{aligned}
\left\{\eth f\right\}_{\ell,m}
&=
\int \left\{\eth f(R)\right\}\, {}_{s+1}\bar{Y}_{\ell,m}(R)\, dR \\
&=
\int \left\{\eth \sum_{\ell',m'}f_{\ell',m'}\, {}_{s}Y_{\ell',m'}(R)\right\}\, {}_{s+1}\bar{Y}_{\ell,m}(R)\, dR \\
&=
\sum_{\ell',m'} f_{\ell',m'}\, \int \left\{\eth {}_{s}Y_{\ell',m'}(R)\right\}\, {}_{s+1}\bar{Y}_{\ell,m}(R)\, dR \\
&=
\sum_{\ell',m'} f_{\ell',m'}\, \sqrt{(\ell'-s)(\ell'+s+1)} \int {}_{s+1}Y_{\ell',m'}(R)\, {}_{s+1}\bar{Y}_{\ell,m}(R)\, dR \\
&=
\sum_{\ell',m'} f_{\ell',m'}\, \sqrt{(\ell'-s)(\ell'+s+1)} \delta_{\ell,\ell'} \delta_{m,m'} \\
&=
f_{\ell,m}\, \sqrt{(\ell-s)(\ell+s+1)}
\end{aligned}
```
Similarly ``\bar{\eth}`` — and ``R_\pm`` of course — obey this more "covariant"
form of transformation.


## Module functions

```@autodocs
Modules = [SphericalFunctions]
Pages   = ["operators.jl"]
```
