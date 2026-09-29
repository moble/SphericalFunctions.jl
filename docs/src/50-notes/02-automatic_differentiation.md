# [Automatic differentiation](@id automatic_differentiation)

```@meta
CurrentModule = SphericalFunctions
```

Wigner's ``𝔇`` matrices and the spin-weighted spherical harmonics are
smooth functions of the rotor everywhere on ``\mathrm{Spin}(3)``, so
it is natural to differentiate them with respect to the rotor by
automatic differentiation.  For [`D`](@ref) and [`sYlm`](@ref) of a
single rotor, and for the functions built on them, this package
supplies rules that give the derivatives directly, in terms of the
values themselves.  These rules are used by `ForwardDiff`,
`ReverseDiff`, `Enzyme`, and `Mooncake`, and by the tools that read
`ChainRules`, such as `Zygote`.  Elsewhere — in the calculators, for
example — automatic differentiation differentiates the algorithm, not
the function, and the algorithm used here passes through intermediate
quantities that are singular at two special sets of rotors: those that
take the ``z`` axis to itself, and those that take it to its opposite
— or those with ``β = 0`` or ``β = π``.  These are the rotors at which
the harmonics are evaluated at the poles of the sphere, so we will
refer to both sets as "poles."  This note describes the rules first,
and then explains where the singularities come from, how the
calculators work around them for ``𝔇`` and the harmonics, and why the
same cannot be done for ``d`` and ``H`` of a rotor.


## Rules for the derivatives

The derivative of ``𝔇`` along a rotation is given by the generators
of rotations, which are the angular-momentum operators.  Along the
path ``𝐑(t) = e^{t𝐮/2}\, 𝐑``, for any vector ``𝐮``,
```math
\frac{d}{dt} 𝔇^{(ℓ)}(𝐑(t)) \bigg|_{t=0}
=
-i\, (𝐮 ⋅ 𝐉)\, 𝔇^{(ℓ)}(𝐑),
```
where ``𝐉`` is the angular momentum acting on the index ``m'``, with
``⟨m'|J_z|m'⟩ = m'`` and ``⟨m'±1|J_±|m'⟩ = \sqrt{(ℓ∓m')(ℓ±m'+1)}``.
The functions are taken to depend on the rotor only through
``𝐑/\|𝐑\|``, which is how they are computed, so a tangent
``\dot{𝐑}`` may be any quaternion.  Writing ``\dot{𝐑} = 𝐪\, 𝐑``,
with ``𝐪 = \dot{𝐑}\, \bar{𝐑} / \|𝐑\|^2``, the scalar part of
``𝐪`` changes only the norm of ``𝐑``, and drops out, while its
vector part ``𝐯`` gives ``𝐮 = 2𝐯``.  So, with ``w = v_x + i v_y``,
```math
\dot{𝔇}^{(ℓ)}_{m',m}
=
-i \left[
2 v_z\, m'\, 𝔇^{(ℓ)}_{m',m}
+ \bar{w} \sqrt{(ℓ-m'+1)(ℓ+m')}\, 𝔇^{(ℓ)}_{m'-1,m}
+ w \sqrt{(ℓ+m'+1)(ℓ-m')}\, 𝔇^{(ℓ)}_{m'+1,m}
\right].
```
The harmonics are a conjugated row of ``𝔇``, so the same derivative,
conjugated, applies to them, and couples each harmonic only to those
of the same ``ℓ`` and spin weight with ``m ± 1``.  Differentiating
from the left in this way is what keeps the spin weight fixed; from
the right, the derivative would couple the harmonics of weight ``s``
to those of weights ``s ± 1``.  For a block of ``𝔇`` restricted in
``m'``, the derivatives need the values one row beyond each limit,
which the rules compute along with the block.  The reverse-mode rules
apply the adjoint of this linear map, and return a cotangent that is
orthogonal to ``𝐑``, as the cotangent of a function of ``𝐑/\|𝐑\|``
must be.

The derivative in every direction is therefore a combination of values
of the same ``ℓ``, and is as accurate as the values are, at every
rotor, the poles included.  Because the rules compute those values by
calling the same function again, a tool that nests its derivatives, as
`ForwardDiff` does for a Hessian, reaches the rules once at each
level, and every order of derivative is exact.  The recurrence itself
is never differentiated.

The rules are supplied by package extensions, which are loaded along
with the tool: for `ChainRulesCore`, `EnzymeCore`, `ForwardDiff`,
`Mooncake`, and `ReverseDiff`.  They apply to [`D`](@ref) of a rotor
or of Euler angles, and to [`sYlm`](@ref) and [`Ylm`](@ref) of a rotor
or of spherical coordinates, since each of these reaches the same
underlying function of a single rotor.  The calculators, the forms
that take a vector of rotors, [`sYlm_matrix`](@ref), and the
transforms are differentiated through the algorithm, as described in
the rest of this note.


## The singularity in the recurrence

We write a rotor as ``𝐑 = W + X𝐢 + Y𝐣 + Z𝐤``.  The [``H``
recurrence](@ref "Algorithm for computing ``H``") does not see the
rotor directly.  Instead, it works with the half-angles
```math
\cos\frac{β}{2} = \frac{\sqrt{W^2 + Z^2}}{\|𝐑\|}
\qquad \text{and} \qquad
\sin\frac{β}{2} = \frac{\sqrt{X^2 + Y^2}}{\|𝐑\|},
```
and with the phases
```math
z_+ = e^{i(α+γ)/2} = \frac{W + iZ}{\sqrt{W^2 + Z^2}}
\qquad \text{and} \qquad
z_- = e^{i(α-γ)/2} = \frac{Y - iX}{\sqrt{X^2 + Y^2}}.
```
The half-angles determine ``d``, and the phases supply the factors
that convert ``d`` into ``𝔇``, because ``e^{i(m'α + mγ)} =
z_+^{m'+m}\, z_-^{m'-m}``, in which the exponents are integers even
when the indices are half-integers.

At ``β = 0``, where ``X = Y = 0``, the square root in ``\sin(β/2)`` is
taken of zero, and the derivative of the square root is infinite
there; the phase ``z_-`` is not even defined.  Automatic
differentiation therefore returns `NaN` for every derivative of ``𝔇``
computed this way, even though ``𝔇`` itself is perfectly smooth
there.  The same happens at ``β = π`` with ``\cos(β/2)`` and ``z_+``.

Near a pole the problem is less obvious, but almost as bad.  At a
distance ``r`` from the pole — by which we mean ``r = \sin(β/2)`` near
``β = 0``, or ``r = \cos(β/2)`` near ``β = π`` — the ``k``-th
derivatives of the half-angle and the phase are of order ``r^{-k}``.
They cancel when the two are multiplied together, but the cancellation
leaves an error of order ``ε\, r^{-k}`` relative to the size of the
``k``-th derivative of ``𝔇``, where ``ε`` is the machine epsilon.
The values themselves are unaffected; only the derivatives suffer.


## The expansion about a pole

The quantities that are smooth everywhere are the products
```math
σ = \cos\frac{β}{2}\, z_+ = \frac{W + iZ}{\|𝐑\|}
\qquad \text{and} \qquad
ρ = \sin\frac{β}{2}\, z_- = \frac{Y - iX}{\|𝐑\|},
```
which are the rotor's normalized Cayley–Klein parameters.  Wigner's
formula for ``d``, with the phases of our convention absorbed into
``σ`` and ``ρ``, expresses ``𝔇`` as a polynomial in these parameters
and their conjugates:
```math
𝔇^{(ℓ)}_{m',m}(𝐑)
=
\sum_s (-1)^{k+s} C_s\,
\bar{σ}^{ℓ+m-s}\, σ^{ℓ-m'-s}\, \bar{ρ}^{k+s}\, ρ^s,
\qquad
k = m' - m,
```
where
```math
C_s^2
=
\binom{ℓ+m}{s} \binom{ℓ-m'}{s} \binom{ℓ+m'}{k+s} \binom{ℓ-m}{k+s},
```
and the sum runs over every ``s`` for which all of the exponents are
non-negative.  Every exponent is an integer, even for half-integer
indices.  This polynomial is of no use as a general-purpose algorithm:
away from the poles its terms cancel badly, and for large ``ℓ`` its
coefficients overflow.  Near a pole, however, it is exactly what is
needed.

At ``β = 0`` the parameter ``ρ`` vanishes, and at ``β = π`` the
parameter ``σ`` vanishes.  A term of degree ``n`` in the vanishing
parameter and its conjugate vanishes at the pole together with all of
its derivatives of order less than ``n``.  So every derivative of
order ``N`` or less at the pole is given *exactly* by the terms of
degree at most ``N``.  There are at most ``\lfloor N/2 \rfloor + 1``
such terms in any element, and they appear only in the elements with
``|m'-m| ≤ N`` near ``β = 0``, or ``|m'+m| ≤ N`` near ``β = π``; every
other element vanishes to that order.  The remaining factor is a power
of the parameter that does not vanish, which is written as its modulus
times one of the phases ``z_\pm``, both of which are smooth near that
pole.  Near the pole, the truncated sum is accurate to the size of the
first term omitted, which is about ``\left((ℓ+1)\, r\right)^{N+1-k}``
relative to the size of a ``k``-th derivative.

This package keeps the terms up to degree ``N = 8``.  A rotor within a
distance
```math
r_s = \frac{ε^{1/(N+1)}}{ℓₘₐₓ+1}
```
of either pole is evaluated from the truncated expansion, in place of
the recurrence.  In `Float64` this corresponds to an angle of about
``0.04/(ℓₘₐₓ+1)``.  The radius ``r_s`` is the largest at which the
omitted terms cannot change the values, so the values are as accurate
as the recurrence's.  At the pole itself, every derivative up to
eighth order is exact, and only a derivative of higher order would be
wrong.  Just outside ``r_s``, where the recurrence is used again, the
recurrence's ``k``-th derivatives are wrong by at most about
``ε^{1-k/(N+1)}\, (ℓₘₐₓ+1)^k`` relative to their size.  In `Float64`
with ``ℓₘₐₓ = 32``, the errors measured there are of order
``10^{-13}`` for first derivatives and ``10^{-10}`` for second
derivatives.  A Hessian needs only ``N = 2``, of course, but the order
also sets the radius, and a larger order widens the neighborhood in
which the recurrence's inaccurate derivatives are replaced.  The cost
is negligible, since at most seventeen elements of any column are
nonzero in the expansion, each with at most five terms.  The radius is
set by ``ℓₘₐₓ`` rather than by each ``ℓ``, so that a given rotor is
treated in the same way at every ``ℓ``.

The test for whether a rotor is near a pole is applied to every number
type, not just to dual numbers.  This is necessary because tools like
`Enzyme` differentiate the ordinary floating-point code, and cannot be
distinguished from an ordinary evaluation by the type of the input.
Also, a rotor that is exactly at a pole is given to the recurrence as
the exact constants ``e^{iβ} = ±1`` and the corresponding half-angles,
rather than through the square root of zero.  The recurrence's results
for that rotor are overwritten by the expansion anyway, but a
reverse-mode tool like `ReverseDiff` runs its reverse pass through
every operation it recorded, including those whose results were later
discarded, and would otherwise encounter ``0 × ∞``.

None of this is needed when the recurrence is given an angle rather
than a rotor — as it is for ``d`` of an angle ``β`` or of a phase
``e^{iβ}``, and for [the real harmonics](@ref
interface_real_harmonics) ``{}_{s}λ_{ℓ,m}(θ)`` — because the
recurrence is a smooth function of that angle at the poles.


## ``d`` and ``H`` of a rotor

The expansion is not used for ``d`` and ``H`` of a rotor, and their
derivatives at the poles are `NaN`.  This is not a limitation of the
algorithm; for some of the elements, those derivatives simply do not
exist.  Although ``d`` is a smooth function of ``β``, ``β`` is *not* a
smooth function of the rotor at the poles.  Near the identity, for
example, ``\sin(β/2) = \sqrt{X^2 + Y^2} / \|𝐑\|``, which is a cone
over the ``XY`` plane — the two-dimensional analog of ``|x|``.  In
particular, the rotations by angles ``ε`` and ``-ε`` about the ``x``
axis both have ``β = |ε|``.

Whether a given element of ``d`` survives this depends on its parity
in ``β``.  Each element satisfies
```math
d^{(ℓ)}_{m',m}(-β) = (-1)^{m'-m}\, d^{(ℓ)}_{m',m}(β).
```
If ``m'-m`` is even, the element is an even function of ``β``, and is
therefore a smooth function of ``\sin^2(β/2) = (X^2 + Y^2) /
\|𝐑\|^2``, which is itself smooth in the rotor — just as ``\cos|x|``
is simply ``\cos x``.  If ``m'-m`` is odd, on the other hand, the
element is ``\sin(β/2)`` times a smooth function of ``\sin^2(β/2)``,
and so inherits the cone.  For example,
```math
d^{(1)}_{1,0}(β) = -\frac{\sin β}{\sqrt{2}},
```
which, along the rotations about the ``x`` axis, behaves like
``-|ε|/\sqrt{2}``.  Its one-sided derivatives in opposite directions
have opposite signs, so it has no derivative at the identity, and
`NaN` is the correct answer.  More generally, an element with
``|m'-m| = 3`` behaves like ``|ε|^3``, which has a first derivative
but no second; and so on.  The same thing happens at ``β = π`` with
``\cos(β/2) = \sqrt{W^2 + Z^2} / \|𝐑\|`` and the parity of ``m'+m``.

``𝔇`` escapes this problem because it multiplies the same elements by
phases that are also singular at the pole, and the two singularities
cancel.  The odd power of ``\sin(β/2)`` combines with a power of the
phase ``z_-`` to form a power of ``ρ`` or ``\bar{ρ}``, which is linear
in the rotor's components.  Seen this way, ``d`` is essentially ``𝔇``
with its phases removed, and removing the phases replaces ``ρ =
(Y - iX)/\|𝐑\|`` with its modulus ``|ρ|``, which is precisely the
cone.

To differentiate ``d`` at the poles, therefore, we must give it the
angle ``β`` or the phase ``e^{iβ}``, with respect to which it is
smooth everywhere.
