# [Automatic differentiation](@id automatic_differentiation)

```@meta
CurrentModule = SphericalFunctions
```

Wigner's ``𝔇`` matrices and ``ₛY_{ℓ,m}`` are smooth functions of
their arguments everywhere on ``\mathrm{Spin}(3)``, so it is natural
to differentiate them with respect to the rotor by automatic
differentiation.  But naive automatic differentiation isn't
necessarily the most efficient way, and in fact it can fail at
important points.  Both *functions* are smooth everywhere, with finite
derivatives of all orders.  However, automatic differentiation
necessarily acts on *the algorithms we use to compute the functions*,
rather than the functions themselves.  This package uses an extremely
efficient and accurate recurrence that essentially computes the
``d(β)`` function — or more precisely its extension to rotors ``d(R)``
— and then multiplies by the appropriate complex phases.  While this
is ideal for most purposes, it is not ideal for automatic
differentiation with very specific arguments.

To see why, note that ``d`` is essentially the complex magnitude of
``𝔇``.  And just as the derivative of ``|z|`` is undefined at ``z =
0``, the derivative of ``d(R)`` can also also be undefined at certain
values equivalent to ``β = 0`` or ``β = π``.  More precisely, ``d(β)``
is perfectly differentiable at those values, but the rotor extension
``d(R)`` is not differentiable when ``R𝐳R^{-1} = ±1``.

On the other hand, ``|z|^k`` is perfectly differentiable for ``k≠1``.
In the exact same way, it turns out that ``d^{(ℓ)}_{m',m}`` is
differentiable at ``β = 0`` for ``|m'-m| ≠ 1``, and at ``β = π`` for
``|m'+m| ≠ 1``.  Now, because ``𝔇^{(ℓ)}_{m',m}(R)`` and
``ₛY_{ℓ,m}(R)`` are computed *via* ``d^{(ℓ)}_{m',m}(β)``, the
derivatives of these functions with respect to the rotor are also
undefined at these points.  This is an unavoidable consequence of the
choice of algorithms used to compute ``𝔇``.

To solve this problem, we have a few options.  One option would be to
switch to a different algorithm for problematic values of ``R``.  For
example, we know how to express ``𝔇`` as a polynomial in the
quaternion components, which is easily differentiated, and actually
not such a bad method for ``R`` near the "poles".  However, this is
inefficient and requires fine-tuning to determine when to switch
algorithms.

A better option is to use our knowledge of the derivatives of ``𝔇``
to compute them explicitly in terms of the values of ``𝔇``.  This is
the approach taken by this package, as it is faster and more accurate
than automatic differentiation of the algorithm.  The only downside is
that it requires some extra work to implement the rules for the
variety of automatic differentiation tools that exist in the Julia
ecosystem: `ChainRules`, `ForwardDiff`, `Enzyme`, `Mooncake`, and
`ReverseDiff`.

This note describes the rules first, and then goes into more detail
about the problematic values noted above.


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
This is the derivative from the left, which couples each element to
its neighbors in the same column.  Writing instead ``\dot{𝐑} = 𝐑\,
𝐪'``, with ``𝐪' = \bar{𝐑}\, \dot{𝐑} / \|𝐑\|^2`` and ``𝐯'`` its
vector part, gives the derivative from the right, which couples each
element to its neighbors in the same row:
```math
\dot{𝔇}^{(ℓ)}_{m',m}
=
-i \left[
2 v'_z\, m\, 𝔇^{(ℓ)}_{m',m}
+ w' \sqrt{(ℓ-m+1)(ℓ+m)}\, 𝔇^{(ℓ)}_{m',m-1}
+ \bar{w}' \sqrt{(ℓ+m+1)(ℓ-m)}\, 𝔇^{(ℓ)}_{m',m+1}
\right].
```
A block of ``𝔇`` whose rows include all of ``-ℓ:ℓ`` is differentiated
from the left, and one whose columns do, from the right, so that every
neighbor needed is in the block already.  Only a block restricted in
both ``m'`` and ``m`` depends on values beyond its limits.  A
calculator of such blocks stores only its block, but its wedge reaches
one more row or column on each side, along whichever axis is the
wider, and those values are computed when a derivative needs them: by
the calculator of values inside a calculator whose rotors hold
derivatives, and by `derivative_values` in the rules of the other
tools.  The recurrence, whose cost is set by the narrower axis, is
widened only when the two are equally wide.  The harmonics are a
conjugated row of ``𝔇``, so the derivative from the left, conjugated,
applies to them, and couples each harmonic only to those of the same
``ℓ`` and spin weight with ``m ± 1``; from the right, the derivative
would couple the harmonics of weight ``s`` to those of weights ``s ±
1``.  The reverse-mode rules apply the adjoint of these linear maps,
and give a cotangent that is orthogonal to ``𝐑``, as the cotangent of
a function of ``𝐑/\|𝐑\|`` must be.

The rotations that ``d`` and ``{}_{s}λ_{ℓ,m}`` describe are those
about the ``y`` axis, ``𝐑(β) = e^{β𝐣/2}``, which commute with their
generator, so that a tangent ``\dot{β}`` is the vector ``𝐯 = (0,
\dot{β}/2, 0)`` from the left and from the right alike.  The formulas
above then become real:
```math
\dot{d}^{(ℓ)}_{m',m}
=
\frac{\dot{β}}{2} \left[
b_{m'}\, d^{(ℓ)}_{m'+1,m} - a_{m'}\, d^{(ℓ)}_{m'-1,m}
\right]
=
\frac{\dot{β}}{2} \left[
a_m\, d^{(ℓ)}_{m',m-1} - b_m\, d^{(ℓ)}_{m',m+1}
\right],
\qquad
{}_{s}\dot{λ}_{ℓ,m}
=
\frac{\dot{θ}}{2} \left[
b_m\, {}_{s}λ_{ℓ,m+1} - a_m\, {}_{s}λ_{ℓ,m-1}
\right],
```
with ``a_m = \sqrt{(ℓ-m+1)(ℓ+m)}`` and ``b_m = \sqrt{(ℓ+m+1)(ℓ-m)}``.
The two forms of ``\dot{d}`` are equal, so the side is chosen, as for
``𝔇``, by the neighbors that a block holds, and otherwise for speed.
Each of these derivatives is a difference of two products that nearly
cancel near the poles, and a gradient sums many of them with weights
of either sign.  So the coefficients are used to about twice the
working precision, since the rounding error of a coefficient is the
same in every element that it multiplies, and the reverse-mode rules
add their terms with compensation.  A phase ``e^{iβ}`` is
differentiated through its argument, so that only the part of its
tangent along the unit circle has an effect, as only the part of a
rotor's tangent orthogonal to the rotor does; a phase must have unit
modulus in any case, or the calculator refuses it.  A rotor is
differentiated through its ``β``, as described below.

The derivative of a block in every direction is therefore a
combination of values of that same block, and is as accurate as the
values are, at every rotor, the poles included.  This is what lets a
calculator produce the derivatives of each block as it produces the
block, one ``ℓ`` at a time, in the loop over blocks that the
calculators are designed for.  A calculator whose rotors or angles are
`ForwardDiff`'s dual numbers runs the recurrence on the values of
those rotors, and writes each value together with its derivatives into
its blocks of dual numbers; stepping it allocates nothing.  Under
nested differentiation, as for a Hessian, those values are themselves
lifted from a calculator of their own values, so every order of
derivative is exact.  `Enzyme` and `Mooncake` differentiate a
calculator of floats: the calculator keeps a copy of its rotors, or of
its angles, whose tangents or cotangents those tools follow, and the
rules for each step give the block's derivatives from them, or add the
block's cotangents into them.  A calculator that `Enzyme` is to
differentiate must be `Duplicated`, as any mutable workspace must be,
which it is automatically when it is created within the function being
differentiated.  A calculator of `ReverseDiff`'s tracked numbers, like
one of dual numbers, runs the recurrence on the values of its rotors
or angles, and it records each step as one instruction, whose pullback
adds the cotangents of the block into those of the rotors or angles.
In every case the recurrence itself is never differentiated.  A tape
that `ReverseDiff` records, and then compiles or replays at other
points, gives the derivatives of [`D`](@ref), [`d`](@ref),
[`sYlm`](@ref), and [`sYlm_matrix`](@ref) correctly at each of them,
and so does a tape that loops over a calculator, since each step's
instruction computes its block from the rotors or angles it is given
when it is replayed.  On Julia 1.10, `Enzyme`'s reverse mode fails to
differentiate a loop over a calculator when bounds checking is forced
on, as it is by `Pkg.test`, because its type analysis cannot deduce
the type of an integer in that loop; its forward mode, and every other
tool, is unaffected.

The rules are supplied by package extensions, which are loaded along
with the tool: for `ChainRulesCore`, `EnzymeCore`, `ForwardDiff`,
`Mooncake`, and `ReverseDiff`.  `Zygote` cannot follow the mutation of
a calculator's buffers, so it uses rules for the arrays of
[`D`](@ref), [`d`](@ref), [`sYlm`](@ref), and [`sYlm_matrix`](@ref)
instead, which also serve [`Ylm`](@ref) and the forms of those
functions that take Euler angles or spherical coordinates; it cannot
differentiate a loop over a calculator at all.  `ReverseDiff` uses
rules for those arrays too, which record one instruction for a whole
array rather than one for each block.


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

The rules avoid both problems, because they never differentiate this
split: the derivatives of a block are formed from its values, which
are as accurate at and near a pole as anywhere else.  None of this
arises when the recurrence is given an angle rather than a rotor — as
it is for ``d`` of an angle ``β`` or of a phase ``e^{iβ}``, and for
[the real harmonics](@ref interface_real_harmonics)
``{}_{s}λ_{ℓ,m}(θ)`` — because the recurrence is then a smooth
function of that angle at the poles.  Even so, the rules for ``d`` and
``{}_{s}λ_{ℓ,m}`` never differentiate the recurrence either.


## ``d`` and ``H`` of a rotor

``d`` of a rotor is differentiated through the angle ``β`` of the
rotor's Euler decomposition, ``2\operatorname{atan}(\sqrt{X^2+Y^2},
\sqrt{W^2+Z^2})``, by the same rules as ``d`` of an angle, and is as
accurate away from the poles.  ``H`` of a rotor has no rules.  At the
poles, some of the derivatives of ``d`` of a rotor are `NaN`.  This is
not a limitation of the algorithm; for those elements, the derivatives
simply do not exist.  Although ``d`` is a smooth function of ``β``,
``β`` is *not* a smooth function of the rotor at the poles.  Near the
identity, for example, ``\sin(β/2) = \sqrt{X^2 + Y^2} / \|𝐑\|``,
which is a cone over the ``XY`` plane — the two-dimensional analog of
``|x|``.  In particular, the rotations by angles ``ε`` and ``-ε``
about the ``x`` axis both have ``β = |ε|``.

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

At a pole, the rules give each element of ``d`` of a rotor the
derivative it has there as a function of the rotor: zero, unless
``m'-m`` (at ``β = 0``) or ``m'+m`` (at ``β = π``) is ``±1``, in which
case it has none, and is `NaN`.

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
