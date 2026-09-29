# [Automatic differentiation](@id automatic_differentiation)

```@meta
CurrentModule = SphericalFunctions
```

Wigner's ``𝔇`` matrices and the spin-weighted spherical harmonics are
smooth functions of the rotor everywhere on ``\mathrm{Spin}(3)``, so
it is natural to differentiate them with respect to the rotor by
automatic differentiation.  This package supplies rules that give the
derivatives directly, in terms of the values themselves.
`ForwardDiff`, `Enzyme`, `Mooncake`, and `ReverseDiff` use them for
every step of the calculators, and so for [`D`](@ref), [`sYlm`](@ref),
and [`sYlm_matrix`](@ref), which the calculators compute; the tools
that read `ChainRules`, such as `Zygote`, use them for those three
functions.  A rotor may be given as a `Rotor` or as any other
`Quaternion`, which denotes the rotation of its normalization, so that
the derivatives may be taken with respect to the four components of an
unnormalized quaternion directly.  The rules are needed for accuracy,
and not only for speed.  Where no rule applies, automatic
differentiation differentiates the algorithm, not the function, and
the algorithm used here passes through intermediate quantities that
are singular at two special sets of rotors: those that take the ``z``
axis to itself, and those that take it to its opposite — or those with
``β = 0`` or ``β = π``.  These are the rotors at which the harmonics
are evaluated at the poles of the sphere, so we will refer to both
sets as "poles."  This note describes the rules first, and then
explains where the singularities come from, and why ``d`` and ``H`` of
a rotor cannot be protected from them in the same way.


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
both ``m'`` and ``m`` needs values beyond its limits, and a calculator
of such blocks computes one more row or column on each side, along
whichever axis is the wider, so that the recurrence, whose cost is set
by the narrower, is widened only when the two are equally wide.  The
harmonics are a conjugated row of ``𝔇``, so the derivative from the
left, conjugated, applies to them, and couples each harmonic only to
those of the same ``ℓ`` and spin weight with ``m ± 1``; from the
right, the derivative would couple the harmonics of weight ``s`` to
those of weights ``s ± 1``.  The reverse-mode rules apply the adjoint
of these linear maps, and give a cotangent that is orthogonal to
``𝐑``, as the cotangent of a function of ``𝐑/\|𝐑\|`` must be.

The derivative of a block in every direction is therefore a
combination of values of that same block, and is as accurate as the
values are, at every rotor, the poles included.  This is what lets a
calculator produce the derivatives of each block as it produces the
block, one ``ℓ`` at a time, in the loop over blocks that the
calculators are designed for.  A calculator whose rotors are
`ForwardDiff`'s dual numbers runs the recurrence on the values of
those rotors, and writes each value together with its derivatives into
its blocks of dual numbers; stepping it allocates nothing.  Under
nested differentiation, as for a Hessian, those values are themselves
lifted from a calculator of their own values, so every order of
derivative is exact.  `Enzyme` and `Mooncake` differentiate a
calculator of floats: the calculator keeps a copy of its rotors, whose
tangents or cotangents those tools follow, and the rules for each step
give the block's derivatives from them, or add the block's cotangents
into them.  A calculator that `Enzyme` is to differentiate must be
`Duplicated`, as any mutable workspace must be, which it is
automatically when it is created within the function being
differentiated.  A calculator of `ReverseDiff`'s tracked numbers, like
one of dual numbers, runs the recurrence on the values of its rotors,
and it records each step as one instruction for each rotor, whose
pullback adds the cotangents of that rotor's block into the rotor's.
In every case the recurrence itself is never differentiated.
`ReverseDiff` replays a recorded tape by running its instructions
again, but keeps the pullbacks of the first run, so a tape recorded at
one rotor and replayed at another gives the derivatives at the first;
as for any rule defined with `ReverseDiff`'s `@grad`, each gradient
should record its own tape.  On Julia 1.10, `Enzyme`'s reverse mode
fails to differentiate a loop over a calculator when bounds checking
is forced on, as it is by `Pkg.test`, because its type analysis cannot
deduce the type of an integer in that loop; its forward mode, and
every other tool, is unaffected.

The rules are supplied by package extensions, which are loaded along
with the tool: for `ChainRulesCore`, `EnzymeCore`, `ForwardDiff`,
`Mooncake`, and `ReverseDiff`.  `Zygote` cannot follow the mutation of
a calculator's buffers, so it uses rules for the arrays of
[`D`](@ref), [`sYlm`](@ref), and [`sYlm_matrix`](@ref) instead, which
also serve [`Ylm`](@ref) and the forms of those functions that take
Euler angles or spherical coordinates; it cannot differentiate a loop
over a calculator at all.  `ReverseDiff` uses rules for those arrays
too, which record one instruction for a whole array rather than one
for each block.


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
function of that angle at the poles.


## ``d`` and ``H`` of a rotor

The rules do not apply to ``d`` and ``H`` of a rotor, and their
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
