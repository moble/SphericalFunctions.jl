# [Differential operators](@id background_differential_operators)

On [the previous page](@ref background_domain), we saw that the domain
on which our most important functions — spherical harmonics and Wigner
𝔇 matrices — are defined is most naturally regarded as the group of
unit quaternions, ``\mathrm{Spin}(3) \cong \mathrm{SU}(2)``.  As a
result, we can define a variety of differential operators acting on
these functions, relating to infinitesimal motions in this group,
acting either from the left or the right on their arguments.  Right or
left matters because this group is non-commutative.  The left
derivative operator we define here turns out to be exactly the
angular-momentum operator familiar from physics.  In the coming pages,
we will use these operators to define the spin-weighted spherical
harmonics and Wigner's 𝔇 matrices as eigenfunctions.

## Definitions and properties

In general, the *left* Lie derivative of a function ``f(𝐐)`` over the
unit quaternions with respect to a generator of rotation ``𝐠`` is
defined as
```math
L_𝐠(f)\{𝐐\} := -\frac{i}{2}
    \left. \frac{df\left(e^{ϵ\,𝐠}\, 𝐐\right)}{d ϵ} \right|_{ϵ=0}.
```
Note that the exponential multiplies ``𝐐`` *on the left* — hence the
name.  We will see below that this agrees with the usual definition of
the angular-momentum from physics, except that in *quantum* physics a
factor of ``\hbar`` is usually included.

So, for example, a rotation about the ``𝐳`` axis has the quaternion
``𝐤`` as its generator of rotation, and ``L_𝐤`` defined in this way
agrees with [the usual angular-momentum
operator](https://en.wikipedia.org/wiki/Angular_momentum_operator)
``L_z`` familiar from spherical-harmonic theory, and reduces to it
when the function has spin weight 0, but also applies to functions of
general spin weight.  Similarly, we can compute ``L_x`` and ``L_y``,
and take appropriate combinations to find [the usual raising and
lowering (ladder)
operators](https://en.wikipedia.org/wiki/Ladder_operator#Angular_momentum)
``L_+`` and ``L_-``.

In just the same way, we can define the *right* Lie derivative of a
function ``f(𝐐)`` over the unit quaternions with respect to a
generator of rotation ``𝐠`` as
```math
R_𝐠(f)\{𝐐\} := -\frac{i}{2}
    \left. \frac{df\left(𝐐\, e^{-ϵ\,𝐠}\right)}{d ϵ} \right|_{ϵ=0}.
```
Note that the exponential multiplies ``𝐐`` *on the right* — hence the
name.

This operator is less common in physics, because it represents the
dependence of the function on the choice of frame (or coordinate
system), which is not usually interesting.  Multiplication on the left
represents a rotation of the physical system, while rotation on the
right represents a rotation of the coordinate system.  However, this
dependence on coordinate system is precisely what defines the *spin
weight* of a function, so this class of operators is relevant in
discussions of spin-weighted spherical functions.  In particular,
``R_z`` is the spin-weight operator — meaning that when it acts on a
spin-weighted spherical harmonic of spin weight ``s``, it returns
``s`` times that same function.  Moreover, the operators ``R_\pm``
correspond (up to a sign) to the spin-raising and -lowering operators
``\eth`` and ``\bar{\eth}`` originally introduced by
[Newman_1966](@citet), as explained in greater detail by
[Boyle_2016](@citet).

Note that these definitions are *extremely* general, in that they can
be used for *any* Lie group, and for any complex-valued function on
that group.  And in full generality, we have the useful properties of
linearity:
```math
L_{s𝐚} = sL_{𝐚}
\qquad \text{and} \qquad
R_{s𝐚} = sR_{𝐚},
```
and
```math
L_{𝐚+𝐛} = L_{𝐚} + L_{𝐛}
\qquad \text{and} \qquad
R_{𝐚+𝐛} = R_{𝐚} + R_{𝐛},
```
for any scalar ``s`` and any elements of the Lie algebra
``𝐚`` and ``𝐛``.  In particular, if the Lie algebra
has a basis ``𝐞_{(j)}``, we use the shorthand ``L_j`` and
``R_j`` for ``L_{𝐞_{(j)}}`` and ``R_{𝐞_{(j)}}``,
respectively, and we can expand any operator in terms of these basis
operators:
```math
L_{𝐚} = \sum_{j} a_j L_j
\qquad \text{and} \qquad
R_{𝐚} = \sum_{j} a_j R_j.
```


## Commutators and angular momentum

In general, for generators ``𝐚`` and ``𝐛``, we have the commutator
relations
```math
\left[ L_𝐚, L_𝐛 \right] = \frac{i}{2} L_{[𝐚,𝐛]}
\qquad
\left[ R_𝐚, R_𝐛 \right] = \frac{i}{2} R_{[𝐚,𝐛]},
```
where ``[𝐚,𝐛]`` is the commutator of the two generators, which can
be obtained directly as the commutator of the corresponding
quaternions.  Note that these two equations have the same signs.  The
factors of ``i/2`` are inherited directly from the definitions of
``L_𝐠`` and ``R_𝐠`` given above.  Note the subtle sign difference in
the exponents in those definitions.  The fact that these two
commutator relations have the same sign results from the fact that the
quaternions are multiplied in opposite orders in the two cases.  There
are overall arbitrary sign choices; we choose these purely for
conventional reasons, to reproduce the standard ``L`` operator and to
produce similar commutators for ``R``.[^1]

[^1]:
    In fact, we can define the left and right Lie derivative operators
    quite generally, for functions on *any* Lie group and for the
    corresponding Lie algebra.  And in all cases (at least for
    finite-dimensional Lie algebras) we obtain the same commutator
    relations. The only potential difference is that it may not make
    sense to use the coefficient ``i/2`` in general; it was chosen
    here for consistency with the standard angular-momentum operators.
    If that coefficient is changed in the definitions of the Lie
    derivatives, the only change to the commutator relations would the
    substitution of that coefficient.  The presence of an ``i`` is
    important to ensure that the operators are Hermitian when acting on
    appropriate function spaces.

Again, these results are valid for general (finite-dimensional) Lie
groups, but a particularly interesting case is in application to the
three-dimensional rotation group.  In the following, we will apply our
results to this group.

The commutator relations for ``L`` are consistent — except for the
differing use of ``\hbar`` — with the usual relations from quantum
mechanics:
```math
\left[ L_j, L_k \right] = i \hbar \sum_{l=1}^{3} ε_{jkl} L_l.
```
Here, ``j``, ``k``, and ``l`` are indices that run from 1 to 3, and
index the set of basis vectors ``(\hat{x}, \hat{y}, \hat{z})``.  If we
represent an arbitrary basis vector as ``\hat{e}_j``, then the
quaternion commutator ``[𝐚,𝐛]`` in the expression for ``[L_𝐚, L_𝐛]``
becomes
```math
[\hat{e}_j, \hat{e}_k] = 2 \sum_{l=1}^{3} ε_{jkl} \hat{e}_l.
```
That is, the commutator ``[𝐚,𝐛]`` is essentially twice the cross
product of the corresponding vectors.  Plugging this into the general
expression ``[L_𝐚, L_𝐛] = \frac{i}{2} L_{[𝐚,𝐛]}``, we obtain
(except for that factor of ``\hbar``) the version used in quantum
physics.

The raising and lowering operators relative to ``L_z`` and ``R_z``
satisfy — *by definition of raising and lowering operators* — the
relations
```math
[L_z, L_\pm] = \pm L_\pm
\qquad
[R_z, R_\pm] = \pm R_\pm.
```
These allow us to solve, up to an overall factor, for those operators
in terms of the basic generators:
```math
L_\pm = L_x \pm i L_y
\qquad
R_\pm = R_x \pm i R_y.
```
(Interestingly, the solution process also shows that rasing and
lowering operators can only exist if the factor in front of the
derivatives in the definitions of ``L_g`` and ``R_g`` are pure
imaginary numbers.) In particular, this results in the commutator
relations
```math
[L_+, L_-] = 2L_z
\qquad
[R_+, R_-] = 2R_z.
```

In the functions [listed below](#Module-functions), these operators
are returned as matrices acting on vectors of mode weights.  As such,
we can actually evaluate these commutators as given to cross-validate
the expressions and those functions.
