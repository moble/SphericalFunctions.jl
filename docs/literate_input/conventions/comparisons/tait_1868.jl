md"""
# Tait (1867 and 1868)

[Tait_1867](@cite) provided an elementary introduction to quaternions, which was full of
examples and applications.  He was quite disinterested in mathematics *per se*; his focus
was on applying mathematics in general — and quaternions specifically — to physics.  Recall
that, at this time, vector algebra had not yet been developed, so quaternions were one of
the few available tools for dealing with three-dimensional quantities.  (The tragic
quaternion-vector wars would play out over the next few decades, even though Clifford really
unified the two notions.)

But, as such, this book was a broadly influential introduction to quaternions for many
years.

In it, he defines the quaternion basis elements ``(i, j, k)`` with multiplication rules
```math
\begin{gathered}
i^2 = j^2 = k^2 = -1,
\qquad
ij = -ji = k,
\qquad
jk = -kj = i,
\qquad
ki = -ik = j,
\\
ijk = -1.
\end{gathered}
```
Note that Tait specifically writes expressions like ``xi + yj + zk`` and ``\xi i + \eta j +
\zeta k``; these suggest that Tait does think of ``i`` as corresponding to the ``x`` axis,
``j`` to the ``y`` axis, and ``k`` to the ``z`` axis — evidently unlike Hamilton himself.
Moreover, Tait even says

> Suppose ``i`` to be drawn eastwards, ``j`` northwards, and ``k`` upwards.

This is consistent with our standard right-handed system, with the ``k`` axis pointing out
of the page.



[Tait_1868](@cite) discusses the decomposition of a general rotation into three rotations
about successive axes, which we now call "Euler angles".  Tait uses the symbols ``(ψ, θ,
ϕ)`` for the three angles.  He actually does so just to express that set of rotations in
terms of a single quaternion.

> Here the vectors ``i``, ``j``, ``k`` in the original position of the body correspond to
``\overline{OA}``, ``\overline{OB}``, ``\overline{OC}`` respectively, at time ``t``.  The
transposition is effected by — *first*, a rotation ``ψ`` about ``k``; *second*, a
rotation ``θ`` about the new position of the line originally coinciding with ``j``;
*third*, a rotation ``ϕ`` about the final position of the line at first coinciding with
``k'``.

So this is what would probably now be called the ``z-y'-z''`` convention for Euler angles
``(ψ, θ, ϕ)``: each rotation is about an axis of the body, as carried along by the
preceding rotations.  [Our own convention](@ref summary_euler_angles) is stated the other
way round, as rotations about *fixed* axes — ``γ`` about ``𝐳``, then ``β`` about ``𝐲``,
then ``α`` about ``𝐳`` — but the two descriptions give the same quaternion product
``e^{α𝐤/2}\, e^{β𝐣/2}\, e^{γ𝐤/2}``.  Tait's angles are therefore ours *in the same
order*,
```math
(α, β, γ) = (ψ, θ, ϕ),
```
and it is only when his rotations are re-expressed about the fixed axes that the list
reverses, to ``ϕ``, ``θ``, ``ψ``.

Indeed, Tait goes on to derive (somewhat laboriously) the expression for the quaternion:

```math
q = \cos \frac{ϕ + ψ}{2} \cos \frac{θ}{2}
  + i \sin \frac{ϕ - ψ}{2} \sin \frac{θ}{2}
  + j \cos \frac{ϕ - ψ}{2} \sin \frac{θ}{2}
  + k \sin \frac{ϕ + ψ}{2} \cos \frac{θ}{2},
```

which is exactly the same as our expression from `from_euler_angles(ψ, θ, ϕ)`.

"""
