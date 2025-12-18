md"""
# Hamilton (1844 and 1853)

[Hamilton_1844](@cite) is the foundational work on quaternions, where William Rowan Hamilton
first introduced them in 1844.  On page 492, he introduces the quaternion, denoted as ``a +
i b + j c + k d``, or ``(a, b, c, d)``.  He quickly establishes the multiplication rules for
the basis elements ``(i, j, k)`` as
```math
\begin{gathered}
i^2 = j^2 = k^2 = -1;
\quad
ij = -ji = k;
\quad
jk = -kj = i;
\quad
ki = -ik = j.
\end{gathered}
```
He doesn't include the common ``ijk = -1`` explicitly, but it follows from the above.

He also introduces a parameterization:
```math
a = μ \cos ϱ,
b = μ \sin ϱ \cos ϕ,
c = μ \sin ϱ \sin ϕ \cos ψ,
d = μ \sin ϱ \sin ϕ \sin ψ.
```
He then says (emphasis in original)

> I would call ``ϱ`` the *amplitude* of the quaternion; ``ϕ`` its *colatitude*; and ``ψ``
> its *longitude*.  The *modulus* is ``μ``.

(Note that Hamilton is using the archaic meaning of the word "amplitude", which was the
standard term from complex analysis in the early 1800s, but which we would now call the
"argument" or "angle"; he uses "modulus" where physicists and engineers would now frequently
use "amplitude" or "magnitude".)  This is similar to the axis-angle representation of
rotations, where ``ϱ`` is the (whole, not half) angle of rotation, and ``(ϕ, ψ)`` are
spherical coordinates for the axis of rotation — except that this would imply that ``i``
corresponds to the ``z`` axis if we take the name "colatitude" seriously.  And the angle
``ψ`` is measured from the ``j`` axis, so we would probably identify ``j`` with the ``x``
axis.  So ``i, j, k`` correspond to the axes ``z, x, y``.

There is no mention of "Euler angles" in Hamilton's introduction of quaternions.


[Hamilton_1853](@cite) is Hamilton's later monograph on quaternions, published in 1853.  In
the table of contents — which is really a collection of summaries of each section — we find
in Lecture II §x this statement (emphasis in original):

> the  symbols ``i, j, k`` come to denote here *three rectangular vector-units* (supposed
> usually, in these Lectures, to be in the directions of *south*, *west*, and *up*)

This is a *left*-handed system.  One sometimes sees claims that the quaternions are somehow
"inherently" left-handed, but this is simply not true; Hamilton could have just as easily
chosen a right-handed system, and later authors generally did.


"""
