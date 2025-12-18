md"""
# Clifford (1878)

[Clifford_1878](@citet) introduced Clifford algebras in 1878.  The paper starts off with an
introduction to quaternions — though interestingly from a "projective" perspective.[^1] This
first section is a bit muddled from a modern perspective, but Clifford is using it to
develop his ideas.  Eventually, he gets to a definition that says

> the symbols ``ι₁, ι₂, ι₃`` may be taken to mean unit vectors along the axes.

[^1]:

    That is, in addition to what we would recognize as basis vectors ``ι₁, ι₂, ι₃``, he
    includes an element ``ι₀`` that is "at an infinite distance" from the others.  This
    appears to be because he recognizes the importance of what we now call the even
    subalgebra of the Clifford algebra, and wants to be able to include the pseudoscalar in
    that subalgebra — which cannot happen with a three-dimensional vector space alone.

But he actually settles on *negative* squares for these basis elements:

> We are therefore obliged to write ``ι₂² = -1``, and in a similar way we may find ``ι₁² =
> ι₃² = -1``.

This differs from our modern treatments that usually assume that the basis vectors square to
``+1``.  To compensate for this, he goes on to identify Hamilton's quaternions with products
of pairs of ``ι₁, ι₂, ι₃`` that happen to be exactly the *negative* of modern
definitions — or at least the ones we use here:

```math
i = ι₂ ι₃,
\quad
j = ι₃ ι₁,
\quad
k = ι₁ ι₂.
```

Again, only in later sections does the exposition become clearer and more modern.  There, he
returns to quaternions briefly, and manages to describe them in exactly the way we would
today, in terms of the even subalgebra of the Clifford algebra of three-dimensional space.
(He even notes that they are isomorphic to the *full* Clifford algebra of two-dimensional
space.)

Clifford — like Hamilton and Tait — failed to recognize the vital importance of conjugation
(or "sandwiching") when using quaternions to represent rotations; presumably by analogy with
complex numbers, he only considered left-multiplication by a quaternion as the operation 

"""
