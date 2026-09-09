# Outline and working notes

!!! warning "Scratch page"
    This page is not part of the rendered navigation.  It holds the
    original outline for the conventions pages and assorted working
    notes that were *not* promoted to the [Summary](@ref "Summary") or
    [Details](@ref "Details") pages.  Anything here that contradicts
    those pages is stale; the settled conventions live there.

* Three-dimensional Euclidean space
  - Cartesian coordinates ``(x, y, z)`` => ℝ³
  - Cartesian basis vectors ``(𝐱, 𝐲, 𝐳,)``
  - Euclidean norm => Euclidean metric
  - Spherical coordinates
    - Specifically give transformation to/from ``(x, y, z)``
    - Derive metric in these coordinates from transformation
  - Integration / measure on two-sphere
    - Derive as restriction of full metric, in both coordinate systems
* Four-dimensional Euclidean space
  - Eight-dimensional Clifford algebra over the tangent *vector space* ``Tℝ³``
  - Four-dimensional even sub-algebra => ℝ⁴
  - Coordinates ``(W, X, Y, Z)``
  - Basis vectors ``(𝟏, 𝐢, 𝐣, 𝐤)``, but we usually just omit ``𝟏``
    - Show a few essential formulas establishing the product and its conventions
  - Unit quaternions are isomorphic to ``\mathbf{Spin}(3) =
    \mathbf{SU}(2)``; double covers ``\mathbf{SO}(3)``
    - Be explicit about the mapping between vector in ℝ³ and quaternions
    - Show how a unit quaternion can be used to rotate a vector
  - Spherical coordinates (hyperspherical / Euler)
    - Specifically give transformation to/from ``(W, X, Y, Z)``
    - Derive metric in these coordinates from transformation
    - Express unit quaternion in Euler angles
  - Integration / measure / Haar measure on three-sphere
    - Derive as restriction of full metric, in both coordinate systems
* Angular momentum operators / functional analysis
  - Express angular momentum operators in terms of quaternion components
  - Express angular momentum operators in terms of Euler angles
  - Show for both the three- and two-spheres
  - Show how they act on functions on the three-sphere
* Representation theory / harmonic analysis
  - Representations show up in Fourier analysis on groups
  - Peter-Weyl theorem
    - Generalizes Fourier analysis to compact groups
    - A basis of functions on the group is given by matrix elements of
      group representations
  - Representation theory of ``\mathbf{Spin}(3)``
    - Show how the Lie algebra is represented by the angular-momentum operators
    - Show how the Lie group is represented by the Wigner D-matrices
    - Demonstrate that ``𝔇`` is a representation
    - Demonstrate its behavior under left and right rotation
    - Demonstrate orthonormality
  - Representation theory of ``\mathbf{SO}(3)``
    - There are several places in [Folland](@cite Folland_2016) (e.g.,
      above corollary 5.48) where he mentions that representations of
      a quotient group are just representations that are trivial
      (evidently meaning mapping everything to the identity matrix) on
      the factor.  I can't find anywhere that he explains this
      explicitly, but it seems easy enough to show.  He might do it
      using characters.
    - For ``\mathbf{Spin}(3)`` and ``\mathbf{SO}(3)``, the factor
      group is just ``\{1, -1\}``.  Presumably, every representation
      acting on ``1`` will give the identity matrix, so that's
      trivial.  So we just need a criterion for when a representation
      is trivial on ``-1``.  Noting that ``\exp(π \vec{v}) = -1``
      for any ``\vec{v}``, I think we can show that this requires
      ``m \in \mathbb{Z}``.
    - Basically, the point is that the representations of
      ``\mathbf{SO}(3)`` are just the integer representations of
      ``\mathbf{Spin}(3)``.
  - Restrict to homogeneous space (S³ -> S²)
    - The circle group is a closed (normal?) subgroup of
      ``\mathbf{Spin}(3)``, which we might implement as initial
      multiplication about a particular axis.
    - In Eq. (2.47) [Folland (2016)](@cite Folland_2016) defines a
      functional taking a function on the group to a function on the
      homogeneous space by integrating over the factor (the circle
      group).  This gives you the spherical harmonics, but *not* the
      spin-weighted spherical harmonics — because the spin-weighted
      spherical harmonics cannot be defined on the 2-sphere.
    - Spin weight comes from Fourier analysis on the subgroup.
    - Representation matrices transfer to the homogeneous space, with
      sparsity patterns

 

# Notes

Spherical harmonics as functions on homogeneous space.
https://www.youtube.com/watch?v=TnFvOa9v7do gives some nice
discussion; maybe the paper has better references.

Theorem 2.16 of [Hanson-Yakovlev](@cite HansonYakovlev_2002) says that
an orthonormal basis of a product of ``L^2`` spaces is given by the
product of the orthonormal bases of the individual spaces.
Furthermore, on page 354, they point out that ``\{(1/\sqrt{2π})
e^{imϕ}\}`` is an orthonormal basis of ``L^2(0,2π)``, while the
set ``\{1/c_{n,m} P_n^m(\cos θ)\}`` is an orthonormal basis of
``L^2(0, π)`` in the ``θ`` coordinate.  Therefore, the product
of these two sets is an orthonormal basis of the product space
``L^2\left((0,2π) \times (0, π)\right)``, which forms a coordinate
space for ``𝕊²``.  I would probably modify this to point out that
``(0,2π)`` is really ``𝕊¹``, and then we could extend it to point
out that you can throw on another factor of ``𝕊¹`` to cover ``𝕊³``,
which happens to give us the Wigner D-matrices.

The derivation of the transformation law and of the relation between
``{}_sY`` and ``𝔇`` that used to be sketched here is now carried out
in full on the [Details](@ref conv_wigner_D) page; see in particular
the [rotation law](@ref conv_rotation_law) and the [definition of the
spin-weighted spherical harmonics](@ref conv_swsh).

## collapsible markdown?

```@raw html
<details><summary>CLICK ME</summary>
```
#### yes, even hidden code blocks!

```julia
println("hello world!")
```
```@raw html
</details>
```


# More notes

## Angular-momentum operators

* First, a couple points about ``-i\hbar``:
  - The finite transformations look like ``\exp[-i θ L_j]``, but
    the factor of ``i`` introduced here just cancels the one in the
    ``L_j``, and the sign is just chosen to make the result consistent
    with our notion of active or passive transformations.
  - Any factors of ``\hbar`` are included *purely* for the sake of
     convenience.
  - The factor ``i`` comes from plain functional analysis: We need a
    self-adjoint operator, and ``\partial_x`` by itself is
    anti-self-adjoint (as can be verified by evaluating on ``\langle
    x' | x \rangle = δ(x-x')``, which switches sign based on
    which is being differentiated).  We want self-adjoint operators so
    that we get purely real eigenvalues.  [Van Neerven](@cite
    vanNeerven_2022) cites this in a more rigorous context in his
    Example (10.40) (page 331), with more explanation around Eq.
    (15.17) (page 592).  The "self-adjoint ``\iff`` real eigenvalues"
    condition is item (1) in his Corollary 9.18.

Wigner's ``𝔇`` matrices are defined as matrix elements of a rotation in
the basis of spherical harmonics.  That rotation is defined in terms
of the generators of rotation, which are expressed in terms of the
angular-momentum operators.  Therefore, to really understand
conventions for the ``𝔇`` matrices, we need to understand conventions
for the angular-momentum operators.

There is universal agreement that the angular momentum is defined as
``𝐋 = 𝐱 \times 𝐩``, where ``𝐱`` is
the position vector and ``𝐩`` is the momentum vector.  In
quantum mechanics, there is further agreement that the momentum
operator becomes ``-i\hbar\nabla``.  Thus, in operator form, the
angular momentum can be decomposed as
```math
\begin{aligned}
L_x &= -i\hbar \left( y \frac{\partial}{\partial z} - z \frac{\partial}{\partial y} \right), \\
L_y &= -i\hbar \left( z \frac{\partial}{\partial x} - x \frac{\partial}{\partial z} \right), \\
L_z &= -i\hbar \left( x \frac{\partial}{\partial y} - y \frac{\partial}{\partial x} \right).
\end{aligned}
```
We can transform these to use spherical coordinates and obtain
```math
\begin{aligned}
L_x &= i\hbar \left( \sin ϕ \frac{\partial}{\partial θ} + \cot θ \cos ϕ \frac{\partial}{\partial ϕ} \right), \\
L_y &= -i\hbar \left( \cos ϕ \frac{\partial}{\partial θ} - \cot θ \sin ϕ \frac{\partial}{\partial ϕ} \right), \\
L_z &= -i\hbar \frac{\partial}{\partial ϕ}.
\end{aligned}
```
The conventions we choose *must* be chosen to agree with these —
modulo factors of ``\hbar``, which are nonstandard in mathematics.  We
will have to check this, and the Condon-Shortley requirement that when
applied to spherical harmonics they produce real and positive
coefficients.

I defined these in Eqs. (42) and (43) of [Boyle (2016)](@cite Boyle_2016) as
```math
\begin{aligned}
L_{j} f(𝐑) &\colonequals -z \left. \frac{\partial}{\partial θ}
f\left(e^{θ 𝐞_j / 2} 𝐑 \right) \right|_{θ=0}, \\
K_{j} f(𝐑) &\colonequals -z \left. \frac{\partial}{\partial θ}
f\left(𝐑 e^{θ 𝐞_j / 2}\right) \right|_{θ=0},
\end{aligned}
```
where ``𝐞_j`` is the unit vector in the ``j`` direction.
Surprisingly, I found that [Edmonds](@cite Edmonds_2016) expresses
essentially the same thing in the equations following his Eq. (4.1.5).

Condon and Shortley's Eq. (1) of section 4³ (page 50) defines
```math
L_z = -i \hbar \frac{\partial}{\partial ϕ},
```
while Eq. (8) on the following page defines
```math
\begin{aligned}
L_x + i L_y &= \hbar e^{iϕ} \left( \frac{\partial}{\partial θ} + i \cot θ \frac{\partial}{\partial ϕ} \right), \\
L_x - i L_y &= \hbar e^{-iϕ} \left(-\frac{\partial}{\partial θ} + i \cot θ \frac{\partial}{\partial ϕ} \right).
\end{aligned}
```
Note that one is not the conjugate of the other!  This is because of
the factors of ``-i`` in the definitions of ``L_x`` and ``L_y``.

[Edmonds](@cite Edmonds_2016) gives the *total* angular-momentum
operator for a rigid body in Eq. (2.2.2) as
```math
\begin{aligned}
L_x &= -i\hbar \left(-\cos α \cot β \frac{\partial}{\partial α} - \sin α \frac{\partial}{\partial β} + \frac{\cos α}{\sin β} \frac{\partial}{\partial γ} \right), \\
L_y &= -i\hbar \left(-\sin α \cot β \frac{\partial}{\partial α} + \cos α \frac{\partial}{\partial β} + \frac{\sin α}{\sin β} \frac{\partial}{\partial γ} \right), \\
L_z &= -i\hbar \frac{\partial}{\partial α}.
\end{aligned}
```

!!! note
    The definitions from Boyle (2016) quoted just above use the
    *opposite* sign for the right operator (there called ``K``) from
    the one settled on the [Details](@ref conv_L_R_definitions) page,
    where ``R_𝐮 f(𝐑) = -i\, d/dϵ\, f(𝐑 e^{-ϵ𝐮/2})``.  The
    Condon–Shortley check of the ladder coefficients that used to
    follow here is now derived there as well.

## Rotor scraps

```math
\begin{gather}
R = \cos ϵ + \sin ϵ\, \hat{𝔯} \\
R𝐯 = \cos ϵ 𝐯 + \sin ϵ\, \hat{𝔯}𝐯 \\
R𝐯R^{-1} = (𝐯\cos ϵ + \sin ϵ\, \hat{𝔯}𝐯)(\cos ϵ - \sin ϵ\, \hat{𝔯}) \\
R𝐯R^{-1} = 𝐯\cos^2ϵ + \sin^2ϵ\, \hat{𝔯}𝐯\hat{𝔯}^{-1} + \sin ϵ \cos ϵ\, (\hat{𝔯}𝐯 - 𝐯\hat{𝔯}) \\
R𝐯R^{-1} = \begin{cases}
𝐯 & 𝐯 \hat{𝔯} = \hat{𝔯}𝐯 \\
𝐯(\cos^2ϵ - \sin^2ϵ) + 2 \sin ϵ \cos ϵ\, \frac{[\hat{𝔯}, 𝐯]}{2} & 𝐯 \hat{𝔯} = -\hat{𝔯}𝐯 \\
\end{cases} \\
R𝐯R^{-1} = \begin{cases}
𝐯 & 𝐯 \hat{𝔯} = \hat{𝔯}𝐯 \\
\cos2ϵ 𝐯 + \sin2ϵ \frac{[\hat{𝔯}, 𝐯]}{2} & 𝐯 \hat{𝔯} = -\hat{𝔯}𝐯 \\
\end{cases} \\
\end{gather}
```




Using techniques from geometric algebra, we can easily prove that the
result is another vector, so we can measure its (squared) norm just by
multiplying it by itself:
```math
\begin{aligned}
\| 𝐑\, 𝐯\, 𝐑^{-1} \|^2
&= 𝐑\, 𝐯\, 𝐑^{-1}\, 𝐑\, 𝐯\, 𝐑^{-1} \\
&= 𝐑\, 𝐯\, 𝐯\, 𝐑^{-1} \\
&= \|𝐯\|^2\, 𝐑\, 𝐑^{-1} \\
&= \|𝐯\|^2
\end{aligned}
```
That is, ``𝐯' = 𝐑\, 𝐯\, 𝐑^{-1}`` has the same norm as ``𝐯``,
which means that ``𝐯'`` is a rotation of ``𝐯``.  Given the constraint
on the norm of ``𝐑``, we can rewrite it as

## Representation theory / harmonic analysis (moved from Details)

  - Representations show up in Fourier analysis on groups
  - Peter-Weyl theorem
    - Generalizes Fourier analysis to compact groups
    - Has three parts, [as given by Wikipedia](https://en.wikipedia.org/wiki/Peter%E2%80%93Weyl_theorem):
      1. "The matrix coefficients of irreducible representations of
         ``G`` are dense in the space ``C(G)`` of continuous
         complex-valued functions on ``G``, and thus also in the space
         ``L^2(G)`` of square-integrable functions."
      2. Unitary representations of ``G`` are completely reducible.
      3. "The regular representation of ``G`` on ``L^2(G)`` decomposes
         as the direct sum of all irreducible unitary representations.
         Moreover, the matrix coefficients of the irreducible unitary
         representations form an orthonormal basis of ``L^2(G)``."
  - Representation theory of ``\mathbf{Spin}(3)``
    - Show how the Lie algebra is represented by the angular-momentum operators
    - Show how the Lie group is represented by the Wigner D-matrices
    - Demonstrate that ``𝔇`` is a representation
    - Demonstrate its behavior under left and right rotation
    - Demonstrate orthonormality
  - Representation theory of ``\mathbf{SO}(3)``
    - There are several places in [Folland](@cite Folland_2016) (e.g.,
      above corollary 5.48) where he mentions that representations of
      a quotient group are just representations that are trivial
      (evidently meaning mapping everything to the identity matrix) on
      the factor.  I can't find anywhere that he explains this
      explicitly, but it seems easy enough to show.  He might do it
      using characters.
    - For ``\mathbf{Spin}(3)`` and ``\mathbf{SO}(3)``, the factor
      group is just ``\{1, -1\}``.  Presumably, every representation
      acting on ``1`` will give the identity matrix, so that's
      trivial.  So we just need a criterion for when a representation
      is trivial on ``-1``.  Noting that ``\exp(π \vec{v}) = -1``
      for any ``\vec{v}``, I think we can show that this requires
      ``m \in \mathbb{Z}``.
    - Basically, the point is that the representations of
      ``\mathbf{SO}(3)`` are just the integer representations of
      ``\mathbf{Spin}(3)``.
  - Restrict to homogeneous space (S³ -> S²)
    - The circle group is a closed (normal?) subgroup of
      ``\mathbf{Spin}(3)``, which we might implement as initial
      multiplication about a particular axis.
    - In Eq. (2.47) [Folland (2016)](@cite Folland_2016) defines a
      functional taking a function on the group to a function on the
      homogeneous space by integrating over the factor (the circle
      group).  This gives you the spherical harmonics, but *not* the
      spin-weighted spherical harmonics — because the spin-weighted
      spherical harmonics cannot be defined on the 2-sphere.
    - Spin weight comes from Fourier analysis on the subgroup.
    - Representation matrices transfer to the homogeneous space, with
      sparsity patterns

Theorem 2.16 of [Hanson-Yakovlev](@cite HansonYakovlev_2002) says that
an orthonormal basis of a product of ``L^2`` spaces is given by the
product of the orthonormal bases of the individual spaces.
Furthermore, on page 354, they point out that ``\{(1/\sqrt{2π})
e^{imϕ}\}`` is an orthonormal basis of ``L^2(0,2π)``, while the
set ``\{1/c_{n,m} P_n^m(\cos θ)\}`` is an orthonormal basis of
``L^2(0, π)`` in the ``θ`` coordinate.  Therefore, the product
of these two sets is an orthonormal basis of the product space
``L^2\left((0,2π) \times (0, π)\right)``, which forms a coordinate
space for ``𝕊²``.  I would probably modify this to point out that
``(0,2π)`` is really ``𝕊¹``, and then we could extend it to point
out that you can throw on another factor of ``𝕊¹`` to cover ``𝕊³``,
which happens to give us the Wigner D-matrices.

## Recursion relations (moved from Details)

[Gumerov and Duraiswami (2001)](@cite Gumerov_2001) derive their
recursion relations by differentiating solutions of the Helmholtz
equation ``\nabla^2 ψ + k^2 ψ = 0`` as ``\tfrac{1}{k} \nabla
ψ``.  More precisely, they differentiate both sides of the equation
relating one solution to its rotated form — which naturally involves
Wigner's ``𝔇`` matrix.  Using orthogonal basis functions
for the solution, this allows them to equate terms on the two sides
proportional to a given basis function, which leaves them with
expressions involving sums of only the ``𝔇`` matrices and
some coefficients depending on the indices of the basis functions (and
hence of ``𝔇``) on both sides of the equation.  Since
``\nabla`` is a 3-vector operator, this gives them three relations.

This, of course, is happening in 3-D space, since ``ψ`` is a
function of location in the Helmholtz equation.  It seems likely to
me, however, that we could use the 4-D (quaternionic) version of the
functions.  Note that G&D use ``\partial_z`` and ``\partial_x \pm i
\partial_y`` as their operators to differentiate the functions — that
is, the derivatives are with respect to Cartesian coordinates, which
may be more similar to the right-derivative defined above.  However, I
don't know that we'll necessarily be able to achieve the same results
with just angular-momentum operators, since their operators do involve
moving off of the sphere.  Maybe we'd need to move off of the sphere
in 4-D space to get comparable results.  Or maybe just use something
like ``𝐫 ∧ L``, which should also have 3 degrees of freedom.

The SWSHs/``𝔇`` functions can be naturally promoted to
functions not just on the 3-sphere, but also in 4-D space just by
allowing the quaternions to be non-unit quaternions.
