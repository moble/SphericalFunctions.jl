# Outline

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
      is trivial on ``-1``.  Noting that ``\exp(\pi \vec{v}) = -1``
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
Furthermore, on page 354, they point out that ``\{(1/\sqrt{2\pi})
e^{imϕ}\}`` is an orthonormal basis of ``L^2(0,2\pi)``, while the
set ``\{1/c_{n,m} P_n^m(\cos θ)`` is an orthonormal basis of
``L^2(0, \pi)`` in the ``θ`` coordinate.  Therefore, the product
of these two sets is an orthonormal basis of the product space
``L^2\left((0,2\pi) \times (0, \pi)\right)``, which forms a coordinate
space for ``𝕊²``.  I would probably modify this to point out that
``(0,2\pi)`` is really ``𝕊¹``, and then we could extend it to point
out that you can throw on another factor of ``𝕊¹`` to cover ``𝕊³``,
which happens to give us the Wigner D-matrices.

We first define the rotor that takes ``(\hat{x}, \hat{y}, \hat{z})``
onto ``(\hat{θ}, \hat{ϕ}, \hat{r})``.  Then, we can invert
that, so that given a rotor that specifies such a rotation exactly, we
can get the spherical coordinates — or specifically ``\sin θ``,
``\cos θ``, and ``\exp(iϕ)``.

Then, with the universally agreed-upon ``Y`` as given in terms of
spherical coordinates, we can rewrite it directly to work with
quaternion components, and then it immediately applies to general
rotations, which allows us to figure out where the ``s`` should go.
That is, we can essentially derive ``{}_sY`` from the universal
formula for ``Y``.

Then, we can simply follow Wigner around Eq. (15.21) to derived a
transformation law in the form
```math
{}_sY_{ℓ,m'}(R_{θ', ϕ'}) = \sum_m M_{m',m}(R)
{}_sY_{ℓ,m}(R_{θ, ϕ}),
```
for some matrix ``M``.  Note that I have written this as if the
``{}_sY`` functions are column vectors.  The reason this happens is
because I want to write ``R_{θ', ϕ'} = R\, R_{θ, ϕ}``,
rather than swapping the order of the rotations on the right-hand
side.

The big problem here is that Wigner, in his Eq. (15.21) defines the
transformation matrix as if the eigenfunctions formed a row vector
instead of a column vector, which means that his matrix is transposed
compared to what I want to write.  I suppose maybe other authors then
just consider the inverse rotation, so that they can work with the
conjugate transpose, which is why we see the relative conjugate.

* Since ``Y`` is universal, let's start with that as non-negotiable,
  and see if we can derive the relationship to ``𝔇``.
* ``R_{θ, ϕ}`` is a unit quaternion that rotates the point
  described by Cartesian coordinates (0,0,1) onto the point described
  by spherical coordinates ``(θ, ϕ)``.
* Just textually, it makes the most sense to write
  ```math
  R_{θ', ϕ'} = R\, R_{θ, ϕ}
  ```
  for some rotation ``R``.  Now, we just need to interpret ``R``.
* Again, just textually, it makes the most sense to write
  ```math
  Y_{ℓ,m'}(θ', ϕ') = \sum_m 𝔇^{(ℓ)}_{m',m}(R)
  Y_{ℓ,m}(θ, ϕ),
  ```
  or, generalizing to spin-weighted spherical harmonics
  ```math
  {}_{s}Y_{ℓ,m'}(R_{θ', ϕ'}) = \sum_m 𝔇^{(ℓ)}_{m',m}(R)
  {}_{s}Y_{ℓ,m}(R_{θ, ϕ}).
  ```
* We also have that ``𝔇`` obeys the representation
  property, so
  ```math
  𝔇^{(ℓ)}_{m',m''}(R_{θ', ϕ'})
  = \sum_{m} 𝔇^{(ℓ)}_{m',m}(R)
  𝔇^{(ℓ)}_{m,m''}(R_{θ, ϕ}).
  ```
  - There is no reason that I can see to introduce a conjugation
  - The fact that ``m''`` appears on both sides of the equation means
    that it must correspond to ``s`` — though we have to check the
    behavior under final rotation to determine the sign.

```math
{}_{s}Y_{ℓ,m}(R_{θ, ϕ})
\propto
𝔇^{(ℓ)}_{m,\propto s}(R_{θ, ϕ})
```


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


## Spherical harmonics

Fortunately, there does not seem to be any disagreement in the physics
literature about the definition of the spherical harmonics; everyone
uses the Condon-Shortley convention.  Or at least, they say they do.
The problem arises when people define the spherical harmonics in terms
of the Legendre polynomials, for which there is a sign ambiguity.
Therefore, to ensure that we are using the same conventions, we need
to go back to the original definition of the spherical harmonics by
Condon and Shortley.

### Condon-Shortley phase

The [Condon-Shortley](@cite CondonShortley_1935) phase convention is a
choice of phase factors in the definition of the spherical harmonics
that requires the coefficients in
```math
L_{\pm} |ℓ,m\rangle = α^{\pm}_{ℓ,m} |ℓ, m \pm 1\rangle
```
to be real and positive.  The reasoning behind this choice is
explained more fully in Section 2 of [Ufford and Shortley
(1932)](@cite UffordShortley_1932).  As a more practical matter, the
Condon-Shortley phase describes signs chosen in the expression for
spherical harmonics.  The key expression is Eq. (15) of section 4³
(page 52) of [Condon-Shortley](@cite CondonShortley_1935):
```math
\Theta(ℓ, m) = (-1)^ℓ \sqrt{\frac{2ℓ+1}{2} \frac{(ℓ+m)!}{(ℓ-m)!}}
\frac{1}{2^ℓ ℓ!} \frac{1}{\sin^mθ}
\frac{d^{ℓ-m}}{d(\cos θ)^{ℓ-m}} \sin^{2ℓ}θ.
```
When multiplied by Eq. (5) ``\Phi(m) = e^{imϕ} / \sqrt{2\pi}``,
this gives the spherical harmonic function.  The right-hand side of
the expression above is usually immediately replaced by a simpler
expression using Legendre polynomials, but this just shifts sign
ambiguity into the definition of the Legendre polynomials.  Instead,
we can expand the above expression directly for the first few ``ℓ``
values and/or use automatic differentiation to actually test their
original expression as such against the function implemented in this
package.  The first few values are given in a footnote to Condon and
Shortley's Eq. (15) (and have been verified separately by hand and by
computation with SymPy):
```math
\begin{aligned}
\Theta(0,0) &= \sqrt{\frac{1}{2}} \\
\Theta(1,0) &= \sqrt{\frac{3}{2}} \cos θ &
\Theta(1,\pm1) &= \mp \sqrt{\frac{3}{4}} \sin θ \\
\Theta(2,0) &= \sqrt{\frac{5}{8}} (2\cos^2θ - \sin^2θ) &
\Theta(2,\pm1) &= \mp \sqrt{\frac{15}{4}} \cos θ \sin θ &
\Theta(2,\pm2) &= \sqrt{\frac{15}{16}} \sin^2θ \\
\Theta(3,0) &= \sqrt{\frac{7}{8}} (2\cos^3θ - 3\cos θ\sin^2θ) &
\Theta(3,\pm1) &= \mp \sqrt{\frac{21}{32}} (4\cos^2θ\sin θ - \sin^3θ) &
\Theta(3,\pm2) &= \sqrt{\frac{105}{16}} \cos θ \sin^2θ &
\Theta(3,\pm3) &= \mp \sqrt{\frac{35}{32}} \sin^3θ
\end{aligned}
```
These are tested, along with the results from automatic
differentiation, every time this package is updated.  The result is
perfect agreement, so that we can definitively say that ***the
spherical-harmonic functions provided by this package obey the
Condon-Shortley phase convention.***

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
L_x &= -i\hbar \left( \sin ϕ \frac{\partial}{\partial θ} + \cot θ \cos ϕ \frac{\partial}{\partial ϕ} \right), \\
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


## Wigner ``𝔇`` and ``d`` matrices

Wigner's Eqs. (11.18) and (11.19) define the real orthogonal
transformation ``𝐑`` by
```math
x'_i = R_{ij} x_j
```
and the operator ``𝐏_{𝐑}`` to act on a function
``f`` such that
```math
𝐏_{𝐑} f(x'_1, \ldots) = f(x_1, \ldots).
```
Then, his Eq. (15.5) presumably implies
```math
Y_{ℓ,m}(ϑ', φ')
= 𝐏_{\{α, β, γ\}} Y_{ℓ,m}(ϑ, φ)
= \sum_{m'} 𝔇^{(ℓ)}(\{α, β, γ\})_{m',m}
  Y_{ℓ,m'}(ϑ, φ),
```
where ``\{α, β, γ\}`` takes ``(ϑ, φ)`` to
``(ϑ', φ')``.  In any case, we can now leave behind this
``𝐏`` notation and just look at the beginning and end of the
equation above as the critical relationship in Wigner's notation.


Eq. (44b) of [Boyle (2016)](@cite Boyle_2016) says
```math
L_{\pm} 𝔇^{(ℓ)}_{m',m}(𝐑)
= \sqrt{(ℓ \mp m')(ℓ \pm m' + 1)} 𝔇^{(ℓ)}_{m' \pm 1, m}(𝐑).
```
while Eq. (21) relates the Wigner D-matrix to the spin-weighted spherical harmonics as
```math
{}_{s}Y_{ℓ,m}(𝐑)
= (-1)^s \sqrt{\frac{2ℓ+1}{4\pi}} 𝔇^{(ℓ)}_{m,-s}(𝐑).
```
Plugging the latter into the former, we get
```math
L_{\pm} {}_{s}Y_{ℓ,m}(𝐑)
= \sqrt{(ℓ \mp m)(ℓ \pm m + 1)} {}_{s}Y_{ℓ,m \pm 1}(𝐑).
```
That is, in our conventions we have
```math
α^{\pm}_{ℓ,m} = \sqrt{(ℓ \mp m)(ℓ \pm m + 1)},
```
which is always real and positive, and thus consistent with the Condon-Shortley phase
convention.


### Properties

* ``D^j_{m'm}(α,β,γ) = (-1)^{m'-m} D^j_{-m',-m}(α,β,γ)^*``
* ``(-1)^{m'-m}D^{j}_{mm'}(α,β,γ)=D^{j}_{m'm}(γ,β,α)``
* ``d_{m',m}^{j}=(-1)^{m-m'}d_{m,m'}^{j}=d_{-m,-m'}^{j}``

```math
\begin{aligned}
d_{m',m}^{j}(\pi)        &= (-1)^{j-m}  δ_{m',-m} \\[6pt]
d_{m',m}^{j}(\pi-β)  &= (-1)^{j+m'}  d_{m',-m}^{j}(β)\\[6pt]
d_{m',m}^{j}(\pi+β)  &= (-1)^{j-m}  d_{m',-m}^{j}(β)\\[6pt]
d_{m',m}^{j}(2\pi+β) &= (-1)^{2j}    d_{m',m}^{j}(β)\\[6pt]
d_{m',m}^{j}(-β)     &= d_{m,m'}^{j}(β) = (-1)^{m'-m} d_{m',m}^{j}(β)
\end{aligned}
```







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
