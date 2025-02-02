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
    - Demonstrate that ``\mathfrak{D}`` is a representation
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
e^{im\phi}\}`` is an orthonormal basis of ``L^2(0,2\pi)``, while the
set ``\{1/c_{n,m} P_n^m(\cos\theta)`` is an orthonormal basis of
``L^2(0, \pi)`` in the ``\theta`` coordinate.  Therefore, the product
of these two sets is an orthonormal basis of the product space
``L^2\left((0,2\pi) \times (0, \pi)\right)``, which forms a coordinate
space for ``S^2``.  I would probably modify this to point out that
``(0,2\pi)`` is really ``S^1``, and then we could extend it to point
out that you can throw on another factor of ``S^1`` to cover ``S^3``,
which happens to give us the Wigner D-matrices.

We first define the rotor that takes ``(\hat{x}, \hat{y}, \hat{z})``
onto ``(\hat{\theta}, \hat{\phi}, \hat{r})``.  Then, we can invert
that, so that given a rotor that specifies such a rotation exactly, we
can get the spherical coordinates — or specifically ``\sin\theta``,
``\cos\theta``, and ``\exp(i\phi)``.

Then, with the universally agreed-upon ``Y`` as given in terms of
spherical coordinates, we can rewrite it directly to work with
quaternion components, and then it immediately applies to general
rotations, which allows us to figure out where the ``s`` should go.
That is, we can essentially derive ``{}_sY`` from the universal
formula for ``Y``.

Then, we can simply follow Wigner around Eq. (15.21) to derived a
transformation law in the form
```math
{}_sY_{\ell,m'}(R_{\theta', \phi'}) = \sum_m M_{m',m}(R)
{}_sY_{\ell,m}(R_{\theta, \phi}),
```
for some matrix ``M``.  Note that I have written this as if the
``{}_sY`` functions are column vectors.  The reason this happens is
because I want to write ``R_{\theta', \phi'} = R\, R_{\theta, \phi}``,
rather than swapping the order of the rotations on the right-hand
side.

The big problem here is that Wigner, in his Eq. (15.21) defines the
transformation matrix as if the eigenfunctions formed a row vector
instead of a column vector, which means that his matrix is transposed
compared to what I want to write.  I suppose maybe other authors then
just consider the inverse rotation, so that they can work with the
conjugate transpose, which is why we see the relative conjugate.

* Since ``Y`` is universal, let's start with that as non-negotiable,
  and see if we can derive the relationship to ``\mathfrak{D}``.
* ``R_{\theta, \phi}`` is a unit quaternion that rotates the point
  described by Cartesian coordinates (0,0,1) onto the point described
  by spherical coordinates ``(\theta, \phi)``.
* Just textually, it makes the most sense to write
  ```math
  R_{\theta', \phi'} = R\, R_{\theta, \phi}
  ```
  for some rotation ``R``.  Now, we just need to interpret ``R``.
* Again, just textually, it makes the most sense to write
  ```math
  Y_{\ell,m'}(\theta', \phi') = \sum_m \mathfrak{D}^{(\ell)}_{m',m}(R)
  Y_{\ell,m}(\theta, \phi),
  ```
  or, generalizing to spin-weighted spherical harmonics
  ```math
  {}_{s}Y_{\ell,m'}(R_{\theta', \phi'}) = \sum_m \mathfrak{D}^{(\ell)}_{m',m}(R)
  {}_{s}Y_{\ell,m}(R_{\theta, \phi}).
  ```
* We also have that ``\mathfrak{D}`` obeys the representation
  property, so
  ```math
  \mathfrak{D}^{(\ell)}_{m',m''}(R_{\theta', \phi'})
  = \sum_{m} \mathfrak{D}^{(\ell)}_{m',m}(R)
  \mathfrak{D}^{(\ell)}_{m,m''}(R_{\theta, \phi}).
  ```
  - There is no reason that I can see to introduce a conjugation
  - The fact that ``m''`` appears on both sides of the equation means
    that it must correspond to ``s`` — though we have to check the
    behavior under final rotation to determine the sign.

```math
{}_{s}Y_{\ell,m}(R_{\theta, \phi})
\propto
\mathfrak{D}^{(\ell)}_{m,\propto s}(R_{\theta, \phi})
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
L_{\pm} |\ell,m\rangle = \alpha^{\pm}_{\ell,m} |\ell, m \pm 1\rangle
```
to be real and positive.  The reasoning behind this choice is
explained more fully in Section 2 of [Ufford and Shortley
(1932)](@cite UffordShortley_1932).  As a more practical matter, the
Condon-Shortley phase describes signs chosen in the expression for
spherical harmonics.  The key expression is Eq. (15) of section 4³
(page 52) of [Condon-Shortley](@cite CondonShortley_1935):
```math
\Theta(\ell, m) = (-1)^\ell \sqrt{\frac{2\ell+1}{2} \frac{(\ell+m)!}{(\ell-m)!}}
\frac{1}{2^\ell \ell!} \frac{1}{\sin^m\theta}
\frac{d^{\ell-m}}{d(\cos\theta)^{\ell-m}} \sin^{2\ell}\theta.
```
When multiplied by Eq. (5) ``\Phi(m) = e^{im\phi} / \sqrt{2\pi}``,
this gives the spherical harmonic function.  The right-hand side of
the expression above is usually immediately replaced by a simpler
expression using Legendre polynomials, but this just shifts sign
ambiguity into the definition of the Legendre polynomials.  Instead,
we can expand the above expression directly for the first few ``\ell``
values and/or use automatic differentiation to actually test their
original expression as such against the function implemented in this
package.  The first few values are given in a footnote to Condon and
Shortley's Eq. (15) (and have been verified separately by hand and by
computation with SymPy):
```math
\begin{aligned}
\Theta(0,0) &= \sqrt{\frac{1}{2}} \\
\Theta(1,0) &= \sqrt{\frac{3}{2}} \cos\theta &
\Theta(1,\pm1) &= \mp \sqrt{\frac{3}{4}} \sin\theta \\
\Theta(2,0) &= \sqrt{\frac{5}{8}} (2\cos^2\theta - \sin^2\theta) &
\Theta(2,\pm1) &= \mp \sqrt{\frac{15}{4}} \cos\theta \sin\theta &
\Theta(2,\pm2) &= \sqrt{\frac{15}{16}} \sin^2\theta \\
\Theta(3,0) &= \sqrt{\frac{7}{8}} (2\cos^3\theta - 3\cos\theta\sin^2\theta) &
\Theta(3,\pm1) &= \mp \sqrt{\frac{21}{32}} (4\cos^2\theta\sin\theta - \sin^3\theta) &
\Theta(3,\pm2) &= \sqrt{\frac{105}{16}} \cos\theta \sin^2\theta &
\Theta(3,\pm3) &= \mp \sqrt{\frac{35}{32}} \sin^3\theta
\end{aligned}
```
These are tested, along with the results from automatic
differentiation, every time this package is updated.  The result is
perfect agreement, so that we can definitively say that ***the
spherical-harmonic functions provided by this package obey the
Condon-Shortley phase convention.***

## Angular-momentum operators

* First, a couple points about ``-i\hbar``:
  - The finite transformations look like ``\exp[-i \theta L_j]``, but
    the factor of ``i`` introduced here just cancels the one in the
    ``L_j``, and the sign is just chosen to make the result consistent
    with our notion of active or passive transformations.
  - Any factors of ``\hbar`` are included *purely* for the sake of
     convenience.
  - The factor ``i`` comes from plain functional analysis: We need a
    self-adjoint operator, and ``\partial_x`` by itself is
    anti-self-adjoint (as can be verified by evaluating on ``\langle
    x' | x \rangle = \delta(x-x')``, which switches sign based on
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
``\mathbf{L} = \mathbf{x} \times \mathbf{p}``, where ``\mathbf{x}`` is
the position vector and ``\mathbf{p}`` is the momentum vector.  In
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
L_x &= -i\hbar \left( \sin\phi \frac{\partial}{\partial\theta} + \cot\theta \cos\phi \frac{\partial}{\partial\phi} \right), \\
L_y &= -i\hbar \left( \cos\phi \frac{\partial}{\partial\theta} - \cot\theta \sin\phi \frac{\partial}{\partial\phi} \right), \\
L_z &= -i\hbar \frac{\partial}{\partial\phi}.
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
L_{j} f(\mathbf{R}) &\colonequals -z \left. \frac{\partial}{\partial \theta}
f\left(e^{\theta \mathbf{e}_j / 2} \mathbf{R} \right) \right|_{\theta=0}, \\
K_{j} f(\mathbf{R}) &\colonequals -z \left. \frac{\partial}{\partial \theta}
f\left(\mathbf{R} e^{\theta \mathbf{e}_j / 2}\right) \right|_{\theta=0},
\end{aligned}
```
where ``\mathbf{e}_j`` is the unit vector in the ``j`` direction.
Surprisingly, I found that [Edmonds](@cite Edmonds_2016) expresses
essentially the same thing in the equations following his Eq. (4.1.5).

Condon and Shortley's Eq. (1) of section 4³ (page 50) defines
```math
L_z = -i \hbar \frac{\partial}{\partial \phi},
```
while Eq. (8) on the following page defines
```math
\begin{aligned}
L_x + i L_y &= \hbar e^{i\phi} \left( \frac{\partial}{\partial \theta} + i \cot\theta \frac{\partial}{\partial \phi} \right), \\
L_x - i L_y &= \hbar e^{-i\phi} \left(-\frac{\partial}{\partial \theta} + i \cot\theta \frac{\partial}{\partial \phi} \right).
\end{aligned}
```
Note that one is not the conjugate of the other!  This is because of
the factors of ``-i`` in the definitions of ``L_x`` and ``L_y``.

[Edmonds](@cite Edmonds_2016) gives the *total* angular-momentum
operator for a rigid body in Eq. (2.2.2) as
```math
\begin{aligned}
L_x &= -i\hbar \left(-\cos \alpha \cot\beta \frac{\partial}{\partial\alpha} - \sin\alpha \frac{\partial}{\partial\beta} + \frac{\cos\alpha}{\sin\beta} \frac{\partial}{\partial\gamma} \right), \\
L_y &= -i\hbar \left(-\sin\alpha \cot\beta \frac{\partial}{\partial \alpha} + \cos\alpha \frac{\partial}{\partial\beta} + \frac{\sin\alpha}{\sin\beta} \frac{\partial}{\partial\gamma} \right), \\
L_z &= -i\hbar \frac{\partial}{\partial\alpha}.
\end{aligned}
```


## Wigner ``𝔇`` and ``d`` matrices

Wigner's Eqs. (11.18) and (11.19) define the real orthogonal
transformation ``\mathbf{R}`` by
```math
x'_i = R_{ij} x_j
```
and the operator ``\mathbf{P}_{\mathbf{R}}`` to act on a function
``f`` such that
```math
\mathbf{P}_{\mathbf{R}} f(x'_1, \ldots) = f(x_1, \ldots).
```
Then, his Eq. (15.5) presumably implies
```math
Y_{\ell,m}(\vartheta', \varphi')
= \mathbf{P}_{\{\alpha, \beta, \gamma\}} Y_{\ell,m}(\vartheta, \varphi)
= \sum_{m'} \mathfrak{D}^{(\ell)}(\{\alpha, \beta, \gamma\})_{m',m}
  Y_{\ell,m'}(\vartheta, \varphi),
```
where ``\{\alpha, \beta, \gamma\}`` takes ``(\vartheta, \varphi)`` to
``(\vartheta', \varphi')``.  In any case, we can now leave behind this
``\mathbf{P}`` notation and just look at the beginning and end of the
equation above as the critical relationship in Wigner's notation.


Eq. (44b) of [Boyle (2016)](@cite Boyle_2016) says
```math
L_{\pm} \mathfrak{D}^{(\ell)}_{m',m}(\mathbf{R})
= \sqrt{(\ell \mp m')(\ell \pm m' + 1)} \mathfrak{D}^{(\ell)}_{m' \pm 1, m}(\mathbf{R}).
```
while Eq. (21) relates the Wigner D-matrix to the spin-weighted spherical harmonics as
```math
{}_{s}Y_{\ell,m}(\mathbf{R})
= (-1)^s \sqrt{\frac{2\ell+1}{4\pi}} \mathfrak{D}^{(\ell)}_{m,-s}(\mathbf{R}).
```
Plugging the latter into the former, we get
```math
L_{\pm} {}_{s}Y_{\ell,m}(\mathbf{R})
= \sqrt{(\ell \mp m)(\ell \pm m + 1)} {}_{s}Y_{\ell,m \pm 1}(\mathbf{R}).
```
That is, in our conventions we have
```math
\alpha^{\pm}_{\ell,m} = \sqrt{(\ell \mp m)(\ell \pm m + 1)},
```
which is always real and positive, and thus consistent with the Condon-Shortley phase
convention.


### Properties

* ``D^j_{m'm}(\alpha,\beta,\gamma) = (-1)^{m'-m} D^j_{-m',-m}(\alpha,\beta,\gamma)^*``
* ``(-1)^{m'-m}D^{j}_{mm'}(\alpha,\beta,\gamma)=D^{j}_{m'm}(\gamma,\beta,\alpha)``
* ``d_{m',m}^{j}=(-1)^{m-m'}d_{m,m'}^{j}=d_{-m,-m'}^{j}``

```math
\begin{aligned}
d_{m',m}^{j}(\pi)        &= (-1)^{j-m}  \delta_{m',-m} \\[6pt]
d_{m',m}^{j}(\pi-\beta)  &= (-1)^{j+m'}  d_{m',-m}^{j}(\beta)\\[6pt]
d_{m',m}^{j}(\pi+\beta)  &= (-1)^{j-m}  d_{m',-m}^{j}(\beta)\\[6pt]
d_{m',m}^{j}(2\pi+\beta) &= (-1)^{2j}    d_{m',m}^{j}(\beta)\\[6pt]
d_{m',m}^{j}(-\beta)     &= d_{m,m'}^{j}(\beta) = (-1)^{m'-m} d_{m',m}^{j}(\beta)
\end{aligned}
```







```math
\begin{gather}
R = \cos\epsilon + \sin\epsilon\, \hat{\mathfrak{r}} \\
R𝐯 = \cos\epsilon 𝐯 + \sin\epsilon\, \hat{\mathfrak{r}}𝐯 \\
R𝐯R^{-1} = (𝐯\cos\epsilon + \sin\epsilon\, \hat{\mathfrak{r}}𝐯)(\cos\epsilon - \sin\epsilon\, \hat{\mathfrak{r}}) \\
R𝐯R^{-1} = 𝐯\cos^2\epsilon + \sin^2\epsilon\, \hat{\mathfrak{r}}𝐯\hat{\mathfrak{r}}^{-1} + \sin\epsilon \cos\epsilon\, (\hat{\mathfrak{r}}𝐯 - 𝐯\hat{\mathfrak{r}}) \\
R𝐯R^{-1} = \begin{cases}
𝐯 & 𝐯 \hat{\mathfrak{r}} = \hat{\mathfrak{r}}𝐯 \\
𝐯(\cos^2\epsilon - \sin^2\epsilon) + 2 \sin\epsilon \cos\epsilon\, \frac{[\hat{\mathfrak{r}}, 𝐯]}{2} & 𝐯 \hat{\mathfrak{r}} = -\hat{\mathfrak{r}}𝐯 \\
\end{cases} \\
R𝐯R^{-1} = \begin{cases}
𝐯 & 𝐯 \hat{\mathfrak{r}} = \hat{\mathfrak{r}}𝐯 \\
\cos2\epsilon 𝐯 + \sin2\epsilon \frac{[\hat{\mathfrak{r}}, 𝐯]}{2} & 𝐯 \hat{\mathfrak{r}} = -\hat{\mathfrak{r}}𝐯 \\
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
