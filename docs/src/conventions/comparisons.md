# Comparisons

Here, we compare our conventions to other sources, including
references in the literature as well as other software that implements
some of these.  Each of these comparisons is also performed explicitly
in [this package's test
suite](https://github.com/moble/SphericalFunctions.jl/tree/main/test/conventions).

Among the items that would be good to compare are the following, when
actually used by any of these sources:
* Quaternions
  - Order of components
  - Basis and multiplication table
  - Operation as rotations
* Euler angles
* Spherical coordinates
* Angular momentum operators
  - Fundamental definitions
  - Expression in terms of spherical coordinates
  - Expression in terms of Euler angles
  - Right-derivative form
* Spherical harmonics
  - Condon-Shortley phase
  - Formula
* Spin-weighted spherical harmonics
  - Behavior under rotation
* Wigner D-matrices
  - Representation à la $\langle ℓ, m' | e^{-i α J_z} e^{-i β J_y} e^{-i γ J_z} | ℓ, m \rangle$
  - Rotation of spherical harmonics
  - Order of indices
  - Conjugation
  - Function of rotation or inverse rotation
  - Formula

One major result of this is that almost everyone since 1935 has used
the same exact expression for the (scalar) spherical harmonics.

When choosing conventions for this package, I intend to prioritize
consistency (to the extent that any of these references actually have
anything to say about the above items) with the following sources, in
order:

1. LALSuite
2. NINJA
3. Newman-Penrose
4. Goldberg
5. Thorne / MTW
6. Wikipedia
7. Sakurai
8. Shankar
9. Zettili

I think that should be sufficient to find a consensus on conventions
for each of the above — with the possible exception of quaternions,
for which I have my own strong opinions.


## Cohen-Tannoudji (1991)

(moved)


## Condon-Shortley (1935)

(moved)

## Edmonds (1960)

[Edmonds_2016](@citet) is a standard reference for the theory of
angular momentum.

In Sec. 1.3 he actually does a fair job of defining the Euler angles.
The upshot is that his definition agrees with ours, though he uses the
"active" definition style.  That is, the rotations are to be performed
successively in order:

> 1. A rotation ``α(0 \leq α < 2π)`` about the ``z``-axis,
>    bringing the frame of axes from the initial position ``S`` into
>    the position ``S'``.  The axis of this rotation is commonly
>    called the *vertical*.
>
> 2. A rotation ``β(0 \leq β < π)`` about the ``y``-axis of
>    the frame ``S'``, called the *line of nodes*.  Note that its
>    position is in general different from the initial position of the
>    ``y``-axis of the frame ``S``. The resulting position of the
>    frame of axes is symbolized by ``S''``.
>
> 3. A rotation ``γ(0 \leq γ < 2π)`` about the ``z``-axis
>    of the frame of axes ``S''``, called the *figure axis*; the
>    position of this axis depends on the previous rotations
>    ``α`` and ``β``.  The final position of the frame is
>    symbolized by ``S'''``.

I would simply write the "``y``-axis of the frame ``S'``" as ``y'``,
and so on.  In quaternionic language, I would write these rotations as
``\exp[γ 𝐤''/2]\, \exp[β 𝐣'/2]\, \exp[α 𝐤/2]``.  But
we also have
```math
\exp[β 𝐣'/2] = \exp[α 𝐤/2]\, \exp[β 𝐣/2]\, \exp[-α 𝐤/2]
```
so we can just swap the ``α`` rotation with the ``β``
rotation while dropping the prime from ``𝐣'``.  We can do a similar
trick swapping the ``α`` and ``β`` rotations with the
``γ`` rotation while dropping the double prime from ``𝐤''``.
That is, an easy calculation shows that
```math
\exp[γ 𝐤''/2]\, \exp[β 𝐣'/2]\, \exp[α 𝐤/2]
=
\exp[α 𝐤/2]\, \exp[β 𝐣/2]\, \exp[γ 𝐤/2],
```
which is precisely our definition.

The spherical coordinates are implicitly defined by this statement:

> It should be noted that the polar coordinates ``φ, θ``
> with respect to the original frame ``S`` of the ``z``-axis in its
> final position are identical with the Euler angles ``α, β``
> respectively.

Again, this agrees with our definition.

His expression for the angular-momentum operator in Euler angles — Eq.
(2.2.2) — agrees with ours:
```math
\begin{aligned}
L_x &= -i \hbar \left\{
    -\frac{\cos α}{\tan β} \frac{\partial} {\partial α}
    - \sin α \frac{\partial} {\partial β}
    + \frac{\cos α}{\sin β} \frac{\partial} {\partial γ}
\right\},
\\
L_y &= -i \hbar \left\{
    -\frac{\sin α}{\tan β} \frac{\partial} {\partial α}
    + \cos α \frac{\partial} {\partial β}
    +\frac{\sin α}{\sin β} \frac{\partial} {\partial γ}
\right\},
\\
L_z &= -i \hbar \frac{\partial} {\partial α}.
\end{aligned}
```
(The corresponding restriction to spherical coordinates also precisely
agrees with our results, with the extra factor of ``\hbar``.)

Unfortunately, there is disagreement over the definition of the
Wigner D-matrices.  In Eq. (4.1.12) he defines
```math
𝒟_{α β γ} =
\exp\big( \frac{iα}{\hbar} J_z\big)
\exp\big( \frac{iβ}{\hbar} J_y\big)
\exp\big( \frac{iγ}{\hbar} J_z\big),
```
which is the *conjugate* of most other definitions.


## Goldberg et al. (1967)

[GoldbergEtAl_1967](@citet) presented the first paper specifically
about spin-weighted spherical harmonics (after [Newman_1966](@citet)
introduced them), and the first to relate them to the Wigner
D-matrices.

If we relate two vectors by a rotation matrix as ``x'^k = R^{kl}
x^l``, then Goldberg et al. define ``D`` by its action on spherical
harmonics [Eq. (3.3)]:
```math
Y_{ℓ,m}(x') = \sum_{m'} Y_{ℓ,m'}(x) D^{ℓ}_{m',m}\left( R^{-1} \right).
```
They then define the Euler angles as we do, and write [Eq. (3.4)]
```math
D^{ℓ}_{m', m}(α, β, γ)
\equiv
D^{ℓ}_{m', m}\left( R(α β γ)^{-1} \right)
=
e^{i m' γ} d^{ℓ}_{m', m}(β) e^{i m α}.
```
Finally, they derive [Eq. (3.9)]
```math
D^{j}_{m', m}(α, β, γ)
=
\left[\frac{(j+m)!(j-m)!}{(j+m')!(j-m')!}\right]^{1/2}
(\sin \frac{1}{2}β)^{2j}
\sum_r \binom{j+m'}{r} \binom{j-m'}{r-m-m'}
(-1)^{j+m'-r}
e^{imα}
(\cot \tfrac{1}{2}β)^{2r-m-m'}
e^{im'γ}.
```

Equation (3.11) naturally extends to
```math
  {}_sY_{ℓ, m}(θ, ϕ, γ)
  =
  \left[ \left(2ℓ+1\right) / 4π \right]^{1/2}
  D^{ℓ}_{-s,m}(ϕ, θ, γ),
```
where Eq. (3.4) also shows that ``D^{ℓ}_{m', m}(α, β,
γ) = D^{ℓ}_{m', m}(α, β, 0) e^{i m' γ}``,
so we have
```math
  {}_sY_{ℓ, m}(θ, ϕ, γ)
  =
  {}_sY_{ℓ, m}(θ, ϕ)\, e^{-i s γ}.
```
This is the most natural extension of the standard spin-weighted
spherical harmonics to ``\mathrm{Spin}(3)``.  In particular, the
spin-weight operator is ``i \partial_γ``, which suggests that it
will be most natural to choose the sign of ``R_𝐮`` so that ``R_z = i
\partial_γ``.


## Griffiths (1995)

Griffiths' ["Introduction to Quantum Mechanics"](@cite Griffiths_1995)
is probably the most common introductory text used in undergraduate
physics programs, so it would be useful to compare.

Equation (4.27) gives the associated Legendre function as
```math
P_{ℓ}^{m}(x)
=
(1-x^2)^{|m|/2} \left(\frac{d}{dx}\right)^{|m|} P_{ℓ}(x),
```
and (4.28) gives the Legendre polynomial as
```math
P_{ℓ}(x)
=
\frac{1}{2^ℓ ℓ!} \left(\frac{d}{dx}\right)^ℓ (x^2-1)^ℓ.
```
Then, (4.32) gives the spherical harmonics as
```math
Y_{ℓ}^{m}(θ, ϕ)
=
ϵ
\sqrt{\frac{2ℓ+1}{4π} \frac{(ℓ-|m|)!}{(ℓ+|m|)!}}
e^{imϕ} P_{ℓ}^{m}(\cos θ),
```
where ``ϵ = (-1)^m`` for ``m\geq 0`` and ``ϵ = 1`` for
``m\leq 0``.  In Table 4.2, he explicitly lists the first few
spherical harmonics:
```math
\begin{aligned}
  Y_{0}^{0} &= \left(\frac{1}{4π}\right)^{1/2},\\
  Y_{1}^{0} &= \left(\frac{3}{4π}\right)^{1/2} \cos θ,\\
  Y_{1}^{\pm 1} &= \mp \left(\frac{3}{8π}\right)^{1/2} \sin θ e^{\pm iϕ},\\
  Y_{2}^{0} &= \left(\frac{5}{16π}\right)^{1/2} \left(3\cos^2θ - 1\right),\\
  Y_{2}^{\pm 1} &= \mp \left(\frac{15}{8π}\right)^{1/2} \sin θ \cos θ e^{\pm iϕ},\\
  Y_{2}^{\pm 2} &= \left(\frac{15}{32π}\right)^{1/2} \sin^2θ e^{\pm 2iϕ},\\
  Y_{3}^{0} &= \left(\frac{7}{16π}\right)^{1/2} \left(5\cos^3θ - 3\cos θ\right),\\
  Y_{3}^{\pm 1} &= \mp \left(\frac{21}{64π}\right)^{1/2} \sin θ \left(5\cos^2θ - 1\right) e^{\pm iϕ},\\
  Y_{3}^{\pm 2} &= \left(\frac{105}{32π}\right)^{1/2} \sin^2θ \cos θ e^{\pm 2iϕ},\\
  Y_{3}^{\pm 3} &= \mp \left(\frac{35}{64π}\right)^{1/2} \sin^3θ e^{\pm 3iϕ}.
\end{aligned}
```
In Eqs. (4.127)—(4.129), he gives the angular-momentum operators in
terms of spherical coordinates:
```math
\begin{aligned}
L_x &= \frac{\hbar}{i} \left(
    -\sin ϕ \frac{\partial} {\partial θ}
    - \cos ϕ \cot θ \frac{\partial} {\partial ϕ}
\right), \\
L_y &= \frac{\hbar}{i} \left(
    \cos ϕ \frac{\partial} {\partial θ}
    - \sin ϕ \cot θ \frac{\partial} {\partial ϕ}
\right), \\
L_z &= -i \hbar \frac{\partial} {\partial ϕ}.
\end{aligned}
```


## LALSuite

(moved)

## Le Bellac (2006)

[LeBellac_2006](@citet) (with Foreword by Cohen-Tannoudji) takes an
odd approach, defining [Eq. (10.32)]
```math
D^{(j)}_{m', m} \left[ ℛ(θ, ϕ) \right]
=
\langle j, m' | e^{-iϕ J_z} e^{-iθ J_y} | j, m \rangle,
```
but later allowing that ``e^{-i ψ J_z}`` usually goes on the
right-hand side of the others, in which case ``D^{(j)}(θ, ϕ)
\to D^{(j)}(ϕ, θ, ψ)``.  Figure 10.1 shows that the
spherical coordinates are standard (physicist's) coordinates.

Equation (10.65) shows the rotation law:
```math
Y_{ℓ}^{m}\left( ℛ^{-1} \hat{r} \right)
=
\sum_{m'} D^{(ℓ)}_{m', m}(ℛ) Y_{ℓ}^{m'}(\hat{r}),
```
and Eq. (10.66) relates the spherical harmonics to the Wigner
D-matrices:
```math
D^{(ℓ)}_{m, 0}(θ, ϕ)
=
\sqrt{\frac{4π}{2ℓ+1}} \left[Y_{ℓ}^{m}(θ, ϕ)\right]^\ast.
```


## Mathematica

The Euler angles are defined generally such that

> `EulerMatrix[{α,β,γ},{a,b,c}]` is equivalent to ``R_{α,a} R_{β,b} R_{γ,c}``, where
> ``R_{α,a}``=`RotationMatrix[α,UnitVector[3,a]]`, etc.

and

> `EulerMatrix[{α,β,γ}]` is equivalent to `EulerMatrix[{α,β,γ},{3,2,3}]`

(representing the ``z-y-z`` convention).

Finally, we find that they say that `EulerMatrix` corresponds to three rotations:

```mathematica
rα = RotationMatrix[α, {0, 0, 1}];
rβ = RotationMatrix[β, {0, 1, 0}];
rγ = RotationMatrix[γ, {0, 0, 1}];

Simplify[rα . rβ . rγ == EulerMatrix[{α, β, γ}]]
```

This agrees with the conventions used in this package, so we can directly compare
expressions in terms of Euler angles.


We can find conventions at [this
page](https://reference.wolfram.com/language/ref/WignerD.html).

> The Wolfram Language uses phase conventions where ``D^j_{m_1, m_2}(ψ, θ, ϕ) = \exp(i m_1 ψ + i m_2 ϕ) D^j_{m_1, m_2}(0, θ, 0)``.

> `WignerD[{1, 0, 1}, ψ, θ, ϕ]` = ``-\sqrt{2} e^{i ϕ} \cos\frac{θ}{2}
> \sin\frac{θ}{2}``

> `WignerD[{𝓁, 0, m}, θ, ϕ] == Sqrt[(4 π)/(2 𝓁 + 1)] SphericalHarmonicY[𝓁, m, θ, ϕ]`

> `WignerD[{j, m1, m2},ψ, θ, ϕ] == (-1)^(m1 - m2) Conjugate[WignerD[{j, -m1, -m2}, ψ, θ,
> ϕ]]`


> For ``ℓ \geq 0``, ``Y_ℓ^m = \sqrt{(2ℓ+1)/(4π)} \sqrt{(ℓ-m)! / (ℓ+m)!}  P_ℓ^m(\cos θ) e^{imϕ}`` where ``P_ℓ^m`` is the associated Legendre function.

> The associated Legendre polynomials are defined by ``P_n^m(x) = (-1)^m (1-x^2)^{m/2}(d^m/dx^m)P_n(x)`` where ``P_n(x)`` is the Legendre polynomial.

[NIST (14.7.13)](https://dlmf.nist.gov/14.7#E13) gives the Legendre polynomial for nonnegative integer ``n`` as
```math
P_n(x) = \frac{1}{2^n n!} \frac{d^n}{dx^n} (x^2 - 1)^n.
```


## Newman-Penrose

In their 1966 paper, [Newman_1966](@citet), Newman and Penrose first
introduced the spin-weighted spherical harmonics, ``{}_sY_{ℓ m}``.
They use the standard (physicists') convention for spherical
coordinates and introduce the stereographic coordinate ``\zeta =
e^{iϕ} \cot\frac{θ}{2}``.  They define the spin-raising
operator ``\eth`` acting on a function of spin weight ``s`` as
```math
\eth \eta
=
-\left(\sin θ\right)^s
\left\{
    \frac{\partial}{\partial θ}
    + \frac{i}{\sin θ} \frac{\partial}{\partial ϕ}
\right\} \left\{\left(\sin θ\right)^{-s} \eta\right\},
```
They then compute
```math
{}_sY_{ℓ, m}
\propto
\frac{1}{\left[(ℓ-s)! (ℓ+s)!\right]^{1/2}}
\left(1 + \zeta \bar{\zeta}\right)^{-ℓ}
\sum_p \zeta^p (-\bar{\zeta})^{p+s-m}
\binom{ℓ-s}{p} \binom{ℓ+s}{p+s-m},
```
where the sum is over all integers ``p`` such that the factorials are
nonzero.

They are a little ambiguous about the relationship of the complex
basis vector ``m^\mu`` to the coordinates.  

> The vectors ``\Re(m^\mu)`` and ``\Im(m^\mu)`` may be regarded as
> orthogonal tangent vectors (of length ``2^{-1/2}``) at each point of
> the surface. [...] If spherical polar coordinates are used, a
> natural choice for ``m^\mu`` is to make ``\Re(m^\mu)`` and
> ``\Im(m^\mu)`` tangential, respectively, to the curves ``ϕ =
> \mathrm{const}`` and ``θ = \mathrm{const}.``

The ambiguity is in the sign implied by "tangential", but the natural
choice is to assume they mean that the components are *positive*
multiples of ``\partial_θ`` and ``\partial_ϕ`` respectively,
in which case we have
```math
m^\mu = \frac{1}{\sqrt{2}}
\left[ \partial_θ + i \csc θ \partial_ϕ \right].
```
They define the spin weight in terms of behavior of a quantity under
rotation of ``m^\mu`` in its own plane as
```math
(m^\mu)'
=
e^{iψ} m^\mu
=
\frac{1}{\sqrt{2}}
\left[ 
  \left(\cos ψ\partial_θ - \sin ψ\csc θ \partial_ϕ\right)
  + i \left(\cos ψ\csc θ \partial_ϕ + \sin ψ\partial_θ\right)
\right].
```
Raising the spherical coordinates ``(θ, ϕ)`` to Euler angles
``(ϕ, θ, -ψ)``, we see that the rotor ``R_{ϕ, θ,
-ψ}`` rotates the ``𝐳`` basis vector to the point ``(θ,
ϕ)``, and it rotates ``(𝐱 + i 𝐲) / \sqrt{2}`` onto ``(m^\mu)'``.
Under this rotation, a quantity ``\eta`` has spin weight ``s`` if it
transforms as
```math
\eta' = e^{i s ψ} \eta.
```
Now, supposing that these quantities are functions of Euler angles, we
can write
```math
\eta(ϕ, θ, -ψ) = e^{i s ψ} \eta(ϕ, θ, 0),
```
or
```math
\eta(ϕ, θ, γ) = e^{-i s γ} \eta(ϕ, θ, 0).
```
Thus, the operator with eigenvalue ``s`` is ``i \partial_γ``.


## NINJA

(moved)


## Sakurai (1994)




## Scipy


[`scipy.special.sph_harm_y`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_harm_y.html)


## Shankar (1994)

[Shankar_1994](@citet) writes in Eq. (12.5.35) the spherical harmonics
as
```math
Y_{ℓ}^{m}(θ, ϕ)
=
(-1)^ℓ
\left[ \frac{(2ℓ+1)!}{4π} \right]^{1/2}
\frac{1}{2^ℓ ℓ!}
\left[ \frac{(ℓ+m)!}{(2ℓ)!(ℓ-m)!} \right]^{1/2}
e^{i m ϕ}
(\sin θ)^{-m}
\frac{d^{ℓ-m}}{d(\cos θ)^{ℓ-m}}
(\sin θ)^{2ℓ}
```
for ``m \geq 0``, with (12.5.40) giving the expression
```math
Y_{ℓ}^{-m}(θ, ϕ)
=
(-1)^m \left( Y_{ℓ}^{m}(θ, ϕ) \right)^\ast.
```
The angular-momentum operators are given below (12.5.27) as
```math
\begin{aligned}
L_x &= i \hbar \left(
    \sin ϕ \frac{\partial} {\partial θ}
    + \cos ϕ \cot θ \frac{\partial} {\partial ϕ}
\right),
\\
L_y &= i \hbar \left(
    -\cos ϕ \frac{\partial} {\partial θ}
    + \sin ϕ \cot θ \frac{\partial} {\partial ϕ}
\right),
\\
L_z &= -i \hbar \frac{\partial} {\partial ϕ}.
\end{aligned}
```
In Exercise 12.5.7, the rotation operator is defined by
```math
U\left[ R(α, β, γ) \right]
=
e^{-i α J_z/\hbar}
e^{-i β J_y/\hbar}
e^{-i γ J_z/\hbar},
```
That ``U`` becomes a ``D^{(j)}`` when the operator is acting on the
states ``|j, m\rangle`` for a given ``j``.  Thus, while Shankar never
actually uses notation like ``D^{(j)}_{m', m}``, he does talk about
``\langle j, m' | D^{(j)}\left[ R(α, β, γ) \right] | j,
m \rangle``.


## SymPy

There is no specific Euler angle convention in SymPy, however it is
informative to see what the `sympy.algebras.Quaternion.from_euler`
class method does.  You can specify 

SymPy uses what I would consider just a wrong expression for ``D``.
Specifically, the
[source](https://github.com/sympy/sympy/blob/b4ce69ad5d40e4e545614b6c76ca9b0be0b98f0b/sympy/physics/wigner.py#L1136-L1191)
cites [Edmonds_2016](@citet) when defining
```math
𝒟_{α β γ} =
\exp\big( \frac{iα}{\hbar} J_z\big)
\exp\big( \frac{iβ}{\hbar} J_y\big)
\exp\big( \frac{iγ}{\hbar} J_z\big).
```
But that is an incorrect copy of Edmonds' Eq. (4.1.9), in which the
``α`` and ``γ`` on the right-hand side are swapped.  The
code also implements D in the `wigner_d` function as (essentially)
```python
exp(I*mprime*alpha)*d[i, j]*exp(I*m*gamma)
```
even though the actual equation Eq. (4.1.12) of Edmonds says
```math
\mathscr{D}^{(j)}_{m' m}(α β γ) =
\exp i m' γ d^{(j)}_{m' m}(α, β) \exp(i m α).
```
The ``d`` matrix appears to be implemented consistently with Edmonds,
and thus not affected.

Basically, it appears that SymPy just swapped the order of the Euler
angles relative to Edmonds, who already introduced a conjugate to the
definition of the D matrix.

## Thorne

## Torres del Castillo (2003)

[TorresDelCastillo_2003](@citet) starts by defining a rotation
``ℛ`` as transforming a point ``x_i`` into another point
with coordinates ``x_i' = a_{ij}x_j``.  Under that rotation, any
scalar function ``f`` transforms into another function ``f' =
ℛ f`` defined by [Eq. (2.43)]
```math
f'\big(x_i\big) = f\big( a^{-1}_{ij} x_j \big).
```
In particular, ``f'(x'_i) = f(x_i)``.  He then defines Wigner's
D-matrix to satisfy [Eq. (2.45)]
```math
ℛ Y_{l,m} = \sum_{m} D^l_{m',m}(ℛ) Y_{l,m'}.
```
Including the arguments to the spherical harmonics, this becomes
```math
Y_{l,m}\big(ℛ^{-1} R_{θ, ϕ}\big)
=
\sum_{m} D^l_{m',m}(ℛ) Y_{l,m'}\big(R_{θ, ϕ}\big).
```
In this form, we have [Eq. (2.46)]
```math
D^l_{m'',m}(ℛ_1 ℛ_2)
=
\sum_{m'} D^l_{m'',m'}(ℛ_1) D^l_{m',m}(ℛ_2).
```
He computes [Eq. (2.53)]
```math
D^l_{m',m}(ϕ, θ, \chi)
=
e^{-i m' ϕ} d^l_{m',m}(θ) e^{-i m \chi},
```
where the ``d`` matrix is given by
```math
d^l_{m',m}(θ)
=
\sqrt{(l+m)!(l-m)!(l+m')!(l-m')!}
\sum_{k} \frac{
  (-1)^k
  (\sin \tfrac{1}{2} θ)^{m-m'+2k}
  (\cos \tfrac{1}{2} θ)^{2l-m+m'-2k}
} {
  k!(l+m'-k)!(l-m-k)!(m-m'+k)!
},
```
and the spin-weighted spherical harmonic is related to ``D`` by
```math
{}_{s}Y_{j,m}(θ, ϕ)
=
(-1)^m
\sqrt{\frac{2j+1}{4π}}
d^j_{-m,s}(θ)
e^{i m ϕ}.
```


## Varshalovich et al. (1988)

[Varshalovich_1988](@citet) has a fairly decent comparison of
definitions related to the rotation matrix by previous authors.  

Eq. 1.4.(31) defines the operator
```math
\hat{D}(α, β, γ)
=
e^{-iα \hat{J}_z}
e^{-iβ \hat{J}_y}
e^{-iγ \hat{J}_z},
```
where the ``\hat{J}`` operators are defined in 

> In quantum mechanics the total angular momentum operator ``\hat{J}``
> is defined as an operator which generates transformations of wave
> functions (state vectors) and quantum operators under infinitesimal
> rotations of the coordinate system (see Eqs. 2.1.(1) and 2.1.(2)).
> 
> A transformation of an arbitrary wave function ``\Psi`` under
> rotation of the coordinate system through an infinitesimal angle
> ``δ \omega`` about an axis ``𝐧`` may be written as
> ```math
> \Psi \to \Psi' = \left(1 - i δ \omega 𝐧 \cdot \hat{J} \right)\Psi,
> ```
> where ``\hat{J}`` is the total angular momentum operator.

Eq. 4.1.(1) defines the Wigner D-functions according to
```math
\langle J M | \hat{D}(α, β, γ) | J' M' \rangle
=
δ_{J J'} D^J_{M M'}(α, β, γ).
```
Eq. 4.3.(1) states
```math
D^J_{M M'}(α, β, γ)
=
e^{-i M α}
d^J_{M M'}(β)
e^{-i M' γ}
```


Page 155 has a table of values for ``ℓ \leq 5``

[Varshalovich_1988](@citet) distinguish in Sec. 1.1.3 between
*covariant* and *contravariant* spherical coordinates and the
corresponding basis vectors, which they define as
```math
\begin{aligned}
  𝐞_{+1} &= - \frac{1}{\sqrt{2}} \left( 𝐞_x + i 𝐞_y\right)
  &&&
  𝐞^{+1} &= - \frac{1}{\sqrt{2}} \left( 𝐞_x - i 𝐞_y\right) \\
  𝐞_{0} &= 𝐞_z &&& 𝐞^{0} &= 𝐞_z \\
  𝐞_{-1} &= \frac{1}{\sqrt{2}} \left( 𝐞_x - i 𝐞_y\right)
  &&&
  𝐞^{-1} &= \frac{1}{\sqrt{2}} \left( 𝐞_x + i 𝐞_y\right).
\end{aligned}
```
Then, in Sec. 4.2 they define ``\hat{𝐉}`` as the operator of
angular momentum of the rigid symmetric top.  They then give in Eq.
(6) the "covariant spherical coordinates of ``\hat{𝐉}`` in the
non-rotating (lab-fixed) system" as
```math
\begin{gather}
  \hat{J}_{\pm 1} = \frac{i}{\sqrt{2}} e^{\pm i α} \left[
    \mp \cot β \frac{\partial}{\partial α}
    + i \frac{\partial}{\partial β}
    \pm \frac{1}{\sin β} \frac{\partial}{\partial γ}
  \right] \\
  \hat{J}_0 = - i \frac{\partial}{\partial α},
\end{gather}
```
and "contravariant components of ``\hat{𝐉}`` in the rotating
(body-fixed) system" as
```math
\begin{gather}
  \hat{J}'^{\pm 1} = \frac{i}{\sqrt{2}} e^{\mp i γ} \left[
    \pm \cot β \frac{\partial}{\partial γ}
    + i \frac{\partial}{\partial β}
    \mp \frac{1}{\sin β} \frac{\partial}{\partial α}
  \right] \\
  \hat{J}'^0 = - i \frac{\partial}{\partial γ}.
\end{gather}
```
(Note the prime in the last two equations.)  We can expand these in
Cartesian components to compare to our expressions.  First the
covariant components:
```math
\begin{aligned}
  \hat{J}_{x}
  &= -\frac{1}{\sqrt{2}} \left( \hat{J}_{+1} - \hat{J}_{-1} \right) \\
  % &= -\frac{1}{\sqrt{2}} \left( 
  %   \frac{i}{\sqrt{2}} e^{i α} \left[
  %     - \cot β \frac{\partial}{\partial α}
  %     + i \frac{\partial}{\partial β}
  %     + \frac{1}{\sin β} \frac{\partial}{\partial γ}
  %   \right]
  %   -
  %   \frac{i}{\sqrt{2}} e^{-i α} \left[
  %     + \cot β \frac{\partial}{\partial α}
  %     + i \frac{\partial}{\partial β}
  %     - \frac{1}{\sin β} \frac{\partial}{\partial γ}
  %   \right]
  % \right) \\
  &= i\left[ 
      \frac{\cos α}{\tan β} \frac{\partial}{\partial α}
      + \sin α \frac{\partial}{\partial β}
      - \frac{\cos α}{\sin β} \frac{\partial}{\partial γ}
  \right] \\
  \hat{J}_{y}
  &= -\frac{1}{i\sqrt{2}} \left( \hat{J}_{+1} + \hat{J}_{-1} \right) \\
  % &= -\frac{1}{i\sqrt{2}} \left( 
  %   \frac{i}{\sqrt{2}} e^{i α} \left[
  %     - \cot β \frac{\partial}{\partial α}
  %     + i \frac{\partial}{\partial β}
  %     + \frac{1}{\sin β} \frac{\partial}{\partial γ}
  %   \right]
  %   +
  %   \frac{i}{\sqrt{2}} e^{-i α} \left[
  %     + \cot β \frac{\partial}{\partial α}
  %     + i \frac{\partial}{\partial β}
  %     - \frac{1}{\sin β} \frac{\partial}{\partial γ}
  %   \right]
  % \right) \\
  &= i \left[
      \frac{\sin α}{\tan β} \frac{\partial}{\partial α}
      - \cos α \frac{\partial}{\partial β}
      - \frac{\sin α}{\sin β} \frac{\partial}{\partial γ}
  \right] \\
  \hat{J}_{z}
  &= \hat{J}_{0} \\
  &= -i \frac{\partial}{\partial α}
\end{aligned}
```
We can compare these to the [Full expressions on ``𝕊³``]() `@ref`, and find
that they are precisely equivalent to expressions for ``L_j`` computed in
this package's conventions.

Next, the contravariant components:
```math
\begin{aligned}
  \hat{J}'_{x}
  &= -\frac{1}{\sqrt{2}} \left( \hat{J}'^{+1} - \hat{J}'^{-1} \right) \\
  % &= -\frac{1}{\sqrt{2}} \left(
  %   \frac{i}{\sqrt{2}} e^{- i γ} \left[
  %     + \cot β \frac{\partial}{\partial γ}
  %     + i \frac{\partial}{\partial β}
  %     - \frac{1}{\sin β} \frac{\partial}{\partial α}
  %   \right]
  %   -
  %   \frac{i}{\sqrt{2}} e^{+ i γ} \left[
  %     - \cot β \frac{\partial}{\partial γ}
  %     + i \frac{\partial}{\partial β}
  %     + \frac{1}{\sin β} \frac{\partial}{\partial α}
  %   \right]
  % \right) \\
  &= -i \left(
      \frac{\cos γ}{\tan β} \frac{\partial}{\partial γ}
      + \sin γ \frac{\partial}{\partial β}
      - \frac{\cos γ}{\sin β} \frac{\partial}{\partial α}
  \right) \\
  \hat{J}'_{y}
  &= \frac{1}{i\sqrt{2}} \left( \hat{J}'^{+1} + \hat{J}'^{-1} \right) \\
  % &= \frac{1}{i\sqrt{2}} \left(
  %   \frac{i}{\sqrt{2}} e^{-i γ} \left[
  %     + \cot β \frac{\partial}{\partial γ}  
  %     + i \frac{\partial}{\partial β}
  %     - \frac{1}{\sin β} \frac{\partial}{\partial α}
  %   \right]
  %   +
  %   \frac{i}{\sqrt{2}} e^{+ i γ} \left[
  %     - \cot β \frac{\partial}{\partial γ}
  %     + i \frac{\partial}{\partial β}
  %     + \frac{1}{\sin β} \frac{\partial}{\partial α}
  %   \right]
  % \right) \\
  &= -i \left(
      \frac{\sin γ}{\tan β} \frac{\partial}{\partial γ}
      - \cos γ \frac{\partial}{\partial β}
      - \frac{\sin γ}{\sin β} \frac{\partial}{\partial α}
  \right) \\
  \hat{J}'_{z}
  &= \hat{J}'^{0} \\
  &= -i \frac{\partial}{\partial γ}
\end{aligned}
```
Unfortunately, while we have agreement on ``\hat{J}'^{y} = R_y``, we
also have disagreement on ``\hat{J}'^{x} = -R_x`` and ``\hat{J}'^{z} =
-R_z``, as they have relative minus signs.

It's very easy to check, for example, that ``[\hat{J}'^{z},
\hat{J}'^{x}] = i \hat{J}'^{y}``, as expected from the general
expression in their Eq. (12).  So these expressions are — at least —
consistent with the claims of Varshalovich et al.  I wonder if there
is some subtlety involving the order of operations and passing to the
"body-fixed" frame.  I'm confident that my definitions are internally
consistent, and fit in nicely with the spin-weighted function
literature; maybe Varshalovich et al. are just doing something
different.

## Wikipedia

Euler angles

Angular-momentum operators

Spherical harmonics

Spin-weighted spherical harmonics

Defining the operator
```math
ℛ(α,β,γ) = e^{-iα J_z}e^{-iβ J_y}e^{-iγ J_z},
```
[Wikipedia expresses the Wigner D-matrix](https://en.wikipedia.org/wiki/Wigner_D-matrix#Definition_of_the_Wigner_D-matrix) as
```math
D^j_{m'm}(α,β,γ) \equiv \langle jm' | ℛ(α,β,γ)| jm \rangle =e^{-im'α } d^j_{m'm}(β)e^{-i mγ}.
```


## Wigner




## Zettili (2009)

[Zettili_2009](@citet) is a relatively recent textbook that seems to
be gaining popularity.  (Note that there is a 3rd edition from 2022,
but I do not have access to it; all the references here are to the 2nd
edition from 2009.)

In Appendix B.1, we find that the spherical coordinates are related to
Cartesian coordinates in the usual (physicist's) way.  Equation
(5.132) gives the angular-momentum operator
```math
\hat{L}_z = -i \hbar \frac{\partial}{\partial φ},
```
which agrees with [our expression](@ref "``L`` operators in spherical
coordinates").  This is followed by equation (5.134):
```math
\hat{L}_{\pm}
=
\hat{L}_x \pm i \hat{L}_y
=
\pm \hbar e^{\pm iφ} \left(
  \frac{\partial}{\partial θ}
  \pm i \cot θ \frac{\partial}{\partial φ}
\right),
```
which also agrees with [our results.](@ref "``L_{\pm}`` operators in
spherical coordinates")

Equation (5.180) gives the spherical harmonics as
```math
Y_{l, m}(θ, φ)
=
\frac{(-1)^l}{2^l l!}
\sqrt{\frac{2l+1}{4π} \frac{(l+m)!}{(l-m)!}}
e^{imφ}
\frac{1}{\sin^m θ}
\frac{d^{l-m}}{d(\cos θ)^{l-m}}
(\sin θ)^{2l}.
```

Section 7.2.1 denotes by ``\hat{R}_z(δ ϕ)`` the

> rotation of the coordinates of a *spinless* particle over an
> *infinitesimal* angle ``δ ϕ`` about the ``z``-axis

and shows its action [Eq. (7.16)]
```math
\hat{R}_z (δ ϕ) ψ(r, θ, ϕ)
=
ψ(r, θ, ϕ - δ ϕ).
```

> We may generalize this relation to a rotation of angle ``δ
> ϕ`` about an arbitrary axis whose direction is given by the unit
> vector ``\vec{n}``:

```math
\hat{R}(δ ϕ)
=
1 - \frac{i}{\hbar} δ ϕ \vec{n} \cdot \hat{\vec{L}}.
```
This extends to finite rotation by defining the operator [Eq. (7.48)]
```math
\hat{R}(α, β, γ)
=
e^{-iα J_z / \hbar} e^{-iβ J_y / \hbar} e^{-iγ J_z / \hbar}.
```
Equation (7.52) then defines
```math
D^{(j)}_{m', m}(α, β, γ)
=
\langle j, m' | \hat{R}(α, β, γ) | j, m \rangle,
```
So that [Eq. (7.54)]
```math
D^{(j)}_{m', m}(α, β, γ)
=
e^{-i (m' α + m γ)} d^{(j)}_{m', m}(β),
```
where [Eq. (7.55)]
```math
d^{(j)}_{m', m}(β)
=
\langle j, m' | e^{-iβ J_y / \hbar} | j, m \rangle.
```
The explicit expression for ``d`` is [Eq. (7.56)]
```math
d^{(j)}_{m', m}(β)
=
\sum_k (-1)^{k+m'-m}
\frac{\sqrt{(j+m)!(j-m)!(j+m')!(j-m')!}}
{(j-m'-k)!(j+m-k)!(k+m'-m)!k!}
\left(\cos\frac{β}{2}\right)^{2j+m-m'-2k}
\left(\sin\frac{β}{2}\right)^{m'-m+2k}.
```
In Sec. 7.2.6, we find that if the operator ``\hat{R}(α, β,
γ)`` rotates a vector pointing in the ``(θ, ϕ)`` to a
vector pointing in the ``(θ', ϕ')`` direction, then the
spherical harmonics transform as [Eq. (7.70)]
```math
Y_{ℓ, m}^\ast (θ', ϕ')
=
\sum_{m'} D^{(ℓ)}_{m, m'}(α, β, γ) Y_{ℓ, m'}^\ast (θ, ϕ).
```
