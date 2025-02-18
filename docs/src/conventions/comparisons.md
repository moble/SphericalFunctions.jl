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
  - Representation à la $\langle \ell, m' | e^{-i \alpha J_z} e^{-i \beta J_y} e^{-i \gamma J_z} | \ell, m \rangle$
  - Rotation of spherical harmonics
  - Order of indices
  - Conjugation
  - Function of rotation or inverse rotation
  - Formula

One major result of this is that almost everyone since 1935 has used
the same exact expression for the (scalar) spherical harmonics.

When choosing my conventions, I intend to prioritize consistency (to
the extent that any of these references actually have anything to say
about the above items) with the following sources, in order:

1. LALSuite
2. NINJA
3. Thorne / MTW
4. Goldberg
5. Newman-Penrose
6. Wikipedia
7. Sakurai
8. Shankar
9. Zettili

I think that should be sufficient to find a consensus on conventions
for each of the above — with the possible exception of quaternions,
for which I have my own strong opinions.


## Cohen-Tannoudji (1991)

[CohenTannoudji_1991](@citet) defines spherical coordinates in the
usual (physicist's) way in Chapter VI.  He then computes the
angular-momentum operators as [Eqs. (D-5)]
```math
\begin{aligned}
L_x &= i \hbar \left(
    \sin\phi \frac{\partial} {\partial \theta}
    + \frac{\cos\phi}{\tan\theta} \frac{\partial} {\partial \phi}
\right),
\\
L_y &= i \hbar \left(
    -\cos\phi \frac{\partial} {\partial \theta}
    + \frac{\sin\phi}{\tan\theta} \frac{\partial} {\partial \phi}
\right),
\\
L_z &= \frac{\hbar}{i} \frac{\partial} {\partial \phi}.
\end{aligned}
```


He derives the spherical harmonics in two
ways and gets two different, but equivalent, expressions in Complement
``\mathrm{A}_{\mathrm{VI}}``.  The first is Eq. (26)
```math
Y_{l}^{m}(\theta, \phi)
=
\frac{(-1)^l}{2^l l!} \sqrt{\frac{(2l+1)}{4\pi} \frac{(l+m)!}{(l-m)!}}
e^{i m \phi} (\sin \theta)^{-m}
\frac{d^{l-m}}{d(\cos \theta)^{l-m}} (\sin \theta)^{2l},
```
while the second is Eq. (30)
```math
Y_{l}^{m}(\theta, \phi)
=
\frac{(-1)^{l+m}}{2^l l!} \sqrt{\frac{(2l+1)}{4\pi} \frac{(l-m)!}{(l+m)!}}
e^{i m \phi} (\sin \theta)^m
\frac{d^{l+m}}{d(\cos \theta)^{l+m}} (\sin \theta)^{2l}.
```

In Complement ``\mathrm{B}_{\mathrm{VI}}`` he defines a rotation
operator ``R`` as acting on a state such that [Eq. (21)]
```math
\langle \mathbf{r} | R | \psi \rangle
=
\langle \mathscr{R}^{-1} \mathbf{r} | \psi \rangle.
```
For an infinitesimal rotation through angle ``d\alpha`` about the axis
``\mathbf{u}``, he shows [Eq. (49)]
```math
R_{\mathbf{u}}(d\alpha) = 1 - \frac{i}{\hbar} d\alpha \mathbf{L}.\mathbf{u}.
```

Cohen-Tannoudji does not appear to define the Wigner D-matrices.


## Condon-Shortley

## Edmonds (1960)

[Edmonds_2016](@citet) is a standard reference for the theory of
angular momentum.

In Sec. 1.3 he actually does a fair job of defining the Euler angles.
The upshot is that his definition agrees with ours, though he uses the
"active" definition style.  That is, the rotations are to be performed
successively in order:

> 1. A rotation ``\alpha(0 \leq \alpha < 2\pi)`` about the ``z``-axis,
>    bringing the frame of axes from the initial position ``S`` into
>    the position ``S'``.  The axis of this rotation is commonly
>    called the *vertical*.
>
> 2. A rotation ``\beta(0 \leq \beta < \pi)`` about the ``y``-axis of
>    the frame ``S'``, called the *line of nodes*.  Note that its
>    position is in general different from the initial position of the
>    ``y``-axis of the frame ``S``. The resulting position of the
>    frame of axes is symbolized by ``S''``.
>
> 3. A rotation ``\gamma(0 \leq \gamma < 2\pi)`` about the ``z``-axis
>    of the frame of axes ``S''``, called the *figure axis*; the
>    position of this axis depends on the previous rotations
>    ``\alpha`` and ``\beta``.  The final position of the frame is
>    symbolized by ``S'''``.

I would simply write the "``y``-axis of the frame ``S'``" as ``y'``,
and so on.  In quaternionic language, I would write these rotations as
``\exp[\gamma 𝐤''/2]\, \exp[\beta 𝐣'/2]\, \exp[\alpha 𝐤/2]``.  But
we also have
```math
\exp[\beta 𝐣'/2] = \exp[\alpha 𝐤/2]\, \exp[\beta 𝐣/2]\, \exp[-\alpha 𝐤/2]
```
so we can just swap the ``\alpha`` rotation with the ``\beta``
rotation while dropping the prime from ``𝐣'``.  We can do a similar
trick swapping the ``\alpha`` and ``\beta`` rotations with the
``\gamma`` rotation while dropping the double prime from ``𝐤''``.
That is, an easy calculation shows that
```math
\exp[\gamma 𝐤''/2]\, \exp[\beta 𝐣'/2]\, \exp[\alpha 𝐤/2]
=
\exp[\alpha 𝐤/2]\, \exp[\beta 𝐣/2]\, \exp[\gamma 𝐤/2],
```
which is precisely our definition.

The spherical coordinates are implicitly defined by

> It should be noted that the polar coordinates ``\varphi, \theta``
> with respect to the original frame ``S`` of the ``z``-axis in its
> final position are identical with the Euler angles ``\alpha, \beta``
> respectively.

Again, this agrees with our definition.

His expression for the angular-momentum operator in Euler angles — Eq.
(2.2.2) — agrees with ours:
```math
\begin{aligned}
L_x &= -i \hbar \left\{
    -\frac{\cos\alpha}{\tan\beta} \frac{\partial} {\partial \alpha}
    - \sin\alpha \frac{\partial} {\partial \beta}
    + \frac{\cos\alpha}{\sin\beta} \frac{\partial} {\partial \gamma}
\right\},
\\
L_y &= -i \hbar \left\{
    -\frac{\sin\alpha}{\tan\beta} \frac{\partial} {\partial \alpha}
    + \cos\alpha \frac{\partial} {\partial \beta}
    +\frac{\sin\alpha}{\sin\beta} \frac{\partial} {\partial \gamma}
\right\},
\\
L_z &= -i \hbar \frac{\partial} {\partial \alpha}.
\end{aligned}
```
(The corresponding restriction to spherical coordinates also precisely
agrees with our results, with the extra factor of ``\hbar``.)

Unfortunately, there is disagreement over the definition of the
Wigner D-matrices.  In Eq. (4.1.12) he defines
```math
\mathcal{D}_{\alpha \beta \gamma} =
\exp\big( \frac{i\alpha}{\hbar} J_z\big)
\exp\big( \frac{i\beta}{\hbar} J_y\big)
\exp\big( \frac{i\gamma}{\hbar} J_z\big),
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
Y_{ell,m}(x') = \sum_{m'} Y_{ell,m'}(x) D^{ell}_{m',m}\left( R^{-1} \right).
```
They then define the Euler angles as we do, and write [Eq. (3.4)]
```math
D^{\ell}_{m', m}(\alpha, \beta, \gamma)
\equiv
D^{\ell}_{m', m}\left( R(\alpha \beta \gamma)^{-1} \right)
=
e^{i m' \gamma} d^{\ell}_{m', m}(\beta) e^{i m \alpha}.
```
Finally, they derive [Eq. (3.9)]
```math
D^{j}_{m', m}(\alpha, \beta, \gamma)
=
\left[\frac{(j+m)!(j-m)!}{(j+m')!(j-m')!}\right]^{1/2}
(\sin \frac{1}{2}\beta)^{2j}
\sum_r \binom{j+m'}{r} \binom{j-m'}{r-m-m'}
(-1)^{j+m'-r}
e^{im\alpha}
(\cot \tfrac{1}{2}\beta)^{2r-m-m'}
e^{im'\gamma}.
```

Equation (3.11) naturally extends to
```math
  {}_sY_{\ell, m}(\theta, \phi, \gamma)
  =
  \left[ \left(2\ell+1\right) / 4\pi \right]^{1/2}
  D^{\ell}_{-s,m}(\phi, \theta, \gamma),
```
where Eq. (3.4) also shows that ``D^{\ell}_{m', m}(\alpha, \beta,
\gamma) = D^{\ell}_{m', m}(\alpha, \beta, 0) e^{i m' \gamma}``,
so we have
```math
  {}_sY_{\ell, m}(\theta, \phi, \gamma)
  =
  {}_sY_{\ell, m}(\theta, \phi)\, e^{-i s \gamma}.
```
This is the most natural extension of the standard spin-weighted
spherical harmonics to ``\mathrm{Spin}(3)``.  In particular, the
spin-weight operator is ``i \partial_\gamma``, which suggests that it
will be most natural to choose the sign of ``R_𝐮`` so that ``R_z = i
\partial_\gamma``.


## Griffiths (1995)

Griffiths' ["Introduction to Quantum Mechanics"](@cite Griffiths_1995)
is probably the most common introductory text used in undergraduate
physics programs, so it would be useful to compare.

Equation (4.27) gives the associated Legendre function as
```math
P_{\ell}^{m}(x)
=
(1-x^2)^{|m|/2} \left(\frac{d}{dx}\right)^{|m|} P_{\ell}(x),
```
and (4.28) gives the Legendre polynomial as
```math
P_{\ell}(x)
=
\frac{1}{2^\ell \ell!} \left(\frac{d}{dx}\right)^\ell (x^2-1)^\ell.
```
Then, (4.32) gives the spherical harmonics as
```math
Y_{\ell}^{m}(\theta, \phi)
=
\epsilon
\sqrt{\frac{2\ell+1}{4\pi} \frac{(\ell-|m|)!}{(\ell+|m|)!}}
e^{im\phi} P_{\ell}^{m}(\cos\theta),
```
where ``\epsilon = (-1)^m`` for ``m\geq 0`` and ``\epsilon = 1`` for
``m\leq 0``.  In Table 4.2, he explicitly lists the first few
spherical harmonics:
```math
\begin{aligned}
  Y_{0}^{0} &= \left(\frac{1}{4\pi}\right)^{1/2},\\
  Y_{1}^{0} &= \left(\frac{3}{4\pi}\right)^{1/2} \cos\theta,\\
  Y_{1}^{\pm 1} &= \mp \left(\frac{3}{8\pi}\right)^{1/2} \sin\theta e^{\pm i\phi},\\
  Y_{2}^{0} &= \left(\frac{5}{16\pi}\right)^{1/2} \left(3\cos^2\theta - 1\right),\\
  Y_{2}^{\pm 1} &= \mp \left(\frac{15}{8\pi}\right)^{1/2} \sin\theta \cos\theta e^{\pm i\phi},\\
  Y_{2}^{\pm 2} &= \left(\frac{15}{32\pi}\right)^{1/2} \sin^2\theta e^{\pm 2i\phi}.
  Y_{3}^{0} &= \left(\frac{7}{16\pi}\right)^{1/2} \left(5\cos^3\theta - 3\cos\theta\right),\\
  Y_{3}^{\pm 1} &= \mp \left(\frac{21}{64\pi}\right)^{1/2} \sin\theta \left(5\cos^2\theta - 1\right) e^{\pm i\phi},\\
  Y_{3}^{\pm 2} &= \left(\frac{105}{32\pi}\right)^{1/2} \sin^2\theta \cos\theta e^{\pm 2i\phi},\\
  Y_{3}^{\pm 3} &= \mp \left(\frac{35}{64\pi}\right)^{1/2} \sin^3\theta e^{\pm 3i\phi}.
\end{aligned}
```
In Eqs. (4.127)—(4.129), he gives the angular-momentum operators in
terms of spherical coordinates:
```math
\begin{aligned}
L_x &= \frac{\hbar}{i} \left(
    -\sin\phi \frac{\partial} {\partial \theta}
    - \cos\phi \cot\theta \frac{\partial} {\partial \phi}
\right), \\
L_y &= \frac{\hbar}{i} \left(
    \cos\phi \frac{\partial} {\partial \theta}
    - \sin\phi \cot\theta \frac{\partial} {\partial \phi}
\right), \\
L_z &= -i \hbar \frac{\partial} {\partial \phi}.
\end{aligned}
```


## LALSuite

[LALSuite (LSC Algorithm Library Suite)](@cite LALSuite_2018) is a
collection of software routines, comprising the primary official
software used by the LIGO-Virgo-KAGRA Collaboration to detect and
characterize gravitational waves.  As far as I can tell, the ultimate
source for all spin-weighted spherical harmonic values used in
LALSuite is the function
[`XLALSpinWeightedSphericalHarmonic`](https://git.ligo.org/lscsoft/lalsuite/-/blob/6e653c91b6e8a6728c4475729c4f967c9e09f020/lal/lib/utilities/SphericalHarmonics.c),
which cites the NINJA paper [Ajith_2007](@citet) as its source.
Unfortunately, it cites version *1*, which contained a serious error,
using ``\tfrac{\cos\iota}{2}`` instead of ``\cos \tfrac{\iota}{2}``
and similarly for ``\sin``.  This error was corrected in version 2,
but the citation was not updated.  Nonetheless, it appears that the
actual code is consistent with the corrected versions of the NINJA
paper; the equivalence is
[tested](https://github.com/moble/SphericalFunctions.jl/blob/0f57c77e65da85e4996f0969fe0a931b460135ac/test/wigner_matrices/sYlm.jl#L59)
in this package's test suite.

## Le Bellac (2006)

[LeBellac_2006](@citet) (with Foreword by Cohen-Tannoudji) takes an
odd approach, defining [Eq. (10.32)]
```math
D^{(j)}_{m', m} \left[ \mathcal{R}(\theta, \phi) \right]
=
\langle j, m' | e^{-i\phi J_z} e^{-i\theta J_y} | j, m \rangle,
```
but later allowing that ``e^{-i \psi J_z}`` usually goes on the
right-hand side of the others, in which case ``D^{(j)}(\theta, \phi)
\to D^{(j)}(\phi, \theta, \psi)``.  Figure 10.1 shows that the
spherical coordinates are standard (physicist's) coordinates.

Equation (10.65) shows the rotation law:
```math
Y_{\ell}^{m}\left( \mathcal{R}^{-1} \hat{r} \right)
=
\sum_{m'} D^{(\ell)}_{m', m}(\mathcal{R}) Y_{\ell}^{m'}(\hat{r}),
```
and Eq. (10.66) relates the spherical harmonics to the Wigner
D-matrices:
```math
D^{(\ell)}_{m, 0}(\theta, \phi)
=
\sqrt{\frac{4\pi}{2\ell+1}} \left[Y_{\ell}^{m}(\theta, \phi)\right]^\ast.
```


## Mathematica

The Euler angles are defined generally such that

> `EulerMatrix[{α,β,γ},{a,b,c}]` is equivalent to ``R_{α,a} R_{β,b} R_{γ,c}``, where
> ``R_{α,a}``=`RotationMatrix[α,UnitVector[3,a]]`, etc.

and

> `EulerMatrix[{α,β,γ}]` is equivalent to `EulerMatrix[{α,β,γ},{3,2,3}]`

(representing the ``z-y-z`` convention).

Finally, we find that they say that `EulerMatrix`` corresponds to three rotations:

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

> The Wolfram Language uses phase conventions where ``D^j_{m_1, m_2}(\psi, \theta, \phi) = \exp(i m_1 \psi + i m_2 \phi) D^j_{m_1, m_2}(0, \theta, 0)``.

> `WignerD[{1, 0, 1}, ψ, θ, ϕ]` = ``-\sqrt{2} e^{i \phi} \cos\frac{\theta}{2}
> \sin\frac{\theta}{2}``

> `WignerD[{𝓁, 0, m}, θ, ϕ] == Sqrt[(4 π)/(2 𝓁 + 1)] SphericalHarmonicY[𝓁, m, θ, ϕ]`

> `WignerD[{j, m1, m2},ψ, θ, ϕ] == (-1)^(m1 - m2) Conjugate[WignerD[{j, -m1, -m2}, ψ, θ,
> ϕ]]`


> For ``\ell \geq 0``, ``Y_\ell^m = \sqrt{(2\ell+1)/(4\pi)} \sqrt{(\ell-m)! / (\ell+m)!}  P_\ell^m(\cos \theta) e^{im\phi}`` where ``P_\ell^m`` is the associated Legendre function.

> The associated Legendre polynomials are defined by ``P_n^m(x) = (-1)^m (1-x^2)^{m/2}(d^m/dx^m)P_n(x)`` where ``P_n(x)`` is the Legendre polynomial.

[NIST (14.7.13)](https://dlmf.nist.gov/14.7#E13) gives the Legendre polynomial for nonnegative integer ``n`` as
```math
P_n(x) = \frac{1}{2^n n!} \frac{d^n}{dx^n} (x^2 - 1)^n.
```


## Newman-Penrose

In their 1966 paper, [Newman_1966](@citet), Newman and Penrose first
introduced the spin-weighted spherical harmonics, ``{}_sY_{\ell m}``.
They use the standard (physicists') convention for spherical
coordinates and introduce the stereographic coordinate ``\zeta =
e^{i\phi} \cot\frac{\theta}{2}``.  They define the spin-raising
operator ``\eth`` acting on a function of spin weight ``s`` as
```math
\eth \eta
=
-\left(\sin\theta\right)^s
\left\{
    \frac{\partial}{\partial\theta}
    + \frac{i}{\sin\theta} \frac{\partial}{\partial\phi}
\right\} \left\{\left(\sin\theta\right)^{-s} \eta\right\},
```
They then compute
```math
{}_sY_{\ell, m}
\propto
\frac{1}{\left[(\ell-s)! (\ell+s)!\right]^{1/2}}
\left(1 + \zeta \bar{\zeta}\right)^{-\ell}
\sum_p \zeta^p (-\bar{\zeta})^{p+s-m}
\binom{\ell-s}{p} \binom{\ell+s}{p+s-m},
```
where the sum is over all integers ``p`` such that the factorials are
nonzero.

They are a little ambiguous about the relationship of the complex
basis vector ``m^\mu`` to the coordinates.  

> The vectors ``\Re(m^\mu)`` and ``\Im(m^\mu)`` may be regarded as
> orthogonal tangent vectors (of length ``2^{-1/2}``) at each point of
> the surface. [...] If spherical polar coordinates are used, a
> natural choice for ``m^\mu`` is to make ``\Re(m^\mu)`` and
> ``\Im(m^\mu)`` tangential, respectively, to the curves ``\phi =
> \mathrm{const}`` and ``\theta = \mathrm{const}.``

The ambiguity is in the sign implied by "tangential", but the natural
choice is to assume they mean that the components are *positive*
multiples of ``\partial_\theta`` and ``\partial_\phi`` respectively,
in which case we have
```math
m^\mu = \frac{1}{\sqrt{2}}
\left[ \partial_\theta + i \csc\theta \partial_\phi \right].
```
They define the spin weight in terms of behavior of a quantity under
rotation of ``m^\mu`` in its own plane as
```math
(m^\mu)'
=
e^{i\psi} m^\mu
=
\frac{1}{\sqrt{2}}
\left[ 
  \left(\cos\psi\partial_\theta - \sin\psi\csc\theta \partial_\phi\right)
  + i \left(\cos\psi\csc\theta \partial_\phi + \sin\psi\partial_\theta\right)
\right].
```
Raising the spherical coordinates ``(\theta, \phi)`` to Euler angles
``(\phi, \theta, -\psi)``, we see that the rotor ``R_{\phi, \theta,
-\psi}`` rotates the ``𝐳`` basis vector to the point ``(\theta,
\phi)``, and it rotates ``(𝐱 + i 𝐲) / \sqrt{2}`` onto ``(m^\mu)'``.
Under this rotation, a quantity ``\eta`` has spin weight ``s`` if it
transforms as
```math
\eta' = e^{i s \psi} \eta.
```
Now, supposing that these quantities are functions of Euler angles, we
can write
```math
\eta(\phi, \theta, -\psi) = e^{i s \psi} \eta(\phi, \theta, 0),
```
or
```math
\eta(\phi, \theta, \gamma) = e^{-i s \gamma} \eta(\phi, \theta, 0).
```
Thus, the operator with eigenvalue ``s`` is ``i \partial_\gamma``.


## NINJA

[Ajith_2007](@citet) was prepared by a broad cross-section of
researchers (including the author of this package) involved in
modeling gravitational waves with the intent of providing a shared set
of conventions.  The spherical coordinates are standard physicist's
coordinates, except that the polar angle is denoted ``\iota``.
Equation (II.7) is
```math
  {}^{-s}Y_{l,m} = (-1)^s\sqrt{\frac{2\ell+1}{4\pi}}
  d^{\ell}_{m,s}(\iota)e^{im\phi},
```
where
```math
  d^{\ell}_{m,s}(\iota)
  =
  \sum_{k = k_1}^{k_2}
  \frac{(-1)^k[(\ell+m)!(\ell-m)!(\ell+s)!(\ell-s)!]^{1/2}}
  {(\ell+m-k)!(\ell-s-k)!k!(k+s-m)!}
  \left(\cos\left(\frac{\iota}{2}\right)\right)^{2\ell+m-s-2k}
  \left(\sin\left(\frac{\iota}{2}\right)\right)^{2k+s-m}
```
with ``k_1 = \textrm{max}(0, m-s)`` and ``k_2=\textrm{min}(\ell+m,
\ell-s)``. For reference, they provide several values [Eqs.
(II.9)--(II.13)]:
```math
\begin{aligned}
  {}^{-2}Y_{2,2} &= \sqrt{\frac{5}{64\pi}}(1+\cos\iota)^2e^{2i\phi},\\
  {}^{-2}Y_{2,1} &= \sqrt{\frac{5}{16\pi}}  \sin\iota( 1 + \cos\iota )e^{i\phi},\\
  {}^{-2}Y_{2,0} &= \sqrt{\frac{15}{32\pi}} \sin^2\iota,\\
  {}^{-2}Y_{2,-1} &= \sqrt{\frac{5}{16\pi}}  \sin\iota( 1 - \cos\iota
  )e^{-i\phi},\\
  {}^{-2}Y_{2,-2} &=& \sqrt{\frac{5}{64\pi}}(1-\cos\iota)^2e^{-2i\phi}.
\end{aligned}
```
Note that most of the above was copied directly from the TeX source of
the paper.  Also note the annoying negative sign on the left-hand side
of the first expression.  Getting rid of it and combining the first
two expressions, we have the full formula for the spin-weighted
spherical harmonics in this convention:
```math
\begin{aligned}
  {}_{s}Y_{lm}
  &=
  (-1)^s\sqrt{\frac{2\ell+1}{4\pi}} e^{im\phi}
  \sum_{k = k_1}^{k_2}
  \frac{(-1)^k[(\ell+m)!(\ell-m)!(\ell-s)!(\ell+s)!]^{1/2}}
  {(\ell+m-k)!(\ell+s-k)!k!(k-s-m)!}
  \\ &\qquad \times
  \left(\cos\left(\frac{\iota}{2}\right)\right)^{2\ell+m+s-2k}
  \left(\sin\left(\frac{\iota}{2}\right)\right)^{2k-s-m}
\end{aligned}
```
where ``k_1 = \textrm{max}(0, m+s)`` and ``k_2=\textrm{min}(\ell+m,
\ell+s)``.


## Sakurai (1994)




## Scipy


[`scipy.special.sph_harm_y`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_harm_y.html)


## Shankar (1994)

[Shankar_1994](@citet) writes in Eq. (12.5.35) the spherical harmonics
as
```math
Y_{\ell}^{m}(\theta, \phi)
=
(-1)^\ell
\left[ \frac{(2\ell+1)!}{4\pi} \right]^{1/2}
\frac{1}{2^\ell \ell!}
\left[ \frac{(\ell+m)!}{(2\ell)!(\ell-m)!} \right]^{1/2}
e^{i m \phi}
(\sin \theta)^{-m}
\frac{d^{\ell-m}}{d(\cos\theta)^{\ell-m}}
(\sin\theta)^{2\ell}
```
for ``m \geq 0``, with (12.5.40) giving the expression
```math
Y_{\ell}^{-m}(\theta, \phi)
=
(-1)^m \left( Y_{\ell}^{m}(\theta, \phi) \right)^\ast.
```
The angular-momentum operators are given below (12.5.27) as
```math
\begin{aligned}
L_x &= i \hbar \left(
    \sin\phi \frac{\partial} {\partial \theta}
    + \cos\phi \cot\theta \frac{\partial} {\partial \phi}
\right),
\\
L_y &= i \hbar \left(
    -\cos\phi \frac{\partial} {\partial \theta}
    + \sin\phi \cot\theta \frac{\partial} {\partial \phi}
\right),
\\
L_z &= -i \hbar \frac{\partial} {\partial \phi}.
\end{aligned}
```
In Exercise 12.5.7, the rotation operator is defined by
```math
U\left[ R(\alpha, \beta, \gamma) \right]
=
e^{-i \alpha J_z/\hbar}
e^{-i \beta J_y/\hbar}
e^{-i \gamma J_z/\hbar},
```
That ``U`` becomes a ``D^{(j)}`` when the operator is acting on the
states ``|j, m\rangle`` for a given ``j``.  Thus, while Shankar never
actually uses notation like ``D^{(j)}_{m', m}``, he does talk about
``\langle j, m' | D^{(j)}\left[ R(\alpha, \beta, \gamma) \right] | j,
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
\mathcal{D}_{\alpha \beta \gamma} =
\exp\big( \frac{i\alpha}{\hbar} J_z\big)
\exp\big( \frac{i\beta}{\hbar} J_y\big)
\exp\big( \frac{i\gamma}{\hbar} J_z\big).
```
But that is an incorrect copy of Edmonds' Eq. (4.1.9), in which the
``\alpha`` and ``\gamma`` on the right-hand side are swapped.  The
code also implements D in the `wigner_d` function as (essentially)
```python
exp(I*mprime*alpha)*d[i, j]*exp(I*m*gamma)
```
even though the actual equation Eq. (4.1.12) says
```math
\mathscr{D}^{(j)}_{m' m}(\alpha \beta \gamma) =
\exp i m' \gamma d^{(j)}_{m' m}(\alpha, \beta) \exp(i m \alpha).
```
The ``d`` matrix appears to be implemented consistently with Edmonds,
and thus not affected.

Basically, it appears that SymPy just swapped the order of the Euler
angles relative to Edmonds, who already introduced a conjugate to the
definition of the D matrix.

## Thorne

## Torres del Castillo (2003)

[TorresDelCastillo_2003](@citet) starts by defining a rotation
``\mathcal{R}`` as transforming a point ``x_i`` into another point
with coordinates ``x_i' = a_{ij}x_j``.  Under that rotation, any
scalar function ``f`` transforms into another function ``f' =
\mathcal{R} f`` defined by [Eq. (2.43)]
```math
f'\big(x_i\big) = f\big( a^{-1}_{ij} x_j \big).
```
In particular, ``f'(x'_i) = f(x_i)``.  He then defines Wigner's
D-matrix to satisfy [Eq. (2.45)]
```math
\mathcal{R} Y_{l,m} = \sum_{m} D^l_{m',m}(\mathcal{R}) Y_{l,m'}.
```
Including the arguments to the spherical harmonics, this becomes
```math
Y_{l,m}\big(\mathcal{R}^{-1} R_{\theta, \phi}\big)
=
\sum_{m} D^l_{m',m}(\mathcal{R}) Y_{l,m'}\big(R_{\theta, \phi}\big).
```
In this form, we have [Eq. (2.46)]
```math
D^l_{m'',m}(\mathcal{R}_1 \mathcal{R}_2)
=
\sum_{m'} D^l_{m'',m'}(\mathcal{R}_1) D^l_{m',m}(\mathcal{R}_2).
```
He computes [Eq. (2.53)]
```math
D^l_{m',m}(\phi, \theta, \chi)
=
e^{-i m' \phi} d^l_{m',m}(\theta) e^{-i m \chi},
```
where the ``d`` matrix is given by
```math
d^l_{m',m}(\theta)
=
\sqrt{(l+m)!(l-m)!(l+m')!(l-m')!}
\sum_{k} \frac{
  (-1)^k
  (\sin \tfrac{1}{2} \theta)^{m-m'+2k}
  (\cos \tfrac{1}{2} \theta)^{2l-m+m'-2k}
} {
  k!(l+m'-k)!(l-m-k)!(m-m'+k)!
},
```
and the spin-weighted spherical harmonic is related to ``D`` by
```math
{}_{s}Y_{j,m}(\theta, \phi)
=
(-1)^m
\sqrt{\frac{2j+1}{4\pi}}
d^j_{-m,s}(\theta)
e^{i m \phi}.
```


## Varshalovich et al. (1988)

[Varshalovich_1988](@citet) has a fairly decent comparison of
definitions related to the rotation matrix by previous authors.  

Eq. 1.4.(31) defines the operator
```math
\hat{D}(\alpha, \beta, \gamma)
=
e^{-i\alpha \hat{J}_z}
e^{-i\beta \hat{J}_y}
e^{-i\gamma \hat{J}_z},
```
where the ``\hat{J}`` operators are defined in 

> In quantum mechanics the total angular momentum operator ``\hat{J}``
> is defined as an operator which generates transformations of wave
> functions (state vectors) and quantum operators under infinitesimal
> rotations of the coordinate system (see Eqs. 2.1.(1) and 2.1.(2)).
> 
> A transformation of an arbitrary wave function ``\Psi`` under
> rotation of the coordinate system through an infinitesimal angle
> ``\delta \omega`` about an axis ``\mathbf{n}`` may be written as
> ```math
> \Psi \to \Psi' = \left(1 - i \delta \omega \mathbf{n} \cdot \hat{J} \right)\Psi,
> ```
> where ``\hat{J}`` is the total angular momentum operator.

Eq. 4.1.(1) defines the Wigner D-functions according to
```math
\langle J M | \hat{D}(\alpha, \beta, \gamma) | J' M' \rangle
=
\delta_{J J'} D^J_{M M'}(\alpha, \beta, \gamma).
```
Eq. 4.3.(1) states
```math
D^J_{M M'}(\alpha, \beta, \gamma)
=
e^{-i M \alpha}
d^J_{M M'}(\beta)
e^{-i M' \gamma}
```


Page 155 has a table of values for ``\ell \leq 5``

[Varshalovich_1988](@citet) distinguish in Sec. 1.1.3 between
*covariant* and *contravariant* spherical coordinates and the
corresponding basis vectors, which they define as
```math
\begin{aligned}
  \mathbf{e}_{+1} &= - \frac{1}{\sqrt{2}} \left( \mathbf{e}_x + i \mathbf{e}_y\right)
  &&&
  \mathbf{e}^{+1} &= - \frac{1}{\sqrt{2}} \left( \mathbf{e}_x - i \mathbf{e}_y\right) \\
  \mathbf{e}_{0} &= \mathbf{e}_z &&& \mathbf{e}^{0} &= \mathbf{e}_z \\
  \mathbf{e}_{-1} &= \frac{1}{\sqrt{2}} \left( \mathbf{e}_x - i \mathbf{e}_y\right)
  &&&
  \mathbf{e}^{-1} &= \frac{1}{\sqrt{2}} \left( \mathbf{e}_x + i \mathbf{e}_y\right).
\end{aligned}
```
Then, in Sec. 4.2 they define ``\hat{\mathbf{J}}`` as the operator of
angular momentum of the rigid symmetric top.  They then give in Eq.
(6) the "covariant spherical coordinates of ``\hat{\mathbf{J}}`` in the
non-rotating (lab-fixed) system" as
```math
\begin{gather}
  \hat{J}_{\pm 1} = \frac{i}{\sqrt{2}} e^{\pm i \alpha} \left[
    \mp \cot\beta \frac{\partial}{\partial \alpha}
    + i \frac{\partial}{\partial \beta}
    \pm \frac{1}{\sin\beta} \frac{\partial}{\partial \gamma}
  \right] \\
  \hat{J}_0 = - i \frac{\partial}{\partial \alpha},
\end{gather}
```
and "contravariant components of ``\hat{\mathbf{J}}`` in the rotating
(body-fixed) system" as
```math
\begin{gather}
  \hat{J}'^{\pm 1} = \frac{i}{\sqrt{2}} e^{\mp i \gamma} \left[
    \pm \cot\beta \frac{\partial}{\partial \gamma}
    + i \frac{\partial}{\partial \beta}
    \mp \frac{1}{\sin\beta} \frac{\partial}{\partial \alpha}
  \right] \\
  \hat{J}'^0 = - i \frac{\partial}{\partial \gamma}.
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
  %   \frac{i}{\sqrt{2}} e^{i \alpha} \left[
  %     - \cot\beta \frac{\partial}{\partial \alpha}
  %     + i \frac{\partial}{\partial \beta}
  %     + \frac{1}{\sin\beta} \frac{\partial}{\partial \gamma}
  %   \right]
  %   -
  %   \frac{i}{\sqrt{2}} e^{-i \alpha} \left[
  %     + \cot\beta \frac{\partial}{\partial \alpha}
  %     + i \frac{\partial}{\partial \beta}
  %     - \frac{1}{\sin\beta} \frac{\partial}{\partial \gamma}
  %   \right]
  % \right) \\
  &= i\left[ 
      \frac{\cos\alpha}{\tan\beta} \frac{\partial}{\partial \alpha}
      + \sin\alpha \frac{\partial}{\partial \beta}
      - \frac{\cos\alpha}{\sin\beta} \frac{\partial}{\partial \gamma}
  \right] \\
  \hat{J}_{y}
  &= -\frac{1}{i\sqrt{2}} \left( \hat{J}_{+1} + \hat{J}_{-1} \right) \\
  % &= -\frac{1}{i\sqrt{2}} \left( 
  %   \frac{i}{\sqrt{2}} e^{i \alpha} \left[
  %     - \cot\beta \frac{\partial}{\partial \alpha}
  %     + i \frac{\partial}{\partial \beta}
  %     + \frac{1}{\sin\beta} \frac{\partial}{\partial \gamma}
  %   \right]
  %   +
  %   \frac{i}{\sqrt{2}} e^{-i \alpha} \left[
  %     + \cot\beta \frac{\partial}{\partial \alpha}
  %     + i \frac{\partial}{\partial \beta}
  %     - \frac{1}{\sin\beta} \frac{\partial}{\partial \gamma}
  %   \right]
  % \right) \\
  &= i \left[
      \frac{\sin\alpha}{\tan\beta} \frac{\partial}{\partial \alpha}
      - \cos\alpha \frac{\partial}{\partial \beta}
      - \frac{\sin\alpha}{\sin\beta} \frac{\partial}{\partial \gamma}
  \right] \\
  \hat{J}_{z}
  &= \hat{J}_{0} \\
  &= -i \frac{\partial}{\partial \alpha}
\end{aligned}
```
We can compare these to the [Full expressions on ``S^3``](@ref), and find
that they are precisely equivalent to expressions for ``L_j`` computed in
this package's conventions.

Next, the contravariant components:
```math
\begin{aligned}
  \hat{J}'_{x}
  &= -\frac{1}{\sqrt{2}} \left( \hat{J}'^{+1} - \hat{J}'^{-1} \right) \\
  % &= -\frac{1}{\sqrt{2}} \left(
  %   \frac{i}{\sqrt{2}} e^{- i \gamma} \left[
  %     + \cot\beta \frac{\partial}{\partial \gamma}
  %     + i \frac{\partial}{\partial \beta}
  %     - \frac{1}{\sin\beta} \frac{\partial}{\partial \alpha}
  %   \right]
  %   -
  %   \frac{i}{\sqrt{2}} e^{+ i \gamma} \left[
  %     - \cot\beta \frac{\partial}{\partial \gamma}
  %     + i \frac{\partial}{\partial \beta}
  %     + \frac{1}{\sin\beta} \frac{\partial}{\partial \alpha}
  %   \right]
  % \right) \\
  &= -i \left(
      \frac{\cos\gamma}{\tan\beta} \frac{\partial}{\partial \gamma}
      + \sin\gamma \frac{\partial}{\partial \beta}
      - \frac{\cos\gamma}{\sin\beta} \frac{\partial}{\partial \alpha}
  \right) \\
  \hat{J}'_{y}
  &= \frac{1}{i\sqrt{2}} \left( \hat{J}'^{+1} + \hat{J}'^{-1} \right) \\
  % &= \frac{1}{i\sqrt{2}} \left(
  %   \frac{i}{\sqrt{2}} e^{-i \gamma} \left[
  %     + \cot\beta \frac{\partial}{\partial \gamma}  
  %     + i \frac{\partial}{\partial \beta}
  %     - \frac{1}{\sin\beta} \frac{\partial}{\partial \alpha}
  %   \right]
  %   +
  %   \frac{i}{\sqrt{2}} e^{+ i \gamma} \left[
  %     - \cot\beta \frac{\partial}{\partial \gamma}
  %     + i \frac{\partial}{\partial \beta}
  %     + \frac{1}{\sin\beta} \frac{\partial}{\partial \alpha}
  %   \right]
  % \right) \\
  &= -i \left(
      \frac{\sin\gamma}{\tan\beta} \frac{\partial}{\partial \gamma}
      - \cos\gamma \frac{\partial}{\partial \beta}
      - \frac{\sin\gamma}{\sin\beta} \frac{\partial}{\partial \alpha}
  \right) \\
  \hat{J}'_{z}
  &= \hat{J}'^{0} \\
  &= -i \frac{\partial}{\partial \gamma}
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

Defining the operator
```math
\mathcal{R}(\alpha,\beta,\gamma) = e^{-i\alpha J_z}e^{-i\beta J_y}e^{-i\gamma J_z},
```
[Wikipedia defines the Wigner D-matrix](https://en.wikipedia.org/wiki/Wigner_D-matrix#Definition_of_the_Wigner_D-matrix) as
```math
D^j_{m'm}(\alpha,\beta,\gamma) \equiv \langle jm' | \mathcal{R}(\alpha,\beta,\gamma)| jm \rangle =e^{-im'\alpha } d^j_{m'm}(\beta)e^{-i m\gamma}.
```


## Wigner




## Zettili (2009)

[Zettili_2009](@citet) denotes by ``\hat{R}_z(\delta \phi)`` the

> rotation of the coordinates of a *spinless* particle over an
> *infinitesimal* angle ``\delta \phi`` about the ``z``-axis

and shows its action [Eq. (7.16)]
```math
\hat{R}_z (\delta \phi) \psi(r, \theta, \phi)
=
\psi(r, \theta, \phi - \delta \phi).
```

> We may generalize this relation to a rotation of angle ``\delta
> \phi`` about an arbitrary axis whose direction is given by the unit
> vector ``\vec{n}``:

```math
\hat{R}(\delta \phi)
=
1 - \frac{i}{\hbar} \delta \phi \vec{n} \cdot \hat{\vec{L}}.
```
This extends to finite rotation by defining the operator [Eq. (7.48)]
```math
\hat{R}(\alpha, \beta, \gamma)
=
e^{-i\alpha J_z / \hbar} e^{-i\beta J_y / \hbar} e^{-i\gamma J_z / \hbar}.
```
Equation (7.52) then defines
```math
D^{(j)}_{m', m}(\alpha, \beta, \gamma)
=
\langle j, m' | \hat{R}(\alpha, \beta, \gamma) | j, m \rangle,
```
So that [Eq. (7.54)]
```math
D^{(j)}_{m', m}(\alpha, \beta, \gamma)
=
e^{-i (m' \alpha + m \gamma)} d^{(j)}_{m', m}(\beta),
```
where [Eq. (7.55)]
```math
d^{(j)}_{m', m}(\beta)
=
\langle j, m' | e^{-i\beta J_y / \hbar} | j, m \rangle.
```
The explicit expression for ``d`` is [Eq. (7.56)]
```math
d^{(j)}_{m', m}(\beta)
=
\sum_k (-1)^{k+m'-m}
\frac{\sqrt{(j+m)!(j-m)!(j+m')!(j-m')!}}
{(j-m'-k)!(j+m-k)!(k+m'-m)!k!}
\left(\cos\frac{\beta}{2}\right)^{2j+m-m'-2k}
\left(\sin\frac{\beta}{2}\right)^{m'-m+2k}.
```
In Sec. 7.2.6, we find that if the operator ``\hat{R}(\alpha, \beta,
\gamma)`` rotates a vector pointing in the ``(\theta, \phi)`` to a
vector pointing in the ``(\theta', \phi')`` direction, then the
spherical harmonics transform as [Eq. (7.70)]
```math
Y_{\ell, m}^\ast (\theta', \phi')
=
\sum_{m'} D^{(\ell)}_{m, m'}(\alpha, \beta, \gamma) Y_{\ell, m'}^\ast (\theta, \phi).
```

In Appendix B.1, we find that the spherical coordinates are related to
Cartesian coordinates in the usual (physicist's) way, and Eqs.
(B.25)—(B.27) give the components of the angular-momentum operator in
spherical coordinates as
```math
\begin{aligned}
L_x &= i \hbar \left(
    \sin\phi \frac{\partial}{\partial \theta}
    + \cot\theta \cos\phi \frac{\partial}{\partial \phi}
\right),
\\
L_y &= i \hbar \left(
    -\cos\phi \frac{\partial}{\partial \theta}
    + \cot\theta \sin\phi \frac{\partial}{\partial \phi}
\right),
\\
L_z &= -i \hbar \frac{\partial}{\partial \phi}.
\end{aligned}
```