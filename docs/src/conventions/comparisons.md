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
  - Basis
  - Operation as rotations
* Euler angles
* Spherical coordinates
* Spherical harmonics
  - Condon-Shortley phase
  - Formula
* Spin-weighted spherical harmonics
  - Behavior under rotation
* Wigner D-matrices
  - Order of indices
  - Conjugation
  - Function of rotation or inverse rotation
  - Formula

One major result of this is that almost everyone since 1935 has used
the same exact expression for the (scalar) spherical harmonics.

## Condon-Shortley

## Wigner

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


## Goldberg

Eq. (3.11) of [GoldbergEtAl_1967](@citet) naturally extends to
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

## Wikipedia

## Mathematica

## SymPy

## Sakurai

## Thorne

## Torres del Castillo

## NINJA

Combining Eqs. (II.7) and (II.8) of [Ajith_2007](@citet), we have
```math
\begin{align}
  {}_{-s}Y_{lm}
  &=
  (-1)^s\sqrt{\frac{2\ell+1}{4\pi}} e^{im\phi}
  \sum_{k = k_1}^{k_2}
  \frac{(-1)^k[(\ell+m)!(\ell-m)!(\ell+s)!(\ell-s)!]^{1/2}}
  {(\ell+m-k)!(\ell-s-k)!k!(k+s-m)!}
  \\ &\qquad \times
  \left(\cos\left(\frac{\iota}{2}\right)\right)^{2\ell+m-s-2k}
  \left(\sin\left(\frac{\iota}{2}\right)\right)^{2k+s-m}
\end{align}
```
with ``k_1 = \textrm{max}(0, m-s)`` and ``k_2=\textrm{min}(\ell+m,
\ell-s)``.  Note that most of the above was copied directly from the
TeX source of the paper, but the two equations were trivially combined
into one.  Also note the annoying negative sign on the left-hand side.
That's so annoying that I'm going to duplicate the expression just to
get rid of it:
```math
\begin{align}
  {}_{s}Y_{lm}
  &=
  (-1)^s\sqrt{\frac{2\ell+1}{4\pi}} e^{im\phi}
  \sum_{k = k_1}^{k_2}
  \frac{(-1)^k[(\ell+m)!(\ell-m)!(\ell-s)!(\ell+s)!]^{1/2}}
  {(\ell+m-k)!(\ell+s-k)!k!(k-s-m)!}
  \\ &\qquad \times
  \left(\cos\left(\frac{\iota}{2}\right)\right)^{2\ell+m+s-2k}
  \left(\sin\left(\frac{\iota}{2}\right)\right)^{2k-s-m}
\end{align}
```
where ``k_1 = \textrm{max}(0, m+s)`` and ``k_2=\textrm{min}(\ell+m,
\ell+s)``.


## LALSuite

## Varshalovich et al.

Page 155 has a table of values for ``\ell \leq 5``

[Varshalovich_1988](@citet) distinguish in Sec. 1.1.3 between
*covariant* and *contravariant* spherical coordinates and the
corresponding basis vectors, which they define as
```math
\begin{align}
  \mathbf{e}_{+1} &= - \frac{1}{\sqrt{2}} \left( \mathbf{e}_x + i \mathbf{e}_y\right)
  &&&
  \mathbf{e}^{+1} &= - \frac{1}{\sqrt{2}} \left( \mathbf{e}_x - i \mathbf{e}_y\right) \\
  \mathbf{e}_{0} &= \mathbf{e}_z &&& \mathbf{e}^{0} &= \mathbf{e}_z \\
  \mathbf{e}_{-1} &= \frac{1}{\sqrt{2}} \left( \mathbf{e}_x - i \mathbf{e}_y\right)
  &&&
  \mathbf{e}^{-1} &= \frac{1}{\sqrt{2}} \left( \mathbf{e}_x + i \mathbf{e}_y\right).
\end{align}
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
\begin{align}
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
\end{align}
```
We can compare these to the [Full expressions on ``S^3``](@ref), and find
that they are precisely equivalent to expressions for ``L_j`` computed in
this package's conventions.

Next, the contravariant components:
```math
\begin{align}
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
\end{align}
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
