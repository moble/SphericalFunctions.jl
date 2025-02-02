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

## Goldberg

## Wikipedia

## Mathematica

## SymPy

## Sakurai

## Thorne

## Torres del Castillo

## NINJA

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
Unfortunately, while ``\hat{J}'^{x} = R_x`` and ``\hat{J}'^{z} =
R_z``, we have ``\hat{J}'^{y} = -R_y`` with an unexplained relative
minus sign.

```math
\begin{align}
  \hat{J}'_{x}
  &= -i \left(
      \frac{\cos\gamma}{\tan\beta} \frac{\partial}{\partial \gamma}
      + \sin\gamma \frac{\partial}{\partial \beta}
      - \frac{\cos\gamma}{\sin\beta} \frac{\partial}{\partial \alpha}
  \right) \\
  \hat{J}'_{y}
  &= -i \left(
      \frac{\sin\gamma}{\tan\beta} \frac{\partial}{\partial \gamma}
      - \cos\gamma \frac{\partial}{\partial \beta}
      - \frac{\sin\gamma}{\sin\beta} \frac{\partial}{\partial \alpha}
  \right)
\end{align}
```
We can write out these two operators acting on a function `f` in SymPy:
```python
from sympy import symbols, sin, cos, tan, diff, I
α, β, γ = symbols("α, β, γ", real=True)
f = symbols("f", cls=Function)

def Jx(f, α, β, γ):
  return -I * (
      cos(γ) / tan(β) * diff(f(α, β, γ), γ)
      + sin(γ) * diff(f(α, β, γ), β)
      - cos(γ) / sin(β) * diff(f(α, β, γ), α)
  )
def Jy(f, α, β, γ):
  return -I * (
    sin(γ) / tan(β) * diff(f(α, β, γ), γ)
    - cos(γ) * diff(f(α, β, γ), β)
    - sin(γ) / sin(β) * diff(f(α, β, γ), α)
  )
(
    Jx(lambda α, β, γ: Jy(f, α, β, γ), α, β, γ)
    - Jy(lambda α, β, γ: Jx(f, α, β, γ), α, β, γ)
).expand().simplify()
```

Just to check that we have the right expression, let's check
``[\hat{J}'_{x}, \hat{J}'_{y}] = i \hat{J}'_{z}`` (Eq. 12).  We can
skip any contributions where a derivative on the left will act on a
derivative on the right, and ``\partial_\alpha`` on the left will have
no other effect, so we just have six terms subtracted from six terms.

```math
\begin{align}
  [\hat{J}'_{x}, \hat{J}'_{y}]
  &= \left[
    -i \left(
      \frac{\cos\gamma}{\tan\beta} \frac{\partial}{\partial \gamma}
      + \sin\gamma \frac{\partial}{\partial \beta}
      - \frac{\cos\gamma}{\sin\beta} \frac{\partial}{\partial \alpha}
    \right),
    -i \left(
      \frac{\sin\gamma}{\tan\beta} \frac{\partial}{\partial \gamma}
      - \cos\gamma \frac{\partial}{\partial \beta}
      - \frac{\sin\gamma}{\sin\beta} \frac{\partial}{\partial \alpha}
    \right)
  \right] \\
  &= -\left[
    \frac{\cos\gamma}{\tan\beta} \frac{\partial}{\partial \gamma}
    \left(
      \frac{\sin\gamma}{\tan\beta} \frac{\partial}{\partial \gamma}
      - \cos\gamma \frac{\partial}{\partial \beta}
      - \frac{\sin\gamma}{\sin\beta} \frac{\partial}{\partial \alpha}
    \right)
    + \sin\gamma \frac{\partial}{\partial \beta} \left(
      \frac{\sin\gamma}{\tan\beta} \frac{\partial}{\partial \gamma}
      - \cos\gamma \frac{\partial}{\partial \beta}
      - \frac{\sin\gamma}{\sin\beta} \frac{\partial}{\partial \alpha}
    \right)
    -
    \frac{\sin\gamma}{\tan\beta} \frac{\partial}{\partial \gamma} \left(
      \frac{\cos\gamma}{\tan\beta} \frac{\partial}{\partial \gamma}
      + \sin\gamma \frac{\partial}{\partial \beta}
      - \frac{\cos\gamma}{\sin\beta} \frac{\partial}{\partial \alpha}
    \right)
    + \cos\gamma \frac{\partial}{\partial \beta} \left(
      \frac{\cos\gamma}{\tan\beta} \frac{\partial}{\partial \gamma}
      + \sin\gamma \frac{\partial}{\partial \beta}
      - \frac{\cos\gamma}{\sin\beta} \frac{\partial}{\partial \alpha}
    \right)
  \right] \\
  &= -\left[
    \frac{\cos\gamma}{\tan\beta} 
    \left(
      \frac{\partial}{\partial \gamma}\frac{\sin\gamma}{\tan\beta} \frac{\partial}{\partial \gamma}
      - \frac{\partial}{\partial \gamma}\cos\gamma \frac{\partial}{\partial \beta}
      - \frac{\partial}{\partial \gamma}\frac{\sin\gamma}{\sin\beta} \frac{\partial}{\partial \alpha}
    \right)
    + \sin\gamma \left(
      \frac{\partial}{\partial \beta} \frac{\sin\gamma}{\tan\beta} \frac{\partial}{\partial \gamma}
      - \cos\gamma \frac{\partial}{\partial \beta}
      - \frac{\partial}{\partial \beta} \frac{\sin\gamma}{\sin\beta} \frac{\partial}{\partial \alpha}
    \right)
    -
    \frac{\sin\gamma}{\tan\beta} \left(
      \frac{\partial}{\partial \gamma} \frac{\cos\gamma}{\tan\beta} \frac{\partial}{\partial \gamma}
      + \frac{\partial}{\partial \gamma} \sin\gamma \frac{\partial}{\partial \beta}
      - \frac{\partial}{\partial \gamma} \frac{\cos\gamma}{\sin\beta} \frac{\partial}{\partial \alpha}
    \right)
    + \cos\gamma \left(
       \frac{\partial}{\partial \beta}\frac{\cos\gamma}{\tan\beta} \frac{\partial}{\partial \gamma}
      + \sin\gamma \frac{\partial}{\partial \beta}
      -  \frac{\partial}{\partial \beta} \frac{\cos\gamma}{\sin\beta} \frac{\partial}{\partial \alpha}
    \right)
  \right] \\
  &= -\left[
    \frac{\cos\gamma}{\tan\beta} 
    \left(
      \frac{\cos\gamma}{\tan\beta} \frac{\partial}{\partial \gamma}
      + \sin\gamma \frac{\partial}{\partial \beta}
      - \frac{\cos\gamma}{\sin\beta} \frac{\partial}{\partial \alpha}
    \right)
    + \sin\gamma \left(
      \frac{\sin\gamma}{-\sin^2\beta} \frac{\partial}{\partial \gamma}
      - \frac{\sin\gamma\cos\beta}{-\sin^2\beta} \frac{\partial}{\partial \alpha}
    \right)
    -
    \frac{\sin\gamma}{\tan\beta} \left(
      \frac{-\sin\gamma}{\tan\beta} \frac{\partial}{\partial \gamma}
      + \cos\gamma \frac{\partial}{\partial \beta}
      - \frac{-\sin\gamma}{\sin\beta} \frac{\partial}{\partial \alpha}
    \right)
    + \cos\gamma \left(
       \frac{\cos\gamma}{-\sin^2\beta} \frac{\partial}{\partial \gamma}
      + \sin\gamma \frac{\partial}{\partial \beta}
      -  \frac{\cos\gamma\cos\beta}{-\sin^2\beta} \frac{\partial}{\partial \alpha}
    \right)
  \right] \\
  &= -\left[
    \left(
      + \sin\gamma\cos\gamma \frac{\partial}{\partial \beta}
      - \frac{\cos^2\gamma\cos\beta}{\sin^2\beta} \frac{\partial}{\partial \alpha}
      + \frac{-\sin^2\gamma\cos\beta}{\sin^2\beta} \frac{\partial}{\partial \alpha}
      +  \frac{\cos\beta}{\sin^2\beta} \frac{\partial}{\partial \alpha}
      + \frac{\cos^2\gamma}{\tan^2\beta} \frac{\partial}{\partial \gamma}
      + \frac{\sin^2\gamma}{-\sin^2\beta} \frac{\partial}{\partial \gamma}
      - \frac{-\sin^2\gamma}{\tan^2\beta} \frac{\partial}{\partial \gamma}
      + \frac{\cos^2\gamma}{-\sin^2\beta} \frac{\partial}{\partial \gamma}
    \right)
  \right]
\end{align}
```

