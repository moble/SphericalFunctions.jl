md"""
# Varshalovich et al. (1988)

!!! info "Summary"
    Varshalovich et al.'s Euler angles, Wigner ``D`` and ``d`` functions, and lab-fixed
    angular-momentum operators agree with the definitions used in the `SphericalFunctions`
    package: their ``D^J_{MM'}(α, β, γ)`` is our ``𝔇^{(J)}_{m',m}(α, β, γ)`` with ``(M, M')
    = (m', m)``.  Their *body-fixed* operators ``\hat{J}'`` are related to our ``R``
    operators by ``\hat{J}'_x = -R_x``, ``\hat{J}'_y = R_y``, ``\hat{J}'_z = -R_z``.  Their
    closed-form ``d`` function and tables are also valid for half-integer ``J``, where they
    agree with the complex conjugate of the half-integer ``𝔇`` of [Boyle (2016)](@ref
    "Boyle (2016)"), providing an independent reference for half-integer indices.

[Varshalovich_1988](@citet) is the encyclopedic reference on the quantum theory of angular
momentum, and has a fairly decent comparison of definitions related to the rotation matrix by
previous authors.  Varshalovich et al. define their Euler angles (scheme B, page 22) in the
same way we do, except that they specify that this describes the rotation *of the coordinate
system*.

## Angular-momentum operators and the rotation operator

Varshalovich et al. define the ``\hat{J}`` operators as follows:

> In quantum mechanics the total angular momentum operator ``\hat{J}`` is defined as an
> operator which generates transformations of wave functions (state vectors) and quantum
> operators under infinitesimal rotations of the coordinate system (see Eqs. 2.1.(1) and
> 2.1.(2)).
>
> A transformation of an arbitrary wave function ``\Psi`` under rotation of the coordinate
> system through an infinitesimal angle ``δ \omega`` about an axis ``𝐧`` may be written as
> ```math
> \Psi \to \Psi' = \left(1 - i δ \omega 𝐧 \cdot \hat{J} \right)\Psi,
> ```
> where ``\hat{J}`` is the total angular momentum operator.

Eq. 1.4.(31) defines the operator
```math
\hat{D}(α, β, γ)
=
e^{-iα \hat{J}_z}
e^{-iβ \hat{J}_y}
e^{-iγ \hat{J}_z},
```
which is [our ``U(𝐑_{α,β,γ})``](@ref summary_wigner_D).  Eq. 4.1.(1) defines the Wigner
D-functions according to
```math
\langle J M | \hat{D}(α, β, γ) | J' M' \rangle
=
δ_{J J'} D^J_{M M'}(α, β, γ),
```
and Eq. 4.3.(1) states
```math
D^J_{M M'}(α, β, γ)
=
e^{-i M α}
d^J_{M M'}(β)
e^{-i M' γ},
```
with the ``d`` function given in Eq. 4.3.1(2) as
```math
d^J_{MM'}(β)
=
(-1)^{J-M'}
\sqrt{(J+M)!\,(J-M)!\,(J+M')!\,(J-M')!}
\sum_k (-1)^k
\frac{\left(\cos\frac{β}{2}\right)^{M+M'+2k}
      \left(\sin\frac{β}{2}\right)^{2J-M-M'-2k}}
     {k!\,(J-M-k)!\,(J-M'-k)!\,(M+M'+k)!}.
```
Note that Varshalovich et al. label the indices ``M`` and ``M'``, in the opposite order to
our ``m'`` and ``m``; with ``(M, M') = (m', m)`` these are precisely [our
definitions](@ref summary_wigner_D), and we expect exact agreement.  Table 4.3 (page 119)
gives the spin-``1/2`` matrix explicitly:
```math
\begin{aligned}
D^{1/2}_{1/2, 1/2} &= e^{-iα/2} \cos\tfrac{β}{2}\, e^{-iγ/2}, &
D^{1/2}_{1/2, -1/2} &= -e^{-iα/2} \sin\tfrac{β}{2}\, e^{iγ/2}, \\
D^{1/2}_{-1/2, 1/2} &= e^{iα/2} \sin\tfrac{β}{2}\, e^{-iγ/2}, &
D^{1/2}_{-1/2, -1/2} &= e^{iα/2} \cos\tfrac{β}{2}\, e^{iγ/2}.
\end{aligned}
```

## Half-integer indices

Sec. 4.8.2 (page 92) relates elements with half-integer indices to those with the
neighboring integer indices, by way of the Clebsch–Gordan series for the product of a
spin-``1/2`` matrix with an integer-``J`` matrix.  Specifically, Eqs. 4.8.2(14) and (15) read
```math
\begin{aligned}
D^J_{M M'}
&=
\sqrt{\frac{J-M}{J-M'}}\, \cos\tfrac{β}{2}\, e^{i(α+γ)/2}\, D^{J-1/2}_{M+1/2, M'+1/2}
-
\sqrt{\frac{J+M}{J-M'}}\, \sin\tfrac{β}{2}\, e^{-i(α-γ)/2}\, D^{J-1/2}_{M-1/2, M'+1/2},
\qquad (M' \neq J), \\
D^J_{M M'}
&=
\sqrt{\frac{J-M}{J+M'}}\, \sin\tfrac{β}{2}\, e^{i(α-γ)/2}\, D^{J-1/2}_{M+1/2, M'-1/2}
+
\sqrt{\frac{J+M}{J+M'}}\, \cos\tfrac{β}{2}\, e^{-i(α+γ)/2}\, D^{J-1/2}_{M-1/2, M'-1/2},
\qquad (M' \neq -J),
\end{aligned}
```
where all ``D`` on the right-hand sides have the same Euler angles as the left.  (The forms
given here are the ones that we have verified numerically against the closed-form
expression; earlier transcriptions of these equations in the author's notes contained sign
errors.)  Tables 4.3–4.12 list the ``d^J_{MM'}`` explicitly for ``J \leq 9/2``; we
transcribe below the entries with ``M \geq 1/2``, and obtain the rows with ``M < 0`` from
``d^J_{MM'} = (-1)^{M-M'} d^J_{-M,-M'}``.  The closed form 4.3.1(2), the recursions,
and the tables are all valid for half-integer ``J``, which makes Varshalovich et al. an
independent reference for the half-integer matrices computed by the algorithm of [Boyle
(2016)](@ref "Boyle (2016)"), whose matrices are the complex conjugates of ours.

Finally, the spin-weighted spherical harmonics of half-integer spin weight — defined on
[our summary page](@ref summary_swsh) by ``{}_sY_{ℓ,m} = (-1)^s \sqrt{(2ℓ+1)/4π}\,
\overline{𝔇_{m,-s}}`` with ``(-1)^s \equiv e^{iπs}`` — can be built from Varshalovich's
``D``, and we use them to check the anchor ``{}_sY_{ℓ,-s}(𝟏) = i^{2s}\sqrt{(2ℓ+1)/4π}`` and
the conjugation relation ``\overline{{}_sY_{ℓ,m}} = (-1)^{m+s}\, {}_{-s}Y_{ℓ,-m}`` for
half-integer indices.

## Body-fixed operators

Varshalovich et al. distinguish in Sec. 1.1.3 between *covariant* and *contravariant*
spherical components and the corresponding basis vectors, which they define as
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
Then, in Sec. 4.2 they define ``\hat{𝐉}`` as the operator of angular momentum of the rigid
symmetric top.  They then give in Eq. (6) the "covariant spherical coordinates of
``\hat{𝐉}`` in the non-rotating (lab-fixed) system" as
```math
\begin{gathered}
  \hat{J}_{\pm 1} = \frac{i}{\sqrt{2}} e^{\pm i α} \left[
    \mp \cot β \frac{\partial}{\partial α}
    + i \frac{\partial}{\partial β}
    \pm \frac{1}{\sin β} \frac{\partial}{\partial γ}
  \right] \\
  \hat{J}_0 = - i \frac{\partial}{\partial α},
\end{gathered}
```
and in Eq. (7) the "contravariant components of ``\hat{𝐉}`` in the rotating (body-fixed)
system" as
```math
\begin{gathered}
  \hat{J}'^{\pm 1} = \frac{i}{\sqrt{2}} e^{\mp i γ} \left[
    \pm \cot β \frac{\partial}{\partial γ}
    + i \frac{\partial}{\partial β}
    \mp \frac{1}{\sin β} \frac{\partial}{\partial α}
  \right] \\
  \hat{J}'^0 = - i \frac{\partial}{\partial γ}.
\end{gathered}
```
(Note the prime in the last two equations.)  We can expand these in Cartesian components to
compare to our expressions.  First the covariant components:
```math
\begin{aligned}
  \hat{J}_{x}
  &= -\frac{1}{\sqrt{2}} \left( \hat{J}_{+1} - \hat{J}_{-1} \right)
  = i\left[
      \frac{\cos α}{\tan β} \frac{\partial}{\partial α}
      + \sin α \frac{\partial}{\partial β}
      - \frac{\cos α}{\sin β} \frac{\partial}{\partial γ}
  \right], \\
  \hat{J}_{y}
  &= -\frac{1}{i\sqrt{2}} \left( \hat{J}_{+1} + \hat{J}_{-1} \right)
  = i \left[
      \frac{\sin α}{\tan β} \frac{\partial}{\partial α}
      - \cos α \frac{\partial}{\partial β}
      - \frac{\sin α}{\sin β} \frac{\partial}{\partial γ}
  \right], \\
  \hat{J}_{z}
  &= \hat{J}_{0}
  = -i \frac{\partial}{\partial α}.
\end{aligned}
```
We can compare these to the [full expressions on ``𝕊³``](@ref euler_full_S3), and find that
they are precisely our ``L_x``, ``L_y``, and ``L_z``.  Next, the contravariant components:
```math
\begin{aligned}
  \hat{J}'_{x}
  &= -\frac{1}{\sqrt{2}} \left( \hat{J}'^{+1} - \hat{J}'^{-1} \right)
  = -i \left(
      \frac{\cos γ}{\tan β} \frac{\partial}{\partial γ}
      + \sin γ \frac{\partial}{\partial β}
      - \frac{\cos γ}{\sin β} \frac{\partial}{\partial α}
  \right), \\
  \hat{J}'_{y}
  &= \frac{1}{i\sqrt{2}} \left( \hat{J}'^{+1} + \hat{J}'^{-1} \right)
  = -i \left(
      \frac{\sin γ}{\tan β} \frac{\partial}{\partial γ}
      - \cos γ \frac{\partial}{\partial β}
      - \frac{\sin γ}{\sin β} \frac{\partial}{\partial α}
  \right), \\
  \hat{J}'_{z}
  &= \hat{J}'^{0}
  = -i \frac{\partial}{\partial γ}.
\end{aligned}
```
Comparing with [our ``R`` operators](@ref euler_R_S3), we have ``\hat{J}'_y = R_y``, but
``\hat{J}'_x = -R_x`` and ``\hat{J}'_z = -R_z``.  It's very easy to check that
``[\hat{J}'_{z}, \hat{J}'_{x}] = i \hat{J}'_{y}``, as expected from the general expression in
their Eq. (12), so these expressions are — at least — consistent with the claims of
Varshalovich et al.; the sign difference arises from their treatment of the body-fixed frame.
Below we verify all six relations numerically, by applying both sets of operators to the
``D`` functions with automatic differentiation.

## Implementing formulas

We begin by writing code that implements the formulas from Varshalovich et al.  Because the
closed form and the tables are also needed on the [Boyle (2016)](@ref "Boyle (2016)") page,
we define them in a test module that both pages can use.
"""

# TODO: Confirm the section/page of the quoted definition of Ĵ (Sec. 1.4?) and of Eq. (12) for the commutators.  #src
using TestItems: @testmodule, @testitem  #hide
@testmodule Varshalovich begin  #hide

import ForwardDiff

const 𝒾 = im
#+

# Factorials of integers and of integer-valued rationals (half-integer arithmetic produces
# the latter), computed exactly, in the postfix form `(n)❗` used throughout these pages:
struct Factorial end
Base.:*(n::Integer, ::Factorial) = factorial(big(n))
Base.:*(n::Rational, ::Factorial) = factorial(big(Int(n)))
const ❗ = Factorial()
#+

# Eq. 4.3.1(2).  The sum runs over all ``k`` for which the factorials have non-negative
# arguments, ``\max(0, -(M+M')) \leq k \leq \min(J-M, J-M')``.  The formula is valid for
# integer and half-integer ``J``; note that ``J - M'`` and ``M + M' + 2k`` are always
# integers.
function d(J, M, M′, β::T) where {T<:Real}
    if abs(M) > J || abs(M′) > J
        return zero(T)  # convenient when applying Eqs. 4.8.2(14)-(15) at the edges
    end
    (-1)^Int(J-M′) * √T((J+M)❗ * (J-M)❗ * (J+M′)❗ * (J-M′)❗) *
    sum(
        (-1)^k * cos(β/2)^Int(M+M′+2k) * sin(β/2)^Int(2J-M-M′-2k)
        / T((k)❗ * (J-M-k)❗ * (J-M′-k)❗ * (M+M′+k)❗)
        for k ∈ Int(max(0, -(M+M′))):Int(min(J-M, J-M′));
        init=zero(T)
    )
end
#+

# Eq. 4.3.(1):
function D(J, M, M′, α, β, γ)
    exp(-𝒾 * M * α) * d(J, M, M′, β) * exp(-𝒾 * M′ * γ)
end
#+

# Table 4.3, the spin-1/2 matrix:
function D½(M, M′, α, β, γ)
    if (M, M′) == (1//2, 1//2)
        exp(-𝒾*α/2) * cos(β/2) * exp(-𝒾*γ/2)
    elseif (M, M′) == (1//2, -1//2)
        -exp(-𝒾*α/2) * sin(β/2) * exp(𝒾*γ/2)
    elseif (M, M′) == (-1//2, 1//2)
        exp(𝒾*α/2) * sin(β/2) * exp(-𝒾*γ/2)
    elseif (M, M′) == (-1//2, -1//2)
        exp(𝒾*α/2) * cos(β/2) * exp(𝒾*γ/2)
    end
end
#+

# Eqs. 4.8.2(14) and (15), expressing a half-integer-``J`` element in terms of elements with
# ``J - 1/2``:
function D_recursion(J, M, M′, α, β, γ)
    if M′ ≠ J  # Eq. 4.8.2(14)
        (
            √((J-M)/(J-M′)) * cos(β/2) * exp(𝒾*(α+γ)/2) * D(J-1//2, M+1//2, M′+1//2, α, β, γ)
            -
            √((J+M)/(J-M′)) * sin(β/2) * exp(-𝒾*(α-γ)/2) * D(J-1//2, M-1//2, M′+1//2, α, β, γ)
        )
    else  # Eq. 4.8.2(15)
        (
            √((J-M)/(J+M′)) * sin(β/2) * exp(𝒾*(α-γ)/2) * D(J-1//2, M+1//2, M′-1//2, α, β, γ)
            +
            √((J+M)/(J+M′)) * cos(β/2) * exp(-𝒾*(α+γ)/2) * D(J-1//2, M-1//2, M′-1//2, α, β, γ)
        )
    end
end
#+

# Explicit values for the half-integer ``d`` functions, as given in Tables 4.3—4.12 for ``J
# ∈ [1/2, 9/2]``.  The tables list only the rows ``M ≥ 1/2``, and for each such row only the
# entries not obtainable from the ones already given by the symmetries; we return `nothing`
# for entries the tables omit, and use ``d^J_{MM'} = (-1)^{M-M'} d^J_{-M,-M'}`` for ``M <
# 0``.
function d_½_explicit(J::Rational{Int}, M::Rational{Int}, M′::Rational{Int}, β::T) where T
    if denominator(J) != 2 || denominator(M) != 2 || denominator(M′) != 2
        error("Only half-integer J, M, M′ are supported")
    end
    if J < 1//2 || J > 9//2
        error("Only J = 1/2, 3/2, 5/2, 7/2, 9/2 are supported")
    end
    if abs(M) > J || abs(M′) > J
        error("abs(M) and abs(M′) must be ≤ J")
    end
    if M < 0
        r = d_½_explicit(J, -M, -M′, β)
        return r === nothing ? nothing : (-1)^Int(M-M′) * r
    else
        let √ = (x -> √T(x))
            if (J, M, M′) == (1//2, 1//2,-1//2)
                -sin(β/2)
            elseif (J, M, M′) == (1//2, 1//2, 1//2)
                cos(β/2)

            elseif (J, M, M′) == (3//2, 1//2,-3//2)
                √3 * sin(β/2)^2 * cos(β/2)
            elseif (J, M, M′) == (3//2, 1//2,-1//2)
                sin(β/2) * (3 * sin(β/2)^2 - 2)
            elseif (J, M, M′) == (3//2, 1//2, 1//2)
                cos(β/2) * (3 * cos(β/2)^2 - 2)
            elseif (J, M, M′) == (3//2, 1//2, 3//2)
                √3 * sin(β/2) * cos(β/2)^2
            elseif (J, M, M′) == (3//2, 3//2,-3//2)
                -sin(β/2)^3
            elseif (J, M, M′) == (3//2, 3//2,-1//2)
                √3 * sin(β/2)^2 * cos(β/2)
            elseif (J, M, M′) == (3//2, 3//2, 1//2)
                -√3 * sin(β/2) * cos(β/2)^2
            elseif (J, M, M′) == (3//2, 3//2, 3//2)
                cos(β/2)^3

            elseif (J, M, M′) == (5//2, 5//2, 5//2)
                cos(β/2)^5
            elseif (J, M, M′) == (5//2, 5//2, 3//2)
                -√5 * sin(β/2) * cos(β/2)^4
            elseif (J, M, M′) == (5//2, 5//2, 1//2)
                √10 * sin(β/2)^2 * cos(β/2)^3
            elseif (J, M, M′) == (5//2, 5//2,-1//2)
                -√10 * sin(β/2)^3 * cos(β/2)^2
            elseif (J, M, M′) == (5//2, 5//2,-3//2)
                √5 * sin(β/2)^4 * cos(β/2)
            elseif (J, M, M′) == (5//2, 5//2,-5//2)
                -sin(β/2)^5
            elseif (J, M, M′) == (5//2, 3//2, 3//2)
                cos(β/2)^3 * (1 - 5 * sin(β/2)^2)
            elseif (J, M, M′) == (5//2, 3//2, 1//2)
                -√2 * sin(β/2) * cos(β/2)^2 * (2 - 5 * sin(β/2)^2)
            elseif (J, M, M′) == (5//2, 3//2,-1//2)
                -√2 * sin(β/2)^2 * cos(β/2) * (2 - 5 * cos(β/2)^2)
            elseif (J, M, M′) == (5//2, 3//2,-3//2)
                sin(β/2)^3 * (1 - 5 * cos(β/2)^2)
            elseif (J, M, M′) == (5//2, 1//2, 1//2)
                cos(β/2) * (3 - 12 * cos(β/2)^2 + 10 * cos(β/2)^4)
            elseif (J, M, M′) == (5//2, 1//2,-1//2)
                -sin(β/2) * (3 - 12 * sin(β/2)^2 + 10 * sin(β/2)^4)

            elseif (J, M, M′) == (7//2, 7//2, 7//2)
                cos(β/2)^7
            elseif (J, M, M′) == (7//2, 7//2, 5//2)
                -√7 * cos(β/2)^6 * sin(β/2)
            elseif (J, M, M′) == (7//2, 7//2, 3//2)
                √21 * cos(β/2)^5 * sin(β/2)^2
            elseif (J, M, M′) == (7//2, 7//2, 1//2)
                -√35 * cos(β/2)^4 * sin(β/2)^3
            elseif (J, M, M′) == (7//2, 7//2,-1//2)
                √35 * cos(β/2)^3 * sin(β/2)^4
            elseif (J, M, M′) == (7//2, 7//2,-3//2)
                -√21 * cos(β/2)^2 * sin(β/2)^5
            elseif (J, M, M′) == (7//2, 7//2,-5//2)
                √7 * cos(β/2) * sin(β/2)^6
            elseif (J, M, M′) == (7//2, 7//2,-7//2)
                -sin(β/2)^7
            elseif (J, M, M′) == (7//2, 5//2, 5//2)
                cos(β/2)^5 * (1 - 7 * sin(β/2)^2)
            elseif (J, M, M′) == (7//2, 5//2, 3//2)
                -√3 * cos(β/2)^4 * sin(β/2) * (2 - 7 * sin(β/2)^2)
            elseif (J, M, M′) == (7//2, 5//2, 1//2)
                √5 * cos(β/2)^3 * sin(β/2)^2 * (3 - 7 * sin(β/2)^2)
            elseif (J, M, M′) == (7//2, 5//2,-1//2)
                √5 * cos(β/2)^2 * sin(β/2)^3 * (3 - 7 * cos(β/2)^2)
            elseif (J, M, M′) == (7//2, 5//2,-3//2)
                -√3 * cos(β/2) * sin(β/2)^4 * (2 - 7 * cos(β/2)^2)
            elseif (J, M, M′) == (7//2, 5//2,-5//2)
                sin(β/2)^5 * (1 - 7 * cos(β/2)^2)
            elseif (J, M, M′) == (7//2, 3//2, 3//2)
                cos(β/2)^3 * (10 - 30 * cos(β/2)^2 + 21 * cos(β/2)^4)
            elseif (J, M, M′) == (7//2, 3//2, 1//2)
                -√15 * cos(β/2)^2 * sin(β/2) * (2 - 8 * cos(β/2)^2 + 7 * cos(β/2)^4)
            elseif (J, M, M′) == (7//2, 3//2,-1//2)
                √15 * cos(β/2) * sin(β/2)^2 * (2 - 8 * sin(β/2)^2 + 7 * sin(β/2)^4)
            elseif (J, M, M′) == (7//2, 3//2,-3//2)
                -sin(β/2)^3 * (10 - 30 * sin(β/2)^2 + 21 * sin(β/2)^4)
            elseif (J, M, M′) == (7//2, 1//2, 1//2)
                -cos(β/2) * (4 - 30 * cos(β/2)^2 + 60 * cos(β/2)^4 - 35 * cos(β/2)^6)
            elseif (J, M, M′) == (7//2, 1//2,-1//2)
                -sin(β/2) * (4 - 30 * sin(β/2)^2 + 60 * sin(β/2)^4 - 35 * sin(β/2)^6)

            elseif (J, M, M′) == (9//2, 9//2, 9//2)
                cos(β/2)^9
            elseif (J, M, M′) == (9//2, 9//2, 7//2)
                -3 * cos(β/2)^8 * sin(β/2)
            elseif (J, M, M′) == (9//2, 9//2, 5//2)
                6 * cos(β/2)^7 * sin(β/2)^2
            elseif (J, M, M′) == (9//2, 9//2, 3//2)
                -2 * √21 * cos(β/2)^6 * sin(β/2)^3
            elseif (J, M, M′) == (9//2, 9//2, 1//2)
                3 * √14 * cos(β/2)^5 * sin(β/2)^4
            elseif (J, M, M′) == (9//2, 9//2,-1//2)
                -3 * √14 * cos(β/2)^4 * sin(β/2)^5
            elseif (J, M, M′) == (9//2, 9//2,-3//2)
                2 * √21 * cos(β/2)^3 * sin(β/2)^6
            elseif (J, M, M′) == (9//2, 9//2,-5//2)
                -6 * cos(β/2)^2 * sin(β/2)^7
            elseif (J, M, M′) == (9//2, 9//2,-7//2)
                3 * cos(β/2) * sin(β/2)^8
            elseif (J, M, M′) == (9//2, 9//2,-9//2)
                -sin(β/2)^9
            elseif (J, M, M′) == (9//2, 7//2, 7//2)
                cos(β/2)^7 * (1 - 9 * sin(β/2)^2)
            elseif (J, M, M′) == (9//2, 7//2, 5//2)
                -2 * cos(β/2)^6 * sin(β/2) * (2 - 9 * sin(β/2)^2)
            elseif (J, M, M′) == (9//2, 7//2, 3//2)
                2 * √21 * cos(β/2)^5 * sin(β/2)^2 * (1 - 3 * sin(β/2)^2)
            elseif (J, M, M′) == (9//2, 7//2, 1//2)
                -√14 * cos(β/2)^4 * sin(β/2)^3 * (4 - 9 * sin(β/2)^2)
            elseif (J, M, M′) == (9//2, 7//2,-1//2)
                -√14 * cos(β/2)^3 * sin(β/2)^4 * (4 - 9 * cos(β/2)^2)
            elseif (J, M, M′) == (9//2, 7//2,-3//2)
                2 * √21 * cos(β/2)^2 * sin(β/2)^5 * (1 - 3 * cos(β/2)^2)
            elseif (J, M, M′) == (9//2, 7//2,-5//2)
                -2 * cos(β/2) * sin(β/2)^6 * (2 - 9 * cos(β/2)^2)
            elseif (J, M, M′) == (9//2, 7//2,-7//2)
                sin(β/2)^7 * (1 - 9 * cos(β/2)^2)
            elseif (J, M, M′) == (9//2, 5//2, 5//2)
                cos(β/2)^5 * (21 - 56 * cos(β/2)^2 + 36 * cos(β/2)^4)
            elseif (J, M, M′) == (9//2, 5//2, 3//2)
                -√21 * cos(β/2)^4 * sin(β/2) * (5 - 16 * cos(β/2)^2 + 12 * cos(β/2)^4)
            elseif (J, M, M′) == (9//2, 5//2, 1//2)
                √14 * cos(β/2)^3 * sin(β/2)^2 * (5 - 20 * cos(β/2)^2 + 18 * cos(β/2)^4)
            elseif (J, M, M′) == (9//2, 5//2,-1//2)
                -√14 * cos(β/2)^2 * sin(β/2)^3 * (5 - 20 * sin(β/2)^2 + 18 * sin(β/2)^4)
            elseif (J, M, M′) == (9//2, 5//2,-3//2)
                √21 * cos(β/2) * sin(β/2)^4 * (5 - 16 * sin(β/2)^2 + 12 * sin(β/2)^4)
            elseif (J, M, M′) == (9//2, 5//2,-5//2)
                -sin(β/2)^5 * (21 - 56 * sin(β/2)^2 + 36 * sin(β/2)^4)
            elseif (J, M, M′) == (9//2, 3//2, 3//2)
                -cos(β/2)^3 * (20 - 105 * cos(β/2)^2 + 168 * cos(β/2)^4 - 84 * cos(β/2)^6)
            elseif (J, M, M′) == (9//2, 3//2, 1//2)
                √6 * cos(β/2)^2 * sin(β/2) * (5 - 35 * cos(β/2)^2 + 70 * cos(β/2)^4 - 42 * cos(β/2)^6)
            elseif (J, M, M′) == (9//2, 3//2,-1//2)
                √6 * cos(β/2) * sin(β/2)^2 * (5 - 35 * sin(β/2)^2 + 70 * sin(β/2)^4 - 42 * sin(β/2)^6)
            elseif (J, M, M′) == (9//2, 3//2,-3//2)
                -sin(β/2)^3 * (20 - 105 * sin(β/2)^2 + 168 * sin(β/2)^4 - 84 * sin(β/2)^6)
            elseif (J, M, M′) == (9//2, 1//2, 1//2)
                cos(β/2) * (5 - 60 * cos(β/2)^2 + 210 * cos(β/2)^4 - 280 * cos(β/2)^6 + 126 * cos(β/2)^8)
            elseif (J, M, M′) == (9//2, 1//2,-1//2)
                -sin(β/2) * (5 - 60 * sin(β/2)^2 + 210 * sin(β/2)^4 - 280 * sin(β/2)^6 + 126 * sin(β/2)^8)
            end
        end
    end
end
#+

# The angular-momentum operators of Sec. 4.2, Eqs. (6) and (7), in their Cartesian forms
# as expanded above.  Each takes a function `f(α, β, γ)` and returns a new function, with the
# derivatives evaluated by forward-mode automatic differentiation.
∂α(f) = (α, β, γ) -> ForwardDiff.derivative(α′ -> f(α′, β, γ), α)
∂β(f) = (α, β, γ) -> ForwardDiff.derivative(β′ -> f(α, β′, γ), β)
∂γ(f) = (α, β, γ) -> ForwardDiff.derivative(γ′ -> f(α, β, γ′), γ)
Ĵx(f) = (α, β, γ) -> 𝒾 * (cos(α)/tan(β) * ∂α(f)(α, β, γ) + sin(α) * ∂β(f)(α, β, γ) - cos(α)/sin(β) * ∂γ(f)(α, β, γ))
Ĵy(f) = (α, β, γ) -> 𝒾 * (sin(α)/tan(β) * ∂α(f)(α, β, γ) - cos(α) * ∂β(f)(α, β, γ) - sin(α)/sin(β) * ∂γ(f)(α, β, γ))
Ĵz(f) = (α, β, γ) -> -𝒾 * ∂α(f)(α, β, γ)
Ĵ′x(f) = (α, β, γ) -> -𝒾 * (cos(γ)/tan(β) * ∂γ(f)(α, β, γ) + sin(γ) * ∂β(f)(α, β, γ) - cos(γ)/sin(β) * ∂α(f)(α, β, γ))
Ĵ′y(f) = (α, β, γ) -> -𝒾 * (sin(γ)/tan(β) * ∂γ(f)(α, β, γ) - cos(γ) * ∂β(f)(α, β, γ) - sin(γ)/sin(β) * ∂α(f)(α, β, γ))
Ĵ′z(f) = (α, β, γ) -> -𝒾 * ∂γ(f)(α, β, γ)
#+

end  #hide

md"""
## Tests

We can now test the functions against the equivalent functions from the `SphericalFunctions`
package.  For the operator comparison we also need our own ``L`` and ``R`` operators in
Euler angles, which we transcribe from the [summary page](@ref summary_L_R_euler).
"""

@testitem "Varshalovich conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities, Varshalovich, Boyle2016] begin  #hide
import ForwardDiff
import Quaternionic: from_euler_angles
const 𝒾 = im
#+

# Our operators, from the summary page, in the same form as Varshalovich's above:
∂α(f) = (α, β, γ) -> ForwardDiff.derivative(α′ -> f(α′, β, γ), α)
∂β(f) = (α, β, γ) -> ForwardDiff.derivative(β′ -> f(α, β′, γ), β)
∂γ(f) = (α, β, γ) -> ForwardDiff.derivative(γ′ -> f(α, β, γ′), γ)
Lx(f) = (α, β, γ) -> 𝒾 * (cos(α)/tan(β) * ∂α(f)(α, β, γ) + sin(α) * ∂β(f)(α, β, γ) - cos(α)/sin(β) * ∂γ(f)(α, β, γ))
Ly(f) = (α, β, γ) -> 𝒾 * (sin(α)/tan(β) * ∂α(f)(α, β, γ) - cos(α) * ∂β(f)(α, β, γ) - sin(α)/sin(β) * ∂γ(f)(α, β, γ))
Lz(f) = (α, β, γ) -> -𝒾 * ∂α(f)(α, β, γ)
Rx(f) = (α, β, γ) -> 𝒾 * (-cos(γ)/sin(β) * ∂α(f)(α, β, γ) + sin(γ) * ∂β(f)(α, β, γ) + cos(γ)/tan(β) * ∂γ(f)(α, β, γ))
Ry(f) = (α, β, γ) -> 𝒾 * (sin(γ)/sin(β) * ∂α(f)(α, β, γ) + cos(γ) * ∂β(f)(α, β, γ) - sin(γ)/tan(β) * ∂γ(f)(α, β, γ))
Rz(f) = (α, β, γ) -> 𝒾 * ∂γ(f)(α, β, γ)
#+

# We will need to test approximate floating-point equality, so we set absolute and relative
# tolerances (respectively) in terms of the machine epsilon:
ϵₐ = 100eps()
ϵᵣ = 1000eps()
#+

# The closed-form expression for ``d`` overflows `Float64` for ``J \gtrsim 8``, so we test
# up to
Jₘₐₓ = 5
#+
# for integers and ``15/2`` for half-integers.  The closed form is slow (it uses exact
# integer arithmetic for the factorials), so we use a modest grid of Euler angles for the
# integer tests and a smaller one for the half-integer and operator tests:
αβγs = αβγrange(Float64, 5)
αβγs_small = αβγrange(Float64, 1)
#+

# First, the integer case: with ``(M, M') = (m', m)``, Varshalovich's ``d`` and ``D`` are
# ours:
for β ∈ βrange()
    for (J, m′, m) ∈ ℓm′mrange(Jₘₐₓ)
        @test Varshalovich.d(J, m′, m, β) ≈ ConventionsUtilities.d(J, m′, m, β) atol=ϵₐ rtol=ϵᵣ
    end
end
for (α, β, γ) ∈ αβγs
    for (J, m′, m) ∈ ℓm′mrange(Jₘₐₓ)
        @test Varshalovich.D(J, m′, m, α, β, γ) ≈ ConventionsUtilities.D(J, m′, m, α, β, γ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Table 4.3 agrees with the closed form:
for (α, β, γ) ∈ αβγs
    for M ∈ (-1//2, 1//2), M′ ∈ (-1//2, 1//2)
        @test Varshalovich.D½(M, M′, α, β, γ) ≈ Varshalovich.D(1//2, M, M′, α, β, γ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# The recursion relations 4.8.2(14)–(15) reproduce the closed form for half-integer ``J``:
for (α, β, γ) ∈ αβγs_small
    for J ∈ 1//2:15//2, M ∈ -J:J, M′ ∈ -J:J
        @test Varshalovich.D_recursion(J, M, M′, α, β, γ) ≈ Varshalovich.D(J, M, M′, α, β, γ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Every entry transcribed from Tables 4.3–4.12 agrees with the closed form (entries the
# tables omit return `nothing` and are skipped):
for β ∈ βrange()
    for J ∈ 1//2:9//2, M ∈ -J:J, M′ ∈ -J:J
        dₜ = Varshalovich.d_½_explicit(J, M, M′, β)
        dₜ === nothing && continue
        @test dₜ ≈ Varshalovich.d(J, M, M′, β) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# For half-integer ``J``, Varshalovich's closed form and the quaternion algorithm of Boyle
# (2016) are independent references; they agree up to the complex conjugation that
# distinguishes the 2016 convention from the present one, for ``J \leq 15/2``:
for (α, β, γ) ∈ αβγs_small
    R = from_euler_angles(α, β, γ)
    for J ∈ 1//2:15//2, M ∈ -J:J, M′ ∈ -J:J
        @test Varshalovich.D(J, M, M′, α, β, γ) ≈ conj(Boyle2016.WignerDElement(R, J, M, M′)) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Half-integer spin-weighted spherical harmonics built from Varshalovich's ``D``, following
# our definition with ``(-1)^s ≡ e^{iπs}``, satisfy the anchor condition and the conjugation
# relation from the summary page:
ₛYₗₘ(s, ℓ, m, α, β, γ) = exp(𝒾 * π * s) * √((2ℓ+1) / (4π)) * conj(Varshalovich.D(ℓ, m, -s, α, β, γ))
for ℓ ∈ 1//2:7//2, s ∈ -ℓ:ℓ
    @test ₛYₗₘ(s, ℓ, -s, 0.0, 0.0, 0.0) ≈ (1.0𝒾)^Int(2s) * √((2ℓ+1) / (4π)) atol=ϵₐ rtol=ϵᵣ
end
for (α, β, γ) ∈ αβγs_small
    for ℓ ∈ 1//2:7//2, s ∈ -ℓ:ℓ, m ∈ -ℓ:ℓ
        @test conj(ₛYₗₘ(s, ℓ, m, α, β, γ)) ≈ (-1)^Int(m+s) * ₛYₗₘ(-s, ℓ, -m, α, β, γ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Finally, the operators.  Applying Varshalovich's lab-fixed operators and our ``L``
# operators to the ``D`` functions gives identical results; applying the body-fixed operators
# gives ``\hat{J}'_x = -R_x``, ``\hat{J}'_y = R_y``, and ``\hat{J}'_z = -R_z``.  The
# operators involve ``1/\sin β``, so we avoid the poles; and because each application
# evaluates the closed form many times under automatic differentiation, we use just a few
# generic Euler-angle triples.
for (α, β, γ) ∈ [(0.7, 1.1, 2.3), (2.9, 0.4, 5.1), (4.0, 2.2, 0.3), (1.3, 2.9, 4.7), (5.5, 1.7, 1.9)]
    for (J, M, M′) ∈ ℓm′mrange(2)
        f(α, β, γ) = Varshalovich.D(J, M, M′, α, β, γ)
        @test Varshalovich.Ĵx(f)(α, β, γ) ≈ Lx(f)(α, β, γ) atol=ϵₐ rtol=ϵᵣ
        @test Varshalovich.Ĵy(f)(α, β, γ) ≈ Ly(f)(α, β, γ) atol=ϵₐ rtol=ϵᵣ
        @test Varshalovich.Ĵz(f)(α, β, γ) ≈ Lz(f)(α, β, γ) atol=ϵₐ rtol=ϵᵣ
        @test Varshalovich.Ĵ′x(f)(α, β, γ) ≈ -Rx(f)(α, β, γ) atol=ϵₐ rtol=ϵᵣ
        @test Varshalovich.Ĵ′y(f)(α, β, γ) ≈ Ry(f)(α, β, γ) atol=ϵₐ rtol=ϵᵣ
        @test Varshalovich.Ĵ′z(f)(α, β, γ) ≈ -Rz(f)(α, β, γ) atol=ϵₐ rtol=ϵᵣ
        ## and, as a consequence, the eigenvalue relations Ĵz D = -M D and Ĵ′z D = -M′ D
        @test Varshalovich.Ĵz(f)(α, β, γ) ≈ -M * f(α, β, γ) atol=ϵₐ rtol=ϵᵣ
        @test Varshalovich.Ĵ′z(f)(α, β, γ) ≈ -M′ * f(α, β, γ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# These successful tests show that Varshalovich et al.'s ``D`` and ``d`` functions agree
# with ours (with the indices read in their order), for integer and half-integer ``J``; that
# their lab-fixed operators are our ``L`` operators; and that their body-fixed operators are
# ``(-R_x, R_y, -R_z)``.

end  #hide
