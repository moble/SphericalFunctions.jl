# [``{}_{s}Y_{ℓ,m}`` and ``𝔇^{(ℓ)}_{m', m}``](@id sYlm_and_Dlmpm)

The [previous page](@ref background_differential_operators) introduced
the left and right Lie derivative operators acting on functions of a
quaternion argument.  Here, we show how the spin-weighted spherical
harmonics arise naturally as eigenfunctions of these operators.

The familiar way of arriving at the standard (scalar, ``s=0``)
spherical harmonics is to consider solutions to the Laplace equation
in three-dimensional space, then to separate variables in spherical
coordinates, and finally to identify the angular part of the solution
as the spherical harmonics.  Since scalar spherical harmonics really
are just functions on the sphere ``𝕊²``, rather than the full group
``\mathrm{Spin}(3)``, they can be sensibly written as functions of
spherical coordinates ``(θ, ϕ)`` alone.  We'll suspend for a moment
what, exactly, the angular-momentum operators mean for functions on
the sphere, and just accept the usual definition of ``L`` as the
relevant operator in this case.  Recall that we can only
simultaneously diagonalize operators (find nontrivial functions that
are eigenfunctions of both operators at the same time) if they
commute.  In this case, the only two commuting operators are ``L^2 =
L_x^2 + L_y^2 + L_z^2`` and any given component of ``L`` —
conventionally chosen to be ``L_z``.  Specifically, the spherical
harmonics are defined to satisfy
```math
\begin{aligned}
L² \left\{ Y_{ℓ,m} \right\}(θ, ϕ)
&= ℓ(ℓ+1) \left\{ Y_{ℓ,m} \right\}(θ, ϕ),
\\
L_z \left\{ Y_{ℓ,m} \right\}(θ, ϕ)
&= m \left\{ Y_{ℓ,m} \right\}(θ, ϕ).
\end{aligned}
```
(As usual, we do not include the factors of ``\hbar`` that would
normally be included in quantum mechanics texts.)

But now, we consider functions defined on the full group
``\mathrm{Spin}(3)``.  In that case, we have both left and right Lie
derivative operators available.  As before, we can choose
eigenfunctions of ``L²`` and ``L_z``.  It turns out that ``[L_g, R_h]
= 0``, which means that we could also simultaneously diagonalize with
respect to any given component of ``R`` — which again is
conventionally chosen to be ``R_z``.  We might also expect to be able
to diagonalize with respect to ``R²``, but that turns out to equal
``L²``, so while it is true that the spin-weighted spherical harmonics
are also eigenvalues of ``R²``, that statement contains no additional
information.  Thus, we select the spin-weighted spherical harmonics as
functions satisfying
```math
\begin{aligned}
L² \left\{ {}_{s}Y_{ℓ,m} \right\}(Q)
&= ℓ(ℓ+1) \left\{ {}_{s}Y_{ℓ,m} \right\}(Q),
\\
L_z \left\{ {}_{s}Y_{ℓ,m} \right\}(Q)
&= m \left\{ {}_{s}Y_{ℓ,m} \right\}(Q),
\\
R_z \left\{ {}_{s}Y_{ℓ,m} \right\}(Q)
&= s \left\{ {}_{s}Y_{ℓ,m} \right\}(Q).
\end{aligned}
```
These conditions only determine the functions up to a normalization
factor, which is conventionally chosen so that the functions are
orthonormal with respect to the natural measure on
``\mathrm{Spin}(3)`` — which we can implement using [extended Euler
angles](@ref Quaternions-and-Euler-angles).  This still leaves an
overall complex phase freedom, which is conventionally fixed by
requiring that ``{}_{s}Y_{ℓ,m}(I)``, where ``I`` is the identity
rotor, be real and nonnegative.

Now, we can return to the idea of functions on ``𝕊²``.  We can
identify each point in ``\mathrm{Spin}(3)`` with a point on ``𝕊²`` by
considering the direction that a given rotor takes the ``z`` axis.
This loses information regarding rotation *about* that point.  A
function on ``\mathrm{Spin}(3)`` can then be converted into a function
on ``𝕊²`` if and only if the function does not depend on the rotation
about that point, which is described by ``R_z``.  Thus, only functions
with ``s=0`` can be considered as functions on ``𝕊²``.  In fact, the
scalar spherical harmonics are precisely the spin-weighted spherical
harmonics with ``s=0``.

If we alter notation slightly and write the spherical harmonics as
functions of a unit vector to a point rather than the spherical
coordinates describing the same point, we can define scalar spherical
harmonics as functions on ``\mathrm{Spin}(3)`` as well:
```math
Y_{ℓ,m}(𝐐) = Y_{ℓ,m}\left(𝐐\, 𝐳\, 𝐐⁻¹\right)
```
where ``𝐳`` is the unit vector in the ``z`` direction.  The
right-hand side is just the usual scalar spherical harmonic evaluated
at the point on ``𝕊²`` corresponding to the rotation of the ``z``
axis by the rotor ``𝐐``.  Now, we can explicitly write
```math
Y_{ℓ,m}(𝐐) = {}_0Y_{ℓ,m}(𝐐);
```
the scalar spherical harmonics are *precisely* the spin-weighted
spherical harmonics with spin weight ``s=0``.

Finally, we can consider the Wigner ``𝔇`` matrices.  These are
usually defined as
```math
𝔇^{(ℓ)}_{m', m}(α, β, γ)
= \big\langle ℓ, m' \big|
  e^{-iL_z α} e^{-iL_y β} e^{-iL_z γ}
  \big| ℓ, m \big\rangle.
```
We can rewrite this as
```math
𝔇^{(ℓ)}_{m', m}(𝐐)
= \int_{\mathrm{Spin}(3)}
  \bar{Y}_{ℓ,m'}(𝐐')\,
  Y_{ℓ,m}\left(𝐐⁻¹ 𝐐'\right)\,
  d𝐐',
```
where the integral is taken over ``\mathrm{Spin}(3)`` with the
appropriate measure to ensure orthonormality.  We can evaluate the
action of the differential operators on this function pretty easily,
noting that the derivative ``d/dϵ`` can pass through the integral sign
by the Leibniz integral rule.  The result is that ``𝔇`` satisfies
```math
\begin{aligned}
L² \left\{ 𝔇^{(ℓ)}_{m', m} \right\}(𝐐)
&= ℓ(ℓ+1) \left\{ 𝔇^{(ℓ)}_{m', m} \right\}(𝐐),
\\
L_z \left\{ 𝔇^{(ℓ)}_{m', m} \right\}(𝐐)
&= m' \left\{ 𝔇^{(ℓ)}_{m', m} \right\}(𝐐),
\\
R_z \left\{ 𝔇^{(ℓ)}_{m', m} \right\}(𝐐)
&= m \left\{ 𝔇^{(ℓ)}_{m', m} \right\}(𝐐).
\end{aligned}
```
That is, Wigner's ``𝔇`` matrices are proportional to the
spin-weighted spherical harmonics.  We can get the proportionality
factor from ``𝔇^{(ℓ)}_{m', m}(1) = \delta_{m', m}``.

