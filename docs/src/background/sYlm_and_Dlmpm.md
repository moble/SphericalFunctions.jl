# [``{}_{s}Y_{ℓ,m}`` and ``𝔇^{(ℓ)}_{m', m}``](@id sYlm_and_Dlmpm)

The [previous page](@ref background_differential_operators) introduced
the left and right Lie derivative operators acting on functions of a
quaternion argument.  Here, we show how the spin-weighted spherical
harmonics arise naturally as eigenfunctions of these operators.

## ``{}_{s}Y_{ℓ,m}`` as eigenfunctions

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
normally be included in quantum mechanics texts.)  Familiar arguments
— either from spectra of operators, or continuity of the functions —
tell us that the allowed values of ``ℓ`` are non-negative integers,
and for each ``ℓ``, the allowed values of ``m`` are integers
satisfying ``-ℓ ≤ m ≤ ℓ``.

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
L² \left\{ {}_{s}Y_{ℓ,m} \right\}(𝐐)
&= ℓ(ℓ+1) \left\{ {}_{s}Y_{ℓ,m} \right\}(𝐐),
\\
L_z \left\{ {}_{s}Y_{ℓ,m} \right\}(𝐐)
&= m \left\{ {}_{s}Y_{ℓ,m} \right\}(𝐐),
\\
R_z \left\{ {}_{s}Y_{ℓ,m} \right\}(𝐐)
&= s \left\{ {}_{s}Y_{ℓ,m} \right\}(𝐐).
\end{aligned}
```
In this case, similar arguments show that the allowed values of ``ℓ``
are non-negative integers or half-integers, and for each ``ℓ``, the
allowed values of ``m`` and ``s`` are correspondingly integers or
half-integers satisfying ``-ℓ ≤ m ≤ ℓ`` and ``-ℓ ≤ s ≤ ℓ``.  These
conditions only determine the functions up to a normalization factor,
which we discuss in the next section.

First, though, we take a moment to consider the effect of the raising
and lowering operators.  We know the commutators from the previous
page, and we now know the eigenvalues, so we can compute just like we
do in elementary treatments of scalar spherical harmonics.  ``L_±``
operates just the same as usual, modifying the ``m`` index, while
``R_±`` in an exactly analogous way, but modifies the ``s`` index.
Thus, we can follow the standard derivation to find the standard
ladder relations:
```math
\begin{aligned}
L_± \left\{ {}_{s}Y_{ℓ,m} \right\}(𝐐)
&= \sqrt{(ℓ ∓ m)(ℓ ± m + 1)}\, \left\{ {}_{s}Y_{ℓ,m±1} \right\}(𝐐),
\\
R_± \left\{ {}_{s}Y_{ℓ,m} \right\}(𝐐)
&= \sqrt{(ℓ ∓ s)(ℓ ± s + 1)}\, \left\{ {}_{s±1}Y_{ℓ,m} \right\}(𝐐).
\end{aligned}
```
As usual, we have chosen the coefficients to be real and positive.

## Integration and normalization

It is natural to define an inner product on the space of functions on
``\mathrm{Spin}(3)`` using integration over the group itself.  That is,
for two functions ``f`` and ``g``, we define[^1]
```math
⟨f|g⟩_{\mathrm{Spin}(3)} = \int_{\mathrm{Spin}(3)} f̄(𝐐)\, g(𝐐)\, d𝐐,
```
where the integration can be implemented using [extended Euler
angles](@ref Quaternions-and-Euler-angles).  This induces a norm on
the space of functions on ``\mathrm{Spin}(3)``:
```math
\left\| f \right\|²_{\mathrm{Spin}(3)}
=⟨f|f⟩_{\mathrm{Spin}(3)}
=\int_{\mathrm{Spin}(3)} |f(𝐐)|²\, d𝐐.
```
One important point is the norm of the constant function ``f=1``,
which gives the total volume of ``\mathrm{Spin}(3)``.  We can compute
it explicitly using extended Euler angles [and find](@ref
Quaternions-and-Euler-angles) that the result is ``2π²``.  Recalling
that ``\mathrm{Spin}(3)`` considered as a subset of ``ℝ⁴`` is just the
unit three-sphere ``𝕊³``, which has volume ``2π²``, this makes sense.

[^1]: Note that we are using the physicists' bra-ket notation here
    [Dirac_1939](@cite).  Specifically, the inner product is
    conjugate-linear in its first argument and linear in its second
    argument: for complex numbers ``a`` and ``b``, we have
    ``⟨af|bg⟩=ā⟨f|g⟩b``.  This convention is common in physics,
    engineering, and computer science.  But note that mathematicians
    more commonly write ``⟨f, g⟩`` (with a comma instead of a vertical
    bar) and define the inner product with opposite linearity: linear
    in its first argument and conjugate-linear in its second argument,
    so that ``⟨af|bg⟩=a⟨f|g⟩b̄``.

Unfortunately, these are not actually the conventional inner product
and norm used in the literature.  The difference is not too difficult
to understand or deal with, however.  Note that ``|f(𝐐)|²`` is
actually a field of spin weight 0, no matter what the spin weight of
``f`` itself is.[^2]  Therefore, we can push it forward to a
nontrivial function on ``𝕊²`` as described [previously](@ref
Pushing-forward-to-𝕊), and then integrate over ``𝕊²`` instead of
``\mathrm{Spin}(3)``.  Thus, we can define a distinct second norm
```math
\left\| f \right\|²_{𝕊²} = \int_{𝕊²} |f|²\, dΩ.
```
This is actually the conventional norm used in the literature on
spin-weighted spherical functions.  Again, we can compute the norm of
the constant function ``f=1`` and find that it is just the area of the
unit 2-sphere, ``4π``.  We can relate the two norms with a simple
constant factor:
```math
\left\| f \right\|²_{𝕊²}
=\frac{2}{π} \left\| f \right\|²_{\mathrm{Spin}(3)}.
```
Clearly, we must choose one or the other of these norms when
normalizing the spin-weighted spherical harmonics.  To agree with
standard scalar spherical harmonics, the conventional choice is to use
the ``𝕊²`` norm.

[^2]: It is possible for a function to have no definite spin weight;
    to simply not be an eigenfunction of ``R_z``.  In that case,
    ``|f(𝐐)|²`` can also have indefinite spin weight.  Nonetheless,
    its integral over ``\mathrm{Spin}(3)`` will only pick up
    contributions from the part of that function with spin weight 0.

!!! danger "#TODO"

    Finish this section, noting that SWSHs with different 
    spin weights are not orthogonal under
    the ``𝕊²`` inner product.


## Defining ``𝔇^{(ℓ)}_{m', m}``

Finally, we can consider the Wigner ``𝔇`` matrices.  These are
usually defined as
```math
𝔇^{(ℓ)}_{m', m}(α, β, γ)
= \big\langle ℓ, m' \big|
  e^{-iL_z α} e^{-iL_y β} e^{-iL_z γ}
  \big| ℓ, m \big\rangle.
```
Note that the bra-ket notation usually represents integration over
``𝕊²``.  We can extend that to obtain the corresponding norm, while
integrating over ``\mathrm{Spin}(3)`` instead.  Thus, we can rewrite
this as
```math
𝔇^{(ℓ)}_{m', m}(𝐐)
= \frac{2}{π} \int_{\mathrm{Spin}(3)}
  \bar{Y}_{ℓ,m'}(𝐏)\,
  Y_{ℓ,m}\left(𝐐⁻¹ 𝐏\right)\,
  d𝐏.
```
We can evaluate the action of the differential operators on this
function pretty easily, noting that the derivative ``d/dϵ`` can pass
through the integral sign by the Leibniz integral rule.  The result is
that ``𝔇`` satisfies
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
factor by applying the definition above with ``𝐐=𝟏``, in which case
the integral simplifies to the orthonormality condition for the
spin-weighted spherical harmonics:
```math
𝔇^{(ℓ)}_{m', m}(𝟏) = \delta_{m', m}.
```

!!! danger "#TODO"

    Check the signs on the eigenvalues above, and refer to 
    the Notes page on normalization to relate Y to D.
