md"""
# Wigner (1959)

!!! info "Summary"
    Wigner's representation matrices ``𝔇^{(j)}(α, β, γ)_{\mu'\mu}`` are related to ours
    by
    ```math
    𝔇^{(j)}(α, β, γ)_{\mu'\mu}\big|_{\text{Wigner}}
    =
    (-1)^{\mu'-\mu}\, \overline{𝔇^{(j)}_{\mu',\mu}(α, β, γ)}
    =
    𝔇^{(j)}_{-\mu',-\mu}(α, β, γ).
    ```
    That is, they are the complex conjugates of ours with an additional factor of
    ``(-1)^{\mu'-\mu}`` — or, equivalently, ours with the signs of both indices reversed.

[Wigner_1959](@citet) is the English translation of Wigner's 1931 book, which introduced
the representation matrices of the rotation group that now bear his name.  As with most
early sources, the conventions take some care to extract.

## Euler angles

Figure 2 on page 59 shows Wigner's Euler angles ``\{α, β, γ\}``.  In that figure, the
position of the ``z'`` axis is independent of ``α``, which appears to represent a final
rotation about the ``z'`` axis; in our convention, this rotation would be described by the
*final* Euler angle, ``γ``, so the labels appear to be swapped relative to ours.  On the
other hand, on page 156, if ``𝔇^{(ℓ)}`` obeys the representation-composition property, then
``\{α, β, γ\}`` represents the rotation ``\{α, 0, 0\} \circ \{0, β, 0\} \circ \{0, 0, γ\}``,
which is the same as our convention.  Wigner is most explicit about his Euler angles in
Appendix A: Eq. (A.2) gives the rotation matrix in terms of the Euler angles, and
multiplying it on the left by the column vector ``(0, 0, 1)`` shows where the point on the
``z`` axis is rotated.  It is independent of ``γ`` and depends on ``α`` (and ``β``), which
is inconsistent with Fig. 2 but consistent with our convention; we conclude that the labels
in Fig. 2 have simply been swapped, and that Wigner's Euler angles are otherwise ours.

## Representation matrices

Equation (15.27) gives the explicit expression
```math
𝔇^{(j)}(α, β, γ)_{\mu'\mu}
=
\sum_κ (-1)^κ
\frac{\sqrt{(j+\mu)!\,(j-\mu)!\,(j+\mu')!\,(j-\mu')!}}
     {(j-\mu'-κ)!\,(j+\mu-κ)!\,κ!\,(κ+\mu'-\mu)!}\,
e^{i\mu'α}
\left(\cos\frac{β}{2}\right)^{2j+\mu-\mu'-2κ}
\left(\sin\frac{β}{2}\right)^{2κ+\mu'-\mu}
e^{i\muγ},
```
which we implement below.  The phases ``e^{+i\mu'α}`` and ``e^{+i\muγ}`` are the conjugates
of [ours](@ref summary_wigner_D), and the ``β``-dependent factor differs from Wigner's own
``d`` formula as quoted by later authors by the sign ``(-1)^{\mu'-\mu}``.  Using [the
symmetry](@ref summary_wigner_D) ``\overline{𝔇_{\mu',\mu}} = (-1)^{\mu'-\mu}
𝔇_{-\mu',-\mu}`` of our matrices, the relation in the summary box follows, and is tested
below.

Wigner's relation to the spherical harmonics, Eq. (A.11), can be written as
```math
Y_{ℓ}^{m}(θ, ϕ)
=
c\, (-1)^m\, 𝔇^{(ℓ)}(ϕ, θ, 0)_{m 0},
```
for a positive normalization constant ``c`` — note the factor ``(-1)^m``, which does not
appear in [our relation](@ref summary_spherical_harmonics) ``Y_{ℓ,m} = \sqrt{(2ℓ+1)/4π}\,
\overline{𝔇_{m,0}}``.  In fact, with the relation between Wigner's ``𝔇`` and ours given
above, ``(-1)^m 𝔇^{(ℓ)}(ϕ, θ, 0)_{m0}|_{\text{Wigner}} = \overline{𝔇_{m,0}(ϕ, θ, 0)}``, so
the two relations are the same and Wigner's spherical harmonics are the standard ones.  We
verify this as well.

## Implementing formulas

We begin by writing code that implements the formulas from Wigner.  We encapsulate the
formulas in a module so that we can test them against the `SphericalFunctions` package.
"""

using TestItems: @testitem  #hide
@testitem "Wigner conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide

module Wigner
#+

# We'll use some predefined utilities to make the code look more like the equations.
import ..ConventionsUtilities: 𝒾, ❗
#+

# Equation (15.27).  The sum runs over all ``κ`` for which the factorials have non-negative
# arguments, ``\max(0, \mu-\mu') \leq κ \leq \min(j-\mu', j+\mu)``.
function 𝔇(j, μ′, μ, α::T, β::T, γ::T) where {T<:Real}
    sum(
        (-1)^κ
        * T(√((j+μ)❗ * (j-μ)❗ * (j+μ′)❗ * (j-μ′)❗)
            / ((j-μ′-κ)❗ * (j+μ-κ)❗ * (κ)❗ * (κ+μ′-μ)❗))
        * exp(𝒾 * μ′ * α) * cos(β/2)^(2j+μ-μ′-2κ) * sin(β/2)^(2κ+μ′-μ) * exp(𝒾 * μ * γ)
        for κ ∈ max(0, μ-μ′):min(j-μ′, j+μ);
        init=zero(Complex{T})
    )
end
#+

end  # module Wigner
#+

# ## Tests
#
# We can now test the functions against the equivalent functions from the
# `SphericalFunctions` package.  We will need to test approximate floating-point equality,
# so we set absolute and relative tolerances (respectively) in terms of the machine epsilon:
ϵₐ = 100eps()
ϵᵣ = 1000eps()
#+

# We only test up to
ℓₘₐₓ = 4
#+
# because the formulas are slow, and this will be sufficient to sort out any sign or
# normalization differences, which are the most likely source of error.  For the same reason
# we use a modest grid of Euler angles.
αβγs = αβγrange(Float64, 5)
#+

# Wigner's ``𝔇`` is ``(-1)^{\mu'-\mu}`` times the complex conjugate of ours, or equivalently
# ours with both indices negated:
for (α, β, γ) ∈ αβγs
    for (j, μ′, μ) ∈ ℓm′mrange(ℓₘₐₓ)
        @test Wigner.𝔇(j, μ′, μ, α, β, γ) ≈
            (-1)^(μ′-μ) * conj(ConventionsUtilities.D(j, μ′, μ, α, β, γ)) atol=ϵₐ rtol=ϵᵣ
        @test Wigner.𝔇(j, μ′, μ, α, β, γ) ≈
            ConventionsUtilities.D(j, -μ′, -μ, α, β, γ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# And Wigner's Eq. (A.11) reproduces the standard spherical harmonics (with ``c =
# \sqrt{(2ℓ+1)/4π}``):
for (θ, ϕ) ∈ θϕrange()
    for (ℓ, m) ∈ ℓmrange(ℓₘₐₓ)
        @test √((2ℓ+1) / (4π)) * (-1)^m * Wigner.𝔇(ℓ, m, 0, ϕ, θ, zero(θ)) ≈
            ConventionsUtilities.Y(ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# These successful tests show that Wigner's representation matrices are ``(-1)^{\mu'-\mu}``
# times the complex conjugates of the ``𝔇`` matrices defined by the `SphericalFunctions`
# package, and that his spherical harmonics are ours.

end  #hide
