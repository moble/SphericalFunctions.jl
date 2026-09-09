md"""
# Mathematica (2026)

!!! info "Summary"
    Mathematica's `SphericalHarmonicY` agrees with the spherical harmonics used in the
    `SphericalFunctions` package, and its `EulerMatrix` agrees with our Euler angles.  Its
    `WignerD` function is documented only through a set of identities; those identities are
    all satisfied by
    ```math
    \mathtt{WignerD[\{j, m_1, m_2\}, ψ, θ, ϕ]}
    = 𝔇^{(j)}_{-m_1,-m_2}(ψ, θ, ϕ)
    = (-1)^{m_1 - m_2}\, \overline{𝔇^{(j)}_{m_1,m_2}(ψ, θ, ϕ)},
    ```
    which is [Wigner's](@ref "Wigner (1959)") convention, and are *not* all satisfied by
    the plain complex conjugate ``\overline{𝔇^{(j)}_{m_1,m_2}}``.

We cannot run Mathematica as part of this package's test suite, so this page relies on the
statements in the [Wolfram Language documentation](@cite Mathematica_WignerD), accessed
2026-09-08.  (`WignerD` was introduced in Version 8.0, 2010; `SphericalHarmonicY` in
Version 1.0, 1988.)

## Euler angles

The Euler angles are defined generally such that

> `EulerMatrix[{α,β,γ},{a,b,c}]` is equivalent to ``R_{α,a} R_{β,b} R_{γ,c}``, where
> ``R_{α,a}``=`RotationMatrix[α,UnitVector[3,a]]`, etc.

and

> `EulerMatrix[{α,β,γ}]` is equivalent to `EulerMatrix[{α,β,γ},{3,2,3}]`

(representing the ``z-y-z`` convention).  Finally, we find that they say that `EulerMatrix`
corresponds to three rotations:

```mathematica
rα = RotationMatrix[α, {0, 0, 1}];
rβ = RotationMatrix[β, {0, 1, 0}];
rγ = RotationMatrix[γ, {0, 0, 1}];

Simplify[rα . rβ . rγ == EulerMatrix[{α, β, γ}]]
```

This agrees with [the conventions used in this package](@ref summary_euler_angles), so we
can directly compare expressions in terms of Euler angles.

## Spherical harmonics

The [`SphericalHarmonicY`](@cite Mathematica_SphericalHarmonicY) documentation states

> For ``ℓ \geq 0``, ``Y_ℓ^m(θ, ϕ) = \sqrt{(2ℓ+1)/(4π)} \sqrt{(ℓ-m)! / (ℓ+m)!}
> P_ℓ^m(\cos θ) e^{imϕ}`` where ``P_ℓ^m`` is the associated Legendre function,

and the [`LegendreP`](@cite Mathematica_LegendreP) documentation defines that function
(for its default "type 1") as

> The associated Legendre polynomials are defined by ``P_n^m(x) = (-1)^m (1-x^2)^{m/2}
> (d^m/dx^m) P_n(x)`` where ``P_n(x)`` is the Legendre polynomial.

The Condon–Shortley phase is thus included in `LegendreP`, and we expect agreement with
[our spherical harmonics](@ref summary_spherical_harmonics).  For negative ``m``, the
documentation does not spell out the continuation of ``P_n^m``; we use the standard
relation for Legendre functions with this phase convention (which is what Mathematica
evaluates to), [DLMF 14.9.3](https://dlmf.nist.gov/14.9#E3),
```math
P_n^{-m}(x) = (-1)^m \frac{(n-m)!}{(n+m)!} P_n^m(x).
```
The Legendre polynomial itself is given by Rodrigues' formula, [DLMF
14.7.13](https://dlmf.nist.gov/14.7#E13), ``P_n(x) = \frac{1}{2^n n!} \frac{d^n}{dx^n} (x^2 -
1)^n``.

## Wigner ``D`` function

The [`WignerD`](@cite Mathematica_WignerD) documentation does not give a general formula.
Instead, it states the following identities:

> The Wolfram Language uses phase conventions where ``D^j_{m_1, m_2}(ψ, θ, ϕ) = \exp(i m_1
> ψ + i m_2 ϕ) D^j_{m_1, m_2}(0, θ, 0)``.

> `WignerD[{j, m1, m2}, ψ, θ, ϕ] == (-1)^(m1 - m2) Conjugate[WignerD[{j, -m1, -m2}, ψ, θ,
> ϕ]]`

> `WignerD[{j, m1, m2}, ψ, θ, ϕ] == (-1)^(m1 - m2) WignerD[{j, m2, m1}, ϕ, θ, ψ]`

> `WignerD[{𝓁, 0, m}, θ, ϕ] == Sqrt[(4 π)/(2 𝓁 + 1)] SphericalHarmonicY[𝓁, m, θ, ϕ]`

> `WignerD[{1, 0, 1}, ψ, θ, ϕ]` = ``-\sqrt{2} e^{i ϕ} \cos\frac{θ}{2} \sin\frac{θ}{2}``

(The two-argument form `WignerD[{j, m1, m2}, θ, ϕ]` is `WignerD[{j, m1, m2}, 0, θ, ϕ]`.)
The phases ``e^{+im_1ψ}`` and ``e^{+im_2ϕ}`` are the complex conjugates of
[ours](@ref summary_wigner_D), so Mathematica's ``D`` is certainly not equal to our ``𝔇``.
There are two natural candidates for the relation: the plain complex conjugate,
``\overline{𝔇^{(j)}_{m_1,m_2}}``, and ``(-1)^{m_1-m_2} \overline{𝔇^{(j)}_{m_1,m_2}} =
𝔇^{(j)}_{-m_1,-m_2}`` (which is [Wigner's own convention](@ref "Wigner (1959)")).  Both
satisfy the first three identities.  The last two identities distinguish them: for the
spherical-harmonic relation, our conventions give ``\sqrt{4π/(2ℓ+1)}\, Y_{ℓ,m}(θ, ϕ) =
\overline{𝔇_{m,0}(ϕ, θ, 0)} = e^{imϕ} d_{m,0}(θ)``, whereas
``\overline{𝔇_{0,m}(0, θ, ϕ)} = e^{imϕ} d_{0,m}(θ) = (-1)^m e^{imϕ} d_{m,0}(θ)``; and the
explicit example gives ``-\sin θ / \sqrt{2}`` where ``d_{0,1}(θ) = +\sin θ/\sqrt{2}`` but
``d_{1,0}(θ) = -\sin θ/\sqrt{2}``.  So the plain conjugate fails both, and the second
candidate satisfies both.  We therefore take
```math
\mathtt{WignerD[\{j, m_1, m_2\}, ψ, θ, ϕ]} = 𝔇^{(j)}_{-m_1,-m_2}(ψ, θ, ϕ)
```
as the relation, and verify below that it satisfies every documented identity.  Note that
this is as far as the documentation allows us to go: the identities determine the ``m_1 = 0``
row and ``m_2 = 0`` column of Mathematica's ``D`` completely, but not the general element.

## Implementing formulas

We encapsulate the formulas in a module so that we can test them against the
`SphericalFunctions` package.
"""

using TestItems: @testitem  #hide
@testitem "Mathematica conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide

module Mathematica
#+

# We'll use some predefined utilities to make the code look more like the equations,
# including `∂ⁿ`, which computes derivatives symbolically so that we can transcribe the
# Legendre formulas literally, and our reference ``𝔇`` for the candidate `WignerD`.
import ..ConventionsUtilities: 𝒾, ❗, ∂ⁿ, D
#+

# The Legendre polynomial (Rodrigues' formula) and Mathematica's `LegendreP[n, m, x]`, with
# DLMF 14.9.3 for negative ``m``:
function LegendreP(n, x)
    1 / (2^n * (n)❗) * ∂ⁿ(x -> (x^2 - 1)^n, n)(x)
end
function LegendreP(n, m, x)
    if m < 0
        return (-1)^m * (n+m)❗ / (n-m)❗ * LegendreP(n, -m, x)
    end
    (-1)^m * (1 - x^2)^(m/2) * ∂ⁿ(x -> LegendreP(n, x), m)(x)
end
#+

# `SphericalHarmonicY[l, m, θ, ϕ]`.  We capture the floating-point type `T` to ensure that
# we don't lose precision when converting π and the factorials to floating-point numbers.
function SphericalHarmonicY(ℓ, m, θ::T, ϕ::T) where {T<:Real}
    √T((2ℓ+1) / (4big(π))) * √T((ℓ-m)❗ / (ℓ+m)❗) * T(LegendreP(ℓ, m, cos(θ))) * exp(𝒾 * m * ϕ)
end
#+

# The candidate for `WignerD[{j, m1, m2}, ψ, θ, ϕ]`, and its two-argument form:
WignerD(j, m₁, m₂, ψ, θ, ϕ) = D(j, -m₁, -m₂, ψ, θ, ϕ)
WignerD(j, m₁, m₂, θ, ϕ) = WignerD(j, m₁, m₂, zero(θ), θ, ϕ)
#+

end  # module Mathematica
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

# We will also need the imaginary unit from the utilities module.
import .ConventionsUtilities: 𝒾
#+

# First, `SphericalHarmonicY` agrees with ours:
for (θ, ϕ) ∈ θϕrange(Float64, 7)
    for (ℓ, m) ∈ ℓmrange(ℓₘₐₓ)
        @test Mathematica.SphericalHarmonicY(ℓ, m, θ, ϕ) ≈ ConventionsUtilities.Y(ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Now the documented identities for `WignerD`, in the order quoted above.  The phase
# convention:
for (ψ, θ, ϕ) ∈ αβγs
    for (j, m₁, m₂) ∈ ℓm′mrange(ℓₘₐₓ)
        @test Mathematica.WignerD(j, m₁, m₂, ψ, θ, ϕ) ≈
            exp(𝒾 * m₁ * ψ + 𝒾 * m₂ * ϕ) * Mathematica.WignerD(j, m₁, m₂, zero(θ), θ, zero(θ)) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# The two symmetries:
for (ψ, θ, ϕ) ∈ αβγs
    for (j, m₁, m₂) ∈ ℓm′mrange(ℓₘₐₓ)
        @test Mathematica.WignerD(j, m₁, m₂, ψ, θ, ϕ) ≈
            (-1)^(m₁-m₂) * conj(Mathematica.WignerD(j, -m₁, -m₂, ψ, θ, ϕ)) atol=ϵₐ rtol=ϵᵣ
        @test Mathematica.WignerD(j, m₁, m₂, ψ, θ, ϕ) ≈
            (-1)^(m₁-m₂) * Mathematica.WignerD(j, m₂, m₁, ϕ, θ, ψ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# The relation to `SphericalHarmonicY`, and the explicit example — which, as explained
# above, are the identities that rule out the plain complex conjugate of our ``𝔇``:
for (θ, ϕ) ∈ θϕrange(Float64, 7)
    for (ℓ, m) ∈ ℓmrange(ℓₘₐₓ)
        @test Mathematica.WignerD(ℓ, 0, m, θ, ϕ) ≈
            √(4π / (2ℓ+1)) * Mathematica.SphericalHarmonicY(ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
    for ψ ∈ (0.0, 0.7, 2.9)
        @test Mathematica.WignerD(1, 0, 1, ψ, θ, ϕ) ≈
            -√2 * exp(𝒾 * ϕ) * cos(θ/2) * sin(θ/2) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# These successful tests show that Mathematica's `SphericalHarmonicY` agrees with our
# spherical harmonics, and that every identity its documentation states for `WignerD` is
# satisfied by ``𝔇^{(j)}_{-m_1,-m_2}(ψ, θ, ϕ) = (-1)^{m_1-m_2}\,\overline{𝔇^{(j)}_{m_1,m_2}(ψ,
# θ, ϕ)}``.

end  #hide
