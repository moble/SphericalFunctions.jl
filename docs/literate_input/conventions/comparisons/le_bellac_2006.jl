md"""
# Le Bellac (2006)

!!! info "Summary"
    Le Bellac's definition of the Wigner ``D`` matrix, his rotation law for the spherical
    harmonics, and his relation between the spherical harmonics and ``D`` all agree with the
    conventions used in the `SphericalFunctions` package.

[LeBellac_2006](@citet) (with Foreword by Cohen-Tannoudji) is a graduate textbook on quantum
mechanics.  Figure 10.1 shows that the spherical coordinates are the standard (physicist's)
coordinates.

## Wigner's ``D`` matrix

Le Bellac takes an odd approach, defining [Eq. (10.32)]
```math
D^{(j)}_{m', m} \left[ ℛ(θ, ϕ) \right]
=
\langle j, m' | e^{-iϕ J_z} e^{-iθ J_y} | j, m \rangle,
```
but later allowing that ``e^{-i ψ J_z}`` usually goes on the right-hand side of the others,
in which case ``D^{(j)}(θ, ϕ) \to D^{(j)}(ϕ, θ, ψ)``.  That is, the general rotation
operator is ``e^{-iϕ J_z} e^{-iθ J_y} e^{-iψ J_z}``, whose matrix elements are [our
``𝔇``](@ref summary_wigner_D) with Euler angles ``(α, β, γ) = (ϕ, θ, ψ)``.

Equation (10.65) shows the rotation law:
```math
Y_{ℓ}^{m}\left( ℛ^{-1} \hat{r} \right)
=
\sum_{m'} D^{(ℓ)}_{m', m}(ℛ) Y_{ℓ}^{m'}(\hat{r}),
```
which is precisely [our canonical form](@ref summary_spherical_harmonics), and Eq. (10.66)
relates the spherical harmonics to the Wigner D-matrices:
```math
D^{(ℓ)}_{m, 0}(θ, ϕ)
=
\sqrt{\frac{4π}{2ℓ+1}} \left[Y_{ℓ}^{m}(θ, ϕ)\right]^\ast,
```
again in agreement with [ours](@ref summary_spherical_harmonics).

Rather than transcribing an explicit formula for the ``d`` matrix from Le Bellac, we test
his definition directly: we construct the matrices of ``J_z`` and ``J_y = (J_+ - J_-)/2i`` in
the ``|j, m\rangle`` basis from the standard ladder relations, and exponentiate them
numerically.

## Implementing formulas

We begin by writing code that implements the formulas from Le Bellac.  We encapsulate the
formulas in a module so that we can test them against the `SphericalFunctions` package.
"""

# TODO: Confirm whether Le Bellac gives an explicit d-matrix formula; if so, transcribe and test it.  #src
using TestItems: @testitem  #hide
@testitem "Le Bellac conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide

module LeBellac
#+

# We'll use some predefined utilities to make the code look more like the equations, and
# the matrix exponential from `LinearAlgebra` for the rotation operator.
import ..ConventionsUtilities: 𝒾
import LinearAlgebra: exp, Diagonal
#+

# The matrices of ``J_z``, ``J_\pm``, and ``J_y`` in the ``|j, m\rangle`` basis, ordered so
# that row and column `k` correspond to ``m = -j + k - 1``:
Jz(j) = Diagonal([m for m ∈ -j:j])
function J₊(j)
    M = zeros(2j+1, 2j+1)
    for (k, m) ∈ enumerate(-j:j-1)
        M[k+1, k] = √((j-m) * (j+m+1))  # ⟨j, m+1| J₊ |j, m⟩
    end
    M
end
J₋(j) = J₊(j)'
Jy(j) = (J₊(j) - J₋(j)) / 2𝒾
#+

# The matrix element of Eq. (10.32), extended by ``e^{-iψJ_z}`` on the right:
function D(j, m′, m, ϕ, θ, ψ)
    U = exp(-𝒾 * ϕ * Matrix(Jz(j))) * exp(-𝒾 * θ * Matrix(Jy(j))) * exp(-𝒾 * ψ * Matrix(Jz(j)))
    U[m′ + j + 1, m + j + 1]
end
D(j, m′, m, θ, ϕ) = D(j, m′, m, ϕ, θ, zero(θ))
#+

end  # module LeBellac
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

# First, Le Bellac's ``D`` agrees with ours:
for (ϕ, θ, ψ) ∈ αβγs
    for (j, m′, m) ∈ ℓm′mrange(ℓₘₐₓ)
        @test LeBellac.D(j, m′, m, ϕ, θ, ψ) ≈ ConventionsUtilities.D(j, m′, m, ϕ, θ, ψ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Equation (10.66) holds with our spherical harmonics:
for (θ, ϕ) ∈ θϕrange(Float64, 7)
    for (ℓ, m) ∈ ℓmrange(ℓₘₐₓ)
        @test LeBellac.D(ℓ, m, 0, θ, ϕ) ≈ √(4π/(2ℓ+1)) * conj(ConventionsUtilities.Y(ℓ, m, θ, ϕ)) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Finally, the rotation law of Eq. (10.65).  We rotate the unit vector ``\hat{r}(θ, ϕ)`` by
# the *inverse* of the rotor ``𝐑_{α,β,γ}``, compute the spherical coordinates ``(θ', ϕ')`` of
# the result, and check that our spherical harmonics there are given by the stated
# combination of the harmonics at the original point.  We avoid the poles so that ``ϕ'`` is
# well defined.
import Quaternionic: from_euler_angles, from_spherical_coordinates, imz, components
for (α, β, γ) ∈ αβγrange(Float64, 3)
    R = from_euler_angles(α, β, γ)
    for (θ, ϕ) ∈ θϕrange(Float64, 3; avoid_poles=1e-3)
        r̂ = from_spherical_coordinates(θ, ϕ) * imz * conj(from_spherical_coordinates(θ, ϕ))
        r̂′ = conj(R) * r̂ * R  # ℛ⁻¹ r̂
        _, x, y, z = components(r̂′)
        θ′, ϕ′ = atan(hypot(x, y), z), atan(y, x)  # well conditioned near the poles
        for (ℓ, m) ∈ ℓmrange(ℓₘₐₓ)
            @test ConventionsUtilities.Y(ℓ, m, θ′, ϕ′) ≈ sum(
                LeBellac.D(ℓ, m′, m, α, β, γ) * ConventionsUtilities.Y(ℓ, m′, θ, ϕ)
                for m′ ∈ -ℓ:ℓ
            ) atol=ϵₐ rtol=ϵᵣ
        end
    end
end
#+

# These successful tests show that Le Bellac's ``D`` matrix, rotation law, and relation
# between the spherical harmonics and ``D`` agree with the conventions of the
# `SphericalFunctions` package.

end  #hide
