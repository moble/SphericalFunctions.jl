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
his definition directly: we construct the matrices of ``J_z`` and ``J_y = (J_+ - J_-)/2i``
in the ``|j, m\rangle`` basis from the standard ladder relations, and exponentiate them
numerically.

## Implementing formulas

We begin by writing code that implements the formulas from Le Bellac.  We encapsulate the
formulas in a module so that we can test them against the `SphericalFunctions` package.
"""

using TestItems: @testitem  #hide
@testitem "Le Bellac conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide
import .Utilities: αβγrange, θϕrange  #hide

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

# The matrix of the rotation operator ``e^{-iϕ J_z} e^{-iθ J_y} e^{-iψ J_z}`` in the same
# basis.  The first and last tests below use every element of each matrix, so each matrix is
# computed only once for each set of arguments, and is kept in `U_matrices` for the later
# elements; the three matrix exponentials cost far more than looking the result up.
const U_matrices = Dict{Any, Any}()
function U(j, ϕ, θ, ψ)
    get!(U_matrices, (j, ϕ, θ, ψ)) do
        exp(-𝒾 * ϕ * Matrix(Jz(j))) * exp(-𝒾 * θ * Matrix(Jy(j))) * exp(-𝒾 * ψ * Matrix(Jz(j)))
    end
end
#+

# The matrix element of Eq. (10.32), extended by ``e^{-iψJ_z}`` on the right:
D(j, m′, m, ϕ, θ, ψ) = U(j, ϕ, θ, ψ)[m′ + j + 1, m + j + 1]
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
αβγs = αβγrange(rng, Float64, 5)
#+

# First, Le Bellac's ``D`` agrees with ours:
for (ϕ, θ, ψ) ∈ αβγs
    for (j, 𝔇ʲ) ∈ SphericalFunctions.DCalculator(ϕ, θ, ψ, ℓₘₐₓ), m′ ∈ -j:j, m ∈ -j:j
        @test LeBellac.D(j, m′, m, ϕ, θ, ψ) ≈ 𝔇ʲ[m′, m] atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Equation (10.66) holds with our spherical harmonics:
for (θ, ϕ) ∈ θϕrange(rng, Float64, 7)
    for (ℓ, Yˡ) ∈ SphericalFunctions.YlmCalculator(θ, ϕ, ℓₘₐₓ), m ∈ -ℓ:ℓ
        @test LeBellac.D(ℓ, m, 0, θ, ϕ) ≈ √(4π/(2ℓ+1)) * conj(Yˡ[m]) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Finally, the rotation law of Eq. (10.65).  We rotate the unit vector ``\hat{r}(θ, ϕ)`` by
# the *inverse* of the rotor ``𝐑_{α,β,γ}``, compute the spherical coordinates ``(θ', ϕ')``
# of the result, and check that our spherical harmonics there are given by the stated
# combination of the harmonics at the original point.  The original points avoid the poles
# only to keep the grid generic; the two-argument `atan` below computes ``(θ', ϕ')`` stably
# wherever the rotated point lands.  Two calculators, one at each point, are iterated in
# step, so that the harmonics of each ``ℓ`` at both points are at hand together.
import Quaternionic: from_euler_angles, from_spherical_coordinates, imz, components
for (α, β, γ) ∈ αβγrange(rng, Float64, 3)
    R = from_euler_angles(α, β, γ)
    for (θ, ϕ) ∈ θϕrange(rng, Float64, 3; avoid_poles=1e-3)
        r̂ = from_spherical_coordinates(θ, ϕ) * imz * conj(from_spherical_coordinates(θ, ϕ))
        r̂′ = conj(R) * r̂ * R  # ℛ⁻¹ r̂
        _, x, y, z = components(r̂′)
        θ′, ϕ′ = atan(hypot(x, y), z), atan(y, x)  # well conditioned near the poles
        Y = SphericalFunctions.YlmCalculator(θ, ϕ, ℓₘₐₓ)
        Y′ = SphericalFunctions.YlmCalculator(θ′, ϕ′, ℓₘₐₓ)
        for ((ℓ, Yˡ), (_, Y′ˡ)) ∈ zip(Y, Y′), m ∈ -ℓ:ℓ
            @test Y′ˡ[m] ≈ sum(
                LeBellac.D(ℓ, m′, m, α, β, γ) * Yˡ[m′]
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
