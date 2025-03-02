md"""
# NINJA (2011)

!!! info "Summary"
    The NINJA collaboration's definitions of the spherical harmonics and Wigner's ``d``
    functions agrees with the definitions used in the `SphericalFunctions` package.

Motivated by the need for a shared set of conventions in the NINJA project, a broad
cross-section of researchers involved in modeling gravitational waves (including the author
of this package) prepared Ref. [AjithEtAl_2011](@cite).  It is worth noting that the first
version posted to the arXiv included an unfortunate typo in the definition of the
spin-weighted spherical harmonics.  This was corrected in the second version, and remains
correct in the final — third — version.

The spherical coordinates are standard physicists' coordinates, except that the polar angle
is denoted ``\iota``:

> we define standard spherical coordinates ``(r, ι, φ)`` where ``ι`` is the inclination
> angle from the z-axis and ``φ`` is the phase angle.


## Implementing formulas

We begin by writing code that implements the formulas from Ref. [AjithEtAl_2011](@cite).  We
encapsulate the formulas in a module so that we can test them against the SphericalFunctions
package.

"""
using TestItems: @testitem  #hide
@testitem "NINJA conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide

module NINJA

import ..ConventionsUtilities: 𝒾, ❗
#+

# The spin-weighted spherical harmonics are defined in Eq. (II.7) as
# ```math
#   {}^{-s}Y_{l,m} = (-1)^s \sqrt{\frac{2l+1}{4\pi}}
#   d^{l}_{m,s}(\iota) e^{im\phi}.
# ```
# Just for convenience, we eliminate the negative sign on the left-hand side:
# ```math
#   {}^{s}Y_{l,m} = (-1)^{-s} \sqrt{\frac{2l+1}{4\pi}}
#   d^{l}_{m,-s}(\iota) e^{im\phi}.
# ```
function ₛYₗₘ(s, l, m, ι::T, ϕ::T) where {T<:Real}
    (-1)^(-s) * √((2l + 1) / (4T(π))) * d(l, m, -s, ι) * exp(𝒾 * m * ϕ)
end
#+

# Immediately following that, in Eq. (II.8), we find the definition of Wigner's ``d``
# function (again, noting that this expression was incorrect in version 1 of the paper, but
# correct in versions 2 and 3):
# ```math
#   d^{l}_{m,s}(\iota)
#   =
#   \sum_{k = k_1}^{k_2}
#   \frac{(-1)^k [(l+m)!(l-m)!(l+s)!(l-s)!]^{1/2}}
#   {(l+m-k)!(l-s-k)!k!(k+s-m)!}
#   \left(\cos\left(\frac{\iota}{2}\right)\right)^{2l+m-s-2k}
#   \left(\sin\left(\frac{\iota}{2}\right)\right)^{2k+s-m},
# ```
# with ``k_1 = \textrm{max}(0, m-s)`` and ``k_2=\textrm{min}(l+m, l-s)``.
function d(l, m, s, ι::T) where {T<:Real}
    k₁ = max(0, m - s)
    k₂ = min(l + m, l - s)
    sum(
        (-1)^k
        * T(
            ((l + m)❗ * (l - m)❗ * (l + s)❗ * (l - s)❗)^(1//2)
           / ((l + m - k)❗ * (l - s - k)❗ * (k)❗ * (k + s - m)❗)
        )
         * cos(ι / 2) ^ (2l + m - s - 2k)
         * sin(ι / 2) ^ (2k + s - m)
        for k in k₁:k₂;
        init=zero(T)
    )
end
#+

# For reference, several explicit formulas are also provided in Eqs.  (II.9)--(II.13):
# ```math
# \begin{aligned}
#   {}^{-2}Y_{2,2} &= \sqrt{\frac{5}{64\pi}} (1+\cos\iota)^2 e^{2i\phi},\\
#   {}^{-2}Y_{2,1} &= \sqrt{\frac{5}{16\pi}} \sin\iota (1 + \cos\iota) e^{i\phi},\\
#   {}^{-2}Y_{2,0} &= \sqrt{\frac{15}{32\pi}} \sin^2\iota,\\
#   {}^{-2}Y_{2,-1} &= \sqrt{\frac{5}{16\pi}} \sin\iota (1 - \cos\iota) e^{-i\phi},\\
#   {}^{-2}Y_{2,-2} &= \sqrt{\frac{5}{64\pi}} (1-\cos\iota)^2 e^{-2i\phi}.
# \end{aligned}
# ```
₋₂Y₂₂(ι::T, ϕ::T) where {T<:Real} = √(5 / (64T(π))) * (1 + cos(ι))^2 * exp(2𝒾*ϕ)
₋₂Y₂₁(ι::T, ϕ::T) where {T<:Real} = √(5 / (16T(π))) * sin(ι) * (1 + cos(ι)) * exp(𝒾*ϕ)
₋₂Y₂₀(ι::T, ϕ::T) where {T<:Real} = √(15 / (32T(π))) * sin(ι)^2
₋₂Y₂₋₁(ι::T, ϕ::T) where {T<:Real} = √(5 / (16T(π))) * sin(ι) * (1 - cos(ι)) * exp(-𝒾*ϕ)
₋₂Y₂₋₂(ι::T, ϕ::T) where {T<:Real} = √(5 / (64T(π))) * (1 - cos(ι))^2 * exp(-2𝒾*ϕ)
#+

# The paper did not give an expression for the Wigner D-matrices, but the definition of the
# spin-weighted spherical harmonics is probably most relevant, so this will suffice.

end  # module NINJA
#+

# ## Tests
#
# We can now test the functions against the equivalent functions from the SphericalFunctions
# package.  We will need to test approximate floating-point equality, so we set absolute and
# relative tolerances (respectively) in terms of the machine epsilon:
ϵₐ = 10eps()
ϵᵣ = 10eps()
#+

# First, we compare the explicit formulas to the general formulas.
for (ι, ϕ) ∈ θϕrange(Float64, 1)
    @test NINJA.ₛYₗₘ(-2, 2, 2, ι, ϕ) ≈ NINJA.₋₂Y₂₂(ι, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test NINJA.ₛYₗₘ(-2, 2, 1, ι, ϕ) ≈ NINJA.₋₂Y₂₁(ι, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test NINJA.ₛYₗₘ(-2, 2, 0, ι, ϕ) ≈ NINJA.₋₂Y₂₀(ι, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test NINJA.ₛYₗₘ(-2, 2, -1, ι, ϕ) ≈ NINJA.₋₂Y₂₋₁(ι, ϕ) atol=ϵₐ rtol=ϵᵣ
    @test NINJA.ₛYₗₘ(-2, 2, -2, ι, ϕ) ≈ NINJA.₋₂Y₂₋₂(ι, ϕ) atol=ϵₐ rtol=ϵᵣ
end
#+

# Next, we compare the general formulas to the SphericalFunctions package.
# We will only test up to
ℓₘₐₓ = 4
#+
# and
sₘₐₓ = 2
#+
# because the formulas are very slow, and this will be sufficient to sort out any sign or
# normalization differences, which are the most likely source of error.
for (θ, ϕ) ∈ θϕrange()
    for (s, ℓ, m) ∈ sℓmrange(ℓₘₐₓ, sₘₐₓ)
        @test NINJA.ₛYₗₘ(s, ℓ, m, θ, ϕ) ≈ SphericalFunctions.Y(s, ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# Finally, we compare the Wigner ``d`` matrix to the SphericalFunctions package.
for ι ∈ θrange()
    for (ℓ, m′, m) ∈ ℓm′mrange(ℓₘₐₓ)
        @test NINJA.d(ℓ, m′, m, ι) ≈ SphericalFunctions.d(ℓ, m′, m, ι) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# These successful tests show that both the spin-weighted spherical harmonics and the Wigner
# ``d`` matrix defined by the NINJA collaboration agree with the corresponding functions
# defined by the `SphericalFunctions` package.

end  #hide
