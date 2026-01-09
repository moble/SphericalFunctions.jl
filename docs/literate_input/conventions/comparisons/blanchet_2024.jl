md"""
# Blanchet (2024)

!!! info "Summary"
    Blanchet's  definition of the spherical harmonics agrees with the definition used in the
    `SphericalFunctions` package.

Luc Blanchet is one of the pre-eminent researchers in post-Newtonian approximations, and has
written a "living" review article on the subject [Blanchet_2024](@cite), which he has kept
up-to-date with the latest developments.

The spherical coordinates are standard physicists' coordinates, implicitly defined by the
direction vector below Eq. (188b):
```math
  N_i = \left(\sin θ\cos ϕ, \sin θ\sin ϕ, \cos θ\right).
```


## Implementing formulas

We begin by writing code that implements the formulas from Ref. [Blanchet_2024](@cite).  We
encapsulate the formulas in a module so that we can test them against the
`SphericalFunctions` package.

"""
using TestItems: @testitem  #hide
@testitem "Blanchet conventions" setup=[ConventionsUtilities, ConventionsSetup, Utilities] begin  #hide

module Blanchet
#+

# We'll also use some predefined utilities to make the code look more like the equations.
import ..ConventionsUtilities: 𝒾, ❗
#+

# The ``s=-2`` spin-weighted spherical harmonics are defined in Eq. (184a) as
# ```math
#   Y^{l,m}_{-2} = \sqrt{\frac{2l+1}{4π}} d^{ℓ m}(θ) e^{imϕ}.
# ```
function Yˡᵐ₋₂(l, m, θ::T, ϕ::T) where {T<:Real}
    √((2l + 1) / (4T(π))) * d(l, m, θ) * exp(𝒾 * m * ϕ)
end
#+

# Immediately following that, in Eq. (184b), we find the definition of the ``d`` function:
# ```math
#   d^{ℓ m}
#   =
#   \sum_{k = k_1}^{k_2}
#   \frac{(-)^k}{k!}
#   e_k^{ℓ m}
#   \left(\cos\frac{θ}{2}\right)^{2ℓ+m-2k-2}
#   \left(\sin\frac{θ}{2}\right)^{2k-m+2},
# ```
# with ``k_1 = \textrm{max}(0, m-2)`` and ``k_2=\textrm{min}(l+m, l-2)``.
function d(l, m, θ::T) where {T<:Real}
    k₁ = max(0, m - 2)
    k₂ = min(l + m, l - 2)
    sum(
        T((-1)^k / (k)❗ * eₖˡᵐ(k, l, m))
         * cos(θ / 2) ^ (2l + m - 2k - 2)
         * sin(θ / 2) ^ (2k - m + 2)
        for k in k₁:k₂;
        init=zero(T)
    )
end
#+

# The ``e_k^{ℓ m}`` symbol is defined in Eq. (184c) as
# ```math
#   e_k^{ℓ m} = \frac{
#     \sqrt{(ℓ+m)!(ℓ-m)!(ℓ+2)!(ℓ-2)!}
#   }{
#     (k-m+2)!(ℓ+m-k)!(ℓ-k-2)!
#   }.
# ```
function eₖˡᵐ(k, l, m)
    (
        √((l + m)❗ * (l - m)❗ * (l + 2)❗ * (l - 2)❗)
        / ((k - m + 2)❗ * (l + m - k)❗ * (l - k - 2)❗)
    )
end
#+

# The paper did not give an expression for the Wigner D-matrices, but the definition of the
# spin-weighted spherical harmonics is probably most relevant, so this will suffice.

end  # module Blanchet
#+

# ## Tests
#
# We can now test the functions against the equivalent functions from the
# `SphericalFunctions` package.  We will test up to
ℓₘₐₓ = 8
#+

# because that's the maximum ``ℓ`` used for PN results — and that's roughly the limit to
# which I'd trust these expressions anyway.  We will also only test the
s = -2
#+

# case, which is the only one defined in the paper.
# We will need to test approximate floating-point equality,
# so we set absolute and relative tolerances (respectively) in terms of the machine epsilon:
ϵₐ = 30eps()
ϵᵣ = 1500eps()
#+

# This loose relative tolerance is necessary because the numerical errors in Blanchet's
# explicit expressions grow rapidly with ``ℓ``.
for (θ, ϕ) ∈ θϕrange()
    for (ℓ, m) ∈ ℓmrange(abs(s), ℓₘₐₓ)
        @test Blanchet.Yˡᵐ₋₂(ℓ, m, θ, ϕ) ≈ SphericalFunctions.Deprecated.Y(s, ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# These successful tests show that Blanchet's expression agrees with ours.


end  #hide
