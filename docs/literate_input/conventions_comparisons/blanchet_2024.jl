md"""
# Blanchet (2024)

!!! info "Summary"
    TODO

Luc Blanchet is one of the pre-eminent researchers in post-Newtonian approximations, and has
written a "living" review article on the subject [Blanchet_2024](@cite), which he has kept
up-to-date with the latest developments.

The spherical coordinates are standard physicists' coordinates, except that the polar angle
is denoted ``\iota``:

> we define standard spherical coordinates ``(r, ι, φ)`` where ``ι`` is the inclination
> angle from the z-axis and ``φ`` is the phase angle.


## Implementing formulas

We begin by writing code that implements the formulas from Ref. [AjithEtAl_2011](@cite).  We
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
#   Y^{l,m}_{-2} = \sqrt{\frac{2l+1}{4\pi}} d^{\ell m}(\theta) e^{im\phi}.
# ```
function Yˡᵐ₋₂(l, m, θ::T, ϕ::T) where {T<:Real}
    √((2l + 1) / (4T(π))) * d(l, m, θ) * exp(𝒾 * m * ϕ)
end
#+

# Immediately following that, in Eq. (184b), we find the definition of the ``d`` function:
# ```math
#   d^{\ell m}(\theta)
#   =
#   \sum_{k = k_1}^{k_2}
#   \frac{(-1)^k}{k!}
#   e_k^{\ell m}
#   \left(\cos\frac{\theta}{2}\right)^{2\ell+m-2k-2}
#   \left(\sin\frac{\theta}{2}\right)^{2k-m+2},
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

# Note that he seems to have flipped the sign of ``s=-2`` in that equation, so that he has
# evidently given formulas for the ``s=2`` harmonics instead.  His notation seems consistent
# with the NINJA paper [AjithEtAl_2011](@cite), which unfortunately included a confusing
# negative sign on the left-hand side of the definition of the spin-weighted spherical
# harmonics.
#
# The ``e_k^{\ell m}`` symbol is defined in Eq. (184c) as
# ```math
#   e_k^{\ell m} = \frac{
#     \sqrt{(\ell+m)!(\ell-m)!(\ell+2)!(\ell-2)!}
#   }{
#     (k-m+2)!(\ell+m-k)!(\ell-k-2)!
#   }.
# ```
function eₖˡᵐ(k, l, m)
    (
        √((l + m)❗ * (l - m)❗ * (l + 2)❗ * (l - 2)❗)
        / ((k - m + 2)❗ * (l + m - k)❗ * (l - k - 2)❗)
    )
end


end  # module Blanchet
#+

# ## Tests
#
# We can now test the functions against the equivalent functions from the
# `SphericalFunctions` package.  We will need to test approximate floating-point equality,
# so we set absolute and relative tolerances (respectively) in terms of the machine epsilon:
ϵₐ = 10eps()
ϵᵣ = 10eps()
#+

# We will only test up to
ℓₘₐₓ = 2
#+
# because the formulas are very slow, and this will be sufficient to sort out any sign or
# normalization differences, which are the most likely source of error.
for (θ, ϕ) ∈ θϕrange(Float64, 1)
    for (ℓ, m) ∈ ℓmrange(ℓₘₐₓ)
        @test Blanchet.Yˡᵐ₋₂(ℓ, m, θ, ϕ) ≈ SphericalFunctions.Y(2, ℓ, m, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
    end
end
#+

# These successful tests show that TODO: finish this sentence


end  #hide
