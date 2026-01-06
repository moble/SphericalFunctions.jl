raw"""
Formulas and conventions from [Wigner's "Group Theory and Its Applications to the Quantum
Mechanics of Atomic Spectra"](@cite Wigner_1959).

The conclusion here is that Wigner's ``Dʲₘₚ,ₘ`` includes a factor of ``(-1)^{m'-m}``
relative to ours.
"""
@testmodule Wigner begin

const 𝒾 = im

include("../utilities/naive_factorial.jl")
import .NaiveFactorials: ❗


raw"""
Figure 2 on page 59 shows that Wigner's Euler angles notated {α, β, γ} are simply swapped
with respect to ours.  For example, note that the position of ``z'`` is independent of α,
which appears to represent a final rotation about the ``z'`` axis.  In our convention, this
rotation would be described by the final Euler angle, γ.

On the other hand, on page 156, if ``𝔇^{(\ell)}`` obeys the representation-composition
property, then {α, β, γ} represents the rotation {α, 0, 0}∘{0, β, 0}∘{0, 0, γ}, which is the
same as our convention.

Wigner is most explicit about his Euler angles in Appendix A.  Specifically, note Eq. (A.2),
where he gives the rotation matrix in terms of the Euler angles.  If we multiply on the left
by the column vector (0, 0, 1), we can see where the point on the ``z`` axis should be
rotated in Fig. 2.  Specifically, it is independent of γ, and depends on α (and β).  This is
inconsistent with Fig. 2, indicating that the labels for α and γ have just been swapped.

Eq. (15.8) says

    D(ℓ, m′, m, α, β, γ) = exp(𝒾*m′*α) * d(ℓ, m′, m, β) * exp(𝒾*m*γ)

while Eq. (A.11) says

    Y(ℓ, m, θ, ϕ) = c * (-1)^m * exp(𝒾*m*ϕ) * d(ℓ, m, 0, θ)
                  = c * (-1)^m * D(ℓ, m, 0, ϕ, θ, 0)

Note the factor of (-1)^m, which I do not include.  But also note that there is no
conjugation of the D function, which is consistent with our convention.

"""


# @doc raw"""
#     D(ℓ, m′, m, α, β, γ)

# Eq. (15.8) of [Wigner](@cite Wigner_1959), implementing
# ```math
#     D^\ell_{m',m}(\alpha, \beta, \gamma).
# ```
# """
# function D(ℓ, m′, m, α, β, γ)
#     exp(𝒾*m′*α) * d(ℓ, m′, m, β) * exp(𝒾*m*γ)
# end


@doc raw"""
    D(ℓ, m′, m, α, β, γ)

Eq. (15.27) of [Wigner](@cite Wigner_1959), implementing
```math
    D^{(j)}(\alpha, \beta, \gamma)_{\mu',\mu}.
```
"""
function D(j, μ′, μ, α, β, γ)
    if j < 0
        throw(DomainError("The degree j=$j must be non-negative."))
    end
    if abs(μ′) > j
        throw(DomainError("The order abs(μ′)=$μ′ must be ≤ j=$j."))
    end
    if abs(μ) > j
        throw(DomainError("The order abs(μ)=$μ must be ≤ j=$j."))
    end

    α, β, γ = promote(α, β, γ)

    # The summation index `κ` ranges over all values for which the factorial arguments are
    # valid.
    κₘᵢₙ = max(0, μ-μ′)
    κₘₐₓ = min(j-μ′, j+μ)

    sinβ╱2, cosβ╱2 = sincos(β/2)
    T = typeof(sinβ╱2)

    sum(
        κ -> (
            (-1)^(κ)
            * T(√((j+μ)❗ * (j-μ)❗ * (j+μ′)❗ * (j-μ′)❗)
                / ((j-μ′-κ)❗ * (j+μ-κ)❗ * (κ)❗ * (κ+μ′-μ)❗))
            * exp(𝒾*μ′*α) * cosβ╱2^(2j+μ-μ′-2κ) * sinβ╱2^(2κ+μ′-μ) * exp(𝒾*μ*γ)
        ),
        κₘᵢₙ:κₘₐₓ,
        init=complex(zero(T))
    )
end

end # @testmodule Wigner


@testitem "Wigner conventions" setup=[Utilities, Wigner] begin
    using Random
    using Quaternionic: from_spherical_coordinates
    using SphericalFunctions: Deprecated

    Random.seed!(1234)
    const T = Float64
    const ℓₘₐₓ = 5
    ϵₐ = 2eps(T)
    ϵᵣ = 50eps(T)

    # Tests for D(j, m′, m, α, β, γ)
    let ϵₐ=√ϵᵣ, ϵᵣ=√ϵᵣ, 𝒟=Wigner.D
        for α ∈ αrange(T)
            for β ∈ βrange(T)
                for γ ∈ γrange(T)
                    D = Deprecated.D_matrices(α, β, γ, ℓₘₐₓ)
                    i = 1
                    for j in 0:ℓₘₐₓ
                        for m′ in -j:j
                            for m in -j:j
                                @test (-1)^(m′-m) * 𝒟(j, m′, m, α, β, γ) ≈ D[i] atol=ϵₐ rtol=ϵᵣ
                                i += 1
                            end
                        end
                    end
                end
            end
        end
    end
end
