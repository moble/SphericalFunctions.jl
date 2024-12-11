raw"""

Formulas and conventions from [Thorne (1980)](@cite Thorne_1980).

"""
@testmodule Thorne begin

const 𝒾 = im

include("../utilities/naive_factorial.jl")
import .NaiveFactorials: ❗


# Eq. (2.8) upper
function C(ℓ, m, T)
    let π=convert(T, π), √=sqrt∘T
        (-1)^m * √T(
            ((2ℓ+1) * (ℓ-m)❗)
            / (4π * (ℓ+m)❗)
        )
    end
end


# Eq. (2.8) lower
function a(ℓ, m, j, T)
    T(((-1)^j / (2^ℓ * (j)❗ * (ℓ-j)❗)) * ((2ℓ-2j)❗ / (ℓ-m-2j)❗))
end


@doc raw"""
    Y(ℓ, m, θ, ϕ)

Eqs. (2.7) of [Thorne](@cite Thorne_1980), implementing
```math
    Y^{\ell,m}(\theta, \phi).
```
"""
function Y(ℓ, m, θ, ϕ)
    if m < 0
        return (-1)^m * conj(Y(ℓ, abs(m), θ, ϕ))
    end
    θ, ϕ = promote(θ, ϕ)
    sinθ, cosθ = sincos(θ)
    T = typeof(sinθ)
    C(ℓ, m, T) * (exp(𝒾*ϕ) * sinθ)^m * sum(
        j -> a(ℓ, m, j, T) * (cosθ)^(ℓ-m-2j),
        0:floor((ℓ-m)÷2),
        init=zero(T)
    )
end

end  # @testmodule Thorne


@testitem "Thorne conventions" setup=[Utilities, Thorne] begin
    using Random
    using Quaternionic: from_spherical_coordinates

    Random.seed!(1234)
    const T = Float64
    const ℓₘₐₓ = 5
    ϵₐ = nextfloat(T(0), 4)
    ϵᵣ = 20eps(T)

    # Tests for Y(ℓ, m, θ, ϕ)
    for θ ∈ βrange(T)
        for ϕ ∈ αrange(T)

            # Test Thorne's Eq. (2.9b)
            for ℓ in 0:ℓₘₐₓ
                for m in -ℓ:-1
                    @test conj(Thorne.Y(ℓ, m, θ, ϕ)) ≈ (-1)^-m * Thorne.Y(ℓ, -m, θ, ϕ)
                end
            end

            # Compare to SphericalFunctions
            let s=0
                Y = ₛ𝐘(s, ℓₘₐₓ, T, [from_spherical_coordinates(θ, ϕ)])
                i = 1
                for ℓ in 0:ℓₘₐₓ
                    for m in -ℓ:ℓ
                        @test Thorne.Y(ℓ, m, θ, ϕ) ≈ Y[i]
                        i += 1
                    end
                end
            end
        end
    end
end
