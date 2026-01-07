raw"""
Formulas and conventions from [Edmonds' "Angular Momentum in Quantum Mechanics"](@cite
Edmonds_2016).

Note that Edmonds explains on page 8 that his Euler angles agree with ours.  His spherical
harmonics agree also, but his ``𝔇`` is transposed.  Alternatively, we could think of his
``𝔇`` being conjugated — just like other modern conventions — but taking the inverse
rotation as argument.

TODO: Figure out the meaning of those rotations.

"""
@testmodule Edmonds begin

import FastDifferentiation

const 𝒾 = im

include("../utilities/naive_factorial.jl")
import .NaiveFactorials: ❗


@doc raw"""
    Y(ℓ, m, θ, φ)

Eq. (2.5.5) of [Edmonds](@cite Edmonds_2016), implementing
```math
    Yₗₘ(θ, φ).
```
"""
function Y(ℓ, m, θ::T, φ::T)::Complex{T} where {T}
    (-1)^(ℓ+m) / (2^ℓ * (ℓ)❗) * √((2ℓ+1)*(ℓ-m)❗/(4big(π) * (ℓ+m)❗)) *
    (sin(θ)^T(m)) * dʲsin²ᵏθdcosθʲ(ℓ+m, ℓ, θ) * exp(𝒾*m*φ)
end


@doc raw"""
    𝒟(j, m′, m, α, β, γ)

Eqs. (4.1.12) of [Edmonds](@cite Edmonds_2016), implementing
```math
    𝒟^{(j)}_{m',m}(α, β, γ).
```

See also [`d`](@ref) for Edmonds' version the Wigner d-function.
"""
function 𝒟(j, m′, m, α, β, γ)
    exp(𝒾*m′*γ) * d(j, m′, m, β) * exp(𝒾*m*α)
end


@doc raw"""
    d(j, m′, m, β)

Eqs. (4.1.15) of [Edmonds](@cite Edmonds_2016), implementing
```math
    d^{(j)}_{m',m}(β).
```

See also [`𝒟`](@ref) for Edmonds' version the Wigner D-function.
"""
function d(j, m′, m, β)
    if j < 0
        throw(DomainError("j=$j must be non-negative"))
    end
    if abs(m′) > j || abs(m) > j
        throw(DomainError("abs(m′=$m′) and abs(m=$m) must be ≤ j=$j"))
    end
    if j ≥ 8
        throw(DomainError("j=$j≥8 will lead to overflow errors"))
    end

    # The summation index `k` ranges over all values for which the factorials are
    # non-negative.
    σₘᵢₙ = 0
    σₘₐₓ = j - m′

    T = typeof(β)

    # Note that Edmonds' actual formula is reproduced here, even though it leads to overflow
    # errors for `j ≥ 8`, which could be eliminated by other means.
    return √T((j+m′)❗ * (j-m′)❗ / ((j+m)❗ * (j-m)❗)) *
    sum(
        σ -> (
            binomial(j+m, j-m′-σ) * binomial(j-m, σ) *
            (-1)^(j-m′-σ) * cos(β/2)^(2σ+m′+m) * sin(β/2)^(2j-2σ-m′-m)
        ),
        σₘᵢₙ:σₘₐₓ,
        init=zero(T)
    )
end


@doc raw"""
    dʲsin²ᵏθdcosθʲ(j, k, θ)

Compute the ``j``th derivative of the function ``\sin^{2k}(θ)`` with respect to ``\cos(θ)``.
Note that ``\sin^{2k}(θ) = (1 - \cos^2(θ))^k``, so this is equivalent to evaluating the
``j``th derivative of ``(1-x^2)^k`` with respect to ``x``, evaluated at ``x = \cos(θ)``.
"""
function dʲsin²ᵏθdcosθʲ(j, k, θ)
    if j < 0
        throw(ArgumentError("j=$j must be non-negative"))
    end
    if j == 0
        return sin(θ)^(2k)
    end
    x = FastDifferentiation.make_variables(:x)[1]
    ∂ₓʲfᵏ = FastDifferentiation.derivative((1 - x^2)^k, (x for _ ∈ 1:j)...)
    return FastDifferentiation.make_function([∂ₓʲfᵏ,], [x,])(cos(θ))[1]
end


end  # @testmodule Edmonds


@testitem "Edmonds conventions" setup=[Utilities, Edmonds] begin
    using Random
    using Quaternionic: from_spherical_coordinates
    using SphericalFunctions: Deprecated

    Random.seed!(1234)
    const T = Float64
    const ℓₘₐₓ = 3
    ϵₐ = 8eps(T)
    ϵᵣ = 20eps(T)

    # Tests for Y(ℓ, m, θ, ϕ)
    for θ ∈ βrange(T, 3)
        if abs(sin(θ)) ≤ eps(T)
            continue
        end

        for ϕ ∈ αrange(T, 3)
            # Test Edmonds' Eq. (2.5.5)
            let Y = Edmonds.Y
                for ℓ in 0:ℓₘₐₓ
                    for m in -ℓ:0
                        @test Y(ℓ, -m, θ, ϕ) ≈ (-1)^-m * conj(Y(ℓ, m, θ, ϕ)) atol=ϵₐ rtol=ϵᵣ
                    end
                end
            end

            # Compare to SphericalFunctions
            let s=0
                Y = Deprecated.ₛ𝐘(s, ℓₘₐₓ, T, [from_spherical_coordinates(θ, ϕ)])
                i = 1
                for ℓ in 0:ℓₘₐₓ
                    for m in -ℓ:ℓ
                        @test Edmonds.Y(ℓ, m, θ, ϕ) ≈ Y[i] atol=ϵₐ rtol=ϵᵣ
                        i += 1
                    end
                end
            end
        end
    end

    # Tests for 𝒟(j, m′, m, α, β, γ)
    let ϵₐ=√ϵᵣ, ϵᵣ=√ϵᵣ, 𝒟=Edmonds.𝒟
        for α ∈ αrange(T)
            for β ∈ βrange(T)
                if abs(sin(β)) ≤ eps(T)
                    continue
                end

                for γ ∈ γrange(T)
                    D = Deprecated.D_matrices(α, β, γ, ℓₘₐₓ)
                    i = 1
                    for j in 0:ℓₘₐₓ
                        for m′ in -j:j
                            for m in -j:j
                                #@test 𝒟(j, m, m′, α, β, γ) ≈ D[i] atol=ϵₐ rtol=ϵᵣ
                                @test 𝒟(j, m′, m, -γ, -β, -α) ≈ conj(D[i]) atol=ϵₐ rtol=ϵᵣ
                                i += 1
                            end
                        end
                    end
                end
            end
        end
    end

end
