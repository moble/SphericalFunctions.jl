raw"""
Formulas and conventions from [Torres del Castillo's "3-D spinors, spin-weighted functions
and their applications"](@cite TorresDelCastillo_2003).

The conclusion here is that del Castillo's ``ₛYₗₘ(θ, ϕ)`` is identical to ours, while
his ``Dʲₘₚ,ₘ`` is conjugated relative to ours.

"""
@testmodule TorresDelCastillo begin

const 𝒾 = im

include("../utilities/naive_factorial.jl")
import .NaiveFactorials: ❗


@doc raw"""
    D(j, m′, m, ϕ, θ, χ)

Eq. (2.52) of [Torres del Castillo](@cite TorresDelCastillo_2003), implementing
```math
    D^l_{m',m}(ϕ, θ, \chi).
```
"""
function D(l, m′, m, ϕ, θ, χ)
    if l < 0
        throw(DomainError("The degree l=$l must be non-negative."))
    end
    if abs(m′) > l
        throw(DomainError("The order abs(m′)=$m′ must be ≤ l=$l."))
    end
    if abs(m) > l
        throw(DomainError("The order abs(m)=$m must be ≤ l=$l."))
    end

    ₘYₗ₋ₘₚ = Y(m, l, -m′, θ, ϕ)

    let π = oftype(ₘYₗ₋ₘₚ, π)
        (-1)^m′ * √((4π)/(2l+1)) * ₘYₗ₋ₘₚ * exp(-𝒾*m*χ)
    end
end


@doc raw"""
    Y(s, j, m, θ, ϕ)

The equation following Eq. (2.53) of [Torres del Castillo](@cite TorresDelCastillo_2003),
implementing
```math
    {}_sY_{j,m}(θ, ϕ).
```
"""
function Y(s, j, m, θ, ϕ)
    if j < 0
        throw(DomainError("The degree j=$j must be non-negative."))
    end
    if abs(m) > j
        throw(DomainError("The order abs(m)=$m must be ≤ j=$j."))
    end
    if abs(s) > j
        throw(DomainError("The spin abs(s)=$s must be ≤ j=$j."))
    end

    dʲ₋ₘₛ = d(j, -m, s, θ)

    let π = oftype(dʲ₋ₘₛ, π)
        (-1)^m * √((2j+1)/(4π)) * dʲ₋ₘₛ * exp(𝒾*m*ϕ)
    end
end


@doc raw"""
    d(l, m′, m, θ)

Second equation below Eq. (2.53) of [Torres del Castillo](@cite TorresDelCastillo_2003),
implementing
```math
    d^l_{m',m}(θ).
```
"""
function d(l, m′, m, θ)
    if l < 0
        throw(DomainError("The degree l=$l must be non-negative."))
    end
    if abs(m′) > l
        throw(DomainError("The order abs(m′)=$m′ must be ≤ l=$l."))
    end
    if abs(m) > l
        throw(DomainError("The order abs(m)=$m must be ≤ l=$l."))
    end

    # The summation index `k` ranges over all values for which the factorial arguments are
    # valid.
    kₘᵢₙ = max(0, -m+m′)
    kₘₐₓ = min(l+m′, l-m)

    sinθ╱2, cosθ╱2 = sincos(θ/2)
    T = typeof(sinθ╱2)

    √T((l+m)❗ * (l-m)❗ * (l+m′)❗ * (l-m′)❗) *
    sum(
        k -> (
            (-1)^(k) * sinθ╱2^(m-m′+2k) * cosθ╱2^(2l-m+m′-2k) /
            T((k)❗ * (l+m′-k)❗ * (l-m-k)❗ * (m-m′+k)❗)
        ),
        kₘᵢₙ:kₘₐₓ,
        init=complex(zero(T))
    )
end


end  # @testmodule TorresDelCastillo


@testitem "TorresDelCastillo conventions" setup=[Utilities, TorresDelCastillo] begin
    using Random
    using Quaternionic: from_spherical_coordinates
    using SphericalFunctions: Deprecated

    Random.seed!(1234)
    const T = Float64
    const ℓₘₐₓ = 5
    ϵₐ = 2eps(T)
    ϵᵣ = 50eps(T)

    # Tests for Y(ℓ, m, θ, ϕ)
    const Y = TorresDelCastillo.Y
    for θ ∈ βrange(T)
        for ϕ ∈ αrange(T)
            for s ∈ -ℓₘₐₓ:ℓₘₐₓ
                Y₁ = Deprecated.ₛ𝐘(s, ℓₘₐₓ, T, [from_spherical_coordinates(θ, ϕ)])[1,:]
                Y₂ = [Y(s, ℓ, m, θ, ϕ) for ℓ ∈ abs(s):ℓₘₐₓ for m ∈ -ℓ:ℓ]
                @test Y₁ ≈ Y₂ atol=ϵₐ rtol=ϵᵣ
            end
        end
    end

    # Tests for D(j, m′, m, α, β, γ)
    let ϵₐ=√ϵᵣ, ϵᵣ=√ϵᵣ, 𝒟=TorresDelCastillo.D
        for α ∈ αrange(T)
            for β ∈ βrange(T)
                for γ ∈ γrange(T)
                    D = Deprecated.D_matrices(α, β, γ, ℓₘₐₓ)
                    i = 1
                    for j in 0:ℓₘₐₓ
                        for m′ in -j:j
                            for m in -j:j
                                @test conj(𝒟(j, m′, m, α, β, γ)) ≈ D[i] atol=ϵₐ rtol=ϵᵣ
                                i += 1
                            end
                        end
                    end
                end
            end
        end
    end
end
