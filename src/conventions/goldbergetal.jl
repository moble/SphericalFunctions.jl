@doc raw"""
Formulas and conventions from [Goldberg et al.'s "Spin-``s`` Spherical Harmonics and
``\eth``"](@cite GoldbergEtAl_1967).

The conclusion here is that Goldberg et al.'s ``ₛYₗₘ(θ, ϕ)`` differs from ours by a
factor of ``(-1)^m``, while Goldberg et al.'s ``Dʲₘₚ,ₘ`` differs from ours by a transpose
and a factor of ``(-1)^{m+m'}``.

"""
module GoldbergEtAl

const 𝒾 = im


@doc raw"""
    D(j, m′, m, α, β, γ)

Eq. (3.9) of [Goldberg et al.](@cite GoldbergEtAl_1967),
implementing
```math
    D^j_{m',m}(\alpha, \beta, \gamma).
```
"""
function D(j, m′, m, α, β, γ)
    if j < 0
        throw(DomainError("The degree j=$j must be non-negative."))
    end
    if abs(m′) > j
        throw(DomainError("The order abs(m′)=$m′ must be ≤ j=$j."))
    end
    if abs(m) > j
        throw(DomainError("The order abs(m)=$m must be ≤ j=$j."))
    end

    α, β, γ = promote(α, β, γ)

    # The summation index `r` ranges over all values for which the binomials are
    # positive.
    rₘᵢₙ = max(0, m+m′)
    rₘₐₓ = min(j+m′, j+m)

    sinβ╱2 = sin(β/2)
    T = typeof(complex(sinβ╱2))

    √(factorial(j+m) * factorial(j-m) / T(factorial(j+m′) * factorial(j-m′))) *
    sum(
        r -> (
            binomial(j+m′, r) * binomial(j-m′, r-m-m′) * (-1)^(j+m′-r)
            * exp(𝒾*m*α) * cos(β/2)^(2r-m-m′) * sinβ╱2^(2j-2r+m+m′) * exp(𝒾*m′*γ)
        ),
        rₘᵢₙ:rₘₐₓ,
        init=zero(T)
    )
end


@doc raw"""
    Y(s, ℓ, m, θ, ϕ)

Eq. (3.1) of [Goldberg et al.](@cite GoldbergEtAl_1967),
implementing
```math
    {}_sY_{\ell,m}(\theta, \phi).
```
"""
function Y(s, ℓ, m, θ, ϕ)
    if ℓ < 0
        throw(DomainError("The degree ℓ=$ℓ must be non-negative."))
    end
    if abs(m) > ℓ
        throw(DomainError("The order abs(m)=$m must be ≤ ℓ=$ℓ."))
    end
    if abs(s) > ℓ
        throw(DomainError("The spin abs(s)=$s must be ≤ ℓ=$ℓ."))
    end

    θ, ϕ = promote(θ, ϕ)

    # The summation index `r` ranges over all values for which the binomials are
    # positive.
    rₘᵢₙ = max(0, m-s)
    rₘₐₓ = min(ℓ-s, ℓ+m)

    sinθ╱2 = sin(θ/2)
    T = typeof(complex(sinθ╱2))

    √(factorial(ℓ+m) * factorial(ℓ-m) * (2ℓ+1) / (factorial(ℓ+s) * factorial(ℓ-s) * 4T(π))) *
    sum(
        r -> (
            binomial(ℓ-s, r) * binomial(ℓ+s, r+s-m) * (-1)^(ℓ-r-s)
            * exp(𝒾*m*ϕ) * cos(θ/2)^(2r+s-m) * sinθ╱2^(2ℓ-2r-s+m)
        ),
        rₘᵢₙ:rₘₐₓ,
        init=zero(T)
    )
end


end # module GoldbergEtAl


@testitem "GoldbergEtAl conventions" setup=[Utilities] begin
    using Random
    import SphericalFunctions: GoldbergEtAl
    using Quaternionic: from_spherical_coordinates

    Random.seed!(1234)
    const T = Float64
    const ℓₘₐₓ = 8
    ϵₐ = nextfloat(T(0), 4)
    ϵᵣ = 50eps(T)

    # Tests for Y(ℓ, m, θ, ϕ)
    const Y = GoldbergEtAl.Y
    for θ ∈ βrange(T)
        for ϕ ∈ αrange(T)

            # Test Eq. (2.6) of [Goldberg et al.](@cite GoldbergEtAl_1967)
            for ℓ ∈ 0:ℓₘₐₓ
                for s ∈ -ℓ:-1
                    for m ∈ -ℓ:ℓ
                        @test conj(Y(s, ℓ, m, θ, ϕ)) ≈ (-1)^(m+s) * Y(-s, ℓ, -m, θ, ϕ)
                    end
                end
            end

            # Compare to SphericalHarmonics Y
            for s ∈ -ℓₘₐₓ:ℓₘₐₓ
                Y₁ = ₛ𝐘(s, ℓₘₐₓ, T, [from_spherical_coordinates(θ, ϕ)])[1,:]
                Y₂ = [(-1)^m * Y(s, ℓ, m, θ, ϕ) for ℓ ∈ abs(s):ℓₘₐₓ for m ∈ -ℓ:ℓ]
                @test Y₁ ≈ Y₂ atol=ϵₐ rtol=ϵᵣ
            end
        end
    end

    # Tests for D(j, m′, m, α, β, γ)
    let ℓₘₐₓ=6, ϵₐ=√ϵᵣ, ϵᵣ=√ϵᵣ, 𝒟=GoldbergEtAl.D
        for α ∈ αrange(T)
            for β ∈ βrange(T)
                for γ ∈ γrange(T)
                    D = D_matrices(α, β, γ, ℓₘₐₓ)
                    i = 1
                    for j in 0:ℓₘₐₓ
                        for m′ in -j:j
                            for m in -j:j
                                @test (-1)^(m+m′) * 𝒟(j, m, m′, α, β, γ) ≈ D[i] atol=ϵₐ rtol=ϵᵣ
                                i += 1
                            end
                        end
                    end
                end
            end
        end
    end
end
