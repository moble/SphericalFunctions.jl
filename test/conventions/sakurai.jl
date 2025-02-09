raw"""
Formulas and conventions from [Sakurai's "Modern Quantum Mechanics"](@cite Sakurai_1994).

The conclusion here is that Sakurai's Yₗᵐ(θ, ϕ) is the same as ours, but his
𝒟ʲₘₚ,ₘ is conjugated relative to ours.

- On p. 154 he says that "a rotation operation affects the physical system itself, ...,
  while the coordinate axes remain *unchanged*."
- On p. 156 he poses "``|\alpha\rangle_R = \mathcal{D}(R) | \alpha \rangle``, where
  ``|\alpha\rangle_R`` and ``|\alpha \rangle`` stand for the kets of the rotated and
  original system, respectively."
- On p. 157 he says "``\mathcal{D}(\hat{\mathbf{n}}, d\phi) = 1 - i\left( \frac{\mathbf{J}
  \cdot \hat{\mathbf{n}}} {\hbar} \right) d\phi``"
- On p. 173 he defines his Euler angles in the same way as Quaternionic.
- On p. 192 he defines "``\mathcal{D}^{(j)}_{m',m}(R) =  \langle j,m'| \exp \left(
  \frac{-i\mathbf{J} \cdot \hat{\mathbf{n}} \phi} {\hbar} \right) |j, m\rangle``".
- On p. 194 he gives the expression in terms of Euler angles.
- On p. 223 he gives an explicit formula for ``d``.
- On p. 203 he relates ``\mathcal{D} to Y_{\ell}^m$ (note the upper index of ``m``).


Below (1.6.14), we find the translation operator acts as
``\mathscr{T}_{dx'} \alpha(x') = \alpha(x' - dx')``.  Then Eq.
(1.6.32)
```math
\mathscr{T}_{dx'} = 1 - i p\, dx',
```
for infinitesimal ``dx'``.  Eq. (1.7.17) gives the momentum operator
as ``p \alpha(x') = -i \partial_{x'} \alpha(x')``.  Combining these,
we can verify consistency:
```math
\mathscr{T}_{dx'} \alpha(x')
=
\alpha(x' - dx')
=
\alpha(x') - \partial_{x'}\, \alpha(x')\, dx',
```
which is exactly what we expect from Taylor expanding ``\alpha(x' -
dx')``.


```math
\begin{aligned}
f\left(𝐑\right)
&\to
f\left(e^{-\epsilon 𝐮/2}𝐑\right) \\
&\approx
f\left(𝐑\right) + \epsilon \left. \frac{d}{d\epsilon} \right|_{\epsilon=0}
f\left(e^{-\epsilon 𝐮/2}𝐑\right) \\
&=
f\left(𝐑\right) - i \epsilon L_𝐮 f\left(𝐑\right)
```

"""
@testmodule Sakurai begin

const 𝒾 = im

include("../utilities/naive_factorial.jl")
import .NaiveFactorials: ❗


@doc raw"""
    𝒟(j, m′, m, α, β, γ)

Eqs. (3.5.50)-(3.5.51) of [Sakurai](@cite Sakurai_1994), p. 194,
implementing
```math
    \mathcal{D}^{(j)}_{m',m}(\alpha, \beta, \gamma).
```

See also [`d`](@ref) for Sakurai's version the Wigner d-function.
"""
function 𝒟(j, m′, m, α, β, γ)
    exp(-𝒾*(m′*α + m*γ)) * d(j, m′, m, β)
end

@doc raw"""
    d(j, m′, m, β)

Eqs. (3.5.50)-(3.5.51) of [Sakurai](@cite Sakurai_1994), p. 194
(or Eq. (3.8.33), p. 223), implementing
```math
    d^{(j)}_{m',m}(\beta).
```

See also [`𝒟`](@ref) for Sakurai's version the Wigner D-function.
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
    kₘᵢₙ = max(0, m-m′)
    kₘₐₓ = min(j-m′, j+m)

    T = typeof(β)

    # Note that Sakurai's actual formula is reproduced here, even though it leads to
    # overflow errors for `j ≥ 8`, which could be eliminated by other means.
    return sum(
        k -> (
            (-1)^(k-m+m′) * T(
                √((j+m)❗ * (j-m)❗ * (j+m′)❗ * (j-m′)❗)
                / ((j+m-k)❗ * (k)❗ * (j-k-m′)❗ * (k-m+m′)❗)
            )
            * cos(β/2)^(2j-2k+m-m′) * sin(β/2)^(2k-m+m′)
        ),
        kₘᵢₙ:kₘₐₓ,
        init=zero(T)
    )
end

@doc raw"""
    Y(ℓ, m, θ, ϕ)

Eqs. (3.6.51) of [Sakurai](@cite Sakurai_1994), p. 203,
implementing
```math
    Y_{\ell}^m(\theta, \phi).
```
"""
function Y(ℓ, m, θ, ϕ)
    conj(√((2ℓ+1)/(4π)) * 𝒟(ℓ, m, 0, ϕ, θ, 0))
end

# Sakurai's explicit formulas from Eq. (A.5.7e-i) of [Sakurai](@cite Sakurai_1994), p. 451
Y₀⁰(θ, ϕ) = 1 / √(4π)
Y₁⁻¹(θ, ϕ) = +√(3/(8π)) * sin(θ) * exp(-𝒾*ϕ)
Y₁⁰(θ, ϕ) = √(3/(4π)) * cos(θ)
Y₁⁺¹(θ, ϕ) = -√(3/(8π)) * sin(θ) * exp(+𝒾*ϕ)
Y₂⁻²(θ, ϕ) = √(15/(32π)) * sin(θ)^2 * exp(-2𝒾*ϕ)
Y₂⁻¹(θ, ϕ) = +√(15/(8π)) * sin(θ) * cos(θ) * exp(-𝒾*ϕ)
Y₂⁰(θ, ϕ) = √(5/(16π)) * (3cos(θ)^2 - 1)
Y₂⁺¹(θ, ϕ) = -√(15/(8π)) * sin(θ) * cos(θ) * exp(+𝒾*ϕ)
Y₂⁺²(θ, ϕ) = √(15/(32π)) * sin(θ)^2 * exp(+2𝒾*ϕ)

end  # @testmodule Sakurai


@testitem "Sakurai conventions" setup=[Utilities, Sakurai] begin
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
            # Test Sakurai's own explicit formulas
            @test Sakurai.Y₀⁰(θ, ϕ) ≈ Sakurai.Y(0, 0, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
            @test Sakurai.Y₁⁻¹(θ, ϕ) ≈ Sakurai.Y(1, -1, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
            @test Sakurai.Y₁⁰(θ, ϕ) ≈ Sakurai.Y(1, 0, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
            @test Sakurai.Y₁⁺¹(θ, ϕ) ≈ Sakurai.Y(1, 1, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
            @test Sakurai.Y₂⁻²(θ, ϕ) ≈ Sakurai.Y(2, -2, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
            @test Sakurai.Y₂⁻¹(θ, ϕ) ≈ Sakurai.Y(2, -1, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
            @test Sakurai.Y₂⁰(θ, ϕ) ≈ Sakurai.Y(2, 0, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
            @test Sakurai.Y₂⁺¹(θ, ϕ) ≈ Sakurai.Y(2, 1, θ, ϕ) atol=ϵₐ rtol=ϵᵣ
            @test Sakurai.Y₂⁺²(θ, ϕ) ≈ Sakurai.Y(2, 2, θ, ϕ) atol=ϵₐ rtol=ϵᵣ

            # Test Sakurai's Eq. (A.5.7b) of [Sakurai](@cite Sakurai_1994), p. 451
            for ℓ in 0:ℓₘₐₓ
                for m in -ℓ:-1
                    @test Sakurai.Y(ℓ, m, θ, ϕ) ≈ (-1)^-m * conj(Sakurai.Y(ℓ, -m, θ, ϕ))
                end
            end

            # Compare to SphericalFunctions
            let s=0
                Y = ₛ𝐘(s, ℓₘₐₓ, T, [from_spherical_coordinates(θ, ϕ)])
                i = 1
                for ℓ in 0:ℓₘₐₓ
                    for m in -ℓ:ℓ
                        @test Sakurai.Y(ℓ, m, θ, ϕ) ≈ Y[i]
                        i += 1
                    end
                end
            end
        end
    end

    # Tests for 𝒟(j, m′, m, α, β, γ)
    let ϵₐ=√ϵᵣ, ϵᵣ=√ϵᵣ, 𝒟=Sakurai.𝒟
        for α ∈ αrange(T)
            for β ∈ βrange(T)
                for γ ∈ γrange(T)
                    D = D_matrices(α, β, γ, ℓₘₐₓ)
                    i = 1
                    for j in 0:ℓₘₐₓ
                        for m′ in -j:j
                            for m in -j:j
                                @test 𝒟(j, m′, m, α, β, γ) ≈ conj(D[i]) atol=ϵₐ rtol=ϵᵣ
                                i += 1
                            end
                        end
                    end
                end
            end
        end
    end

end
