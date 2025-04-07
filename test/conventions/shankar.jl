@testmodule Shankar begin

const 𝒾 = im

include("../utilities/naive_factorial.jl")
import .NaiveFactorials: ❗


# Shankar's explicit formulas from Eq. (12.5.39)
Y₀⁰(θ, ϕ) = 1 / √(4π)
Y₁⁻¹(θ, ϕ) = +√(3/(8π)) * sin(θ) * exp(-𝒾*ϕ)
Y₁⁰(θ, ϕ) = √(3/(4π)) * cos(θ)
Y₁⁺¹(θ, ϕ) = -√(3/(8π)) * sin(θ) * exp(+𝒾*ϕ)
Y₂⁻²(θ, ϕ) = √(15/(32π)) * sin(θ)^2 * exp(-2𝒾*ϕ)
Y₂⁻¹(θ, ϕ) = +√(15/(8π)) * sin(θ) * cos(θ) * exp(-𝒾*ϕ)
Y₂⁰(θ, ϕ) = √(5/(16π)) * (3cos(θ)^2 - 1)
Y₂⁺¹(θ, ϕ) = -√(15/(8π)) * sin(θ) * cos(θ) * exp(+𝒾*ϕ)
Y₂⁺²(θ, ϕ) = √(15/(32π)) * sin(θ)^2 * exp(+2𝒾*ϕ)

end  # @testmodule Shankar
