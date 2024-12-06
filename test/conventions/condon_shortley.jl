raw"""
Formulas and conventions from [Condon and Shortley's "The Theory Of Atomic Spectra"](@cite
CondonShortley_1935).

The method we use here is as direct and explicit as possible.  In particular, Condon and
Shortley provide a formula for the φ=0 part in terms of iterated derivatives of a power of
sin(θ).  Rather than expressing these derivatives in terms of the Legendre polynomials —
which would subject us to another round of ambiguity — the functions in this module use
automatic differentiation to compute the derivatives explicitly.

The result is that the original Condon-Shortley spherical harmonics agree perfectly with the
ones computed by this package.

(Condon and Shortley do not give an expression for the Wigner D-matrices.)

"""
@testmodule CondonShortley begin

import FastDifferentiation

const 𝒾 = im

include("../utilities/naive_factorial.jl")
import .NaiveFactorials: ❗


"""
    Θ(ℓ, m, θ)

Equation (15) of section 4³ (page 52) of [Condon-Shortley](@cite CondonShortley_1935),
implementing
```math
    Θ(ℓ, m),
```
which is implicitly a function of the spherical coordinate ``θ``.
"""
function Θ(ℓ, m, θ::T) where {T}
    (-1)^ℓ * T(√(((2ℓ+1) * (ℓ+m)❗) / (2 * (ℓ - m)❗)) * (1 / (2^ℓ * (ℓ)❗))) *
    (1 / sin(θ)^T(m)) * dʲsin²ᵏθdcosθʲ(ℓ-m, ℓ, θ)
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


"""
    Φ(mₗ, φ)

Equation (5) of section 4³ (page 50) of [Condon-Shortley](@cite CondonShortley_1935),
implementing
```math
    Φ(mₗ),
```
which is implicitly a function of the spherical coordinate ``φ``.
"""
function Φ(mₗ, φ::T) where {T}
    1 / √(2T(π)) * exp(𝒾 * mₗ * φ)
end


"""
    ϕ(ℓ, m, θ, φ)

Spherical harmonics.  This is defined as such below Eq. (5) of section 5⁵ (page 127) of
[Condon-Shortley](@cite CondonShortley_1935), implementing
```math
    ϕ(ℓ, mₗ),
```
which is implicitly a function of the spherical coordinates ``θ`` and ``φ``.
"""
function ϕ(ℓ, mₗ, θ, φ)
    Θ(ℓ, mₗ, θ) * Φ(mₗ, φ)
end

@doc raw"""
    ϴ(ℓ, m, θ)

Explicit formulas for the first few spherical harmonics as given by Condon-Shortley in the
footnote to Eq. (15) of Sec. 4³ (page 52).

Note that the name of this function is `\varTheta`, as opposed to the `\Theta` function
that implements Condon-Shortley's general form.
"""
ϴ(ℓ, m, θ) = ϴ(Val(ℓ), Val(m), θ)
ϴ(::Val{0}, ::Val{0}, θ) = √(1/2)
ϴ(::Val{1}, ::Val{0}, θ) = √(3/2) * cos(θ)
ϴ(::Val{2}, ::Val{0}, θ) = √(5/8) * (2cos(θ)^2 - sin(θ)^2)
ϴ(::Val{3}, ::Val{0}, θ) = √(7/8) * (2cos(θ)^3 - 3cos(θ)sin(θ)^2)
ϴ(::Val{1}, ::Val{+1}, θ) = -√(3/4) * sin(θ)
ϴ(::Val{1}, ::Val{-1}, θ) = +√(3/4) * sin(θ)
ϴ(::Val{2}, ::Val{+1}, θ) = -√(15/4) * cos(θ) * sin(θ)
ϴ(::Val{2}, ::Val{-1}, θ) = +√(15/4) * cos(θ) * sin(θ)
ϴ(::Val{3}, ::Val{+1}, θ) = -√(21/32) * (4cos(θ)^2*sin(θ) - sin(θ)^3)
ϴ(::Val{3}, ::Val{-1}, θ) = +√(21/32) * (4cos(θ)^2*sin(θ) - sin(θ)^3)
ϴ(::Val{2}, ::Val{+2}, θ) = √(15/16) * sin(θ)^2
ϴ(::Val{2}, ::Val{-2}, θ) = √(15/16) * sin(θ)^2
ϴ(::Val{3}, ::Val{+2}, θ) = √(105/16) * cos(θ) * sin(θ)^2
ϴ(::Val{3}, ::Val{-2}, θ) = √(105/16) * cos(θ) * sin(θ)^2
ϴ(::Val{3}, ::Val{+3}, θ) = -√(35/32) * sin(θ)^3
ϴ(::Val{3}, ::Val{-3}, θ) = +√(35/32) * sin(θ)^3

end  # @testmodule CondonShortley


@testitem "Condon-Shortley conventions" setup=[Utilities, CondonShortley] begin
    using Random
    using Quaternionic: from_spherical_coordinates
    #const check = NaNChecker.NaNCheck

    Random.seed!(1234)
    const T = Float64
    const ℓₘₐₓ = 4
    ϵₐ = 4eps(T)
    ϵᵣ = 1000eps(T)

    # Tests for Y(ℓ, m, θ, ϕ)
    let Y=CondonShortley.ϕ, Θ=CondonShortley.Θ, ϴ=CondonShortley.ϴ, ϕ=zero(T)
        for θ ∈ βrange(T)
            if abs(sin(θ)) < ϵₐ
                continue
            end

            # # Find where NaNs are coming from
            # for ℓ ∈ 0:ℓₘₐₓ
            #     for m ∈ -ℓ:ℓ
            #         Θ(ℓ,  m, check(θ))
            #     end
            # end

            # Test footnote to Eq. (15) of Sec. 4³ of Condon-Shortley
            let Y = ₛ𝐘(0, 3, T, [from_spherical_coordinates(θ, ϕ)])[1,:]
                for ℓ ∈ 0:3
                    for m ∈ -ℓ:ℓ
                        @test ϴ(ℓ, m, θ) / √(2π) ≈ Y[Yindex(ℓ, m)] atol=ϵₐ rtol=ϵᵣ
                    end
                end
            end

            # Test Eq. (18) of Sec. 4³ of Condon-Shortley
            for ℓ ∈ 0:ℓₘₐₓ
                for m ∈ -ℓ:ℓ
                    @test Θ(ℓ, m, θ) ≈ (-1)^(m) * Θ(ℓ, -m, θ) atol=ϵₐ rtol=ϵᵣ
                end
            end

            # Compare to SphericalHarmonics Y
            let s = 0
                Y₁ = ₛ𝐘(s, ℓₘₐₓ, T, [from_spherical_coordinates(θ, ϕ)])[1,:]
                Y₂ = [Y(ℓ, m, θ, ϕ) for ℓ ∈ abs(s):ℓₘₐₓ for m ∈ -ℓ:ℓ]
                @test Y₁ ≈ Y₂ atol=ϵₐ rtol=ϵᵣ
            end
        end
    end

end
