@testmodule Boyle2016 begin

    using Quaternionic

    """
        WignerDElement(R, ℓ, m′, m)

    Compute a single Wigner-D matrix element for half-integer or integer (ℓ, m′, m).

    `R` is a `Rotor`, and `ℓ`, `m′`, and `m` are the indices of the Wigner-D matrix element.
    The indices must all be integers or all be `Rational` with denominators of 2.

    """
    function WignerDElement(R::Rotor{T}, ℓ::I, m′::I, m::I) where {T, I}
        # If `I` is Rational, check that the denominators are 2
        if I <: Rational
            if (denominator(ℓ) != 2) || (denominator(m′) != 2) || (denominator(m) != 2)
                error("The indices ℓ, m′, and m must all be integers or all be half-integers")
            end
        end

        # Convert to twice the input values for half-integer support
        L  = Int(2ℓ)
        M′ = Int(2m′)
        M  = Int(2m)

        if L > 16
            error(
                "The maximum supported ℓ for this function is 8; " *
                "larger numbers become numerically unstable.\n" *
                "Consider using the `WignerD` function instead."
            )
        end

        # Simple helper for 0+0im
        zeroCT = zero(Complex{T})

        if abs(M′) > L || abs(M) > L
            return zeroCT
        end

        let π = T(π)
            # Split input `R` into its two complex components and extract magnitude and phase
            Rₛ = Complex(R[1], R[4])
            Rₐ = Complex(R[3], R[2])
            rₛ = abs(Rₛ)
            rₐ = abs(Rₐ)
            ϕₛ = angle(Rₛ)
            ϕₐ = angle(Rₐ)

            # Check simple limiting cases
            if rₐ ≤ 4eps(rₛ)
                if M′ != M
                    return zeroCT
                else
                    return cis(M * ϕₛ)
                end

            elseif rₛ ≤ 4eps(rₐ)
                if -M′ != M
                    return zeroCT
                else
                    return cis(M * ϕₐ) * (((L - M) % 4 == 0) ? 1 : -1)
                end

            elseif rₐ ≤ rₛ
                λ = -(rₐ/rₛ)^2
                ρₘᵢₙ = max(0, (M′ - M)÷2)
                κ = √T(
                        (factorial((L + M)÷2) * factorial((L - M)÷2))
                        / (factorial((L + M′)÷2) * factorial((L - M′)÷2))
                    ) *
                    binomial((L + M′)÷2, ρₘᵢₙ) * binomial((L - M′)÷2, (L - M)÷2 - ρₘᵢₙ)
                if (ρₘᵢₙ % 2) != 0
                    κ = -κ
                end
                ρₘₐₓ = min((L + M′)÷2, (L - M)÷2)
                N₁ = L + M′ + 2
                N₂ = L - M + 2
                N₃  = M - M′

                total = one(T)
                for P in reverse(2ρₘᵢₙ+2:2:2ρₘₐₓ)
                    total *= λ * ((N₁ - P)*(N₂ - P)) / (P*(N₃ + P))
                    total += one(T)
                end
                return κ *
                    (rₛ ^ (L - (M - M′)÷2 - 2ρₘᵢₙ)) *
                    (rₐ ^ ((M - M′)÷2 + 2ρₘᵢₙ)) *
                    cis((M + M′)÷2 * ϕₛ + (M - M′)÷2 * ϕₐ) *
                    total

                # return κ *
                #     (rₛ ^ (L - (M - M′)÷2 - 2ρₘᵢₙ)) * (rₐ ^ ((M - M′)÷2 + 2ρₘᵢₙ)) *
                #     cis((M + M′)÷2 * ϕₛ + (M - M′)÷2 * ϕₐ) *
                #     foldl(
                #         (acc, P) -> acc * λ * ((N₁ - P)*(N₂ - P)) / (P*(N₃ + P)) + one(T),
                #         reverse(2ρₘᵢₙ+2:2:2ρₘₐₓ),
                #         init=one(T)
                #     )

            else # rₛ < rₐ
                λ = -(rₛ/rₐ)^2
                ρₘᵢₙ = max(0, -(M′ + M)÷2)
                κ = √T(
                        (factorial((L + M)÷2) * factorial((L - M)÷2))
                        / (factorial((L + M′)÷2) * factorial((L - M′)÷2))
                    ) *
                    binomial((L + M′)÷2, (L - M)÷2 - ρₘᵢₙ) * binomial((L - M′)÷2, ρₘᵢₙ)
                if (((L - M)÷2 - ρₘᵢₙ) % 2) != 0
                    κ = -κ
                end
                ρₘₐₓ = min((L - M′)÷2, (L - M)÷2)
                N₁ = L - M′ + 2
                N₂ = L - M + 2
                N₃  = M + M′

                total = one(T)
                for P in reverse(2ρₘᵢₙ+2:2:2ρₘₐₓ)
                    total *= λ * ((N₁ - P)*(N₂ - P)) / (P*(N₃ + P))
                    total += one(T)
                end
                return κ *
                    (rₐ ^ (L - (M + M′)÷2 - 2ρₘᵢₙ)) *
                    (rₛ ^ ((M + M′)÷2 + 2ρₘᵢₙ)) *
                    cis((M - M′)÷2 * ϕₐ + (M + M′)÷2 * ϕₛ) *
                    total

                # return κ *
                #     (rₐ ^ (L - (M + M′)÷2 - 2ρₘᵢₙ)) *
                #     (rₛ ^ ((M + M′)÷2 + 2ρₘᵢₙ)) *
                #     cis((M - M′)÷2 * ϕₐ + (M + M′)÷2 * ϕₛ) *
                #     foldl(
                #         (acc, P) -> acc * λ * ((N₁ - P) * (N₂ - P)) / (P * (N₃ + P)) + one(T),
                #         reverse(2ρₘᵢₙ+2:2:2ρₘₐₓ),
                #         init=one(T)
                #     )

            end
        end
    end

end # @testmodule Boyle2016


@testitem "WignerDElement" setup=[Boyle2016] begin
    using Quaternionic
    using Random
    Random.seed!(1234)
    const T = Float64
    const ℓₘₐₓ = 8
    ϵₐ = 8eps(T)
    ϵᵣ = 20eps(T)

    Rs = [
        exp(3eps(T)*imx);
        exp(3eps(T)*imy);
        exp(3eps(T)*imz);
        Rotor{T}(imx)*exp(3eps(T)*imx);
        Rotor{T}(imx)*exp(3eps(T)*imy);
        Rotor{T}(imx)*exp(3eps(T)*imz);
        Rotor{T}(imy)*exp(3eps(T)*imx);
        Rotor{T}(imy)*exp(3eps(T)*imy);
        Rotor{T}(imy)*exp(3eps(T)*imz);
        Rotor{T}(imz)*exp(3eps(T)*imx);
        Rotor{T}(imz)*exp(3eps(T)*imy);
        Rotor{T}(imz)*exp(3eps(T)*imz);
        randn(Rotor{T}, 20)
    ]

    for R′ ∈ Rs
        for R ∈ (R′, conj(R′))
            𝔇1 = [
                Boyle2016.WignerDElement(R, ℓ, m′, m)
                for (ℓ, m′, m) ∈ eachrow(SphericalFunctions.WignerDrange(ℓₘₐₓ))
            ]
            𝔇2 = D_matrices(R, ℓₘₐₓ)
            @test 𝔇1 ≈ 𝔇2
        end
    end

end  # @testitem "WignerDElement"
