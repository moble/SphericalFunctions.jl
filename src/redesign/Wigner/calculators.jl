struct WignerCalculator{IT, RT<:Real, NT<:Union{RT,Complex{RT}}, ST}
    H::WignerHCalculator{IT, RT, ST}
    e⁻ⁱᵐ′ᵅ::FixedSizeArrayDefault{NT}
    e⁻ⁱᵐᵞ::FixedSizeArrayDefault{NT}

    function WignerCalculator(
        ::Type{NT},
        rotors::AbstractVector{AbstractQuaternion{RT}},
        ℓₘₐₓ::IT, m′ₘₐₓ::IT=ℓₘₐₓ, m′ₘᵢₙ::IT=-ℓₘₐₓ
    ) where {IT, RT<:Real, NT<:Union{RT,Complex{RT}}}
        eⁱᵝ = FixedSizeVectorDefault{NT}(undef, length(rotors))
        e⁻ⁱᵐ′ᵅ = FixedSizeArrayDefault{NT}(undef, length(rotors), Int(ℓₘₐₓ-ℓₘᵢₙ(ℓₘₐₓ))+1)
        e⁻ⁱᵐᵞ  = FixedSizeArrayDefault{NT}(undef, length(rotors), Int(ℓₘₐₓ-ℓₘᵢₙ(ℓₘₐₓ))+1)
        @inbounds for (i, j) ∈ eachindex(rotors, eⁱᵝ)
            eⁱᵅ⁰, eⁱᵝ⁰, eⁱᵞ⁰ = Quaternionic.to_euler_phases(rotors[i])
            eⁱᵝ[j] = eⁱᵝ⁰
            conjexpiℓₘᵢₙα╱2, conjexpiℓₘᵢₙγ╱2 = if IT <: Rational
                conj(√eⁱᵅ⁰), conj(√eⁱᵞ⁰)
            else
                1, 1
            end
            α = ComplexPowers(eⁱᵅ⁰)
            γ = ComplexPowers(eⁱᵞ⁰)
            for (k, eⁱᵏᵅ, eⁱᵏᵞ) ∈ zip(0:Int(ℓₘₐₓ-ℓₘᵢₙ(ℓₘₐₓ)), α, γ)
                e⁻ⁱᵐ′ᵅ[j, k+1] = conj(eⁱᵏᵅ) * conjexpiℓₘᵢₙα╱2
                e⁻ⁱᵐᵞ[ j, k+1] = conj(eⁱᵏᵞ) * conjexpiℓₘᵢₙγ╱2
            end
        end
        H = WignerHCalculator(eⁱᵝ, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ)
        ST = typeof(parent(Hˡ(H)))
        new{IT, RT, NT, ST}(H, e⁻ⁱᵐ′ᵅ, e⁻ⁱᵐᵞ)
    end
end
