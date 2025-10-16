
struct WignerCalculator{IT, RT<:Real, NT<:Union{RT, Complex{RT}}}
    ℓₘₐₓ::IT
    m′ₘₐₓ::IT
    m′ₘᵢₙ::IT
    mₘₐₓ::IT
    mₘᵢₙ::IT
    Hᵃ::Matrix{RT}
    Hᵇ::Matrix{RT}
    Wˡ::Matrix{NT}
    swapH::Base.RefValue{Bool}  # wc(ℓ) returns (Wˡ, Hᵇ, Hᵃ) if `true`, otherwise (Wˡ, Hᵃ, Hᵇ)
    function WignerCalculator(
        ℓₘₐₓ::IT, rt::Type{RT}, ::Type{NT}=rt;
        m′ₘₐₓ::IT=ℓₘₐₓ, m′ₘᵢₙ::IT=-ℓₘₐₓ, mₘₐₓ::IT=ℓₘₐₓ, mₘᵢₙ::IT=-ℓₘₐₓ
    ) where {IT, RT, NT}
        if real(NT) ≠ RT
            error("RT=$RT is supposed to be the real type of NT=$NT.")
        end
        validate_index_ranges(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
        # One of the H matrices will (eventually) be required to store all the coefficients
        # for Hˡ⁺¹₀ₘ with non-negative `m` (and that will be strictly necessary), so we give
        # it one extra column.  Since we may not know which one that will be, we give both
        # of them that extra column.
        Hᵃp = Matrix{RT}(undef, 1, Int(mₘₐₓ-ℓₘᵢₙ(ℓₘₐₓ))+2)
        Hᵇp = Matrix{RT}(undef, 1, Int(mₘₐₓ-ℓₘᵢₙ(ℓₘₐₓ))+2)
        Wˡp = Matrix{NT}(undef, Int(m′ₘₐₓ-m′ₘᵢₙ)+1, Int(mₘₐₓ-mₘᵢₙ)+1)
        new{IT, RT, NT}(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, Hᵃp, Hᵇp, Wˡp, Ref(false))
    end
end

ℓₘᵢₙ(wc::WignerCalculator{IT}) where {IT} = ℓₘᵢₙ(wc.ℓₘₐₓ)
ℓₘₐₓ(wc::WignerCalculator{IT}) where {IT} = wc.ℓₘₐₓ
m′ₘₐₓ(wc::WignerCalculator{IT}) where {IT} = wc.m′ₘₐₓ
m′ₘᵢₙ(wc::WignerCalculator{IT}) where {IT} = wc.m′ₘᵢₙ
mₘₐₓ(wc::WignerCalculator{IT}) where {IT} = wc.mₘₐₓ
mₘᵢₙ(wc::WignerCalculator{IT}) where {IT} = wc.mₘᵢₙ

Hᵃ(wc::WignerCalculator) = swapH(wc) ? wc.Hᵇ : wc.Hᵃ
Hᵇ(wc::WignerCalculator) = swapH(wc) ? wc.Hᵃ : wc.Hᵇ
Wˡ(wc::WignerCalculator) = wc.Wˡ

function swapH(wc::WignerCalculator)
    wc.swapH[]
end

function swapH!(wc::WignerCalculator)
    wc.swapH[] = !wc.swapH[]
    wc
end

function fillW!(wc::WignerCalculator{IT}, ℓ::IT) where {IT}
    Wˡ, Hˡ, Hˡ⁺¹ = wc(ℓ)
    Wˡ[0:0, 0:ℓ] .= Hˡ[0:0, 0:ℓ]
    wc
end

function (wc::WignerCalculator{IT, RT, NT})(ℓ::IT) where {IT, RT, NT}
    if IT <: Rational
        if denominator(ℓ) ≠ 2
            error("For IT=$IT <: Rational, ℓ=$ℓ must have denominator 2")
        end
    end
    if ℓ < ℓₘᵢₙ(wc) || ℓ > ℓₘₐₓ(wc)
        error(
            "ℓ=$ℓ is out of range for this WignerCalculator " *
            "(which has ℓₘᵢₙ=$(ℓₘᵢₙ(wc)) and ℓₘₐₓ=$(ℓₘₐₓ(wc)))."
        )
    end
    let ℓₘᵢₙ=ℓₘᵢₙ(wc), m′ₘₐₓ=min(ℓ, m′ₘₐₓ(wc)), m′ₘᵢₙ=max(-ℓ, m′ₘᵢₙ(wc)),
        mₘₐₓ=min(ℓ, mₘₐₓ(wc)), mₘᵢₙ=max(-ℓ, mₘᵢₙ(wc)),
        Wˡ=WignerMatrix(Wˡ(wc), ℓ; m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ),
        Hˡ=WignerMatrix(Hᵃ(wc), ℓ; m′ₘₐₓ=ℓₘᵢₙ, m′ₘᵢₙ=ℓₘᵢₙ, mₘₐₓ, mₘᵢₙ=ℓₘᵢₙ),
        Hˡ⁺¹=WignerMatrix(Hᵇ(wc), ℓ+1; m′ₘₐₓ=ℓₘᵢₙ, m′ₘᵢₙ=ℓₘᵢₙ, mₘₐₓ=mₘₐₓ+1, mₘᵢₙ=ℓₘᵢₙ)

        Wˡ, Hˡ, Hˡ⁺¹
    end
end

function WignerDCalculator(
    ℓₘₐₓ::IT, ::Type{RT};
    mp_max::IT=ℓₘₐₓ, mp_min::IT=-ℓₘₐₓ, m_max::IT=ℓₘₐₓ, m_min::IT=-ℓₘₐₓ,
    m′ₘₐₓ::IT=mp_max, m′ₘᵢₙ::IT=mp_min, mₘₐₓ::IT=m_max, mₘᵢₙ::IT=m_min
) where {IT, RT<:Real}
    NT = complex(RT)
    WignerCalculator(ℓₘₐₓ, RT, NT; m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
end

function WignerdCalculator(
    ℓₘₐₓ::IT, ::Type{RT};
    mp_max::IT=ℓₘₐₓ, mp_min::IT=-ℓₘₐₓ, m_max::IT=ℓₘₐₓ, m_min::IT=-ℓₘₐₓ,
    m′ₘₐₓ::IT=mp_max, m′ₘᵢₙ::IT=mp_min, mₘₐₓ::IT=m_max, mₘᵢₙ::IT=m_min
) where {IT, RT<:Real}
    WignerCalculator(ℓₘₐₓ, RT, RT; m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
end

function recurrence_step1!(w::WignerCalculator{IT}) where {IT<:Signed}
    W⁰, H⁰, H¹ = w(0)
    recurrence_step1!(H⁰)
    W⁰[0, 0] = H⁰[0, 0]
    w
end

function recurrence_step2!(w::WignerCalculator{IT}, eⁱᵝ, ℓ) where {IT<:Signed}
    Wˡ⁻¹, Hˡ⁻¹, Hˡ = w(ℓ-1)
    cosβ, sinβ = reim(eⁱᵝ)
    recurrence_step2!(Hˡ, Hˡ⁻¹, sinβ, cosβ)
    w
end

function recurrence_step3!(w::WignerCalculator{IT}, eⁱᵝ, ℓ) where {IT<:Signed}
    Wˡ, Hˡ, Hˡ⁺¹ = w(ℓ)
    cosβ, sinβ = reim(eⁱᵝ)
    recurrence_step3!(Wˡ, Hˡ⁺¹, sinβ, cosβ)
    w
end

function recurrence_step4!(w::WignerCalculator{IT}, eⁱᵝ, ℓ) where {IT<:Signed}
    Wˡ, Hˡ, Hˡ⁺¹ = w(ℓ)
    cosβ, sinβ = reim(eⁱᵝ)
    recurrence_step4!(Wˡ, sinβ, cosβ)
    w
end

function recurrence_step5!(w::WignerCalculator{IT}, eⁱᵝ, ℓ) where {IT<:Signed}
    Wˡ, Hˡ, Hˡ⁺¹ = w(ℓ)
    cosβ, sinβ = reim(eⁱᵝ)
    recurrence_step5!(Wˡ, sinβ, cosβ)
    w
end

function recurrence_step6!(w::WignerCalculator{IT}, ℓ) where {IT<:Signed}
    Wˡ, Hˡ, Hˡ⁺¹ = w(ℓ)
    recurrence_step6!(Wˡ)
    w
end
    
function recurrence!(
    w::WignerCalculator{IT, RT, NT}, α::RT, β::RT, γ::RT, ℓ::IT
) where {IT<:Signed, RT, NT}
    eⁱᵅ, eⁱᵝ, eⁱᵞ = cis(α), cis(β), cis(γ)

    # H⁰₀₀ = 1
    recurrence_step1!(w)

    if ℓ == ℓₘᵢₙ(w)
        fillW!(w, ℓₘᵢₙ(w))
        Wˡ, Hˡ, Hˡ⁺¹ = w(ℓ)
        Wˡ
    else
        for ℓ′ in ℓₘᵢₙ(w)+1:ℓ+1
            # Hˡ⁻¹₀ₘ -> Hˡ₀ₘ
            recurrence_step2!(w, eⁱᵝ, ℓ′)
            swapH!(w)
        end
        swapH!(w)

        # Copy Hˡ₀ₘ to Wˡ₀ₘ
        fillW!(w, ℓ)

        # Hˡ⁺¹₀ₘ -> Hˡ₁ₘ
        recurrence_step3!(w, eⁱᵝ, ℓ)

        # Hˡₘ′ₘ₋₁, Hˡₘ′₋₁ₘ, Hˡₘ′ₘ₊₁ -> Hˡₘ′₊₁ₘ
        recurrence_step4!(w, eⁱᵝ, ℓ)

        # Hˡₘ′ₘ₋₁, Hˡₘ′₊₁ₘ, Hˡₘ′ₘ₊₁ -> Hˡₘ′₋₁ₘ
        recurrence_step5!(w, eⁱᵝ, ℓ)

        # Impose symmetries
        recurrence_step6!(w, ℓ)

        # Finish conversion to Dˡ
        Wˡ, Hˡ, Hˡ⁺¹ = w(ℓ)
        convert_H_to_D!(Wˡ, eⁱᵅ, eⁱᵞ)
        Wˡ
    end

end
