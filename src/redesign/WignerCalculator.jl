
struct WignerCalculator{IT, RT<:Real, NT<:Union{RT, Complex{RT}}}
    ℓₘₐₓ::IT
    m′ₘₐₓ::IT
    m′ₘᵢₙ::IT
    mₘₐₓ::IT
    mₘᵢₙ::IT
    Hᵃ::Matrix{RT}
    Hᵇ::Matrix{RT}
    Wˡ::Matrix{NT}
    swapH::Base.RefValue{Bool}  # w(ℓ) returns (Wˡ, Hᵇ, Hᵃ) if `true`, otherwise (Wˡ, Hᵃ, Hᵇ)
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

ℓₘᵢₙ(w::WignerCalculator{IT}) where {IT} = ℓₘᵢₙ(w.ℓₘₐₓ)
ℓₘₐₓ(w::WignerCalculator{IT}) where {IT} = w.ℓₘₐₓ
m′ₘₐₓ(w::WignerCalculator{IT}) where {IT} = w.m′ₘₐₓ
m′ₘᵢₙ(w::WignerCalculator{IT}) where {IT} = w.m′ₘᵢₙ
mₘₐₓ(w::WignerCalculator{IT}) where {IT} = w.mₘₐₓ
mₘᵢₙ(w::WignerCalculator{IT}) where {IT} = w.mₘᵢₙ

Hᵃ(w::WignerCalculator) = swapH(w) ? w.Hᵇ : w.Hᵃ
Hᵇ(w::WignerCalculator) = swapH(w) ? w.Hᵃ : w.Hᵇ
Wˡ(w::WignerCalculator) = w.Wˡ

function swapH(w::WignerCalculator)
    w.swapH[]
end

function swapH!(w::WignerCalculator)
    w.swapH[] = !w.swapH[]
    w
end

function fillW!(w::WignerCalculator{IT}, ℓ::IT) where {IT}
    Wˡ, Hˡ, Hˡ⁺¹ = w(ℓ)
    Wˡ[0:0, 0:ℓ] .= Hˡ[0:0, 0:ℓ]
    w
end

function (w::WignerCalculator{IT, RT, NT})(ℓ::IT) where {IT, RT, NT}
    if IT <: Rational
        if denominator(ℓ) ≠ 2
            error("For IT=$IT <: Rational, ℓ=$ℓ must have denominator 2")
        end
    end
    if ℓ < ℓₘᵢₙ(w) || ℓ > ℓₘₐₓ(w)
        error(
            "ℓ=$ℓ is out of range for this WignerCalculator " *
            "(which has ℓₘᵢₙ=$(ℓₘᵢₙ(w)) and ℓₘₐₓ=$(ℓₘₐₓ(w)))."
        )
    end
    let ℓₘᵢₙ=ℓₘᵢₙ(w), m′ₘₐₓ=min(ℓ, m′ₘₐₓ(w)), m′ₘᵢₙ=max(-ℓ, m′ₘᵢₙ(w)),
        mₘₐₓ=min(ℓ, mₘₐₓ(w)), mₘᵢₙ=max(-ℓ, mₘᵢₙ(w)),
        Wˡ=WignerMatrix(Wˡ(w), ℓ; m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ),
        Hˡ=WignerMatrix(Hᵃ(w), ℓ; m′ₘₐₓ=ℓₘᵢₙ, m′ₘᵢₙ=ℓₘᵢₙ, mₘₐₓ, mₘᵢₙ=ℓₘᵢₙ),
        Hˡ⁺¹=WignerMatrix(Hᵇ(w), ℓ+1; m′ₘₐₓ=ℓₘᵢₙ, m′ₘᵢₙ=ℓₘᵢₙ, mₘₐₓ=mₘₐₓ+1, mₘᵢₙ=ℓₘᵢₙ)

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
    w::WignerCalculator{IT, RT, NT}, α::RT, β::RT, γ::RT, ℓ::IT,
    skip_ℓ_recurrence::Bool=false
) where {IT<:Signed, RT, NT}
    eⁱᵅ, eⁱᵝ, eⁱᵞ = cis(α), cis(β), cis(γ)

    # NOTE: In the comments explaining the recurrence steps below, we use notation with
    # ℓₘᵢₙ=0 for simplicity, but the code should work for ℓₘᵢₙ=1//2 as well.

    if ℓ == ℓₘᵢₙ(w)
        # H⁰₀₀ = 1
        recurrence_step1!(w)

        # Record the result in Wˡ
        fillW!(w, ℓ)

        # Now, to leave `w` in a good state for the next ℓ, compute H¹₀ₘ and swap
        # H⁰₀ₘ -> H¹₀ₘ
        recurrence_step2!(w, eⁱᵝ, ℓ+1)
        swapH!(w)

        Wˡ, Hˡ, Hˡ⁺¹ = w(ℓ)
        Wˡ
    else
        if !skip_ℓ_recurrence
            # H⁰₀₀ = 1
            recurrence_step1!(w)

            for ℓ′ in ℓₘᵢₙ(w)+1:ℓ
                # Hˡ⁻¹₀ₘ -> Hˡ₀ₘ
                recurrence_step2!(w, eⁱᵝ, ℓ′)
                swapH!(w)
            end
        end

        let ℓ′ = ℓ+1
            # Hˡ₀ₘ -> Hˡ⁺¹₀ₘ
            recurrence_step2!(w, eⁱᵝ, ℓ′)
        end

        let
            Wˡ, Hˡ, Hˡ⁺¹ = w(ℓ)
            @info "recurrence! called with skip_ℓ_recurrence=$skip_ℓ_recurrence" Hˡ Hˡ⁺¹
        end

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

        # Swap the H matrices once more so that the current Hˡ⁺¹ is the next loop's Hˡ
        swapH!(w)

        Wˡ
    end

end
