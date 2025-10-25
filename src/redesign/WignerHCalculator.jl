struct WignerHCalculator{IT, RT<:Real, ST}
    ℓₘₐₓ::IT
    m′ₘₐₓ::IT
    m′ₘᵢₙ::IT
    Hᵃ::HAxis{IT, RT}
    Hᵇ::HAxis{IT, RT}
    Wˡ::HWedge{IT, RT, ST}
    swapH::Base.RefValue{Bool}  # w(ℓ) returns (Wˡ, Hᵇ, Hᵃ) if `true`, otherwise (Wˡ, Hᵃ, Hᵇ)
    # function WignerHCalculator(
    #     ℓₘₐₓ::IT, rt::Type{RT};
    #     m′ₘₐₓ::IT=ℓₘₐₓ, m′ₘᵢₙ::IT=-ℓₘₐₓ
    # ) where {IT, RT}
    #     if real(NT) ≠ RT
    #         error("RT=$RT is supposed to be the real type of NT=$NT.")
    #     end
    #     validate_index_ranges(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    #     # One of the H matrices will (eventually) be required to store all the coefficients
    #     # for Hˡ⁺¹₀ₘ with non-negative `m` (and that will be strictly necessary), so we give
    #     # it one extra column.  Since we may not know which one that will be, we give both
    #     # of them that extra column.
    #     Hᵃp = Matrix{RT}(undef, 1, Int(mₘₐₓ-ℓₘᵢₙ(ℓₘₐₓ))+2)
    #     Hᵇp = Matrix{RT}(undef, 1, Int(mₘₐₓ-ℓₘᵢₙ(ℓₘₐₓ))+2)
    #     Wˡp = Matrix{NT}(undef, Int(m′ₘₐₓ-m′ₘᵢₙ)+1, Int(mₘₐₓ-mₘᵢₙ)+1)
    #     new{IT, RT, NT}(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, Hᵃp, Hᵇp, Wˡp, Ref(false))
    # end
end

ℓₘᵢₙ(w::WignerHCalculator{IT}) where {IT} = ℓₘᵢₙ(w.ℓₘₐₓ)
ℓₘₐₓ(w::WignerHCalculator{IT}) where {IT} = w.ℓₘₐₓ
m′ₘₐₓ(w::WignerHCalculator{IT}) where {IT} = w.m′ₘₐₓ
m′ₘᵢₙ(w::WignerHCalculator{IT}) where {IT} = w.m′ₘᵢₙ

Hᵃ(w::WignerHCalculator) = swapH(w) ? w.Hᵇ : w.Hᵃ
Hᵇ(w::WignerHCalculator) = swapH(w) ? w.Hᵃ : w.Hᵇ
Wˡ(w::WignerHCalculator) = w.Wˡ

function Base.fill!(w::WignerHCalculator{IT, RT}, v::Real) where {IT, RT}
    let Wˡ = Wˡ(w), Hᵃ = Hᵃ(w), Hᵇ = Hᵇ(w)
        fill!(Wˡ, eltype(Wˡ)(v))
        fill!(Hᵃ, eltype(Hᵃ)(v))
        fill!(Hᵇ, eltype(Hᵇ)(v))
    end
    w
end

function swapH(w::WignerHCalculator)
    w.swapH[]
end

function swapH!(w::WignerHCalculator)
    w.swapH[] = !w.swapH[]
    w
end

function fillW!(w::WignerHCalculator{IT}, ℓ::IT) where {IT}
    Wˡ, Hˡ, Hˡ⁺¹ = w(ℓ)
    @views copyto!(Wˡ[0:0, 0:ℓ], Hˡ[0:0, 0:ℓ])
    w
end

function (w::WignerHCalculator{IT, RT, NT})(ℓ::IT) where {IT, RT, NT}
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

# function WignerDCalculator(
#     ℓₘₐₓ::IT, ::Type{RT};
#     mp_max::IT=ℓₘₐₓ, mp_min::IT=-ℓₘₐₓ, m_max::IT=ℓₘₐₓ, m_min::IT=-ℓₘₐₓ,
#     m′ₘₐₓ::IT=mp_max, m′ₘᵢₙ::IT=mp_min, mₘₐₓ::IT=m_max, mₘᵢₙ::IT=m_min
# ) where {IT, RT<:Real}
#     NT = complex(RT)
#     WignerHCalculator(ℓₘₐₓ, RT, NT; m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
# end

# function WignerdCalculator(
#     ℓₘₐₓ::IT, ::Type{RT};
#     mp_max::IT=ℓₘₐₓ, mp_min::IT=-ℓₘₐₓ, m_max::IT=ℓₘₐₓ, m_min::IT=-ℓₘₐₓ,
#     m′ₘₐₓ::IT=mp_max, m′ₘᵢₙ::IT=mp_min, mₘₐₓ::IT=m_max, mₘᵢₙ::IT=m_min
# ) where {IT, RT<:Real}
#     WignerHCalculator(ℓₘₐₓ, RT, RT; m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
# end

function recurrence_step1!(w::WignerHCalculator{IT}) where {IT<:Signed}
    ℓ = ℓₘᵢₙ(w)
    W⁰, H⁰, H¹ = w(ℓ)
    recurrence_step1!(H⁰)
    fillW!(w, ℓ)
    w
end

function recurrence_step2!(w::WignerHCalculator{IT}, eⁱᵝ, ℓ) where {IT<:Signed}
    Wˡ⁻¹, Hˡ⁻¹, Hˡ = w(ℓ-1)
    cosβ, sinβ = reim(eⁱᵝ)
    recurrence_step2!(Hˡ, Hˡ⁻¹, sinβ, cosβ)
    w
end

function recurrence_step3!(w::WignerHCalculator{IT}, eⁱᵝ, ℓ) where {IT<:Signed}
    Wˡ, Hˡ, Hˡ⁺¹ = w(ℓ)
    cosβ, sinβ = reim(eⁱᵝ)
    recurrence_step3!(Wˡ, Hˡ⁺¹, sinβ, cosβ)
    w
end

function recurrence_step4!(w::WignerHCalculator{IT}, eⁱᵝ, ℓ) where {IT<:Signed}
    Wˡ, Hˡ, Hˡ⁺¹ = w(ℓ)
    cosβ, sinβ = reim(eⁱᵝ)
    recurrence_step4!(Wˡ, sinβ, cosβ)
    w
end

function recurrence_step5!(w::WignerHCalculator{IT}, eⁱᵝ, ℓ) where {IT<:Signed}
    Wˡ, Hˡ, Hˡ⁺¹ = w(ℓ)
    cosβ, sinβ = reim(eⁱᵝ)
    recurrence_step5!(Wˡ, sinβ, cosβ)
    w
end

function recurrence_step6!(w::WignerHCalculator{IT}, ℓ) where {IT<:Signed}
    Wˡ, Hˡ, Hˡ⁺¹ = w(ℓ)
    recurrence_step6!(Wˡ)
    w
end
    
function recurrence!(
    w::WignerHCalculator{IT, RT, NT}, α::RT, β::RT, γ::RT, ℓ::IT,
    skip_ℓ_recurrence::Bool=false
) where {IT<:Signed, RT, NT<:Complex}
    eⁱᵅ, eⁱᵝ, eⁱᵞ = cis(α), cis(β), cis(γ)
    recurrence!(w, eⁱᵅ, eⁱᵝ, eⁱᵞ, ℓ, skip_ℓ_recurrence)
end

function recurrence!(
    w::WignerHCalculator{IT, RT, NT}, β::RT, ℓ::IT,
    skip_ℓ_recurrence::Bool=false
) where {IT<:Signed, RT, NT<:Real}
    eⁱᵝ = cis(β)
    recurrence!(w, eⁱᵝ, ℓ, skip_ℓ_recurrence)
end

function recurrence!(
    w::WignerHCalculator{IT, RT, NT}, eⁱᵅ::Complex{RT}, eⁱᵝ::Complex{RT}, eⁱᵞ::Complex{RT},
    ℓ::IT, skip_ℓ_recurrence::Bool=false
) where {IT<:Signed, RT, NT<:Complex}
    _recurrence!(w, eⁱᵝ, ℓ, skip_ℓ_recurrence)
    Wˡ, Hˡ, Hˡ⁺¹ = w(ℓ)
    convert_H_to_D!(Wˡ, eⁱᵅ, eⁱᵞ)
    Wˡ
end

function recurrence!(
    w::WignerHCalculator{IT, RT, NT}, eⁱᵝ::Complex{RT}, ℓ::IT,
    skip_ℓ_recurrence::Bool=false
) where {IT<:Signed, RT, NT<:Real}
    _recurrence!(w, eⁱᵝ, ℓ, skip_ℓ_recurrence)
    Wˡ, Hˡ, Hˡ⁺¹ = w(ℓ)
    convert_H_to_d!(Wˡ)
    Wˡ
end

function _recurrence!(
    w::WignerHCalculator{IT, RT, NT}, eⁱᵝ::Complex{RT}, ℓ::IT,
    skip_ℓ_recurrence::Bool=false
) where {IT<:Signed, RT, NT}
    # NOTE: In the comments explaining the recurrence steps below, we use notation with
    # ℓₘᵢₙ=0 for simplicity, but this sequence may work for ℓₘᵢₙ=1//2 as well.

    if ℓ == ℓₘᵢₙ(w)
        recurrence_step1!(w)  # H⁰₀₀ = 1
        fillW!(w, ℓ)  # Record the result in Wˡ

        # Now, to leave `w` in a good state for the next ℓ, compute H¹₀ₘ and swap.
        recurrence_step2!(w, eⁱᵝ, ℓ+1)  # H⁰₀ₘ -> H¹₀ₘ
        swapH!(w)
    else
        if !skip_ℓ_recurrence
            recurrence_step1!(w)  # H⁰₀₀ = 1

            for ℓ′ in ℓₘᵢₙ(w)+1:ℓ
                recurrence_step2!(w, eⁱᵝ, ℓ′)  # Hˡ⁻¹₀ₘ -> Hˡ₀ₘ
                swapH!(w)  # Prepare for the next iteration of ℓ′
            end
        end

        # Do one more step of the recurrence to get Hˡ⁺¹₀ₘ, regardless of whether or not we
        # asked to skip the ℓ recurrence.  If we did, the user is responsible for having
        # already set Hˡ₀ₘ correctly; if we didn't, we just set it in the loop above.
        let ℓ′ = ℓ+1
            recurrence_step2!(w, eⁱᵝ, ℓ′)  # Hˡ₀ₘ -> Hˡ⁺¹₀ₘ
        end

        fillW!(w, ℓ)  # Copy Hˡ₀ₘ to Wˡ₀ₘ
        recurrence_step3!(w, eⁱᵝ, ℓ)  # Hˡ⁺¹₀ₘ -> Hˡ₁ₘ
        recurrence_step4!(w, eⁱᵝ, ℓ)  # Hˡₘ′ₘ₋₁, Hˡₘ′₋₁ₘ, Hˡₘ′ₘ₊₁ -> Hˡₘ′₊₁ₘ
        recurrence_step5!(w, eⁱᵝ, ℓ)  # Hˡₘ′ₘ₋₁, Hˡₘ′₊₁ₘ, Hˡₘ′ₘ₊₁ -> Hˡₘ′₋₁ₘ

        # Impose symmetries
        recurrence_step6!(w, ℓ)

        # Swap the H matrices once more so that the current Hˡ⁺¹ is the next loop's Hˡ
        swapH!(w)
    end

    Wˡ, Hˡ, Hˡ⁺¹ = w(ℓ)
    Wˡ

end
