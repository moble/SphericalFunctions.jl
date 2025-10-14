function validate_index_ranges(ℓₘₐₓ::IT, m′ₘᵢₙ::IT, m′ₘₐₓ::IT, mₘᵢₙ::IT, mₘₐₓ::IT) where
    {IT<:Union{Signed, Rational}}
    if IT <: Rational
        if (
            denominator(ℓₘₐₓ) ≠ 2 ||
            denominator(m′ₘᵢₙ) ≠ 2 || denominator(m′ₘₐₓ) ≠ 2 ||
            denominator(mₘᵢₙ) ≠ 2 || denominator(mₘₐₓ) ≠ 2
        )
            error(
                "For IT=$IT <: Rational, indices must have denominator 2:\n"
                * "\tℓₘₐₓ=$ℓₘₐₓ, m′ₘᵢₙ=$m′ₘᵢₙ, m′ₘₐₓ=$m′ₘₐₓ, mₘᵢₙ=$mₘᵢₙ, mₘₐₓ=$mₘₐₓ."
            )
        end
    end

    if ℓₘₐₓ < ℓₘᵢₙ(ℓₘₐₓ)
        error("ℓₘₐₓ=$ℓₘₐₓ must be non-negative.")
    end

    # The m′ and m values must bracket ℓₘᵢₙ
    if m′ₘᵢₙ > ℓₘᵢₙ(ℓₘₐₓ)
        error("m′ₘᵢₙ=$m′ₘᵢₙ is too large for this index type.")
    end
    if m′ₘₐₓ < ℓₘᵢₙ(ℓₘₐₓ)
        error("m′ₘₐₓ=$m′ₘₐₓ is too small for this index type.")
    end
    if mₘᵢₙ > ℓₘᵢₙ(ℓₘₐₓ)
        error("mₘᵢₙ=$mₘᵢₙ is too large for this index type.")
    end
    if mₘₐₓ < ℓₘᵢₙ(ℓₘₐₓ)
        error("mₘₐₓ=$mₘₐₓ is too small for this index type.")
    end

    # The m′ and m values must be in range for ℓₘₐₓ
    if abs(m′ₘᵢₙ) > ℓₘₐₓ
        error("|m′ₘᵢₙ|=|$m′ₘᵢₙ| is too large for ℓₘₐₓ=$ℓₘₐₓ.")
    end
    if abs(m′ₘₐₓ) > ℓₘₐₓ
        error("|m′ₘₐₓ|=|$m′ₘₐₓ| is too large for ℓₘₐₓ=$ℓₘₐₓ.")
    end
    if abs(mₘᵢₙ) > ℓₘₐₓ
        error("|mₘᵢₙ|=|$mₘᵢₙ| is too large for ℓₘₐₓ=$ℓₘₐₓ.")
    end
    if abs(mₘₐₓ) > ℓₘₐₓ
        error("|mₘₐₓ|=|$mₘₐₓ| is too large for ℓₘₐₓ=$ℓₘₐₓ.")
    end

end

struct WignerComputer{IT, RT<:Real, NT<:Union{RT, Complex{RT}}}
    ℓₘₐₓ::IT
    m′ₘᵢₙ::IT
    m′ₘₐₓ::IT
    mₘᵢₙ::IT
    mₘₐₓ::IT
    H⁻::Matrix{RT}
    H⁺::Matrix{RT}
    Wˡ::Matrix{NT}
    function WignerComputer(
        ℓₘₐₓ::IT, rt::Type{RT}, ::Type{NT}=rt;
        m′ₘᵢₙ::IT=-ℓₘₐₓ, m′ₘₐₓ::IT=ℓₘₐₓ, mₘᵢₙ::IT=-ℓₘₐₓ, mₘₐₓ::IT=ℓₘₐₓ
    ) where {IT, RT, NT}
        if real(NT) ≠ RT
            error("RT=$RT is supposed to be the real type of NT=$NT.")
        end
        validate_index_ranges(ℓₘₐₓ, m′ₘᵢₙ, m′ₘₐₓ, mₘᵢₙ, mₘₐₓ)
        # `H⁻p` may (eventually) be required to store all the coefficients for Hˡ₀ₘ with
        # non-negative `m`; even though that won't strictly be necessary for `ℓₘₐₓ`, it is
        # just one extra `RT`, and may simplify the coding significantly.  `H⁺p` will
        # (eventually) be required to store all the coefficients for Hˡ⁺¹₀ₘ with
        # non-negative `m` (and that will be strictly necessary), so we give it one extra
        # column.
        H⁻p = Matrix{RT}(undef, 1, Int(ℓₘₐₓ-ℓₘᵢₙ(ℓₘₐₓ))+1)
        H⁺p = Matrix{RT}(undef, 1, Int(ℓₘₐₓ-ℓₘᵢₙ(ℓₘₐₓ))+2)
        Wˡp = Matrix{NT}(undef, Int(2ℓₘₐₓ)+1, Int(2ℓₘₐₓ)+1)
        new{IT, RT, NT}(ℓₘₐₓ, m′ₘᵢₙ, m′ₘₐₓ, mₘᵢₙ, mₘₐₓ, H⁻p, H⁺p, Wˡp)
    end
end

function (wc::WignerComputer{IT, RT, NT})(ℓ::IT) where {IT, RT, NT}
    if IT <: Rational
        if denominator(ℓ) ≠ 2
            error("For IT=$IT <: Rational, ℓ=$ℓ must have denominator 2")
        end
    end
    if ℓ < ℓₘᵢₙ(wc) || ℓ > ℓₘₐₓ(wc)
        error(
            "ℓ=$ℓ is out of range for this WignerComputer "
            * "(which has ℓₘᵢₙ=$(ℓₘᵢₙ(wc)) and ℓₘₐₓ=$(ℓₘₐₓ(wc))).")
    end
    Wˡ = WignerMatrix{IT, NT}(wc.Wˡ, ℓ, ℓ)
    H⁻ = Hˡrow{IT, RT}(wc.H⁻, ℓ, 0)
    H⁺ = Hˡrow{IT, RT}(wc.H⁺, ℓ+1, 0)
    return Wˡ, H⁻, H⁺
end

function WignerDComputer(
    ℓₘₐₓ::IT, ::Type{RT};
    m′ₘᵢₙ::IT=-ℓₘₐₓ, m′ₘₐₓ::IT=ℓₘₐₓ, mₘᵢₙ::IT=-ℓₘₐₓ, mₘₐₓ::IT=ℓₘₐₓ
) where {IT, RT<:Real}
    NT = complex(RT)
    WignerComputer(ℓₘₐₓ, RT, NT; m′ₘᵢₙ, m′ₘₐₓ, mₘᵢₙ, mₘₐₓ)
end

function WignerdComputer(
    ℓₘₐₓ::IT, ::Type{RT};
    m′ₘᵢₙ::IT=-ℓₘₐₓ, m′ₘₐₓ::IT=ℓₘₐₓ, mₘᵢₙ::IT=-ℓₘₐₓ, mₘₐₓ::IT=ℓₘₐₓ
) where {IT, RT<:Real}
    WignerComputer(ℓₘₐₓ, RT, RT; m′ₘᵢₙ, m′ₘₐₓ, mₘᵢₙ, mₘₐₓ)
end

function recurrence_step1!(w::WignerComputer{IT}) where {IT<:Signed}
    W⁰, H⁰, H¹ = w(0)
    initialize!(H⁰)
    w
end

function recurrence_step2!(w::WignerComputer{IT}, eⁱᵝ, ℓ) where {IT<:Signed}
    Wˡ⁻¹, Hˡ⁻¹, Hˡ = w(ℓ-1)
    sinβ, cosβ = reim(eⁱᵝ)
    recurrence_step2!(Hˡ, Hˡ⁻¹, sinβ, cosβ)
    Wˡ, Hˡ, Hˡ⁺¹ = w(ℓ)
    Wˡ[0:0, 0:ℓ] .= Hˡ[0:0, 0:ℓ]
    recurrence_step2!(Hˡ⁺¹, Hˡ, sinβ, cosβ)
    w
end

function recurrence_step3!(w::WignerComputer{IT}, eⁱᵝ, ℓ) where {IT<:Signed}
    Wˡ, Hˡ, Hˡ⁺¹ = w(ℓ)
    sinβ, cosβ = reim(eⁱᵝ)
    recurrence_step3!(Wˡ, Hˡ⁺¹, sinβ, cosβ)
    w
end

function recurrence_step4!(w::WignerComputer{IT}, eⁱᵝ, ℓ) where {IT<:Signed}
    Wˡ, Hˡ, Hˡ⁺¹ = w(ℓ)
    sinβ, cosβ = reim(eⁱᵝ)
    recurrence_step4!(Wˡ, sinβ, cosβ)
    w
end

function recurrence_step5!(w::WignerComputer{IT}, eⁱᵝ, ℓ) where {IT<:Signed}
    Wˡ, Hˡ, Hˡ⁺¹ = w(ℓ)
    sinβ, cosβ = reim(eⁱᵝ)
    recurrence_step5!(Wˡ, sinβ, cosβ)
    w
end

function recurrence_step6!(w::WignerComputer{IT}, ℓ) where {IT<:Signed}
    Wˡ, Hˡ, Hˡ⁺¹ = w(ℓ)
    recurrence_step6!(Wˡ)
    w
end

function recurrence!(w::WignerComputer{IT, RT}, α::RT, β::RT, γ::RT, ℓ::IT) where {IT<:Signed, RT}
    eⁱᵅ, eⁱᵝ, eⁱᵞ = cis(α), cis(β), cis(γ)
    recurrence_step1!(w)
    for ℓ′ in 1:ℓ
        recurrence_step2!(w, eⁱᵝ, ℓ′)
        recurrence_step3!(w, eⁱᵝ, ℓ′)
        recurrence_step4!(w, eⁱᵝ, ℓ′)
        recurrence_step5!(w, eⁱᵝ, ℓ′)
        recurrence_step6!(w, ℓ′)
    end
    Wˡ, Hˡ, Hˡ⁺¹ = w(ℓ)
    convert_H_to_D!(Wˡ, eⁱᵅ, eⁱᵞ)
    Wˡ
end
