

struct WignerComputer{IT, RT, NT}
    ℓₘₐₓ::IT
    H⁻::Matrix{RT}
    H⁺::Matrix{RT}
    Wˡ::Matrix{NT}
    function WignerComputer(ℓₘₐₓ::IT, rt::Type{RT}, ::Type{NT}=rt) where {IT, RT, NT}
        if real(NT) ≠ RT
            error("RT=$RT is supposed to be the real type of NT=$NT.")
        end
        if IT <: Rational
            if denominator(ℓₘₐₓ) ≠ 2
                error("For IT=$IT <: Rational, ℓₘₐₓ must have denominator 2.")
            end
        end
        if ℓₘₐₓ < ℓₘᵢₙ(ℓₘₐₓ)
            error("ℓₘₐₓ=$ℓₘₐₓ must be non-negative.")
        end
        # `H⁻p` may (eventually) be required to store all the coefficients for Hˡ₀ₘ with
        # non-negative `m`; even though that won't strictly be necessary for `ℓₘₐₓ`, it is
        # just one extra `RT`, and may simplify the coding significantly.  `H⁺p` will
        # (eventually) be required to store all the coefficients for Hˡ⁺¹₀ₘ with
        # non-negative `m` (and that will be strictly necessary), so we give it one extra
        # column.
        H⁻p = Matrix{RT}(undef, 1, Int(ℓₘₐₓ-ℓₘᵢₙ(ℓₘₐₓ))+1)
        H⁺p = Matrix{RT}(undef, 1, Int(ℓₘₐₓ-ℓₘᵢₙ(ℓₘₐₓ))+2)
        Wˡp = Matrix{NT}(undef, Int(2ℓₘₐₓ)+1, Int(2ℓₘₐₓ)+1)
        new{IT, RT, NT}(ℓₘₐₓ, H⁻p, H⁺p, Wˡp)
    end
end
function (wc::WignerComputer{IT, RT, NT})(ℓ::IT) where {IT, RT, NT}
    if ℓ < ℓₘᵢₙ(wc) || ℓ > ℓₘₐₓ(wc)
        error(
            "ℓ=$ℓ is out of range for this WignerComputer "
            * "(which has ℓₘᵢₙ=$(ℓₘᵢₙ(wc)) and ℓₘₐₓ=$(ℓₘₐₓ(wc))).")
    end
    Wˡ = Redesign.WignerMatrix{IT, NT}(wc.Wˡparent, ℓ, ℓ)
    H⁻ = Redesign.Hˡrow{IT, RT}(wc.H⁻parent, ℓ, 0)
    H⁺ = Redesign.Hˡrow{IT, RT}(wc.H⁺parent, ℓ+1, 0)
    return Wˡ, H⁻, H⁺
end

function initialize!(computer::WignerComputer{IT, RT}, sinβ::RT, cosβ::RT) where {IT, RT}
    H⁰ = Redesign.Hˡrow(computer.H⁻parent, 0, 0)
    H¹ = Redesign.Hˡrow(computer.H⁺parent, 1, 0)
    Redesign.initialize!(H⁰, sinβ, cosβ)
    Redesign.initialize!(H¹, sinβ, cosβ)
    Redesign.recurrence_0_m!(H¹, H⁰, sinβ, cosβ)
    computer
end

"""
    incrementH!(computer, ℓ, sinβ, cosβ)


"""
function incrementH!(computer::DˡComputer{IT, NT}, ℓ::IT, sinβ::NT, cosβ::NT) where {IT, NT}
    Hˡ = Redesign.Hˡrow(computer.H⁻parent, ℓ, 0)
    Hˡ⁺¹ = Redesign.Hˡrow(computer.H⁺parent, ℓ+1, 0)
    Hˡ[0:0, 0:ℓ] .= Hˡ⁺¹[0:0, 0:ℓ]
    Redesign.recurrence_0_m!(Hˡ, Hˡ⁻¹, sinβ, cosβ)
    computer
end

function recurrence!(computer::DˡComputer{IT, NT}, α::NT, β::NT, γ::NT) where {IT, NT}
    H⁺, H⁻, Dˡ = computer.H⁺, computer.H⁻, computer.Dˡ
    sinβ, cosβ = sincos(β)
    eⁱᵅ, eⁱᵞ = cis(α), cis(γ)
    parent(Dˡ) .= 0
    Redesign.initialize!(H⁻, sinβ, cosβ)
    Redesign.initialize!(H⁺, sinβ, cosβ)
    Redesign.recurrence_0_m!(H⁻, H⁺, sinβ, cosβ)


    Redesign.recurrence_0_m!(Hˡ⁺¹, Dˡ, sinβ, cosβ)
    Redesign.recurrence_1_m!(Dˡ, Hˡ⁺¹, sinβ, cosβ)
    Redesign.recurrence_m′₊!(Dˡ, sinβ, cosβ)
    Redesign.recurrence_m′₋!(Dˡ, sinβ, cosβ)
    Redesign.impose_symmetries!(Dˡ)
    Redesign.convert_H_to_D!(Dˡ, eⁱᵅ, eⁱᵞ)
    Dˡ
end
