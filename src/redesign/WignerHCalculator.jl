struct WignerHCalculator{IT, RT<:Real, ST}
    ℓₘₐₓ::IT
    m′ₘₐₓ::IT
    m′ₘᵢₙ::IT
    h⃗ᵃ::HAxis{IT, RT}
    h⃗ᵇ::HAxis{IT, RT}
    Hˡ::HWedge{IT, RT, ST}
    eⁱᵝ::FixedSizeVectorDefault{Complex{RT}}
    swapH::Base.RefValue{Bool}  # h⃗ˡ(w) returns h⃗ᵃ if `false`, otherwise h⃗ᵇ; and vice versa for h⃗ˡ⁺¹(w)
    function WignerHCalculator(
        eⁱᵝ::AbstractVector{Complex{RT}}, ℓₘₐₓ::IT, m′ₘₐₓ::IT=ℓₘₐₓ, m′ₘᵢₙ::IT=-ℓₘₐₓ
    ) where {IT, RT<:Real}
        eⁱᵝ = FixedSizeVector(eⁱᵝ)
        Nᵣ = length(eⁱᵝ)
        # One of the H matrices will (eventually) be required to store all the coefficients
        # for Hˡ⁺¹₀ₘ with non-negative `m` (and that will be strictly necessary), so we give
        # it one extra column.  Since we may not know which one that will be, we give both
        # of them that extra column.
        h⃗ᵃ = HAxis(RT, Nᵣ, ℓₘₐₓ+1)
        h⃗ᵇ = HAxis(RT, Nᵣ, ℓₘₐₓ+1)
        Hˡ = HWedge(RT, Nᵣ, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ)
        new{IT, RT, NT}(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, h⃗ᵃ, h⃗ᵇ, Hˡ, eⁱᵝ, Ref(false))
    end
end

consistent_ℓ(w::WignerHCalculator{IT}) where {IT} = ℓ(h⃗ˡ(w)) == ℓ(Hˡ(w)) == ℓ(h⃗ˡ⁺¹(w)) - 1

function ℓ(w::WignerHCalculator{IT}) where {IT}
    if !consistent_ℓ(w)
        error(
            "Inconsistent ℓ values in WignerHCalculator:\n"
            * "    ℓ(h⃗ˡ)=$(ℓ(h⃗ˡ(w))), ℓ(h⃗ˡ⁺¹)=$(ℓ(h⃗ˡ⁺¹(w))), ℓ(Hˡ)=$(ℓ(Hˡ(w)))."
        )
    end
    ℓ(Hˡ(w))
end
ℓₘᵢₙ(w::WignerHCalculator{IT}) where {IT} = ℓₘᵢₙ(w.ℓₘₐₓ)
ℓₘₐₓ(w::WignerHCalculator{IT}) where {IT} = w.ℓₘₐₓ
m′ₘₐₓ(w::WignerHCalculator{IT}) where {IT} = w.m′ₘₐₓ
m′ₘᵢₙ(w::WignerHCalculator{IT}) where {IT} = w.m′ₘᵢₙ

h⃗ˡ(w::WignerHCalculator) = swapH(w) ? w.h⃗ᵇ : w.h⃗ᵃ
h⃗ˡ⁺¹(w::WignerHCalculator) = swapH(w) ? w.h⃗ᵃ : w.h⃗ᵇ
Hˡ(w::WignerHCalculator) = w.Hˡ
eⁱᵝ(w::WignerHCalculator) = w.eⁱᵝ

function Base.fill!(w::WignerHCalculator{IT, RT}, v::Real) where {IT, RT}
    let h⃗ˡ = h⃗ˡ(w), h⃗ˡ⁺¹ = h⃗ˡ⁺¹(w), Hˡ = Hˡ(w), v = convert(RT, v)
        fill!(parent(h⃗ˡ), v)
        fill!(parent(h⃗ˡ⁺¹), v)
        fill!(parent(Hˡ), v)
    end
    w
end

function swapH(w::WignerHCalculator)
    w.swapH[]
end

function increment_axes!(w::WignerHCalculator)
    # The data that is now stored as h⃗ˡ(w) will get swapped below so that it will be
    # returned by h⃗ˡ⁺¹(w), so we need to increment its ℓ value twice.
    h⃗ˡ = h⃗ˡ(w)
    h⃗ˡ.ℓ = ℓ(h⃗ˡ) + 2
    # The data that is now stored as h⃗ˡ⁺¹(w) will get swapped below so that it will be
    # returned by h⃗ˡ(w), which will already be correct for the next ℓ value.
    w.swapH[] = !w.swapH[]
    w
end

function increment_ℓ!(w::WignerHCalculator)
    increment_axes!(w)
    Hˡ = Hˡ(w)
    Hˡ.ℓ = ℓ(Hˡ) + 1
    w
end

function fillHˡ₀ₘ!(w::WignerHCalculator{IT}) where {IT}
    let h⃗ˡ = h⃗ˡ(w), Hˡ = Hˡ(w)
        if ℓ(h⃗ˡ) != ℓ(Hˡ)
            error("Cannot fill Hˡ₀ₘ for ℓ=$(ℓ(Hˡ)) from h⃗ˡ for ℓ=$(ℓ(h⃗ˡ)).")
        end
        # Get the index to the start of the central row, m′ = 0 or 1//2
        iˡ₀₀ = row_index(w)[Int(ℓₘᵢₙ(w) - m′ₘᵢₙ(w)) + 1]
        # Figure out how many entries to copy
        N = Nᵣ(Hˡ) * (Int(ℓ(w) - ℓₘᵢₙ(w)) + 1)
        # Now just copy that many entries from h⃗ˡ₀ₘ into row Hˡ₀ₘ
        copyto!(parent(Hˡ), iˡ₀₀, parent(h⃗ˡ), 1, N)
    end
    w
end

function Base.setproperty!(w::WignerHCalculator{IT}, s::Symbol, ℓ::IIT) where {IT, IIT}
    if s === :ℓ
        h⃗ˡ(w).ℓ = ℓ
        h⃗ˡ⁺¹(w).ℓ = ℓ+1
        Hˡ(w).ℓ = ℓ
        ℓ
    else
        error("Cannot set property `$s` on HWedge; only `ℓ` is allowed to be changed.")
    end
end

function recurrence_step1!(w::WignerHCalculator{IT}) where {IT<:Signed}
    let h⃗⁰ = h⃗ˡ(w), ℓ = ℓ(h⃗⁰)
        if ℓ ≠ ℓₘᵢₙ(w)
            error("recurrence_step1! can only be called for ℓ=$(ℓₘᵢₙ(w)); current ℓ=$ℓ.")
        end
        parent(h⃗⁰)[:, 1] .= 1
    end
    w
end

function recurrence_step2!(w::WignerHCalculator{IT}) where {IT<:Signed}
    let h⃗ˡ = h⃗ˡ(w), h⃗ˡ⁺¹ = h⃗ˡ⁺¹(w), ℓ = ℓ(h⃗ˡ), eⁱᵝ = eⁱᵝ(w)
        if ℓ < ℓₘᵢₙ(IT)
            error(
                "recurrence_step2! can only be called for ℓ≥ℓₘᵢₙ=$(ℓₘᵢₙ(IT)); current ℓ=$ℓ."
            )
        end
        recurrence_step2!(h⃗ˡ⁺¹, h⃗ˡ, eⁱᵝ)
    end
    w
end

function recurrence_step3!(w::WignerHCalculator{IT}) where {IT<:Signed}
    let Hˡ = Hˡ(w), h⃗ˡ⁺¹ = h⃗ˡ⁺¹(w), eⁱᵝ = eⁱᵝ(w)
        recurrence_step3!(Hˡ, h⃗ˡ⁺¹, eⁱᵝ)
    end
    w
end

function recurrence_step4!(w::WignerHCalculator{IT}) where {IT<:Signed}
    let Hˡ = Hˡ(w), eⁱᵝ = eⁱᵝ(w)
        recurrence_step4!(Hˡ, eⁱᵝ)
    end
    w
end

function recurrence_step5!(w::WignerHCalculator{IT}, eⁱᵝ) where {IT<:Signed}
    let Hˡ = Hˡ(w), eⁱᵝ = eⁱᵝ(w)
        recurrence_step5!(Hˡ, eⁱᵝ)
    end
    w
end

function recurrence_step6!(w::WignerHCalculator{IT}, ℓ) where {IT<:Signed}
    let Hˡ = Hˡ(w)
        recurrence_step6!(Hˡ)
    end
    w
end
    
function recurrence!(
    w::WignerHCalculator{IT, RT, NT}, α::RT, β::RT, γ::RT, ℓ::IT
) where {IT<:Signed, RT, NT<:Complex}
    eⁱᵅ, eⁱᵝ, eⁱᵞ = cis(α), cis(β), cis(γ)
    recurrence!(w, eⁱᵅ, eⁱᵝ, eⁱᵞ, ℓ)
end

function recurrence!(
    w::WignerHCalculator{IT, RT, NT}, β::RT, ℓ::IT
) where {IT<:Signed, RT, NT<:Real}
    eⁱᵝ = cis(β)
    recurrence!(w, eⁱᵝ, ℓ)
end

function recurrence!(
    w::WignerHCalculator{IT, RT, NT}, eⁱᵅ::Complex{RT}, eⁱᵝ::Complex{RT}, eⁱᵞ::Complex{RT},
    ℓ::IT
) where {IT<:Signed, RT, NT<:Complex}
    _recurrence!(w, eⁱᵝ, ℓ)
    let Hˡ = Hˡ(w)
        convert_H_to_D!(Hˡ, eⁱᵅ, eⁱᵞ)
        Hˡ
    end
end

function recurrence!(
    w::WignerHCalculator{IT, RT, NT}, eⁱᵝ::Complex{RT}, ℓ::IT
) where {IT<:Signed, RT, NT<:Real}
    _recurrence!(w, eⁱᵝ, ℓ)
    let Hˡ = Hˡ(w)
        convert_H_to_d!(Hˡ)
        Hˡ
    end
end

function _recurrence!(
    w::WignerHCalculator{IT, RT, NT}, eⁱᵝ::Complex{RT}, ℓ::IT
) where {IT<:Signed, RT, NT}
    # NOTE: In the comments explaining the recurrence steps below, we use notation with
    # ℓₘᵢₙ=0 for simplicity, but this sequence may work for ℓₘᵢₙ=1//2 as well.

    if ℓ == ℓₘᵢₙ(w)
        w.ℓ = ℓ
        recurrence_step1!(w)  # H⁰₀₀ = 1
        fillHˡ₀ₘ!(w)  # Record the result in Hˡ (it was computed in h⃗ˡ)

        # Now, to leave `w` in a good state for the next ℓ, compute H¹₀ₘ.
        recurrence_step2!(w, eⁱᵝ)  # H⁰₀ₘ -> H¹₀ₘ
    else
        if !consistent_ℓ(w) || ℓ(w) ≠ ℓ-1
            # We need to start over, but will only be using the axes, so we only reset them.
            h⃗ˡ(w).ℓ = ℓₘᵢₙ(w)
            h⃗ˡ⁺¹(w).ℓ = ℓₘᵢₙ(w)+1

            recurrence_step1!(w)  # H⁰₀₀ = 1
            recurrence_step2!(w, eⁱᵝ)  # H⁰₀ₘ -> H¹₀ₘ

            for ℓ′ in ℓₘᵢₙ(w)+1:ℓ-1
                increment_axes!(w)
                recurrence_step2!(w, eⁱᵝ)  # Hˡ₀ₘ -> Hˡ⁺¹₀ₘ
            end
        end

        # Do one more step of the recurrence to get Hˡ⁺¹₀ₘ, regardless of whether or not the
        # `if` block above was executed.  If not, it was set in a previous call to this
        # function, but is currently stored in h⃗ˡ⁺¹ so we have to swap and increment the
        # axes first.  If so, that block just got us to the same point.
        increment_axes!(w)
        recurrence_step2!(w, eⁱᵝ, ℓ)  # Hˡ₀ₘ -> Hˡ⁺¹₀ₘ

        # So far, we've only used h⃗ˡ and h⃗ˡ⁺¹.  Now we're ready to start using Hˡ.
        Hˡ(w).ℓ = ℓ

        # The h⃗ˡ and h⃗ˡ⁺¹ are what we need to start the recurrence for Hˡ
        fillHˡ₀ₘ!(w)  # Copy h⃗ˡ₀ₘ to Hˡ₀ₘ
        recurrence_step3!(w, eⁱᵝ)  # Hˡ⁺¹₀ₘ -> Hˡ₁ₘ
        recurrence_step4!(w, eⁱᵝ)  # Hˡₘ′ₘ₋₁, Hˡₘ′₋₁ₘ, Hˡₘ′ₘ₊₁ -> Hˡₘ′₊₁ₘ
        recurrence_step5!(w, eⁱᵝ)  # Hˡₘ′ₘ₋₁, Hˡₘ′₊₁ₘ, Hˡₘ′ₘ₊₁ -> Hˡₘ′₋₁ₘ

        # # Impose symmetries
        # recurrence_step6!(w, ℓ)

        # # Swap the H matrices once more so that the current Hˡ⁺¹ is the next loop's Hˡ
        # increment_axes!(w)
    end
    Hˡ(w)
end
