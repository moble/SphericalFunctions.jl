function initialize!(w::WignerMatrix{NT, IT}, sinβ, cosβ) where {NT, IT}
    T = real(NT)
    let √ = sqrt ∘ T
        if ℓ(w) == 0
            w[0, 0] = 0
        elseif ℓ(w) == 1
            w[0, 0] = cosβ
            w[0, 1] = sinβ / √2
        end
    end
end
