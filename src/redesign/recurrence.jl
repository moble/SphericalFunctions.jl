# Eq. (44) in Gumerov and Duraiswami (2015).  Note that they define `sgn` as follows, which
# is different from the usual definition, including from Julia's `sign` function, at 0:
sgn(m) = m ≥ 0 ? 1 : -1

# Eq. (7) in Gumerov and Duraiswami (2015)
ϵ(m) = (m ≥ 0 ? (-1)^m : 1)


@doc raw"""
    initialize!(Hˡ, sinβ, cosβ)

Step 1 of the computation of ``H``: Initialize the Wigner matrix `Hˡ` for the recurrence
relations.  This only sets the values ``H^0_{0,0}=0``, ``H^1_{0,0}=cosβ``, and
``H^1_{0,1}=sinβ/√2``.

Note that `Hˡ` can be any `WignerMatrix` with integer indices.  In particular, it can be a
`D` matrix or a `d` matrix.
"""
function initialize!(Hˡ::WignerMatrix{NT, IT}, sinβ::T, cosβ::T) where {NT, IT<:Signed, T}
    @inbounds let √=sqrt∘T, ℓ=ℓ(Hˡ)
        if ℓ == 0
            Hˡ[0, 0] = 0
        elseif ℓ == 1
            Hˡ[0, 0] = cosβ
            Hˡ[0, 1] = sinβ / √2
        end
    end
end

@doc raw"""
    recurrence_0_m!(Hˡ, Hˡ⁻¹, sinβ, cosβ)

Step 2 of the computation of ``H``: Given ``H^{\ell-1}`` with its ``(0, m)`` entries,
compute ``H^{\ell}_{0,m}`` for all ``m \geq 0``.

"""
function recurrence_0_m!(
    Hˡ::WignerMatrix{NT, IT}, Hˡ⁻¹::WignerMatrix{NT, IT}, sinβ::T, cosβ::T
) where {NT, IT<:Signed, T}
    @assert ℓ(Hˡ⁻¹) == ℓ(Hˡ) - 1
    # Note that in this step only, we use notation derived from Xing et al., denoting the
    # coefficients as b̄ₗ, c̄ₗₘ, d̄ₗₘ, ēₗₘ.  In the following steps, we will use notation
    # from Gumerov and Duraiswami, who denote their different coefficients aₗᵐ, etc.
    @inbounds let √=sqrt∘T, ℓ=ℓ(Hˡ)
        if ℓ > 1
            b̄ₗ = √(T(ℓ-1)/ℓ)
            Hˡ[0, 0] = cosβ * Hˡ⁻¹[0, 0] - b̄ₗ * sinβ * Hˡ⁻¹[0, 1]
            for m ∈ 1:ℓ-2
                c̄ₙₘ = √((ℓ+m)*(ℓ-m)) / ℓ
                d̄ₙₘ = √((ℓ-m)*(ℓ-m-1)) / 2ℓ
                ēₙₘ = √((ℓ+m)*(ℓ+m-1)) / 2ℓ
                Hˡ[0, m] = (
                    c̄ₗₘ * cosβ * Hˡ⁻¹[0, m]
                    - sinβ * (d̄ₗₘ * Hˡ⁻¹[0, m+1] - ēₗₘ * Hˡ⁻¹[0, m-1])
                )
            end
            let m = ℓ-1
                c̄ₙₘ = √((ℓ+m)*(ℓ-m)) / ℓ
                ēₙₘ = √((ℓ+m)*(ℓ+m-1)) / 2ℓ
                Hˡ[0, m] = (
                    c̄ₗₘ * cosβ * Hˡ⁻¹[0, m]
                    - sinβ * (- ēₗₘ * Hˡ⁻¹[0, m-1])
                )
            end
            let m = ℓ
                ēₙₘ = √((ℓ+m)*(ℓ+m-1)) / 2ℓ
                Hˡ[0, m] = (
                    - sinβ * (- ēₗₘ * Hˡ⁻¹[0, m-1])
                )
            end
        end
    end
end

@doc raw"""
    recurrence_1_m!(Hˡ, Hˡ⁺¹, sinβ, cosβ)

Step 3 of the computation of ``H``: Given ``H^{\ell+1}`` with all its ``(0, m)`` entries,
compute ``H^{\ell}_{1,m}`` for all ``m \geq 1``.

"""
function recurrence_1_m!(
    Hˡ::WignerMatrix{NT, IT}, Hˡ⁺¹::WignerMatrix{NT, IT}, sinβ::T, cosβ::T
) where {NT, IT<:Signed, T}
    @assert ℓ(Hˡ⁺¹) == ℓ(Hˡ) + 1
    @inbounds let √=sqrt∘T, ℓ=ℓ(Hˡ), m′ₘₐₓ=m′ₘₐₓ(Hˡ)
        if ℓ > 0 && m′ₘₐₓ ≥ 1
            c = 1 / √(ℓ*(ℓ+1))
            for m ∈ 0:ℓ
                āₗᵐ = √((ℓ+m+1)*(ℓ-m+1))
                b̄ₗ₊₁ᵐ⁻¹ = √((ℓ-m+1)*(ℓ-m+2))
                b̄ₗ₊₁⁻ᵐ⁻¹ = √((ℓ+m+1)*(ℓ+m+2))
                Hˡ[1, m] = -c * (
                    b̄ₗ₊₁⁻ᵐ⁻¹ * (1 - cosβ) / 2 * Hˡ⁺¹[0, m+1]
                    + b̄ₗ₊₁ᵐ⁻¹ * (1 + cosβ) / 2 * Hˡ⁺¹[0, m-1]
                    + āₗᵐ * sinβ * Hˡ⁺¹[0, m]
                )
            end
        end
    end
end

@doc raw"""
    recurrence_m′₊!(Hˡ, sinβ, cosβ)

Step 4 of the computation of ``H``: Given ``H^{\ell}`` with all its ``(0, m)`` and ``(1,
m)`` entries, compute ``H^{\ell}_{m′+1,m}`` for all ``m′ \geq 1`` and ``m \geq 1``.

"""
function recurrence_m′₊!(
    Hˡ::WignerMatrix{NT, IT}, sinβ::T, cosβ::T
) where {NT, IT<:Signed, T}
    @inbounds let √=sqrt∘T, ℓ=ℓ(Hˡ), m′ₘₐₓ=m′ₘₐₓ(Hˡ)
        for m′ ∈ 1:min(ℓ, m′ₘₐₓ)-1
            # Note that the signs of m′ and m are always +1, so we leave them out of the
            # calculations of d̄ in this function.
            d̄ₗᵐ′ = √((ℓ-m′)*(ℓ+m′+1))
            d̄ₗᵐ′⁻¹ = √((ℓ-m′+1)*(ℓ+m′))
            for m ∈ m′:ℓ-1
                d̄ₗᵐ⁻¹ = √((ℓ-m+1)*(ℓ+m))
                d̄ₗᵐ = √((ℓ-m)*(ℓ+m+1))
                Hˡ[m′+1, m] = (
                    d̄ₗᵐ′⁻¹ * Hˡ[m′-1, m]
                    - d̄ₗᵐ⁻¹ * Hˡ[m′, m-1]
                    + d̄ₗᵐ * Hˡ[m′, m+1]
                ) / d̄ₗᵐ′
            end
            let m = ℓ
                d̄ₗᵐ⁻¹ = √((ℓ-m+1)*(ℓ+m))
                Hˡ[m′+1, m] = (
                    d̄ₗᵐ′⁻¹ * Hˡ[m′-1, m]
                    - d̄ₗᵐ⁻¹ * Hˡ[m′, m-1]
                ) / d̄ₗᵐ′
            end
        end
    end
end

@doc raw"""
    recurrence_m′₋!(Hˡ, sinβ, cosβ)

Step 5 of the computation of ``H``: Given ``H^{\ell}`` with all its ``(m′, m)`` entries for
`m′ ≥ 0`, compute ``H^{\ell}_{m′-1,m}`` for all `m′ ≤ -1` and `m`.

"""
function recurrence_m′₋!(
    Hˡ::WignerMatrix{NT, IT}, sinβ::T, cosβ::T
) where {NT, IT<:Signed, T}
    @inbounds let √=sqrt∘T, ℓ=ℓ(Hˡ), m′ₘₐₓ=m′ₘₐₓ(Hˡ)
        for m′ ∈ 0:-1:-min(ℓ, m′ₘₐₓ)+1
            d̄ₗᵐ′ = sgn(m′) * √((ℓ-m′)*(ℓ+m′+1))
            d̄ₗᵐ′⁻¹ = sgn(m′-1) * √((ℓ-m′+1)*(ℓ+m′))
            for m ∈ -m′:ℓ-1
                d̄ₗᵐ = sgn(m) * √((ℓ-m)*(ℓ+m+1))
                d̄ₗᵐ⁻¹ = sgn(m-1) * √((ℓ-m+1)*(ℓ+m))
                Hˡ[m′-1, m] = (
                    d̄ₗᵐ′ * Hˡ[m′+1, m]
                    + d̄ₗᵐ⁻¹ * Hˡ[m′, m-1]
                    - d̄ₗᵐ * Hˡ[m′, m+1]
                ) / d̄ₗᵐ′⁻¹
            end
            let m = ℓ
                d̄ₗᵐ⁻¹ = sgn(m-1) * √((ℓ-m+1)*(ℓ+m))
                Hˡ[m′-1, m] = (
                    d̄ₗᵐ′ * Hˡ[m′+1, m]
                    + d̄ₗᵐ⁻¹ * Hˡ[m′, m-1]
                ) / d̄ₗᵐ′⁻¹
            end
        end
    end
end

@doc raw"""
    impose_symmetries!(Hˡ)

Assuming that `Hˡ` has already been computed as much as possible by the recurrence
relations, this function imposes the symmetries, rather than recalculating terms.
Specifically, the recurrence relations will calculate the terms for all `m`, and
`m′ ≥ abs(m)`, and this function will complete the calculations using the symmetries
```math
\begin{aligned}
H^ℓ_{m′, m} &= H^ℓ_{m, m′}, \\
H^ℓ_{m′, m} &= H^ℓ_{-m′, -m}.
\end{aligned}
```

"""
function impose_symmetries!(Hˡ::WignerMatrix{NT, IT}) where {NT, IT<:Signed}
    @inbounds let ℓ=ℓ(Hˡ), m′ₘₐₓ=m′ₘₐₓ(Hˡ)
        for m ∈ -ℓ:ℓ
            for m′ ∈ abs(m):m′ₘₐₓ
                Hˡ[m, m′] = Hˡ[-m, -m′] = Hˡ[-m′, -m] = Hˡ[m′, m]
            end
        end
    end
end


"""
    convert_H_to_d!(Hˡ)

Convert the Wigner matrix `Hˡ` to the d matrix `dˡ`, which just involves multiplying by
signs related to the `m′` and `m` indices.

"""
function convert_H_to_d!(Hˡ::WignerMatrix{NT, IT}) where {NT, IT<:Signed}
    @inbounds let ℓ=ℓ(Hˡ), m′ₘₐₓ=m′ₘₐₓ(Hˡ)
        for m ∈ -ℓ:ℓ
            for m′ ∈ -m′ₘₐₓ:m′ₘₐₓ
                Hˡ[m′, m] *= ϵ(m′) * ϵ(-m)
            end
        end
    end
end


"""
    convert_H_to_D!(Hˡ)

Convert the Wigner matrix `Hˡ` to the D matrix `Dˡ`, which just involves multiplying by
complex phases related to the `m′` and `m` indices.

"""
function convert_H_to_D!(Hˡ::WignerMatrix{NT, IT}, eⁱᵅ::NT, eⁱᵞ::NT) where {NT, IT<:Signed}
    # NOTE: This function will have to be modified to work for Rational indices because the
    # phases will not be integer powers; we'll have to incorporate √eⁱᵅ and √eⁱᵞ.
    @inbounds let ℓ=ℓ(Hˡ), ℓₘᵢₙ=ℓₘᵢₙ(Hˡ), m′ₘₐₓ=m′ₘₐₓ(Hˡ)
        ϕᵞ = ComplexPowers(eⁱᵞ)
        ϕᵅ = ComplexPowers(eⁱᵅ)
        for (m, eⁱᵐᵞ) ∈ zip(ℓₘᵢₙ:ℓ, ϕᵞ)
            for (m′, eⁱᵐ′ᵅ) ∈ zip(ℓₘᵢₙ:m′ₘₐₓ, ϕᵅ)
                Hˡ[m′, m] *= ϵ(m′) * ϵ(-m) * conj(eⁱᵐ′ᵅ) * conj(eⁱᵐᵞ)
                if m′ ≠ 0
                    Hˡ[-m′, m] *= ϵ(-m′) * ϵ(-m) * eⁱᵐ′ᵅ * conj(eⁱᵐᵞ)
                    if m ≠ 0
                        Hˡ[-m′, -m] *= ϵ(-m′) * ϵ(m) * eⁱᵐ′ᵅ * eⁱᵐᵞ
                    end
                end
                if m ≠ 0
                    Hˡ[m′, -m] *= ϵ(m′) * ϵ(m) * conj(eⁱᵐ′ᵅ) * eⁱᵐᵞ
                end
            end
        end
    end
end
