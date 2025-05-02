# Eq. (44) in Gumerov and Duraiswami (2015).  Note that they define `sgn` as follows, which
# is different from the usual definition, including from Julia's `sign` function, at 0:
sgn(m) = m ≥ 0 ? 1 : -1

# Eq. (7) in Gumerov and Duraiswami (2015)
ϵ(m) = (m ≥ 0 ? (-1)^m : 1)


function initialize!(Hˡ::WignerMatrix{NT, IT}, sinβ::T, cosβ::T) where {NT, IT<:Integer, T}
    @inbounds let √=sqrt∘T, ℓ=ℓ(Hˡ)
        if ℓ == 0
            Hˡ[0, 0] = 0
        elseif ℓ == 1
            Hˡ[0, 0] = cosβ
            Hˡ[0, 1] = sinβ / √2
        end
    end
end

function recurrence_0_m!(
    Hˡ::WignerMatrix{NT, IT}, Hˡ⁻¹::WignerMatrix{NT, IT}, sinβ::T, cosβ::T
) where {NT, IT<:Integer, T}
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

function recurrence_1_m!(
    Hˡ::WignerMatrix{NT, IT}, Hˡ⁺¹::WignerMatrix{NT, IT}, sinβ::T, cosβ::T
) where {NT, IT<:Integer, T}
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

function recurrence_m′₊!(
    Hˡ::WignerMatrix{NT, IT}, sinβ::T, cosβ::T
) where {NT, IT<:Integer, T}
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

function recurrence_m′₋!(
    Hˡ::WignerMatrix{NT, IT}, sinβ::T, cosβ::T
) where {NT, IT<:Integer, T}
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
function impose_symmetries!(Hˡ::WignerMatrix{NT, IT}, cisα::NT, cisβ::NT) where {NT, IT<:Integer}
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
function convert_H_to_d!(Hˡ::WignerMatrix{NT, IT}) where {NT, IT<:Integer}
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
function convert_H_to_D!(Hˡ::WignerMatrix{NT, IT}) where {NT, IT<:Integer}
    @inbounds let ℓ=ℓ(Hˡ), m′ₘₐₓ=m′ₘₐₓ(Hˡ)
        throw(ErrorException("convert_H_to_D! is not yet implemented"))
        for m ∈ -ℓ:ℓ
            for m′ ∈ -m′ₘₐₓ:m′ₘₐₓ
                Hˡ[m′, m] *= ϵ(m′) * ϵ(-m)
            end
        end
    end
end
