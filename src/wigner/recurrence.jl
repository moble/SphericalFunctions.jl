# `sgn`, `ϵ` and `δ²` are defined in `wigner_H.jl`.


@doc raw"""
    recurrence_step1!(H⁰)

Initialize the Wigner matrix `H⁰` for the recurrence relations.  This only sets the values
`H⁰[0,0]=1`.

Note that `H⁰` can be any `AbstractWignerMatrix` with integer indices.  In particular, it
can be a `D` matrix or a `d` matrix.
"""
function recurrence_step1!(H⁰::AbstractWignerMatrix{IT, NT}) where {IT<:Signed, NT}
    @inbounds let ℓ=ℓ(H⁰)
        if ℓ == 0
            H⁰[0, 0] = 1
        else
            error("Trying to initialize ℓ=$ℓ; only ℓ=0 is supported.")
        end
    end
    H⁰
end

@doc raw"""
    recurrence_step2!(Hˡ, Hˡ⁻¹, sinβ, cosβ)

Compute the values of ``H^{ℓ}_{0,m}``, from the values of ``H^{ℓ-1}_{0,m}`` for all
``m \geq 0``.

"""
function recurrence_step2!(
    Hˡ::AbstractWignerMatrix{IT, NT}, Hˡ⁻¹::AbstractWignerMatrix{IT, NT2}, sinβ::T, cosβ::T
) where {IT<:Signed, NT, NT2, T}
    @assert ℓ(Hˡ⁻¹) == ℓ(Hˡ) - 1
    # Note that in this step only, we use notation derived from Xing et al., denoting the
    # coefficients as b̄ₗ, c̄ₗₘ, d̄ₗₘ, ēₗₘ.  In the following steps, we will use notation
    # from Gumerov and Duraiswami, who denote their different coefficients aₗᵐ, etc.
    @inbounds let √=sqrt∘T, ℓ=ℓ(Hˡ)
        if ℓ == 1
            # The ℓ>1 branch would try to access invalid indices of H⁰; if we treat those
            # elements as zero, we can simplify that branch to just the following much
            # simpler code anyway.  So fundamentally, this branch is the same as the other
            # branch.
            Hˡ[0, 0] = cosβ
            Hˡ[0, 1] = sinβ / √2
        elseif ℓ > 1
            b̄ₗ = √(T(ℓ-1)/ℓ)
            Hˡ[0, 0] = cosβ * Hˡ⁻¹[0, 0] - b̄ₗ * sinβ * Hˡ⁻¹[0, 1]
            for m ∈ 1:ℓ-2
                c̄ₗₘ = √((ℓ+m)*(ℓ-m)) / ℓ
                d̄ₗₘ = √((ℓ-m)*(ℓ-m-1)) / 2ℓ
                ēₗₘ = √((ℓ+m)*(ℓ+m-1)) / 2ℓ
                Hˡ[0, m] = (
                    c̄ₗₘ * cosβ * Hˡ⁻¹[0, m]
                    - sinβ * (d̄ₗₘ * Hˡ⁻¹[0, m+1] - ēₗₘ * Hˡ⁻¹[0, m-1])
                )
            end
            let m = ℓ-1
                c̄ₗₘ = √((ℓ+m)*(ℓ-m)) / ℓ
                ēₗₘ = √((ℓ+m)*(ℓ+m-1)) / 2ℓ
                Hˡ[0, m] = (
                    c̄ₗₘ * cosβ * Hˡ⁻¹[0, m]
                    - sinβ * (- ēₗₘ * Hˡ⁻¹[0, m-1])
                )
            end
            let m = ℓ
                ēₗₘ = √((ℓ+m)*(ℓ+m-1)) / 2ℓ
                Hˡ[0, m] = (
                    - sinβ * (- ēₗₘ * Hˡ⁻¹[0, m-1])
                )
            end
        else
            error("Tried to recurse with ℓ=$ℓ; only integer ℓ ≥ 1 is supported.")
        end
    end
    Hˡ
end

@doc raw"""
    recurrence_step3!(Hˡ, Hˡ⁺¹, sinβ, cosβ)

Compute the values of ``H^{ℓ}_{1,m}``, from the values of ``H^{ℓ+1}_{0,m}`` for all
``m \geq 0``.

"""
function recurrence_step3!(
    Hˡ::AbstractWignerMatrix{IT, NT}, Hˡ⁺¹::AbstractWignerMatrix{IT, NT2}, sinβ::T, cosβ::T
) where {IT<:Signed, NT, NT2, T}
    @assert ℓ(Hˡ⁺¹) == ℓ(Hˡ) + 1
    @inbounds let √=sqrt∘T, ℓ=ℓ(Hˡ), m′ₘₐₓ=m′ₘₐₓ(Hˡ)
        if ℓ > 0 && m′ₘₐₓ ≥ 1
            c = 1 / √(ℓ*(ℓ+1))
            for m ∈ 1:ℓ
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
    Hˡ
end

@doc raw"""
    recurrence_step4!(Hˡ, sinβ, cosβ)

Compute the values of ``H^{ℓ}_{m'+1,m}``, from the values of ``H^{ℓ}_{m',m-1}``,
``H^{ℓ}_{m'-1,m}``, and  ``H^{ℓ}_{m',m+1}``, for all ``m' > 1`` and ``m \geq m'``.

"""
function recurrence_step4!(
    Hˡ::AbstractWignerMatrix{IT, NT}, sinβ::T, cosβ::T
) where {IT<:Signed, NT, T}
    @inbounds let √=sqrt∘T, ℓ=ℓ(Hˡ), m′ₘₐₓ=m′ₘₐₓ(Hˡ)
        for m′ ∈ 1:min(ℓ, m′ₘₐₓ)-1
            # Note that the signs of m′ and m are always +1 for *integer* indices, so we
            # leave them out of the calculations of d̄ in this function.  They are not for
            # half-integer indices, where sgn(m′-1) = -1 at m′ = 1/2 (see the v3 design
            # memo, §5.2); this function is integer-only, and the batched engine's
            # `recurrence_step4!` applies the sign.
            d̄ₗᵐ′ = √((ℓ-m′)*(ℓ+m′+1))
            d̄ₗᵐ′⁻¹ = √((ℓ-m′+1)*(ℓ+m′))
            for m ∈ (m′+1):ℓ-1
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
    Hˡ
end

@doc raw"""
    recurrence_step5!(Hˡ, sinβ, cosβ)

Compute the values of ``H^{ℓ}_{m'-1,m}``, from the values of ``H^{ℓ}_{m',m-1}``,
``H^{ℓ}_{m'+1,m}``, and  ``H^{ℓ}_{m',m+1}``, for all ``m' \leq 0`` and ``m > -m'``.

"""
function recurrence_step5!(
    Hˡ::AbstractWignerMatrix{IT, NT}, sinβ::T, cosβ::T
) where {IT<:Signed, NT, T}
    @inbounds let √=sqrt∘T, ℓ=ℓ(Hˡ), m′ₘᵢₙ=m′ₘᵢₙ(Hˡ)
        for m′ ∈ 0:-1:max(-ℓ, m′ₘᵢₙ)+1
            d̄ₗᵐ′ = sgn(m′) * √((ℓ-m′)*(ℓ+m′+1))
            d̄ₗᵐ′⁻¹ = sgn(m′-1) * √((ℓ-m′+1)*(ℓ+m′))
            for m ∈ -(m′-1):ℓ-1
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
    Hˡ
end

@doc raw"""
    recurrence_step6!(Hˡ)

Impose the symmetries of the Wigner matrix `Hˡ` to fill in all the values that have not yet
been computed.

Assuming that `Hˡ` has already been computed as much as possible by the recurrence
relations, this function imposes the symmetries, rather than recalculating terms.
Specifically, steps 1–5 fill the wedge `m ≥ abs(m′)` for `abs(m′) ≤ m′ₘₐₓ`, and this
function completes the requested block using the symmetries
```math
\begin{aligned}
H^ℓ_{m′, m} &= H^ℓ_{m, m′}, \\
H^ℓ_{m′, m} &= H^ℓ_{-m′, -m}.
\end{aligned}
```

!!! note "Half-integer indices"

    Both of those symmetries acquire the sign ``σ = \mathrm{sgn}(m)\,\mathrm{sgn}(m')`` for
    half-integer indices (see [`transpose_sign`](@ref) and the notes on the
    [``H`` recursion](@ref "Algorithm for computing ``H`` (redesigned)")); only
    ``H^ℓ_{m′, m} = H^ℓ_{-m, -m′}`` is sign-free.  This
    function is therefore restricted to integer indices.  The batched engine never runs
    step 6 at all: every out-of-wedge element is read through
    [`wedge_source`](@ref), which already accounts for ``σ``.

"""
function recurrence_step6!(Hˡ::AbstractWignerMatrix{IT, NT}) where {IT<:Signed, NT}
    @inbounds let ℓ=ℓ(Hˡ), m′ₘₐₓ=m′ₘₐₓ(Hˡ)
        # The idea here is to impose
        #   Hˡ[m, m′] = Hˡ[-m, -m′] = Hˡ[-m′, -m] = Hˡ[m′, m]
        # without double-counting any entries, and accounting for m′ₘₐₓ.
        for m ∈ 1:ℓ
            for m′ ∈ -min(m′ₘₐₓ, m):min(m′ₘₐₓ, m)
                Hˡ[-m′, -m] = Hˡ[m′, m]
            end
            # Rows ±m of the transposed region exist only when |m| ≤ m′ₘₐₓ.  Without this
            # guard, a matrix with 0 < m′ₘₐₓ < ℓ writes past the end of its own block — and
            # the enclosing `@inbounds` turns that into silent corruption rather than a
            # `BoundsError` (the v3 design memo, bug B2, which was fixed in the batched
            # engine by always computing the symmetric wedge).
            if m ≤ m′ₘₐₓ
                for m′ ∈ -min(m′ₘₐₓ, m-1):min(m′ₘₐₓ, m-1)
                    Hˡ[m, m′] = Hˡ[-m, -m′] = Hˡ[m′, m]
                end
            end
        end
    end
    Hˡ
end


"""
    convert_H_to_d!(Hˡ)

Convert the Wigner matrix `Hˡ` to the d matrix `dˡ`, which just involves multiplying by
signs related to the `m′` and `m` indices.

"""
function convert_H_to_d!(Hˡ::AbstractWignerMatrix{IT, NT}) where {IT<:Signed, NT<:Real}
    @inbounds let ℓ=ℓ(Hˡ), m′ₘₐₓ=m′ₘₐₓ(Hˡ)
        for m ∈ -ℓ:ℓ
            for m′ ∈ -m′ₘₐₓ:m′ₘₐₓ
                Hˡ[m′, m] *= ϵ(m′) * ϵ(-m)
            end
        end
    end
    Hˡ
end


"""
    convert_H_to_D!(Hˡ)

Convert the Wigner matrix `Hˡ` to the D matrix `Dˡ`, which just involves multiplying by
complex phases related to the `m′` and `m` indices.

"""
function convert_H_to_D!(Hˡ::AbstractWignerMatrix{IT, NT}, eⁱᵅ::NT, eⁱᵞ::NT) where {IT<:Signed, NT<:Complex}
    # For half-integer indices this form does not apply, because e^{-im′α} and e^{-imγ} are
    # not integer powers of eⁱᵅ and eⁱᵞ.  No square roots are needed to fix that, though:
    # m′ ± m *are* integers, so e^{i(m′α+mγ)} = z₊^{m′+m} z₋^{m′-m} with
    # z₊ = e^{i(α+γ)/2}, z₋ = e^{i(α-γ)/2} (the v3 design memo, §5.4).  That is what the
    # batched `materialize!` implements, for both index types at once.
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
    Hˡ
end
