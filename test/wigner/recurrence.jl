# A dense form of the `H` recurrence, in which one `WignerMatrix` holds the whole `Hˡ` of
# one rotor, and the test that compares it with the engine that the package runs.
#
# The module `DenseRecurrence` defines `recurrence_step1!` … `recurrence_step6!`,
# `convert_H_to_d!`, and `convert_H_to_D!`, which take the steps of the notes on the `H`
# recursion one at a time, for integer indices.  The engine (`HCalculator`) runs steps 1
# to 5 as the package's functions of the same names in `src/recurrence/h_calculator.jl`,
# which work on a batched quarter-wedge; it never applies step 6, and the calculators'
# `materialize!` applies the phases of step 7.  The functions here are functions of the test
# module, not methods of the package's.  The dense form is an independent second
# implementation of the same recurrence (different loop structure, different storage, no
# batching), so that an error in either shows up as a disagreement between them.

@testmodule DenseRecurrence begin

import SphericalFunctions: WignerMatrix, ℓ, ℓₘᵢₙ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, sgn, ϵ,
    ComplexPowers

# The functions in this module loop over every m of a block's ℓ, and take m′ₘᵢₙ to be
# -m′ₘₐₓ, under `@inbounds`, and the `WignerMatrix` accessors are `@propagate_inbounds`, so
# the container's own bounds check is elided along with the storage's.  Each function
# therefore refuses, before it writes anything, a block that does not have the full range of
# m and a symmetric range of m′.  (The batched engine in `src/recurrence/h_calculator.jl`,
# which is what the package itself runs, has no such restriction.)
function check_dense_block(Hˡ::WignerMatrix, name::Symbol)
    if !(mₘᵢₙ(Hˡ) == -ℓ(Hˡ) && mₘₐₓ(Hˡ) == ℓ(Hˡ) && m′ₘᵢₙ(Hˡ) == -m′ₘₐₓ(Hˡ))
        throw(ArgumentError(
            "`$name` needs a block with the full range m ∈ -ℓ:ℓ and a symmetric range of "
            * "m′; this one has ℓ=$(ℓ(Hˡ)), m′ ∈ $(m′ₘᵢₙ(Hˡ)):$(m′ₘₐₓ(Hˡ)) and "
            * "m ∈ $(mₘᵢₙ(Hˡ)):$(mₘₐₓ(Hˡ))."
        ))
    end
    nothing
end

# Steps 2 and 3 combine blocks of two consecutive orders, the lower one first.
function check_consecutive_orders(lower::WignerMatrix, upper::WignerMatrix, name::Symbol)
    if ℓ(upper) != ℓ(lower) + 1
        throw(ArgumentError(
            "`$name` needs blocks of consecutive orders ℓ and ℓ+1; these have ℓ=$(ℓ(lower)) "
            * "and ℓ=$(ℓ(upper))."
        ))
    end
    nothing
end


@doc raw"""
    recurrence_step1!(H⁰)

Initialize the Wigner matrix `H⁰` for the recurrence relations.  This only sets the values
`H⁰[0,0]=1`.

Note that `H⁰` can be any `WignerMatrix` with integer indices — the only container indexed
by `(m′, m)`.  In particular, it can be a `D` matrix or a `d` matrix.  As for the other
dense reference functions (`recurrence_step2!` to `recurrence_step6!`, `convert_H_to_d!`,
and `convert_H_to_D!`), the block must have the full range of ``m`` and a symmetric range of
``m′``.
"""
function recurrence_step1!(H⁰::WignerMatrix{IT, NT}) where {IT<:Signed, NT}
    check_dense_block(H⁰, :recurrence_step1!)
    @inbounds let ℓ=ℓ(H⁰)
        if ℓ == 0
            H⁰[0, 0] = 1
        else
            throw(ArgumentError("Trying to initialize ℓ=$ℓ; only ℓ=0 is supported."))
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
    Hˡ::WignerMatrix{IT, NT}, Hˡ⁻¹::WignerMatrix{IT, NT2}, sinβ::T, cosβ::T
) where {IT<:Signed, NT, NT2, T}
    check_dense_block(Hˡ, :recurrence_step2!)
    check_dense_block(Hˡ⁻¹, :recurrence_step2!)
    check_consecutive_orders(Hˡ⁻¹, Hˡ, :recurrence_step2!)
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
            throw(ArgumentError(
                "Tried to recurse with ℓ=$ℓ; only integer ℓ ≥ 1 is supported."
            ))
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
    Hˡ::WignerMatrix{IT, NT}, Hˡ⁺¹::WignerMatrix{IT, NT2}, sinβ::T, cosβ::T
) where {IT<:Signed, NT, NT2, T}
    check_dense_block(Hˡ, :recurrence_step3!)
    check_dense_block(Hˡ⁺¹, :recurrence_step3!)
    check_consecutive_orders(Hˡ, Hˡ⁺¹, :recurrence_step3!)
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
``H^{ℓ}_{m'-1,m}``, and ``H^{ℓ}_{m',m+1}``, for all ``1 \leq m' < m'_{\mathrm{max}}`` and
``m \geq m'+1``.

"""
function recurrence_step4!(
    Hˡ::WignerMatrix{IT, NT}, sinβ::T, cosβ::T
) where {IT<:Signed, NT, T}
    check_dense_block(Hˡ, :recurrence_step4!)
    @inbounds let √=sqrt∘T, ℓ=ℓ(Hˡ), m′ₘₐₓ=m′ₘₐₓ(Hˡ)
        for m′ ∈ 1:min(ℓ, m′ₘₐₓ)-1
            # Note that the signs of m′ and m are always +1 for *integer* indices, so we
            # leave them out of the calculations of d̄ in this function.  They are not for
            # half-integer indices, where sgn(m′-1) = -1 at m′ = 1/2 (see step 4 of the
            # notes on the H recursion); this function is integer-only, and the batched
            # engine's `recurrence_step4!` applies the sign.
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
``H^{ℓ}_{m'+1,m}``, and ``H^{ℓ}_{m',m+1}``, for all ``m' \leq 0`` and ``m > -m'``.

"""
function recurrence_step5!(
    Hˡ::WignerMatrix{IT, NT}, sinβ::T, cosβ::T
) where {IT<:Signed, NT, T}
    check_dense_block(Hˡ, :recurrence_step5!)
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
    [``H`` recursion](@ref "Algorithm for computing ``H``")); only
    ``H^ℓ_{m′, m} = H^ℓ_{-m, -m′}`` is sign-free.  This
    function is therefore restricted to integer indices.  The batched engine never runs
    step 6 at all: every out-of-wedge element is read through
    [`wedge_source`](@ref), which already accounts for ``σ``.

"""
function recurrence_step6!(Hˡ::WignerMatrix{IT, NT}) where {IT<:Signed, NT}
    check_dense_block(Hˡ, :recurrence_step6!)
    @inbounds let ℓ=ℓ(Hˡ), m′ₘₐₓ=m′ₘₐₓ(Hˡ)
        # The idea here is to impose
        #   Hˡ[m, m′] = Hˡ[-m, -m′] = Hˡ[-m′, -m] = Hˡ[m′, m]
        # without double-counting any entries, and accounting for m′ₘₐₓ.
        for m ∈ 1:ℓ
            for m′ ∈ -min(m′ₘₐₓ, m):min(m′ₘₐₓ, m)
                Hˡ[-m′, -m] = Hˡ[m′, m]
            end
            # Rows ±m of the transposed region exist only when |m| ≤ m′ₘₐₓ.  Without this
            # guard, a matrix with 0 < m′ₘₐₓ < ℓ would write past the end of its own block,
            # and the enclosing `@inbounds` would turn that into silent corruption rather
            # than a `BoundsError`.  (The batched engine never runs this step: it reads
            # every element outside its symmetric wedge through `wedge_source`.)
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
function convert_H_to_d!(Hˡ::WignerMatrix{IT, NT}) where {IT<:Signed, NT<:Real}
    check_dense_block(Hˡ, :convert_H_to_d!)
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
    convert_H_to_D!(Hˡ, eⁱᵅ, eⁱᵞ)

Convert the Wigner matrix `Hˡ` to the D matrix `Dˡ`, which just involves multiplying by the
complex phases ``e^{-im′α}`` and ``e^{-imγ}``, given the phases `eⁱᵅ` and `eⁱᵞ`.

"""
function convert_H_to_D!(Hˡ::WignerMatrix{IT, NT}, eⁱᵅ::NT, eⁱᵞ::NT) where {IT<:Signed, NT<:Complex}
    # For half-integer indices this form does not apply, because e^{-im′α} and e^{-imγ} are
    # not integer powers of eⁱᵅ and eⁱᵞ.  No square roots are needed to fix that, though: m′
    # ± m *are* integers, so e^{i(m′α+mγ)} = z₊^{m′+m} z₋^{m′-m} with z₊ = e^{i(α+γ)/2}, z₋
    # = e^{i(α-γ)/2} (see step 7 of the notes on the H recursion).  That is what the batched
    # `materialize!` implements, for both index types at once.
    check_dense_block(Hˡ, :convert_H_to_D!)
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

end  # @testmodule DenseRecurrence


@testitem "Dense H recurrence vs the batched engine" setup=[RefusalChecks, DenseRecurrence] begin
    import SphericalFunctions as SF
    import SphericalFunctions: WignerMatrix, D, d, spinor_phases
    import .DenseRecurrence:
        recurrence_step1!, recurrence_step2!, recurrence_step3!,
        recurrence_step4!, recurrence_step5!, recurrence_step6!,
        convert_H_to_d!, convert_H_to_D!
    using Quaternionic: Rotor, 𝐢, 𝐣, 𝐤
    using Random: Xoshiro

    # Drive the six steps by hand, exactly as the notes on the `H` recursion describe them,
    # and return the filled `Hˡ`.  `axes_[n+1]` holds the m′=0 axis of `Hⁿ`; step 3 needs
    # the axis one order *above* the target ℓ.
    function dense_H(::Type{NT}, ℓ, cosβ, sinβ; m′ₘₐₓ=ℓ) where {NT}
        axes_ = [WignerMatrix(zeros(NT, 2n + 1, 2n + 1), n) for n ∈ 0:ℓ+1]
        recurrence_step1!(axes_[1])                               # H⁰₀₀ = 1
        for n ∈ 1:ℓ+1
            recurrence_step2!(axes_[n+1], axes_[n], sinβ, cosβ)   # Hⁿ⁻¹₀ₘ -> Hⁿ₀ₘ
        end
        Hˡ = WignerMatrix(
            zeros(NT, 2m′ₘₐₓ + 1, 2ℓ + 1), ℓ; m′ₘₐₓ=m′ₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ
        )
        for m ∈ 0:ℓ
            Hˡ[0, m] = axes_[ℓ+1][0, m]
        end
        recurrence_step3!(Hˡ, axes_[ℓ+2], sinβ, cosβ)             # Hˡ⁺¹₀ₘ -> Hˡ₁ₘ
        recurrence_step4!(Hˡ, sinβ, cosβ)                         # ... -> Hˡₘ′₊₁ₘ
        recurrence_step5!(Hˡ, sinβ, cosβ)                         # ... -> Hˡₘ′₋₁ₘ
        recurrence_step6!(Hˡ)                                     # the symmetries
        Hˡ
    end

    # The poles (where the recurrence's special branches live), both signs of the double
    # cover, and a few seeded random rotors.
    rotors(::Type{T}) where {T} = [
        Rotor{T}(1); Rotor{T}(𝐢); Rotor{T}(𝐣); Rotor{T}(𝐤);
        -Rotor{T}(1); -Rotor{T}(𝐣);
        randn(Xoshiro(1729), Rotor{T}, 4)
    ]

    @testset "$T" for T ∈ (Float64, BigFloat)
        ℓₘₐₓ = 6
        # Measured worst case over everything below: 3.3e-16 (d) and 4.9e-16 (𝔇) in
        # Float64, i.e. about 2 eps; the two implementations differ only in rounding.
        atol = 20 * eps(T)
        for R ∈ rotors(T)
            eⁱᵝ, z₊, z₋, _, _ = spinor_phases(R)
            cosβ, sinβ = reim(eⁱᵝ)
            # e^{iα} = z₊ z₋ and e^{iγ} = z₊ conj(z₋), from z₊ = e^{i(α+γ)/2},
            # z₋ = e^{i(α-γ)/2}; taken this way, no Euler angle is ever extracted.
            eⁱᵅ, eⁱᵞ = z₊ * z₋, z₊ * conj(z₋)
            for ℓ ∈ 0:ℓₘₐₓ
                dref = d(eⁱᵝ, ℓ)[ℓ]
                𝔇ref = D(R, ℓ)[ℓ]
                # Errors are accumulated and asserted once per (rotor, ℓ): a per-element
                # `@test` would be tens of thousands of assertions.
                errd = zero(T)
                err𝔇 = zero(T)
                for m′ₘₐₓ ∈ 0:ℓ
                    Hd = dense_H(T, ℓ, cosβ, sinβ; m′ₘₐₓ)
                    convert_H_to_d!(Hd)
                    H𝔇 = dense_H(Complex{T}, ℓ, cosβ, sinβ; m′ₘₐₓ)
                    convert_H_to_D!(H𝔇, eⁱᵅ, eⁱᵞ)
                    for m′ ∈ -m′ₘₐₓ:m′ₘₐₓ, m ∈ -ℓ:ℓ
                        errd = max(errd, abs(Hd[m′, m] - dref[m′, m]))
                        err𝔇 = max(err𝔇, abs(H𝔇[m′, m] - 𝔇ref[m′, m]))
                    end
                end
                @test errd ≤ atol
                @test err𝔇 ≤ atol
            end
        end
    end

    # `H` itself, before the phases: the symmetries step 6 imposes must hold exactly, since
    # they are assignments from one stored element to another.
    let T = Float64, ℓ = 5
        eⁱᵝ, = spinor_phases(randn(Xoshiro(11), Rotor{T}))
        cosβ, sinβ = reim(eⁱᵝ)
        H = dense_H(T, ℓ, cosβ, sinβ)
        for m′ ∈ -ℓ:ℓ, m ∈ -ℓ:ℓ
            @test H[m′, m] == H[m, m′]
            @test H[m′, m] == H[-m′, -m]
        end
    end

    # Step 1 only initializes ℓ=0, and steps 2 and 3 combine blocks of consecutive orders.
    @test refuses(
        () -> recurrence_step1!(WignerMatrix(zeros(3, 3), 1)), ArgumentError,
        "only ℓ=0 is supported"
    )
    let H⁰ = WignerMatrix(zeros(1, 1), 0), H² = WignerMatrix(zeros(5, 5), 2)
        @test refuses(
            () -> recurrence_step2!(H², H⁰, 0.1, 0.9), ArgumentError, "consecutive orders"
        )
        @test refuses(
            () -> recurrence_step3!(H⁰, H², 0.1, 0.9), ArgumentError, "consecutive orders"
        )
    end
end
