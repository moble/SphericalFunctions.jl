# Stage-1 oracle for the half-integer Wigner tests (design memo §6).
#
# Two *independent* reference implementations of the Wigner 𝔇 and d functions that are
# valid for half-integer indices, plus the transcribed tables from Varshalovich et al.
# Both are copied verbatim from the Literate comparison pages under
# `docs/literate_input/conventions/comparisons/`, which are finished, reviewed work; the
# copies exist because `@testmodule` cannot depend on another `@testmodule`, so the
# `Boyle2016` and `Varshalovich` modules defined on those pages cannot be imported here.
# Each copy names its source page and the line range it was taken from.  The
# "half-integer oracle: references agree with each other" test item in
# `test/wigner/half_integer.jl` re-derives, against these copies, the cross-checks that
# the comparison pages make, so that a drift between copy and original shows up as a
# failure here rather than as a silently wrong oracle.
#
# Conventions.  Varshalovich's `d(J, M, M′, β)` is this package's `dᴶ_{m′m}(β)` with
# `(M, M′) = (m′, m)`, and his `D(J, M, M′, α, β, γ)` is this package's `𝔇ᴶ_{m′m}`; the
# Boyle (2016) `WignerDElement(R, ℓ, m′, m)` is the *complex conjugate* of this package's
# `𝔇ˡ_{m′m}(R)`.  Both relations are established (for integer indices, and for
# half-integers between the two references) by the test items on the two comparison
# pages; `d_oracle` and `D_oracle` below package them up.

@testmodule HalfIntegerOracle begin

using Quaternionic: Rotor, from_euler_angles
import Random


# ---------------------------------------------------------------------------------------
# Boyle (2016), Eq. (35) and Appendix A.
# Copied verbatim from docs/literate_input/conventions/comparisons/boyle_2016.jl:59-167.
# ---------------------------------------------------------------------------------------

function WignerDElement(R::Rotor{T}, ℓ::I, m′::I, m::I) where {T, I}
    ## If `I` is Rational, check that the denominators are 2
    if I <: Rational
        if (denominator(ℓ) != 2) || (denominator(m′) != 2) || (denominator(m) != 2)
            error("The indices ℓ, m′, and m must all be integers or all be half-integers")
        end
    end

    ## Convert to twice the input values for half-integer support
    L  = Int(2ℓ)
    M′ = Int(2m′)
    M  = Int(2m)

    if L > 16
        error(
            "The maximum supported ℓ for this function is 8; " *
            "larger numbers become numerically unstable.\n" *
            "Consider using the `WignerD` function instead."
        )
    end

    ## Simple helper for 0+0im
    zeroCT = zero(Complex{T})

    if abs(M′) > L || abs(M) > L
        return zeroCT
    end

    let π = T(π)
        ## Split input `R` into its two complex components and extract magnitude and phase
        Rₛ = Complex(R[1], R[4])
        Rₐ = Complex(R[3], R[2])
        rₛ = abs(Rₛ)
        rₐ = abs(Rₐ)
        ϕₛ = angle(Rₛ)
        ϕₐ = angle(Rₐ)

        ## Check simple limiting cases
        if rₐ ≤ 4eps(rₛ)
            if M′ != M
                return zeroCT
            else
                return cis(M * ϕₛ)
            end

        elseif rₛ ≤ 4eps(rₐ)
            if -M′ != M
                return zeroCT
            else
                return cis(M * ϕₐ) * (((L - M) % 4 == 0) ? 1 : -1)
            end

        elseif rₐ ≤ rₛ
            λ = -(rₐ/rₛ)^2
            ρₘᵢₙ = max(0, (M′ - M)÷2)
            κ = √T(
                    (factorial((L + M)÷2) * factorial((L - M)÷2))
                    / (factorial((L + M′)÷2) * factorial((L - M′)÷2))
                ) *
                binomial((L + M′)÷2, ρₘᵢₙ) * binomial((L - M′)÷2, (L - M)÷2 - ρₘᵢₙ)
            if (ρₘᵢₙ % 2) != 0
                κ = -κ
            end
            ρₘₐₓ = min((L + M′)÷2, (L - M)÷2)
            N₁ = L + M′ + 2
            N₂ = L - M + 2
            N₃  = M - M′

            total = one(T)
            for P in reverse(2ρₘᵢₙ+2:2:2ρₘₐₓ)
                total *= λ * ((N₁ - P)*(N₂ - P)) / (P*(N₃ + P))
                total += one(T)
            end
            return κ *
                (rₛ ^ (L - (M - M′)÷2 - 2ρₘᵢₙ)) *
                (rₐ ^ ((M - M′)÷2 + 2ρₘᵢₙ)) *
                cis((M + M′)÷2 * ϕₛ + (M - M′)÷2 * ϕₐ) *
                total

        else # rₛ < rₐ
            λ = -(rₛ/rₐ)^2
            ρₘᵢₙ = max(0, -(M′ + M)÷2)
            κ = √T(
                    (factorial((L + M)÷2) * factorial((L - M)÷2))
                    / (factorial((L + M′)÷2) * factorial((L - M′)÷2))
                ) *
                binomial((L + M′)÷2, (L - M)÷2 - ρₘᵢₙ) * binomial((L - M′)÷2, ρₘᵢₙ)
            if (((L - M)÷2 - ρₘᵢₙ) % 2) != 0
                κ = -κ
            end
            ρₘₐₓ = min((L - M′)÷2, (L - M)÷2)
            N₁ = L - M′ + 2
            N₂ = L - M + 2
            N₃  = M + M′

            total = one(T)
            for P in reverse(2ρₘᵢₙ+2:2:2ρₘₐₓ)
                total *= λ * ((N₁ - P)*(N₂ - P)) / (P*(N₃ + P))
                total += one(T)
            end
            return κ *
                (rₐ ^ (L - (M + M′)÷2 - 2ρₘᵢₙ)) *
                (rₛ ^ ((M + M′)÷2 + 2ρₘᵢₙ)) *
                cis((M - M′)÷2 * ϕₐ + (M + M′)÷2 * ϕₛ) *
                total

        end
    end
end


# ---------------------------------------------------------------------------------------
# Varshalovich, Moskalev & Khersonskii (1988).
# Copied verbatim from docs/literate_input/conventions/comparisons/varshalovich_1988.jl:
# the `Factorial` helper (lines 231-236) and Eq. 4.3.1(2) (lines 243-254).
# ---------------------------------------------------------------------------------------

# Factorials of integers and of integer-valued rationals (half-integer arithmetic produces
# the latter), computed exactly, in the postfix form `(n)❗` used throughout these pages:
struct Factorial end
Base.:*(n::Integer, ::Factorial) = factorial(big(n))
Base.:*(n::Rational, ::Factorial) = factorial(big(Int(n)))
const ❗ = Factorial()

function d(J, M, M′, β::T) where {T<:Real}
    if abs(M) > J || abs(M′) > J
        return zero(T)  # convenient when applying Eqs. 4.8.2(14)-(15) at the edges
    end
    (-1)^Int(J-M′) * √T((J+M)❗ * (J-M)❗ * (J+M′)❗ * (J-M′)❗) *
    sum(
        (-1)^k * cos(β/2)^Int(M+M′+2k) * sin(β/2)^Int(2J-M-M′-2k)
        / T((k)❗ * (J-M-k)❗ * (J-M′-k)❗ * (M+M′+k)❗)
        for k ∈ Int(max(0, -(M+M′))):Int(min(J-M, J-M′));
        init=zero(T)
    )
end

# Eq. 4.3.(1):
function D(J, M, M′, α, β, γ)
    exp(-im * M * α) * d(J, M, M′, β) * exp(-im * M′ * γ)
end


# ---------------------------------------------------------------------------------------
# Varshalovich Tables 4.3-4.12, the explicit half-integer `d` functions for J ≤ 9/2.
# Copied verbatim from docs/literate_input/conventions/comparisons/varshalovich_1988.jl:
# 300-467.  The tables list only the rows M ≥ 1/2, and for each such row only the entries
# not obtainable by symmetry from the ones already given; the transcription returns
# `nothing` for the entries the book omits, so every caller must guard on that.
# ---------------------------------------------------------------------------------------

# 0``.
function d_½_explicit(J::Rational{Int}, M::Rational{Int}, M′::Rational{Int}, β::T) where T
    if denominator(J) != 2 || denominator(M) != 2 || denominator(M′) != 2
        error("Only half-integer J, M, M′ are supported")
    end
    if J < 1//2 || J > 9//2
        error("Only J = 1/2, 3/2, 5/2, 7/2, 9/2 are supported")
    end
    if abs(M) > J || abs(M′) > J
        error("abs(M) and abs(M′) must be ≤ J")
    end
    if M < 0
        r = d_½_explicit(J, -M, -M′, β)
        return r === nothing ? nothing : (-1)^Int(M-M′) * r
    else
        let √ = (x -> √T(x))
            if (J, M, M′) == (1//2, 1//2,-1//2)
                -sin(β/2)
            elseif (J, M, M′) == (1//2, 1//2, 1//2)
                cos(β/2)

            elseif (J, M, M′) == (3//2, 1//2,-3//2)
                √3 * sin(β/2)^2 * cos(β/2)
            elseif (J, M, M′) == (3//2, 1//2,-1//2)
                sin(β/2) * (3 * sin(β/2)^2 - 2)
            elseif (J, M, M′) == (3//2, 1//2, 1//2)
                cos(β/2) * (3 * cos(β/2)^2 - 2)
            elseif (J, M, M′) == (3//2, 1//2, 3//2)
                √3 * sin(β/2) * cos(β/2)^2
            elseif (J, M, M′) == (3//2, 3//2,-3//2)
                -sin(β/2)^3
            elseif (J, M, M′) == (3//2, 3//2,-1//2)
                √3 * sin(β/2)^2 * cos(β/2)
            elseif (J, M, M′) == (3//2, 3//2, 1//2)
                -√3 * sin(β/2) * cos(β/2)^2
            elseif (J, M, M′) == (3//2, 3//2, 3//2)
                cos(β/2)^3

            elseif (J, M, M′) == (5//2, 5//2, 5//2)
                cos(β/2)^5
            elseif (J, M, M′) == (5//2, 5//2, 3//2)
                -√5 * sin(β/2) * cos(β/2)^4
            elseif (J, M, M′) == (5//2, 5//2, 1//2)
                √10 * sin(β/2)^2 * cos(β/2)^3
            elseif (J, M, M′) == (5//2, 5//2,-1//2)
                -√10 * sin(β/2)^3 * cos(β/2)^2
            elseif (J, M, M′) == (5//2, 5//2,-3//2)
                √5 * sin(β/2)^4 * cos(β/2)
            elseif (J, M, M′) == (5//2, 5//2,-5//2)
                -sin(β/2)^5
            elseif (J, M, M′) == (5//2, 3//2, 3//2)
                cos(β/2)^3 * (1 - 5 * sin(β/2)^2)
            elseif (J, M, M′) == (5//2, 3//2, 1//2)
                -√2 * sin(β/2) * cos(β/2)^2 * (2 - 5 * sin(β/2)^2)
            elseif (J, M, M′) == (5//2, 3//2,-1//2)
                -√2 * sin(β/2)^2 * cos(β/2) * (2 - 5 * cos(β/2)^2)
            elseif (J, M, M′) == (5//2, 3//2,-3//2)
                sin(β/2)^3 * (1 - 5 * cos(β/2)^2)
            elseif (J, M, M′) == (5//2, 1//2, 1//2)
                cos(β/2) * (3 - 12 * cos(β/2)^2 + 10 * cos(β/2)^4)
            elseif (J, M, M′) == (5//2, 1//2,-1//2)
                -sin(β/2) * (3 - 12 * sin(β/2)^2 + 10 * sin(β/2)^4)

            elseif (J, M, M′) == (7//2, 7//2, 7//2)
                cos(β/2)^7
            elseif (J, M, M′) == (7//2, 7//2, 5//2)
                -√7 * cos(β/2)^6 * sin(β/2)
            elseif (J, M, M′) == (7//2, 7//2, 3//2)
                √21 * cos(β/2)^5 * sin(β/2)^2
            elseif (J, M, M′) == (7//2, 7//2, 1//2)
                -√35 * cos(β/2)^4 * sin(β/2)^3
            elseif (J, M, M′) == (7//2, 7//2,-1//2)
                √35 * cos(β/2)^3 * sin(β/2)^4
            elseif (J, M, M′) == (7//2, 7//2,-3//2)
                -√21 * cos(β/2)^2 * sin(β/2)^5
            elseif (J, M, M′) == (7//2, 7//2,-5//2)
                √7 * cos(β/2) * sin(β/2)^6
            elseif (J, M, M′) == (7//2, 7//2,-7//2)
                -sin(β/2)^7
            elseif (J, M, M′) == (7//2, 5//2, 5//2)
                cos(β/2)^5 * (1 - 7 * sin(β/2)^2)
            elseif (J, M, M′) == (7//2, 5//2, 3//2)
                -√3 * cos(β/2)^4 * sin(β/2) * (2 - 7 * sin(β/2)^2)
            elseif (J, M, M′) == (7//2, 5//2, 1//2)
                √5 * cos(β/2)^3 * sin(β/2)^2 * (3 - 7 * sin(β/2)^2)
            elseif (J, M, M′) == (7//2, 5//2,-1//2)
                √5 * cos(β/2)^2 * sin(β/2)^3 * (3 - 7 * cos(β/2)^2)
            elseif (J, M, M′) == (7//2, 5//2,-3//2)
                -√3 * cos(β/2) * sin(β/2)^4 * (2 - 7 * cos(β/2)^2)
            elseif (J, M, M′) == (7//2, 5//2,-5//2)
                sin(β/2)^5 * (1 - 7 * cos(β/2)^2)
            elseif (J, M, M′) == (7//2, 3//2, 3//2)
                cos(β/2)^3 * (10 - 30 * cos(β/2)^2 + 21 * cos(β/2)^4)
            elseif (J, M, M′) == (7//2, 3//2, 1//2)
                -√15 * cos(β/2)^2 * sin(β/2) * (2 - 8 * cos(β/2)^2 + 7 * cos(β/2)^4)
            elseif (J, M, M′) == (7//2, 3//2,-1//2)
                √15 * cos(β/2) * sin(β/2)^2 * (2 - 8 * sin(β/2)^2 + 7 * sin(β/2)^4)
            elseif (J, M, M′) == (7//2, 3//2,-3//2)
                -sin(β/2)^3 * (10 - 30 * sin(β/2)^2 + 21 * sin(β/2)^4)
            elseif (J, M, M′) == (7//2, 1//2, 1//2)
                -cos(β/2) * (4 - 30 * cos(β/2)^2 + 60 * cos(β/2)^4 - 35 * cos(β/2)^6)
            elseif (J, M, M′) == (7//2, 1//2,-1//2)
                -sin(β/2) * (4 - 30 * sin(β/2)^2 + 60 * sin(β/2)^4 - 35 * sin(β/2)^6)

            elseif (J, M, M′) == (9//2, 9//2, 9//2)
                cos(β/2)^9
            elseif (J, M, M′) == (9//2, 9//2, 7//2)
                -3 * cos(β/2)^8 * sin(β/2)
            elseif (J, M, M′) == (9//2, 9//2, 5//2)
                6 * cos(β/2)^7 * sin(β/2)^2
            elseif (J, M, M′) == (9//2, 9//2, 3//2)
                -2 * √21 * cos(β/2)^6 * sin(β/2)^3
            elseif (J, M, M′) == (9//2, 9//2, 1//2)
                3 * √14 * cos(β/2)^5 * sin(β/2)^4
            elseif (J, M, M′) == (9//2, 9//2,-1//2)
                -3 * √14 * cos(β/2)^4 * sin(β/2)^5
            elseif (J, M, M′) == (9//2, 9//2,-3//2)
                2 * √21 * cos(β/2)^3 * sin(β/2)^6
            elseif (J, M, M′) == (9//2, 9//2,-5//2)
                -6 * cos(β/2)^2 * sin(β/2)^7
            elseif (J, M, M′) == (9//2, 9//2,-7//2)
                3 * cos(β/2) * sin(β/2)^8
            elseif (J, M, M′) == (9//2, 9//2,-9//2)
                -sin(β/2)^9
            elseif (J, M, M′) == (9//2, 7//2, 7//2)
                cos(β/2)^7 * (1 - 9 * sin(β/2)^2)
            elseif (J, M, M′) == (9//2, 7//2, 5//2)
                -2 * cos(β/2)^6 * sin(β/2) * (2 - 9 * sin(β/2)^2)
            elseif (J, M, M′) == (9//2, 7//2, 3//2)
                2 * √21 * cos(β/2)^5 * sin(β/2)^2 * (1 - 3 * sin(β/2)^2)
            elseif (J, M, M′) == (9//2, 7//2, 1//2)
                -√14 * cos(β/2)^4 * sin(β/2)^3 * (4 - 9 * sin(β/2)^2)
            elseif (J, M, M′) == (9//2, 7//2,-1//2)
                -√14 * cos(β/2)^3 * sin(β/2)^4 * (4 - 9 * cos(β/2)^2)
            elseif (J, M, M′) == (9//2, 7//2,-3//2)
                2 * √21 * cos(β/2)^2 * sin(β/2)^5 * (1 - 3 * cos(β/2)^2)
            elseif (J, M, M′) == (9//2, 7//2,-5//2)
                -2 * cos(β/2) * sin(β/2)^6 * (2 - 9 * cos(β/2)^2)
            elseif (J, M, M′) == (9//2, 7//2,-7//2)
                sin(β/2)^7 * (1 - 9 * cos(β/2)^2)
            elseif (J, M, M′) == (9//2, 5//2, 5//2)
                cos(β/2)^5 * (21 - 56 * cos(β/2)^2 + 36 * cos(β/2)^4)
            elseif (J, M, M′) == (9//2, 5//2, 3//2)
                -√21 * cos(β/2)^4 * sin(β/2) * (5 - 16 * cos(β/2)^2 + 12 * cos(β/2)^4)
            elseif (J, M, M′) == (9//2, 5//2, 1//2)
                √14 * cos(β/2)^3 * sin(β/2)^2 * (5 - 20 * cos(β/2)^2 + 18 * cos(β/2)^4)
            elseif (J, M, M′) == (9//2, 5//2,-1//2)
                -√14 * cos(β/2)^2 * sin(β/2)^3 * (5 - 20 * sin(β/2)^2 + 18 * sin(β/2)^4)
            elseif (J, M, M′) == (9//2, 5//2,-3//2)
                √21 * cos(β/2) * sin(β/2)^4 * (5 - 16 * sin(β/2)^2 + 12 * sin(β/2)^4)
            elseif (J, M, M′) == (9//2, 5//2,-5//2)
                -sin(β/2)^5 * (21 - 56 * sin(β/2)^2 + 36 * sin(β/2)^4)
            elseif (J, M, M′) == (9//2, 3//2, 3//2)
                -cos(β/2)^3 * (20 - 105 * cos(β/2)^2 + 168 * cos(β/2)^4 - 84 * cos(β/2)^6)
            elseif (J, M, M′) == (9//2, 3//2, 1//2)
                √6 * cos(β/2)^2 * sin(β/2) * (5 - 35 * cos(β/2)^2 + 70 * cos(β/2)^4 - 42 * cos(β/2)^6)
            elseif (J, M, M′) == (9//2, 3//2,-1//2)
                √6 * cos(β/2) * sin(β/2)^2 * (5 - 35 * sin(β/2)^2 + 70 * sin(β/2)^4 - 42 * sin(β/2)^6)
            elseif (J, M, M′) == (9//2, 3//2,-3//2)
                -sin(β/2)^3 * (20 - 105 * sin(β/2)^2 + 168 * sin(β/2)^4 - 84 * sin(β/2)^6)
            elseif (J, M, M′) == (9//2, 1//2, 1//2)
                cos(β/2) * (5 - 60 * cos(β/2)^2 + 210 * cos(β/2)^4 - 280 * cos(β/2)^6 + 126 * cos(β/2)^8)
            elseif (J, M, M′) == (9//2, 1//2,-1//2)
                -sin(β/2) * (5 - 60 * sin(β/2)^2 + 210 * sin(β/2)^4 - 280 * sin(β/2)^6 + 126 * sin(β/2)^8)
            end
        end
    end
end


# ---------------------------------------------------------------------------------------
# The oracle proper: the two references above, expressed in *this package's* conventions.
# ---------------------------------------------------------------------------------------

"""
    d_oracle(J, m′, m, β)

`dᴶ_{m′m}(β)` in this package's convention, evaluated from Varshalovich Eq. 4.3.1(2).
Varshalovich's index order is ours, so this is just his closed form.  The factorials are
computed exactly (`BigInt`) and only the final combination is in the precision of `β`.
"""
d_oracle(J, m′, m, β::T) where {T<:Real} = d(J, m′, m, β)

"""
    D_oracle(R, J, m′, m)

`𝔇ᴶ_{m′m}(R)` in this package's convention, evaluated from the quaternionic form of Boyle
(2016).  The 2016 convention is the complex conjugate of the present one.
"""
D_oracle(R::Rotor, J, m′, m) = conj(WignerDElement(R, J, m′, m))

"""
    d_table(J, M, M′, β)

The entry of Varshalovich's Tables 4.3-4.12 for `dᴶ_{MM′}(β)`, or `nothing` where the book
does not print it.  Only `1//2 ≤ J ≤ 9//2` is transcribed.
"""
d_table(J, M, M′, β) = d_½_explicit(J, M, M′, β)

"Whole `dᴶ` block from the Varshalovich closed form, as a plain `Matrix` indexed `1:2J+1`."
d_oracle_block(J, β::T) where {T<:Real} = T[d_oracle(J, m′, m, β) for m′ ∈ -J:J, m ∈ -J:J]

"Whole `𝔇ᴶ` block from the Boyle (2016) reference, as a plain `Matrix` indexed `1:2J+1`."
D_oracle_block(R::Rotor{T}, J) where {T} =
    Complex{T}[D_oracle(R, J, m′, m) for m′ ∈ -J:J, m ∈ -J:J]

"""
    maxabsdiff(block, oracle_block)

Largest absolute difference between a half-integer container (or any object supporting
`collect`) and a plain matrix of oracle values in the same `(m′, m)` order.  Accumulating
this and asserting once keeps a failing item from printing thousands of separate failures.
"""
function maxabsdiff(block, oracle_block)
    A = collect(block)
    size(A) == size(oracle_block) || return Inf
    maximum(abs, A .- oracle_block; init=0.0)
end


# ---------------------------------------------------------------------------------------
# Shared, reproducible inputs, so that every item below samples the same angles and rotors.
# ---------------------------------------------------------------------------------------

"""
Nine values of β: both poles, points a whisker away from each of them, and five generic
interior points.  The poles matter because the half-angle pair `(cos(β/2), sin(β/2))` that
seeds the half-integer recurrence degenerates there.
"""
const βvalues = [0.0, 1.0e-3, 0.3, 1.1, 1.5, 2.0, 2.9, π - 1.0e-9, Float64(π)]

"""
    rotors(n=5; seed=1234)

`n + 5` rotors: the identity, the three basis rotors, one rotor with a negative scalar part
(so the double cover is exercised), and `n` random ones drawn from a private RNG, so the
list does not depend on how much of the global RNG stream earlier code has consumed.
"""
function rotors(n::Int=5; seed::Int=1234)
    rng = Random.Xoshiro(seed)
    [
        Rotor{Float64}(1.0);
        Rotor{Float64}(0.0, 1.0, 0.0, 0.0);
        Rotor{Float64}(0.0, 0.0, 1.0, 0.0);
        Rotor{Float64}(0.0, 0.0, 0.0, 1.0);
        -Rotor{Float64}(0.6, 0.8, 0.0, 0.0);
        [randn(rng, Rotor{Float64}) for _ ∈ 1:n]
    ]
end

end  # HalfIntegerOracle
