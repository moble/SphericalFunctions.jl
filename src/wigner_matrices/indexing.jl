"""
    WignerHsize(ℓₘₐₓ, m′ₘₐₓ=ℓₘₐₓ)

Total size of array of wedges of width m′ₘₐₓ up to ℓₘₐₓ.  If m′ₘₐₓ is not
given, it defaults to ℓₘₐₓ.

See also [`WignerHrange`](@ref) and [`WignerHindex`](@ref).

# Notes
Here, it is assumed that only data with m≥|m′| are stored, and only
corresponding values are passed.  We also assume |m|≤ℓ and |m′|≤ℓ.  Neither
of these are checked.  The wedge array that this function indexes is ordered as

    [
        H(ℓ, m′, m) for ℓ in 0:ℓₘₐₓ
        for m′ in -min(ℓ, m′ₘₐₓ):min(ℓ, m′ₘₐₓ)
        for m in abs(m′):ℓ
    ]

"""
function WignerHsize(ℓₘₐₓ)
    if ℓₘₐₓ < 0
        return 0
    else
        return (ℓₘₐₓ+1) * (ℓₘₐₓ+2) * (2*ℓₘₐₓ+3) ÷ 6
    end
end
function WignerHsize(ℓₘₐₓ, m′ₘₐₓ)
    if ℓₘₐₓ < 0
        return 0
    elseif m′ₘₐₓ >= ℓₘₐₓ
        return ((ℓₘₐₓ+1) * (ℓₘₐₓ+2) * (2*ℓₘₐₓ+3)) ÷ 6
    else
        return (
            (ℓₘₐₓ+1) * (ℓₘₐₓ+2) * (2*ℓₘₐₓ+3)
            - 2*(ℓₘₐₓ-m′ₘₐₓ)*(ℓₘₐₓ-m′ₘₐₓ+1)*(ℓₘₐₓ-m′ₘₐₓ+2)
        ) ÷ 6
    end
end


"""
    WignerHrange(ℓₘₐₓ, m′ₘₐₓ=ℓₘₐₓ)

Create an array of (ℓ, m', m) indices as in H array

See also [`WignerHsize`](@ref) and [`WignerHindex`](@ref)

# Notes

Here, it is assumed that only data with m≥|m'| are stored, and only
corresponding values are passed.  We also assume |m|≤ℓ and |m'|≤ℓ.  Neither of
these are checked.  The wedge array that this function indexes is ordered as

    [
        H(ℓ, m′, m) for ℓ in range(ℓₘₐₓ+1)
        for m′ in range(-min(ℓ, m′ₘₐₓ), min(ℓ, m′ₘₐₓ)+1)
        for m in range(abs(m′), ℓ+1)
    ]

"""
function WignerHrange(ℓₘₐₓ, m′ₘₐₓ=ℓₘₐₓ)
    r = zeros(typeof(m′ₘₐₓ), (WignerHsize(ℓₘₐₓ, m′ₘₐₓ), 3))
    i = 1
    for ℓ in 0:ℓₘₐₓ
        for m′ in -min(ℓ, m′ₘₐₓ):min(ℓ, m′ₘₐₓ)
            for m in abs(m′):ℓ
                r[i, 1] = ℓ
                r[i, 2] = m′
                r[i, 3] = m
                i += 1
            end
        end
    end
    return r
end


"""Helper function for `WignerHindex`"""
function _WignerHindex(ℓ, m′, m, m′ₘₐₓ)
    m′ₘₐₓ = min(m′ₘₐₓ, ℓ)
    i = WignerHsize(ℓ-1, m′ₘₐₓ)  # total size of everything with smaller ℓ
    if m′<1
        i += (m′ₘₐₓ + m′) * (2*ℓ - m′ₘₐₓ + m′ + 1) ÷ 2  # size of wedge to the left of m'
    else
        i += (m′ₘₐₓ + 1) * (2*ℓ - m′ₘₐₓ + 2) ÷ 2  # size of entire left half of wedge
        i += (m′ - 1) * (2*ℓ - m′ + 2) ÷ 2  # size of right half of wedge to the left of m'
    end
    i += m - abs(m′)  # size of column in wedge between m and |m'|
    return i + 1
end


"""
    WignerHindex(ℓ, m′, m, m′ₘₐₓ)

Index to "wedge" arrays

See also [`WignerHsize`](@ref) and [`WignerHrange`](@ref).

# Notes

Here, it is assumed that only data with m≥|m'| are stored, and only
corresponding values are passed.  We also assume |m|≤ℓ and |m'|≤ℓ.  Neither
of these are checked.  The wedge array that this function indexes is ordered as

    [
        H(ℓ, m′, m) for ℓ in range(ℓₘₐₓ+1)
        for m′ in range(-min(ℓ, m′ₘₐₓ), min(ℓ, m′ₘₐₓ)+1)
        for m in range(abs(m′), ℓ+1)
    ]

"""
function WignerHindex(ℓ, m′, m, m′ₘₐₓ=ℓ)
    if ℓ == 0
        return 1
    end
    m′max = min(m′ₘₐₓ, ℓ)
    if m < -m′
        if m < m′
            return _WignerHindex(ℓ, -m′, -m, m′max)
        else
            return _WignerHindex(ℓ, -m, -m′, m′max)
        end
    else
        if m < m′
            return _WignerHindex(ℓ, m, m′, m′max)
        else
            return _WignerHindex(ℓ, m′, m, m′max)
        end
    end
end


"""
    WignerDsize(ℓₘₐₓ, m′ₘₐₓ=ℓₘₐₓ)

Compute total size of Wigner 𝔇 matrix

See also [`WignerDrange`](@ref) and [`WignerDindex`](@ref).

# Notes

This assumes that the Wigner 𝔇 matrix is arranged as

    [
        𝔇(ℓ, m′, m)
        for ℓ in ℓₘᵢₙ:ℓₘₐₓ
        for m′ in -min(ℓ, m′ₘₐₓ):min(ℓ, m′ₘₐₓ)
        for m in -ℓ:ℓ
    ]

"""
function WignerDsize(ℓₘₐₓ)
    (ℓₘₐₓ * (ℓₘₐₓ * (4 * ℓₘₐₓ + 12) + 11) + 3) ÷ 3
end
function WignerDsize(ℓₘₐₓ, m′ₘₐₓ)
    # from sympy import symbols, summation, horner
    # from sympy.printing.pycode import pycode
    # ℓ,m′,m,ℓₘᵢₙ,ℓₘₐₓ,m′ₘₐₓ = symbols('ℓ,m′,m,ℓₘᵢₙ,ℓₘₐₓ,m′ₘₐₓ', integer=True)
    #
    # def nice(expr)
    #     return horner(expr.expand().simplify(), (m′ₘₐₓ, ℓₘᵢₙ, ℓₘₐₓ))
    #
    # # Assuming ℓₘᵢₙ <= ℓₘₐₓ <= m′ₘₐₓ:
    # WignerDsize_ℓmin_ℓmax_m′max = horner(
    #     (
    #         # sum over all m′=0 elements
    #         summation(summation(summation(1, (m, -ℓ, ℓ)), (m′, 0, 0)), (ℓ, ℓₘᵢₙ, ℓₘₐₓ))
    #         # sum over all ℓₘᵢₙ <= ℓ <= ℓₘₐₓ elements
    #         + summation(summation(summation(2, (m, -ℓ, ℓ)), (m′, 1, ℓ)), (ℓ, ℓₘᵢₙ, ℓₘₐₓ))
    #     ).expand().simplify(),
    #     (ℓₘᵢₙ, ℓₘₐₓ)
    # )
    # print(f"({pycode(nice(3*WignerDsize_ℓmin_ℓmax_m′max.subs(ℓₘₐₓ, ℓ-1)))}) ÷ 3")
    #
    # # Assuming ℓₘᵢₙ <= m′ₘₐₓ <= ℓₘₐₓ:
    # WignerDsize_ℓmin_m′max_ℓmax = horner(
    #     (
    #         # sum over all m′=0 elements
    #         summation(summation(summation(1, (m, -ℓ, ℓ)), (m′, 0, 0)), (ℓ, ℓₘᵢₙ, ℓₘₐₓ))
    #         # sum over all ℓ <= m′ₘₐₓ elements
    #         + summation(summation(summation(2, (m, -ℓ, ℓ)), (m′, 1, ℓ)), (ℓ, ℓₘᵢₙ, m′ₘₐₓ))
    #         # sum over all ℓ >= m′ₘₐₓ elements
    #         + summation(summation(summation(2, (m, -ℓ, ℓ)), (m′, 1, m′ₘₐₓ)), (ℓ, m′ₘₐₓ+1, ℓₘₐₓ))
    #     ).expand().simplify(),
    #     (m′ₘₐₓ, ℓₘᵢₙ, ℓₘₐₓ)
    # )
    # print(f"({pycode(nice(3*WignerDsize_ℓmin_m′max_ℓmax.subs(ℓₘₐₓ, ℓ-1)))}) ÷ 3")
    #
    # # Assuming m′ₘₐₓ <= ℓₘᵢₙ <= ℓₘₐₓ:
    # WignerDsize_m′max_ℓmin_ℓmax = horner(
    #     (
    #         # sum over all m′=0 elements
    #         summation(summation(summation(1, (m, -ℓ, ℓ)), (m′, 0, 0)), (ℓ, ℓₘᵢₙ, ℓₘₐₓ))
    #         # sum over all remaining |m′| <= m′ₘₐₓ elements
    #         + summation(summation(summation(2, (m, -ℓ, ℓ)), (m′, 1, m′ₘₐₓ)), (ℓ, ℓₘᵢₙ, ℓₘₐₓ))
    #     ).expand().simplify(),
    #     (m′ₘₐₓ, ℓₘᵢₙ, ℓₘₐₓ)
    # )
    # print(f"{pycode(nice(WignerDsize_m′max_ℓmin_ℓmax.subs(ℓₘₐₓ, ℓ-1)).factor())}")
    if m′ₘₐₓ >= ℓₘₐₓ
        return (ℓₘₐₓ * (ℓₘₐₓ * (4 * ℓₘₐₓ + 12) + 11) + 3) ÷ 3
    elseif m′ₘₐₓ > 0
        return (
            3 * ℓₘₐₓ * (ℓₘₐₓ + 2)
            + m′ₘₐₓ * (
                3 * ℓₘₐₓ * (2 * ℓₘₐₓ + 4)
                + m′ₘₐₓ * (-2 * m′ₘₐₓ - 3) + 5
            )
            + 3
        ) ÷ 3
    else
        return (ℓₘₐₓ * (ℓₘₐₓ + 2)) * (1 + 2 * m′ₘₐₓ) + 2 * m′ₘₐₓ + 1
    end
end


"""
    WignerDrange(ℓₘₐₓ, m′ₘₐₓ=ℓₘₐₓ)

Create an array of (ℓ, m', m) indices as in 𝔇 array

See also [`WignerDsize`](@ref) and [`WignerDindex`](@ref).

# Notes

This assumes that the Wigner 𝔇 matrix is arranged as

    [
        𝔇(ℓ, m′, m)
        for ℓ in range(ℓₘᵢₙ, ℓₘₐₓ+1)
        for m′ in range(-min(ℓ, m′ₘₐₓ), min(ℓ, m′ₘₐₓ)+1)
        for m in range(-ℓ, ℓ+1)
    ]

"""
function WignerDrange(ℓₘₐₓ, m′ₘₐₓ=ℓₘₐₓ)
    r = zeros(typeof(ℓₘₐₓ), (WignerDsize(ℓₘₐₓ, m′ₘₐₓ), 3))
    i = 1
    for ℓ in 0:ℓₘₐₓ
        for m′ in -min(ℓ, m′ₘₐₓ):min(ℓ, m′ₘₐₓ)
            for m in -ℓ:ℓ
                r[i, 1] = ℓ
                r[i, 2] = m′
                r[i, 3] = m
                i += 1
            end
        end
    end
    return r
end


"""
    WignerDindex(ℓ, m′, m, m′ₘₐₓ=ℓ)

Compute index into Wigner 𝔇 matrix

See also [`WignerDrange`](@ref) and [`WignerDsize`](@ref).

# Notes

This assumes that the Wigner 𝔇 matrix is arranged as

    [
        𝔇(ℓ, m′, m)
        for ℓ in range(ℓₘᵢₙ, ℓₘₐₓ+1)
        for m′ in range(-min(ℓ, m′ₘₐₓ), min(ℓ, m′ₘₐₓ)+1)
        for m in range(-ℓ, ℓ+1)
    ]

"""
function WignerDindex(ℓ, m′, m, m′ₘₐₓ=ℓ)
    WignerDsize(ℓ-1, m′ₘₐₓ) + (m′ + min(m′ₘₐₓ, ℓ)) * (2 * ℓ + 1) + m + ℓ + 1
end
