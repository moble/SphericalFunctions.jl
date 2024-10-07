"""
    Ysize(ℓₘₐₓ)

Compute total size of array of mode weights

See Also
--------
Yrange : Array of (ℓ, m) indices corresponding to this array
Yindex : Index of a particular element of the mode weight array

Notes
-----
This assumes that the modes are arranged (with fixed s value) as

    [
        Y(s, ℓ, m)
        for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ
        for m ∈ -ℓ:ℓ
    ]

"""
function Ysize(ℓₘₐₓ)
    return ℓₘₐₓ * (ℓₘₐₓ + 2) + 1
end
function Ysize(ℓₘᵢₙ, ℓₘₐₓ)
    # from sympy import symbols, summation, horner
    # from sympy.printing.pycode import pycode
    # ℓ,m,ℓₘᵢₙ,ℓₘₐₓ = symbols('ℓ,m,ℓₘᵢₙ,ℓₘₐₓ', integer=True)
    #
    # def nice(expr)
    #     return horner(expr.expand().simplify(), (m′ₘₐₓ, ℓₘᵢₙ, ℓₘₐₓ))
    #
    # Ysize = horner(
    #     summation(summation(1, (m, -ℓ, ℓ)), (ℓ, ℓₘᵢₙ, ℓₘₐₓ)).expand().simplify(),
    #     (ℓₘᵢₙ, ℓₘₐₓ)
    # )
    return ℓₘₐₓ * (ℓₘₐₓ + 2) - ℓₘᵢₙ^2 + 1
end


"""
    deduce_limits(ysize, [ℓmin])

Deduce the value of `(ℓmin, ℓmax)` that produces ``Y`` arrays of the given size.

If `ℓmin` is not given, it is assumed to be 0.  If it is set to `nothing`, the
smallest possible value of `ℓmin` will be used.  However, note that this is not
a well-posed problem; multiple combinations of `(ℓmin, ℓmax)` can give rise to
``Y`` arrays of the same size.

See also [`Ysize`](@ref)
"""
function deduce_limits(ysize, ℓmin=0)
    if ℓmin === nothing
        maxℓmax = (ysize-1) ÷ 2
        ℓmin_range = 0:maxℓmax
    else
        ℓmin_range = [ℓmin]
    end
    for ℓmin ∈ ℓmin_range
        ℓmax = round(Int, √(ysize + ℓmin^2) - 1, RoundDown)
        if ℓmax * (ℓmax + 2) - ℓmin^2 + 1 == ysize
            return (ℓmin, ℓmax)
        end
    end
    error("Input size $ysize does not correspond to a possible array of modes with ℓmin=$ℓmin")
end


"""
    Yrange(ℓₘᵢₙ, ℓₘₐₓ)

Create an array of (ℓ, m) indices as in Y array

Parameters
----------
ℓₘᵢₙ : int
    Integer satisfying 0 <= ℓₘᵢₙ <= ℓₘₐₓ
ℓₘₐₓ : int
    Integer satisfying 0 <= ℓₘᵢₙ <= ℓₘₐₓ

Returns
-------
i : int
    Total size of array of mode weights arranged as described below

See Also
--------
Ysize : Total size of array of mode weights
Yindex : Index of a particular element of the mode weight array

Notes
-----
This assumes that the modes are arranged (with fixed s value) as

    [
        Y(s, ℓ, m)
        for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ
        for m ∈ -ℓ:ℓ
    ]

"""
Yrange(ℓₘₐₓ) = Yrange(0, ℓₘₐₓ)
function Yrange(ℓₘᵢₙ, ℓₘₐₓ)
    r = zeros(typeof(ℓₘₐₓ), (Ysize(ℓₘᵢₙ, ℓₘₐₓ), 2))
    i = 1
    for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ
        for m ∈ -ℓ:ℓ
            r[i, 1] = ℓ
            r[i, 2] = m
            i += 1
        end
    end
    return r
end


"""
    Yindex(ℓ, m, ℓₘᵢₙ=0)

Compute index into array of mode weights

Parameters
----------
ℓ : int
    Integer satisfying ℓₘᵢₙ <= ℓ <= ℓₘₐₓ
m : int
    Integer satisfying -ℓ <= m <= ℓ
ℓₘᵢₙ : int, optional
    Integer satisfying 0 <= ℓₘᵢₙ.  Defaults to 0.

Returns
-------
i : int
    Index of a particular element of the mode-weight array as described below

See Also
--------
Ysize : Total size of array of mode weights
Yrange : Array of (ℓ, m) indices corresponding to this array

Notes
-----
This assumes that the modes are arranged (with fixed s value) as

    [
        Y(s, ℓ, m)
        for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ
        for m ∈ -ℓ:ℓ
    ]

"""
@inline function Yindex(ℓ, m, ℓₘᵢₙ=0)
    # from sympy import symbols, summation, horner
    # from sympy.printing.pycode import pycode
    # ℓ,m,m′,ℓₘᵢₙ, = symbols('ℓ,m,m′,ℓₘᵢₙ', integer=True)
    #
    # def nice(expr)
    #     return horner(expr.expand().simplify(), (ℓₘᵢₙ, ℓ, m))
    #
    # Yindex = horner(
    #     (Ysize.subs(ℓₘₐₓ, ℓ-1) + summation(1, (m′, -ℓ, m)) - 1).expand().simplify(),
    #     (ℓₘₐₓ, ℓ, m)
    # )
    return ℓ*(ℓ + 1) - ℓₘᵢₙ^2 + m + 1
end


"""
    theta_phi(Nθ, Nϕ, [T=Float64])

Construct (theta, phi) grid in `spinsfast` order.

Note that this order is different from the one assumed by this package;
use [`phi_theta`](@ref) for the opposite ordering.

"""
function theta_phi(nθ, nϕ, ::Type{T}=Float64) where T
    [
        [θ, ϕ][i]
        for θ ∈ range(0, T(π), length=nθ),
            ϕ ∈ range(0, 2*T(π), length=nϕ+1)[begin:end-1],
            i ∈ 1:2
    ]
end

"""
    phi_theta(Nϕ, Nθ, [T=Float64])

Construct (phi, theta) grid in order expected by this package.

See also [`theta_phi`](@ref) for the opposite ordering.

"""
function phi_theta(nϕ, nθ, ::Type{T}=Float64) where T
    [
        [ϕ, θ][i]
        for ϕ ∈ range(0, 2*T(π), length=nϕ+1)[begin:end-1],
            θ ∈ range(0, T(π), length=nθ),
            i ∈ 1:2
    ]
end


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
        H(ℓ, m′, m) for ℓ ∈ 0:ℓₘₐₓ
        for m′ ∈ -min(ℓ, m′ₘₐₓ):min(ℓ, m′ₘₐₓ)
        for m ∈ abs(m′):ℓ
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
        H(ℓ, m′, m) for ℓ ∈ 0:ℓₘₐₓ
        for m′ ∈ -min(ℓ, m′ₘₐₓ):min(ℓ, m′ₘₐₓ)
        for m ∈ abs(m′):ℓ
    ]

"""
function WignerHrange(ℓₘₐₓ, m′ₘₐₓ=ℓₘₐₓ)
    r = zeros(typeof(m′ₘₐₓ), (WignerHsize(ℓₘₐₓ, m′ₘₐₓ), 3))
    i = 1
    for ℓ ∈ 0:ℓₘₐₓ
        for m′ ∈ -min(ℓ, m′ₘₐₓ):min(ℓ, m′ₘₐₓ)
            for m ∈ abs(m′):ℓ
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
WignerHindex(ℓ, m′, m, m′ₘₐₓ)

Helper function for [`WignerHindex`](@ref), with more constraints.

This function assumes that `m ≥ |m′|`.  The main [`WignerHindex`](@ref) function uses
symmetries of the H array to account for cases that violate this assumption.  (But note that
both that function and this one assume that `|m| ≤ ℓ` and `|m′| ≤ ℓ`.)

"""
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

Index to "wedge" arrays.

See also [`WignerHsize`](@ref) and [`WignerHrange`](@ref).

# Notes

Here, it is assumed that only data with m≥|m'| are stored, and only
corresponding values are passed.  We also assume |m|≤ℓ and |m'|≤ℓ.  Neither
of these are checked.  The wedge array that this function indexes is ordered as

    [
        H(ℓ, m′, m) for ℓ ∈ 0:ℓₘₐₓ
        for m′ ∈ -min(ℓ, m′ₘₐₓ):min(ℓ, m′ₘₐₓ)
        for m ∈ abs(m′):ℓ
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
        for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ
        for m′ ∈ -min(ℓ, m′ₘₐₓ):min(ℓ, m′ₘₐₓ)
        for m ∈ -ℓ:ℓ
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
        for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ
        for m′ ∈ -min(ℓ, m′ₘₐₓ):min(ℓ, m′ₘₐₓ)
        for m ∈ -ℓ:ℓ
    ]

"""
function WignerDrange(ℓₘₐₓ, m′ₘₐₓ=ℓₘₐₓ)
    r = zeros(typeof(ℓₘₐₓ), (WignerDsize(ℓₘₐₓ, m′ₘₐₓ), 3))
    i = 1
    for ℓ ∈ 0:ℓₘₐₓ
        for m′ ∈ -min(ℓ, m′ₘₐₓ):min(ℓ, m′ₘₐₓ)
            for m ∈ -ℓ:ℓ
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
        for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ
        for m′ ∈ -min(ℓ, m′ₘₐₓ):min(ℓ, m′ₘₐₓ)
        for m ∈ -ℓ:ℓ
    ]

"""
function WignerDindex(ℓ, m′, m, m′ₘₐₓ=ℓ)
    WignerDsize(ℓ-1, m′ₘₐₓ) + (m′ + min(m′ₘₐₓ, ℓ)) * (2 * ℓ + 1) + m + ℓ + 1
end
