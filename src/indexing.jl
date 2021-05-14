"""
    WignerHsize(m′ₘₐₓ, ℓₘₐₓ=-2)

Total size of array of wedges of width m′ₘₐₓ up to ℓₘₐₓ

Parameters
----------
ℓₘₐₓ : int
m′ₘₐₓ : int, optional
    If nothing, it is assumed to be at least ℓ

See Also
--------
WignerHrange : Array of (ℓ, m', m) indices corresponding to this wedge
WignerHindex : Index inside these wedges

Notes
-----
Here, it is assumed that only data with m≥|m'| are stored, and only
corresponding values are passed.  We also assume |m|≤ℓ and |m'|≤ℓ.  Neither
of these are checked.  The wedge array that this function indexes is ordered as

    [
        H(ℓ, m′, m) for ℓ in range(ℓₘₐₓ+1)
        for m′ in range(-min(ℓ, m′ₘₐₓ), min(ℓ, m′ₘₐₓ)+1)
        for m in range(abs(m′), ℓ+1)
    ]

"""
function WignerHsize(m′ₘₐₓ, ℓₘₐₓ=-2)
    if ℓₘₐₓ == -2
        ℓₘₐₓ = m′ₘₐₓ
    elseif ℓₘₐₓ < 0
        return 0
    end
    if m′ₘₐₓ === nothing || m′ₘₐₓ >= ℓₘₐₓ
        return (ℓₘₐₓ+1) * (ℓₘₐₓ+2) * (2*ℓₘₐₓ+3) ÷ 6
    else
        return ((ℓₘₐₓ+1) * (ℓₘₐₓ+2) * (2*ℓₘₐₓ+3) - 2*(ℓₘₐₓ-m′ₘₐₓ)*(ℓₘₐₓ-m′ₘₐₓ+1)*(ℓₘₐₓ-m′ₘₐₓ+2)) ÷ 6
    end
end


"""
    WignerHrange(m′ₘₐₓ, ℓₘₐₓ=-1)

Create an array of (ℓ, m', m) indices as in H array

Parameters
----------
ℓₘₐₓ : int
m′ₘₐₓ : int, optional
    If nothing, it is assumed to be at least ℓ

See Also
--------
WignerHsize : Total size of wedge array
WignerHindex : Index inside these wedges

Notes
-----
Here, it is assumed that only data with m≥|m'| are stored, and only
corresponding values are passed.  We also assume |m|≤ℓ and |m'|≤ℓ.  Neither of
these are checked.  The wedge array that this function indexes is ordered as

    [
        H(ℓ, m′, m) for ℓ in range(ℓₘₐₓ+1)
        for m′ in range(-min(ℓ, m′ₘₐₓ), min(ℓ, m′ₘₐₓ)+1)
        for m in range(abs(m′), ℓ+1)
    ]

"""
function WignerHrange(m′ₘₐₓ, ℓₘₐₓ=-1)
    if ℓₘₐₓ < 0
        ℓₘₐₓ = m′ₘₐₓ
    end
    r = zeros(typeof(m′ₘₐₓ), (WignerHsize(m′ₘₐₓ, ℓₘₐₓ), 3))
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


function _WignerHindex(ℓ, m′, m, m′ₘₐₓ)
    """Helper function for `WignerHindex`"""
    m′ₘₐₓ = min(m′ₘₐₓ, ℓ)
    i = WignerHsize(m′ₘₐₓ, ℓ-1)  # total size of everything with smaller ℓ
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

Parameters
----------
ℓ : int
m′ : int
m : int
m′ₘₐₓ : int, optional
    If nothing, it is assumed to be at least ℓ

See Also
--------
WignerHsize : Total size of wedge array
WignerHrange : Array of (ℓ, m', m) indices corresponding to this wedge

Notes
-----
Here, it is assumed that only data with m≥|m'| are stored, and only
corresponding values are passed.  We also assume |m|≤ℓ and |m'|≤ℓ.  Neither
of these are checked.  The wedge array that this function indexes is ordered as

    [
        H(ℓ, m′, m) for ℓ in range(ℓₘₐₓ+1)
        for m′ in range(-min(ℓ, m′ₘₐₓ), min(ℓ, m′ₘₐₓ)+1)
        for m in range(abs(m′), ℓ+1)
    ]

"""
function WignerHindex(ℓ, m′, m, m′ₘₐₓ=nothing)
    if ℓ == 0
        return 1
    end
    m′max = ℓ
    if m′ₘₐₓ !== nothing
        m′max = min(m′ₘₐₓ, m′max)
    end
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


# function WignerDsize(ℓₘᵢₙ, m′ₘₐₓ, ℓₘₐₓ=-1)
#     """Compute total size of Wigner 𝔇 matrix

#     Parameters
#     ----------
#     ℓₘᵢₙ : int
#         Integer satisfying 0 <= ℓₘᵢₙ <= ℓₘₐₓ
#     m′ₘₐₓ : int, optional
#         Integer satisfying 0 <= m′ₘₐₓ.  Defaults to ℓₘₐₓ.
#     ℓₘₐₓ : int
#         Integer satisfying 0 <= ℓₘᵢₙ <= ℓₘₐₓ

#     Returns
#     -------
#     i : int
#         Total size of Wigner 𝔇 matrix arranged as described below

#     See Also
#     --------
#     WignerDrange : Array of (ℓ, m', m) indices corresponding to the 𝔇 matrix
#     WignerDindex : Index of a particular element of the 𝔇 matrix

#     Notes
#     -----
#     This assumes that the Wigner 𝔇 matrix is arranged as

#         [
#             𝔇(ℓ, m′, m)
#             for ℓ in range(ℓₘᵢₙ, ℓₘₐₓ+1)
#             for m′ in range(-min(ℓ, m′ₘₐₓ), min(ℓ, m′ₘₐₓ)+1)
#             for m in range(-ℓ, ℓ+1)
#         ]

#     """
#     # from sympy import symbols, summation, horner
#     # from sympy.printing.pycode import pycode
#     # ℓ,m′,m,ℓₘᵢₙ,ℓₘₐₓ,m′ₘₐₓ = symbols('ℓ,m′,m,ℓₘᵢₙ,ℓₘₐₓ,m′ₘₐₓ', integer=True)
#     # 
#     # def nice(expr)
#     #     return horner(expr.expand().simplify(), (m′ₘₐₓ, ℓₘᵢₙ, ℓₘₐₓ))
#     #
#     # # Assuming ℓₘᵢₙ <= ℓₘₐₓ <= m′ₘₐₓ:
#     # WignerDsize_ℓmin_ℓmax_m′max = horner(
#     #     (
#     #         # sum over all m′=0 elements
#     #         summation(summation(summation(1, (m, -ℓ, ℓ)), (m′, 0, 0)), (ℓ, ℓₘᵢₙ, ℓₘₐₓ))
#     #         # sum over all ℓₘᵢₙ <= ℓ <= ℓₘₐₓ elements
#     #         + summation(summation(summation(2, (m, -ℓ, ℓ)), (m′, 1, ℓ)), (ℓ, ℓₘᵢₙ, ℓₘₐₓ))
#     #     ).expand().simplify(),
#     #     (ℓₘᵢₙ, ℓₘₐₓ)
#     # )
#     # print(f"({pycode(nice(3*WignerDsize_ℓmin_ℓmax_m′max.subs(ℓₘₐₓ, ℓ-1)))}) ÷ 3")
#     # 
#     # # Assuming ℓₘᵢₙ <= m′ₘₐₓ <= ℓₘₐₓ:
#     # WignerDsize_ℓmin_m′max_ℓmax = horner(
#     #     (
#     #         # sum over all m′=0 elements
#     #         summation(summation(summation(1, (m, -ℓ, ℓ)), (m′, 0, 0)), (ℓ, ℓₘᵢₙ, ℓₘₐₓ))
#     #         # sum over all ℓ <= m′ₘₐₓ elements
#     #         + summation(summation(summation(2, (m, -ℓ, ℓ)), (m′, 1, ℓ)), (ℓ, ℓₘᵢₙ, m′ₘₐₓ))
#     #         # sum over all ℓ >= m′ₘₐₓ elements
#     #         + summation(summation(summation(2, (m, -ℓ, ℓ)), (m′, 1, m′ₘₐₓ)), (ℓ, m′ₘₐₓ+1, ℓₘₐₓ))
#     #     ).expand().simplify(),
#     #     (m′ₘₐₓ, ℓₘᵢₙ, ℓₘₐₓ)
#     # )
#     # print(f"({pycode(nice(3*WignerDsize_ℓmin_m′max_ℓmax.subs(ℓₘₐₓ, ℓ-1)))}) ÷ 3")
#     #
#     # # Assuming m′ₘₐₓ <= ℓₘᵢₙ <= ℓₘₐₓ:
#     # WignerDsize_m′max_ℓmin_ℓmax = horner(
#     #     (
#     #         # sum over all m′=0 elements
#     #         summation(summation(summation(1, (m, -ℓ, ℓ)), (m′, 0, 0)), (ℓ, ℓₘᵢₙ, ℓₘₐₓ))
#     #         # sum over all remaining |m′| <= m′ₘₐₓ elements
#     #         + summation(summation(summation(2, (m, -ℓ, ℓ)), (m′, 1, m′ₘₐₓ)), (ℓ, ℓₘᵢₙ, ℓₘₐₓ))
#     #     ).expand().simplify(),
#     #     (m′ₘₐₓ, ℓₘᵢₙ, ℓₘₐₓ)
#     # )
#     # print(f"{pycode(nice(WignerDsize_m′max_ℓmin_ℓmax.subs(ℓₘₐₓ, ℓ-1)).factor())}")
#     if ℓₘₐₓ < 0
#         ℓₘₐₓ = m′ₘₐₓ
#     end
#     if m′ₘₐₓ >= ℓₘₐₓ
#         return (
#             ℓₘₐₓ * (ℓₘₐₓ * (4 * ℓₘₐₓ + 12) + 11)
#             + ℓₘᵢₙ * (1 - 4 * ℓₘᵢₙ^2)
#             + 3
#         ) ÷ 3
#     end
#     if m′ₘₐₓ > ℓₘᵢₙ
#         return (
#             3 * ℓₘₐₓ * (ℓₘₐₓ + 2)
#             + ℓₘᵢₙ * (1 - 4 * ℓₘᵢₙ^2)
#             + m′ₘₐₓ * (
#                 3 * ℓₘₐₓ * (2 * ℓₘₐₓ + 4)
#                 + m′ₘₐₓ * (-2 * m′ₘₐₓ - 3) + 5
#             )
#             + 3
#         ) ÷ 3
#     else
#         return (ℓₘₐₓ * (ℓₘₐₓ + 2) - ℓₘᵢₙ^2) * (1 + 2 * m′ₘₐₓ) + 2 * m′ₘₐₓ + 1
#     end
# end


# function WignerDrange(ℓₘᵢₙ, m′ₘₐₓ, ℓₘₐₓ=-1)
#     """Create an array of (ℓ, m', m) indices as in 𝔇 array

#     Parameters
#     ----------
#     ℓₘᵢₙ : int
#         Integer satisfying 0 <= ℓₘᵢₙ <= ℓₘₐₓ
#     m′ₘₐₓ : int, optional
#         Integer satisfying 0 <= m′ₘₐₓ.  Default is ℓₘₐₓ.
#     ℓₘₐₓ : int
#         Integer satisfying 0 <= ℓₘᵢₙ <= ℓₘₐₓ

#     See Also
#     --------
#     WignerDsize : Total size of 𝔇 array
#     WignerDindex : Index inside these wedges

#     Notes
#     -----
#     This assumes that the Wigner 𝔇 matrix is arranged as

#         [
#             𝔇(ℓ, m′, m)
#             for ℓ in range(ℓₘᵢₙ, ℓₘₐₓ+1)
#             for m′ in range(-min(ℓ, m′ₘₐₓ), min(ℓ, m′ₘₐₓ)+1)
#             for m in range(-ℓ, ℓ+1)
#         ]

#     """
#     if ℓₘₐₓ < 0
#         ℓₘₐₓ = m′ₘₐₓ
#     end
#     r = np.zeros((WignerDsize(ℓₘᵢₙ, m′ₘₐₓ, ℓₘₐₓ), 3), dtype=np.int64)
#     i = 0
#     for ℓ in range(ℓₘᵢₙ, ℓₘₐₓ+1)
#         for m′ in range(-min(ℓ, m′ₘₐₓ), min(ℓ, m′ₘₐₓ)+1)
#             for m in range(-ℓ, ℓ+1)
#                 r[i, 0] = ℓ
#                 r[i, 1] = m′
#                 r[i, 2] = m
#                 i += 1
#             end
#         end
#     end
#     return r
# end


# @jit
# def WignerDindex(ℓ, m′, m, ℓₘᵢₙ=0, m′ₘₐₓ=-1)
#     """Compute index into Wigner 𝔇 matrix

#     Parameters
#     ----------
#     ℓ : int
#         Integer satisfying ℓₘᵢₙ <= ℓ <= ℓₘₐₓ
#     m′ : int
#         Integer satisfying -min(ℓₘₐₓ, m′ₘₐₓ) <= m′ <= min(ℓₘₐₓ, m′ₘₐₓ)
#     m : int
#         Integer satisfying -ℓ <= m <= ℓ
#     ℓₘᵢₙ : int, optional
#         Integer satisfying 0 <= ℓₘᵢₙ <= ℓₘₐₓ.  Defaults to 0.
#     m′ₘₐₓ : int, optional
#         Integer satisfying 0 <= m′ₘₐₓ.  Defaults to ℓ.

#     Returns
#     -------
#     i : int
#         Index into Wigner 𝔇 matrix arranged as described below

#     See Also
#     --------
#     WignerDsize : Total size of the 𝔇 matrix
#     WignerDrange : Array of (ℓ, m', m) indices corresponding to the 𝔇 matrix

#     Notes
#     -----
#     This assumes that the Wigner 𝔇 matrix is arranged as

#         [
#             𝔇(ℓ, m′, m)
#             for ℓ in range(ℓₘᵢₙ, ℓₘₐₓ+1)
#             for m′ in range(-min(ℓ, m′ₘₐₓ), min(ℓ, m′ₘₐₓ)+1)
#             for m in range(-ℓ, ℓ+1)
#         ]

#     """
#     if m′ₘₐₓ < 0
#         m′ₘₐₓ = ℓ
#     i = (m′ + min(m′ₘₐₓ, ℓ)) * (2 * ℓ + 1) + m + ℓ
#     if ℓ > ℓₘᵢₙ
#         i += WignerDsize(ℓₘᵢₙ, m′ₘₐₓ, ℓ-1)
#     return i


# @jit
# def Ysize(ℓₘᵢₙ, ℓₘₐₓ)
#     """Compute total size of array of mode weights

#     Parameters
#     ----------
#     ℓₘᵢₙ : int
#         Integer satisfying 0 <= ℓₘᵢₙ <= ℓₘₐₓ
#     ℓₘₐₓ : int
#         Integer satisfying 0 <= ℓₘᵢₙ <= ℓₘₐₓ

#     Returns
#     -------
#     i : int
#         Total size of array of mode weights arranged as described below

#     See Also
#     --------
#     Yrange : Array of (ℓ, m) indices corresponding to this array
#     Yindex : Index of a particular element of the mode weight array

#     Notes
#     -----
#     This assumes that the modes are arranged (with fixed s value) as

#         [
#             Y(s, ℓ, m)
#             for ℓ in range(ℓₘᵢₙ, ℓₘₐₓ+1)
#             for m in range(-ℓ, ℓ+1)
#         ]

#     """
#     # from sympy import symbols, summation, horner
#     # from sympy.printing.pycode import pycode
#     # ℓ,m,ℓₘᵢₙ,ℓₘₐₓ = symbols('ℓ,m,ℓₘᵢₙ,ℓₘₐₓ', integer=True)
#     # 
#     # def nice(expr)
#     #     return horner(expr.expand().simplify(), (m′ₘₐₓ, ℓₘᵢₙ, ℓₘₐₓ))
#     #
#     # Ysize = horner(
#     #     summation(summation(1, (m, -ℓ, ℓ)), (ℓ, ℓₘᵢₙ, ℓₘₐₓ)).expand().simplify(),
#     #     (ℓₘᵢₙ, ℓₘₐₓ)
#     # )
#     return ℓₘₐₓ * (ℓₘₐₓ + 2) - ℓₘᵢₙ^2 + 1


# @jit
# def Yrange(ℓₘᵢₙ, ℓₘₐₓ)
#     """Create an array of (ℓ, m) indices as in Y array

#     Parameters
#     ----------
#     ℓₘᵢₙ : int
#         Integer satisfying 0 <= ℓₘᵢₙ <= ℓₘₐₓ
#     ℓₘₐₓ : int
#         Integer satisfying 0 <= ℓₘᵢₙ <= ℓₘₐₓ

#     Returns
#     -------
#     i : int
#         Total size of array of mode weights arranged as described below

#     See Also
#     --------
#     Ysize : Total size of array of mode weights
#     Yindex : Index of a particular element of the mode weight array

#     Notes
#     -----
#     This assumes that the modes are arranged (with fixed s value) as

#         [
#             Y(s, ℓ, m)
#             for ℓ in range(ℓₘᵢₙ, ℓₘₐₓ+1)
#             for m in range(-ℓ, ℓ+1)
#         ]

#     """
#     r = np.zeros((Ysize(ℓₘᵢₙ, ℓₘₐₓ), 2), dtype=np.int64)
#     i = 0
#     for ℓ in range(ℓₘᵢₙ, ℓₘₐₓ+1)
#         for m in range(-ℓ, ℓ+1)
#             r[i, 0] = ℓ
#             r[i, 1] = m
#             i += 1
#     return r


# @jit
# def Yindex(ℓ, m, ℓₘᵢₙ=0)
#     """Compute index into array of mode weights

#     Parameters
#     ----------
#     ℓ : int
#         Integer satisfying ℓₘᵢₙ <= ℓ <= ℓₘₐₓ
#     m : int
#         Integer satisfying -ℓ <= m <= ℓ
#     ℓₘᵢₙ : int, optional
#         Integer satisfying 0 <= ℓₘᵢₙ.  Defaults to 0.

#     Returns
#     -------
#     i : int
#         Index of a particular element of the mode-weight array as described below

#     See Also
#     --------
#     Ysize : Total size of array of mode weights
#     Yrange : Array of (ℓ, m) indices corresponding to this array

#     Notes
#     -----
#     This assumes that the modes are arranged (with fixed s value) as

#         [
#             Y(s, ℓ, m)
#             for ℓ in range(ℓₘᵢₙ, ℓₘₐₓ+1)
#             for m in range(-ℓ, ℓ+1)
#         ]

#     """
#     # from sympy import symbols, summation, horner
#     # from sympy.printing.pycode import pycode
#     # ℓ,m,m′,ℓₘᵢₙ, = symbols('ℓ,m,m′,ℓₘᵢₙ', integer=True)
#     # 
#     # def nice(expr)
#     #     return horner(expr.expand().simplify(), (ℓₘᵢₙ, ℓ, m))
#     #
#     # Yindex = horner(
#     #     (Ysize.subs(ℓₘₐₓ, ℓ-1) + summation(1, (m′, -ℓ, m)) - 1).expand().simplify(),
#     #     (ℓₘₐₓ, ℓ, m)
#     # )
#     if ℓ > ℓₘᵢₙ
#         return ℓ*(ℓ + 1) - ℓₘᵢₙ^2 + m
#     else
#         return m + ℓ


# def theta_phi(n_theta, n_phi)
#     """Construct (theta, phi) grid

#     This grid is in the order expected by spinsfast

#     Parameters
#     ----------
#     n_theta : int
#         Number of points in the theta direction
#     n_phi : int
#         Number of points in the phi direction

#     Returns
#     -------
#     theta_phi_grid : ndarray
#         Array of pairs of floats giving the respective [theta, phi] pairs.  The
#         shape of this array is (n_theta, n_phi, 2).

#     Notes
#     -----
#     The array looks like

#         [
#             [θ, ϕ]
#             for ϕ ∈ [0, 2π)
#             for θ ∈ [0, π]
#         ]

#     (note the open and closed endpoints, respectively), where ϕ and θ are uniformly
#     sampled in their respective ranges.

#     """
#     return np.array([[[theta, phi]
#                       for phi in np.linspace(0.0, 2*π, num=n_phi, endpoint=False)]
#                      for theta in np.linspace(0.0, π, num=n_theta, endpoint=True)])

