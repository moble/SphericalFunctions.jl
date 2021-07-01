"""
    Ysize(ℓₘᵢₙ, ℓₘₐₓ)

Compute total size of array of mode weights

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
Yrange : Array of (ℓ, m) indices corresponding to this array
Yindex : Index of a particular element of the mode weight array

Notes
-----
This assumes that the modes are arranged (with fixed s value) as

    [
        Y(s, ℓ, m)
        for ℓ in range(ℓₘᵢₙ, ℓₘₐₓ+1)
        for m in range(-ℓ, ℓ+1)
    ]

"""
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
    for ℓmin in ℓmin_range
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
        for ℓ in range(ℓₘᵢₙ, ℓₘₐₓ+1)
        for m in range(-ℓ, ℓ+1)
    ]

"""
function Yrange(ℓₘᵢₙ, ℓₘₐₓ)
    r = zeros(typeof(ℓₘᵢₙ), (Ysize(ℓₘᵢₙ, ℓₘₐₓ), 2))
    i = 1
    for ℓ in ℓₘᵢₙ:ℓₘₐₓ
        for m in -ℓ:ℓ
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
        for ℓ in range(ℓₘᵢₙ, ℓₘₐₓ+1)
        for m in range(-ℓ, ℓ+1)
    ]

"""
function Yindex(ℓ, m, ℓₘᵢₙ=0)
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


"""Construct (theta, phi) grid

This grid is in the order expected by spinsfast

Parameters
----------
n_theta : int
    Number of points in the theta direction
n_phi : int
    Number of points in the phi direction

Returns
-------
theta_phi_grid : ndarray
    Array of pairs of floats giving the respective [theta, phi] pairs.  The
    shape of this array is (n_theta, n_phi, 2).

Notes
-----
The array looks like

    [
        [θ, ϕ]
        for ϕ ∈ [0, 2π)
        for θ ∈ [0, π]
    ]

(note the open and closed endpoints, respectively), where ϕ and θ are uniformly
sampled in their respective ranges.

"""
function theta_phi(nθ, nϕ)
    [[θ, ϕ][i] for θ in range(0, π, length=nθ), ϕ in range(0, 2π, length=nϕ+1)[begin:end-1], i in 1:2]
end
