### TODO:
### 1. Test speeds without caching a, b, d; maybe switch
### 2. Separate ALF computation to a different module
### 3. Compute H directly inside D array
### 4. Test skipping all the complicated indexing tricks; use fancy indexing
### 5. Allow specifying the extent of recursion / iterating over ℓ matrices


### NOTES:
### 1. Caching coefficients provides a ~11% speedup at low ℓ, but that falls off
###    as ℓ increases, reaching ~4% around ℓ=512.
### 2. This is probably a much more significant advantage for ALFs.



export WignerMatrixCalculator, H!, d!, D!, Y!, ϵ, WignerHindex, _WignerHindex

using ..Spherical: complex_powers!
using Quaternionic: AbstractQuaternion, to_euler_phases!

include("indexing.jl")
include("calculator.jl")
include("Hrecursions.jl")


@inline ϵ(m) = (m <= 0 ? 1 : (isodd(m) ? -1 : 1))


"""Return sign of input, with sign(0)=1"""
@inline sign(m) = (m < 0 ? -1 : 1)


"""
    d!(d, wigner, expiβ)

Compute Wigner's d matrix dˡₘₚ,ₘ(β)

# Parameters

* expiβ : array_like
    Values of expi(i*β) on which to evaluate the d matrix.
* out : array_like, optional
    Array into which the d values should be written.  It should be an array of
    floats, with size `self.dsize`.  If not present, the array will be created.
    In either case, the array will also be returned.
* workspace : array_like, optional
    A working array like the one returned by Wigner.new_workspace().  If not
    present, this object's default workspace will be used.  Note that it is not
    safe to use the same workspace on multiple threads.

# Returns

* d : array
    This is a 1-dimensional array of floats; see below.

# See Also

* H : Compute a portion of the H matrix
* D : Compute the full Wigner 𝔇 matrix
* rotate : Avoid computing the full 𝔇 matrix and rotate modes directly
* evaluate : Avoid computing the full 𝔇 matrix and evaluate modes directly

# Notes

This function is the preferred method of computing the d matrix for large ell
values.  In particular, above ell≈32 standard formulas become completely
unusable because of numerical instabilities and overflow.  This function uses
stable recursion methods instead, and should be usable beyond ell≈1000.

The result is returned in a 1-dimensional array ordered as

    [
        d(ell, mp, m, β)
        for ell in range(ell_max+1)
        for mp in range(-min(ℓ, mp_max), min(ℓ, mp_max)+1)
        for m in range(-ell, ell+1)
    ]

"""
function d!(d, w::WignerMatrixCalculator, expiβ::Complex)
    ell_min = ℓₘᵢₙ(w)
    ell_max = ℓₘₐₓ(w)
    mp_max = m′ₘₐₓ(w)
    if mp_max < ell_max
        throw(DomainError("ℓₘₐₓ = $ℓₘₐₓ",
                          "Cannot compute full d matrix up to ℓₘₐₓ with m′ₘₐₓ only $(mp_max)"
        ))
    end

    H!(w, expiβ)

    for ℓ in ell_min:ell_max
        for m′ in -ℓ:ℓ
            for m in -ℓ:ℓ
                i_d = WignerDindex(ℓ, m′, m, ell_min)
                i_H = WignerHindex(ℓ, m′, m, mp_max)
                d[i_d] = ϵ(m′) * ϵ(-m) * w.Hwedge[i_H]
            end
        end
    end

    d
end


function d!(w::WignerMatrixCalculator, β::Real)
    d = zeros(T(w), Wignerdsize(w))
    d!(w, expiβ, out)
end


"""
    D!(𝔇, w, R)
    D!(w, R)

Compute Wigner's 𝔇 matrix

# Parameters

* 𝔇 : array_like, optional
    Array into which the 𝔇 values should be written.  It should be an array of
    complex, with size `self.Dsize`.  If not present, the array will be
    created.  In either case, the array will also be returned.
* workspace : optional
    A working array like the one returned by Wigner.new_workspace().  If not
    present, this object's default workspace will be used.  Note that it is not
    safe to use the same workspace on multiple threads.
* R : Quaternion
    Array to be interpreted as a quaternionic array (thus its final dimension
    must have size 4), representing the rotations on which the 𝔇 matrix will be
    evaluated.

# Returns

* D : array
    This is a 1-dimensional array of complex; see below.

# See Also

* H : Compute a portion of the H matrix
* d : Compute the full Wigner d matrix
* rotate : Avoid computing the full 𝔇 matrix and rotate modes directly
* evaluate : Avoid computing the full 𝔇 matrix and evaluate modes directly

# Notes

This function is the preferred method of computing the 𝔇 matrix for large ell
values.  In particular, above ell≈32 standard formulas become completely
unusable because of numerical instabilities and overflow.  This function uses
stable recursion methods instead, and should be usable beyond ell≈1000.

This function computes 𝔇ˡₘₚ,ₘ(R).  The result is returned in a 1-dimensional
array ordered as

    [
        𝔇(ell, mp, m, R)
        for ell in range(ell_max+1)
        for mp in range(-min(ℓ, mp_max), min(ℓ, mp_max)+1)
        for m in range(-ell, ell+1)
    ]

"""
function D!(𝔇, w::WignerMatrixCalculator, R::AbstractQuaternion)
    ell_min = ℓₘᵢₙ(w)
    ell_max = ℓₘₐₓ(w)
    mp_max = m′ₘₐₓ(w)
    if mp_max < ell_max
        throw(DomainError("ℓₘₐₓ = $(ell_max)",
                          "Cannot compute full d matrix up to ℓₘₐₓ with m′ₘₐₓ only $(mp_max)"
        ))
    end

    to_euler_phases!(w.z, R)
    H!(w, w.z[2])
    complex_powers!(w.zₐpowers, w.z[1])
    complex_powers!(w.zᵧpowers, w.z[3])

    # 𝔇ˡₘₚ,ₘ(R) = dˡₘₚ,ₘ(R) exp[iϕₐ(m-mp)+iϕₛ(m+mp)] = dˡₘₚ,ₘ(R) exp[i(ϕₛ+ϕₐ)m+i(ϕₛ-ϕₐ)mp]
    # exp[iϕₛ] = R̂ₛ = hat(R[0] + 1j * R[3]) = zp
    # exp[iϕₐ] = R̂ₐ = hat(R[2] + 1j * R[1]) = zm.conjugate()
    # exp[i(ϕₛ+ϕₐ)] = zp * zm.conjugate() = z[2] = zᵧ
    # exp[i(ϕₛ-ϕₐ)] = zp * zm = z[0] = zₐ
    for ell in ell_min:ell_max
        for mp in -ell:-1
            i_D = WignerDindex(ell, mp, -ell, ell_min)
            for m in -ell:-1
                i_H = WignerHindex(ell, mp, m, mp_max)
                𝔇[i_D] = ϵ(mp) * ϵ(-m) * w.Hwedge[i_H] * conj(w.zᵧpowers[-m+1]) * conj(w.zₐpowers[-mp+1])
                i_D += 1
            end
            for m in 0:ell
                i_H = WignerHindex(ell, mp, m, mp_max)
                𝔇[i_D] = ϵ(mp) * ϵ(-m) * w.Hwedge[i_H] * w.zᵧpowers[m+1] * conj(w.zₐpowers[-mp+1])
                i_D += 1
            end
        end
        for mp in 0:ell
            i_D = WignerDindex(ell, mp, -ell, ell_min)
            for m in -ell:-1
                i_H = WignerHindex(ell, mp, m, mp_max)
                𝔇[i_D] = ϵ(mp) * ϵ(-m) * w.Hwedge[i_H] * conj(w.zᵧpowers[-m+1]) * w.zₐpowers[mp+1]
                i_D += 1
            end
            for m in 0:ell
                i_H = WignerHindex(ell, mp, m, mp_max)
                𝔇[i_D] = ϵ(mp) * ϵ(-m) * w.Hwedge[i_H] * w.zᵧpowers[m+1] * w.zₐpowers[mp+1]
                i_D += 1
            end
        end
    end
    𝔇
end


function D!(w::WignerMatrixCalculator, R::AbstractQuaternion)
    𝔇 = zeros(Complex{T(w)}, WignerDsize(w))
    D!(𝔇, w, R)
end


@doc raw"""
    Y!(Y, wigner, s, R)
    Y!(wigner, s, R)

Evaluate (and write into `Y`, if present) the values of ``{}_{s}Y_{\ell,
m}(R)`` for the input value of `s`, for all ``(\ell, m)`` throughout the range
specified by `wigner`.  `R` is assumed to be a unit quaternion (which may be
`Rotor`, or simply a `Quaternion`).  If `R` does not have unit magnitude, the
output elements will be too large by a factor ``|R|^{\ell}``.  If `Y` is not
present, a new array will be created.

The spherical harmonics of spin weight ``s`` are related to Wigner's
``\mathfrak{D}`` matrix as
```math
\begin{aligned}
{}_{s}Y_{\ell, m}(R)
  &= (-1)^s \sqrt{\frac{2\ell+1}{4\pi}} \mathfrak{D}^{(\ell)}_{m, -s}(R) \\
  &= (-1)^s \sqrt{\frac{2\ell+1}{4\pi}} \bar{\mathfrak{D}}^{(\ell)}_{-s, m}(\bar{R}).
\end{aligned}
```
"""
function Y!(Y, w::WignerMatrixCalculator, s::Int, R::AbstractQuaternion)
    if length(Y) < Ysize(w)
        error("Input `Y` has length $(length(Y)); it should be at least $(Ysize(w))")
    end
    ell_min = ℓₘᵢₙ(w)
    ell_max = ℓₘₐₓ(w)
    mp_max = m′ₘₐₓ(w)
    if mp_max < abs(s)
        throw(DomainError("ℓₘₐₓ = $(ell_max)",
                          "Cannot compute sYlm for spin weight $s with m′ₘₐₓ only $(mp_max)"
        ))
    end

    to_euler_phases!(w.z, R)
    H!(w, w.z[2])
    complex_powers!(w.zₐpowers, w.z[1])
    complex_powers!(w.zᵧpowers, w.z[3])

    # Yˡₘₚ,ₘ(R) = dˡₘₚ,ₘ(R) exp[iϕₐ(m-mp)+iϕₛ(m+mp)] = dˡₘₚ,ₘ(R) exp[i(ϕₛ+ϕₐ)m+i(ϕₛ-ϕₐ)mp]
    # exp[iϕₛ] = R̂ₛ = hat(R[0] + 1j * R[3]) = zp
    # exp[iϕₐ] = R̂ₐ = hat(R[2] + 1j * R[1]) = zm.conjugate()
    # exp[i(ϕₛ+ϕₐ)] = zp * zm.conjugate() = z[2] = zᵧ
    # exp[i(ϕₛ-ϕₐ)] = zp * zm = z[0] = zₐ
    i_D = 1
    @inbounds for ell in ell_min:ell_max
        if ell < abs(s)
            for mp in -ell:ell
                Y[i_D] = 0
                i_D += 1
            end
        else
            factor = (isodd(s) ? -1 : 1) * √((2ell+1)/(4T(w)(π)))
            for mp in -ell:-1
                i_H = WignerHindex(ell, mp, -s, mp_max)
                if -s < 0
                    Y[i_D] = factor * ϵ(mp) * ϵ(s) * w.Hwedge[i_H] * conj(w.zᵧpowers[s+1]) * conj(w.zₐpowers[-mp+1])
                    # println((ell, mp, i_D, factor, ϵ(mp), ϵ(s), w.Hwedge[i_H], conj(w.zᵧpowers[s+1]), conj(w.zₐpowers[-mp+1])))
                else
                    Y[i_D] = factor * ϵ(mp) * ϵ(s) * w.Hwedge[i_H] * w.zᵧpowers[-s+1] * conj(w.zₐpowers[-mp+1])
                    # println((ell, mp, i_D, factor, ϵ(mp), ϵ(s), w.Hwedge[i_H], w.zᵧpowers[-s+1], conj(w.zₐpowers[-mp+1])))
                end
                i_D += 1
            end
            for mp in 0:ell
                i_H = WignerHindex(ell, mp, -s, mp_max)
                if -s < 0
                    Y[i_D] = factor * ϵ(mp) * ϵ(s) * w.Hwedge[i_H] * conj(w.zᵧpowers[s+1]) * w.zₐpowers[mp+1]
                    # println((ell, mp, i_D, factor, ϵ(mp), ϵ(s), w.Hwedge[i_H], conj(w.zᵧpowers[s+1]), w.zₐpowers[mp+1]))
                else
                    Y[i_D] = factor * ϵ(mp) * ϵ(s) * w.Hwedge[i_H] * w.zᵧpowers[-s+1] * w.zₐpowers[mp+1]
                    # println((ell, mp, i_D, factor, ϵ(mp), ϵ(s), w.Hwedge[i_H], w.zᵧpowers[-s+1], w.zₐpowers[mp+1]))
                end
                i_D += 1
            end
        end
    end
    Y
end


function Y!(w::WignerMatrixCalculator, s::Int, R::AbstractQuaternion)
    Y = zeros(Complex{T(w)}, Ysize(w))
    Y!(Y, w, s, R)
end
