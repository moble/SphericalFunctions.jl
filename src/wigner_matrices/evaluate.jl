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



#export WignerMatrixCalculator,
export H!, abd
export WignerDsize, WignerHsize, WignerDindex, WignerHindex, _WignerHindex
#export d!, D!, Y!, ϵ

using ..SphericalFunctions: complex_powers!
using Quaternionic: AbstractQuaternion, to_euler_phases!

include("indexing.jl")
# include("calculator.jl")
#include("Hrecursions.jl")
include("Hrecursor.jl")


@inline ϵ(m) = (m <= 0 ? 1 : (isodd(m) ? -1 : 1))


"""Return sign of input, with sign(0)=1"""
@inline sign(m) = (m < 0 ? -1 : 1)


"""
    d!(d, expiβ, ℓₘₐₓ, abd)
    d!(d, expiβ, ℓₘₐₓ)
    d!(d, β, ℓₘₐₓ, abd)
    d!(d, β, ℓₘₐₓ)

Compute Wigner's d matrix dˡₘₚ,ₘ(β)

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
function d!(d, expiβ::Complex, ℓₘₐₓ, (avals,bvals,dvals))
    H!(d, expiβ, ℓₘₐₓ, ℓₘₐₓ, (avals,bvals,dvals), WignerDindex)

    @inbounds for ℓ in 0:ℓₘₐₓ
        i0 = WignerDindex(ℓ, -ℓ, -ℓ)
        for m in -ℓ:-1
            oddm = isodd(m)
            for m′ in -ℓ:m
                i1 = i0 + (ℓ + m′) * (2ℓ + 1) + m + ℓ
                i2 = i0 + (ℓ - m) * (2ℓ + 1) - m′ + ℓ
                d[i1] = ifelse(oddm, -d[i2], d[i2])
            end
            for m′ in m+1:0
                i1 = i0 + (ℓ + m′) * (2ℓ + 1) + m + ℓ
                i2 = i0 + (ℓ - m′) * (2ℓ + 1) - m + ℓ
                d[i1] = ifelse(oddm, -d[i2], d[i2])
            end
            for m′ in 1:-m
                i1 = i0 + (ℓ + m′) * (2ℓ + 1) + m + ℓ
                i2 = i0 + (ℓ - m′) * (2ℓ + 1) - m + ℓ
                d[i1] = ifelse(isodd(m′)⊻oddm, -d[i2], d[i2])
            end
            for m′ in 1-m:ℓ
                i1 = i0 + (ℓ + m′) * (2ℓ + 1) + m + ℓ
                i2 = i0 + (ℓ + m) * (2ℓ + 1) + m′ + ℓ
                d[i1] = ifelse(isodd(m′)⊻oddm, -d[i2], d[i2])
            end
        end
        for m in 0:ℓ
            for m′ in -ℓ:-m-1
                i1 = i0 + (ℓ + m′) * (2ℓ + 1) + m + ℓ
                i2 = i0 + (ℓ - m) * (2ℓ + 1) - m′ + ℓ
                d[i1] = d[i2]
            end
            for m′ in m+1:ℓ
                i1 = i0 + (ℓ + m′) * (2ℓ + 1) + m + ℓ
                i2 = i0 + (ℓ + m) * (2ℓ + 1) + m′ + ℓ
                d[i1] = ifelse(isodd(m′), -d[i2], d[i2])
            end
        end
        for m′ in 1:2:ℓ
            i1 = i0 + (ℓ + m′) * (2ℓ + 1) + ℓ
            for m in abs(m′):ℓ
                d[i1+m] *= -1
            end
        end
    end
    d
end
d!(d, expiβ::Complex, ℓₘₐₓ) = d!(d, expiβ, ℓₘₐₓ, abd(ℓₘₐₓ, real(typeof(β))))
d!(d, β::Real, ℓₘₐₓ, (avals,bvals,dvals)) = d!(d, exp(im*β), ℓₘₐₓ, (avals,bvals,dvals))
d!(d, β::Real, ℓₘₐₓ) = d!(d, exp(im*β), ℓₘₐₓ, abd(ℓₘₐₓ, typeof(β)))


"""
    D!(𝔇, R, ℓₘₐₓ, (avals,bvals,dvals), expimα, expimγ)

Compute Wigner's 𝔇 matrix

This function implements the preferred method of computing the 𝔇 matrix for large ell
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
function D!(𝔇, R::AbstractQuaternion, ℓₘₐₓ, (avals,bvals,dvals), expimα, expimγ)
    expiα, expiβ, expiγ = to_euler_phases(R)
    H!(𝔇, expiβ, ℓₘₐₓ, ℓₘₐₓ, (avals,bvals,dvals), WignerDindex)
    complex_powers!(expimα, expiα)
    complex_powers!(expimγ, expiγ)

    # 𝔇ˡₘₚ,ₘ(R) = dˡₘₚ,ₘ(R) exp[iϕₐ(m-mp)+iϕₛ(m+mp)] = dˡₘₚ,ₘ(R) exp[i(ϕₛ+ϕₐ)m+i(ϕₛ-ϕₐ)mp]
    # exp[iϕₛ] = R̂ₛ = hat(R[0] + 1j * R[3]) = zp
    # exp[iϕₐ] = R̂ₐ = hat(R[2] + 1j * R[1]) = zm.conjugate()
    # exp[i(ϕₛ+ϕₐ)] = zp * zm.conjugate() = z[2] = zᵧ
    # exp[i(ϕₛ-ϕₐ)] = zp * zm = z[0] = zₐ
    @inbounds for ℓ in 0:ℓₘₐₓ
        i0 = WignerDindex(ℓ, -ℓ, -ℓ)
        for m in -ℓ:-1
            oddm_factor = ifelse(isodd(m), -1, 1)
            for m′ in -ℓ:m
                i1 = i0 + (ℓ + m′) * (2ℓ + 1) + m + ℓ
                i2 = i0 + (ℓ - m) * (2ℓ + 1) - m′ + ℓ
                𝔇[i1] = oddm_factor * 𝔇[i2] * conj(expimγ[-m+1] * expimα[-m′+1])
            end
            for m′ in m+1:0
                i1 = i0 + (ℓ + m′) * (2ℓ + 1) + m + ℓ
                i2 = i0 + (ℓ - m′) * (2ℓ + 1) - m + ℓ
                𝔇[i1] = oddm_factor * 𝔇[i2] * conj(expimγ[-m+1] * expimα[-m′+1])
            end
            for m′ in 1:-m
                i1 = i0 + (ℓ + m′) * (2ℓ + 1) + m + ℓ
                i2 = i0 + (ℓ - m′) * (2ℓ + 1) - m + ℓ
                𝔇[i1] = ifelse(isodd(m′), -1, 1) * oddm_factor * 𝔇[i2] * conj(expimγ[-m+1]) * expimα[m′+1]
            end
            for m′ in 1-m:ℓ
                i1 = i0 + (ℓ + m′) * (2ℓ + 1) + m + ℓ
                i2 = i0 + (ℓ + m) * (2ℓ + 1) + m′ + ℓ
                𝔇[i1] = ifelse(isodd(m′), -1, 1) * oddm_factor * 𝔇[i2] * conj(expimγ[-m+1]) * expimα[m′+1]
            end
        end
        for m in 0:ℓ
            for m′ in -ℓ:-m-1
                i1 = i0 + (ℓ + m′) * (2ℓ + 1) + m + ℓ
                i2 = i0 + (ℓ - m) * (2ℓ + 1) - m′ + ℓ
                𝔇[i1] = 𝔇[i2] * expimγ[m+1] * conj(expimα[-m′+1])
            end
            for m′ in m+1:ℓ
                i1 = i0 + (ℓ + m′) * (2ℓ + 1) + m + ℓ
                i2 = i0 + (ℓ + m) * (2ℓ + 1) + m′ + ℓ
                𝔇[i1] = ifelse(isodd(m′), -𝔇[i2], 𝔇[i2]) * expimγ[m+1] * expimα[m′+1]
            end
        end
        for m′ in -ℓ:0
            i1 = i0 + (ℓ + m′) * (2ℓ + 1) + ℓ
            for m in abs(m′):ℓ
                𝔇[i1+m] *= expimγ[m+1] * conj(expimα[-m′+1])
            end
        end
        for m′ in 1:ℓ
            i1 = i0 + (ℓ + m′) * (2ℓ + 1) + ℓ
            for m in abs(m′):ℓ
                𝔇[i1+m] *= ifelse(isodd(m′), -1, 1) * expimγ[m+1] * expimα[m′+1]
            end
        end
    end
    𝔇
end


# @doc raw"""
#     Y!(Y, wigner, s, R)
#     Y!(wigner, s, R)

# Evaluate (and write into `Y`, if present) the values of ``{}_{s}Y_{\ell,
# m}(R)`` for the input value of `s`, for all ``(\ell, m)`` throughout the range
# specified by `wigner`.  `R` is assumed to be a unit quaternion (which may be
# `Rotor`, or simply a `Quaternion`).  If `R` does not have unit magnitude, the
# output elements will be too large by a factor ``|R|^{\ell}``.  If `Y` is not
# present, a new array will be created.

# The spherical harmonics of spin weight ``s`` are related to Wigner's
# ``\mathfrak{D}`` matrix as
# ```math
# \begin{aligned}
# {}_{s}Y_{\ell, m}(R)
#   &= (-1)^s \sqrt{\frac{2\ell+1}{4\pi}} \mathfrak{D}^{(\ell)}_{m, -s}(R) \\
#   &= (-1)^s \sqrt{\frac{2\ell+1}{4\pi}} \bar{\mathfrak{D}}^{(\ell)}_{-s, m}(\bar{R}).
# \end{aligned}
# ```
# """
# function Y!(Y, w::WignerMatrixCalculator, s::Int, R::AbstractQuaternion)
#     if length(Y) < Ysize(w)
#         error("Input `Y` has length $(length(Y)); it should be at least $(Ysize(w))")
#     end
#     ell_min = ℓₘᵢₙ(w)
#     ell_max = ℓₘₐₓ(w)
#     mp_max = m′ₘₐₓ(w)
#     if mp_max < abs(s)
#         throw(DomainError("ℓₘₐₓ = $(ell_max)",
#                           "Cannot compute sYlm for spin weight $s with m′ₘₐₓ only $(mp_max)"
#         ))
#     end

#     to_euler_phases!(w.z, R)
#     H!(w, w.z[2])
#     complex_powers!(w.zₐpowers, w.z[1])
#     complex_powers!(w.zᵧpowers, w.z[3])

#     # Yˡₘₚ,ₘ(R) = dˡₘₚ,ₘ(R) exp[iϕₐ(m-mp)+iϕₛ(m+mp)] = dˡₘₚ,ₘ(R) exp[i(ϕₛ+ϕₐ)m+i(ϕₛ-ϕₐ)mp]
#     # exp[iϕₛ] = R̂ₛ = hat(R[0] + 1j * R[3]) = zp
#     # exp[iϕₐ] = R̂ₐ = hat(R[2] + 1j * R[1]) = zm.conjugate()
#     # exp[i(ϕₛ+ϕₐ)] = zp * zm.conjugate() = z[2] = zᵧ
#     # exp[i(ϕₛ-ϕₐ)] = zp * zm = z[0] = zₐ
#     i_D = 1
#     @inbounds for ell in ell_min:ell_max
#         if ell < abs(s)
#             for mp in -ell:ell
#                 Y[i_D] = 0
#                 i_D += 1
#             end
#         else
#             factor = (isodd(s) ? -1 : 1) * √((2ell+1)/(4T(w)(π)))
#             for mp in -ell:-1
#                 i_H = WignerHindex(ell, mp, -s, mp_max)
#                 if -s < 0
#                     Y[i_D] = factor * ϵ(mp) * ϵ(s) * w.Hwedge[i_H] * conj(w.zᵧpowers[s+1]) * conj(w.zₐpowers[-mp+1])
#                     # println((ell, mp, i_D, factor, ϵ(mp), ϵ(s), w.Hwedge[i_H], conj(w.zᵧpowers[s+1]), conj(w.zₐpowers[-mp+1])))
#                 else
#                     Y[i_D] = factor * ϵ(mp) * ϵ(s) * w.Hwedge[i_H] * w.zᵧpowers[-s+1] * conj(w.zₐpowers[-mp+1])
#                     # println((ell, mp, i_D, factor, ϵ(mp), ϵ(s), w.Hwedge[i_H], w.zᵧpowers[-s+1], conj(w.zₐpowers[-mp+1])))
#                 end
#                 i_D += 1
#             end
#             for mp in 0:ell
#                 i_H = WignerHindex(ell, mp, -s, mp_max)
#                 if -s < 0
#                     Y[i_D] = factor * ϵ(mp) * ϵ(s) * w.Hwedge[i_H] * conj(w.zᵧpowers[s+1]) * w.zₐpowers[mp+1]
#                     # println((ell, mp, i_D, factor, ϵ(mp), ϵ(s), w.Hwedge[i_H], conj(w.zᵧpowers[s+1]), w.zₐpowers[mp+1]))
#                 else
#                     Y[i_D] = factor * ϵ(mp) * ϵ(s) * w.Hwedge[i_H] * w.zᵧpowers[-s+1] * w.zₐpowers[mp+1]
#                     # println((ell, mp, i_D, factor, ϵ(mp), ϵ(s), w.Hwedge[i_H], w.zᵧpowers[-s+1], w.zₐpowers[mp+1]))
#                 end
#                 i_D += 1
#             end
#         end
#     end
#     Y
# end


# function Y!(w::WignerMatrixCalculator, s::Int, R::AbstractQuaternion)
#     Y = zeros(Complex{T(w)}, Ysize(w))
#     Y!(Y, w, s, R)
# end
