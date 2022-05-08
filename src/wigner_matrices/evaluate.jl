### TODO:
### 2. Separate ALF computation to a different module
### 4. Test skipping all the complicated indexing tricks; use fancy indexing
### 5. Allow specifying the extent of recursion / iterating over ℓ matrices


### NOTES:
### 1. Caching coefficients provides a ~11% speedup at low ℓ, but that falls off
###    as ℓ increases, reaching ~4% around ℓ=512.
### 2. This is probably a much more significant advantage for ALFs.



export H!, H_recursion_coefficients
export WignerDsize, WignerHsize, WignerDindex, WignerHindex, _WignerHindex
export d!, d, D!, Y!
export dprep, dstorage, Dprep, Dstorage, Yprep, Ystorage#, ϵ

using ..SphericalFunctions: complex_powers!
using Quaternionic: AbstractQuaternion, to_euler_phases!

include("indexing.jl")
include("Hrecursor.jl")


@inline ϵ(m) = ifelse(m > 0 && isodd(m), -1, 1)


"""
    d!(d, expiβ, ℓₘₐₓ, H_rec_coeffs)
    d!(d, expiβ, ℓₘₐₓ)
    d!(d, β, ℓₘₐₓ, H_rec_coeffs)
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
function d!(d, expiβ::Complex, ℓₘₐₓ, H_rec_coeffs)
    H!(d, expiβ, ℓₘₐₓ, ℓₘₐₓ, H_rec_coeffs, WignerDindex)

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
function d!(d, expiβ::Complex{T}, ℓₘₐₓ) where {T<:Real}
    d!(d, expiβ, ℓₘₐₓ, H_recursion_coefficients(ℓₘₐₓ, T))
end
function d!(d, β::T, ℓₘₐₓ, H_rec_coeffs) where {T<:Real}
    d!(d, cis(β), ℓₘₐₓ, H_rec_coeffs)
end
function d!(d, β::T, ℓₘₐₓ) where {T<:Real}
    d!(d, cis(β), ℓₘₐₓ, H_recursion_coefficients(ℓₘₐₓ, T))
end
function d(expiβ::Complex{T}, ℓₘₐₓ) where {T<:Real}
    𝔡 = Array{T}(undef, WignerDsize(ℓₘₐₓ, ℓₘₐₓ))
    d!(𝔡, expiβ, ℓₘₐₓ, H_recursion_coefficients(ℓₘₐₓ, T))
end
d(β::T, ℓₘₐₓ) where {T<:Real} = d(cis(β), ℓₘₐₓ)

function dstorage(ℓₘₐₓ, ::Type{T}) where {T<:Real}
    Vector{T}(undef, WignerDsize(ℓₘₐₓ))
end

function dprep(ℓₘₐₓ, ::Type{T}) where {T<:Real}
    d = dstorage(ℓₘₐₓ, T)
    H_rec_coeffs = H_recursion_coefficients(ℓₘₐₓ, T)
    d, H_rec_coeffs
end


"""
    D!(𝔇, R, ℓₘₐₓ, (aₙᵐ,bₙᵐ,dₙᵐ), expimα, expimγ)

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
function D!(𝔇, R::AbstractQuaternion, ℓₘₐₓ, H_rec_coeffs, expimα, expimγ)
    expiα, expiβ, expiγ = to_euler_phases(R)
    H!(𝔇, expiβ, ℓₘₐₓ, ℓₘₐₓ, H_rec_coeffs, WignerDindex)
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

function Dstorage(ℓₘₐₓ, ::Type{T}) where {T<:Real}
    Vector{Complex{T}}(undef, WignerDsize(ℓₘₐₓ))
end

function Dprep(ℓₘₐₓ, ::Type{T}) where {T<:Real}
    𝔇 = Dstorage(ℓₘₐₓ, T)
    H_rec_coeffs = H_recursion_coefficients(ℓₘₐₓ, T)
    expimα = Vector{Complex{T}}(undef, ℓₘₐₓ+1)
    expimγ = Vector{Complex{T}}(undef, ℓₘₐₓ+1)
    𝔇, H_rec_coeffs, expimα, expimγ
end


@doc raw"""
    Y(s, R, ℓₘₐₓ)
    Y!(Y, s, R, ℓₘₐₓ)

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
function Y!(Y, R, ℓₘₐₓ, spin, H_rec_coeffs, Hwedge, expimϕ, ℓₘᵢₙ=0)
    if length(Y) < Ysize(ℓₘᵢₙ, ℓₘₐₓ)
        error("Input `Y` has length $(length(Y)); which is not enough for ℓₘₐₓ=$ℓₘₐₓ")
    end
    if length(Hwedge) < WignerHsize(ℓₘₐₓ, abs(spin))
        error(
            "Input `Hwedge` has length $(length(Hwedge)); "
            *"which is not enough for ℓₘₐₓ=$ℓₘₐₓ with spin=$spin"
        )
    end

    expiϕ, expiθ, expiγ = to_euler_phases(R)
    H!(Hwedge, expiθ, ℓₘₐₓ, abs(spin), H_rec_coeffs)
    complex_powers!(expimϕ, expiϕ)
    spin_factor = (isodd(spin) ? -1 : 1) * ϵ(spin) * expiγ^-spin

    # Yˡₘₚ,ₘ(R) ∝ dˡₘₚ,ₘ(R) exp[iϕₐ(m-mp)+iϕₛ(m+mp)] = dˡₘₚ,ₘ(R) exp[i(ϕₛ+ϕₐ)m+i(ϕₛ-ϕₐ)mp]
    # exp[iϕₛ] = R̂ₛ = hat(R[0] + 1j * R[3]) = zp
    # exp[iϕₐ] = R̂ₐ = hat(R[2] + 1j * R[1]) = zm.conjugate()
    # exp[i(ϕₛ+ϕₐ)] = zp * zm.conjugate() = z[2] = zᵧ
    # exp[i(ϕₛ-ϕₐ)] = zp * zm = z[0] = zₐ
    iᴰ = 1
    @inbounds for ℓ in ℓₘᵢₙ:ℓₘₐₓ
        if ℓ < abs(spin)
            Y[iᴰ:iᴰ+2ℓ] .= 0
            iᴰ += 2ℓ+1
        else
            factor = spin_factor * √((2ℓ+1)/(4eltype(Y)(π)))
            for m in -ℓ:-1
                iᴴ = WignerHindex(ℓ, m, -spin, abs(spin))
                Y[iᴰ] = factor * Hwedge[iᴴ] * conj(expimϕ[-m+1])  # ϵ(m′)==1
                iᴰ += 1
            end
            for m in 0:ℓ
                iᴴ = WignerHindex(ℓ, m, -spin, abs(spin))
                Y[iᴰ] = factor * ϵ(m) * Hwedge[iᴴ] * expimϕ[m+1]
                iᴰ += 1
            end
        end
    end
    Y
end

function Ystorage(ℓₘₐₓ, ::Type{T}, ℓₘᵢₙ=0) where {T<:Real}
    Vector{Complex{T}}(undef, Ysize(ℓₘᵢₙ, ℓₘₐₓ))
end

function Yworkspace(ℓₘₐₓ, sₘₐₓ, ::Type{T}) where {T<:Real}
    Hwedge = Vector{T}(undef, WignerHsize(ℓₘₐₓ, abs(sₘₐₓ)))
    expimϕ = Vector{Complex{T}}(undef, ℓₘₐₓ+1)
    Hwedge, expimϕ
end

"""
    Yprep(ℓₘₐₓ, sₘₐₓ, T, ℓₘᵢₙ)

Prepare the storage, recursion coefficients, and workspace to compute ₛYₗ,ₘ
data up to the maximum sizes given.

Returns a tuple of `Y, H_rec_coeffs, Hwedge, expimϕ`, which can be passed to
the correspondingly named arguments of `Y!`.

Note that the same results of this function can be passed to `Y!`, even if the
value of `ℓₘₐₓ` passed to that function is smaller than the value passed to
this function, or the value of `spin` passed to that function is smaller (in
absolute value) than the `sₘₐₓ` passed to this function.  However, the value
of `ℓₘᵢₙ` passed to that function *must not* be smaller than the value passed
to this function (unless one of the other sizes is sufficiently smaller).

"""
function Yprep(ℓₘₐₓ, sₘₐₓ, ::Type{T}, ℓₘᵢₙ=0) where {T<:Real}
    Y = Ystorage(ℓₘₐₓ, T, ℓₘᵢₙ)
    H_rec_coeffs = H_recursion_coefficients(ℓₘₐₓ, T)
    Hwedge, expimϕ = Yworkspace(ℓₘₐₓ, sₘₐₓ, T)
    Y, H_rec_coeffs, Hwedge, expimϕ
end
