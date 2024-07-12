### TODO: Test skipping all the complicated indexing tricks; use fancy indexing
### NOTES:
### 1. Caching coefficients provides a ~11% speedup at low ℓ, but that falls off
###    as ℓ increases, reaching ~4% around ℓ=512.
### 2. This is probably a much more significant advantage for ALFs.

using .SphericalFunctions: complex_powers!
using Quaternionic: AbstractQuaternion, to_euler_phases!

@inline ϵ(m) = ifelse(m > 0 && isodd(m), -1, 1)


@doc raw"""
    d_matrices(β, ℓₘₐₓ)
    d_matrices(expiβ, ℓₘₐₓ)

Compute Wigner's ``d^{(\ell)}`` matrices with elements ``d^{(\ell)}_{m',m}(\beta)`` for all
``\ell \leq \ell_\mathrm{max}``.  The ``d`` matrices are sometimes called the "reduced"
Wigner matrices, in contrast to the full ``\mathfrak{D}`` matrices.

See [`d_matrices!`](@ref) for details about the input and output values.

This function only appropriate when you need to evaluate the matrices for a single value of
`β` or `expiβ` because it allocates large arrays and performs many calculations that could
be reused.  If you need to evaluate the matrices for many values of `β` or `expiβ`, you
should pre-allocate the storage with [`d_prep`](@ref), and then call [`d_matrices!`](@ref)
with the result instead.

"""
d_matrices(β::Real, ℓₘₐₓ) = d_matrices(cis(β), ℓₘₐₓ)


@doc raw"""
    d_matrices!(d_storage, β)
    d_matrices!(d_storage, expiβ)
    d_matrices!(d, β, ℓₘₐₓ)
    d_matrices!(d, expiβ, ℓₘₐₓ)

Compute Wigner's ``d^{(\ell)}`` matrices with elements ``d^{(\ell)}_{m',m}(\beta)`` for all
``\ell \leq \ell_\mathrm{max}``.  The ``d`` matrices are sometimes called the "reduced"
Wigner matrices, in contrast to the full ``\mathfrak{D}`` matrices.

In all cases, the result is returned in a 1-dimensional array ordered as

    [
        dˡₘₚ,ₘ(β)
        for ℓ ∈ 0:ℓₘₐₓ
        for mp ∈ -ℓ:ℓ
        for m ∈ -ℓ:ℓ
    ]

When the first argument is `d`, it will be modified, so it must be at least as large as that
array.  When the first argument is `d_storage`, it should be the quantity returned by
[`d_prep`](@ref), and the result will be written into the `d` field of that tuple.  Both of
these options — especially the latter — reduce the number of allocations needed on each call
to the corresponding functions, which should increase the speed significantly.

!!! warn
    When using the `d_storage` argument (which is recommended), the returned quantity `d`
    will be an alias of `d_storage[1]`.  If you want to retain that data after the next call
    to [`d_matrices!`](@ref), you should copy it with `copy(d)`.

See also [`d_matrices`](@ref) for a simpler function call when you only need to evaluate the
matrices for a single value of `β` or `expiβ`.

# Examples

```julia
using SphericalFunctions
ℓₘₐₓ = 8
T = Float64
β = T(1)/8
d_storage = d_prep(ℓₘₐₓ, T)
d = d_matrices!(d_storage, β)
```

"""
d_matrices!(d, β::Real, ℓₘₐₓ) = d_matrices!(d, cis(β), ℓₘₐₓ)

d_matrices!(d_storage, β::Real) = d_matrices!(d_storage, cis(β))

function d_matrices(expiβ::Complex{T}, ℓₘₐₓ) where T
    d_storage = d_prep(ℓₘₐₓ, T)
    d = d_matrices!(d_storage, expiβ)
    return d
end

function d_matrices!(d, expiβ::Complex{T}, ℓₘₐₓ) where T
    d!(d, expiβ, ℓₘₐₓ, H_recursion_coefficients(ℓₘₐₓ, T))
    return d
end

function d_matrices!(d_storage, expiβ::Complex{T}) where T
    d, H_rec_coeffs, ℓₘₐₓ = d_storage
    d!(d, expiβ, ℓₘₐₓ, H_rec_coeffs)
    d
end

"""
    d_prep(ℓₘₐₓ, [T=Float64])

Construct space and pre-compute recursion coefficients to compute Wigner's ``d`` matrix in
place.

This returns the `d_storage` arguments needed by [`d_matrices!`](@ref).

"""
function d_prep(ℓₘₐₓ, ::Type{T}=Float64) where {T<:Real}
    d, H_rec_coeffs = dprep(ℓₘₐₓ, T)
    (d, H_rec_coeffs, ℓₘₐₓ)
end

# Legacy API for d_matrices
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
function dprep(ℓₘₐₓ, ::Type{T}) where {T<:Real}
    d = Vector{T}(undef, WignerDsize(ℓₘₐₓ))
    H_rec_coeffs = H_recursion_coefficients(ℓₘₐₓ, T)
    d, H_rec_coeffs
end


@doc raw"""
    D_matrices(R, ℓₘₐₓ)
    D_matrices(α, β, γ, ℓₘₐₓ)

Compute Wigner's 𝔇 matrices ``\mathfrak{D}^{(\ell)}_{m',m}(\beta)`` for all ``\ell \leq
\ell_\mathrm{max}``.

See [`D_matrices!`](@ref) for details about the input and output values.

This function only appropriate when you need to evaluate the matrices for a single value of
`R` or `α, β, γ` because it allocates large arrays and performs many calculations that could
be reused.  If you need to evaluate the matrices for many values of `R` or `α, β, γ`, you
should pre-allocate the storage with [`D_prep`](@ref), and then call [`D_matrices!`](@ref)
with the result instead.

"""
function D_matrices(R, ℓₘₐₓ)
    D_storage = D_prep(ℓₘₐₓ, eltype(R))
    D_matrices!(D_storage, R)
end

function D_matrices(α, β, γ, ℓₘₐₓ)
    T = promote_type(typeof.((α, β, γ))...)
    D_storage = D_prep(ℓₘₐₓ, T)
    D_matrices!(D_storage, R)
end

@doc raw"""
    D_matrices!(D_storage, R)
    D_matrices!(D_storage, α, β, γ)
    D_matrices!(D, R, ℓₘₐₓ)
    D_matrices!(D, α, β, γ, ℓₘₐₓ)

Compute Wigner's 𝔇 matrices ``\mathfrak{D}^{(\ell)}_{m',m}(\beta)`` for all ``\ell \leq
\ell_\mathrm{max}``.

In all cases, the result is returned in a 1-dimensional array ordered as

    [
        𝔇ˡₘₚ,ₘ(R)
        for ℓ ∈ 0:ℓₘₐₓ
        for mp ∈ -ℓ:ℓ
        for m ∈ -ℓ:ℓ
    ]

When the first argument is `D`, it will be modified, so it must be at least as large as that
array. When the first argument is `D_storage`, it should be the quantity returned by
[`D_prep`](@ref), and the result will be written into the `D` field of that tuple.  Both of
these options — especially the latter — reduce the number of allocations needed on each call
to the corresponding functions, which should increase the speed significantly.  Note that
the `D` or `D_storage` arguments must have types compatible with the type of `R` or `α, β,
γ`.

!!! warn
    When using the `D_storage` argument (which is recommended), the returned quantity `D`
    will be an alias of `D_storage[1]`.  If you want to retain that data after the next call
    to [`D_matrices!`](@ref), you should copy it with `copy(D)`.

The `α, β, γ` arguments are Euler angles as described in the documentation of
[`Quaternionic.from_euler_angles`](https://moble.github.io/Quaternionic.jl/dev/manual/#Quaternionic.from_euler_angles-Tuple{Any,%20Any,%20Any}).

See also [`D_matrices`](@ref) for a simpler function call when you only need to evaluate the
matrices for a single value of `R` or `α, β, γ`.

# Examples

```julia
using Quaternionic, SphericalFunctions
ℓₘₐₓ = 8
T = Float64
R = Rotor{T}(1, 2, 3, 4)  # Will be normalized automatically
D_storage = D_prep(ℓₘₐₓ, T)
D = D_matrices!(D_storage, R)
```

"""
function D_matrices!(D, R, ℓₘₐₓ)
    D_storage = (D, Dworkspace(ℓₘₐₓ, eltype(R))...)
    D_matrices!(D_storage, R)
end

function D_matrices!(D, α, β, γ, ℓₘₐₓ)
    T = promote_type(typeof.((α, β, γ))...)
    D_storage = (D, Dworkspace(ℓₘₐₓ, T)...)
    D_matrices!(D_storage, α, β, γ)
end

function D_matrices!(D_storage, R)
    (𝔇, ℓₘₐₓ, H_rec_coeffs, eⁱᵐᵅ, eⁱᵐᵞ) = D_storage
    D!(𝔇, R, ℓₘₐₓ, H_rec_coeffs, eⁱᵐᵅ, eⁱᵐᵞ)
end

function D_matrices!(D_storage, α, β, γ)
    (𝔇, ℓₘₐₓ, H_rec_coeffs, eⁱᵐᵅ, eⁱᵐᵞ) = D_storage
    D!(𝔇, cis(α), cis(β), cis(γ), ℓₘₐₓ, H_rec_coeffs, eⁱᵐᵅ, eⁱᵐᵞ)
end

@doc raw"""
    D_prep(ℓₘₐₓ, [T=Float64])

Construct storage space and pre-compute recursion coefficients to compute Wigner's
``\mathfrak{D}`` matrix in place.

This returns the `D_storage` arguments needed by [`D_matrices!`](@ref).

"""
function D_prep(ℓₘₐₓ, ::Type{T}=Float64) where {T<:Real}
    𝔇, H_rec_coeffs, eⁱᵐᵅ, eⁱᵐᵞ = Dprep(ℓₘₐₓ, T)
    (𝔇, ℓₘₐₓ, H_rec_coeffs, eⁱᵐᵅ, eⁱᵐᵞ)
end

# Legacy API for D_matrices
function D!(𝔇, R::AbstractQuaternion, ℓₘₐₓ, H_rec_coeffs, eⁱᵐᵅ, eⁱᵐᵞ)
    expiα, expiβ, expiγ = to_euler_phases(R)
    D!(𝔇, expiα, expiβ, expiγ, ℓₘₐₓ, H_rec_coeffs, eⁱᵐᵅ, eⁱᵐᵞ)
end
function D!(𝔇, expiα, expiβ, expiγ, ℓₘₐₓ, H_rec_coeffs, eⁱᵐᵅ, eⁱᵐᵞ)
    H!(𝔇, expiβ, ℓₘₐₓ, ℓₘₐₓ, H_rec_coeffs, WignerDindex)
    complex_powers!(eⁱᵐᵅ, expiα)
    complex_powers!(eⁱᵐᵞ, expiγ)

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
                𝔇[i1] = oddm_factor * 𝔇[i2] * conj(eⁱᵐᵞ[-m+1] * eⁱᵐᵅ[-m′+1])
            end
            for m′ in m+1:0
                i1 = i0 + (ℓ + m′) * (2ℓ + 1) + m + ℓ
                i2 = i0 + (ℓ - m′) * (2ℓ + 1) - m + ℓ
                𝔇[i1] = oddm_factor * 𝔇[i2] * conj(eⁱᵐᵞ[-m+1] * eⁱᵐᵅ[-m′+1])
            end
            for m′ in 1:-m
                i1 = i0 + (ℓ + m′) * (2ℓ + 1) + m + ℓ
                i2 = i0 + (ℓ - m′) * (2ℓ + 1) - m + ℓ
                𝔇[i1] = ifelse(isodd(m′), -1, 1) * oddm_factor * 𝔇[i2] * conj(eⁱᵐᵞ[-m+1]) * eⁱᵐᵅ[m′+1]
            end
            for m′ in 1-m:ℓ
                i1 = i0 + (ℓ + m′) * (2ℓ + 1) + m + ℓ
                i2 = i0 + (ℓ + m) * (2ℓ + 1) + m′ + ℓ
                𝔇[i1] = ifelse(isodd(m′), -1, 1) * oddm_factor * 𝔇[i2] * conj(eⁱᵐᵞ[-m+1]) * eⁱᵐᵅ[m′+1]
            end
        end
        for m in 0:ℓ
            for m′ in -ℓ:-m-1
                i1 = i0 + (ℓ + m′) * (2ℓ + 1) + m + ℓ
                i2 = i0 + (ℓ - m) * (2ℓ + 1) - m′ + ℓ
                𝔇[i1] = 𝔇[i2] * eⁱᵐᵞ[m+1] * conj(eⁱᵐᵅ[-m′+1])
            end
            for m′ in m+1:ℓ
                i1 = i0 + (ℓ + m′) * (2ℓ + 1) + m + ℓ
                i2 = i0 + (ℓ + m) * (2ℓ + 1) + m′ + ℓ
                𝔇[i1] = ifelse(isodd(m′), -𝔇[i2], 𝔇[i2]) * eⁱᵐᵞ[m+1] * eⁱᵐᵅ[m′+1]
            end
        end
        for m′ in -ℓ:0
            i1 = i0 + (ℓ + m′) * (2ℓ + 1) + ℓ
            for m in abs(m′):ℓ
                𝔇[i1+m] *= eⁱᵐᵞ[m+1] * conj(eⁱᵐᵅ[-m′+1])
            end
        end
        for m′ in 1:ℓ
            i1 = i0 + (ℓ + m′) * (2ℓ + 1) + ℓ
            for m in abs(m′):ℓ
                𝔇[i1+m] *= ifelse(isodd(m′), -1, 1) * eⁱᵐᵞ[m+1] * eⁱᵐᵅ[m′+1]
            end
        end
    end
    𝔇
end
function Dprep(ℓₘₐₓ, ::Type{T}) where {T<:Real}
    𝔇 = Vector{Complex{T}}(undef, WignerDsize(ℓₘₐₓ))
    H_rec_coeffs, eⁱᵐᵅ, eⁱᵐᵞ = Dworkspace(ℓₘₐₓ, T)
    𝔇, H_rec_coeffs, eⁱᵐᵅ, eⁱᵐᵞ
end
function Dworkspace(ℓₘₐₓ, ::Type{T}) where {T<:Real}
    H_rec_coeffs = H_recursion_coefficients(ℓₘₐₓ, T)
    eⁱᵐᵅ = Vector{Complex{T}}(undef, ℓₘₐₓ+1)
    eⁱᵐᵞ = Vector{Complex{T}}(undef, ℓₘₐₓ+1)
    H_rec_coeffs, eⁱᵐᵅ, eⁱᵐᵞ
end


@doc raw"""
    sYlm_values(R, ℓₘₐₓ, spin)
    sYlm_values(θ, ϕ, ℓₘₐₓ, spin)

Compute values of the spin-weighted spherical harmonic ``{}_{s}Y_{\ell, m}(R)`` for all
``\ell \leq \ell_\mathrm{max}``.

See [`sYlm_values!`](@ref) for details about the input and output values.

This function only appropriate when you need to evaluate the ``{}_{s}Y_{\ell, m}`` for a
single value of `R` or `θ, ϕ` because it allocates large arrays and performs many
calculations that could be reused.  If you need to evaluate the matrices for many values of
`R` or `θ, ϕ`, you should pre-allocate the storage with [`sYlm_prep`](@ref), and then call
[`sYlm_values!`](@ref) with the result instead.

"""
function sYlm_values(R::AbstractQuaternion, ℓₘₐₓ, spin)
    sYlm_storage = sYlm_prep(ℓₘₐₓ, spin, eltype(R), abs(spin))
    sYlm_values!(sYlm_storage, R, spin)
end

function sYlm_values(θ::Tθ, ϕ::Tϕ, ℓₘₐₓ, spin) where {Tθ<:Real, Tϕ<:Real}
    sYlm_storage = sYlm_prep(ℓₘₐₓ, spin, promote_type(Tθ, Tϕ), abs(spin))
    sYlm_values!(sYlm_storage, θ, ϕ, spin)
end

@doc raw"""
    sYlm_values!(sYlm_storage, R, spin)
    sYlm_values!(sYlm_storage, θ, ϕ, spin)
    sYlm_values!(sYlm, R, ℓₘₐₓ, spin)
    sYlm_values!(sYlm, θ, ϕ, ℓₘₐₓ, spin)

Compute values of the spin-weighted spherical harmonic ``{}_{s}Y_{\ell, m}(R)`` for all
``\ell \leq \ell_\mathrm{max}``.

The spherical harmonics of spin weight ``s`` are related to Wigner's ``\mathfrak{D}`` matrix
as
```math
\begin{aligned}
{}_{s}Y_{\ell, m}(R)
  &= (-1)^s \sqrt{\frac{2\ell+1}{4\pi}} \mathfrak{D}^{(\ell)}_{m, -s}(R) \\
  &= (-1)^s \sqrt{\frac{2\ell+1}{4\pi}} \bar{\mathfrak{D}}^{(\ell)}_{-s, m}(\bar{R}).
\end{aligned}
```

In all cases, the result is returned in a 1-dimensional array ordered as

    [
        ₛYₗₘ(R)
        for ℓ ∈ 0:ℓₘₐₓ
        for m ∈ -ℓ:ℓ
    ]

When the first argument is `Y`, it will be modified, so it must be at least as large as that
array. When the first argument is `sYlm_storage`, it should be the quantity returned by
[`sYlm_prep`](@ref), and the result will be written into the `Y` field of that tuple.  Both
of these options — especially the latter — reduce the number of allocations needed on each
call to the corresponding functions, which should increase the speed significantly.  Note
that the `Y` or `sYlm_storage` arguments must have types compatible with the type of `R` or
`θ, ϕ`.

!!! warn
    When using the `sYlm_storage` argument (which is recommended), the returned quantity
    `sYlm` will be an alias of `sYlm_storage[1]`.  If you want to retain that data after the
    next call to [`sYlm_values!`](@ref), you should copy it with `copy(sYlm)`.

The `θ, ϕ` arguments are spherical coordinates as described in the documentation of
[`Quaternionic.from_spherical_coordinates`](https://moble.github.io/Quaternionic.jl/dev/manual/#Quaternionic.from_spherical_coordinates-Tuple{Any,%20Any}).

See also [`sYlm_values`](@ref) for a simpler function call when you only need to evaluate
the ``{}_{s}Y_{\ell, m}`` for a single value of `R` or `θ, ϕ`.

# Examples

```julia
using Quaternionic, SphericalFunctions
spin = -2
ℓₘₐₓ = 8
T = Float64
R = Rotor{T}(1, 2, 3, 4)  # Will be normalized automatically
sYlm_storage = sYlm_prep(ℓₘₐₓ, spin, T)
sYlm = sYlm_values!(sYlm_storage, R, spin)
```

"""
function sYlm_values!(Y, R::AbstractQuaternion, ℓₘₐₓ, spin)
    sYlm_storage = (Y, Y_workspace(ℓₘₐₓ, spin, eltype(R), abs(spin))...)
    sYlm_values!(sYlm_storage, R, spin)
end

function sYlm_values!(Y, θ::Tθ, ϕ::Tϕ, ℓₘₐₓ, spin) where {Tθ<:Real, Tϕ<:Real}
    sYlm_storage = (Y, Y_workspace(ℓₘₐₓ, spin, promote_type(Tθ, Tϕ), abs(spin))...)
    sYlm_values!(sYlm_storage, θ, ϕ, spin)
end

function sYlm_values!(sYlm_storage, R::AbstractQuaternion, spin)
    (Y, ℓₘₐₓ, sₘₐₓ, H_rec_coeffs, Hwedge, expimϕ, ℓₘᵢₙ) = sYlm_storage
    if abs(spin) > abs(sₘₐₓ)
        error(
            "Input `sYlm_storage` was built created for maximum spin `sₘₐₓ`=$(sₘₐₓ), "
            *"but `spin`=$(spin) was requested."
        )
    end
    Y!(Y, R, ℓₘₐₓ, spin, H_rec_coeffs, Hwedge, expimϕ, ℓₘᵢₙ)
end

function sYlm_values!(sYlm_storage, θ::Tθ, ϕ::Tϕ, spin) where {Tθ<:Real, Tϕ<:Real}
    (Y, ℓₘₐₓ, sₘₐₓ, H_rec_coeffs, Hwedge, expimϕ, ℓₘᵢₙ) = sYlm_storage
    if abs(spin) > abs(sₘₐₓ)
        error(
            "Input `sYlm_storage` was created for maximum spin `sₘₐₓ`=$(sₘₐₓ), "
            *"but `spin`=$(spin) was requested."
        )
    end
    expiθ, expiϕ = cis.(promote(θ, ϕ))
    expiγ = one(typeof(expiθ))
    Y!(Y, expiϕ, expiθ, expiγ, ℓₘₐₓ, spin, H_rec_coeffs, Hwedge, expimϕ, ℓₘᵢₙ)
end

function Y_workspace(ℓₘₐₓ, sₘₐₓ, ::Type{T}, ℓₘᵢₙ=0) where {T<:Real}
    H_rec_coeffs = H_recursion_coefficients(ℓₘₐₓ, T)
    Hwedge = Vector{T}(undef, WignerHsize(ℓₘₐₓ, abs(sₘₐₓ)))
    expimϕ = Vector{Complex{T}}(undef, ℓₘₐₓ+1)
    ℓₘₐₓ, sₘₐₓ, H_rec_coeffs, Hwedge, expimϕ, ℓₘᵢₙ
end

@doc raw"""
    sYlm_prep(ℓₘₐₓ, sₘₐₓ, [T=Float64, [ℓₘᵢₙ=0]])

Construct storage space and pre-compute recursion coefficients to compute spin-weighted
spherical-harmonic values ``{}_{s}Y_{\ell, m}`` in place.

This returns the `sYlm_storage` arguments needed by [`sYlm_values!`](@ref).

Note that the result of this function can be passed to `sYlm_values!`, even if the value of
`spin` passed to that function is smaller (in absolute value) than the `sₘₐₓ` passed to this
function.  That is, the `sYlm_storage` returned by this function can be used to compute
``{}_{s}Y_{\ell, m}`` values for numerous values of the spin.

"""
function sYlm_prep(ℓₘₐₓ, sₘₐₓ, ::Type{T}=Float64, ℓₘᵢₙ=0) where {T<:Real}
    Y = Vector{Complex{T}}(undef, Ysize(ℓₘᵢₙ, ℓₘₐₓ))
    ℓₘₐₓ, sₘₐₓ, H_rec_coeffs, Hwedge, expimϕ, ℓₘᵢₙ = Y_workspace(ℓₘₐₓ, sₘₐₓ, T, ℓₘᵢₙ)
    (Y, ℓₘₐₓ, sₘₐₓ, H_rec_coeffs, Hwedge, expimϕ, ℓₘᵢₙ)
end

# Legacy API for sYlm_values
function Y!(Y, R, ℓₘₐₓ, spin, H_rec_coeffs, Hwedge, expimϕ, ℓₘᵢₙ=0)
    expiϕ, expiθ, expiγ = to_euler_phases(R)
    Y!(Y, expiϕ, expiθ, expiγ, ℓₘₐₓ, spin, H_rec_coeffs, Hwedge, expimϕ, ℓₘᵢₙ)
end
function Y!(Y, expiϕ, expiθ, expiγ, ℓₘₐₓ, spin, H_rec_coeffs, Hwedge, expimϕ, ℓₘᵢₙ=0)
    if length(Y) < Ysize(ℓₘᵢₙ, ℓₘₐₓ)
        error("Input `Y` has length $(length(Y)); which is not enough for ℓₘₐₓ=$ℓₘₐₓ")
    end
    if length(Hwedge) < WignerHsize(ℓₘₐₓ, abs(spin))
        error(
            "Input `Hwedge` has length $(length(Hwedge)); "
            *"which is not enough for ℓₘₐₓ=$ℓₘₐₓ with spin=$spin"
        )
    end

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
function Yprep(ℓₘₐₓ, sₘₐₓ, ::Type{T}, ℓₘᵢₙ=0) where {T<:Real}
    sYlm_storage = sYlm_prep(ℓₘₐₓ, sₘₐₓ, T, ℓₘᵢₙ)
    # Y, H_rec_coeffs, Hwedge, expimϕ
    sYlm_storage[1], sYlm_storage[4], sYlm_storage[5], sYlm_storage[6]
end
@doc raw"""
    ₛ𝐘(s, ℓₘₐₓ, [T=Float64], [Rθϕ=golden_ratio_spiral_rotors(s, ℓₘₐₓ, T)])


Construct a matrix of ``ₛYₗₘ(Rθϕ)`` values for the input `s` and all nontrivial ``(\ell,
m)`` up to `ℓₘₐₓ`.

This is a fast and accurate method for mapping between the vector of spin-weighted
spherical-harmonic mode weights ``ₛ𝐟ₗₘ`` and the vector of function values on the sphere
``ₛ𝐟ⱼₖ``, as
```math
ₛ𝐟ⱼₖ = ₛ𝐘\, ₛ𝐟ₗₘ,
```
where the right-hand side represents the matrix-vector product.  As usual, we assume that
the ``ₛ𝐟ₗₘ`` modes are ordered by increasing ``m ∈ [-ℓ:ℓ]``, and ``ℓ ∈ [|s|:ℓₘₐₓ]``.  The
ordering of the ``ₛ𝐟ⱼₖ`` values will be determined by the ordering of the argument `Rθϕ`.

Note that the number of modes need not be the same as the number of points on which the
function is evaluated, which would imply that the output matrix is not square.  To be able
to invert the relationship, however, we need the number of points ``ₛ𝐟ⱼₖ`` to be *at least
as large* as the number of modes ``ₛ𝐟ₗₘ``.

Note that the usefulness of this approach is limited by the fact that the size of this
matrix scales as ℓₘₐₓ⁴.  As such, it is mostly useful only for ℓₘₐₓ of order dozens, rather
than — say — the tens of thousands that CMB astronomy or lensing require, for example.

Direct application and inversion of this matrix are used in the "direct" methods of
``s``-SHT transformations.  See [`SSHTDirect`](@ref) for details about the implementation.

"""
function ₛ𝐘(s, ℓₘₐₓ, ::Type{T}=Float64, Rθϕ=golden_ratio_spiral_rotors(s, ℓₘₐₓ, T)) where T
    Y, H_rec_coeffs, Hwedge, expimϕ = Yprep(ℓₘₐₓ, s, T, abs(s))
    ₛ𝐘 = zeros(Complex{T}, length(Rθϕ), length(Y))
    for (j1,j2) ∈ zip(axes(ₛ𝐘, 1), axes(Rθϕ, 1))
        ₛ𝐘[j1, :] = Y!(Y, Rθϕ[j2], ℓₘₐₓ, s, H_rec_coeffs, Hwedge, expimϕ, abs(s))
    end
    ₛ𝐘
end
