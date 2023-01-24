
"""Storage for EKKM spin-spherical-harmonic transform

The EKKM algorithm was described in [this paper](https://arxiv.org/abs/1809.01321), and allows for
the minimal number of function samples.

"""
struct SSHTEKKM{T<:Real, Inplace} <: SSHT{T}
    """Spin weight"""
    s::Integer

    """Highest ℓ value present in the data"""
    ℓₘₐₓ::Integer

    """OffsetVector of colatitudes of sampling rings

    These are the coordinates of the rings, where ring `j` contains 2j+1
    equally spaced points.  The index ranges over j ∈ abs(s):ℓₘₐₓ.
    """
    θ::OffsetVector

    """OffsetVector of index ranges for each θ ring

    This OffsetVector provides a series of `UnitRange`s indexing each
    successive colatitude ring (as given by `θ`), so that the data along
    the first such ring will be given by `f[θindices[abs(s)]]`, and the
    last `f[θindices[ℓₘₐₓ]]`, where `f` is a standard (non-offset)
    `Vector`.
    """
    θindices::OffsetVector

    # """Spin-weighted spherical harmonic values"""
    # ₛ𝐘

    """OffsetVector of Fourier-transform plans

    These transform physical-space function values to the Fourier domain on
    each ring of colatitude.  The indices range over j ∈ abs(s):ℓₘₐₓ.
    """
    plans::OffsetVector

    """OffsetVector of LU-decomposed matrices to solve for mode weights

    There is one matrix for each value of m = k ∈ -ℓₘₐₓ:ℓₘₐₓ; the rows in each
    matrix correspond to j ∈ max(abs(s),abs(m)):ℓₘₐₓ; the columns correspond
    to ℓ over the same range as j.
    """
    #ₛ𝐝::OffsetVector
    ₛΛ::OffsetVector

    """Preallocated storage for solving linear equations"""
    workspace::Vector

    """Preallocated storage for all FT modes with a given positive ``m`` value"""
    ₛf̃₊ₘ::OffsetVector

    """Preallocated storage for all FT modes with a given negative ``m`` value"""
    ₛf̃₋ₘ::OffsetVector

    """Preallocated storage for all SH modes with a given positive ``m`` value"""
    ₛf̃₊::OffsetVector

    """Preallocated storage for all SH modes with a given negative ``m`` value"""
    ₛf̃₋::OffsetVector
end

function SSHTEKKM(
    s, ℓₘₐₓ;
    T=Float64, θ=sorted_rings(s, ℓₘₐₓ, T),
    plan_fft_flags=FFTW.ESTIMATE, plan_fft_timelimit=Inf,
    inplace=true
)
    @assert length(θ) == ℓₘₐₓ-abs(s)+1 """
        Length of `θ` ($(length(θ))) must equal `ℓₘₐₓ-abs(s)+1` ($(ℓₘₐₓ-abs(s)+1))
    """
    T = eltype(θ)

    θindices = let iθ = cumsum([2j+1 for j ∈ abs(s):ℓₘₐₓ])
        [a:b for (a,b) in eachrow(hcat([1; iθ[begin:end-1].+1], iθ))]
    end

    #s𝐘 = ₛ𝐘(s, ℓₘₐₓ, T, Rθϕ)

    plans = OffsetVector(
        [
            plan_fft!(
                Vector{Complex{T}}(undef, 2j+1),
                flags=plan_fft_flags,
                timelimit=plan_fft_timelimit
            )
            for j ∈ abs(s):ℓₘₐₓ
        ],
        abs(s):ℓₘₐₓ  # This is the range of valid indices
    )

    # Compute ₛ𝐝 as a series of LU-decomposed matrices — one for each m=k value
    ₛΛ = let π=T(π)
        H_rec_coeffs = H_recursion_coefficients(ℓₘₐₓ, T)
        d = dstorage(ℓₘₐₓ, T)
        ₛ𝐝′ = OffsetVector(
            [
                begin
                    r = max(abs(s), abs(m)):ℓₘₐₓ
                    OffsetArray(zeros(T, length(r), length(r)), r, r)
                end
                for m ∈ -ℓₘₐₓ:ℓₘₐₓ
            ],
            -ℓₘₐₓ:ℓₘₐₓ
        )
        for j ∈ abs(s):ℓₘₐₓ
            θⱼ = θ[j-abs(s)+1]
            d!(d, θⱼ, ℓₘₐₓ, H_rec_coeffs)
            for ℓ ∈ abs(s):ℓₘₐₓ
                coefficient = (-1)^s * √(π*T(2ℓ+1))
                for m ∈ -min(j,ℓ):min(j,ℓ)  # Note that k = m
                    ₛ𝐝′[m][j, ℓ] = coefficient * d[WignerDindex(ℓ, m, -s)]
                end
            end
        end
        OffsetVector([LinearAlgebra.lu(ₛ𝐝′[m]) for m ∈ -ℓₘₐₓ:ℓₘₐₓ], -ℓₘₐₓ:ℓₘₐₓ)
    end

    # Pre-allocate the workspace used to solve the linear equations
    workspace = Vector{Complex{T}}(undef, 2ℓₘₐₓ+1)

    ₛf̃₊ₘ = OffsetVector(Vector{Complex{T}}(undef, ℓₘₐₓ+1), 0:ℓₘₐₓ)
    ₛf̃₋ₘ = OffsetVector(Vector{Complex{T}}(undef, ℓₘₐₓ+1), 0:ℓₘₐₓ)
    ₛf̃₊ = OffsetVector(Vector{Complex{T}}(undef, ℓₘₐₓ+1), 0:ℓₘₐₓ)
    ₛf̃₋ = OffsetVector(Vector{Complex{T}}(undef, ℓₘₐₓ+1), 0:ℓₘₐₓ)

    SSHTEKKM{T, inplace}(
        s, ℓₘₐₓ, OffsetVector(θ, abs(s):ℓₘₐₓ), OffsetVector(θindices, abs(s):ℓₘₐₓ),
        #s𝐘,
        plans, ₛΛ, workspace, ₛf̃₊ₘ, ₛf̃₋ₘ, ₛf̃₊, ₛf̃₋
    )
end

function pixels(𝒯::SSHTEKKM{T}) where {T}
    let π = convert(T, π)
        [
            @SVector [𝒯.θ[j], iϕ * 2π / (2j+1)]
            for j ∈ abs(𝒯.s):𝒯.ℓₘₐₓ
            for iϕ ∈ 0:2j
        ]
    end
end

function rotors(𝒯::SSHTEKKM)
    from_spherical_coordinates.(pixels(𝒯))
end

function Base.:*(𝒯::SSHTEKKM, f̃)
    𝒯.ₛ𝐘 * f̃
end

function LinearAlgebra.mul!(f, 𝒯::SSHTEKKM, f̃)
    mul!(f, 𝒯.ₛ𝐘, f̃)
end

function Base.:\(𝒯::SSHTEKKM, f)
    ldiv!(𝒯, copy(f))
end

function Base.:\(𝒯::SSHTEKKM{T, true}, ff̃) where {T}
    ldiv!(𝒯, ff̃)
end

function LinearAlgebra.ldiv!(f̃, 𝒯::SSHTEKKM, f)
    f̃[:] = f
    ldiv!(𝒯, f̃)
end

function LinearAlgebra.ldiv!(𝒯::SSHTEKKM{T}, ff̃) where {T}
    s1 = size(ff̃, 1)
    s2 = Ysize(abs(𝒯.s), 𝒯.ℓₘₐₓ)
    @assert s1==s2 """
        Size of input `ff̃` along first dimension is incorrect for spins `s`=$s and
        `ℓₘₐₓ`=$ℓₘₐₓ; it has size $s1, but should be $s2.
    """
    s = 𝒯.s
    ℓₘₐₓ = 𝒯.ℓₘₐₓ
    mₘₐₓ = ℓₘₐₓ
    θindices = 𝒯.θindices
    ₛΛ = 𝒯.ₛΛ
    ff̃′ = reshape(ff̃, size(ff̃, 1), :)

    @inbounds let π = T(π)
        for ₛf ∈ eachcol(ff̃′)
            # FFT the data in place
            @threads for j ∈ abs(s):ℓₘₐₓ
                jk = θindices[j]
                @views 𝒯.plans[j] * ₛf[jk]
                @. ₛf[jk] *= 2π / (2j+1)
                @views ₛf[jk] .= fftshift(ₛf[jk])
            end

            # TODO: correct this loop for m==0
            for m ∈ ℓₘₐₓ:-1:0  # We do both ±m inside this loop
                Δ = max(abs(s), m)
                # Gather the data for ±m into temporary workspaces
                for j ∈ Δ:ℓₘₐₓ
                    iⱼ₀ = Yindex(j, 0, abs(s))
                    𝒯.ₛf̃₊ₘ[j] = ₛf[iⱼ₀ + m]
                    𝒯.ₛf̃₋ₘ[j] = ₛf[iⱼ₀ - m]
                end
                # Solve for the mode weights from the Fourier components
                # TODO: Use workspace to do linear solves more efficiently
                # TODO: See if I can just do in-place solves
                @views 𝒯.ₛf̃₊[Δ:ℓₘₐₓ] .= ₛΛ[m] \ 𝒯.ₛf̃₊ₘ[Δ:ℓₘₐₓ]
                # ldiv!(ₛΛ[m], 𝒯.ₛf̃₊ₘ[Δ:ℓₘₐₓ])
                @views 𝒯.ₛf̃₋[Δ:ℓₘₐₓ] .= ₛΛ[-m] \ 𝒯.ₛf̃₋ₘ[Δ:ℓₘₐₓ]
                # Scatter the data back into the output
                for ℓ ∈ Δ:ℓₘₐₓ
                    iₗ₀ = Yindex(ℓ, 0, abs(s))
                    ₛf[iₗ₀+mod(ℓ+m, 2ℓ+1)-ℓ] = 𝒯.ₛf̃₊[ℓ]
                    ₛf[iₗ₀+mod(ℓ-m, 2ℓ+1)-ℓ] = 𝒯.ₛf̃₋[ℓ]
                end
                # De-alias remaining Fourier components
                @threads for j′ ∈ m-1:-1:abs(s)
                    α₊ = zero(T)
                    α₋ = zero(T)
                    θ = 𝒯.θ[j′]
                    # TODO: Use λRecursion objects to simplify this
                    cosθ = cos(θ)
                    sin½θ, cos½θ = sincos(θ/2)
                    ₛλₗ₋₁ₘ = zero(T)
                    ₛλₗ₋₁₋ₘ = zero(T)
                    ₛλₗₘ = λ_recursion_initialize(sin½θ, cos½θ, s, Δ, m)
                    ₛλₗ₋ₘ = λ_recursion_initialize(sin½θ, cos½θ, s, Δ, -m)
                    c⁺ₗ₋₁ = zero(T)
                    c⁻ₗ₋₁ = zero(T)
                    for ℓ ∈ Δ:ℓₘₐₓ-1
                        α₊ += 𝒯.ₛf̃₊[ℓ] * ₛλₗₘ
                        α₋ += 𝒯.ₛf̃₋[ℓ] * ₛλₗ₋ₘ
                        # Take another step in the λ recursions
                        c⁺ₗ₊₁, c⁺ₗ = λ_recursion_coefficients(cosθ, s, ℓ, m)
                        c⁻ₗ₊₁, c⁻ₗ = λ_recursion_coefficients(cosθ, s, ℓ, -m)
                        ₛλₗ₊₁ₘ = if ℓ == 0
                            √(3/4π) * cosθ
                        else
                            (c⁺ₗ * ₛλₗₘ + c⁺ₗ₋₁ * ₛλₗ₋₁ₘ) / c⁺ₗ₊₁
                        end
                        ₛλₗ₊₁₋ₘ = if ℓ == 0
                            √(3/4π) * cosθ
                        else
                            (c⁻ₗ * ₛλₗ₋ₘ + c⁻ₗ₋₁ * ₛλₗ₋₁₋ₘ) / c⁻ₗ₊₁
                        end
                        ₛλₗ₋₁ₘ = ₛλₗₘ
                        ₛλₗ₋₁₋ₘ = ₛλₗ₋ₘ
                        ₛλₗₘ = ₛλₗ₊₁ₘ
                        ₛλₗ₋ₘ = ₛλₗ₊₁₋ₘ
                        c⁺ₗ₋₁ = -c⁺ₗ₊₁ * √((2ℓ+1)/T(2ℓ+3))
                        c⁻ₗ₋₁ = -c⁻ₗ₊₁ * √((2ℓ+1)/T(2ℓ+3))
                    end
                    # Do the above once more for ℓ==ℓₘₐₓ, and skip recursion
                    α₊ += 𝒯.ₛf̃₊[ℓₘₐₓ] * ₛλₗₘ
                    α₋ += 𝒯.ₛf̃₋[ℓₘₐₓ] * ₛλₗ₋ₘ
                    # Finally, de-alias
                    iⱼ′₀ = Yindex(j′, 0, abs(s))
                    ₛf[iⱼ′₀ + mod(j′+m, 2j′+1)-j′] -= 2π * α₊
                    ₛf[iⱼ′₀ + mod(j′-m, 2j′+1)-j′] -= 2π * α₋
                end  # j′
            end  # m
        end # ₛf
    end  # π
    ff̃
end
