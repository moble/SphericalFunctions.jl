
"""Storage for Minimal spin-spherical-harmonic transform

The Minimal algorithm was described in [this paper](https://arxiv.org/abs/1809.01321), and allows for
the minimal number of function samples.

"""
struct SSHTMinimal{T<:Real, Inplace} <: SSHT{T}
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
    ₛΛ::OffsetVector

    """Preallocated storage for solving linear equations"""
    workspace::Vector

    """Preallocated storage for all FT modes with a given ``m``"""
    ₛfₘ::OffsetVector

    """Preallocated storage for all SH modes with a given positive ``m`` value"""
    ₛf̃ₘ::OffsetVector
end

function SSHTMinimal(
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

    ₛfₘ = OffsetVector(Vector{Complex{T}}(undef, ℓₘₐₓ+1), 0:ℓₘₐₓ)
    ₛf̃ₘ = OffsetVector(Vector{Complex{T}}(undef, ℓₘₐₓ+1), 0:ℓₘₐₓ)

    SSHTMinimal{T, inplace}(
        s, ℓₘₐₓ, OffsetVector(θ, abs(s):ℓₘₐₓ), OffsetVector(θindices, abs(s):ℓₘₐₓ),
        plans, ₛΛ, workspace, ₛfₘ, ₛf̃ₘ
    )
end

function pixels(𝒯::SSHTMinimal{T}) where {T}
    let π = convert(T, π)
        [
            @SVector [𝒯.θ[j], iϕ * 2π / (2j+1)]
            for j ∈ abs(𝒯.s):𝒯.ℓₘₐₓ
            for iϕ ∈ 0:2j
        ]
    end
end

function rotors(𝒯::SSHTMinimal)
    from_spherical_coordinates.(pixels(𝒯))
end

function Base.:*(𝒯::SSHTMinimal, f̃)
    𝒯.ₛ𝐘 * f̃
end

function LinearAlgebra.mul!(f, 𝒯::SSHTMinimal, f̃)
    mul!(f, 𝒯.ₛ𝐘, f̃)
end

function Base.:\(𝒯::SSHTMinimal, f)
    ldiv!(𝒯, copy(f))
end

function Base.:\(𝒯::SSHTMinimal{T, true}, ff̃) where {T}
    ldiv!(𝒯, ff̃)
end

function LinearAlgebra.ldiv!(f̃, 𝒯::SSHTMinimal, f)
    f̃[:] = f
    ldiv!(𝒯, f̃)
end

function LinearAlgebra.ldiv!(𝒯::SSHTMinimal{T}, ff̃) where {T}
    s1 = size(ff̃, 1)
    s2 = Ysize(abs(𝒯.s), 𝒯.ℓₘₐₓ)
    @assert s1==s2 """
        Size of input `ff̃` along first dimension is incorrect for spins `s`=$s and
        `ℓₘₐₓ`=$ℓₘₐₓ; it has size $s1, but should be $s2.
    """
    s = 𝒯.s
    ℓₘₐₓ = 𝒯.ℓₘₐₓ
    ff̃′ = reshape(ff̃, size(ff̃, 1), :)

    @inbounds let π = T(π)
        for ₛf ∈ eachcol(ff̃′)
            # FFT the data in place
            @threads for j ∈ abs(s):ℓₘₐₓ
                jk = 𝒯.θindices[j]
                @views 𝒯.plans[j] * ₛf[jk]
                @. ₛf[jk] *= 2π / (2j+1)
                @views ₛf[jk] .= fftshift(ₛf[jk])
            end

            for m ∈ AlternatingCountdown(ℓₘₐₓ)
                Δ = max(abs(s), abs(m))

                # Gather the `m` data from each ring into a temporary workspace
                for j ∈ Δ:ℓₘₐₓ
                    𝒯.ₛfₘ[j] = ₛf[Yindex(j, m, abs(s))]
                end

                # Solve for the mode weights from the Fourier components
                # TODO: Use workspace to do linear solves more efficiently
                # TODO: See if I can just do in-place solves
                @views 𝒯.ₛf̃ₘ[Δ:ℓₘₐₓ] .= 𝒯.ₛΛ[m] \ 𝒯.ₛfₘ[Δ:ℓₘₐₓ]
                # ldiv!(ₛΛ[m], 𝒯.ₛf̃ₘ[Δ:ℓₘₐₓ])

                # Scatter the data back into the output
                for ℓ ∈ Δ:ℓₘₐₓ
                    ₛf[Yindex(ℓ, m, abs(s))] = 𝒯.ₛf̃ₘ[ℓ]
                end

                # De-alias remaining Fourier components
                @threads for j′ ∈ abs(s):abs(m)-1
                    m′ = mod(j′+m, 2j′+1)-j′  # `m` aliases into `(j′, m′)`
                    α = 2π * sum(
                        𝒯.ₛf̃ₘ[ℓ] * ₛλₗₘ
                        for (ℓ, ₛλₗₘ) ∈ zip(Δ:ℓₘₐₓ, λiterator(𝒯.θ[j′], s, m))
                    )
                    ₛf[Yindex(j′, m′, abs(s))] -= α
                end  # j′
            end  # m
        end # ₛf
    end  # π
    ff̃
end
