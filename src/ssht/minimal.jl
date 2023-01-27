
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
    bplans::OffsetVector

    """OffsetVector of matrices used to compute mode weights

    There is one matrix for each value of m = k ∈ -ℓₘₐₓ:ℓₘₐₓ.  Indexing should
    look like ₛΛ[m][j, ℓ], where the rows in each matrix correspond to
    j ∈ abs(s):abs(m) and the columns correspond to ℓ ∈ max(abs(s),abs(m)):ℓₘₐₓ.
    In particular, note that j ≤ |m| in these matrices; the opposite is true for
    the `luₛΛ` variable.
    """
    ₛΛ::OffsetVector

    """OffsetVector of LU-decomposed ₛΛ matrices, used to solve for mode weights

    There is one matrix for each value of m = k ∈ -ℓₘₐₓ:ℓₘₐₓ.  Indexing should
    look like ₛΛ[m][j, ℓ], where the rows in each matrix correspond to
    j ∈ max(abs(s),abs(m)):ℓₘₐₓ and the columns correspond to ℓ over the same
    range as j.  In particular, note that j ≥ |m| in these matrices; the opposite
    is true for the `ₛΛ` variable.
    """
    luₛΛ::OffsetVector

    """Preallocated storage for all FT modes with a given ``m``"""
    ₛfₘ::OffsetVector

    """Preallocated storage for all SH modes with a given positive ``m`` value"""
    ₛf̃ₘ::OffsetVector

    """Preallocated storage for function values / Fourier modes on a given ring"""
    ₛf̃ⱼ
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

    ₛfₘ = OffsetVector(Vector{Complex{T}}(undef, ℓₘₐₓ+1), 0:ℓₘₐₓ)
    ₛf̃ₘ = OffsetVector(Vector{Complex{T}}(undef, ℓₘₐₓ+1), 0:ℓₘₐₓ)
    ₛf̃ⱼ = OffsetVector([Vector{Complex{T}}(undef, 2j+1) for j ∈ abs(s):ℓₘₐₓ], abs(s):ℓₘₐₓ)

    plans = OffsetVector(
        if T ∈ [Float64, Float32]  # Only supported types in FFTW
            [
                plan_fft!(
                    ₛf̃ⱼ[j],
                    flags=plan_fft_flags,
                    timelimit=plan_fft_timelimit
                )
                for j ∈ abs(s):ℓₘₐₓ
            ]
        else
            [plan_fft!(ₛf̃ⱼ[j]) for j ∈ abs(s):ℓₘₐₓ]
        end,
        abs(s):ℓₘₐₓ  # This is the range of valid indices
    )
    bplans = OffsetVector(
        if T ∈ [Float64, Float32]  # Only supported types in FFTW
            [
                plan_bfft!(
                    Vector{Complex{T}}(undef, 2j+1),
                    flags=plan_fft_flags,
                    timelimit=plan_fft_timelimit
                )
                for j ∈ abs(s):ℓₘₐₓ
            ]
        else
            [plan_bfft!(Vector{Complex{T}}(undef, 2j+1)) for j ∈ abs(s):ℓₘₐₓ]
        end,
        abs(s):ℓₘₐₓ  # This is the range of valid indices
    )

    # Compute ₛ𝐝 as a series of LU-decomposed matrices — one for each m=k value
    ₛΛ, luₛΛ = let π=T(π)
        H_rec_coeffs = H_recursion_coefficients(ℓₘₐₓ, T)
        d = dstorage(ℓₘₐₓ, T)
        ₛΛ = OffsetVector(
            [
                let J = abs(s):abs(m), Δ = max(abs(s), abs(m)):ℓₘₐₓ
                    OffsetArray(zeros(T, length(J), length(Δ)), J, Δ)
                end
                for m ∈ -ℓₘₐₓ:ℓₘₐₓ
            ],
            -ℓₘₐₓ:ℓₘₐₓ
        )
        luₛΛ = OffsetVector(
            [
                let Δ = max(abs(s), abs(m)):ℓₘₐₓ
                    OffsetArray(zeros(T, length(Δ), length(Δ)), Δ, Δ)
                end
                for m ∈ -ℓₘₐₓ:ℓₘₐₓ
            ],
            -ℓₘₐₓ:ℓₘₐₓ
        )
        for j ∈ abs(s):ℓₘₐₓ
            θⱼ = θ[j-abs(s)+1]
            d!(d, θⱼ, ℓₘₐₓ, H_rec_coeffs)
            for ℓ ∈ abs(s):ℓₘₐₓ
                coefficient = (-1)^s * √((2ℓ+1)/4π)
                for m ∈ [-ℓ:-j... ; j:ℓ...]
                    ₛΛ[m][j, ℓ] = coefficient * d[WignerDindex(ℓ, m, -s)]
                end
                coefficient *= 2π
                for m ∈ -min(j,ℓ):min(j,ℓ)
                    luₛΛ[m][j, ℓ] = coefficient * d[WignerDindex(ℓ, m, -s)]
                end
            end
        end
        ₛΛ, OffsetVector([LinearAlgebra.lu(luₛΛ[m]) for m ∈ -ℓₘₐₓ:ℓₘₐₓ], -ℓₘₐₓ:ℓₘₐₓ)
    end

    SSHTMinimal{T, inplace}(
        s, ℓₘₐₓ, OffsetVector(θ, abs(s):ℓₘₐₓ), OffsetVector(θindices, abs(s):ℓₘₐₓ),
        plans, bplans, ₛΛ, luₛΛ, ₛfₘ, ₛf̃ₘ, ₛf̃ⱼ
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
    mul!(𝒯::SSHTMinimal, copy(f̃))
end

function Base.:*(𝒯::SSHTMinimal{T, true}, f̃) where {T}
    mul!(𝒯::SSHTMinimal, f̃)
end

function LinearAlgebra.mul!(f, 𝒯::SSHTMinimal, f̃)
    ff̃ .= f̃
    mul!(𝒯.ₛ𝐘, ff̃)
end

function LinearAlgebra.mul!(𝒯::SSHTMinimal{T}, ff̃) where {T}
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
        for ₛf̃ ∈ eachcol(ff̃′)
            for m ∈ AlternatingCountup(ℓₘₐₓ)  # Iterate over +m, then -m, up from m=0
                Δ = max(abs(s), abs(m))

                # Iterate over rings, combining contributions for this `m` value
                @threads for j ∈ Δ:ℓₘₐₓ
                    # We will accumulate into 𝒯.ₛfₘ, and write it out at the end of the loop
                    𝒯.ₛfₘ[j] = false

                    # Direct (non-aliased) contributions from m′ == m
                    λ = λiterator(𝒯.θ[j], s, m)
                    for (ℓ, ₛλₗₘ) ∈ zip(Δ:ℓₘₐₓ, λ)
                        𝒯.ₛfₘ[j] += ₛf̃[Yindex(ℓ, m, abs(s))] * ₛλₗₘ
                    end

                    # Aliased contributions from |m′| > j > |m|
                    for ℓ′ ∈ j:ℓₘₐₓ
                        for n ∈ cld(-ℓ′-m, 2j+1):fld(ℓ′-m, 2j+1)
                            m′ = m + n*(2j+1)
                            if abs(m′) > j
                                ₛλₗ′ₘ′ = 𝒯.ₛΛ[m′][j,ℓ′]
                                𝒯.ₛfₘ[j] += ₛf̃[Yindex(ℓ′, m′, abs(s))] * ₛλₗ′ₘ′
                            end
                        end
                    end

                end  # j

                # Distribute the data back into the output
                @threads for j ∈ Δ:ℓₘₐₓ
                    ₛf̃[Yindex(j, m, abs(s))] = 𝒯.ₛfₘ[j]
                end

            end  # m

            # Iterate over rings, doing Fourier decompositions on each
            @threads for j ∈ abs(s):ℓₘₐₓ
                jk = 𝒯.θindices[j]
                @views ifftshift!(𝒯.ₛf̃ⱼ[j], ₛf̃[jk]) # Cycle order of modes in place to match order of FFT elements
                @views 𝒯.bplans[j] * 𝒯.ₛf̃ⱼ[j] # Perform in-place BFFT
                @. ₛf̃[jk] = 𝒯.ₛf̃ⱼ[j]  # Copy data back into main array
            end

        end # ₛf̃
    end  # π
    ff̃
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
            # Iterate over rings, doing Fourier decompositions on each
            for j ∈ abs(s):ℓₘₐₓ
                jk = 𝒯.θindices[j]
                @. 𝒯.ₛf̃ⱼ[j] = ₛf[jk] * 2π / (2j+1) # Copy data from main array and normalize
                @views 𝒯.plans[j] * 𝒯.ₛf̃ⱼ[j] # Perform in-place IFFT
                @views fftshift!(ₛf[jk], 𝒯.ₛf̃ⱼ[j]) # Cycle order of modes in place to match order of FFT elements
            end

            for m ∈ AlternatingCountdown(ℓₘₐₓ)
                Δ = max(abs(s), abs(m))

                # Gather the `m` data from each ring into a temporary workspace
                @threads for j ∈ Δ:ℓₘₐₓ
                    𝒯.ₛfₘ[j] = ₛf[Yindex(j, m, abs(s))]
                end

                # Solve for the mode weights from the Fourier components
                @views ldiv!(𝒯.ₛf̃ₘ.parent[Δ+1:ℓₘₐₓ+1], 𝒯.luₛΛ[m], 𝒯.ₛfₘ.parent[Δ+1:ℓₘₐₓ+1])

                # Distribute the data back into the output
                @threads for ℓ ∈ Δ:ℓₘₐₓ
                    ₛf[Yindex(ℓ, m, abs(s))] = 𝒯.ₛf̃ₘ[ℓ]
                end

                # De-alias Fourier components from rings with values of j < Δ
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
