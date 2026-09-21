"""
    SSHTMinimal(s, ℓₘₐₓ; T=Float64, θ=sorted_rings(s, ℓₘₐₓ, T), plan_fft_flags=FFTW.ESTIMATE, plan_fft_timelimit=Inf, inplace=true)

Construct an ``s``-SHT object that uses the optimal-dimensionality algorithm described by
[Elahi et al.](@cite Elahi_2018), which samples the function at exactly as many points as
there are modes.  This may also be achieved by calling the main [`SSHT`](@ref) function with
the same keywords, along with `method="Minimal"`.

The spin-weighted spherical harmonics are evaluated on a series of "rings" at constant
colatitude, where ring ``j`` (for ``j ∈ |s|:ℓₘₐₓ``) contains ``2j+1`` equally spaced points.
Their colatitudes are given by the `θ` keyword argument, which defaults to
[`sorted_rings(s, ℓₘₐₓ, T)`](@ref sorted_rings); the first element of `θ` is the colatitude
of the smallest ring (with ``2|s|+1`` points), and the last is that of the largest ring (with
``2ℓₘₐₓ+1`` points).  See [`pixels`](@ref) and [`rotors`](@ref) for the sample points.

Whenever `T` is either `Float64` or `Float32`, the keyword arguments `plan_fft_flags` and
`plan_fft_timelimit` may also be useful for obtaining more efficient FFTs.  They default to
`FFTW.ESTIMATE` and `Inf`, respectively, and are passed to
[`AbstractFFTs.plan_fft!`](https://juliamath.github.io/AbstractFFTs.jl/stable/api/#AbstractFFTs.plan_fft).

Because this algorithm achieves optimal dimensionality, the transformation is performed in
place by default: `𝒯 * f̃` overwrites `f̃` with the function values (and returns it), and
`𝒯 \\ f` overwrites `f`.  If this is not desired, pass the keyword argument `inplace=false`,
which makes those operations work on a copy of the input.  See [`SSHT`](@ref).

The values ``{}_sλ_{ℓ,m}(θ_j)`` needed by the algorithm are precomputed at construction (with
one batched [`sλlmCalculator`](@ref)) and stored, which takes ``O(ℓₘₐₓ^3)`` memory.

This method is defined only for integer spin weights.  Its bookkeeping — rings of ``2j+1``
points indexed by ``j``, and the aliasing of ``m`` into rings too small to hold it — is
written for integer indices throughout, and has not been extended; a half-integer spin weight
is refused with a message naming the two methods, `"RS"` and `"Matrix"`, that do accept one.
"""
struct SSHTMinimal{T<:Real, Inplace, P, BP} <: SSHT{T}
    s::Int
    ℓₘₐₓ::Int
    θ::OffsetVector{T, Vector{T}}  # colatitude of ring j, for j ∈ abs(s):ℓₘₐₓ
    θindices::OffsetVector{UnitRange{Int}, Vector{UnitRange{Int}}}  # pixel indices of ring j
    plans::OffsetVector{P, Vector{P}}  # forward in-place FFT plans, ring j
    bplans::OffsetVector{BP, Vector{BP}}  # backward in-place FFT plans, ring j
    # ₛΛ[m][j, ℓ] = ₛλₗₘ(θⱼ) for j ∈ abs(s):abs(m) and ℓ ∈ max(abs(s),abs(m)):ℓₘₐₓ: the
    # values that alias into rings too small to hold mode m.
    ₛΛ::OffsetVector{OffsetMatrix{T, Matrix{T}}, Vector{OffsetMatrix{T, Matrix{T}}}}
    # Λ[m][j, ℓ] = ₛλₗₘ(θⱼ) for j, ℓ ∈ max(abs(s),abs(m)):ℓₘₐₓ: the direct contributions.
    Λ::OffsetVector{OffsetMatrix{T, Matrix{T}}, Vector{OffsetMatrix{T, Matrix{T}}}}
    # LU decompositions of 2π Λ[m], used to solve for the mode weights
    luΛ::OffsetVector{LinearAlgebra.LU{T, Matrix{T}, Vector{Int}}, Vector{LinearAlgebra.LU{T, Matrix{T}, Vector{Int}}}}
    ₛfₘ::OffsetVector{Complex{T}, Vector{Complex{T}}}  # workspace: Fourier modes with a given m, indexed by ring j
    ₛf̃ₘ::OffsetVector{Complex{T}, Vector{Complex{T}}}  # workspace: mode weights with a given m, indexed by ℓ
    ₛf̃ⱼ::OffsetVector{Vector{Complex{T}}, Vector{Vector{Complex{T}}}}  # workspace: values on ring j
end

# Iterate m as 0, 1, -1, 2, -2, …, ℓₘₐₓ, -ℓₘₐₓ (and the reverse)
function alternating_countup(ℓₘₐₓ::Integer)
    ms = Int[0]
    for m ∈ 1:ℓₘₐₓ
        push!(ms, m)
        push!(ms, -m)
    end
    ms
end
alternating_countdown(ℓₘₐₓ::Integer) = reverse(alternating_countup(ℓₘₐₓ))

# The public constructor is the boundary; the indices arrive at the worker as `Int`s, or as
# `HalfOddInteger`s, for which the worker is the refusal described in the docstring.
function SSHTMinimal(s::IndexSpelling, ℓₘₐₓ::IndexSpelling; T::Type{TT}=Float64, kwargs...) where {TT}
    SSHTMinimal(transform_indices(s, ℓₘₐₓ)..., TT; kwargs...)
end
function SSHTMinimal(s::HalfOddInteger, ℓₘₐₓ::HalfOddInteger, ::Type; kwargs...)
    error(
        "The \"Minimal\" s-SHT method is defined only for integer spin weights, but s=$s and "
        * "ℓₘₐₓ=$ℓₘₐₓ are half-odd-integers.  The \"RS\" method (the default) and the "
        * "\"Matrix\" method both accept half-integer indices."
    )
end
function SSHTMinimal(
    s::Int, ℓₘₐₓ::Int, ::Type{TT};
    θ=sorted_rings(s, ℓₘₐₓ, TT),
    plan_fft_flags=FFTW.ESTIMATE, plan_fft_timelimit=Inf,
    inplace=true
) where {TT}
    if abs(s) > ℓₘₐₓ
        error("|s|=$(abs(s)) exceeds ℓₘₐₓ=$ℓₘₐₓ; there are no such modes.")
    end
    if length(θ) != ℓₘₐₓ - abs(s) + 1
        error("Length of θ ($(length(θ))) must equal ℓₘₐₓ-abs(s)+1 ($(ℓₘₐₓ-abs(s)+1)).")
    end
    θ = Vector{TT}(θ)
    J = abs(s):ℓₘₐₓ  # ring indices

    θindices = let stops = cumsum([2j+1 for j ∈ J])
        [(stop - (2j+1) + 1):stop for (j, stop) ∈ zip(J, stops)]
    end

    ₛfₘ = OffsetVector(Vector{Complex{TT}}(undef, ℓₘₐₓ+1), 0:ℓₘₐₓ)
    ₛf̃ₘ = OffsetVector(Vector{Complex{TT}}(undef, ℓₘₐₓ+1), 0:ℓₘₐₓ)
    ₛf̃ⱼ = OffsetVector([Vector{Complex{TT}}(undef, 2j+1) for j ∈ J], J)

    plans, bplans = if TT ∈ (Float64, Float32)  # Only supported types in FFTW
        (
            [plan_fft!(ₛf̃ⱼ[j]; flags=plan_fft_flags, timelimit=plan_fft_timelimit) for j ∈ J],
            [plan_bfft!(ₛf̃ⱼ[j]; flags=plan_fft_flags, timelimit=plan_fft_timelimit) for j ∈ J],
        )
    else
        ([plan_fft!(ₛf̃ⱼ[j]) for j ∈ J], [plan_bfft!(ₛf̃ⱼ[j]) for j ∈ J])
    end

    # Tables of ₛλₗₘ(θⱼ) for every ring, computed with one batched calculator in angle mode
    ₛΛ = OffsetVector(
        [
            let Jm = abs(s):abs(m), L = max(abs(s), abs(m)):ℓₘₐₓ
                OffsetArray(zeros(TT, length(Jm), length(L)), Jm, L)
            end
            for m ∈ -ℓₘₐₓ:ℓₘₐₓ
        ],
        -ℓₘₐₓ:ℓₘₐₓ
    )
    Λ = OffsetVector(
        [
            let L = max(abs(s), abs(m)):ℓₘₐₓ
                OffsetArray(zeros(TT, length(L), length(L)), L, L)
            end
            for m ∈ -ℓₘₐₓ:ℓₘₐₓ
        ],
        -ℓₘₐₓ:ℓₘₐₓ
    )
    λ = sλlmCalculator(θ, ℓₘₐₓ, s)  # θ is already a Vector{TT}, which fixes the type
    iₛ = spin_index(λ, s)
    Λℓ = λ.Yˡ  # [ring index 1:length(J), spin, m+ℓ+1]
    for ℓ ∈ J
        recurrence!(λ, ℓ)
        for m ∈ -ℓ:ℓ
            for (jᵢ, j) ∈ enumerate(J)
                value = Λℓ[jᵢ, iₛ, m + ℓ + 1]
                if abs(m) > j
                    ₛΛ[m][j, ℓ] = value
                else
                    Λ[m][j, ℓ] = value
                end
            end
        end
    end
    luΛ = OffsetVector(
        [LinearAlgebra.lu(2TT(π) * parent(Λ[m])) for m ∈ -ℓₘₐₓ:ℓₘₐₓ], -ℓₘₐₓ:ℓₘₐₓ
    )

    SSHTMinimal{TT, inplace, eltype(plans), eltype(bplans)}(
        s, ℓₘₐₓ, OffsetVector(θ, J), OffsetVector(θindices, J),
        OffsetVector(plans, J), OffsetVector(bplans, J),
        ₛΛ, Λ, luΛ, ₛfₘ, ₛf̃ₘ, ₛf̃ⱼ
    )
end

function pixels(𝒯::SSHTMinimal{T}) where {T}
    let π = T(π)
        [
            @SVector [𝒯.θ[j], iϕ * 2π / (2j+1)]
            for j ∈ abs(𝒯.s):𝒯.ℓₘₐₓ
            for iϕ ∈ 0:2j
        ]
    end
end
rotors(𝒯::SSHTMinimal) = from_spherical_coordinates.(pixels(𝒯))
npixels(𝒯::SSHTMinimal) = nmodes(𝒯)

function Base.:*(𝒯::SSHTMinimal, f̃)
    check_modes(𝒯, f̃)
    mul!(𝒯, copy(array_view(f̃)))
end
function Base.:*(𝒯::SSHTMinimal{T, true}, f̃) where {T}
    check_modes(𝒯, f̃)
    mul!(𝒯, array_view(f̃))
    f̃
end
function LinearAlgebra.mul!(f, 𝒯::SSHTMinimal, f̃)
    check_modes(𝒯, f̃)
    check_pixels(𝒯, f)
    f .= array_view(f̃)
    mul!(𝒯, f)
end

# Synthesis in place: the mode weights in `ff̃` are replaced by the function values
function LinearAlgebra.mul!(𝒯::SSHTMinimal{T}, ff̃) where {T}
    check_modes(𝒯, ff̃)
    s, ℓₘₐₓ = 𝒯.s, 𝒯.ℓₘₐₓ
    ff̃′ = reshape(array_view(ff̃), size(ff̃, 1), :)

    @inbounds for ₛf̃ ∈ eachcol(ff̃′)
        for m ∈ alternating_countup(ℓₘₐₓ)  # Iterate over +m, then -m, up from m=0
            Δ = max(abs(s), abs(m))
            Λm = 𝒯.Λ[m]

            # Iterate over rings, combining contributions for this `m` value
            @threads for j ∈ Δ:ℓₘₐₓ
                # We will accumulate into 𝒯.ₛfₘ, and write it out at the end of the loop
                acc = zero(Complex{T})

                # Direct (non-aliased) contributions from m′ == m
                for ℓ ∈ Δ:ℓₘₐₓ
                    acc += ₛf̃[Yindex(ℓ, m, abs(s))] * Λm[j, ℓ]
                end

                # Aliased contributions from |m′| > j > |m|
                for ℓ′ ∈ j:ℓₘₐₓ
                    for n ∈ cld(-ℓ′-m, 2j+1):fld(ℓ′-m, 2j+1)
                        m′ = m + n*(2j+1)
                        if abs(m′) > j
                            acc += ₛf̃[Yindex(ℓ′, m′, abs(s))] * 𝒯.ₛΛ[m′][j, ℓ′]
                        end
                    end
                end

                𝒯.ₛfₘ[j] = acc
            end  # j

            # Distribute the data back into the output
            @threads for j ∈ Δ:ℓₘₐₓ
                ₛf̃[Yindex(j, m, abs(s))] = 𝒯.ₛfₘ[j]
            end
        end  # m

        # Iterate over rings, doing Fourier synthesis on each
        @threads for j ∈ abs(s):ℓₘₐₓ
            jk = 𝒯.θindices[j]
            @views ifftshift!(𝒯.ₛf̃ⱼ[j], ₛf̃[jk])  # Reorder modes to match FFT element order
            𝒯.bplans[j] * 𝒯.ₛf̃ⱼ[j]  # Perform in-place unnormalized inverse FFT
            @. ₛf̃[jk] = 𝒯.ₛf̃ⱼ[j]  # Copy data back into main array
        end
    end  # ₛf̃
    ff̃
end

function Base.:\(𝒯::SSHTMinimal, f)
    check_pixels(𝒯, f)
    f̃ = ldiv!(𝒯, copy(f))
    ndims(f) == 1 ? ModeWeights(f̃, 𝒯.s, abs(𝒯.s), 𝒯.ℓₘₐₓ) : f̃
end
function Base.:\(𝒯::SSHTMinimal{T, true}, ff̃) where {T}
    check_pixels(𝒯, ff̃)
    ldiv!(𝒯, array_view(ff̃))
    ff̃
end
function LinearAlgebra.ldiv!(f̃, 𝒯::SSHTMinimal, f)
    check_modes(𝒯, f̃)
    check_pixels(𝒯, f)
    array_view(f̃) .= f
    ldiv!(𝒯, array_view(f̃))
    f̃
end

# Analysis in place: the function values in `ff̃` are replaced by the mode weights
function LinearAlgebra.ldiv!(𝒯::SSHTMinimal{T}, ff̃) where {T}
    check_pixels(𝒯, ff̃)
    s, ℓₘₐₓ = 𝒯.s, 𝒯.ℓₘₐₓ
    ff̃′ = reshape(array_view(ff̃), size(ff̃, 1), :)

    @inbounds let π = T(π)
        for ₛf ∈ eachcol(ff̃′)
            # Iterate over rings, doing Fourier analysis on each
            for j ∈ abs(s):ℓₘₐₓ
                jk = 𝒯.θindices[j]
                @. 𝒯.ₛf̃ⱼ[j] = ₛf[jk] * 2π / (2j+1)  # Copy data from main array and normalize
                𝒯.plans[j] * 𝒯.ₛf̃ⱼ[j]  # Perform in-place FFT
                @views fftshift!(ₛf[jk], 𝒯.ₛf̃ⱼ[j])  # Reorder to m ∈ -j:j
            end

            for m ∈ alternating_countdown(ℓₘₐₓ)
                Δ = max(abs(s), abs(m))

                # Gather the `m` data from each ring into a temporary workspace
                @threads for j ∈ Δ:ℓₘₐₓ
                    𝒯.ₛfₘ[j] = ₛf[Yindex(j, m, abs(s))]
                end

                # Solve for the mode weights from the Fourier components
                @views ldiv!(
                    parent(𝒯.ₛf̃ₘ)[Δ+1:ℓₘₐₓ+1], 𝒯.luΛ[m], parent(𝒯.ₛfₘ)[Δ+1:ℓₘₐₓ+1]
                )

                # Distribute the data back into the output
                @threads for ℓ ∈ Δ:ℓₘₐₓ
                    ₛf[Yindex(ℓ, m, abs(s))] = 𝒯.ₛf̃ₘ[ℓ]
                end

                # De-alias Fourier components from rings with values of j < Δ
                @threads for j′ ∈ abs(s):abs(m)-1
                    m′ = mod(j′+m, 2j′+1)-j′  # `m` aliases into `(j′, m′)`
                    α = 2π * sum(𝒯.ₛf̃ₘ[ℓ] * 𝒯.ₛΛ[m][j′, ℓ] for ℓ ∈ Δ:ℓₘₐₓ)
                    ₛf[Yindex(j′, m′, abs(s))] -= α
                end  # j′
            end  # m
        end  # ₛf
    end  # π
    ff̃
end
