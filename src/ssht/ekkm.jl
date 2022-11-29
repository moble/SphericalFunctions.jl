"""Storage for EKKM spin-spherical-harmonic transform

The EKKM algorithm was described in [this paper](https://arxiv.org/abs/1809.01321), and allows for
the minimal number of function samples.

"""
struct SSHTEKKM{Inplace} <: SSHT
    """Spin weight"""
    s::Integer

    """Highest ℓ value present in the data"""
    ℓₘₐₓ::Integer

    """OffsetVector of colatitudes of sampling rings

    These are the coordinates of the rings, where ring `j` contains 2j+1
    equally spaced points.  The index ranges over j ∈ abs(s):ℓₘₐₓ.
    """
    θ::OffsetVector

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
    ₛ𝐝::OffsetVector

    """Preallocated storage for solving linear equations"""
    workspace::Vector
end


function SSHTEKKM(
    s, ℓₘₐₓ;
    T=Float64, θ=sorted_rings(s, ℓₘₐₓ, T),
    plan_fft_flags=FFTW.ESTIMATE, plan_fft_timelimit=Inf,
    inplace=true
)
    @assert length(θ) == ℓₘₐₓ-abs(s)+1 "Length of `θ` ($(length(θ))) must equal `ℓₘₐₓ-abs(s)+1` ($(ℓₘₐₓ-abs(s)+1))"
    T = eltype(θ)

    plans = OffsetVector(
        [
            plan_fft!(Vector{Complex{T}}(undef, 2j+1), flags=plan_fft_flags, timelimit=plan_fft_timelimit)
            for j ∈ abs(s):ℓₘₐₓ
        ],
        abs(s):ℓₘₐₓ
    )

    # Compute ₛ𝐝 as a series of LU-decomposed matrices — one for each m=k value
    ₛ𝐝 = let
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
                coefficient = √T(2ℓ+1)
                for m ∈ -min(j,ℓ):min(j,ℓ)  # Note that k = m
                    ₛ𝐝′[m][j, ℓ] = coefficient * d[WignerDindex(ℓ, m, s)]
                end
            end
        end
        OffsetVector([LinearAlgebra.lu(ₛ𝐝′[m]) for m ∈ -ℓₘₐₓ:ℓₘₐₓ], -ℓₘₐₓ:ℓₘₐₓ)
    end
    
    # Pre-allocate the workspace used to solve the linear equations
    workspace = Vector{Complex{T}}(undef, 2ℓₘₐₓ+1)

    SSHTEKKM{inplace}(s, ℓₘₐₓ, OffsetVector(θ, abs(s):ℓₘₐₓ), plans, ₛ𝐝, workspace)
end


function Base.:\(ssht::SSHTEKKM, f)
    ldiv!(ssht, copy(f))
end

function Base.:\(ssht::SSHTEKKM{true}, f)
    ldiv!(ssht, f)
end

function LinearAlgebra.ldiv!(Y, ssht::SSHTEKKM, f)
    Y[:] = f
    ldiv!(ssht, Y)
end

function LinearAlgebra.ldiv!(ssht::SSHTEKKM, f)
    # s, ℓₘₐₓ, plans, ₛ𝐝, workspace
    i₁ = firstindex(f)
    for j ∈ abs(ssht.s):ssht.ℓₘₐₓ
        i₂ = i₁ + 2j+1
        @views (plans[j] * f[i₁:i₂])  # performs FFT in place
        i₁ = i₂+1
    end
    for k ∈ ssht.ℓₘₐₓ:-1:0
        w = @view ssht.workspace[begin:begin+ssht.ℓₘₐₓ-abs(k)]

        # Copy all harmonics into workspace
        @error "Not implemented"

        # Solve for mode weights
        ldiv!(ssht.ₛ𝐝[k], w)

        # Copy all mode weights back into f
        @error "Not implemented"

        # De-alias lower |k| harmonics
        for ℓ ∈ abs(k):-1:abs(s)
            for n ∈ cld(-ℓ-k, 2j+1):fld(ℓ-k, 2j+1)
                m = k + n * (2j+1)
                @error "Not implemented"
            end
        end
    end
end
