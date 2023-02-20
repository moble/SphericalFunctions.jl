# TODO: Add enough `G` storage to use on multiple threads / SIMD registers (`plan`s are thread-safe)
# TODO: Reorganize to match actual R&S algorithm recommendations
# TODO: Allow SSHTRS to take a single Int, rather than a full vector for Nϕ
# TODO: Optimize: Check allocations, fuse ±m loops, etc.


"""Storage for  spin-spherical-harmonic transform

The algorithm was described in [this paper by Reinecke and
Seljebotn](https://arxiv.org/abs/1303.4945).

"""
struct SSHTRS{T<:Real} <: SSHT{T}
    """Spin weight"""
    s

    """Highest ℓ value present in the data"""
    ℓₘₐₓ

    """Locations of rings along the colatitude (θ) coordinate"""
    θ

    """Quadrature weights for the given θ values"""
    quadrature_weight

    """Number of points along the azimuthal (ϕ) coordinate in each ring"""
    Nϕ

    """Index range of each ring along the colatitude (θ) coordinate"""
    iθ

    """Preallocated storage for FTs of individual rings"""
    G

    """Fourier-transform plan

    This transforms physical-space function values to the Fourier domain on
    each ring of colatitude.
    """
    plan
end

"""
    SSHTRS(s, ℓₘₐₓ; kwargs...)

Construct a `SSHTRS` object directly.  This may also be achieved by calling the
main `SSHT` function with the same keywords, along with `method="RS"`.

This object uses the algorithm described in [this paper by Reinecke and
Seljebotn](https://arxiv.org/abs/1303.4945).

The basic floating-point number type may be adjusted with the keyword
argument `T`, which defaults to `Float64`.

The SSHs are evaluated on a series of "rings" at constant colatitude.  Their
locations are specified by the `θ` keyword argument, which defaults to
`fejer1_rings(2ℓₘₐₓ+1, T)`.  If this is changed, the user should also provide
the corresponding `quadrature_weights` argument — the default being
`fejer1(length(θ), T)`.

On each of these rings, an FFT is performed.  To reach the band limit of
``m = ± ℓₘₐₓ``, the number of points along each ring must therefore be *at
least* ``2ℓₘₐₓ+1``, but may be greater.  For example, if ``2ℓₘₐₓ+1`` does
not factorize neatly into a product of small primes, it may be preferable
to use ``2ℓₘₐₓ+2`` points along each ring.  (In that case, whenever `ℓₘₐₓ`
is 1 less than a power of 2, the number of points will be exactly a power
of 2, which is usually particularly efficient.)  The number of points on
each ring can be modified independently, if given as a vector with the same
length as `θ`, or as a single number which is assumed to be the same for
all rings.

Whenever `T` is either `Float64` or `Float32`, the keyword arguments
`plan_fft_flags` and `plan_fft_timelimit` may also be useful for obtaining
more efficient FFTs.  They default to `FFTW.ESTIMATE` and `Inf`,
respectively.  They are passed to
[`AbstractFFTs.plan_fft`](https://juliamath.github.io/AbstractFFTs.jl/stable/api/#AbstractFFTs.plan_fft).
"""
function SSHTRS(
    s, ℓₘₐₓ; T::Type{TT}=Float64,
    θ=fejer1_rings(2ℓₘₐₓ+1, T),
    quadrature_weights=fejer1(length(θ), T),
    Nϕ=fill(2ℓₘₐₓ+1, length(θ)),
    plan_fft_flags=FFTW.ESTIMATE, plan_fft_timelimit=Inf
) where TT
    @assert size(θ) == size(quadrature_weights) """
        size(θ) should equal size(quadrature_weights)
        size(θ) = $(size(θ))
        size(quadrature_weights) = $(size(quadrature_weights))
    """
    @assert size(θ) == size(Nϕ) """
        size(θ) should equal size(Nϕ)
        size(θ) = $(size(θ))
        size(Nϕ) = $(size(Nϕ))
    """
    # COV_EXCL_START
    if (message = check_threads()) != ""
        @warn """$message
        Computations with SSHTRS can benefit greatly from using all available threads if many
        functions are to be transformed.
        """
    end
    # COV_EXCL_STOP
    iθ = let iθ = cumsum(Nϕ)
        [a:b for (a,b) in eachrow(hcat([1; iθ[begin:end-1].+1], iθ))]
    end
    Gs = [Vector{Complex{TT}}(undef, N) for N ∈ Nϕ]
    plans = if TT ∈ [Float64, Float32]  # Only supported types in FFTW
        [plan_fft!(G, flags=plan_fft_flags, timelimit=plan_fft_timelimit) for G ∈ Gs]
    else
        [plan_fft!(G) for G ∈ Gs]
    end
    SSHTRS{TT}(s, ℓₘₐₓ, θ, quadrature_weights, Nϕ, iθ, Gs, plans)
end

function pixels(𝒯::SSHTRS{T}) where {T}
    let π=T(π)
        [
            @SVector [θ, iϕ * 2π / Nϕ]
            for (θ,Nϕ) ∈ zip(𝒯.θ, 𝒯.Nϕ)
            for iϕ ∈ 0:Nϕ-1
        ]
    end
end

function rotors(𝒯::SSHTRS)
    from_spherical_coordinates.(pixels(𝒯))
end

function Base.:*(𝒯::SSHTRS, f̃)
    s1 = Ysize(abs(𝒯.s), 𝒯.ℓₘₐₓ)
    @assert size(f̃, 1) ≥ s1 """
        Size of input `f̃` along first dimension ($(size(f̃,1))) is insufficient for
        `ℓₘₐₓ` property of input transform `𝒯`; it must be at least $s1.
    """
    s2 = maximum(iθ.stop for iθ ∈ 𝒯.iθ)
    f = similar(f̃, (s2, size(f̃)[2:end]...))
    mul!(f, 𝒯, f̃)
end

function LinearAlgebra.mul!(f, 𝒯::SSHTRS{T}, f̃) where {T}
    s1 = Ysize(abs(𝒯.s), 𝒯.ℓₘₐₓ)
    @assert size(f̃, 1) ≥ s1 """
        Size of input `f̃` along first dimension ($(size(f̃,1))) is insufficient for
        `ℓₘₐₓ` property of input transform `𝒯`; it must be at least $s1.
    """
    s2 = size(f)
    s̃2 = (maximum(iθ.stop for iθ ∈ 𝒯.iθ), size(f̃)[2:end]...)
    @assert s̃2==s2 """
        Size of output `f` is not matched to size of input `f̃`:
        size(f) = $(s2)
        (maximum(iθ.stop for iθ ∈ 𝒯.iθ), size(f̃)[2:end]...) = $(s̃2)
    """

    s = 𝒯.s
    ℓₘₐₓ = 𝒯.ℓₘₐₓ
    mₘₐₓ = ℓₘₐₓ
    f .= false  # Zero out all elements to prepare for accumulation below

    f̃′ = reshape(f̃, size(f̃, 1), :)
    f′ = reshape(f, size(f, 1), :)
    @inbounds let π = T(π)
        for (f̃′ⱼ, f′ⱼ) ∈ zip(eachcol(f̃′), eachcol(f′))
            @threads for m ∈ -mₘₐₓ:mₘₐₓ  # Note: Contrary to R&S, we include negative m
                ℓ₀ = max(abs(s), abs(m))
                for (θ, Fy) ∈ zip(𝒯.θ, 𝒯.G)
                    Fmy = zero(T)
                    λ = λiterator(θ, s, m)
                    for (ℓ, ₛλₗₘ) ∈ zip(ℓ₀:ℓₘₐₓ, λ)
                        Fmy += f̃′ⱼ[Yindex(ℓ, m, abs(s))] * ₛλₗₘ
                    end  # ℓ
                    Fy[1+mod(m, length(Fy))] = Fmy
                end  # (θ, Nϕ, G)
            end  # m
            @threads for y ∈ eachindex(𝒯.G)
                Nϕy = 𝒯.Nϕ[y]
                iθy = 𝒯.iθ[y]
                Fy = 𝒯.G[y]
                plany = 𝒯.plan[y]
                plany \ Fy
                @. f′ⱼ[iθy] = Fy * Nϕy
            end
        end  # (f̃′ⱼ, f′ⱼ)
    end  # π

    f
end

function Base.:\(𝒯::SSHTRS, f)
    f̃ = similar(f, (Ysize(abs(𝒯.s), 𝒯.ℓₘₐₓ), size(f)[2:end]...))
    ldiv!(f̃, 𝒯, f)
end

# Compute `f̃ = 𝒯 \ f`, storing the result in `f̃`
function LinearAlgebra.ldiv!(f̃, 𝒯::SSHTRS{T}, f) where {T}
    s1 = maximum(iθ.stop for iθ ∈ 𝒯.iθ)
    @assert size(f, 1) ≥ s1 """
        Size of input `f` along first dimension ($(size(f,1))) is insufficient for
        `Nϕ` property of input transform `𝒯`; it must be at least $s1.
    """
    s̃2 = size(f̃)
    s2 = (Ysize(abs(𝒯.s), 𝒯.ℓₘₐₓ), size(f)[2:end]...)
    @assert s̃2==s2 """
        Size of output `f̃` is not matched to size of input `f`:
        size(f̃) = $(s̃2)
        (Ysize(abs(𝒯.s), ℓₘₐₓ), size(f)[2:end]...) = $(s2)
    """

    s = 𝒯.s
    ℓₘₐₓ = 𝒯.ℓₘₐₓ
    mₘₐₓ = ℓₘₐₓ
    f̃ .= false  # Zero out all elements to prepare for accumulation below

    f̃′ = reshape(f̃, size(f̃, 1), :)
    f′ = reshape(f, size(f, 1), :)
    @inbounds let π = T(π)
        for (f̃′ⱼ, f′ⱼ) ∈ zip(eachcol(f̃′), eachcol(f′))
            @threads for y ∈ eachindex(𝒯.G)
                wy = 𝒯.quadrature_weight[y]
                Nϕy = 𝒯.Nϕ[y]
                iθy = 𝒯.iθ[y]
                Gy = 𝒯.G[y]
                plany = 𝒯.plan[y]
                @. Gy = f′ⱼ[iθy] * wy * 2π / Nϕy
                plany * Gy
            end
            @threads for m ∈ -mₘₐₓ:mₘₐₓ  # Note: Contrary to R&S, we include negative m
                ℓ₀ = max(abs(s), abs(m))
                for (θ, Gy) ∈ zip(𝒯.θ, 𝒯.G)
                    Gmy = Gy[1+mod(m, length(Gy))]
                    λ = λiterator(θ, s, m)
                    for (ℓ, ₛλₗₘ) ∈ zip(ℓ₀:ℓₘₐₓ, λ)
                        # Be careful of the following when adding threads!!!
                        # We need this element of f̃′ⱼ to be used in only one thread,
                        # and we need G to be used in only one thread at a time.
                        f̃′ⱼ[Yindex(ℓ, m, abs(s))] += Gmy * ₛλₗₘ
                    end  # ℓ
                end  # (θ, Nϕ, G)
            end  # m
        end  # (f̃′ⱼ, f′ⱼ)
    end  # π

    f̃
end
