using Quaternionic: from_spherical_coordinates

# TODO: Move λ recursion to separate file, as way to evaluate ₛYₗₘs
# TODO: Deal with binomial coefficients in `λ_recursion_initialize` better
# TODO: Add enough `G` storage to use on multiple threads / SIMD registers (`plan`s are thread-safe)
# TODO: Reorganize to match actual R&S algorithm recommendations

"""Helper function for [`λ_recursion_initialize`](@ref)"""
binom(T, n, k) = T(binomial(big(n), big(k)))

# # Eq. (10) of Reinecke & Seljebotn https://dx.doi.org/10.1051/0004-6361/201321494
# ₛλₗₘ(ϑ) = (-1)ᵐ √((2ℓ+1)/(4π)) dˡ₋ₘₛ(ϑ)
#
# # Eq. (4.11) of Kostelec & Rockmore https://dx.doi.org/10.1007/s00041-008-9013-5
# # Note that terms with out-of-range indices should be treated as 0.
# ₛλₗ₊₁ₘ = √((2ℓ+3)/(2ℓ+1)) (ℓ+1) (2ℓ+1) / √(((ℓ+1)²-m²) ((ℓ+1)²-s²)) (cosϑ + ms/(ℓ(ℓ+1))) ₛλₗₘ
#          -  √((2ℓ+3)/(2ℓ-1)) (ℓ+1) (2ℓ+1) √((ℓ-m²) (ℓ-s²)) / √(((ℓ+1)²-m²) ((ℓ+1)²-s²)) ((ℓ+1)/ℓ) ₛλₗ₋₁ₘ
#
# # Eqs. (4.7) and (4.6) of Kostelec & Rockmore
# for 0 ≤ s ≤ ℓ
# ₛλₗₗ(ϑ) = (-1)ᵐ √((2ℓ+1)/(4π)) √(((2ℓ)!)/((ℓ+s)!(ℓ-s)!)) cosˡ⁻ˢ ϑ/2 sinˡ⁺ˢ ϑ/2
# ₛλₗ₋ₗ(ϑ) = (-1)ᵐ⁺ˡ⁺ˢ √((2ℓ+1)/(4π)) √(((2ℓ)!)/((ℓ+s)!(ℓ-s)!)) cosˡ⁺ˢ ϑ/2 sinˡ⁻ˢ ϑ/2
#
# # https://en.wikipedia.org/wiki/Wigner_D-matrix#Symmetries_and_special_cases
# dˡ₋ₘₛ(ϑ) = (-1)ˡ⁺ᵐ dˡ₋ₘ₋ₛ(π-ϑ)
#  ₛλₗₘ(ϑ) = (-1)ˡ⁺ᵐ  ₋ₛλₗₘ(π-ϑ)
#
# for -ℓ ≤ s ≤ 0
# ₛλₗₗ(ϑ) = (-1)ˡ √((2ℓ+1)/(4π)) √(((2ℓ)!)/((ℓ+s)!(ℓ-s)!)) cosˡ⁺ˢ (π-ϑ)/2 sinˡ⁻ˢ (π-ϑ)/2
# ₛλₗ₋ₗ(ϑ) = (-1)ˢ √((2ℓ+1)/(4π)) √(((2ℓ)!)/((ℓ+s)!(ℓ-s)!)) cosˡ⁻ˢ (π-ϑ)/2 sinˡ⁺ˢ (π-ϑ)/2

@doc raw"""
    λ_recursion_initialize(cosθ, sin½θ, cos½θ, s, ℓ, m)

This provides initial values for the recursion to find
``{}_{s}\lambda_{\ell,m}`` along indices of increasing ``\ell``, due to
[Kostelec & Rockmore](https://dx.doi.org/10.1007/s00041-008-9013-5).
Specifically, this function computes values with ``\ell=m``.

{}_{s}\lambda_{\ell,m}(\theta) := {}_{s}Y_{\ell,m}(\theta, 0) = (-1)^m\, \sqrt{\frac{2\ell+1}{4\pi}} d^\ell_{-m,s}(\theta)

"""
function λ_recursion_initialize(sin½θ::T, cos½θ::T, s, ℓ, m) where T
    if abs(s) > abs(m)
        λ_recursion_initialize(-sin½θ, cos½θ, m, ℓ, s)
    else
        if abs(m) != ℓ
            error("""Value of m=$m can only be ±ℓ=±$ℓ for this initial-value function.
                s=$s
                sin½θ=$sin½θ
                cos½θ=$cos½θ
            """)
        end
        let π = T(π)
            c = √((2ℓ+1) * binom(T, 2ℓ, ℓ-abs(s)) / (4π))
            if s < 0
                if m == ℓ
                    (-1)^ℓ * c * sin½θ^(ℓ+s) * cos½θ^(ℓ-s)
                else # m == -ℓ
                    (-1)^s * c * sin½θ^(ℓ-s) * cos½θ^(ℓ+s)
                end
            else
                if m == ℓ
                    (-1)^m * c * sin½θ^(ℓ+s) * cos½θ^(ℓ-s)
                else # m == -ℓ
                    (-1)^(ℓ+s+m) * c * sin½θ^(ℓ-s) * cos½θ^(ℓ+s)
                end
            end
        end
    end
end

function λ_recursion_coefficients(cosθ::T, s, ℓ, m) where T
    cₗ₊₁ = √(((ℓ+1)^2-m^2) * ((ℓ+1)^2-s^2) / T(2ℓ+3)) / T((ℓ+1)*(2ℓ+1))
    cₗ = (cosθ + m*s/T(ℓ*(ℓ+1))) / √T(2ℓ+1)
    cₗ₊₁, cₗ
end




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
    s, ℓₘₐₓ; T=Float64,
    θ=fejer1_rings(2ℓₘₐₓ+1, T),
    quadrature_weights=fejer1(length(θ), T),
    Nϕ=fill(2ℓₘₐₓ+1, length(θ)),
    plan_fft_flags=FFTW.ESTIMATE, plan_fft_timelimit=Inf
)
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
    if (message = check_threads()) != ""
        @warn """$message
        Computations with SSHTRS can benefit greatly from using all available threads if many
        functions are to be transformed.
        """
    end
    iθ = let iθ = cumsum(Nϕ)
        [a:b for (a,b) in eachrow(hcat([1; iθ[begin:end-1].+1], iθ))]
    end
    Gs = [Vector{Complex{T}}(undef, N) for N ∈ Nϕ]
    plans = if T ∈ [Float64, Float32]  # Only supported types in FFTW
        [plan_fft(G, flags=plan_fft_flags, timelimit=plan_fft_timelimit) for G ∈ Gs]
    else
        [plan_fft(G) for G ∈ Gs]
    end
    SSHTRS{T}(s, ℓₘₐₓ, θ, quadrature_weights, Nϕ, iθ, Gs, plans)
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
            for (wy, Nϕy, iθy, Gy, plany) ∈ zip(𝒯.quadrature_weight, 𝒯.Nϕ, 𝒯.iθ, 𝒯.G, 𝒯.plan)
                mul!(Gy, plany, f′ⱼ[iθy])
                @. Gy *= wy * 2π / Nϕy
            end
            for m ∈ -mₘₐₓ:mₘₐₓ  # Note: Contrary to R&S, we include negative m
                ℓ₀ = max(abs(s), abs(m))
                for (θ, Gy) ∈ zip(𝒯.θ, 𝒯.G)
                    Gmy = Gy[1+mod(m, length(Gy))]
                    cosθ = cos(θ)
                    sin½θ, cos½θ = sincos(θ/2)
                    ₛλₗ₋₁ₘ = zero(T)
                    ₛλₗₘ = λ_recursion_initialize(sin½θ, cos½θ, s, ℓ₀, m)
                    cₗ₋₁ = zero(T)
                    for ℓ ∈ ℓ₀:ℓₘₐₓ
                        lm = Yindex(ℓ, m, abs(s))
                        # Be careful of the following when adding threads!!!
                        # We need this element of f̃′ⱼ to be used in only one thread,
                        # and we need G to be used in only one thread at a time.
                        f̃′ⱼ[lm] += Gmy * ₛλₗₘ
                        if ℓ < ℓₘₐₓ  # Take another step in the λ recursion
                            cₗ₊₁, cₗ = λ_recursion_coefficients(cosθ, s, ℓ, m)
                            ₛλₗ₊₁ₘ = if ℓ == 0
                                # The only case in which this will ever be used is when
                                # s == m == ℓ == 0.  So we want ₀Y₁₀, which is known:
                                √(3/4π) * cosθ
                            else
                                (cₗ * ₛλₗₘ + cₗ₋₁ * ₛλₗ₋₁ₘ) / cₗ₊₁
                            end
                            ₛλₗ₋₁ₘ = ₛλₗₘ
                            ₛλₗₘ = ₛλₗ₊₁ₘ
                            cₗ₋₁ = -cₗ₊₁ * √((2ℓ+1)/T(2ℓ+3))
                        end
                    end  # ℓ
                end  # (θ, Nϕ, G)
            end  # m
        end  # (f̃′ⱼ, f′ⱼ)
    end  # π

    f̃
end
