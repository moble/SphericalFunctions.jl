
"""Helper function for [`λ_recursion_initialize`](@ref)"""
binom(T, n, k) = T(binomial(big(n), big(k)))

@doc raw"""
    λ_recursion_initialize(cosθ, sin½θ, cos½θ, s, ℓ, m)

This provides initial values for the recursion to find
``{}_{s}\lambda_{\ell,m}`` along indices of increasing ``\ell``, due to
[Kostelec & Rockmore](https://dx.doi.org/10.1007/s00041-008-9013-5).
Specifically, this function computes values with ``\ell=m``.

{}_{s}\lambda_{\ell,m}(\theta) := {}_{s}Y_{\ell,m}(\theta, 0) = (-1)^m\, \sqrt{\frac{2\ell+1}{4\pi}} d^\ell_{-m,s}(\theta)

"""
function λ_recursion_initialize(sin½θ::T, cos½θ::T, s, ℓ, m) where T
    if abs(m) != ℓ
        @error "Value of m=$m can only be ±ℓ=±$ℓ for this initial-value function.  (θ=$θ; s=$s)."
    end
    if abs(s) > abs(m)
        λ_recursion_initialize(-sin½θ, cos½θ, m, ℓ, s)
    else
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
struct SSHTRS <: SSHT
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

    """Fourier-transform plan

    This transforms physical-space function values to the Fourier domain on
    each ring of colatitude.
    """
    plan

    """Preallocated storage for FTs of individual rings"""
    G
end


function SSHTRS(
    s, ℓₘₐₓ; T=Float64,
    θ=clenshaw_curtis_rings(s, ℓₘₐₓ, T),
    quadrature_weights=clenshaw_curtis(length(θ), T),
    Nϕ=fill(2ℓₘₐₓ+2, length(θ)),
    plan_fft_flags=FFTW.ESTIMATE, plan_fft_timelimit=Inf
)
    @assert size(θ) == size(quadrature_weights) "size(θ) should equal size(quadrature_weights)"
    if (message = check_threads()) != ""
        @warn """$message
        Computations with SSHTRS can benefit greatly from using all available threads if many
        functions are to be transformed.
        """
    end
    m′ₘₐₓ = abs(s)
    Gs = [Vector{Complex{T}}(undef, N) for N ∈ Nϕ]
    plans = [plan_fft(G, flags=plan_fft_flags, timelimit=plan_fft_timelimit) for G ∈ Gs]
    SSHTRS(s, ℓₘₐₓ, θ, quadrature_weights, Nϕ, plans, Gs)
end


function LinearAlgebra.ldiv!(𝐟̃, 𝒯::SSHTRS, 𝐟)  # Compute `𝐟̃ = 𝒯 \ 𝐟`, storing the result in `𝐟̃`
    s1 = size(𝐟̃)
    s2 = (Ysize(ℓₘₐₓ), size(𝐟)[3:end]...)
    @assert s1==s2 "size(𝐟̃)=$s1  !=  (Ysize(ℓₘₐₓ), size(𝐟)[3:end]...)=$s2"

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

    s = 𝒯.s
    lmax = 𝒯.ℓₘₐₓ
    mmax = lmax
    𝐟̃[:] .= false  # Zero out all elements to prepare for accumulation below
    𝐟̃′ = reshape(𝐟̃, size(𝐟̃, 1), :)


    # Based loosely on Fig. 2 of Reinecke & Seljebotn
    @threads for (G,y) ∈ zip(𝒯.G, axes(𝐟, ))
        iₜ = threadid()
        for j ∈ jobs
            G!(G[j,iₜ,:], 𝐟[j,y])
        end  # j
    end  # y

    @threads for m ∈ -mmax:mmax  # Note: Contrary to R&S, we include negative m
        ℓ₀ = max(abs(s), abs(m))
        for y ∈ b.rings
            let θ = θ[y]
                cosθ = cos(θ)
                sin½θ, cos½θ = sincos(θ/2)
                ₛλₗ₋₁ₘ = zero(T)
                ₛλₗₘ = λ_recursion_initialize(sin½θ, cos½θ, s, ℓ₀, m)
                cₗ₋₁ = zero(T)
                for ℓ ∈ ℓ₀:lmax
                    @turbo for j ∈ jobs
                        𝐟̃[j,l,m] += G[j,m,y] * ₛλₗₘ
                    end  # j
                    if ℓ < lmax
                        cₗ₊₁, cₗ = λ_recursion_coefficients(cosθ, s, ℓ, m)
                        ₛλₗ₊₁ₘ = if ℓ == 0
                            √(3/4T(π)) * cosθ
                        else
                            (cₗ * ₛλₗₘ + cₗ₋₁ * ₛλₗ₋₁ₘ) / cₗ₊₁
                        end
                        ₛλₗ₋₁ₘ = ₛλₗₘ
                        ₛλₗₘ = ₛλₗ₊₁ₘ
                        cₗ₋₁ = -cₗ₊₁ * √((2ℓ+1)/T(2ℓ+3))
                    end
                end  # l
            end  # θ
        end  # y
    end  # m

    𝐟̃
end
