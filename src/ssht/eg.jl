"""Storage for spin-spherical-harmonic transform on an equiangular grid

"""
struct SSHTEG <: SSHT
    """Spin weight"""
    s

    """Highest ℓ value present in the data"""
    ℓₘₐₓ

    """Locations of rings along the colatitude (θ) coordinate"""
    θ

    """Number of points along the azimuthal (ϕ) coordinate in each ring"""
    Nϕ

    """Quadrature weights for the given θ values"""
    quadrature_weights

    """Fourier-transform plan

    This transforms physical-space function values to the Fourier domain on
    each ring of colatitude.
    """
    plans

    """Preallocated storage for FTs of individual rings"""
    Gs

    """Preallocated storage for computations of Wigner-d values at a given angles"""
    Hwedge

    """Preallocated storage for constants in computations of Wigner-d values"""
    H_rec_coeffs
end


function SSHTEG(
    s, ℓₘₐₓ; T=Float64,
    θ=clenshaw_curtis_rings(s, ℓₘₐₓ, T), Nϕ=2ℓₘₐₓ+2, quadrature_weights=clenshaw_curtis(length(θ), T),
    plan_fft_flags=FFTW.ESTIMATE, plan_fft_timelimit=Inf
)
    @assert size(θ) == size(quadrature_weights) "size(θ) should equal size(quadrature_weights)"
    if (message = check_threads()) != ""
        @warn """$message
        Computations with SSHTEG can benefit greatly from using all available threads if many
        functions are to be transformed.
        """
    end
    m′ₘₐₓ = abs(s)
    Gs = [Vector{Complex{T}}(undef, Nϕ) for _ ∈ 1:nthreads()]
    plans = [plan_fft(G, flags=plan_fft_flags, timelimit=plan_fft_timelimit) for G ∈ Gs]
    Hwedge = Array{T}(undef, WignerHsize(ℓₘₐₓ, m′ₘₐₓ))
    H_rec_coeffs = H_recursion_coefficients(ℓₘₐₓ, T)
    SSHTEG(s, ℓₘₐₓ, θ, Nϕ, quadrature_weights, plans, Gs, Hwedge, H_rec_coeffs)
end


function LinearAlgebra.ldiv!(f̃, ssht::SSHTEG, f)  # Compute `f̃ = ssht \ f`, storing the result in `f̃`
    s1 = size(f̃)
    s2 = (Ysize(ℓₘₐₓ), size(f)[3:end]...)
    @assert s1==s2 "size(f̃)=$s1  !=  (Ysize(ℓₘₐₓ), size(f)[3:end]...)=$s2"

    # # Eq. (10) of Reinecke & Seljebotn https://dx.doi.org/10.1051/0004-6361/201321494
    # ₛλₗₘ(ϑ) = (-1)ᵐ √((2ℓ+1)/(4π)) dˡ₋ₘₛ(ϑ)
    #
    # # Eq. (4.11) of Kostelec & Rockmore https://dx.doi.org/10.1007/s00041-008-9013-5
    # # Note that terms with out-of-range indices should be treated as 0.
    # ₛλₗ₊₁ₘ = √((2ℓ+3)/(2ℓ+1)) (ℓ+1) (2ℓ+1) / √(((ℓ+1)²-m²) ((ℓ+1)²-s²)) (cosϑ + ms/(ℓ(ℓ+1))) ₛλₗₘ
    #          -  √((2ℓ+3)/(2ℓ-1)) (ℓ+1) (2ℓ+1) √((ℓ-m²) (ℓ-s²)) / √(((ℓ+1)²-m²) ((ℓ+1)²-s²)) ((ℓ+1)/ℓ) ₛλₗ₋₁ₘ
    #
    # # Eqs. (4.6) and (4.7) of Kostelec & Rockmore
    #
    # # https://en.wikipedia.org/wiki/Wigner_D-matrix#Symmetries_and_special_cases
    # dˡ₋ₘₛ(ϑ) = (-1)ˡ⁺ᵐ dˡ₋ₘ₋ₛ(π-ϑ)
    #  ₛλₗₘ(ϑ) = (-1)ˡ⁺ᵐ  ₋ₛλₗₘ(π-ϑ)

    s = ssht.s
    ℓₘₐₓ = ssht.ℓₘₐₓ
    abss = abs(s)
    m′ₘₐₓ = abss
    ϵ₋ₛ = ϵ(-s)
    f̃[:] .= false  # Zero out all elements to prepare for accumulation below
    f̃′ = reshape(f̃, size(f̃, 1), :)

    @inbounds for (θ, weight) ∈ zip(ssht.θ, ssht.quadrature_weights)
        H!(ssht.Hwedge, cis(θ), ℓₘₐₓ, m′ₘₐₓ, ssht.H_rec_coeffs, WignerHindex)
        f′ = reshape(selectdim(f, 2, iθ), Nϕ, :)
        @threads for iextra ∈ axes(f̃′, 2)  # for extra ∈ collect(extra_dims)
            # NOTE: We can't thread at a higher level because each thread could access the
            # same element of `f̃` simultaneously below; by threading at this level, we
            # are assured that the `iextra` index used below is unique to each thread.
            iₜ = Threads.threadid()
            G = Gs[iₜ]
            computeG!(G, f′[:, iextra], weight, ssht.plans[iₜ])
            for ℓ ∈ abss:ℓₘₐₓ
                λ_factor = ϵ₋ₛ * √((2ℓ+1)*T(π)) / Nϕ

                i₀ = WignerHindex(ℓ, s, 0, m′ₘₐₓ)

                let m=0
                    f̃′[Yindex(ℓ, m), iextra] += G[m+1] * λ_factor * Hwedge[i₀]
                end

                i₊ = i₀
                i₋ = i₀
                if !signbit(s)
                    for m ∈ 1:min(ℓ, abss)
                        i₊ -= ℓ-m+2
                        i₋ += ℓ-m+1
                        f̃′[Yindex(ℓ, m), iextra] += G[m+1] * ϵ(m) * λ_factor * Hwedge[i₊]
                        f̃′[Yindex(ℓ, -m), iextra] += G[Nϕ-m+1] * λ_factor * Hwedge[i₋]
                    end
                else
                    for m ∈ 1:min(ℓ, abss)
                        i₊ += ℓ-m+1
                        i₋ -= ℓ-m+2
                        f̃′[Yindex(ℓ, m), iextra] += G[m+1] * ϵ(m) * λ_factor * Hwedge[i₊]
                        f̃′[Yindex(ℓ, -m), iextra] += G[Nϕ-m+1] * λ_factor * Hwedge[i₋]
                    end
                end
                for m ∈ abss+1:ℓ
                    i₊ += 1
                    i₋ += 1
                    f̃′[Yindex(ℓ, m), iextra] += G[m+1] * ϵ(m) * λ_factor * Hwedge[i₊]
                    f̃′[Yindex(ℓ, -m), iextra] += G[Nϕ-m+1] * λ_factor * Hwedge[i₋]
                end
            end
        end
    end
    f̃
end
