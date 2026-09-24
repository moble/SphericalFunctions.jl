"""
    SSHTRS(s, ℓₘₐₓ; T=Float64, [θ, quadrature_weights], Nϕ=2ℓₘₐₓ+1, plan_fft_flags=FFTW.ESTIMATE, plan_fft_timelimit=Inf)

Construct an ``s``-SHT object that uses the ring-based algorithm described by [Reinecke and
Seljebotn](@cite Reinecke_2013).  This may also be achieved by calling the main [`SSHT`](@ref)
function with the same keywords, along with `method="RS"` (the default).

The spin-weighted spherical harmonics are evaluated on a series of "rings" at constant
colatitude, whose locations are given by the `θ` keyword argument, and the analysis
integrates over ``θ`` with the `quadrature_weights` of the rule that placed those rings.  When
both are omitted they are the nodes and weights of Fejér's first rule,
`fejer1_rings(2ℓₘₐₓ+1, T)` and `fejer1(2ℓₘₐₓ+1, T)`.  Only the caller knows which rule
placed a given set of rings, so the two must be given together — for example
`θ=clenshaw_curtis_rings(N, T)` with `quadrature_weights=clenshaw_curtis(N, T)` — and either
one without the other is refused.
The analysis is exact for band-limited functions when the quadrature rule integrates
polynomials of degree ``2ℓₘₐₓ`` in ``\\cos θ`` exactly, as the Fejér and Clenshaw–Curtis
rules with at least ``2ℓₘₐₓ+1`` nodes do.  The constructor checks this, and warns when the
rule falls short; synthesis does not use the weights, and is exact on any rings.

On each ring, an FFT is performed.  To reach the band limit of ``m = ±ℓₘₐₓ``, the number of
points along each ring must be *at least* ``2ℓₘₐₓ+1``, but may be greater.  For example, if
``2ℓₘₐₓ+1`` does not factorize neatly into a product of small primes, it may be preferable to
use ``2ℓₘₐₓ+2`` points along each ring.  The number of points on each ring can be modified
independently, if given as a vector with the same length as `θ`, or as a single number which
is used for all rings.

The sample points are ordered ring by ring, with the azimuth ``ϕ_k = 2πk/N_ϕ`` for
``k = 0, …, N_ϕ-1`` varying fastest; see [`pixels`](@ref) and [`rotors`](@ref).

Whenever `T` is either `Float64` or `Float32`, the keyword arguments `plan_fft_flags` and
`plan_fft_timelimit` may also be useful for obtaining more efficient FFTs.  They default to
`FFTW.ESTIMATE` and `Inf`, respectively, and are passed to
[`AbstractFFTs.plan_fft!`](https://juliamath.github.io/AbstractFFTs.jl/stable/api/#AbstractFFTs.plan_fft).

The harmonics on the rings are computed with one batched [`sλlmCalculator`](@ref), one ``ℓ``
at a time, so the cost is ``O(N_θ ℓₘₐₓ^2)`` and the memory ``O(N_θ ℓₘₐₓ)``.  The object holds
that workspace, so it must not be used from several threads at once.

Half-integer `s` and `ℓₘₐₓ` are accepted as `Rational`s with denominator 2; see
[`SSHT`](@ref) for what the function values then mean.  The algorithm is the same one: the
harmonics on a ring are ``i^{2s}`` times real functions of ``θ``, so the ``θ`` stage works
with those real functions and restores the phase once per ring, and ``e^{imϕ}`` with half-odd
``m`` is ``e^{iϕ/2}`` times a Fourier mode of integer frequency ``m - 1/2``, so the FFT runs
at integer frequencies and each sample of a ring is multiplied by ``e^{±iϕ/2}``.  The
defaults — ``2ℓₘₐₓ+1`` rings and ``2ℓₘₐₓ+1`` points per ring — are unchanged, and are even
numbers.
"""
struct SSHTRS{T<:Real, ST, P, BP, B, IT<:IntegerHalf} <: SSHT{T}
    s::IT
    ℓₘₐₓ::IT
    θ::Vector{T}
    quadrature_weights::Vector{T}
    Nϕ::Vector{Int}
    iθ::Vector{UnitRange{Int}}  # index range of each ring in the pixel vector
    λ::sλlmCalculator{IT, T, ST, IT, B}  # ₛλₗₘ(θ) for all rings at once (angle mode)
    F::Matrix{Complex{T}}  # Fourier coefficients [ring, m+ℓₘₐₓ+1] for m ∈ -ℓₘₐₓ:ℓₘₐₓ
    G::Vector{Vector{Complex{T}}}  # per-ring FFT buffers
    plans::Vector{P}  # forward in-place FFT plans (analysis)
    bplans::Vector{BP}  # backward (unnormalized inverse) in-place FFT plans (synthesis)
end

# The public constructor is the boundary: it normalizes the two indices and re-dispatches
# to the worker, whose keyword defaults are then computed from indices of one kind.  It also
# checks the quadrature rule.  `map2salm_plan` and `salm2map` call the worker directly: they
# build the Clenshaw–Curtis rule themselves, and the first states the number of rings that
# rule needs more plainly than the general check could, while the second only synthesizes.
function SSHTRS(s::IndexArgument, ℓₘₐₓ::IndexArgument; T::Type{TT}=Float64, kwargs...) where {TT}
    𝒯 = SSHTRS(transform_indices(s, ℓₘₐₓ)..., TT; kwargs...)
    warn_if_inexact(𝒯)
    𝒯
end
function SSHTRS(
    s::IT, ℓₘₐₓ::IT, ::Type{TT};
    θ=nothing,
    quadrature_weights=nothing,
    Nϕ=2ℓₘₐₓ+1,
    plan_fft_flags=FFTW.ESTIMATE, plan_fft_timelimit=Inf
) where {IT<:IntegerHalf, TT}
    if abs(s) > ℓₘₐₓ
        error("|s|=$(abs(s)) exceeds ℓₘₐₓ=$ℓₘₐₓ; there are no such modes.")
    end
    # The weights belong to the rule that placed the rings, and only the caller knows which
    # rule that was; a default for one of the two would silently pair it with another rule.
    if θ === nothing && quadrature_weights === nothing
        θ = fejer1_rings(2ℓₘₐₓ+1, TT)
        quadrature_weights = fejer1(2ℓₘₐₓ+1, TT)
    elseif quadrature_weights === nothing
        throw(ArgumentError(
            "The rings `θ` were given without their `quadrature_weights`.  The weights belong "
            * "to the rule that placed the rings — for example `clenshaw_curtis(length(θ), T)` "
            * "for `clenshaw_curtis_rings` — and must be given with them."
        ))
    elseif θ === nothing
        throw(ArgumentError(
            "The `quadrature_weights` were given without the rings `θ` of their rule; the two "
            * "must be given together."
        ))
    end
    check_sample_reals(TT, θ, "θ")
    check_sample_reals(TT, quadrature_weights, "quadrature_weights")
    θ = Vector{TT}(θ)
    quadrature_weights = Vector{TT}(quadrature_weights)
    if length(θ) != length(quadrature_weights)
        error(
            "θ and quadrature_weights must have the same length; got $(length(θ)) and "
            * "$(length(quadrature_weights))."
        )
    end
    Nϕ = Nϕ isa Integer ? fill(Int(Nϕ), length(θ)) : Vector{Int}(Nϕ)
    if length(Nϕ) != length(θ)
        error("Nϕ must be a single number or have the same length as θ ($(length(θ))).")
    end
    if any(<(1), Nϕ)
        error("Every ring needs at least one point; got Nϕ=$Nϕ.")
    end
    if any(<(2ℓₘₐₓ+1), Nϕ)
        @warn "Some rings have fewer than 2ℓₘₐₓ+1=$(2ℓₘₐₓ+1) points, so modes with large |m| will alias."
    end
    Nθ = length(θ)
    iθ = let stops = cumsum(Nϕ)
        [(stop - n + 1):stop for (n, stop) ∈ zip(Nϕ, stops)]
    end
    λ = sλlmCalculator(θ, ℓₘₐₓ, s)  # θ is already a Vector{TT}
    F = Matrix{Complex{TT}}(undef, Nθ, 2ℓₘₐₓ + 1)
    G = [Vector{Complex{TT}}(undef, n) for n ∈ Nϕ]
    plans, bplans = if TT ∈ (Float64, Float32)  # Only supported types in FFTW
        (
            [plan_fft!(g; flags=plan_fft_flags, timelimit=plan_fft_timelimit) for g ∈ G],
            [plan_bfft!(g; flags=plan_fft_flags, timelimit=plan_fft_timelimit) for g ∈ G],
        )
    else
        ([plan_fft!(g) for g ∈ G], [plan_bfft!(g) for g ∈ G])
    end
    SSHTRS{TT, typeof(parent(λ.H.Hˡ)), eltype(plans), eltype(bplans), isbatched(λ), IT}(
        s, ℓₘₐₓ, θ, quadrature_weights, Nϕ, iθ, λ, F, G, plans, bplans
    )
end

# The analysis integrates, ring by ring, products ₛλₗₘ ₛλₗ′ₘ of harmonics with the same m,
# and each such product is a polynomial in cos θ of degree ℓ+ℓ′ ≤ 2ℓₘₐₓ, for half-odd indices
# as for integers.  The analysis is therefore exact when the quadrature rule integrates every
# polynomial of that degree exactly, which is checked here with the Legendre moments: the sum
# Σ_y w_y P_k(cos θ_y) must be 2δ_{k0} for each k ∈ 0:2ℓₘₐₓ.  Unlike the monomials, the P_k
# are bounded by 1 on the whole interval, so the moments are well conditioned, and a fixed
# multiple of the rounding error of the sums serves as the tolerance.  The degree is 2ℓₘₐₓ
# itself rather than 2⌊ℓₘₐₓ⌋, so that a rule without the reflection symmetry of the standard
# ones is checked in every degree the analysis needs; a symmetric rule integrates the odd
# degrees exactly in any case, which is why 2ℓₘₐₓ of its rings suffice for a half-odd ℓₘₐₓ.
function warn_if_inexact(𝒯::SSHTRS{T}) where {T}
    degree = 2𝒯.ℓₘₐₓ
    moments = zeros(T, degree + 1)  # moments[k+1] = Σ_y w_y P_k(cos θ_y)
    for (θ, w) ∈ zip(𝒯.θ, 𝒯.quadrature_weights)
        x = cos(θ)
        P₋, P = zero(T), one(T)
        for k ∈ 0:degree
            moments[k+1] += w * P
            P₋, P = P, ((2k + 1) * x * P - k * P₋) / (k + 1)
        end
    end
    moments[1] -= 2
    residual, i = findmax(abs, moments)
    if !(residual ≤ 100 * length(𝒯.θ) * eps(T))  # also catches NaN
        @warn (
            "The quadrature rule given by `θ` and `quadrature_weights` does not integrate "
            * "polynomials of degree 2ℓₘₐₓ=$degree in cos θ exactly: its Legendre moment of "
            * "degree $(i-1) is off by $(round(Float64(residual), sigdigits=2)).  Analysis "
            * "with this transform will therefore not be exact; synthesis is unaffected."
        )
    end
    nothing
end

function pixels(𝒯::SSHTRS{T}) where {T}
    let π = T(π)
        [
            @SVector [θ, iϕ * 2π / Nϕ]
            for (θ, Nϕ) ∈ zip(𝒯.θ, 𝒯.Nϕ)
            for iϕ ∈ 0:Nϕ-1
        ]
    end
end
rotors(𝒯::SSHTRS) = from_spherical_coordinates.(pixels(𝒯))
npixels(𝒯::SSHTRS) = 𝒯.iθ[end].stop

function Base.:*(𝒯::SSHTRS, f̃)
    check_modes(𝒯, f̃)
    mul!(pixel_output(𝒯, f̃), 𝒯, f̃)
end

# Synthesis: f = 𝒯 * f̃
function LinearAlgebra.mul!(f, 𝒯::SSHTRS{T}, f̃) where {T}
    check_modes(𝒯, f̃)
    check_pixels(𝒯, f)
    check_trailing(f, f̃)
    s, ℓₘₐₓ, Nθ = 𝒯.s, 𝒯.ℓₘₐₓ, length(𝒯.θ)
    λ, F, G = 𝒯.λ, 𝒯.F, 𝒯.G
    f̃′ = reshape(array_view(f̃), size(f̃, 1), :)
    f′ = reshape(f, size(f, 1), :)
    @inbounds for (f̃ⱼ, fⱼ) ∈ zip(eachcol(f̃′), eachcol(f′))
        # Fourier coefficients on each ring: F[y, m] = Σ_ℓ f̃ₗₘ ₛλₗₘ(θ_y)
        fill!(F, zero(Complex{T}))
        iₛ = spin_index(λ, s)
        Λ = λ.Yˡ  # [y, spin, m+ℓ+1]
        for ℓ ∈ abs(s):ℓₘₐₓ
            recurrence!(λ, ℓ)
            i₀ = Yindex(ℓ, -ℓ, abs(s)) - 1
            @threads for m ∈ -ℓ:ℓ
                f̃ₗₘ = f̃ⱼ[i₀ + ℓ + m + 1]
                jm = m + ℓₘₐₓ + 1
                jℓ = m + ℓ + 1
                @simd for y ∈ 1:Nθ
                    F[y, jm] += f̃ₗₘ * Λ[y, iₛ, jℓ]
                end
            end
        end
        # Inverse Fourier transform on each ring (aliasing the m values if Nϕ < 2ℓₘₐₓ+1)
        @threads for y ∈ 1:Nθ
            Gy = G[y]
            Nϕy = 𝒯.Nϕ[y]
            fill!(Gy, zero(Complex{T}))
            for m ∈ -ℓₘₐₓ:ℓₘₐₓ
                Gy[1 + mod(fourier_index(m), Nϕy)] += F[y, m + ℓₘₐₓ + 1]
            end
            𝒯.bplans[y] * Gy  # unnormalized inverse FFT: Σₘ Gₘ e^{+imϕₖ}
            ring_values!(view(fⱼ, 𝒯.iθ[y]), Gy, s)
        end
    end
    f
end

function Base.:\(𝒯::SSHTRS, f)
    check_pixels(𝒯, f)
    ldiv!(mode_output(𝒯, f), 𝒯, f)
end

# Analysis: f̃ = 𝒯 \ f
function LinearAlgebra.ldiv!(f̃, 𝒯::SSHTRS{T}, f) where {T}
    f̃ = analysis_output(𝒯, f̃, f)
    check_modes(𝒯, f̃)
    check_pixels(𝒯, f)
    check_trailing(f, f̃)
    s, ℓₘₐₓ, Nθ = 𝒯.s, 𝒯.ℓₘₐₓ, length(𝒯.θ)
    λ, F, G = 𝒯.λ, 𝒯.F, 𝒯.G
    f̃′ = reshape(array_view(f̃), size(f̃, 1), :)
    f′ = reshape(f, size(f, 1), :)
    @inbounds let π = T(π)
        for (f̃ⱼ, fⱼ) ∈ zip(eachcol(f̃′), eachcol(f′))
            # Fourier transform on each ring, including the quadrature weight and the ϕ
            # measure: F[y, m] = w_y (2π/Nϕ) Σₖ f(θ_y, ϕₖ) e^{-imϕₖ}
            @threads for y ∈ 1:Nθ
                Gy = G[y]
                Nϕy = 𝒯.Nϕ[y]
                factor = 𝒯.quadrature_weights[y] * 2π / Nϕy
                ring_samples!(Gy, fⱼ[𝒯.iθ[y]], factor, s)
                𝒯.plans[y] * Gy
                for m ∈ -ℓₘₐₓ:ℓₘₐₓ
                    F[y, m + ℓₘₐₓ + 1] = Gy[1 + mod(fourier_index(m), Nϕy)]
                end
            end
            # Mode weights: f̃ₗₘ = Σ_y F[y, m] ₛλₗₘ(θ_y)
            iₛ = spin_index(λ, s)
            Λ = λ.Yˡ  # [y, spin, m+ℓ+1]
            for ℓ ∈ abs(s):ℓₘₐₓ
                recurrence!(λ, ℓ)
                i₀ = Yindex(ℓ, -ℓ, abs(s)) - 1
                @threads for m ∈ -ℓ:ℓ
                    jm = m + ℓₘₐₓ + 1
                    jℓ = m + ℓ + 1
                    acc = zero(Complex{T})
                    @simd for y ∈ 1:Nθ
                        acc += F[y, jm] * Λ[y, iₛ, jℓ]
                    end
                    f̃ⱼ[i₀ + ℓ + m + 1] = acc
                end
            end
        end
    end
    f̃
end


### Equiangular-grid convenience functions

@doc raw"""
    map2salm(map, s, ℓₘₐₓ)
    map2salm(map, 𝒯::SSHTRS)

Transform function values `map` sampled on an equiangular grid to spin-weighted
spherical-harmonic mode weights ``{}_sa_{ℓ,m}``.

The `map` array must have size ``N_ϕ`` along its first dimension and ``N_θ`` along its
second, with any number of dimensions following, sampled at the Clenshaw–Curtis nodes
[`clenshaw_curtis_rings`](@ref)`(Nθ)` in ``θ`` (which include both poles) and at
``ϕ_k = 2πk/N_ϕ``; [`pixels`](@ref) returns that grid for a constructed
[`SSHTRS`](@ref).  For the analysis to be exact for a band-limited function, one needs
``N_ϕ ≥ 2ℓₘₐₓ+1`` and ``N_θ ≥ 2ℓₘₐₓ+1`` (or ``N_θ ≥ 2ℓₘₐₓ`` for a half-odd ``ℓₘₐₓ``); with
fewer, the result is not exact, and a warning is issued.

The result is a [`ModeWeights`](@ref) for a one-dimensional map (``N_ϕ × N_θ``), or an
array whose first dimension indexes the modes in the canonical ordering `ℓ ∈ abs(s):ℓₘₐₓ,
m ∈ -ℓ:ℓ` otherwise.  Note that, unlike the pre-3.0 version of this function, the output
starts at ``ℓ = |s|`` rather than ``ℓ = 0``.

For repeated use with different `map`s of the same shape, construct the transform once with
`map2salm_plan(map, s, ℓₘₐₓ)` (an [`SSHTRS`](@ref) on the Clenshaw–Curtis rings) and pass it
as the second argument.  See also [`salm2map`](@ref).

The spin weight and ``ℓₘₐₓ`` may be half-integers, passed as `Rational`s with denominator 2.
The map is then antiperiodic in ``ϕ`` — its values are those of the function at the rotors
`from_spherical_coordinates(θ, ϕ)`, and a full circuit of the azimuth reaches the antipodal
rotor — and the requirements ``N_ϕ ≥ 2ℓₘₐₓ+1``, ``N_θ ≥ 2ℓₘₐₓ+1`` are unchanged in form; see
[`SSHT`](@ref).
"""
function map2salm end

function map2salm(map::MapOrModes, s::IndexArgument, ℓₘₐₓ::IndexArgument)
    map2salm(map, map2salm_plan(map, s, ℓₘₐₓ))
end
function map2salm(map::MapOrModes, 𝒯::SSHTRS)
    Nϕ, Nθ = size(map, 1), size(map, 2)
    if 𝒯.Nϕ != fill(Nϕ, Nθ) || length(𝒯.θ) != Nθ
        error("The transform was planned for a different grid than the $(Nϕ)×$(Nθ) map.")
    end
    check_clenshaw_curtis(𝒯, "map2salm")
    f = reshape(map, Nϕ * Nθ, size(map)[3:end]...)
    𝒯 \ f
end

# `map2salm` and `salm2map` promise the Clenshaw–Curtis grid, but the shape alone does not
# establish it: any `SSHTRS` with the right numbers of rings and points has that shape, and in
# particular the default `SSHT(s, ℓₘₐₓ)`, on the Fejér rings, has exactly the shape of the
# smallest such grid.  The rings and weights themselves are therefore compared.  (The rings
# are compared first, so the weights are not computed for a degenerate single ring.)
function check_clenshaw_curtis(𝒯::SSHTRS{T}, name) where {T}
    Nθ = length(𝒯.θ)
    if !(𝒯.θ ≈ clenshaw_curtis_rings(Nθ, T)) || !(𝒯.quadrature_weights ≈ clenshaw_curtis(Nθ, T))
        error(
            "$name works on the Clenshaw–Curtis grid, but this transform has other rings or "
            * "quadrature weights; construct it with `map2salm_plan(map, s, ℓₘₐₓ)`."
        )
    end
end

"""
    map2salm_plan(map, s, ℓₘₐₓ)

Construct the [`SSHTRS`](@ref) transform used by [`map2salm`](@ref) and [`salm2map`](@ref)
for maps of the shape of `map` (``N_ϕ × N_θ × …``) on the Clenshaw–Curtis grid.
"""
function map2salm_plan(map::AbstractArray{Complex{T}}, s::IndexArgument, ℓₘₐₓ::IndexArgument) where {T<:Real}
    Nϕ, Nθ = size(map, 1), size(map, 2)
    𝒯 = SSHTRS(
        transform_indices(s, ℓₘₐₓ)..., T;
        θ=clenshaw_curtis_rings(Nθ, T), quadrature_weights=clenshaw_curtis(Nθ, T), Nϕ
    )
    # The worker warns about too few points on a ring.  The general check of the quadrature
    # rule, which the public `SSHTRS` constructor makes, is replaced here by the condition it
    # amounts to on this fixed grid, stated as a number of rings: 2⌊ℓₘₐₓ⌋+1 — 2ℓₘₐₓ+1 for an
    # integer ℓₘₐₓ, one fewer for a half-odd one (measured: exact there, and wrong by O(1) a
    # few rings below).  Too few rings make the analysis silently wrong, while synthesis is
    # exact on any number, so the warning is here rather than in `salm2map`.
    let needed = 2floor(Int, 𝒯.ℓₘₐₓ) + 1
        if Nθ < needed
            @warn (
                "The map has Nθ=$Nθ rings, but the Clenshaw–Curtis analysis needs at least "
                * "$needed for ℓₘₐₓ=$(𝒯.ℓₘₐₓ); the mode weights will not be exact."
            )
        end
    end
    𝒯
end

@doc raw"""
    salm2map(salm, s, ℓₘₐₓ, Nϕ, Nθ)
    salm2map(salm, 𝒯::SSHTRS)

Evaluate the spin-weighted function with mode weights `salm` on the equiangular ``N_ϕ ×
N_θ`` grid used by [`map2salm`](@ref).  The mode weights must be given in the canonical
ordering `ℓ ∈ abs(s):ℓₘₐₓ, m ∈ -ℓ:ℓ` along the first dimension, as in a
[`ModeWeights`](@ref).  The result has size ``N_ϕ × N_θ`` followed by the trailing
dimensions of `salm`.
"""
function salm2map end

function salm2map(salm::MapOrModes, s::IndexArgument, ℓₘₐₓ::IndexArgument, Nϕ::Integer, Nθ::Integer)
    T = real(eltype(salm))
    # Synthesis is exact on any number of rings, so the worker is called directly, without
    # the check of the quadrature rule that the public constructor makes for the analysis.
    𝒯 = SSHTRS(
        transform_indices(s, ℓₘₐₓ)..., T;
        θ=clenshaw_curtis_rings(Nθ, T), quadrature_weights=clenshaw_curtis(Nθ, T), Nϕ
    )
    salm2map(salm, 𝒯)
end
function salm2map(salm::MapOrModes, 𝒯::SSHTRS)
    Nθ = length(𝒯.θ)
    Nϕ = 𝒯.Nϕ[1]
    if 𝒯.Nϕ != fill(Nϕ, Nθ)
        error("salm2map requires the same number of points on every ring.")
    end
    check_clenshaw_curtis(𝒯, "salm2map")
    f = 𝒯 * salm
    reshape(f, Nϕ, Nθ, size(salm)[2:end]...)
end
