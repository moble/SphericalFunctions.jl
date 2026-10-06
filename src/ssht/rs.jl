"""
    SSHTRS(s, ℓₘₐₓ, [T=Float64]; [θ, quadrature_weights], Nϕ=2ℓₘₐₓ+1, plan_fft_flags=FFTW.ESTIMATE, plan_fft_timelimit=Inf)

Construct an ``s``-SHT object that uses the ring-based algorithm described by [Reinecke and
Seljebotn](@cite Reinecke_2013).  This may also be achieved by calling the main
[`SSHT`](@ref) function with the same keywords, along with `method="RS"` (the default).

The parameters of the type `SSHTRS{T, ST, P, BP, B, IT}` are as follows:
- `T` is the real type the transform works in.
- `ST` is the storage type of the ``H`` wedge in the transform's [`sλlmCalculator`](@ref).
- `P` and `BP` are the types of the forward and backward FFT plans.
- `B` is `true` when that calculator is batched (see [`isbatched`](@ref)).
- `IT` is the index type, `Int` or [`HalfOddInteger`](@ref).

The spin-weighted spherical harmonics are evaluated on a series of "rings" at constant
colatitude, whose locations are given by the `θ` keyword argument, and the analysis
integrates over ``θ`` with the `quadrature_weights` of the rule that placed those rings.
When both are omitted they are the nodes and weights of Fejér's first rule,
`fejer1_rings(2ℓₘₐₓ+1, T)` and `fejer1(2ℓₘₐₓ+1, T)`.  Only the caller knows which rule
placed a given set of rings, so the two must be given together — for example
`θ=clenshaw_curtis_rings(N, T)` with `quadrature_weights=clenshaw_curtis(N, T)` — and either
one without the other is refused.  The analysis is exact for band-limited functions when the
quadrature rule integrates polynomials of degree ``2ℓₘₐₓ`` in ``\\cos θ`` exactly, as the
Fejér and Clenshaw–Curtis rules with at least ``2ℓₘₐₓ+1`` nodes do.  The constructor checks
this, and warns when the rule falls short; synthesis does not use the weights, and is exact
on any rings.

On each ring, an FFT is performed.  To reach the band limit of ``m = ±ℓₘₐₓ``, the number of
points along each ring must be *at least* ``2ℓₘₐₓ+1``, but may be greater.  For example, if
``2ℓₘₐₓ+1`` does not factorize neatly into a product of small primes, it may be preferable
to use ``2ℓₘₐₓ+2`` points along each ring.  The number of points on each ring can be
modified independently, if given as a vector with the same length as `θ`, or as a single
number which is used for all rings.

The sample points are ordered ring by ring, with the azimuth ``ϕ_k = 2πk/N_ϕ`` for ``k = 0,
…, N_ϕ-1`` varying fastest; see [`pixels`](@ref) and [`rotors`](@ref).  The transform works
in the floating-point type `T`, and the colatitudes and weights must be given in that type,
or as integers.  It never acts in place: `𝒯 * f̃` and `𝒯 \\ f` allocate their results, and
real data are accepted; see [`SSHT`](@ref).

Whenever `T` is either `Float64` or `Float32`, the keyword arguments `plan_fft_flags` and
`plan_fft_timelimit` may also be useful for obtaining more efficient FFTs.  They default to
`FFTW.ESTIMATE` and `Inf`, respectively, and are passed to
[`AbstractFFTs.plan_fft!`](https://juliamath.github.io/AbstractFFTs.jl/stable/api/#AbstractFFTs.plan_fft).
One pair of plans is made for each distinct number of points on a ring, and each plan runs
on a single thread, since a ring is too short for FFTW's threads to be worth their cost.
For other element types the ring FFTs are computed by generic code, and for `BigFloat` they
are most of the cost of a transform; one transform still runs on a single thread, so the way
to use several threads is to give each task its own `copy(𝒯)`.

The harmonics on the rings are computed with one batched [`sλlmCalculator`](@ref), one ``ℓ``
at a time, so the cost is ``O(N_θ ℓₘₐₓ^2)`` and the memory ``O(N_θ ℓₘₐₓ)``.  Data of several
columns are transformed up to eight columns at a time, with the harmonics of each ``ℓ``
computed once for all of them, which makes each column about 1.5 to 2 times faster to
transform; the results of each column are exactly those of a transform of that column alone.
The object holds workspace for these computations, so two tasks must not use it at the same
time; `copy(𝒯)` returns another transform, which shares the rings, weights and FFT plans of
`𝒯` and has workspace of its own.  See [`SSHT`](@ref).

Half-integer `s` and `ℓₘₐₓ` are accepted as `Rational`s with denominator 2; see
[`SSHT`](@ref) for what the function values then mean.  The algorithm is the same one: the
harmonics on a ring are ``i^{2s}`` times real functions of ``θ``, so the ``θ`` stage works
with those real functions and restores the phase once per ring, and ``e^{imϕ}`` with
half-odd ``m`` is ``e^{iϕ/2}`` times a Fourier mode of integer frequency ``m - 1/2``, so the
FFT runs at integer frequencies and each sample of a ring is multiplied by ``e^{±iϕ/2}``.
The defaults are ``2ℓₘₐₓ+1`` rings of ``2ℓₘₐₓ+1`` points, as for integer indices, and these
are then even numbers.
"""
struct SSHTRS{T<:Real, ST, P, BP, B, IT<:IntegerHalf} <: SSHT{T}
    s::IT
    ℓₘₐₓ::IT
    θ::Vector{T}
    quadrature_weights::Vector{T}
    Nϕ::Vector{Int}
    ring_ranges::Vector{UnitRange{Int}}  # index range of each ring in the pixel vector
    plans::RingPlans{T, P, BP}  # FFT plans for each distinct ring size (see `RingPlans`)
    synthesis_phases::Vector{Vector{Complex{T}}}  # for each ring size (see `ring_phases`)
    analysis_phases::Vector{Vector{Complex{T}}}
    # Workspace, which `copy` allocates afresh.  `F` holds, for each column of a chunk (see
    # `rs_chunk`), the Fourier coefficients [ring, m+ℓₘₐₓ+1] for m ∈ -ℓₘₐₓ:ℓₘₐₓ.
    λ::HarmonicCalculator{IT, T, T, ST, IT, B, T, Nothing}  # ₛλₗₘ(θ) for all rings at once
    F::Vector{Matrix{Complex{T}}}
    G::Vector{Vector{Complex{T}}}  # an FFT buffer for each distinct ring size
end

# The public constructor checks the quadrature rule and the number of points on each ring,
# once the transform is built.  `map2salm_plan` and `salm2map` build their transforms with
# `rs_transform` directly: they build the Clenshaw–Curtis rule themselves, and the first
# states the number of rings that rule needs more plainly than the general check could,
# while the second only synthesizes, which is exact on any grid.
@index_methods function SSHTRS(
    s::IndexType, ℓₘₐₓ::IndexType, ::Type{TT}=Float64;
    θ=nothing,
    quadrature_weights=nothing,
    Nϕ=2ℓₘₐₓ+1,
    plan_fft_flags=FFTW.ESTIMATE, plan_fft_timelimit=Inf
) where {TT}
    𝒯 = rs_transform(
        s, ℓₘₐₓ, TT; θ, quadrature_weights, Nϕ, plan_fft_flags, plan_fft_timelimit
    )
    warn_if_aliased(𝒯)
    warn_if_inexact(𝒯)
    𝒯
end

function rs_transform(
    s::IT, ℓₘₐₓ::IT, ::Type{TT};
    θ=nothing,
    quadrature_weights=nothing,
    Nϕ=2ℓₘₐₓ+1,
    plan_fft_flags=FFTW.ESTIMATE, plan_fft_timelimit=Inf
) where {IT<:IntegerHalf, TT}
    check_transform_type(TT)
    check_band_limit(s, ℓₘₐₓ)
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
        throw(DimensionMismatch(
            "θ and quadrature_weights must have the same length; got $(length(θ)) and "
            * "$(length(quadrature_weights))."
        ))
    end
    Nϕ = Nϕ isa Integer ? fill(Int(Nϕ), length(θ)) : Vector{Int}(Nϕ)
    if length(Nϕ) != length(θ)
        throw(DimensionMismatch(
            "Nϕ must be a single number or have the same length as θ ($(length(θ)))."
        ))
    end
    if any(<(1), Nϕ)
        throw(ArgumentError("Every ring needs at least one point; got Nϕ=$Nϕ."))
    end
    Nθ = length(θ)
    ring_ranges = let stops = cumsum(Nϕ)
        [(stop - n + 1):stop for (n, stop) ∈ zip(Nϕ, stops)]
    end
    plans = ring_plans(TT, Nϕ; flags=plan_fft_flags, timelimit=plan_fft_timelimit)
    synthesis_phases, analysis_phases = ring_phases(TT, plans.sizes, s)
    λ = sλlmCalculator(θ, ℓₘₐₓ, s)  # θ is already a Vector{TT}
    F = [Matrix{Complex{TT}}(undef, Nθ, 2ℓₘₐₓ + 1)]
    G = [Vector{Complex{TT}}(undef, N) for N ∈ plans.sizes]
    P, BP = eltype(plans.forward), eltype(plans.backward)
    SSHTRS{TT, typeof(parent(λ.engine.H.Hˡ)), P, BP, isbatched(λ), IT}(
        s, ℓₘₐₓ, θ, quadrature_weights, Nϕ, ring_ranges, plans,
        synthesis_phases, analysis_phases, λ, F, G
    )
end

# An independent transform, for use by another task: the rings, weights, plans and phases
# are shared, since no transform modifies them, and the workspace is new.
function Base.copy(𝒯::SSHTRS{T, ST, P, BP, B, IT}) where {T, ST, P, BP, B, IT}
    SSHTRS{T, ST, P, BP, B, IT}(
        𝒯.s, 𝒯.ℓₘₐₓ, 𝒯.θ, 𝒯.quadrature_weights, 𝒯.Nϕ, 𝒯.ring_ranges, 𝒯.plans,
        𝒯.synthesis_phases, 𝒯.analysis_phases,
        similar(𝒯.λ), [similar(𝒯.F[1])], [similar(g) for g ∈ 𝒯.G]
    )
end

# A ring of fewer than 2ℓₘₐₓ+1 points cannot tell the frequencies m = ±ℓₘₐₓ apart, so the
# analysis of a band-limited function on it is not exact.  Synthesis is exact on any ring,
# so `salm2map`, which only synthesizes, does not warn.
function warn_if_aliased(𝒯::SSHTRS)
    if any(<(2𝒯.ℓₘₐₓ+1), 𝒯.Nϕ)
        @warn (
            "Some rings have fewer than 2ℓₘₐₓ+1=$(2𝒯.ℓₘₐₓ+1) points, so modes with large |m| "
            * "will alias, and analysis on these rings will not be exact."
        )
    end
    nothing
end

# The analysis integrates, ring by ring, products ₛλₗₘ ₛλₗ′ₘ of harmonics with the same m,
# and each such product is a polynomial in cos θ of degree ℓ+ℓ′ ≤ 2ℓₘₐₓ, for half-odd
# indices as for integers.  The analysis is therefore exact when the quadrature rule
# integrates every polynomial of that degree exactly, which is checked here with the
# Legendre moments: the sum Σ_y w_y P_k(cos θ_y) must be 2δ_{k0} for each k ∈ 0:2ℓₘₐₓ.
# Unlike the monomials, the P_k are bounded by 1 on the whole interval, so the moments are
# well conditioned, and a fixed multiple of the rounding error of the sums serves as the
# tolerance.  The degree is 2ℓₘₐₓ itself rather than 2⌊ℓₘₐₓ⌋, so that a rule without the
# reflection symmetry of the standard ones is checked in every degree the analysis needs; a
# symmetric rule integrates the odd degrees exactly in any case, which is why 2ℓₘₐₓ of its
# rings suffice for a half-odd ℓₘₐₓ.
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

pixels(𝒯::SSHTRS) = ring_pixels(𝒯.θ, 𝒯.Nϕ)
rotors(𝒯::SSHTRS) = from_spherical_coordinates.(pixels(𝒯))
npixels(𝒯::SSHTRS) = 𝒯.ring_ranges[end].stop

function Base.:*(𝒯::SSHTRS, f̃::SSHTData)
    d = synthesis_modes(𝒯, f̃)
    mul!(pixel_output(𝒯, d), 𝒯, d)
end

# Several columns are transformed in chunks, so that the harmonics ₛλₗₘ of each ℓ are
# computed once for a whole chunk rather than once for each column: the recurrence is about
# a third of the cost of transforming a column.  Each column of a chunk needs a matrix of
# Fourier coefficients, Nθ × (2ℓₘₐₓ+1), and these are allocated when a chunk first needs
# them, and kept.  A chunk has at most 8 columns, and fewer where the matrices are large, so
# that they hold no more than about 2²² complex numbers in all (64 MiB in Float64); in
# particular a single column never allocates.  For each column the arithmetic is exactly
# that of a transform of the column alone, in the same order, so the results are the same to
# the last bit.
function rs_chunk(𝒯::SSHTRS, ncolumns)
    K = max(1, min(ncolumns, 8, 2^22 ÷ max(1, length(𝒯.F[1]))))
    while length(𝒯.F) < K
        push!(𝒯.F, similar(𝒯.F[1]))
    end
    K
end

# Synthesis: f = 𝒯 * f̃
function LinearAlgebra.mul!(f, 𝒯::SSHTRS, f̃)
    check_modes(𝒯, f̃)
    check_pixels(𝒯, f)
    check_trailing(f, f̃)
    check_complex_output(f, "f")
    f̃′ = reshape(array_view(f̃), size(f̃, 1), :)
    f′ = reshape(array_view(f), size(f, 1), :)
    for columns ∈ Iterators.partition(axes(f′, 2), rs_chunk(𝒯, size(f′, 2)))
        rs_synthesis!(view(f′, :, columns), 𝒯, view(f̃′, :, columns))
    end
    f
end

# The Fourier coefficients on each ring, F[y, m] = Σ_ℓ f̃ₗₘ ₛλₗₘ(θ_y), for each column of
# the chunk, and then the inverse Fourier transform on each ring (aliasing the m values if
# Nϕ < 2ℓₘₐₓ+1)
function rs_synthesis!(f, 𝒯::SSHTRS{T}, f̃) where {T}
    s, ℓₘₐₓ, Nθ = 𝒯.s, 𝒯.ℓₘₐₓ, length(𝒯.θ)
    λ, F, G, plans = 𝒯.λ, 𝒯.F, 𝒯.G, 𝒯.plans
    columns = axes(f̃, 2)
    for k ∈ columns
        fill!(F[k], zero(Complex{T}))
    end
    for ℓ ∈ abs(s):ℓₘₐₓ
        Λ = array_view(recurrence!(λ, ℓ))  # [y, m+ℓ+1]
        i₀ = Yindex(ℓ, -ℓ, abs(s)) - 1
        @inbounds for k ∈ columns
            Fₖ = F[k]
            for m ∈ -ℓ:ℓ
                f̃ₗₘ = f̃[i₀ + ℓ + m + 1, k]
                jm = m + ℓₘₐₓ + 1
                jℓ = m + ℓ + 1
                @simd for y ∈ 1:Nθ
                    Fₖ[y, jm] += f̃ₗₘ * Λ[y, jℓ]
                end
            end
        end
    end
    @inbounds for k ∈ columns
        Fₖ = F[k]
        fₖ = view(f, :, k)
        for y ∈ 1:Nθ
            j = plans.index[y]
            Gy, Nϕy = G[j], 𝒯.Nϕ[y]
            fill!(Gy, zero(Complex{T}))
            for m ∈ -ℓₘₐₓ:ℓₘₐₓ
                Gy[1 + mod(floor_int(m), Nϕy)] += Fₖ[y, m + ℓₘₐₓ + 1]
            end
            plans.backward[j] * Gy  # unnormalized inverse FFT: Σₘ Gₘ e^{+imϕₖ}
            ring_values!(view(fₖ, 𝒯.ring_ranges[y]), Gy, 𝒯.synthesis_phases[j], s)
        end
    end
    f
end

function Base.:\(𝒯::SSHTRS, f::SSHTData)
    check_pixels(𝒯, f)
    ldiv!(mode_output(𝒯, f), 𝒯, f)
end

# Analysis: f̃ = 𝒯 \ f
function LinearAlgebra.ldiv!(f̃, 𝒯::SSHTRS, f)
    f̃ = analysis_output(𝒯, f̃, f)
    check_modes(𝒯, f̃)
    check_pixels(𝒯, f)
    check_trailing(f, f̃)
    check_complex_output(f̃, "f̃")
    f̃′ = reshape(array_view(f̃), size(f̃, 1), :)
    f′ = reshape(array_view(f), size(f, 1), :)
    for columns ∈ Iterators.partition(axes(f′, 2), rs_chunk(𝒯, size(f′, 2)))
        rs_analysis!(view(f̃′, :, columns), 𝒯, view(f′, :, columns))
    end
    f̃
end

# The Fourier transform on each ring of each column of the chunk, including the quadrature
# weight and the ϕ measure, F[y, m] = w_y (2π/Nϕ) Σₖ f(θ_y, ϕₖ) e^{-imϕₖ}, and then the mode
# weights f̃ₗₘ = Σ_y F[y, m] ₛλₗₘ(θ_y)
function rs_analysis!(f̃, 𝒯::SSHTRS{T}, f) where {T}
    s, ℓₘₐₓ, Nθ = 𝒯.s, 𝒯.ℓₘₐₓ, length(𝒯.θ)
    λ, F, G, plans = 𝒯.λ, 𝒯.F, 𝒯.G, 𝒯.plans
    columns = axes(f, 2)
    twoπ = 2T(π)
    @inbounds for k ∈ columns
        Fₖ = F[k]
        fₖ = view(f, :, k)
        for y ∈ 1:Nθ
            j = plans.index[y]
            Gy, Nϕy = G[j], 𝒯.Nϕ[y]
            factor = 𝒯.quadrature_weights[y] * twoπ / Nϕy
            ring_samples!(Gy, view(fₖ, 𝒯.ring_ranges[y]), factor, 𝒯.analysis_phases[j], s)
            plans.forward[j] * Gy
            for m ∈ -ℓₘₐₓ:ℓₘₐₓ
                Fₖ[y, m + ℓₘₐₓ + 1] = Gy[1 + mod(floor_int(m), Nϕy)]
            end
        end
    end
    # The sum over rings is formed in order, without `@simd`.  `@simd` would permit the
    # compiler to reassociate it, and whether it does depends on the context in which each
    # specialization of this function happens to be compiled, so that the same data analyzed
    # from real and from complex storage, or into different outputs, could differ in the
    # last bit.  The compiler still vectorizes the loop, with a reduction kept in order.
    for ℓ ∈ abs(s):ℓₘₐₓ
        Λ = array_view(recurrence!(λ, ℓ))  # [y, m+ℓ+1]
        i₀ = Yindex(ℓ, -ℓ, abs(s)) - 1
        @inbounds for k ∈ columns
            Fₖ = F[k]
            for m ∈ -ℓ:ℓ
                jm = m + ℓₘₐₓ + 1
                jℓ = m + ℓ + 1
                acc = zero(Complex{T})
                for y ∈ 1:Nθ
                    acc += Fₖ[y, jm] * Λ[y, jℓ]
                end
                f̃[i₀ + ℓ + m + 1, k] = acc
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

The `map` array, of real or complex numbers, must have size ``N_ϕ`` along its first
dimension and ``N_θ ≥ 2`` along its second, with any number of dimensions following, sampled
at the Clenshaw–Curtis nodes [`clenshaw_curtis_rings`](@ref)`(Nθ)` in ``θ`` (which include
both poles) and at ``ϕ_k = 2πk/N_ϕ``; [`pixels`](@ref) returns that grid for a constructed
[`SSHTRS`](@ref).  For the analysis to be exact for a band-limited function, one needs ``N_ϕ
≥ 2ℓₘₐₓ+1`` points on each ring and ``N_θ ≥ 2⌊ℓₘₐₓ⌋+1`` rings — that is, ``2ℓₘₐₓ+1`` rings
for an integer ``ℓₘₐₓ``, and ``2ℓₘₐₓ`` for a half-odd one; with fewer, the result is not
exact, and a warning is issued.  The transform works in the floating-point type of the map's
elements, `float(real(eltype(map)))`.

The result is a [`ModeWeights`](@ref) for a one-dimensional map (``N_ϕ × N_θ``), or an array
whose first dimension indexes the modes in the canonical ordering `ℓ ∈ abs(s):ℓₘₐₓ, m ∈
-ℓ:ℓ` otherwise.  Note that the output starts at ``ℓ = |s|``, not at ``ℓ = 0``.

For repeated use with different `map`s of the same shape, construct the transform once with
[`SphericalFunctions.map2salm_plan`](@ref)`(map, s, ℓₘₐₓ)` (an [`SSHTRS`](@ref) on the
Clenshaw–Curtis rings) and pass it as the second argument.  See also [`salm2map`](@ref).

The spin weight and ``ℓₘₐₓ`` may be half-integers, passed as `Rational`s with denominator 2.
The map is then antiperiodic in ``ϕ`` — its values are those of the function at the rotors
`from_spherical_coordinates(θ, ϕ)`, and a full circuit of the azimuth reaches the antipodal
rotor — and each ring needs ``N_ϕ ≥ 2ℓₘₐₓ+1`` points, as for an integer ``ℓₘₐₓ``; see
[`SSHT`](@ref).
"""
function map2salm end

@index_methods map2salm(map::MapArray, s::IndexType, ℓₘₐₓ::IndexType) =
    map2salm(map, map2salm_plan(map, s, ℓₘₐₓ))
function map2salm(map::MapArray, 𝒯::SSHTRS)
    Nϕ, Nθ = size(map, 1), size(map, 2)
    if 𝒯.Nϕ != fill(Nϕ, Nθ) || length(𝒯.θ) != Nθ
        throw(DimensionMismatch(
            "The transform was planned for a different grid than the $(Nϕ)×$(Nθ) map."
        ))
    end
    check_clenshaw_curtis(𝒯, "map2salm")
    f = reshape(map, Nϕ * Nθ, size(map)[3:end]...)
    𝒯 \ f
end

# `map2salm` and `salm2map` promise the Clenshaw–Curtis grid, but the shape alone does not
# establish it: any `SSHTRS` with the right numbers of rings and points has that shape, and
# in particular the default `SSHT(s, ℓₘₐₓ)`, on the Fejér rings, has exactly the shape of
# the smallest such grid.  The rings and weights themselves are therefore compared.  (The
# rings are compared first, so the weights are not computed for a degenerate single ring.)
function check_clenshaw_curtis(𝒯::SSHTRS{T}, name) where {T}
    Nθ = length(𝒯.θ)
    if !(𝒯.θ ≈ clenshaw_curtis_rings(Nθ, T)) || !(𝒯.quadrature_weights ≈ clenshaw_curtis(Nθ, T))
        throw(ArgumentError(
            "$name works on the Clenshaw–Curtis grid, but this transform has other rings or "
            * "quadrature weights; construct it with `map2salm_plan(map, s, ℓₘₐₓ)`."
        ))
    end
end

# The Clenshaw–Curtis rings of the equiangular grid include both poles, so there must be two
# of them at least; and a map must have the two dimensions of the grid, which is a matter of
# its size, and so a `DimensionMismatch`, as it is for `map2salm(map, 𝒯)`.
function check_equiangular_grid(name, Nθ, shape=nothing)
    if shape !== nothing && length(shape) < 2
        throw(DimensionMismatch(
            "`$name` takes a map of size Nϕ × Nθ, with any further dimensions following, but "
            * "this map has size $shape."
        ))
    end
    if Nθ < 2
        throw(ArgumentError(
            "`$name` needs at least two rings of the Clenshaw–Curtis grid, one at each pole; "
            * "got Nθ=$Nθ."
        ))
    end
    nothing
end

"""
    map2salm_plan(map, s, ℓₘₐₓ)

Construct the [`SSHTRS`](@ref) transform used by [`map2salm`](@ref) and [`salm2map`](@ref)
for maps of the shape of `map` (``N_ϕ × N_θ × …``) on the Clenshaw–Curtis grid, in the
floating-point type of the map's elements.
"""
@index_methods function map2salm_plan(map::MapArray, s::IndexType, ℓₘₐₓ::IndexType)
    T = float(real(eltype(map)))
    Nϕ, Nθ = size(map, 1), size(map, 2)
    check_equiangular_grid("map2salm", Nθ, size(map))
    𝒯 = rs_transform(
        s, ℓₘₐₓ, T;
        θ=clenshaw_curtis_rings(Nθ, T), quadrature_weights=clenshaw_curtis(Nθ, T), Nϕ
    )
    warn_if_aliased(𝒯)
    # The general check of the quadrature rule, which the public `SSHTRS` constructor makes,
    # is replaced here by the condition it amounts to on this fixed grid, stated as a number
    # of rings: 2⌊ℓₘₐₓ⌋+1 — 2ℓₘₐₓ+1 for an integer ℓₘₐₓ, one fewer for a half-odd one
    # (measured: exact there, and wrong by O(1) a few rings below).  Too few rings make the
    # analysis silently wrong, while synthesis is exact on any number, so the warning is
    # here rather than in `salm2map`.
    let needed = 2floor_int(𝒯.ℓₘₐₓ) + 1
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
    salm2map(w::ModeWeights, Nϕ, Nθ)
    salm2map(salm, 𝒯::SSHTRS)

Evaluate the spin-weighted function with mode weights `salm` on the equiangular ``N_ϕ ×
N_θ`` grid used by [`map2salm`](@ref).  The mode weights must be given in the canonical
ordering `ℓ ∈ abs(s):ℓₘₐₓ, m ∈ -ℓ:ℓ` along the first dimension, as real or complex numbers,
or as a [`ModeWeights`](@ref), which may then cover any range of ``ℓ`` up to ``ℓₘₐₓ`` (see
[`SSHT`](@ref)).  Given a `ModeWeights` `w`, the spin weight and ``ℓₘₐₓ`` may be omitted,
and are then those of `w`.  The result has size ``N_ϕ × N_θ`` followed by the trailing
dimensions of `salm`, and is computed in the floating-point type of the elements of `salm`.

Synthesis is exact on any grid of ``N_θ ≥ 2`` rings, so no condition on ``N_ϕ`` or ``N_θ``
is needed beyond that; but a map with fewer than ``2ℓₘₐₓ+1`` points on each ring, or fewer
rings than [`map2salm`](@ref) needs, cannot be analyzed back to the same mode weights.
"""
function salm2map end

@index_methods function salm2map(
    salm::ModesArray, s::IndexType, ℓₘₐₓ::IndexType, Nϕ::Integer, Nθ::Integer
)
    T = float(real(eltype(salm)))
    check_equiangular_grid("salm2map", Nθ)
    𝒯 = rs_transform(
        s, ℓₘₐₓ, T;
        θ=clenshaw_curtis_rings(Nθ, T), quadrature_weights=clenshaw_curtis(Nθ, T), Nϕ
    )
    salm2map(salm, 𝒯)
end
salm2map(w::ModeWeights, Nϕ::Integer, Nθ::Integer) = salm2map(w, w.s, w.ℓₘₐₓ, Nϕ, Nθ)
function salm2map(salm::ModesArray, 𝒯::SSHTRS)
    Nθ = length(𝒯.θ)
    Nϕ = 𝒯.Nϕ[1]
    if 𝒯.Nϕ != fill(Nϕ, Nθ)
        throw(ArgumentError("salm2map requires the same number of points on every ring."))
    end
    check_clenshaw_curtis(𝒯, "salm2map")
    f = 𝒯 * salm
    reshape(f, Nϕ, Nθ, size(salm)[2:end]...)
end
