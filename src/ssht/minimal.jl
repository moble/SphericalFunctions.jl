# One group of modes solved together: the modes of a set of m values that alias into one
# another (a strongly connected component of the graph described in `minimal_blocks`), and
# the Fourier coefficients — `(ring, index)` pairs, with index `mod(m, Nϕ)+1` — that
# determine them.  There are exactly as many coefficients as modes.  `T` is the number type
# of the LU decomposition of their system.
struct MinimalBlock{T}
    modes::Vector{Int}
    coefficients::Vector{Tuple{Int, Int}}
    lu::LinearAlgebra.LU{T, Matrix{T}, Vector{Int}}
end

"""
    SSHTMinimal(s, ℓₘₐₓ, [T=Float64]; θ=minimal_rings(s, ℓₘₐₓ, T).θ, plan_fft_flags=FFTW.ESTIMATE, plan_fft_timelimit=Inf, inplace=true)

Construct an ``s``-SHT object that uses the optimal-dimensionality algorithm described by
[Elahi et al.](@cite Elahi_2018), which samples the function at exactly as many points as
there are modes.  This may also be achieved by calling the main [`SSHT`](@ref) function with
the same keywords, along with `method="Minimal"`.

The parameters of the type `SSHTMinimal{T, Inplace, P, BP}` are as follows:
- `T` is the real type the transform works in.
- `Inplace` is `true` when the transforms act in place, as set by the `inplace` keyword.
- `P` and `BP` are the types of the forward and backward FFT plans.

!!! warning
    This method is experimental and not very accurate.  The round-trip error grows
    exponentially with ``ℓₘₐₓ``, and the constructor warns when fewer than half the digits
    of `T` would survive.  The `"RS"` method (the default) has no such limitation, but does
    not have optimal dimensionality.  The `"Matrix"` method does have optimal
    dimensionality, but its memory consumption scales poorly.

The function is sampled on ``ℓₘₐₓ-|s|+1`` "rings" at constant colatitude, each holding an odd
number of equally spaced points starting at ``ϕ = 0``.  Their sizes, and the default
colatitudes, are given by [`minimal_rings`](@ref); the `θ` keyword argument may give other
colatitudes, which must be distinct, one for each ring in the order listed there, which is
also the order of the sample points, given in the type `T` that the transform works in, or
as integers.  See [`pixels`](@ref) and [`rotors`](@ref) for the sample points themselves.

For ``s = 0`` the rings have ``1, 3, …, 2ℓₘₐₓ+1`` points.  For any other spin weight that
choice is badly conditioned: near the north pole a function of spin weight ``s`` is
dominated by the modes with ``m`` near ``-s``, and near the south pole by those near ``+s``,
so that the smallest rings — which are placed nearest the poles — see some of the
frequencies they are responsible for more weakly than the higher frequencies that alias onto
them.  The error of a round trip then grows by more than an order of magnitude with each
unit of ℓₘₐₓ, and no choice of colatitudes cures it.  For ``s ≠ 0`` the rings are therefore
arranged so that each polar ring is centered, in frequency, on the modes that dominate near
its pole; see [`minimal_rings`](@ref).  The analysis is then no longer a sequence of solves
for one ``m`` at a time, but for small groups of ``m`` values that alias into one another —
in every case measured, at most ``4|s|-1`` of them, however large ℓₘₐₓ is.

Whenever `T` is either `Float64` or `Float32`, the keyword arguments `plan_fft_flags` and
`plan_fft_timelimit` may also be useful for obtaining more efficient FFTs.  They default to
`FFTW.ESTIMATE` and `Inf`, respectively, and are passed to
[`AbstractFFTs.plan_fft!`](https://juliamath.github.io/AbstractFFTs.jl/stable/api/#AbstractFFTs.plan_fft).
One pair of plans is made for each distinct number of points on a ring, and each plan runs
on a single thread.

Because this algorithm achieves optimal dimensionality, the transformation is performed in
place by default, in both directions: `𝒯 * f̃` overwrites the storage of `f̃` with the
function values (and returns that storage, as a plain array), and `𝒯 \\ f` overwrites `f`
with the mode weights (and returns them, for one-dimensional `f`, as a `ModeWeights`
wrapping that storage).  That storage must hold complex floating-point numbers.  If this is
not desired, pass the keyword argument `inplace=false`, which makes those operations work on
a copy of the input, in `Complex{T}`, and so accept real data.  Whatever the option, the
two-argument `LinearAlgebra.mul!(𝒯, x)` and `LinearAlgebra.ldiv!(𝒯, x)` act in place, the
first returning the function values as a plain array and the second the mode weights as `𝒯
\\ x` would.  See [`SSHT`](@ref).

The values ``{}_sλ_{ℓ,m}(θ_r)`` of every mode on every ring are precomputed at construction
(with one batched [`sλlmCalculator`](@ref)) and stored, which takes ``O(ℓₘₐₓ^3)`` memory, as
are the LU decompositions of the matrices for the groups of ``m`` values.  The object also
holds workspace for the transforms, so two tasks must not use it at the same time;
`copy(𝒯)` returns another transform, which shares those tables and the FFT plans of `𝒯`
and has workspace of its own.  See [`SSHT`](@ref).

Even so, the sample points become increasingly badly conditioned as ℓₘₐₓ grows, for every
spin weight: in `Float64` a round trip loses about 3.5 digits by ℓₘₐₓ = 32, 5.5 by ℓₘₐₓ = 48
and 8 by ℓₘₐₓ = 64 at ``s = 0``, and more at larger ``|s|`` — about 9 by ℓₘₐₓ = 32 at ``s =
2``.  The constructor therefore measures the error of one round trip, and warns when fewer
than half the digits of `T` would survive; the `"RS"` method has no such limitation.

This method is defined only for integer spin weights.  Its bookkeeping — rings of an odd
number of points, and the aliasing of ``m`` into rings too small to hold it — is written for
integer indices throughout, and has not been extended; a half-integer spin weight is
rejected with a message naming the two methods, `"RS"` and `"Matrix"`, that do accept one.
"""
struct SSHTMinimal{T<:Real, Inplace, P, BP} <: SSHT{T}
    s::Int
    ℓₘₐₓ::Int
    θ::Vector{T}  # colatitude of each ring
    Nϕ::Vector{Int}  # number of points on each ring (odd)
    centers::Vector{Int}  # center of each ring's window of frequencies (see `minimal_rings`)
    ring_ranges::Vector{UnitRange{Int}}  # pixel indices of each ring
    plans::RingPlans{T, P, BP}  # FFT plans for each distinct ring size (see `RingPlans`)
    mode_m::Vector{Int}  # m of each mode, in the canonical order
    Λ::Matrix{T}  # Λ[i, r] = ₛλ_{ℓ,m}(θ_r) for the mode (ℓ, m) with index i
    blocks::Vector{MinimalBlock{T}}  # the groups of modes solved together, in solution order
    # Workspace, which `copy` allocates afresh
    F::Vector{Vector{Complex{T}}}  # Fourier coefficients of each ring
    f̃::Vector{Complex{T}}  # a copy of the mode weights, for synthesis
    rhs::Vector{Complex{T}}  # the right-hand side of one block's system
end


# The groups of m values that must be solved together, in an order in which each group's
# aliases into the others' rings are known by the time they are needed.  There is an edge m′
# → m whenever some ring's window includes m but not m′ ≡ m (mod Nϕ), since m′ then aliases
# into the coefficient that measures m, and must be removed from it first.  The groups are
# the strongly connected components of that graph, which Kosaraju's algorithm produces in
# topological order (sources first).  For windows all centered on 0 — in particular for s =
# 0 — every group is a single m, and the order is that of decreasing |m|.
function minimal_blocks(ℓₘₐₓ, Nϕ, centers)
    ms = -ℓₘₐₓ:ℓₘₐₓ
    n = length(ms)
    index(m) = m + ℓₘₐₓ + 1
    successors = [Int[] for _ ∈ 1:n]
    predecessors = [Int[] for _ ∈ 1:n]
    for (N, c) ∈ zip(Nϕ, centers)
        k = N ÷ 2
        for m′ ∈ ms
            if abs(m′ - c) > k
                m = c + mod(m′ - c + k, N) - k  # the frequency in the window that m′ aliases to
                push!(successors[index(m′)], index(m))
                push!(predecessors[index(m)], index(m′))
            end
        end
    end
    # First pass: order the vertices by the time a depth-first search finishes with them
    finished = Int[]
    visited = falses(n)
    for root ∈ 1:n
        visited[root] && continue
        visited[root] = true
        stack = [(root, 1)]
        while !isempty(stack)
            v, i = stack[end]
            if i ≤ length(successors[v])
                stack[end] = (v, i + 1)
                w = successors[v][i]
                if !visited[w]
                    visited[w] = true
                    push!(stack, (w, 1))
                end
            else
                pop!(stack)
                push!(finished, v)
            end
        end
    end
    # Second pass: search the reversed graph in order of decreasing finishing time; each
    # search finds one component, and they are found in topological order
    component = zeros(Int, n)
    groups = Vector{Int}[]
    for root ∈ Iterators.reverse(finished)
        component[root] ≠ 0 && continue
        push!(groups, Int[])
        component[root] = length(groups)
        stack = [root]
        while !isempty(stack)
            v = pop!(stack)
            push!(groups[end], ms[v])
            for w ∈ predecessors[v]
                if component[w] == 0
                    component[w] = length(groups)
                    push!(stack, w)
                end
            end
        end
    end
    groups
end

# Only integer indices reach the body; a half-integer spin weight is refused with the
# message of `@index_methods`, followed by the methods that do accept one.
@index_methods integer_only (
    "The \"RS\" method (the default) and the \"Matrix\" method both accept half-integer "
    * "indices."
) function SSHTMinimal(
    s::IndexType, ℓₘₐₓ::IndexType, ::Type{TT}=Float64;
    θ=nothing,
    plan_fft_flags=FFTW.ESTIMATE, plan_fft_timelimit=Inf,
    inplace=true
) where {TT}
    check_transform_type(TT)
    check_band_limit(s, ℓₘₐₓ)
    check_inplace_option(inplace)
    rings = minimal_rings(s, ℓₘₐₓ, TT)
    nrings = length(rings.Nϕ)
    if θ === nothing
        θ = rings.θ
    elseif length(θ) != nrings
        throw(DimensionMismatch(
            "Length of θ ($(length(θ))) must equal ℓₘₐₓ-abs(s)+1 ($nrings)."
        ))
    end
    check_sample_reals(TT, θ, "θ")
    θ = Vector{TT}(θ)
    # Two rings at one colatitude make the system for the modes they share singular.  The
    # check of each LU decomposition below would catch that too, but not name the cause.
    let θsorted = sort(θ)
        i = findfirst(i -> θsorted[i] == θsorted[i+1], 1:nrings-1)
        if i !== nothing
            throw(ArgumentError(
                "The colatitudes θ must be distinct, but θ = $(θsorted[i])\n"
                * "appears more than once.  Two rings at one colatitude are degenerate."
            ))
        end
    end
    Nϕ, centers = rings.Nϕ, rings.centers

    ring_ranges = let stops = cumsum(Nϕ)
        [(stop - N + 1):stop for (N, stop) ∈ zip(Nϕ, stops)]
    end
    F = [Vector{Complex{TT}}(undef, N) for N ∈ Nϕ]
    plans = ring_plans(TT, Nϕ; flags=plan_fft_flags, timelimit=plan_fft_timelimit)

    # Tables of ₛλₗₘ(θᵣ) for every mode on every ring, computed with one batched calculator
    # in angle mode
    n = Ysize(abs(s), ℓₘₐₓ)
    mode_m = Vector{Int}(undef, n)
    Λ = Matrix{TT}(undef, n, nrings)
    λ = sλlmCalculator(θ, ℓₘₐₓ, s)  # θ is already a Vector{TT}, which fixes the type
    for ℓ ∈ abs(s):ℓₘₐₓ
        Λℓ = array_view(recurrence!(λ, ℓ))  # [ring, m+ℓ+1]
        for m ∈ -ℓ:ℓ
            i = Yindex(ℓ, m, abs(s))
            mode_m[i] = m
            for r ∈ 1:nrings
                Λ[i, r] = Λℓ[r, m + ℓ + 1]
            end
        end
    end

    # The matrix of each group couples its modes to the coefficients of its m values on
    # every ring whose window includes them; a mode contributes to a coefficient of a ring
    # whenever its m is congruent to that coefficient's m modulo the size of the ring.
    blocks = map(minimal_blocks(ℓₘₐₓ, Nϕ, centers)) do group
        modes = [i for m ∈ group for ℓ ∈ max(abs(s), abs(m)):ℓₘₐₓ for i ∈ Yindex(ℓ, m, abs(s))]
        coefficients = [
            (r, mod(m, Nϕ[r]) + 1)
            for m ∈ group for r ∈ 1:nrings if abs(m - centers[r]) ≤ Nϕ[r] ÷ 2
        ]
        if length(coefficients) != length(modes)
            error(  # Cannot happen for the layout of `minimal_rings`; a guard for changes to it
                "Internal error: the m values $group have $(length(modes)) modes but "
                * "$(length(coefficients)) Fourier coefficients."
            )
        end
        M = zeros(TT, length(coefficients), length(modes))
        for (e, (r, q)) ∈ enumerate(coefficients), (u, i) ∈ enumerate(modes)
            if mod(mode_m[i], Nϕ[r]) + 1 == q
                M[e, u] = Λ[i, r]
            end
        end
        factorization = LinearAlgebra.lu(M; check=false)
        if !LinearAlgebra.issuccess(factorization)
            throw(ArgumentError(
                "The colatitudes θ = $θ do not determine the modes with m ∈ $group:\n"
                * "their system is singular.  This happens, for example, when a ring of more "
                * "than one point lies at a pole, where all of its points coincide.  The "
                * "colatitudes given by `minimal_rings` avoid this."
            ))
        end
        MinimalBlock{TT}(modes, coefficients, factorization)
    end
    rhs = Vector{Complex{TT}}(undef, maximum(b -> length(b.modes), blocks))

    𝒯 = SSHTMinimal{TT, inplace, eltype(plans.forward), eltype(plans.backward)}(
        s, ℓₘₐₓ, θ, Nϕ, centers, ring_ranges, plans, mode_m, Λ, blocks, F,
        Vector{Complex{TT}}(undef, n), rhs
    )
    warn_if_inaccurate(𝒯, "Minimal")
    𝒯
end

# An independent transform, for use by another task: the rings, plans, tables and
# decompositions are shared, since no transform modifies them, and the workspace is new.
function Base.copy(𝒯::SSHTMinimal{T, Inplace, P, BP}) where {T, Inplace, P, BP}
    SSHTMinimal{T, Inplace, P, BP}(
        𝒯.s, 𝒯.ℓₘₐₓ, 𝒯.θ, 𝒯.Nϕ, 𝒯.centers, 𝒯.ring_ranges, 𝒯.plans, 𝒯.mode_m, 𝒯.Λ,
        𝒯.blocks, [similar(Fᵣ) for Fᵣ ∈ 𝒯.F], similar(𝒯.f̃), similar(𝒯.rhs)
    )
end

# The sample points become badly conditioned as ℓₘₐₓ grows — in Float64, a round trip loses
# about 3.5 digits by ℓₘₐₓ = 32 and 5.5 by ℓₘₐₓ = 48 at s = 0, and 9 by ℓₘₐₓ = 32 at s = 2 —
# and that conditioning belongs to the points themselves, so "Matrix" on the same points
# fares no better.  The constructor therefore measures it with `warn_if_inaccurate`.

pixels(𝒯::SSHTMinimal) = ring_pixels(𝒯.θ, 𝒯.Nϕ)
rotors(𝒯::SSHTMinimal) = from_spherical_coordinates.(pixels(𝒯))
npixels(𝒯::SSHTMinimal) = nmodes(𝒯)

# Synthesis on a copy, in `Complex{T}`, whatever the element type of the mode weights
function Base.:*(𝒯::SSHTMinimal{T}, f̃::SSHTData) where {T}
    mul!(𝒯, Array{Complex{T}}(synthesis_modes(𝒯, f̃)))
end
# Synthesis in place, in the storage of the mode weights — unless they had to be copied into
# the transform's range of ℓ, when it is the copy that receives the function values
function Base.:*(𝒯::SSHTMinimal{T, true}, f̃::SSHTData) where {T}
    d = synthesis_modes(𝒯, f̃)
    if d === array_view(f̃)
        check_complex_storage(d, "𝒯 * f̃")
        mul!(𝒯, d)
    else
        mul!(𝒯, Array{Complex{T}}(d))
    end
end
function LinearAlgebra.mul!(f, 𝒯::SSHTMinimal, f̃)
    check_modes(𝒯, f̃)
    check_pixels(𝒯, f)
    check_trailing(f, f̃)
    check_complex_output(f, "f")
    F = array_view(f)
    F .= array_view(f̃)
    mul!(𝒯, F)
    f
end

# Synthesis in place: the mode weights in `ff̃` are replaced by the function values.  Each
# ring's Fourier coefficients collect every mode, aliased or not, whose m is congruent to
# the coefficient's frequency modulo the size of the ring; an unnormalized inverse FFT then
# gives the values on the ring.
function LinearAlgebra.mul!(𝒯::SSHTMinimal{T}, ff̃::SSHTData) where {T}
    check_modes(𝒯, ff̃)
    check_complex_storage(ff̃, "mul!(𝒯, f̃)", "𝒯 * f̃")
    ff̃′ = reshape(array_view(ff̃), size(ff̃, 1), :)
    n = nmodes(𝒯)

    plans = 𝒯.plans
    for ₛf̃ ∈ eachcol(ff̃′)
        𝒯.f̃ .= ₛf̃  # every ring needs every mode, so the input is copied before it is overwritten
        @inbounds for r ∈ eachindex(𝒯.Nϕ)
            Fᵣ, N = 𝒯.F[r], 𝒯.Nϕ[r]
            fill!(Fᵣ, zero(Complex{T}))
            for i ∈ 1:n
                Fᵣ[mod(𝒯.mode_m[i], N) + 1] += 𝒯.f̃[i] * 𝒯.Λ[i, r]
            end
            plans.backward[plans.index[r]] * Fᵣ  # unnormalized inverse FFT: Σₘ Fₘ exp(imϕₖ)
            ₛf̃[𝒯.ring_ranges[r]] .= Fᵣ
        end
    end
    array_view(ff̃)
end

# Analysis on a copy, in `Complex{T}`, whatever the element type of the function values; it
# is a `ModeWeights` for one-dimensional data
function Base.:\(𝒯::SSHTMinimal{T}, f::SSHTData) where {T}
    check_pixels(𝒯, f)
    ldiv!(𝒯, Array{Complex{T}}(array_view(f)))
end
function Base.:\(𝒯::SSHTMinimal{T, true}, ff̃::SSHTData) where {T}
    check_pixels(𝒯, ff̃)
    check_complex_storage(ff̃, "𝒯 \\ f")
    ldiv!(𝒯, array_view(ff̃))
    in_place_modes(𝒯, ff̃)
end
function LinearAlgebra.ldiv!(f̃, 𝒯::SSHTMinimal, f)
    f̃ = analysis_output(𝒯, f̃, f)
    check_modes(𝒯, f̃)
    check_pixels(𝒯, f)
    check_trailing(f, f̃)
    check_complex_output(f̃, "f̃")
    array_view(f̃) .= array_view(f)
    ldiv!(𝒯, array_view(f̃))
    f̃
end

# Analysis in place: the function values in `ff̃` are replaced by the mode weights.  After
# an FFT on each ring, the groups of `minimal_blocks` are solved in order; the solution of
# each is then removed from every ring's coefficients, so that by the time a group is
# reached its coefficients hold only its own modes.
function LinearAlgebra.ldiv!(𝒯::SSHTMinimal{T}, ff̃::SSHTData) where {T}
    check_pixels(𝒯, ff̃)
    check_complex_storage(ff̃, "ldiv!(𝒯, f)", "𝒯 \\ f")
    ff̃′ = reshape(array_view(ff̃), size(ff̃, 1), :)

    plans = 𝒯.plans
    for ₛf ∈ eachcol(ff̃′)
        # Fourier coefficients of each ring, normalized as (1/N) Σₖ f(ϕₖ) exp(-imϕₖ)
        @inbounds for r ∈ eachindex(𝒯.Nϕ)
            Fᵣ, N = 𝒯.F[r], 𝒯.Nϕ[r]
            Fᵣ .= view(ₛf, 𝒯.ring_ranges[r]) ./ N
            plans.forward[plans.index[r]] * Fᵣ  # in-place FFT
        end

        # The coefficients are now all in 𝒯.F, so `ₛf` can receive the mode weights
        for block ∈ 𝒯.blocks
            nb = length(block.modes)
            rhs = view(𝒯.rhs, 1:nb)
            @inbounds for (e, (r, q)) ∈ enumerate(block.coefficients)
                rhs[e] = 𝒯.F[r][q]
            end
            ldiv!(block.lu, rhs)
            @inbounds for (u, i) ∈ enumerate(block.modes)
                ₛf[i] = rhs[u]
            end
            # Remove this group's modes from the coefficients of every ring they reach
            @inbounds for r ∈ eachindex(𝒯.Nϕ)
                Fᵣ, N = 𝒯.F[r], 𝒯.Nϕ[r]
                for (u, i) ∈ enumerate(block.modes)
                    Fᵣ[mod(𝒯.mode_m[i], N) + 1] -= rhs[u] * 𝒯.Λ[i, r]
                end
            end
        end
    end
    in_place_modes(𝒯, ff̃)
end
