# One group of modes solved together: the modes of a set of m values that alias into one
# another (a strongly connected component of the graph described in `minimal_blocks`), and
# the Fourier coefficients — `(ring, index)` pairs, with index `mod(m, Nϕ)+1` — that
# determine them.  There are exactly as many coefficients as modes.
struct MinimalBlock{T}
    modes::Vector{Int}
    coefficients::Vector{Tuple{Int, Int}}
    lu::LinearAlgebra.LU{T, Matrix{T}, Vector{Int}}
end

"""
    SSHTMinimal(s, ℓₘₐₓ; T=Float64, θ=minimal_rings(s, ℓₘₐₓ, T).θ, plan_fft_flags=FFTW.ESTIMATE, plan_fft_timelimit=Inf, inplace=true)

Construct an ``s``-SHT object that uses the optimal-dimensionality algorithm described by
[Elahi et al.](@cite Elahi_2018), which samples the function at exactly as many points as
there are modes.  This may also be achieved by calling the main [`SSHT`](@ref) function with
the same keywords, along with `method="Minimal"`.

!!! warning
    This method is experimental and not very accurate.  The round-trip error grows
    exponentially with ``ℓₘₐₓ``, and the constructor warns when fewer than half the digits
    of `T` would survive.  The `"RS"` method (the default) has no such limitation, but does
    not have optimal dimensionality.  The `"Matrix"` method does have optimal
    dimensionality, but its memory consumption scales poorly.

The function is sampled on ``ℓₘₐₓ-|s|+1`` "rings" at constant colatitude, each holding an odd
number of equally spaced points starting at ``ϕ = 0``.  Their sizes, and the default
colatitudes, are given by [`minimal_rings`](@ref); the `θ` keyword argument may give other
colatitudes, one for each ring in the order listed there, which is also the order of the
sample points.  See [`pixels`](@ref) and [`rotors`](@ref) for the sample points themselves.

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

Because this algorithm achieves optimal dimensionality, the transformation is performed in
place by default: `𝒯 * f̃` overwrites the storage of `f̃` with the function values (and
returns that storage), and `𝒯 \\ f` overwrites `f` with the mode weights (and returns them,
for one-dimensional `f`, as a `ModeWeights` wrapping that storage).  If this is not desired,
pass the keyword argument `inplace=false`, which makes those operations work on a copy of
the input.  See [`SSHT`](@ref).

The values ``{}_sλ_{ℓ,m}(θ_r)`` of every mode on every ring are precomputed at construction
(with one batched [`sλlmCalculator`](@ref)) and stored, which takes ``O(ℓₘₐₓ^3)`` memory, as
are the LU decompositions of the matrices for the groups of ``m`` values.  The object holds
workspace for the transforms, so it must not be used from several threads at once.

Even so, the sample points become increasingly badly conditioned as ℓₘₐₓ grows, for every
spin weight: in `Float64` a round trip loses about 5 digits by ℓₘₐₓ = 32 and 10 by ℓₘₐₓ = 48
at ``s = 0``, and more at larger ``|s|`` — about 10 by ℓₘₐₓ = 32 at ``s = 2``.  The
constructor therefore measures the error of one round trip, and warns when fewer than half
the digits of `T` would survive; the `"RS"` method has no such limitation.

This method is defined only for integer spin weights.  Its bookkeeping — rings of an odd
number of points, and the aliasing of ``m`` into rings too small to hold it — is written for
integer indices throughout, and has not been extended; a half-integer spin weight is refused
with a message naming the two methods, `"RS"` and `"Matrix"`, that do accept one.
"""
struct SSHTMinimal{T<:Real, Inplace, P, BP} <: SSHT{T}
    s::Int
    ℓₘₐₓ::Int
    θ::Vector{T}  # colatitude of each ring
    Nϕ::Vector{Int}  # number of points on each ring (odd)
    centers::Vector{Int}  # center of each ring's window of frequencies (see `minimal_rings`)
    θindices::Vector{UnitRange{Int}}  # pixel indices of each ring
    plans::Vector{P}  # forward in-place FFT plans, one per ring
    bplans::Vector{BP}  # backward in-place FFT plans, one per ring
    mode_m::Vector{Int}  # m of each mode, in the canonical order
    Λ::Matrix{T}  # Λ[i, r] = ₛλ_{ℓ,m}(θ_r) for the mode (ℓ, m) with index i
    blocks::Vector{MinimalBlock{T}}  # the groups of modes solved together, in solution order
    F::Vector{Vector{Complex{T}}}  # workspace: Fourier coefficients of each ring
    f̃::Vector{Complex{T}}  # workspace: a copy of the mode weights, for synthesis
    rhs::Vector{Complex{T}}  # workspace: the right-hand side of one block's system
end



@doc raw"""
    minimal_rings(s, ℓₘₐₓ, [T=Float64])

The rings on which [`SSHTMinimal`](@ref) samples a function of spin weight `s` band-limited
at `ℓₘₐₓ`, as a named tuple `(; Nϕ, centers, θ)`: the number of points on each ring, the center
of each ring's window of frequencies, and the default colatitude of each ring.  The rings are
listed in order of size (for rings of equal size, the northern first), and there are
``ℓₘₐₓ-|s|+1`` of them, holding ``(ℓₘₐₓ+1)^2 - s^2`` points in all — exactly the number of
modes.

A ring of ``N = 2k+1`` points cannot distinguish frequencies ``m`` that differ by a multiple of
``N``; the analysis treats its Fourier coefficients as measuring the ``N`` consecutive
frequencies ``|m - c| ≤ k`` of its window, centered on ``c``, and removes the aliases of all
other frequencies from them.  For each ``m`` there must be as many rings whose windows
include ``m`` as there are modes with that ``m``, namely ``ℓₘₐₓ - \max(|m|, |s|) + 1``.  The
windows ``|m| ≤ a`` for ``a ∈ |s|:ℓₘₐₓ``, one ring for each, satisfy this, and are what is
used for ``s = 0``, with the colatitudes of [`sorted_rings`](@ref).

For ``s ≠ 0`` that arrangement is badly conditioned (see [`SSHTMinimal`](@ref)).  Instead, pairs
of those windows are recentered.  The two windows ``|m| ≤ a`` and ``|m| ≤ a+2d`` cover every
``m`` exactly as often as the two windows ``|m + d| ≤ a+d`` and ``|m - d| ≤ a+d`` do, and the
latter become a ring of ``2(a+d)+1`` points in the northern hemisphere, whose window is
centered on ``-d\,\mathrm{sign}(s)``, and one of the same size in the southern hemisphere,
centered on ``+d\,\mathrm{sign}(s)`` — toward the frequencies ``∓s`` that dominate near each
pole.  The pairs are chosen greedily, first with ``d = |s|`` and ``a`` in increasing order,
whenever both windows are still available, and then with successively smaller ``d`` down to
1, which matters when ``ℓₘₐₓ < 3|s|`` and no window has a partner ``2|s|`` larger.  The windows
left unpaired — the largest ones — remain centered on 0, and become the rings nearest the
equator.  The default colatitudes are equally spaced, ``θ = iπ/(n+1)`` for ``i ∈
1:n`` with ``n`` the number of rings; the northern rings take the slots nearest the north
pole, smallest first, and likewise in the south, and the rings centered on 0 take the
remaining slots in the order [`sorted_rings`](@ref) uses.  For ``s = 0`` this reproduces
[`sorted_rings`](@ref) exactly.
"""
function minimal_rings(s::Integer, ℓₘₐₓ::Integer, ::Type{T}=Float64) where {T}
    if abs(s) > ℓₘₐₓ
        error("|s|=$(abs(s)) exceeds ℓₘₐₓ=$ℓₘₐₓ; there are no such modes.")
    end
    s, ℓₘₐₓ = Int(s), Int(ℓₘₐₓ)
    # Pair the windows |m| ≤ a and |m| ≤ a+2d, as described above, with the largest shift d
    # available first
    rings = Tuple{Int, Int}[]  # (k, center), for a ring of 2k+1 points
    used = falses(ℓₘₐₓ + 1)
    for d ∈ abs(s):-1:1, a ∈ abs(s):ℓₘₐₓ
        b = a + 2d
        if !used[a+1] && b ≤ ℓₘₐₓ && !used[b+1]
            used[a+1] = used[b+1] = true
            push!(rings, (a + d, -sign(s) * d), (a + d, sign(s) * d))  # north, then south
        end
    end
    for a ∈ abs(s):ℓₘₐₓ
        used[a+1] || push!(rings, (a, 0))
    end
    north(center) = center * s < 0  # centered on the side of -s
    south(center) = center * s > 0
    sort!(rings, by=((k, c),) -> (k, north(c) ? 0 : south(c) ? 2 : 1))

    # Equally spaced slots: northern rings from the north pole inward, southern rings from
    # the south pole inward, and the rest in the middle, arranged as `sorted_rings` arranges
    # its rings (so that for s = 0 the result is identical to it)
    n = length(rings)
    slots = collect(LinRange{T}(0, π, n + 2))[begin+1:end-1]
    northern = [i for (i, (k, c)) ∈ enumerate(rings) if north(c)]
    southern = [i for (i, (k, c)) ∈ enumerate(rings) if south(c)]
    middle = [i for (i, (k, c)) ∈ enumerate(rings) if c == 0]
    θ = Vector{T}(undef, n)
    for (q, i) ∈ enumerate(northern)  # already in order of size
        θ[i] = slots[q]
    end
    for (q, i) ∈ enumerate(southern)
        θ[i] = slots[n + 1 - q]
    end
    let πo2 = prevfloat(T(π)/2, s), np = length(northern)
        middle_slots = sort(
            slots[np+1:n-np], lt=(x,y)->(abs(x-πo2)<abs(y-πo2)), rev=true
        )
        for (q, i) ∈ enumerate(middle)  # in order of size, smallest farthest from π/2
            θ[i] = middle_slots[q]
        end
    end
    (; Nϕ=[2k+1 for (k, c) ∈ rings], centers=[c for (k, c) ∈ rings], θ)
end

# The groups of m values that must be solved together, in an order in which each group's
# aliases into the others' rings are known by the time they are needed.  There is an edge
# m′ → m whenever some ring's window includes m but not m′ ≡ m (mod Nϕ), since m′ then
# aliases into the coefficient that measures m, and must be removed from it first.  The
# groups are the strongly connected components of that graph, which Kosaraju's algorithm
# produces in topological order (sources first).  For windows all centered on 0 — in
# particular for s = 0 — every group is a single m, and the order is that of decreasing |m|.
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

# The public constructor is the boundary; the indices arrive at the worker as `Int`s, or as
# `HalfOddInteger`s, for which the worker is the refusal described in the docstring.
function SSHTMinimal(s::IndexArgument, ℓₘₐₓ::IndexArgument; T::Type{TT}=Float64, kwargs...) where {TT}
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
    θ=nothing,
    plan_fft_flags=FFTW.ESTIMATE, plan_fft_timelimit=Inf,
    inplace=true
) where {TT}
    if abs(s) > ℓₘₐₓ
        error("|s|=$(abs(s)) exceeds ℓₘₐₓ=$ℓₘₐₓ; there are no such modes.")
    end
    rings = minimal_rings(s, ℓₘₐₓ, TT)
    nrings = length(rings.Nϕ)
    if θ === nothing
        θ = rings.θ
    elseif length(θ) != nrings
        error("Length of θ ($(length(θ))) must equal ℓₘₐₓ-abs(s)+1 ($nrings).")
    end
    check_sample_reals(TT, θ, "θ")
    θ = Vector{TT}(θ)
    Nϕ, centers = rings.Nϕ, rings.centers

    θindices = let stops = cumsum(Nϕ)
        [(stop - N + 1):stop for (N, stop) ∈ zip(Nϕ, stops)]
    end
    F = [Vector{Complex{TT}}(undef, N) for N ∈ Nϕ]
    plans, bplans = if TT ∈ (Float64, Float32)  # Only supported types in FFTW
        (
            [plan_fft!(Fᵣ; flags=plan_fft_flags, timelimit=plan_fft_timelimit) for Fᵣ ∈ F],
            [plan_bfft!(Fᵣ; flags=plan_fft_flags, timelimit=plan_fft_timelimit) for Fᵣ ∈ F],
        )
    else
        ([plan_fft!(Fᵣ) for Fᵣ ∈ F], [plan_bfft!(Fᵣ) for Fᵣ ∈ F])
    end

    # Tables of ₛλₗₘ(θᵣ) for every mode on every ring, computed with one batched calculator
    # in angle mode
    n = Ysize(abs(s), ℓₘₐₓ)
    mode_m = Vector{Int}(undef, n)
    Λ = Matrix{TT}(undef, n, nrings)
    λ = sλlmCalculator(θ, ℓₘₐₓ, s)  # θ is already a Vector{TT}, which fixes the type
    iₛ = spin_index(λ, s)
    Λℓ = λ.Yˡ  # [ring, spin, m+ℓ+1]
    for ℓ ∈ abs(s):ℓₘₐₓ
        recurrence!(λ, ℓ)
        for m ∈ -ℓ:ℓ
            i = Yindex(ℓ, m, abs(s))
            mode_m[i] = m
            for r ∈ 1:nrings
                Λ[i, r] = Λℓ[r, iₛ, m + ℓ + 1]
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
        MinimalBlock{TT}(modes, coefficients, LinearAlgebra.lu(M))
    end
    rhs = Vector{Complex{TT}}(undef, maximum(b -> length(b.modes), blocks))

    𝒯 = SSHTMinimal{TT, inplace, eltype(plans), eltype(bplans)}(
        s, ℓₘₐₓ, θ, Nϕ, centers, θindices, plans, bplans, mode_m, Λ, blocks, F,
        Vector{Complex{TT}}(undef, n), rhs
    )
    warn_if_inaccurate(𝒯)
    𝒯
end

# The sample points become badly conditioned as ℓₘₐₓ grows — in Float64, a round trip loses
# about 5 digits by ℓₘₐₓ = 32 and 10 by ℓₘₐₓ = 48 at s = 0, and 10 by ℓₘₐₓ = 32 at s = 2 —
# and that conditioning belongs to the points themselves, so "Matrix" on the same points
# fares no better.  Nothing in an individual transform reveals it, so the constructor
# measures it once, by synthesizing and analyzing a fixed set of unit weights with
# quasi-random phases — about the cost of one transform — and warns when fewer than half the
# digits of `T` survive.
function warn_if_inaccurate(𝒯::SSHTMinimal{T}) where {T}
    φ = (√5 - 1) / 2
    f̃ = [cis(T(2π) * T(mod(i * φ, 1))) for i ∈ 1:nmodes(𝒯)]
    f = copy(f̃)
    ldiv!(𝒯, mul!(𝒯, f))
    maxerror = maximum(abs, f - f̃)
    if !(maxerror ≤ √eps(T))  # also catches NaN
        @warn (
            "The \"Minimal\" s-SHT with s=$(𝒯.s), ℓₘₐₓ=$(𝒯.ℓₘₐₓ) and T=$T is inaccurate: a "
            * "round trip of unit mode weights has a maximum error of "
            * "$(round(Float64(maxerror), sigdigits=2)).  Its sample points are badly "
            * "conditioned at this ℓₘₐₓ; the \"RS\" method (the default) is accurate here."
        )
    end
    nothing
end

function pixels(𝒯::SSHTMinimal{T}) where {T}
    let π = T(π)
        [
            @SVector [θ, iϕ * 2π / N]
            for (θ, N) ∈ zip(𝒯.θ, 𝒯.Nϕ)
            for iϕ ∈ 0:N-1
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
    in_place_values(f̃)
end
function LinearAlgebra.mul!(f, 𝒯::SSHTMinimal, f̃)
    check_modes(𝒯, f̃)
    check_pixels(𝒯, f)
    f .= array_view(f̃)
    mul!(𝒯, f)
end

# Synthesis in place: the mode weights in `ff̃` are replaced by the function values.  Each
# ring's Fourier coefficients collect every mode, aliased or not, whose m is congruent to the
# coefficient's frequency modulo the size of the ring; an unnormalized inverse FFT then gives
# the values on the ring.
function LinearAlgebra.mul!(𝒯::SSHTMinimal{T}, ff̃) where {T}
    check_modes(𝒯, ff̃)
    ff̃′ = reshape(array_view(ff̃), size(ff̃, 1), :)
    n = nmodes(𝒯)

    for ₛf̃ ∈ eachcol(ff̃′)
        𝒯.f̃ .= ₛf̃  # every ring needs every mode, so the input is copied before it is overwritten
        @threads for r ∈ eachindex(𝒯.Nϕ)
            Fᵣ, N = 𝒯.F[r], 𝒯.Nϕ[r]
            fill!(Fᵣ, zero(Complex{T}))
            @inbounds for i ∈ 1:n
                Fᵣ[mod(𝒯.mode_m[i], N) + 1] += 𝒯.f̃[i] * 𝒯.Λ[i, r]
            end
            𝒯.bplans[r] * Fᵣ  # In-place unnormalized inverse FFT: Σₘ Fₘ exp(imϕₖ)
            @inbounds ₛf̃[𝒯.θindices[r]] .= Fᵣ
        end
    end
    ff̃
end

function Base.:\(𝒯::SSHTMinimal, f)
    check_pixels(𝒯, f)
    ldiv!(𝒯, copy(array_view(f)))  # a `ModeWeights` for one-dimensional data
end
function Base.:\(𝒯::SSHTMinimal{T, true}, ff̃) where {T}
    check_pixels(𝒯, ff̃)
    ldiv!(𝒯, array_view(ff̃))
    in_place_modes(𝒯, ff̃)
end
function LinearAlgebra.ldiv!(f̃, 𝒯::SSHTMinimal, f)
    f̃ = analysis_output(𝒯, f̃, f)
    check_modes(𝒯, f̃)
    check_pixels(𝒯, f)
    array_view(f̃) .= f
    ldiv!(𝒯, array_view(f̃))
    f̃
end

# Analysis in place: the function values in `ff̃` are replaced by the mode weights.  After an
# FFT on each ring, the groups of `minimal_blocks` are solved in order; the solution of each is
# then removed from every ring's coefficients, so that by the time a group is reached its
# coefficients hold only its own modes.
function LinearAlgebra.ldiv!(𝒯::SSHTMinimal{T}, ff̃) where {T}
    check_pixels(𝒯, ff̃)
    ff̃′ = reshape(array_view(ff̃), size(ff̃, 1), :)

    for ₛf ∈ eachcol(ff̃′)
        # Fourier coefficients of each ring, normalized as (1/N) Σₖ f(ϕₖ) exp(-imϕₖ)
        @threads for r ∈ eachindex(𝒯.Nϕ)
            Fᵣ, N = 𝒯.F[r], 𝒯.Nϕ[r]
            @inbounds Fᵣ .= view(ₛf, 𝒯.θindices[r]) ./ N
            𝒯.plans[r] * Fᵣ  # In-place FFT
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
            @threads for r ∈ eachindex(𝒯.Nϕ)
                Fᵣ, N = 𝒯.F[r], 𝒯.Nϕ[r]
                @inbounds for (u, i) ∈ enumerate(block.modes)
                    Fᵣ[mod(𝒯.mode_m[i], N) + 1] -= rhs[u] * 𝒯.Λ[i, r]
                end
            end
        end
    end
    in_place_modes(𝒯, ff̃)
end
