# The rings of the "RS" and "Minimal" transforms: their FFT plans, points, and phases.

# The FFTs along the rings.  One pair of in-place plans is made for each distinct number of
# points on a ring, rather than for each ring: every ring of the default "RS" grid has the
# same size, planning is almost all of the cost of constructing a transform, and FFTW plans
# one transform at a time, under a global lock.  A plan applies to any array of its size,
# stride and alignment, so each is executed on every buffer of that size; the buffers are
# arrays of their own, allocated with the alignment of the one the plan was made for.
#
# The plans of FFTW are made for a single thread.  A ring holds only O(ℓₘₐₓ) points, and a
# plan made with FFTW's default number of threads — that of Julia — spawns tasks and
# allocates on every execution of so short a transform, at a cost that exceeds the transform
# itself many times over.  Other element types are planned by the generic FFT, which takes
# no options.
#
# The fields that describe the plans come first, and are what the plans are remade from when
# a transform is deserialized (see `Serialization.deserialize` below).
#
# `T` is the real type the transform works in, and `P` and `BP` are the types of the forward
# and backward plans.
struct RingPlans{T<:Real, P, BP}
    sizes::Vector{Int}  # the distinct numbers of points on a ring
    index::Vector{Int}  # for each ring, the index of its number of points in `sizes`
    flags::UInt32  # the FFTW planner's flags and time limit, as given to the constructor
    timelimit::Float64
    forward::Vector{P}  # for each size, the in-place FFT, Σₖ xₖ e^{-imϕₖ}
    backward::Vector{BP}  # and the in-place unnormalized inverse FFT, Σₘ xₘ e^{+imϕₖ}
end

function ring_plans(::Type{T}, Nϕ::AbstractVector{Int}; flags, timelimit) where {T}
    sizes = unique(Nϕ)
    position = Dict(N => j for (j, N) ∈ enumerate(sizes))
    ring_plans(T, sizes, [position[N] for N ∈ Nϕ], flags, timelimit)
end
function ring_plans(::Type{T}, sizes, index, flags, timelimit) where {T}
    planned = [ring_fft_plans(Vector{Complex{T}}(undef, N), flags, timelimit) for N ∈ sizes]
    forward, backward = first.(planned), last.(planned)
    RingPlans{T, eltype(forward), eltype(backward)}(
        sizes, index, flags, Float64(timelimit), forward, backward
    )
end
function ring_fft_plans(
    buffer::Vector{Complex{T}}, flags, timelimit
) where {T<:Union{Float32, Float64}}
    (
        plan_fft!(buffer; flags, timelimit, num_threads=1),
        plan_bfft!(buffer; flags, timelimit, num_threads=1),
    )
end
ring_fft_plans(buffer, flags, timelimit) = (plan_fft!(buffer), plan_bfft!(buffer))

# A deep copy shares the plans rather than copying them.  A copy of an FFTW plan object
# would wrap the same pointer to the plan without owning it, and would execute freed memory
# once the original had been garbage-collected.  Sharing is safe, because no plan is
# modified after it is made, and FFTW allows one plan to be executed on different arrays at
# the same time.
Base.deepcopy_internal(p::RingPlans, ::IdDict) = p

# Serialization writes every pointer as a null pointer.  An FFTW plan object wraps a pointer
# to a plan in the memory of the process that made it, so the plans of a transform would
# arrive in another process — or be read back from a file — as plans that crash the process
# when they are executed.  The plans of a `RingPlans` are therefore not reconstructed from
# what was written: they are made again in the receiving process, of the sizes and with the
# planner options that were written.  What is written is the default serialization of the
# struct, its fields in order, and the plans it contains are read and discarded.  (A plan of
# the generic FFT holds no pointer, and would survive, but is remade in the same way.)
function Serialization.deserialize(
    s::Serialization.AbstractSerializer, ::Type{RingPlans{T, P, BP}}
) where {T, P, BP}
    sizes = Serialization.deserialize(s)
    index = Serialization.deserialize(s)
    flags = Serialization.deserialize(s)
    timelimit = Serialization.deserialize(s)
    Serialization.deserialize(s)  # the forward plans of the sending process
    Serialization.deserialize(s)  # and the backward plans
    ring_plans(T, sizes, index, flags, timelimit)::RingPlans{T, P, BP}
end

# The sample points of a transform on rings, ring by ring, with the azimuth ``ϕ_k = 2πk/N``
# for ``k = 0, …, N-1`` on a ring of ``N`` points.
function ring_pixels(θ::Vector{T}, Nϕ) where {T}
    let π = T(π)
        [@SVector [θᵣ, iϕ * 2π / N] for (θᵣ, N) ∈ zip(θ, Nϕ) for iϕ ∈ 0:N-1]
    end
end

# The ring-based algorithm sees the kind of its indices in two places: the Fourier index of
# ``m`` on a ring, which is `floor_int(m)`, and the phases below.  For an integer spin
# weight the Fourier index of ``m`` is ``m``, and a ring's values are its inverse FFT.  For
# a half-odd spin weight ``e^{imϕ}`` is antiperiodic in ``ϕ``, so the FFT runs in the
# integer ``m̂ = m - 1/2 = ⌊m⌋`` and each sample is multiplied by ``e^{±iϕ/2}``, which for
# the ``k``-th point of a ring of ``N`` is ``e^{±iπk/N}``.  The harmonics on the rings are
# ``i^{2s}`` times real functions of ``θ``, which is what the `sλlmCalculator` tabulates, so
# that constant phase is restored once per ring here, and the innermost loops touch real
# numbers only.  The factors depend only on the size of the ring, and are tabulated once
# for each size by `ring_phases`: for synthesis ``i^{2s} e^{iπk/N}``, and for analysis
# ``i^{-2s} e^{-iπk/N}``.  For an integer spin weight the tables are empty.
function ring_phases(::Type{T}, sizes, ::Integer) where {T}
    ([Complex{T}[] for _ ∈ sizes], [Complex{T}[] for _ ∈ sizes])
end
function ring_phases(::Type{T}, sizes, s::HalfOddInteger) where {T}
    phase = im_power(T, 2s)
    (
        [[phase * cispi(T(k) / N) for k ∈ 0:N-1] for N ∈ sizes],
        [[conj(phase) * cispi(-T(k) / N) for k ∈ 0:N-1] for N ∈ sizes],
    )
end
# Synthesis: a ring's function values from its inverse FFT.
ring_values!(dest, Gy, phases, ::Integer) = (dest .= Gy)
function ring_values!(dest, Gy, phases, ::HalfOddInteger)
    @inbounds for k ∈ eachindex(Gy, phases)
        dest[k] = phases[k] * Gy[k]
    end
    dest
end
# Analysis: a ring's FFT input from its function values, with the quadrature factor.
ring_samples!(Gy, src, factor, phases, ::Integer) = (Gy .= src .* factor)
function ring_samples!(Gy, src, factor, phases, ::HalfOddInteger)
    @inbounds for k ∈ eachindex(Gy, phases)
        Gy[k] = phases[k] * (src[k] * factor)
    end
    Gy
end
