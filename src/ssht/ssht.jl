"""
    SSHT{T}

Supertype of the spin-weighted spherical-harmonic transforms.  See [`SSHT`](@ref) (the
constructor), [`SSHTMatrix`](@ref), [`SSHTRS`](@ref) and [`SSHTMinimal`](@ref).
"""
abstract type SSHT{T<:Real} end


# Sample data given to a transform must already be in the type `T` it works in.  Converting it
# silently is how a `QuatVec` becomes a rotation by π about its own direction, an unnormalized
# `Quaternion` scales every harmonic by a power of its norm, and `BigFloat` data is rounded to
# `Float64`; so, as for the calculators (see `check_rotor_type`), anything else is refused and
# the caller converts.  Integer colatitudes or weights convert exactly, and are accepted.
function check_sample_rotors(::Type{T}, Rθϕ) where {T}
    Rθϕ isa NonRotorData && error(not_a_rotor(Rθϕ))
    if !(Rθϕ isa AbstractVector{Rotor{T}})
        error(
            "This transform works in $T, so `Rθϕ` must be a vector of `Rotor{$T}`s, but it is a "
            * "$(typeof(Rθϕ)).  Pass `T` to work in another type, or convert the rotors."
        )
    end
    nothing
end
function check_sample_reals(::Type{T}, x, name) where {T}
    # (`float(Real) === Float64`, so the concreteness check is needed for a `Vector{Real}`.)
    if !(x isa AbstractVector{<:Real} && isconcretetype(eltype(x)) && float(eltype(x)) === T)
        error(
            "This transform works in $T, so `$name` must be a vector of $T, but it is a "
            * "$(typeof(x)).  Pass `T` to work in another type, or convert `$name`."
        )
    end
    nothing
end


@doc raw"""
    SSHT(s, ℓₘₐₓ; [method="RS"], [T=Float64], [kwargs...])

Construct an `SSHT` object to transform between spin-weighted spherical-harmonic mode
weights and function values — performing an ``s``-SHT.

This object behaves similarly to an `AbstractFFTs.Plan` object — specifically in the ability
to use the semantics of algebra to perform transforms.  For example, if the function values
are stored as a vector `f`, the mode weights as `f̃`, and the `SSHT` as `𝒯`, then we can
compute the function values from the mode weights (synthesis) as

    f = 𝒯 * f̃

or solve for the mode weights from the function values (analysis) as

    f̃ = 𝒯 \ f

The first dimension of `f̃` must index the mode weights in the canonical ordering `[f̃(ℓ, m)
for ℓ ∈ abs(s):ℓₘₐₓ for m ∈ -ℓ:ℓ]` (see [`Yindex`](@ref); a [`ModeWeights`](@ref) with `ℓₘᵢₙ
= abs(s)` is accepted), and the first dimension of `f` must index the locations at which the
function is evaluated, in the order given by [`rotors`](@ref) or [`pixels`](@ref).  Any
following dimensions will be broadcast over.

The available `method`s are
- `"RS"` (default): [`SSHTRS`](@ref), the ring-based algorithm of Reinecke and Seljebotn,
  which scales as ``ℓₘₐₓ^3`` and is the choice for large ``ℓₘₐₓ``;
- `"Matrix"`: [`SSHTMatrix`](@ref), the direct dense-matrix method, which is the most
  accurate for small ``ℓₘₐₓ`` and lets the sample points be chosen freely;
- `"Minimal"`: [`SSHTMinimal`](@ref), the optimal-dimensionality algorithm of Elahi et al.,
  which uses exactly as many samples as there are modes; mostly experimental, not very
  accurate.

The remaining keyword arguments are passed to the constructor of the chosen type.  Two of
the types — `"Minimal"` and `"Matrix"` — also have an option to *always* act in place —
meaning that they simply re-use the input storage, even in an expression like `𝒯 \ f`; this
is the `inplace` keyword argument, and is part of the type of the resulting object.  Even
then, the analysis of one-dimensional data returns a `ModeWeights`, as for every other
method, but one that wraps the input's own storage, now holding the mode weights; and
synthesis returns the storage itself, now holding the function values, as a plain array —
not the `ModeWeights` that may have held the input, whose labels no longer describe it.
Regardless of that option, `LinearAlgebra.mul!` and `LinearAlgebra.ldiv!` force operation in
place for every type.  The destination of `ldiv!(f̃, 𝒯, f)` for one-dimensional data may be a
bare vector at least `nmodes(𝒯)` long, and the mode weights then come back as a `ModeWeights`
over its first entries, rather than as the vector; `ldiv!(𝒯, x)` likewise returns a
`ModeWeights` wrapping `x`.

Sample data given as keyword arguments — the rotors `Rθϕ` of `"Matrix"`, and the
colatitudes `θ` and `quadrature_weights` of `"RS"` and `"Minimal"` — must already be in the
type `T` that the transform works in: rotors as `Rotor{T}`, and real numbers whose floating
point type is `T`.  Anything else is refused rather than converted.  In particular a general
`Quaternion` or a `QuatVec` is not taken for a rotor, and `BigFloat` data is not rounded to the
default `T=Float64`; pass `T=BigFloat` to work in that type.

An `SSHT` object holds preallocated workspace, so it must not be used from several threads
at the same time; construct one object per thread.

# Half-integer indices

The spin weight and ``ℓₘₐₓ`` may be half-integers, passed as `Rational`s with denominator 2
— `SSHT(1//2, 7//2)` — for the `"RS"` and `"Matrix"` methods; the `"Minimal"` method is
defined only for integer spin weights, and says so.  The mode weights are then indexed by
half-odd ``ℓ`` and ``m`` in the same canonical ordering, and the function values are, as
before, the values at the points that [`rotors`](@ref) returns.  Those points deserve more
attention than they need for an integer spin weight.  A function of half-integer spin weight
is a function on the double cover of the rotation group rather than on the sphere, so its
value "at a point of the sphere" depends on which of the two rotors above that point is
meant; the rotors used here are `from_spherical_coordinates(θ, ϕ)`, and the function values
are therefore antiperiodic in the azimuth — a full circuit in ``ϕ`` arrives at the antipodal
rotor, and changes the sign.  The sampling requirements are otherwise unchanged in form:
``N_ϕ ≥ 2ℓₘₐₓ+1`` and ``N_θ ≥ 2ℓₘₐₓ+1``, both of which are even numbers when ``ℓₘₐₓ`` is a
half-odd-integer.
"""
function SSHT(s::IndexArgument, ℓₘₐₓ::IndexArgument; method="RS", kwargs...)
    s, ℓₘₐₓ = transform_indices(s, ℓₘₐₓ)
    if method == "RS"
        return SSHTRS(s, ℓₘₐₓ; kwargs...)
    elseif method == "Minimal"
        return SSHTMinimal(s, ℓₘₐₓ; kwargs...)
    elseif method == "Matrix"
        return SSHTMatrix(s, ℓₘₐₓ; kwargs...)
    elseif method == "Direct"
        Base.depwarn(
            "The \"Direct\" s-SHT method has been renamed \"Matrix\"; use method=\"Matrix\".",
            :SSHT
        )
        return SSHTMatrix(s, ℓₘₐₓ; kwargs...)
    else
        error("""Unrecognized s-SHT method "$method"; use "RS", "Minimal", or "Matrix".""")
    end
end


# The transforms store their indices as `Int` or as `HalfOddInteger` — as `Int` rather than
# whatever `Integer` type the caller happened to use, which is what the fields' former `::Int`
# annotations did.  Every public constructor is a boundary that passes its two indices
# through this and re-dispatches to a worker whose signature is `where {IT<:IntegerHalf}`.
stored_index(x::Integer) = Int(x)
stored_index(x::HalfOddInteger) = x
transform_indices(s, ℓₘₐₓ) = map(stored_index, unify_indices(s, ℓₘₐₓ))

# The two places where the ring-based algorithm sees the kind of its indices.  For an integer
# spin weight each is exactly what the code did before half-integer support was added: the
# Fourier index of ``m`` is ``m``, and a ring's values are its inverse FFT.  For a half-odd
# spin weight ``e^{imϕ}`` is antiperiodic in ``ϕ``, so the FFT runs in the integer
# ``m̂ = m - 1/2 = ⌊m⌋`` and each sample is multiplied by ``e^{±iϕ/2}``, which for the ``k``-th
# point of a ring of ``N`` is ``e^{±iπk/N}``.
#
# There used to be a third: the table held ``{}_sY_{ℓ,m}(θ, 0)``, which for a half-odd spin
# weight is ``i^{2s}`` times a real function, and a `λreal` helper divided that constant phase
# out at every read.  The table is now built by an `sλlmCalculator`, which divides it out once
# at the source — so the phase is restored once per ring by `ring_values!` below, and the
# innermost loops touch real numbers directly.
@inline fourier_index(m::Integer) = m
@inline fourier_index(m::HalfOddInteger) = floor(Int, m)
# Synthesis: a ring's function values from its inverse FFT.
ring_values!(dest, Gy, ::Integer) = (dest .= Gy)
function ring_values!(dest, Gy, s::HalfOddInteger)
    T = real(eltype(Gy))
    N = length(Gy)
    phase = im_power(T, 2s)
    @inbounds for k ∈ 0:N-1
        dest[k+1] = (phase * cispi(T(k) / N)) * Gy[k+1]
    end
    dest
end
# Analysis: a ring's FFT input from its function values, with the quadrature factor.
ring_samples!(Gy, src, factor, ::Integer) = (Gy .= src .* factor)
function ring_samples!(Gy, src, factor, s::HalfOddInteger)
    T = real(eltype(Gy))
    N = length(Gy)
    phase = conj(im_power(T, 2s))
    @inbounds for k ∈ 0:N-1
        Gy[k+1] = (phase * cispi(-T(k) / N)) * (src[k+1] * factor)
    end
    Gy
end


"""
    pixels(𝒯)

Return the spherical coordinates `(θ, ϕ)` at which the transform `𝒯` evaluates functions,
in the order used by the first dimension of the function values.  See also [`rotors`](@ref).
"""
function pixels end


"""
    rotors(𝒯)

Return the `Rotor`s at which the transform `𝒯` evaluates functions, in the order used by the
first dimension of the function values.  See also [`pixels`](@ref).
"""
function rotors end

spin(𝒯::SSHT) = 𝒯.s
ℓₘₐₓ(𝒯::SSHT) = 𝒯.ℓₘₐₓ
ℓₘᵢₙ(𝒯::SSHT) = abs(𝒯.s)

"""
    nmodes(𝒯)

Number of mode weights the transform `𝒯` works with, `Ysize(abs(s), ℓₘₐₓ)`.
"""
nmodes(𝒯::SSHT) = Ysize(abs(𝒯.s), 𝒯.ℓₘₐₓ)

"""
    npixels(𝒯)

Number of sample points (function values) the transform `𝒯` works with.
"""
function npixels end

function Base.show(io::IO, 𝒯::SSHT)
    print(io, nameof(typeof(𝒯)), "{", eltype_real(𝒯), "}(s=$(𝒯.s), ℓₘₐₓ=$(𝒯.ℓₘₐₓ))")
end
Base.show(io::IO, ::MIME"text/plain", 𝒯::SSHT) = show(io, 𝒯)
eltype_real(::SSHT{T}) where {T} = T

# The mode-weight vector of an SSHT, as a plain vector of the right length
function check_modes(𝒯::SSHT, f̃)
    n = nmodes(𝒯)
    if size(f̃, 1) != n
        error(
            "The first dimension of the mode weights has length $(size(f̃, 1)), but "
            * "Ysize(abs(s)=$(abs(𝒯.s)), ℓₘₐₓ=$(𝒯.ℓₘₐₓ)) = $n is required."
        )
    end
end
# `ModeWeights` is not an `AbstractVector`, so the entry points that take "a map, or mode
# weights" have to name both.  Everything downstream goes through `array_view`,
# which accepts either and hands back the raw storage.
const MapOrModes = Union{AbstractArray{<:Complex}, ModeWeights}

# The labels are checked, not just the length: weights of spin -s (or 0) have the same length
# as those of spin s, and would otherwise be synthesized as spin s — or, as the output of
# `ldiv!`, be filled with spin-s weights while keeping the wrong label.  The length of the
# storage is then compared with the labels, because a `ModeWeights` wraps its vector without
# copying it, and the vector may have been resized since the labels were checked against it.
function check_modes(𝒯::SSHT, f̃::ModeWeights)
    if spin(f̃) != 𝒯.s
        error(
            "The ModeWeights have spin weight s=$(spin(f̃)), but the transform is for "
            * "s=$(𝒯.s)."
        )
    end
    if f̃.ℓₘᵢₙ != abs(𝒯.s) || f̃.ℓₘₐₓ != 𝒯.ℓₘₐₓ
        error(
            "The ModeWeights have ℓ ∈ $(f̃.ℓₘᵢₙ):$(f̃.ℓₘₐₓ), but the transform requires "
            * "ℓ ∈ $(abs(𝒯.s)):$(𝒯.ℓₘₐₓ)."
        )
    end
    check_storage_length(f̃)
end
function check_pixels(𝒯::SSHT, f)
    n = npixels(𝒯)
    if size(f, 1) != n
        error(
            "The first dimension of the function values has length $(size(f, 1)), but the "
            * "transform has $n sample points."
        )
    end
end
# Every three-argument `mul!` and `ldiv!` transforms each column of the trailing dimensions
# separately, so its input and output must agree in them; broadcasting would otherwise copy
# one input column into every column of the output.
function check_trailing(f, f̃)
    if size(f)[2:end] != size(f̃)[2:end]
        throw(DimensionMismatch(
            "Trailing dimensions of f $(size(f)[2:end]) and f̃ $(size(f̃)[2:end]) differ."
        ))
    end
end

# Output containers
function mode_output(𝒯::SSHT{T}, f) where {T}
    if ndims(f) == 1
        ModeWeights(Vector{Complex{T}}(undef, nmodes(𝒯)), 𝒯.s, abs(𝒯.s), 𝒯.ℓₘₐₓ)
    else
        Array{Complex{T}}(undef, nmodes(𝒯), size(f)[2:end]...)
    end
end
pixel_output(𝒯::SSHT{T}, f̃) where {T} = Array{Complex{T}}(undef, npixels(𝒯), size(f̃)[2:end]...)

# What the in-place `\` and `*` return, once they have overwritten their input's storage.
# Analysis of one-dimensional data returns its mode weights as a `ModeWeights`, as every
# other analysis does, but wrapping that same storage, so that nothing is copied and they are
# still indexed by (ℓ, m); returned as the bare vector, `(𝒯 \ f)[ℓ, m]` would silently read
# a linear index.  Synthesis returns the storage itself, as a plain array, and never a
# `ModeWeights` that held the input, whose labels would no longer describe its contents.
#
# A least-squares `SSHTMatrix` has more points than modes, and its two-argument `ldiv!`
# leaves the solution in the first `nmodes(𝒯)` entries of the longer input; only those are
# labelled, or returned, rather than the whole array with its residual components.
function in_place_modes(𝒯::SSHT, ff̃)
    d = array_view(ff̃)
    n = nmodes(𝒯)
    if ndims(d) == 1
        s, ℓₘᵢₙ, ℓₘₐₓ = 𝒯.s, abs(𝒯.s), 𝒯.ℓₘₐₓ
        length(d) == n ? ModeWeights(d, s, ℓₘᵢₙ, ℓₘₐₓ) : mode_weights_view(d, s, ℓₘᵢₙ, ℓₘₐₓ)
    else
        size(d, 1) == n ? d : view(d, 1:n, ntuple(_ -> Colon(), ndims(d) - 1)...)
    end
end
in_place_values(ff̃) = array_view(ff̃)

# The output of a three-argument analysis, `ldiv!(f̃, 𝒯, f)`.  A `ModeWeights` is used as
# given, and its labels are checked.  For one-dimensional data a bare vector may be longer
# than needed, and is wrapped as a `ModeWeights` over its first `nmodes(𝒯)` entries — which
# is then what comes back, labelled, rather than the vector.  Anything else must have
# exactly the right shape, and comes back as it is.
function analysis_output(𝒯::SSHT, f̃, f)
    if f̃ isa AbstractVector && ndims(array_view(f)) == 1
        mode_weights_view(f̃, 𝒯.s, abs(𝒯.s), 𝒯.ℓₘₐₓ)
    else
        f̃
    end
end


include("matrix.jl")
include("rs.jl")
include("minimal.jl")
