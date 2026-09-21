"""
    SSHT{T}

Supertype of the spin-weighted spherical-harmonic transforms.  See [`SSHT`](@ref) (the
constructor), [`SSHTMatrix`](@ref), [`SSHTRS`](@ref) and [`SSHTMinimal`](@ref).
"""
abstract type SSHT{T<:Real} end


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

The first dimension of `f̃` must index the mode weights in the canonical ordering
`[f̃(ℓ, m) for ℓ ∈ abs(s):ℓₘₐₓ for m ∈ -ℓ:ℓ]` (see [`Yindex`](@ref); a [`ModeWeights`](@ref)
with `ℓₘᵢₙ = abs(s)` is accepted), and the first dimension of `f` must index the locations
at which the function is evaluated, in the order given by [`rotors`](@ref) or
[`pixels`](@ref).  Any following dimensions will be broadcast over.

The available `method`s are
- `"RS"` (default): [`SSHTRS`](@ref), the ring-based algorithm of Reinecke and Seljebotn,
  which scales as ``ℓₘₐₓ^3`` and is the choice for large ``ℓₘₐₓ``;
- `"Minimal"`: [`SSHTMinimal`](@ref), the optimal-dimensionality algorithm of Elahi et al.,
  which uses exactly as many samples as there are modes;
- `"Matrix"`: [`SSHTMatrix`](@ref), the direct dense-matrix method, which is the most
  accurate for small ``ℓₘₐₓ`` and lets the sample points be chosen freely.

The remaining keyword arguments are passed to the constructor of the chosen type.  Certain
types (`"Minimal"` and `"Matrix"`) also have an option to *always* act in place — meaning
that they simply re-use the input storage, even in an expression like `𝒯 \ f`; this is the
`inplace` keyword argument, and is part of the type of the resulting object.  Regardless of
that option, `LinearAlgebra.mul!` and `LinearAlgebra.ldiv!` force operation in place for
every type.

An `SSHT` object holds preallocated workspace, so it must not be used from several threads
at the same time; construct one object per thread.

# Half-integer indices

The spin weight and ``ℓₘₐₓ`` may be half-integers, spelled as `Rational`s with denominator 2
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
function SSHT(s::IndexSpelling, ℓₘₐₓ::IndexSpelling; method="RS", kwargs...)
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
# through this and re-dispatches to a worker whose signature is `where {IT<:HalfInteger}`.
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
# `ModeWeights` stopped being an `AbstractVector` in version 3, so the entry points that take
# "a map, or mode weights" have to name both.  Everything downstream goes through `strided`,
# which accepts either and hands back the raw storage.
const MapOrModes = Union{AbstractArray{<:Complex}, ModeWeights}

function check_modes(𝒯::SSHT, f̃::ModeWeights)
    if f̃.ℓₘᵢₙ != abs(𝒯.s) || f̃.ℓₘₐₓ != 𝒯.ℓₘₐₓ
        error(
            "The ModeWeights have ℓ ∈ $(f̃.ℓₘᵢₙ):$(f̃.ℓₘₐₓ), but the transform requires "
            * "ℓ ∈ $(abs(𝒯.s)):$(𝒯.ℓₘₐₓ)."
        )
    end
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

# Output containers
function mode_output(𝒯::SSHT{T}, f) where {T}
    if ndims(f) == 1
        ModeWeights(Vector{Complex{T}}(undef, nmodes(𝒯)), 𝒯.s, abs(𝒯.s), 𝒯.ℓₘₐₓ)
    else
        Array{Complex{T}}(undef, nmodes(𝒯), size(f)[2:end]...)
    end
end
pixel_output(𝒯::SSHT{T}, f̃) where {T} = Array{Complex{T}}(undef, npixels(𝒯), size(f̃)[2:end]...)


include("matrix.jl")
include("rs.jl")
include("minimal.jl")
