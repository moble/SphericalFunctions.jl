"""
    SSHT{T}

Supertype of the spin-weighted spherical-harmonic transforms.  See [`SSHT`](@ref) (the
constructor), [`SSHTMatrix`](@ref), [`SSHTRS`](@ref) and [`SSHTMinimal`](@ref).  `T` is the
real type the transform works in.
"""
abstract type SSHT{T<:Real} end


# Sample data given to a transform must already be in the type `T` it works in.  Converting
# it silently is how a `QuatVec` becomes a rotation by π about its own direction, and
# `BigFloat` data is rounded to `Float64`; so, as for the calculators (see
# `check_rotor_type`), anything else is refused and the caller converts.  A `Quaternion`
# denotes the rotation of its normalization, as it does throughout, and the transform keeps
# that rotation.  Integer colatitudes or weights convert exactly, and are accepted for every
# `T`.
function check_sample_rotors(::Type{T}, Rθϕ) where {T}
    Rθϕ isa NonRotorData && throw(ArgumentError(not_a_rotor(Rθϕ)))
    if !(Rθϕ isa Union{AbstractVector{Rotor{T}}, AbstractVector{Quaternionic.Quaternion{T}}})
        throw(ArgumentError(
            "This transform works in $T, so `Rθϕ` must be a vector of `Rotor{$T}`s or "
            * "`Quaternion{$T}`s, but it is a $(typeof(Rθϕ)).  Pass `T` to work in another "
            * "type, or convert the rotors."
        ))
    end
    nothing
end
function check_sample_reals(::Type{T}, x, name) where {T}
    # (The concreteness check refuses a `Vector{Real}`, whose elements could be of any type.)
    if !(
        x isa AbstractVector{<:Real} && isconcretetype(eltype(x))
        && (eltype(x) === T || eltype(x) <: Integer)
    )
        throw(ArgumentError(
            "This transform works in $T, so `$name` must be a vector of $T, but it is a "
            * "$(typeof(x)).  Pass `T` to work in another type, or convert `$name`."
        ))
    end
    nothing
end

# A transform computes in a floating-point type; any other would fail deep inside, at
# `T(π)`.
function check_transform_type(::Type{T}) where {T}
    if !(T <: AbstractFloat)
        throw(ArgumentError(
            "A transform works in a floating-point type, such as `Float64` or `BigFloat`, but "
            * "T=$T is not one."
        ))
    end
    nothing
end

# The option that makes "Matrix" and "Minimal" act in place is a type parameter, so it must
# be a `Bool`; anything else would make a type that no method treats as in place.
function check_inplace_option(inplace)
    if !(inplace isa Bool)
        throw(ArgumentError("inplace=$(repr(inplace)) must be `true` or `false`."))
    end
    nothing
end


@doc raw"""
    SSHT(s, ℓₘₐₓ, [T=Float64]; [method="RS"], [kwargs...])

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
for ℓ ∈ abs(s):ℓₘₐₓ for m ∈ -ℓ:ℓ]` (see [`Yindex`](@ref)), and the first dimension of `f`
must index the locations at which the function is evaluated, in the order given by
[`rotors`](@ref) or [`pixels`](@ref).  Any following dimensions will be broadcast over.  The
transform works in the floating-point type `T`, and its results are complex numbers of that
type; the data given to it may be real or complex.

The available `method`s, each given as a string or a `Symbol`, are
- `"RS"` (default): [`SSHTRS`](@ref), the ring-based algorithm of Reinecke and Seljebotn,
  which scales as ``ℓₘₐₓ^3`` and is the choice for large ``ℓₘₐₓ``;
- `"Matrix"`: [`SSHTMatrix`](@ref), the direct dense-matrix method, which is the fastest for
  small ``ℓₘₐₓ`` and lets the sample points be chosen freely, with round-trip errors of
  about ``10^{-14}`` for ``ℓₘₐₓ ≲ 24``, somewhat larger than those of `"RS"`;
- `"Minimal"`: [`SSHTMinimal`](@ref), the optimal-dimensionality algorithm of Elahi et al.,
  which uses exactly as many samples as there are modes; mostly experimental, not very
  accurate.

The remaining keyword arguments are passed to the constructor of the chosen type.

# Mode weights

A [`ModeWeights`](@ref) is accepted wherever mode weights are, and its labels are checked
against the transform: weights of another spin weight are refused, even when they have the
right number of entries.  Synthesis, `𝒯 * w`, accepts weights of the transform's spin
weight with any range of ``ℓ`` up to the transform's ``ℓₘₐₓ``: the modes that `w` lacks are
zero, and its entries with ``ℓ < |s|``, which belong to no function of spin weight ``s``,
are ignored.  The three-argument `mul!(f, 𝒯, w)` and the outputs of analysis need exactly
the range ``|s| ≤ ℓ ≤ ℓₘₐₓ``; `ModeWeights(w; ℓₘᵢₙ, ℓₘₐₓ)` copies mode weights into another
range.  Analysis of one-dimensional data returns a `ModeWeights`, and of more dimensions a
plain array.

# Acting in place

The `"RS"` method never acts in place: `𝒯 * f̃` and `𝒯 \ f` allocate their results and
leave their arguments alone.  The other two methods have an `inplace` keyword argument,
which is part of the type of the resulting object, and which makes them reuse the storage of
their argument.  It is on by default for both — for `"Matrix"` whenever there are exactly as
many sample points as modes, as there are with its default points — so that whether `𝒯 \ f`
preserves `f` depends on the method; pass `inplace=false` for the behavior of `"RS"`.  With
`inplace=true`,
- `"Matrix"` acts in place for analysis only: `𝒯 \ f` overwrites `f` with the mode weights,
  while `𝒯 * f̃` always allocates its result, since a matrix product cannot write over its
  own input;
- `"Minimal"` acts in place in both directions: `𝒯 * f̃` overwrites `f̃` with the function
  values, and `𝒯 \ f` overwrites `f` with the mode weights.  Analysis in place of
  one-dimensional data returns a `ModeWeights` wrapping the argument's own storage, now
  holding the mode weights, and synthesis in place returns the storage itself, now holding
  the function values, as a plain array — not the `ModeWeights` that may have held the
  input, whose labels no longer describe it.  Mode weights whose range of ``ℓ`` differs from
  the transform's are first copied into its range, and it is the copy that synthesis
  overwrites and returns.  Acting in place needs storage that can hold the complex results,
  and real storage is refused with an `ArgumentError`; real data are transformed with
  `inplace=false`.

Whatever the option, the two-argument `LinearAlgebra.ldiv!(𝒯, x)` analyzes in place for
`"Matrix"` and `"Minimal"`, returning the mode weights as `𝒯 \ x` would, and the
two-argument `LinearAlgebra.mul!(𝒯, x)` synthesizes in place for `"Minimal"`, returning the
function values as a plain array.  `"RS"` has neither, because its number of sample points
differs from the number of modes, and `"Matrix"` has no two-argument `mul!`.  The
three-argument `mul!(f, 𝒯, f̃)` and `ldiv!(f̃, 𝒯, f)` write into the given output for
every type, and never touch their input.  The output of `ldiv!(f̃, 𝒯, f)` for
one-dimensional data may be a bare vector at least `nmodes(𝒯)` long, and the mode weights
then come back as a `ModeWeights` over its first entries, rather than as the vector.

# Sample data

Sample data given as keyword arguments — the rotors `Rθϕ` of `"Matrix"`, and the colatitudes
`θ` and `quadrature_weights` of `"RS"` and `"Minimal"` — must already be in the type `T`
that the transform works in: rotors as `Rotor{T}`, and real numbers of type `T`, or
integers, which convert exactly.  Anything else is rejected, rather than converted.  In
particular a general `Quaternion` or a `QuatVec` is not taken for a rotor, and `BigFloat`
data is not rounded to the default `T=Float64`; pass `BigFloat` as `T` to work in that type.

# Use from several tasks

Each transform runs on the thread of the task that calls it, including a transform of many
columns of data at once.  Parallelism is therefore a matter of running transforms in several
tasks at the same time, and the three types differ in what that needs.

The `"RS"` and `"Minimal"` objects hold workspace, which every transform overwrites, so one
object must never be used by two tasks at the same time.  Give each task its own: `copy(𝒯)`
returns an independent transform that shares the read-only tables and FFT plans of `𝒯` and
allocates new workspace, at a small fraction of the cost of constructing one.  A vector of
objects indexed by `Threads.threadid()` is no substitute: a task may move to another thread
whenever it yields, and a second task, then running on the thread it left, is given the same
object.  A pattern that works is to divide the data into chunks and to spawn one task per
chunk, each with its own copy:
```julia
chunks = Iterators.partition(maps, cld(length(maps), Threads.nthreads()))
tasks = [Threads.@spawn(let 𝒯ₖ = copy(𝒯); [𝒯ₖ \ f for f ∈ chunk]; end) for chunk ∈ chunks]
mode_weights = reduce(vcat, fetch.(tasks))
```
The `"Matrix"` object holds no workspace — its transforms only read the matrix and its
decomposition, through BLAS and LAPACK — so one object may be used by any number of tasks at
once, and its `copy` is the object itself.

A `deepcopy` of a transform is also independent, and shares the FFT plans of the original as
a `copy` does, since a plan is never modified once it is made.  A transform may be
serialized, as it is when it is sent to another process with `Distributed`; the FFT plans,
which belong to the process that made them, are then made again where the transform is
deserialized.

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
rotor, and changes the sign.  Each ring needs ``N_ϕ ≥ 2ℓₘₐₓ+1`` points, as for an integer
``ℓₘₐₓ``, and the default ``2ℓₘₐₓ+1`` rings of ``2ℓₘₐₓ+1`` points are even numbers when
``ℓₘₐₓ`` is a half-odd-integer.  A quadrature rule symmetric about the equator, such as
Fejér's or Clenshaw–Curtis, then needs only ``2ℓₘₐₓ`` rings, one fewer than it would for an
integer ``ℓₘₐₓ``.
"""
@index_methods function SSHT(
    s::IndexType, ℓₘₐₓ::IndexType, ::Type{T}=Float64; method="RS", kwargs...
) where {T}
    name = method isa Symbol ? String(method) : method
    if name == "RS"
        return SSHTRS(s, ℓₘₐₓ, T; kwargs...)
    elseif name == "Minimal"
        return SSHTMinimal(s, ℓₘₐₓ, T; kwargs...)
    elseif name == "Matrix"
        return SSHTMatrix(s, ℓₘₐₓ, T; kwargs...)
    elseif name == "Direct"
        Base.depwarn(
            "The \"Direct\" s-SHT method has been renamed \"Matrix\"; use method=\"Matrix\".",
            :SSHT
        )
        return SSHTMatrix(s, ℓₘₐₓ, T; kwargs...)
    else
        throw(ArgumentError(
            "Unrecognized s-SHT method $(repr(method)); use \"RS\", \"Minimal\", or \"Matrix\"."
        ))
    end
end


"""
    pixels(𝒯)

Return the spherical coordinates `(θ, ϕ)` at which the transform `𝒯` evaluates functions,
in the order used by the first dimension of the function values.  See also [`rotors`](@ref).
"""
function pixels end


"""
    rotors(𝒯)

Return the `Rotor`s at which the transform `𝒯` evaluates functions, in the order used by
the first dimension of the function values.  The vector is a new one, which may be modified
without affecting the transform.  See also [`pixels`](@ref).
"""
function rotors end

spin(𝒯::SSHT) = 𝒯.s
floattype(::Type{<:SSHT{T}}) where {T} = T
ℓₘₐₓ(𝒯::SSHT) = 𝒯.ℓₘₐₓ
ℓₘᵢₙ(𝒯::SSHT) = abs(𝒯.s)
nmodes(𝒯::SSHT) = Ysize(abs(𝒯.s), 𝒯.ℓₘₐₓ)

function Base.show(io::IO, 𝒯::SSHT)
    print(io, nameof(typeof(𝒯)), "{", floattype(𝒯), "}(s=$(𝒯.s), ℓₘₐₓ=$(𝒯.ℓₘₐₓ))")
end
Base.show(io::IO, ::MIME"text/plain", 𝒯::SSHT) = show(io, 𝒯)

# The argument of the two-argument `*` and `\`: mode weights or function values, as an array
# or, for either, as a `ModeWeights`.  The methods are typed with this rather than left
# untyped, because a method of `*` or `\` whose second argument is untyped matches every
# call whose second argument is inferred as a number or an index, and loading the package
# would then invalidate the compiled code of `Base` that makes such calls.
const SSHTData = Union{AbstractArray, ModeWeights}

# The data of `map2salm`, a map of real or complex function values, and of `salm2map`, mode
# weights as an array or a `ModeWeights`.
const MapArray = AbstractArray{<:Union{Real, Complex}}
const ModesArray = Union{MapArray, ModeWeights}

# The mode-weight vector of an SSHT, as a plain vector of the right length
function check_modes(𝒯::SSHT, f̃)
    n = nmodes(𝒯)
    if size(f̃, 1) != n
        throw(DimensionMismatch(
            "The first dimension of the mode weights has length $(size(f̃, 1)), but "
            * "Ysize(abs(s)=$(abs(𝒯.s)), ℓₘₐₓ=$(𝒯.ℓₘₐₓ)) = $n is required."
        ))
    end
end

# The labels are checked, not just the length: weights of spin -s (or 0) have the same
# length as those of spin s, and would otherwise be synthesized as spin s — or, as the
# output of `ldiv!`, be filled with spin-s weights while keeping the wrong label.  The
# length of the storage is then compared with the labels, because a `ModeWeights` wraps its
# vector without copying it, and the vector may have been resized since the labels were
# checked against it.
function check_modes(𝒯::SSHT, f̃::ModeWeights)
    check_spin(𝒯, f̃)
    if f̃.ℓₘᵢₙ != abs(𝒯.s) || f̃.ℓₘₐₓ != 𝒯.ℓₘₐₓ
        throw(ArgumentError(
            "The ModeWeights have ℓ ∈ $(f̃.ℓₘᵢₙ):$(f̃.ℓₘₐₓ), but the transform requires "
            * "ℓ ∈ $(abs(𝒯.s)):$(𝒯.ℓₘₐₓ).  " * rerange_advice(𝒯)
        ))
    end
    check_storage_length(f̃)
end
function check_spin(𝒯::SSHT, f̃::ModeWeights)
    if spin(f̃) != 𝒯.s
        throw(ArgumentError(
            "The ModeWeights have spin weight s=$(spin(f̃)), but the transform is for "
            * "s=$(𝒯.s)."
        ))
    end
end
rerange_advice(𝒯::SSHT) = (
    "`ModeWeights(w; ℓₘᵢₙ=$(abs(𝒯.s)), ℓₘₐₓ=$(𝒯.ℓₘₐₓ))` copies mode weights `w` into that "
    * "range, filling the modes it adds with zeros."
)

# The mode weights that synthesis reads, as storage in the transform's own layout.  Mode
# weights of the transform's spin weight may cover any range of ℓ up to the transform's
# ℓₘₐₓ: the modes they lack are zero, and those with ℓ < |s|, which belong to no function of
# spin weight s, are ignored.  They are then copied into the transform's range, as
# `ModeWeights(w; ℓₘᵢₙ, ℓₘₐₓ)` would copy them; when the range is the transform's own, the
# storage itself is returned.  Weights beyond the transform's ℓₘₐₓ are refused rather than
# dropped, since that would change the function without saying so.
synthesis_modes(𝒯::SSHT, f̃) = (check_modes(𝒯, f̃); array_view(f̃))
function synthesis_modes(𝒯::SSHT, f̃::ModeWeights)
    check_spin(𝒯, f̃)
    check_storage_length(f̃)
    ℓ₀, ℓ₁ = abs(𝒯.s), 𝒯.ℓₘₐₓ
    if f̃.ℓₘₐₓ > ℓ₁
        throw(ArgumentError(
            "The ModeWeights have ℓ ∈ $(f̃.ℓₘᵢₙ):$(f̃.ℓₘₐₓ), but the transform synthesizes "
            * "ℓ ≤ $ℓ₁ only.  " * rerange_advice(𝒯)
        ))
    end
    (f̃.ℓₘᵢₙ == ℓ₀ && f̃.ℓₘₐₓ == ℓ₁) && return array_view(f̃)
    d = zeros(eltype(f̃), nmodes(𝒯))
    lo = max(ℓ₀, f̃.ℓₘᵢₙ)
    if lo ≤ f̃.ℓₘₐₓ  # the ℓ that the two ranges share, which are contiguous in both
        i = Yindex(lo, -lo, f̃.ℓₘᵢₙ)
        copyto!(d, Yindex(lo, -lo, ℓ₀), array_view(f̃), i, length(f̃) - i + 1)
    end
    d
end

function check_pixels(𝒯::SSHT, f)
    n = npixels(𝒯)
    if size(f, 1) != n
        throw(DimensionMismatch(
            "The first dimension of the function values has length $(size(f, 1)), but the "
            * "transform has $n sample points."
        ))
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

# The results of a transform are complex, so storage that receives them — the output of a
# three-argument `mul!` or `ldiv!`, or the argument of an operation in place — must be able
# to hold complex numbers.  Real storage would otherwise fail with an `InexactError` from
# deep inside the transform, or part of the way through writing its result.  The type tested
# against is a constant, since a `Complex{<:AbstractFloat}` written in the body would be
# built anew, and allocated, on every call.
const ComplexFloat = Complex{<:AbstractFloat}
function check_complex_output(x, name)
    S = eltype(array_view(x))
    if !(S <: ComplexFloat)
        throw(ArgumentError(
            "The results of the transform are complex, so the output `$name` must hold complex "
            * "floating-point numbers, but its element type is $S."
        ))
    end
end
function check_complex_storage(x, operation, alternative=operation)
    S = eltype(array_view(x))
    if !(S <: ComplexFloat)
        throw(ArgumentError(
            "`$operation` acts in place, storing its complex results in the storage of its "
            * "argument, whose element type is $S.  Pass complex floating-point data, or use "
            * "`$alternative` with a transform constructed with `inplace=false`, which "
            * "accepts real data."
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

# What the in-place `\` and `*` return, once they have overwritten their argument's storage.
# Analysis of one-dimensional data returns its mode weights as a `ModeWeights`, as every
# other analysis does, but wrapping that same storage, so that nothing is copied and they
# are still indexed by (ℓ, m); returned as the bare vector, `(𝒯 \ f)[ℓ, m]` would silently
# read a linear index.  Synthesis returns the storage itself, as a plain array, and never a
# `ModeWeights` that held the input, whose labels would no longer describe its contents.
#
# The storage of an in-place transform holds exactly `nmodes(𝒯)` mode weights, so the
# `ModeWeights` wraps all of it.  A least-squares `SSHTMatrix` has more points than modes,
# and its two-argument `ldiv!` leaves the solution in the first `nmodes(𝒯)` entries of the
# longer argument; only those are labelled, or returned, rather than the whole array with
# its residual components.  Only a transform that is not in place can be of that kind, so
# the type of the result is known from the type of the transform wherever it can be.
function in_place_modes(𝒯::SSHT, ff̃)
    d = array_view(ff̃)
    ndims(d) == 1 ? ModeWeights(d, 𝒯.s, abs(𝒯.s), 𝒯.ℓₘₐₓ) : d
end
function solution_modes(𝒯::SSHT, ff̃)
    d = array_view(ff̃)
    n = nmodes(𝒯)
    if ndims(d) == 1
        s, ℓₘᵢₙ, ℓₘₐₓ = 𝒯.s, abs(𝒯.s), 𝒯.ℓₘₐₓ
        length(d) == n ? ModeWeights(d, s, ℓₘᵢₙ, ℓₘₐₓ) : mode_weights_view(d, s, ℓₘᵢₙ, ℓₘₐₓ)
    else
        size(d, 1) == n ? d : view(d, 1:n, ntuple(_ -> Colon(), ndims(d) - 1)...)
    end
end

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

# The sample points of "Matrix" and "Minimal" can be badly conditioned, and nothing in an
# individual transform reveals it, so their constructors measure it once, by synthesizing
# and analyzing a fixed set of unit weights with quasi-random phases — about the cost of one
# pair of transforms — and warn when fewer than half the digits of `T` survive.  The
# analysis is `𝒯 \ f`, which acts in place only where the transform's option says so,
# rather than the two-argument `ldiv!`, which for "Matrix" solves in place whatever the
# option, and which a `decomposition` given to a transform that does not act in place need
# not support.  The values `f` are scratch, and may be overwritten.
function warn_if_inaccurate(𝒯::SSHT{T}, method) where {T}
    φ = (√5 - 1) / 2
    f̃ = [cis(T(2π) * T(mod(i * φ, 1))) for i ∈ 1:nmodes(𝒯)]
    f = Vector{Complex{T}}(undef, npixels(𝒯))
    mul!(f, 𝒯, f̃)
    maxerror = maximum(abs, array_view(𝒯 \ f) - f̃)
    if !(maxerror ≤ √eps(T))  # also catches NaN
        @warn (
            "The \"$method\" s-SHT with s=$(𝒯.s), ℓₘₐₓ=$(𝒯.ℓₘₐₓ) and T=$T is inaccurate: a "
            * "round trip of unit mode weights has a maximum error of "
            * "$(round(Float64(maxerror), sigdigits=2)).  Its sample points are badly "
            * "conditioned at this ℓₘₐₓ; the \"RS\" method (the default) is accurate here."
        )
    end
    nothing
end


include("matrix.jl")
include("rings.jl")
include("rs.jl")
include("minimal.jl")
