"""
    ModeWeights(data, s=0; ℓₘᵢₙ=abs(s))
    ModeWeights(data, s, ℓₘᵢₙ, ℓₘₐₓ)
    ModeWeights{T}(undef, s, ℓₘᵢₙ, ℓₘₐₓ)
    ModeWeights{T}(undef, s, ℓₘₐₓ)
    ModeWeights(w::ModeWeights; ℓₘᵢₙ=abs(spin(w)), ℓₘₐₓ=ℓₘₐₓ(w))

Vector of mode weights ``f_{ℓ,m}`` of a spin-weighted function ``f = \\sum_{ℓ,m} f_{ℓ,m}\\,
{}_sY_{ℓ,m}``, stored in the canonical ordering `[f(ℓ, m) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ]`
(see [`Yindex`](@ref)), together with the spin weight `s` and the range of ``ℓ``.

A `ModeWeights` is an [`AbstractModeContainer`](@ref
SphericalFunctions.AbstractModeContainer), not an `AbstractVector`; [`array_view`](@ref)
gives the flat 1-based storage, which is what the transforms and the operator matrices take.
It behaves as the vector of the weights in other respects: `w[i]`, for a position `i` in
storage order, reads or writes one weight, `w[:]`, `w[r]` and `w[v]`, for a range `r` or a
vector `v` of positions, read or write several, and iteration, `length`, `keys`, `sum`,
`maximum` and the other reductions see the weights one by one.

The arithmetic of mode weights as the weights of functions keeps the labels: `a + b` and `a
- b` (and their broadcast forms) when the labels agree, `complex.(a, b)` from real and
imaginary parts whose labels agree, and products and quotients with numbers, or elementwise
with a plain vector of factors (a diagonal operator, such as a filter).  So do `zero(w)`,
`fill!(w, x)`, `rmul!(w, a)`, `lmul!(a, w)`, `I * w`, `copyto!(w, v)` from a plain vector of
one weight per mode, and `copyto!`, `copy!`, `axpy!` and `axpby!` between weights whose
labels agree.  Anything else — adding a number to every weight, the elementwise product of
two sets of weights, `conj.(w)`, `abs2.(w)` — would label numbers that are not the mode
weights of any function, and is an error; arithmetic on the raw numbers goes through
`array_view(w)`.  A broadcast whose result is labelled is a `ModeWeights`, and one that
would extend the weights to a longer vector is refused.  `map` returns plain numbers, and so
do `similar(w, n)` and `similar(w, T, dims)`, which give storage of the size asked for,
while `similar(w)` and `similar(w, T)` keep the labels.

`dot(a, b)` is the inner product of the two functions, and refuses two `ModeWeights` whose
labels differ.  `w'` and `transpose(w)` are those of the storage, so `a' * b` is a product
with a plain matrix, which is refused, as described below.  `a == b`, `isequal(a, b)` and `a
≈ b` compare the labels as well as the numbers, and are `false` when the labels differ;
against a plain vector they compare the numbers alone, and so `hash(w)` is the hash of the
numbers.  In addition
- `w[ℓ, m]` reads or writes the weight of mode ``(ℓ, m)``,
- `w[ℓ, :]` is a [`DegreeBlock`](@ref) view of the weights for one ``ℓ``, indexed by
  `m ∈ -ℓ:ℓ`,
- `modes(w)` is the vector of `(ℓ, m)` pairs in storage order,
- `spin(w)`, `ℓₘᵢₙ(w)`, `ℓₘₐₓ(w)` are the parameters, and `parent(w)` is the storage,
- the differential operators [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref),
  [`Lx`](@ref), [`Ly`](@ref), [`R²`](@ref), [`Rz`](@ref), [`R₊`](@ref), [`R₋`](@ref),
  [`ð`](@ref), [`ð̄`](@ref) give a new `ModeWeights` with the spin weight adjusted where
  appropriate and the range of ``ℓ`` of `w`, written either `ð * w` or `ð(w)`.  These build
  no matrix: the operator is applied by a loop, so the only allocation is the result, and
  `mul!(w′, ð, w)` into a correctly labelled destination allocates nothing at all (`w′` may
  also be a bare vector at least as long as the result, which then comes back as a
  `ModeWeights` over it),
- a product with a plain matrix, `A * w`, is refused, because nothing about a matrix of
  numbers says what spin weight or range of ``ℓ`` it was built for, so the labels could not
  be checked: `A * array_view(w)` is the product with the raw numbers, `sYlm(R⃗, ℓₘₐₓ, s) *
  w` evaluates the function at rotors, and `op * w` or `mul!(w′, op, w)` applies a
  differential operator with the labels checked, and
- `w(R)` evaluates the function at the rotor `R` (see [`sYlm`](@ref)).

When constructed from `data` alone, `ℓₘₐₓ` is deduced from `length(data)` and `ℓₘᵢₙ`, which
must match exactly.  The `data` vector is used as storage, not copied, so it must keep its
length for as long as the `ModeWeights` is in use; the operators, `w[ℓ, :]`, and `w[ℓ, m]`
outside `@inbounds` throw a `DimensionMismatch` when they find that it has been resized.
The `undef` forms allocate uninitialized storage of type `T` instead, and in the
three-argument form `ℓₘᵢₙ` defaults to `abs(s)`, as it does when `data` is given.  The
keyword `ℓₘᵢₙ` may also be spelled `ell_min`.  `show` writes the labels, as
`ModeWeights{ComplexF64} with s=2, ℓ ∈ 2:8`, and the text/plain form adds the weights.

# Changing the range of ``ℓ``

The spin-changing operators keep the range of ``ℓ`` of their input, so the result of `ð * w`
for `w` with ``ℓₘᵢₙ = |s|`` starts at ``ℓ = |s|``, not at ``|s + 1|`` as the weights of its
new spin weight usually do.  Weights of one spin weight but different ranges of ``ℓ`` cannot
be added or compared, and the destination of `mul!` must have exactly the labels of the
result.  `ModeWeights(w; ℓₘᵢₙ, ℓₘₐₓ)` copies `w` into another range of ``ℓ``, with the same
spin weight: it fills with zeros every ``ℓ`` that `w` does not hold, drops the entries with
``ℓ < |s|``, which belong to no harmonic, and drops those above the `ℓₘₐₓ` requested.  With
the defaults, it copies `w` into ``|s| ≤ ℓ ≤ ℓₘₐₓ(w)``, the range the transforms' analysis
produces, so it is a copy only when `w` already has that range.  The keywords may also be
spelled `ell_min` and `ell_max`, and must be indices of the kind of `w`'s own.

```julia
w₂ = ModeWeights(ð * (SSHT(1, ℓₘₐₓ) \\ f); ℓₘᵢₙ=2)   # the spin-2 weights on ℓ ∈ 2:ℓₘₐₓ
w₂ ≈ SSHT(2, ℓₘₐₓ) \\ g
```

# Several functions at once

A transform applied to several functions at once, as the columns of a matrix, returns a
plain array whose columns are the mode weights, and a `ModeWeights` holds one vector of
them.  `ModeWeights(view(W, :, j), s)` labels column `j` without copying it, so that an
operator or a rotation can be applied column by column with no allocation:

```julia
for j ∈ axes(W, 2)
    mul!(view(W′, :, j), ð, ModeWeights(view(W, :, j), s))
end
```

The weights are real or complex numbers, and a quaternion element type is refused.  A
`Rotor` is a quaternion, so `R * w` is an error; the rotation of the function by `R` is
`D(R, ℓₘₐₓ(w)) * w` (see [`D`](@ref)).

# Half-integer indices

The spin weight and the range of ``ℓ`` may be half-integers, passed as `Rational`s with
denominator 2 — as in `ModeWeights(data, 1//2)` or `ModeWeights{T}(undef, 1//2, 1//2, 7//2)`
— or as [`HalfOddInteger`](@ref)s, in which case every ``ℓ`` and ``m`` of the ordering is a
half-odd-integer, and `ℓₘᵢₙ` may be as small as `1//2`.  The indices in one call must all be
of one kind, integers of type `Int` or half-odd-integers; a call that mixes them, such as
`ModeWeights(data, 1//2, 0, 7//2)`, is an error, as is an integer index of another type,
such as an `Int32`, which is to be converted with `Int` first.  The parameters are stored as
`HalfOddInteger`s, which is also what `modes(w)` and the axis of `w[ℓ, :]` are made of;
`w[ℓ, m]` accepts either spelling.  For such a `w`, `w[ℓ, :]` is a [`DegreeBlock`](@ref),
indexed by `m ∈ -ℓ:ℓ`, exactly as it is for integer indices.  The natural indices of `w[ℓ,
m]` and `w[ℓ, :]` obey the same rules as those of the constructors: an integer `w` takes
`Int`s, and an index of another kind or integer type is refused with the reason.
"""
struct ModeWeights{T, IT<:IntegerHalf, V<:AbstractVector{T}} <: AbstractModeContainer{T, IT}
    data::V
    s::IT
    ℓₘᵢₙ::IT
    ℓₘₐₓ::IT
    # These checks hold for either kind of index: the floor of ℓₘᵢₙ is 0 for integers and 1/2
    # for half-odd-integers, and `ℓₘᵢₙ < 0` is the right test for both, since no
    # half-odd-integer lies between 0 and 1/2.  The index methods make this the one
    # constructor that every other reaches, with three indices of one kind.
    @index_methods function ModeWeights(
        data::V, s::IT, ℓₘᵢₙ::IT, ℓₘₐₓ::IT
    ) where {T, IT<:IndexType, V<:AbstractVector{T}}
        Base.require_one_based_indexing(data)
        # A `Rotor` is a `Number`, so `R * w` and `R .* w` would otherwise give quaternions
        # under `w`'s labels, although quaternion-valued weights are the weights of no function
        # that this package can evaluate: quaternions do not commute with the complex
        # harmonics.  Refusing the element type here closes every route to them at once.
        # (`Union{}` is a subtype of every type, and is the element type of some empty results.)
        if T !== Union{} && T <: AbstractQuaternion
            throw(ArgumentError(
                "Mode weights are real or complex numbers, not quaternions of type $T.  To "
                * "rotate the function with a rotor `R`, use `D(R, ℓₘₐₓ(w)) * w`."
            ))
        end
        if ℓₘᵢₙ < 0
            throw(ArgumentError("ℓₘᵢₙ=$ℓₘᵢₙ must be non-negative."))
        end
        if ℓₘₐₓ < ℓₘᵢₙ - 1
            throw(ArgumentError("ℓₘₐₓ=$ℓₘₐₓ must be at least ℓₘᵢₙ-1=$(ℓₘᵢₙ-1)."))
        end
        if length(data) != Ysize(ℓₘᵢₙ, ℓₘₐₓ)
            throw(ArgumentError(
                "The data has length $(length(data)), but Ysize(ℓₘᵢₙ=$ℓₘᵢₙ, ℓₘₐₓ=$ℓₘₐₓ) "
                * "= $(Ysize(ℓₘᵢₙ, ℓₘₐₓ))."
            ))
        end
        new{T, IT, V}(data, s, ℓₘᵢₙ, ℓₘₐₓ)
    end
end

# The outer constructors are index methods too, so that each is reached with indices of one
# kind, and a keyword `ℓₘᵢₙ` is normalized against the kind of `s` in either spelling.  The
# default `abs(s)` of the keyword is evaluated where `s` is already an `Int` or a
# `HalfOddInteger`.
@index_methods function ModeWeights(
    data::AbstractVector, s::IndexType=0; ell_min::IndexType=abs(s), ℓₘᵢₙ::IndexType=ell_min
)
    deduced_mode_weights(data, s, ℓₘᵢₙ)
end
@index_methods function ModeWeights{T}(
    ::UndefInitializer, s::IT, ℓₘᵢₙ::IT, ℓₘₐₓ::IT
) where {T, IT<:IndexType}
    ModeWeights(Vector{T}(undef, Ysize(ℓₘᵢₙ, ℓₘₐₓ)), s, ℓₘᵢₙ, ℓₘₐₓ)
end
@index_methods ModeWeights{T}(::UndefInitializer, s::IndexType, ℓₘₐₓ::IndexType) where {T} =
    ModeWeights{T}(undef, s, abs(s), ℓₘₐₓ)

# A matrix of mode weights, such as the result of a transform applied to several functions at
# once, holds one set of weights per column; the constructors label a vector, and this says
# how to label a column.
function ModeWeights(data::AbstractMatrix, args...; kwargs...)
    throw(ArgumentError(
        "A `ModeWeights` holds the mode weights of one function, as a vector, but this is a "
        * "$(join(size(data), "×")) matrix.  To label the weights in column j without copying "
        * "them, use `ModeWeights(view(data, :, j), s)`."
    ))
end

# Deduce ℓₘₐₓ from the length of the data, given `s` and `ℓₘᵢₙ` of one kind.  The kind of the
# indices selects the method, since the relation between the length and ℓₘₐₓ is written
# differently for the two.
function deduced_mode_weights(data::AbstractVector, s::Int, ℓₘᵢₙ::Int)
    # Deduce ℓₘₐₓ from (ℓₘₐₓ+1)² = length + ℓₘᵢₙ²
    N = length(data) + ℓₘᵢₙ^2
    ℓₘₐₓ = isqrt(N) - 1
    if (ℓₘₐₓ + 1)^2 != N
        throw(ArgumentError(
            "The data has length $(length(data)), which is not Ysize(ℓₘᵢₙ=$ℓₘᵢₙ, ℓₘₐₓ) "
            * "for any ℓₘₐₓ."
        ))
    end
    ModeWeights(data, s, ℓₘᵢₙ, ℓₘₐₓ)
end
function deduced_mode_weights(data::AbstractVector, s::HalfOddInteger, ℓₘᵢₙ::HalfOddInteger)
    # Deduce ℓₘₐₓ from (2ℓₘₐₓ+2)² = 4⋅length + (2ℓₘᵢₙ)², which is `Ysize` on the doubled indices.
    # The root must square back exactly.  It is then automatically odd, as 2ℓₘₐₓ+2 must be for
    # a half-odd ℓₘₐₓ, because 4⋅length + (2ℓₘᵢₙ)² is odd whenever 2ℓₘᵢₙ is; a length that
    # corresponds to an integer ℓₘₐₓ has no exact root here and is refused by the one test.
    # Everything here is `Int` arithmetic on numerators, and the result is built from its
    # numerator at the end.
    N = 4length(data) + (2ℓₘᵢₙ)^2
    r = isqrt(N)
    if r^2 != N
        throw(ArgumentError(
            "The data has length $(length(data)), which is not Ysize(ℓₘᵢₙ=$ℓₘᵢₙ, ℓₘₐₓ) "
            * "for any half-odd-integer ℓₘₐₓ."
        ))
    end
    ModeWeights(data, s, ℓₘᵢₙ, unsafe_half_odd_integer(r - 2))
end

# Copying into another range of ℓ.  There is no positional index to dispatch on, so the two
# keywords are normalized against the kind of `w` by hand, with a message that names `w`
# rather than the positional indices of the call.  Each ASCII spelling is normalized under
# its own name before it becomes the default of the Unicode one, so that a refusal names the
# keyword the caller wrote.  The default of `ell_max` reads the field, since the keyword
# `ℓₘₐₓ` shadows the accessor inside the method.
"""
    ModeWeights(w::ModeWeights; ℓₘᵢₙ=abs(spin(w)), ℓₘₐₓ=ℓₘₐₓ(w))

A copy of the mode weights `w` over the range `ℓₘᵢₙ:ℓₘₐₓ` of ``ℓ``, with the same spin
weight, zero where `w` holds no weight.

Every ``ℓ`` of the new range that `w` does not hold is filled with zeros, the entries of `w`
with ``ℓ < |s|``, which belong to no harmonic, are dropped, and so are those above the
`ℓₘₐₓ` requested.  The defaults give the range ``|s| ≤ ℓ ≤ ℓₘₐₓ(w)``, which is what the
transforms' analysis produces.  The result has new storage, of the element type of `w`, even
when the range is that of `w` itself.  The keywords may also be spelled `ell_min` and
`ell_max`; they must be indices of the kind of `w`'s own, integers of type `Int` or
half-odd-integers, each a [`HalfOddInteger`](@ref) or a `Rational{Int}` with denominator 2.

This is the way to bring weights of one spin weight but different ranges of ``ℓ`` together:
the spin-changing operators keep the range of their input, while `+`, `≈`, `dot` and the
destination of `mul!` require the labels to agree exactly.

```julia
w = SSHT(1, 8) \\ f                  # s = 1, ℓ ∈ 1:8
ðw = ModeWeights(ð * w; ℓₘᵢₙ=2)      # s = 2, ℓ ∈ 2:8, as SSHT(2, 8) \\ g is
ð̄w = ModeWeights(ð̄ * w; ℓₘᵢₙ=0)      # s = 0, ℓ ∈ 0:8, with ℓ = 0 filled with zeros
```
"""
function ModeWeights(
    w::ModeWeights;
    ell_min=abs(w.s), ℓₘᵢₙ=range_keyword(w, :ell_min, ell_min),
    ell_max=w.ℓₘₐₓ, ℓₘₐₓ=range_keyword(w, :ell_max, ell_max)
)
    lo = range_keyword(w, :ℓₘᵢₙ, ℓₘᵢₙ)
    hi = range_keyword(w, :ℓₘₐₓ, ℓₘₐₓ)
    check_storage_length(w)
    w′ = ModeWeights{eltype(w)}(undef, w.s, lo, hi)
    fill!(w′.data, zero(eltype(w)))
    # The ℓ ∈ a:b that both hold, leaving out those below |s|, are one contiguous run of modes
    # in each storage, so one copy moves them all.
    a, b = max(lo, w.ℓₘᵢₙ, abs(w.s)), min(hi, w.ℓₘₐₓ)
    if a ≤ b
        i₀, j₀ = Yindex(a, -a, lo), Yindex(a, -a, w.ℓₘᵢₙ)
        copyto!(w′.data, i₀, w.data, j₀, Yindex(b, b, lo) - i₀ + 1)
    end
    w′
end

# One keyword bound of `ModeWeights(w; …)`, as an index of `w`'s own kind.
function range_keyword(w::ModeWeights{T, IT}, name::Symbol, x) where {T, IT}
    if IT === Int
        x isa Int && return x
    else
        x isa HalfOddInteger && return x
        x isa Rational{Int} && denominator(x) == 2 && return half_odd_index(x)
    end
    kind = IT === Int ? "integers of type `Int`, like 3" : (
        "half-odd-integers, each a `HalfOddInteger` or a `Rational{Int}` with denominator 2, "
        * "like 7//2"
    )
    message = (
        "The keyword argument `$name` of `ModeWeights(w; …)` must be an index of the kind of "
        * "`w`'s own, which are $kind; got $name = $(typed_repr(x))."
    )
    index_kind(x) === nothing && (message *= "  " * index_problem(x))
    throw(ArgumentError(message))
end

Base.parent(w::ModeWeights) = w.data

"""
    spin(w)

The spin weight of a [`ModeWeights`](@ref) vector, of an [`SSHT`](@ref) transform, or of an
[`sYlmCalculator`](@ref) built for a single one.  A calculator built for a range of spin
weights has no single value to report, so it has no method here; ask it for [`spins`](@ref
SphericalFunctions.spins) instead, which answers for either kind.

A function of spin weight ``s`` has ``R_z f = s f``, and is expanded in the harmonics
``{}_{s}Y_{ℓ,m}`` with ``ℓ ≥ |s|``.  The spin weight is kept alongside the numbers because
nothing about the numbers themselves reveals it.

```jldoctest
julia> using SphericalFunctions

julia> spin(ModeWeights(zeros(ComplexF64, 21), -2))
-2

julia> spin(SSHT(1, 4))
1
```

See also [`modes`](@ref), [`ModeWeights`](@ref), [`SSHT`](@ref), and
[`spins`](@ref SphericalFunctions.spins).
"""
function spin end

spin(w::ModeWeights) = w.s

"""
    modes(w::ModeWeights)

The `(ℓ, m)` pairs of `w`, in storage order (see [`Yrange`](@ref)).
"""
modes(w::ModeWeights) = Yrange(w.ℓₘᵢₙ, w.ℓₘₐₓ)

# The array-like interface, written out rather than inherited.  A `ModeWeights` is an
# [`AbstractModeContainer`](@ref) like the rest, not an `AbstractVector`, so `op * w` and `w
# .+ 1` do not come for free from the generic machinery; these are the methods that supply
# the useful part of that behavior.  Forgoing the subtyping costs less than it appears to:
# the transforms in `ssht/` reach for the raw storage before every `mul!` and `ldiv!`
# anyway, through [`array_view`](@ref).
Base.size(w::ModeWeights) = size(w.data)
Base.size(w::ModeWeights, d::Integer) = d ≤ 1 ? size(w)[d] : 1
Base.length(w::ModeWeights) = length(w.data)
Base.axes(w::ModeWeights) = axes(w.data)
Base.axes(w::ModeWeights, d::Integer) = d ≤ 1 ? axes(w)[d] : Base.OneTo(1)
Base.ndims(::ModeWeights) = 1
Base.ndims(::Type{<:ModeWeights}) = 1
Base.firstindex(w::ModeWeights) = firstindex(w.data)
Base.lastindex(w::ModeWeights) = lastindex(w.data)
Base.iterate(w::ModeWeights, state...) = iterate(w.data, state...)
Base.keys(w::ModeWeights) = keys(w.data)
Base.eachindex(w::ModeWeights) = eachindex(w.data)
# Linear indexing is by position in the storage, as for any vector, so it takes the
# positions that the storage takes; the natural indices `w[ℓ, m]` below are what obey the
# index rules.  Any index but a single position gives plain numbers, as `map` does.
@propagate_inbounds Base.getindex(w::ModeWeights, i::Integer) = w.data[i]
@propagate_inbounds Base.getindex(w::ModeWeights, I::Union{Colon, AbstractVector{<:Integer}}) =
    w.data[I]
@propagate_inbounds Base.setindex!(w::ModeWeights, v, i::Integer) = (w.data[i] = v)
@propagate_inbounds Base.setindex!(
    w::ModeWeights, v, I::Union{Colon, AbstractVector{<:Integer}}
) = (w.data[I] = v)
# A position may also be a one-dimensional `CartesianIndex`, which is how Julia 1.10's
# broadcasting reads each of its arguments, arrays or not.
@propagate_inbounds Base.getindex(w::ModeWeights, I::CartesianIndex{1}) = w.data[I]
@propagate_inbounds Base.setindex!(w::ModeWeights, v, I::CartesianIndex{1}) = (w.data[I] = v)
Base.similar(w::ModeWeights) = ModeWeights(similar(w.data), w.s, w.ℓₘᵢₙ, w.ℓₘₐₓ)
Base.similar(w::ModeWeights, ::Type{S}) where {S} =
    ModeWeights(similar(w.data, S), w.s, w.ℓₘᵢₙ, w.ℓₘₐₓ)
# A size, even the size of `w`, asks for storage rather than for weights, so these give plain
# vectors and arrays; that keeps the type returned independent of the size asked for.
Base.similar(w::ModeWeights, n::Integer) = similar(w.data, n)
Base.similar(w::ModeWeights, ::Type{S}, n::Integer) where {S} = similar(w.data, S, n)
Base.similar(w::ModeWeights, dims::Dims) = similar(w.data, dims)
Base.similar(w::ModeWeights, ::Type{S}, dims::Dims) where {S} = similar(w.data, S, dims)
Base.collect(w::ModeWeights) = collect(w.data)
Base.Array(w::ModeWeights) = collect(w.data)
Base.Vector(w::ModeWeights) = collect(w.data)
Base.:(==)(a::ModeWeights, b::ModeWeights) = same_labels(a, b) && a.data == b.data
Base.isequal(a::ModeWeights, b::ModeWeights) = same_labels(a, b) && isequal(a.data, b.data)
# `isequal(w, parent(w))` holds, so the hash must be that of the numbers alone; weights under
# different labels then share a hash, which is allowed.
Base.hash(w::ModeWeights, h::UInt) = hash(w.data, h)

# A product with a plain matrix is refused.  A matrix of numbers says nothing about the spin
# weight or the range of ℓ it was built for — `sYlm_matrix(R⃗, ℓₘₐₓ, s)` has the same shape
# for `s` and `-s`, and `ð(s, ℓₘᵢₙ, ℓₘₐₓ)` is a `Diagonal` like `L²(s, ℓₘᵢₙ, ℓₘₐₓ)` — so the
# labels of `w` could not be checked, and weights of the wrong spin weight would give a
# wrong answer with no error.  Each labelled product checks them instead, and `array_view`
# is the explicit route to the raw numbers.  The `AbstractMatrix` methods cover every matrix
# type, the banded ones and the adjoint of a `ModeWeights`, `w'`, included: no method of `*`
# or `\` for a `Diagonal`, `Bidiagonal` or `Tridiagonal` takes an untyped second argument.
# The product in place, `mul!(y, A, w)`, reaches the five-argument form, which is refused in
# the same way.  The messages are shared by the four methods.
function matrix_product_error(product, raw)
    ArgumentError(
        "A plain matrix cannot be checked against the labels of mode weights, so `$product` is "
        * "refused.  Use `$raw` for arithmetic on the raw numbers, `sYlm(R⃗, ℓₘₐₓ, s) * w` to "
        * "evaluate the function at rotors, `op * w` or `mul!(w′, op, w)` to apply a "
        * "differential operator, and `dot(w₁, w₂)` for the inner product of two sets of mode "
        * "weights, all of which check the labels."
    )
end
Base.:*(A::AbstractMatrix, w::ModeWeights) =
    throw(matrix_product_error("A * w", "A * array_view(w)"))
Base.:*(w::ModeWeights, A::AbstractMatrix) =
    throw(matrix_product_error("w * A", "array_view(w) * A"))
Base.:\(A::AbstractMatrix, w::ModeWeights) =
    throw(matrix_product_error("A \\ w", "A \\ array_view(w)"))
LinearAlgebra.mul!(y, A::AbstractMatrix, w::ModeWeights, α::Number, β::Number) =
    throw(matrix_product_error("mul!(y, A, w)", "mul!(y, A, array_view(w))"))
Base.copy(w::ModeWeights) = ModeWeights(copy(w.data), w.s, w.ℓₘᵢₙ, w.ℓₘₐₓ)

# Broadcasting over mode weights is allowed only where the result is again the mode weights
# of a function, with labels that can be stated: sums and differences of weights with the
# same spin weight and range of ℓ, and products and quotients of weights with numbers — or
# with plain vectors, one factor per mode, which is a diagonal operator such as a filter.
# The result then has those labels.  Anything else that involves mode weights is refused,
# because it would put a label on numbers it does not describe: the mode weights of a
# product of two functions are not the product of their mode weights, those of the complex
# conjugate are not the conjugates, and `abs2.(w)` is not the mode weights of anything.
# Such arithmetic on the raw numbers is still available through `array_view(w)`.  A
# `Broadcast.ArrayStyle` would be the usual way to get a wrapped result, but it is available
# only to an `AbstractArray`; a style of this type's own does the same job, given a
# `broadcastable` that hands back the container rather than `collect`ing it, and the `axes`
# and linear `getindex` defined above.  Where an array of another style takes part, such as
# a `StaticArray`, the two styles conflict and the result is a plain array, which has no
# labels to check.
struct ModeWeightsStyle <: Broadcast.AbstractArrayStyle{1} end
ModeWeightsStyle(::Val{0}) = ModeWeightsStyle()
ModeWeightsStyle(::Val{1}) = ModeWeightsStyle()
ModeWeightsStyle(::Val{N}) where {N} = Broadcast.DefaultArrayStyle{N}()
Base.BroadcastStyle(::Type{<:ModeWeights}) = ModeWeightsStyle()
Base.Broadcast.broadcastable(w::ModeWeights) = w
# The result of a broadcast of this style is always labelled, so its type does not depend on
# the lengths involved.  Its length can differ from that of the labels only where weights of
# a single mode, ℓₘₐₓ = ℓₘᵢₙ = 0, are extended to a longer vector of factors, which is
# refused.
function Base.similar(bc::Broadcast.Broadcasted{ModeWeightsStyle}, ::Type{S}) where {S}
    s, ℓₘᵢₙ, ℓₘₐₓ = broadcast_label(bc)
    if length(bc) != Ysize(ℓₘᵢₙ, ℓₘₐₓ)
        throw(DimensionMismatch(
            "This broadcast gives $(length(bc)) values, but the mode weights it involves, "
            * "with s=$s and ℓ ∈ $ℓₘᵢₙ:$ℓₘₐₓ, have $(Ysize(ℓₘᵢₙ, ℓₘₐₓ)) modes.  A plain "
            * "vector combined with mode weights must hold one value per mode; use "
            * "`array_view(w)` for arithmetic on the raw numbers."
        ))
    end
    ModeWeights(similar(Vector{S}, axes(bc)), s, ℓₘᵢₙ, ℓₘₐₓ)
end

# The labels `(s, ℓₘᵢₙ, ℓₘₐₓ)` of what a broadcast expression computes, or `nothing` if it
# involves no mode weights; an expression that combines mode weights in a way that does not
# give mode weights is an error.  See the comment above.
broadcast_label(w::ModeWeights) = (w.s, w.ℓₘᵢₙ, w.ℓₘₐₓ)
broadcast_label(x) = nothing
function broadcast_label(bc::Broadcast.Broadcasted)
    labels = map(broadcast_label, bc.args)
    labelled = filter(!isnothing, labels)
    isempty(labelled) && return nothing
    label, f = first(labelled), bc.f
    if f === (+) || f === (-)
        check_termwise(bc, labels)
    elseif f === (*)
        if length(labelled) > 1
            throw(ArgumentError(
                "The mode weights of a product of functions are not the product of their mode "
                * "weights.  Mode weights may be multiplied by numbers, or elementwise by a "
                * "plain vector of factors; use `array_view(w)` for arithmetic on the raw numbers."
            ))
        end
    elseif f === (/) || f === (\)
        numerator_position = f === (/) ? 1 : 2
        if length(bc.args) != 2 || labels[3 - numerator_position] !== nothing
            throw(ArgumentError(
                "Mode weights may be divided by numbers, or elementwise by a plain vector of "
                * "factors, but nothing may be divided by mode weights; use `array_view(w)` for "
                * "arithmetic on the raw numbers."
            ))
        end
    elseif length(bc.args) == 1 && (
        f === identity || f === float || f === complex || (f isa Type && f <: Number)
    )
        # a copy, or a change of number type, holds the same modes
    elseif length(bc.args) == 2 && (f === complex || (f isa Type && f <: Complex))
        # complex weights from their real and imaginary parts, which are then the weights of
        # one function only if both parts are
        check_termwise(bc, labels)
    else
        throw(ArgumentError(
            "Broadcasting `$f` over mode weights does not give the mode weights of any function "
            * "with the same labels; only sums, differences, complex weights built from real "
            * "and imaginary parts with the same labels, and products and quotients with "
            * "numbers or plain vectors of factors do.  Use `array_view(w)` for arithmetic on "
            * "the raw numbers."
        ))
    end
    label
end

# The rule for the operations that combine mode weights term by term into the weights of one
# function — sums, differences, and complex weights built from their real and imaginary
# parts.  The mode weights among the arguments must have the same labels, and every other
# argument of at most one dimension must hold one value per mode.  That refuses a number,
# however it is written — a literal, a `Ref`, a zero-dimensional array, or a
# zero-dimensional broadcast such as the `a .* b` of `w .+ a .* b` — and an array or tuple
# too short to hold one value per mode, such as `[1.0]` or `(1,)`, which broadcasting would
# extend to every mode just as it does a number.  An argument of two or more dimensions
# makes the result a matrix, such as the outer sum `w .+ transpose(w)`, which is never
# labelled, so it is left alone.  The messages are built only on the branches that throw
# them, so that a broadcast that passes allocates nothing here.
function check_termwise(bc::Broadcast.Broadcasted, labels)
    labelled = filter(!isnothing, labels)
    label = first(labelled)
    sum_or_difference = bc.f === (+) || bc.f === (-)
    if any(!=(label), labelled)
        what = if sum_or_difference
            "added or subtracted"
        else
            "combined as the real and imaginary parts of complex weights"
        end
        throw(ArgumentError(
            "Mode weights can be $what only when their labels agree; got "
            * join(("s=$(l[1]), ℓ ∈ $(l[2]):$(l[3])" for l ∈ labelled), " and ") * "."
        ))
    end
    n = Ysize(label[2], label[3])
    if any(map((x, l) -> l === nothing && extended_to_every_mode(x, n), bc.args, labels))
        what = if sum_or_difference
            "Adding a number to every mode weight"
        else
            "Using one number as the real or imaginary part of every mode weight"
        end
        throw(ArgumentError(
            "$what does not give the mode weights of any function, whether the number is "
            * "written as such, computed in the same broadcast, or given as an array or tuple "
            * "that broadcasting extends to every mode; an array combined with mode weights "
            * "must hold one value per mode.  Use `array_view(w)` for arithmetic on the raw "
            * "numbers."
        ))
    end
    nothing
end
function extended_to_every_mode(x, n)
    ax = axes(x)
    length(ax) == 0 || (length(ax) == 1 && length(only(ax)) != n)
end

# Writing into mode weights with `.=` checks the labels in the same way: whatever the
# right-hand side computes must be what the destination's labels say it holds.  (The other
# containers go through the methods in `array_view.jl`.)
function check_broadcast_destination(dest::ModeWeights, bc)
    label = broadcast_label(bc)
    if label !== nothing && label != broadcast_label(dest)
        throw(ArgumentError(
            "The destination holds s=$(dest.s), ℓ ∈ $(dest.ℓₘᵢₙ):$(dest.ℓₘₐₓ), but the "
            * "right-hand side computes s=$(label[1]), ℓ ∈ $(label[2]):$(label[3])."
        ))
    end
end
@inline function Base.Broadcast.materialize!(dest::ModeWeights, bc)
    check_broadcast_destination(dest, bc)
    Base.Broadcast.materialize!(array_view(dest), bc)
    dest
end
@inline function Base.Broadcast.materialize!(
    dest::ModeWeights, bc::Base.Broadcast.Broadcasted{<:Any}
)
    check_broadcast_destination(dest, bc)
    Base.Broadcast.materialize!(array_view(dest), bc)
    dest
end

# The non-mutating `copy(bc)` builds its destination with the `similar` above and then fills
# it, so a `ModeWeights` destination needs a `copyto!` of its own; `.=` into an existing one
# goes through the two `materialize!` methods above, which check its labels first.
function Base.copyto!(w::ModeWeights, bc::Broadcast.Broadcasted)
    copyto!(w.data, bc)
    w
end

# `map` applies an arbitrary function, whose result there is no way to label, so it returns
# plain numbers; broadcasting is the way to keep the labels, where they still apply.
Base.map(f, w::ModeWeights) = map(f, w.data)

# The linear arithmetic of mode weights, which keeps the labels: sums and differences of
# weights whose labels agree, and products and quotients with numbers.
Base.:-(w::ModeWeights) = ModeWeights(-w.data, w.s, w.ℓₘᵢₙ, w.ℓₘₐₓ)
Base.:+(w::ModeWeights) = w
same_labels(a::ModeWeights, b::ModeWeights) = broadcast_label(a) == broadcast_label(b)
function check_same_labels(a::ModeWeights, b::ModeWeights, what)
    if !same_labels(a, b)
        throw(ArgumentError(
            "Cannot $what mode weights with different labels: s=$(a.s), ℓ ∈ $(a.ℓₘᵢₙ):$(a.ℓₘₐₓ) "
            * "and s=$(b.s), ℓ ∈ $(b.ℓₘᵢₙ):$(b.ℓₘₐₓ)."
        ))
    end
end
function Base.:+(a::ModeWeights, b::ModeWeights)
    check_same_labels(a, b, "add")
    ModeWeights(a.data + b.data, a.s, a.ℓₘᵢₙ, a.ℓₘₐₓ)
end
function Base.:-(a::ModeWeights, b::ModeWeights)
    check_same_labels(a, b, "subtract")
    ModeWeights(a.data - b.data, a.s, a.ℓₘᵢₙ, a.ℓₘₐₓ)
end
Base.:*(x::Number, w::ModeWeights) = ModeWeights(x * w.data, w.s, w.ℓₘᵢₙ, w.ℓₘₐₓ)
Base.:*(w::ModeWeights, x::Number) = ModeWeights(w.data * x, w.s, w.ℓₘᵢₙ, w.ℓₘₐₓ)
Base.:/(w::ModeWeights, x::Number) = ModeWeights(w.data / x, w.s, w.ℓₘᵢₙ, w.ℓₘₐₓ)
Base.:\(x::Number, w::ModeWeights) = ModeWeights(x \ w.data, w.s, w.ℓₘᵢₙ, w.ℓₘₐₓ)
Base.:*(J::LinearAlgebra.UniformScaling, w::ModeWeights) = J.λ * w
Base.:*(w::ModeWeights, J::LinearAlgebra.UniformScaling) = w * J.λ
Base.zero(w::ModeWeights) = ModeWeights(zero(w.data), w.s, w.ℓₘᵢₙ, w.ℓₘₐₓ)
# The in-place forms of the same arithmetic, which is what generic linear-algebra code
# reaches for.  Each returns the weights it wrote.  There is no `w + v` with a plain vector
# `v`: `+` is the sum of two functions, and `w .+ v` is the explicit form, which checks that
# `v` holds one value per mode.
function Base.fill!(w::ModeWeights, x)
    fill!(w.data, x)
    w
end
function LinearAlgebra.rmul!(w::ModeWeights, a::Number)
    LinearAlgebra.rmul!(w.data, a)
    w
end
function LinearAlgebra.lmul!(a::Number, w::ModeWeights)
    LinearAlgebra.lmul!(a, w.data)
    w
end
function Base.copyto!(w::ModeWeights, v::AbstractVector)
    if length(v) != length(w.data)
        throw(DimensionMismatch(
            "The vector has length $(length(v)), but the mode weights, with s=$(w.s) and "
            * "ℓ ∈ $(w.ℓₘᵢₙ):$(w.ℓₘₐₓ), are $(length(w.data)) weights."
        ))
    end
    copyto!(w.data, v)
    w
end
function Base.copyto!(w::ModeWeights, w′::ModeWeights)
    check_same_labels(w, w′, "copy")
    copyto!(w.data, w′.data)
    w
end
function Base.copy!(w::ModeWeights, w′::ModeWeights)
    check_same_labels(w, w′, "copy")
    copy!(w.data, w′.data)
    w
end
function LinearAlgebra.axpy!(a::Number, x::ModeWeights, y::ModeWeights)
    check_same_labels(x, y, "add")
    LinearAlgebra.axpy!(a, x.data, y.data)
    y
end
function LinearAlgebra.axpby!(a::Number, x::ModeWeights, b::Number, y::ModeWeights)
    check_same_labels(x, y, "add")
    LinearAlgebra.axpby!(a, x.data, b, y.data)
    y
end

# Comparison and reduction against plain storage.  Iteration gives `sum`, `maximum` and the
# rest for free; these are the ones that need the two representations to meet.
Base.:(==)(w::ModeWeights, v::AbstractVector) = w.data == v
Base.:(==)(v::AbstractVector, w::ModeWeights) = v == w.data
Base.isequal(w::ModeWeights, v::AbstractVector) = isequal(w.data, v)
Base.isequal(v::AbstractVector, w::ModeWeights) = isequal(v, w.data)
# Between two sets of mode weights the labels count too, as they do for `==`: the same
# numbers under different labels are the weights of different functions.
Base.isapprox(a::ModeWeights, b::ModeWeights; kwargs...) =
    same_labels(a, b) && isapprox(a.data, b.data; kwargs...)
Base.isapprox(a::ModeWeights, b::AbstractVector; kwargs...) = isapprox(a.data, b; kwargs...)
Base.isapprox(a::AbstractVector, b::ModeWeights; kwargs...) = isapprox(a, b.data; kwargs...)
Base.adjoint(w::ModeWeights) = adjoint(w.data)
Base.transpose(w::ModeWeights) = transpose(w.data)
LinearAlgebra.norm(w::ModeWeights, p::Real=2) = LinearAlgebra.norm(w.data, p)
# The inner product of the two functions, which is defined only between weights of one spin
# weight; the ranges of ℓ are required to agree as well, rather than summed over their overlap.
function LinearAlgebra.dot(a::ModeWeights, b::ModeWeights)
    check_same_labels(a, b, "take the inner product of")
    LinearAlgebra.dot(a.data, b.data)
end
LinearAlgebra.dot(a::ModeWeights, b::AbstractVector) = LinearAlgebra.dot(a.data, b)
LinearAlgebra.dot(a::AbstractVector, b::ModeWeights) = LinearAlgebra.dot(a, b.data)

# Natural indexing.
#
# Each of `w[ℓ, m]`, `w[ℓ, m] = v` and `w[ℓ, :]` has a method for each kind of container,
# whose indices are of the container's own type, `Int` or `HalfOddInteger`, and a method
# that accepts an index of any other type and converts it to the container's type with
# `container_index`, under the rules of the index methods, before it re-dispatches.  That
# method is what admits `w[3//2, 1//2]`, and what turns an index of the wrong kind, such as
# an integer applied to a half-integer `w`, or of the wrong type, such as an `Int32`, into
# an explanation rather than an error from deep inside `Yindex`.  It is written by hand
# rather than with `@index_methods`, because the kind of index it must accept is fixed by
# the container rather than by the indices themselves.

# The comparisons here are defined for either kind of index, and between the two kinds, so
# this is one method; a mixed call never reaches it, because the boundary method refuses it.
# A mode within the labels is then looked up in the storage under `@inbounds`, at the
# position the labels give it, so the storage is compared with that position as well: it is
# the caller's vector, not a copy, and may have been resized since the constructor compared
# its length with the labels.
@inline function check_mode(w::ModeWeights, ℓ, m)
    if !(w.ℓₘᵢₙ ≤ ℓ ≤ w.ℓₘₐₓ && -ℓ ≤ m ≤ ℓ)
        throw(BoundsError(w, (ℓ, m)))
    end
    check_storage(w, Yindex(ℓ, m, w.ℓₘᵢₙ))
end
@inline function check_storage(w::ModeWeights, i)
    if i > length(w.data)
        throw(storage_error(w, "an entry at position $i"))
    end
    nothing
end
# The operator kernels index the storage up to the length the labels imply, under
# `@inbounds`, so they compare the whole length with the labels, once per call.
@inline function check_storage_length(w::ModeWeights)
    n = Ysize(w.ℓₘᵢₙ, w.ℓₘₐₓ)
    if length(w.data) != n
        throw(storage_error(w, "length $n"))
    end
    nothing
end
@noinline function storage_error(w::ModeWeights, needed)
    DimensionMismatch(
        "The storage of these mode weights, with s=$(w.s) and ℓ ∈ $(w.ℓₘᵢₙ):$(w.ℓₘₐₓ), has "
        * "length $(length(w.data)), but the labels need $needed.  A `ModeWeights` uses its "
        * "vector as storage without copying it, so the vector must not be resized."
    )
end

# The natural indices of `w` as its own index type, each refused with the reason if it
# cannot be one.
@inline natural_indices(w::ModeWeights{T, IT}, ℓ, m) where {T, IT} =
    (container_index(IT, ℓ, w, "ℓ"), container_index(IT, m, w, "m"))
@inline natural_index(w::ModeWeights{T, IT}, ℓ) where {T, IT} = container_index(IT, ℓ, w, "ℓ")

"""
    w[ℓ, m]

The mode weight of ``(ℓ, m)`` in the [`ModeWeights`](@ref) `w`.  For a `w` with half-integer
indices, `ℓ` and `m` may be passed as `Rational{Int}`s with denominator 2 — `w[3//2, 1//2]` —
or as [`HalfOddInteger`](@ref)s, either way for each; for a `w` with integer indices they
must be `Int`s.  An index of another kind or type is refused with an `ArgumentError` that
says why, as the constructors refuse one.
"""
@propagate_inbounds function Base.getindex(w::ModeWeights{T, Int}, ℓ::Int, m::Int) where {T}
    @boundscheck check_mode(w, ℓ, m)
    @inbounds w.data[Yindex(ℓ, m, w.ℓₘᵢₙ)]
end
@propagate_inbounds function Base.getindex(
    w::ModeWeights{T, HalfOddInteger}, ℓ::HalfOddInteger, m::HalfOddInteger
) where {T}
    @boundscheck check_mode(w, ℓ, m)
    @inbounds w.data[Yindex(ℓ, m, w.ℓₘᵢₙ)]
end
@propagate_inbounds function Base.getindex(w::ModeWeights, ℓ::IndexType, m::IndexType)
    w[natural_indices(w, ℓ, m)...]
end
@propagate_inbounds function Base.setindex!(w::ModeWeights{T, Int}, v, ℓ::Int, m::Int) where {T}
    @boundscheck check_mode(w, ℓ, m)
    @inbounds w.data[Yindex(ℓ, m, w.ℓₘᵢₙ)] = v
end
@propagate_inbounds function Base.setindex!(
    w::ModeWeights{T, HalfOddInteger}, v, ℓ::HalfOddInteger, m::HalfOddInteger
) where {T}
    @boundscheck check_mode(w, ℓ, m)
    @inbounds w.data[Yindex(ℓ, m, w.ℓₘᵢₙ)] = v
end
@propagate_inbounds function Base.setindex!(w::ModeWeights, v, ℓ::IndexType, m::IndexType)
    ℓ′, m′ = natural_indices(w, ℓ, m)
    w[ℓ′, m′] = v
end

# As for `w[ℓ, m]`, an index of another type is converted by `natural_index`, or refused,
# before the method for the container's own type is reached.
"""
    w[ℓ, :]

A view of the mode weights of the [`ModeWeights`](@ref) `w` for the given ``ℓ``, indexed by
`m ∈ -ℓ:ℓ`.  This is a [`DegreeBlock`](@ref) for either kind of index; where the indices are
half-odd-integers they may be passed either as `Rational`s or as [`HalfOddInteger`](@ref)s.
Writing through the view writes into `w`.
"""
Base.getindex(w::ModeWeights{T, Int}, ℓ::Int, ::Colon) where {T} = degree_block(w, ℓ)
Base.getindex(w::ModeWeights{T, HalfOddInteger}, ℓ::HalfOddInteger, ::Colon) where {T} =
    degree_block(w, ℓ)
Base.getindex(w::ModeWeights, ℓ::IndexType, ::Colon) = w[natural_index(w, ℓ), :]
function degree_block(w::ModeWeights{T, IT}, ℓ::IT) where {T, IT}
    if !(w.ℓₘᵢₙ ≤ ℓ ≤ w.ℓₘₐₓ)
        throw(BoundsError(w, (ℓ, :)))
    end
    r = mode_range(w, ℓ)
    check_storage(w, last(r))
    DegreeBlock(view(w.data, r), ℓ)
end

# The labels, which also head the text/plain form, rather than the fields, since no
# constructor takes the fields in their order with the storage printed in full.
function Base.show(io::IO, w::ModeWeights{T}) where {T}
    print(io, "ModeWeights{$T} with s=$(w.s), ℓ ∈ $(w.ℓₘᵢₙ):$(w.ℓₘₐₓ)")
end
function Base.show(io::IO, ::MIME"text/plain", w::ModeWeights)
    show(io, w)
    println(io, ":")
    Base.print_array(io, w.data)
end



### Operators on mode weights
#
# One method covers all twelve: the operator is a value, so it says its own effect on the
# spin weight through `Δspin`, and the container already holds three indices of one kind,
# `Int` or `HalfOddInteger`.  The range of ℓ is unchanged even where the spin weight moves.
# The entries of the result below the new |s| belong to no harmonic; each is the product of
# an entry of the input with a coefficient that vanishes there, so it is zero where the
# input is finite, and is not dropped.  `ModeWeights(w; ℓₘᵢₙ, ℓₘₐₓ)` is what changes the
# range, and drops those entries.
function Base.:*(op::DifferentialOperator, w::ModeWeights{T}) where {T}
    # The result is allocated at the length of the input's storage, so this one check covers
    # both of the vectors that the kernel indexes.
    check_storage_length(w)
    Treal = real(float(T))
    out = similar(w.data, Base.promote_op(*, coefftype(op, Treal), T))
    apply_operator!(out, op, bandstructure(op), w.data, w.s, w.ℓₘᵢₙ, w.ℓₘₐₓ, Treal)
    ModeWeights(out, w.s + Δspin(op), w.ℓₘᵢₙ, w.ℓₘₐₓ)
end
(op::DifferentialOperator)(w::ModeWeights) = op * w

# The in-place form, for a loop over many sets of weights.  Aliasing is refused for the
# banded operators, whose kernels read a neighbor that an in-place write may already have
# clobbered; it would be safe for the diagonal ones, but allowing it there only would be a
# trap.
function LinearAlgebra.mul!(
    w′::ModeWeights, op::DifferentialOperator, w::ModeWeights{T}
) where {T}
    if spin(w′) != w.s + Δspin(op) || ℓₘᵢₙ(w′) != w.ℓₘᵢₙ || ℓₘₐₓ(w′) != w.ℓₘₐₓ
        throw(operator_output_error(w′, op, w))
    end
    check_storage_length(w)
    check_storage_length(w′)
    if Base.mightalias(w′.data, w.data)
        throw(ArgumentError(
            "The output aliases the input.  $(nameof(op)) reads neighboring modes, so it "
            * "cannot be applied in place; pass a separate destination, such as `similar(w)`."
        ))
    end
    apply_operator!(
        w′.data, op, bandstructure(op), w.data, w.s, w.ℓₘᵢₙ, w.ℓₘₐₓ, real(float(T))
    )
    w′
end
# The operators keep the range of ℓ of their input, so a destination allocated with the
# default range of its own spin weight, `abs(s′):ℓₘₐₓ`, is refused whenever |s′| ≠ |s|; the
# message says how to allocate one, and how to change the range of the result afterwards.
@noinline function operator_output_error(w′::ModeWeights, op, w::ModeWeights{T}) where {T}
    s′ = w.s + Δspin(op)
    ArgumentError(
        "The output has s=$(spin(w′)) and ℓ ∈ $(ℓₘᵢₙ(w′)):$(ℓₘₐₓ(w′)), but $(nameof(op)) "
        * "applied to these weights gives s=$s′ and ℓ ∈ $(w.ℓₘᵢₙ):$(w.ℓₘₐₓ), since the "
        * "operators keep the range of ℓ of their input.  Allocate the output with "
        * "`ModeWeights{$T}(undef, $s′, $(w.ℓₘᵢₙ), $(w.ℓₘₐₓ))`, and use "
        * "`ModeWeights(w′; ℓₘᵢₙ, ℓₘₐₓ)` to copy the result into another range of ℓ."
    )
end

# Bare storage, at least as long as the result, is accepted as the output too, and the
# result comes back labelled, as a `ModeWeights` over it (see `mode_weights_view`).
function LinearAlgebra.mul!(w′::AbstractVector, op::DifferentialOperator, w::ModeWeights)
    mul!(mode_weights_view(w′, w.s + Δspin(op), w.ℓₘᵢₙ, w.ℓₘₐₓ), op, w)
end

# The in-place operations that write mode weights accept, as their output, a bare vector at
# least as long as the result, and return the result as a `ModeWeights` over its first
# entries — always a view, even when the length is exact, so that the storage is shared
# rather than copied and the type returned does not depend on the length.
function mode_weights_view(v::AbstractVector, s, ℓₘᵢₙ, ℓₘₐₓ)
    Base.require_one_based_indexing(v)
    n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
    if length(v) < n
        throw(DimensionMismatch(
            "The output has length $(length(v)); at least Ysize($ℓₘᵢₙ, $ℓₘₐₓ) = $n is needed."
        ))
    end
    ModeWeights(view(v, 1:n), s, ℓₘᵢₙ, ℓₘₐₓ)
end
