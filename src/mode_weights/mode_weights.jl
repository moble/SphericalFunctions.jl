"""
    ModeWeights(data, s=0; ℓₘᵢₙ=abs(s))
    ModeWeights(data, s, ℓₘᵢₙ, ℓₘₐₓ)
    ModeWeights{T}(undef, s, ℓₘᵢₙ, ℓₘₐₓ)
    ModeWeights{T}(undef, s, ℓₘₐₓ)

Vector of mode weights ``f_{ℓ,m}`` of a spin-weighted function ``f = \\sum_{ℓ,m} f_{ℓ,m}\\,
{}_sY_{ℓ,m}``, stored in the canonical ordering `[f(ℓ, m) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ]`
(see [`Yindex`](@ref)), together with the spin weight `s` and the range of ``ℓ``.

A `ModeWeights` is an [`AbstractModeContainer`](@ref
SphericalFunctions.AbstractModeContainer), not an `AbstractVector`; [`array_view`](@ref) gives
the flat 1-based storage, which is what the transforms and the operator matrices take.  Linear
indexing and broadcasting still work, and a shape-preserving broadcast keeps the wrapper.  In
addition
- `w[ℓ, m]` reads or writes the weight of mode ``(ℓ, m)``,
- `w[ℓ, :]` is a [`DegreeBlock`](@ref) view of the weights for one ``ℓ``, indexed by
  `m ∈ -ℓ:ℓ`,
- `modes(w)` is the vector of `(ℓ, m)` pairs in storage order,
- `spin(w)`, `ℓₘᵢₙ(w)`, `ℓₘₐₓ(w)` are the parameters, and `parent(w)` is the storage,
- the differential operators [`L²`](@ref), [`Lz`](@ref), [`L₊`](@ref), [`L₋`](@ref),
  [`Lx`](@ref), [`Ly`](@ref), [`R²`](@ref), [`Rz`](@ref), [`R₊`](@ref), [`R₋`](@ref),
  [`ð`](@ref), [`ð̄`](@ref) give a new `ModeWeights` with the spin weight adjusted where
  appropriate, written either `ð * w` or `ð(w)`.  These build no matrix: the operator is
  applied by a loop, so the only allocation is the result, and `mul!(w′, ð, w)` into a
  correctly labelled destination allocates nothing at all,
- multiplying by an operator *matrix* instead — `ð(s, ℓₘᵢₙ, ℓₘₐₓ) * w` — gives a plain
  `Vector`, because a matrix of numbers cannot say what spin weight its result has; use the
  operators themselves when you want the answer labelled, and
- `w(R)` evaluates the function at the rotor `R` (see [`sYlm`](@ref)).

When constructed from `data` alone, `ℓₘₐₓ` is deduced from `length(data)` and `ℓₘᵢₙ`, which
must match exactly.  The `data` vector is used as storage, not copied.  The `undef` forms
allocate uninitialized storage of type `T` instead, and in the three-argument form `ℓₘᵢₙ`
defaults to `abs(s)`, as it does when `data` is given.

# Half-integer indices

The spin weight and the range of ``ℓ`` may be half-integers, passed as `Rational`s with
denominator 2 — as in `ModeWeights(data, 1//2)` or `ModeWeights{T}(undef, 1//2, 1//2, 7//2)`
— in which case every ``ℓ`` and ``m`` of the ordering is a half-odd-integer, and `ℓₘᵢₙ` may be
as small as `1//2`.  The indices in one call must all be of one kind, integers or
half-odd-integers; a call that mixes them, such as `ModeWeights(data, 1//2, 0, 7//2)`, is an
error.  The parameters are stored as [`HalfOddInteger`](@ref)s, which is also what `modes(w)`
and the axis of `w[ℓ, :]` are made of; `w[ℓ, m]` accepts either type.  For such a `w`,
`w[ℓ, :]` is a [`DegreeBlock`](@ref), indexed by `m ∈ -ℓ:ℓ`, exactly as it is for integer
indices.
"""
struct ModeWeights{T, IT<:IntegerHalf, V<:AbstractVector{T}} <: AbstractModeContainer{T, IT}
    data::V
    s::IT
    ℓₘᵢₙ::IT
    ℓₘₐₓ::IT
    # These checks hold for either kind of index: the floor of ℓₘᵢₙ is 0 for integers and 1/2
    # for half-odd-integers, and `ℓₘᵢₙ < 0` is the right test for both, since no
    # half-odd-integer lies between 0 and 1/2.
    function ModeWeights(data::V, s::IT, ℓₘᵢₙ::IT, ℓₘₐₓ::IT) where {T, IT<:IntegerHalf, V<:AbstractVector{T}}
        Base.require_one_based_indexing(data)
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

# The outer constructors are boundary methods: each accepts an index of any permitted type —
# `Integer`, `HalfOddInteger`, or a `Rational` with denominator 2 — and normalizes them with
# `unify_indices`, which is also what refuses a mixture of the two kinds of index with an
# explanation.  The inner constructor above is the only one reached with three indices of one
# type, and is the one that validates them.  `ℓₘᵢₙ` is a keyword argument, and keyword
# arguments take no part in dispatch, so the normalization has to happen here rather than in a
# second method with `ℓₘᵢₙ::IT` in its signature, which would refuse `ℓₘᵢₙ=1//2` with a bare
# `TypeError` before anything could convert it.  Where `ℓₘᵢₙ` defaults to `abs(s)`, the
# spin weight is normalized before `abs` is taken, so that not even that is applied to a
# `Rational`.
function ModeWeights(
    data::AbstractVector, s::IndexArgument=0; ℓₘᵢₙ::IndexArgument=abs(half_integer(s))
)
    deduced_mode_weights(data, unify_indices(s, ℓₘᵢₙ)...)
end
function ModeWeights(
    data::AbstractVector, s::IndexArgument, ℓₘᵢₙ::IndexArgument, ℓₘₐₓ::IndexArgument
)
    ModeWeights(data, unify_indices(s, ℓₘᵢₙ, ℓₘₐₓ)...)
end
function ModeWeights{T}(
    ::UndefInitializer, s::IndexArgument, ℓₘᵢₙ::IndexArgument, ℓₘₐₓ::IndexArgument
) where {T}
    s, ℓₘᵢₙ, ℓₘₐₓ = unify_indices(s, ℓₘᵢₙ, ℓₘₐₓ)
    ModeWeights(Vector{T}(undef, Ysize(ℓₘᵢₙ, ℓₘₐₓ)), s, ℓₘᵢₙ, ℓₘₐₓ)
end
function ModeWeights{T}(::UndefInitializer, s::IndexArgument, ℓₘₐₓ::IndexArgument) where {T}
    s, ℓₘₐₓ = unify_indices(s, ℓₘₐₓ)
    ModeWeights{T}(undef, s, abs(s), ℓₘₐₓ)
end

# Deduce ℓₘₐₓ from the length of the data, given `s` and `ℓₘᵢₙ` of one kind.  The kind of the
# indices selects the method, since the relation between the length and ℓₘₐₓ is written
# differently for the two.
function deduced_mode_weights(data::AbstractVector, s::IT, ℓₘᵢₙ::IT) where {IT<:Integer}
    # Deduce ℓₘₐₓ from (ℓₘₐₓ+1)² = length + ℓₘᵢₙ²
    N = length(data) + ℓₘᵢₙ^2
    ℓₘₐₓ = isqrt(N) - 1
    if (ℓₘₐₓ + 1)^2 != N
        throw(ArgumentError(
            "The data has length $(length(data)), which is not Ysize(ℓₘᵢₙ=$ℓₘᵢₙ, ℓₘₐₓ) "
            * "for any ℓₘₐₓ."
        ))
    end
    # `ℓₘₐₓ` is an `Int` whatever the concrete type of `ℓₘᵢₙ`, because `length` is, so the
    # three indices are promoted to one integer type, as they always have been.
    ModeWeights(data, promote(s, ℓₘᵢₙ, ℓₘₐₓ)...)
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

Base.parent(w::ModeWeights) = w.data

"""
    spin(w)

The spin weight of a [`ModeWeights`](@ref) vector, of an [`SSHT`](@ref) transform, or of an
[`sYlmCalculator`](@ref) built for a single one.  A calculator built for a range of spin
weights has no single value to report, so it has no method here; ask it for
[`spins`](@ref SphericalFunctions.spins) instead, which answers for either kind.

A function of spin weight ``s`` has ``R_z f = s f``, and is expanded in the harmonics
``{}_{s}Y_{ℓ,m}`` with ``ℓ ≥ |s|``.  The spin weight is kept alongside the numbers
because nothing about the numbers themselves reveals it.

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
ℓₘᵢₙ(w::ModeWeights) = w.ℓₘᵢₙ
ℓₘₐₓ(w::ModeWeights) = w.ℓₘₐₓ

"""
    modes(w::ModeWeights)

The `(ℓ, m)` pairs of `w`, in storage order (see [`Yrange`](@ref)).
"""
modes(w::ModeWeights) = Yrange(w.ℓₘᵢₙ, w.ℓₘₐₓ)

# The array-like interface, written out rather than inherited.  A `ModeWeights` was an
# `AbstractVector` before version 3, which made `op * w` and `w .+ 1` work through the generic
# machinery; it is now an [`AbstractModeContainer`](@ref) like the rest, and these are the
# methods that keep the useful part of the old behavior.  Losing the subtyping costs less than
# it appears to: the transforms in `ssht/` never used it, reaching for the raw storage before
# every `mul!` and `ldiv!` (what is now called [`array_view`](@ref)).
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
@propagate_inbounds Base.getindex(w::ModeWeights, i::Int) = w.data[i]
@propagate_inbounds Base.getindex(w::ModeWeights, r::AbstractRange{Int}) = w.data[r]
@propagate_inbounds Base.setindex!(w::ModeWeights, v, i::Int) = (w.data[i] = v)
Base.similar(w::ModeWeights) = ModeWeights(similar(w.data), w.s, w.ℓₘᵢₙ, w.ℓₘₐₓ)
Base.similar(w::ModeWeights, ::Type{S}) where {S} =
    ModeWeights(similar(w.data, S), w.s, w.ℓₘᵢₙ, w.ℓₘₐₓ)
Base.collect(w::ModeWeights) = collect(w.data)
Base.Array(w::ModeWeights) = collect(w.data)
Base.Vector(w::ModeWeights) = collect(w.data)
function Base.:(==)(a::ModeWeights, b::ModeWeights)
    a.s == b.s && a.ℓₘᵢₙ == b.ℓₘᵢₙ && a.ℓₘₐₓ == b.ℓₘₐₓ && a.data == b.data
end

# Multiplying by an operator *matrix* returns plain storage, not a `ModeWeights`.  A matrix
# cannot say what spin weight its result has: `ð(s, ℓₘᵢₙ, ℓₘₐₓ)` raises the spin weight by
# one, `L₊(s, ℓₘᵢₙ, ℓₘₐₓ)` leaves it alone, and the two are both `Diagonal`/`Bidiagonal`
# matrices of numbers with nothing to tell them apart.  Without these methods the generic
# `AbstractVector` machinery would hand the result `w`'s own spin weight, which for the
# spin-changing operators is silently wrong; `ð(w)` is the expression that keeps the label
# right.  These three cover every matrix type the operators in this package return.
Base.:*(A::AbstractMatrix, w::ModeWeights) = A * parent(w)
# ... and on the other side, which is the outer product `w * w'`.
Base.:*(w::ModeWeights, A::AbstractMatrix) = parent(w) * A
Base.:*(A::Diagonal, w::ModeWeights) = A * parent(w)
Base.:*(A::Bidiagonal, w::ModeWeights) = A * parent(w)
Base.:*(A::Tridiagonal, w::ModeWeights) = A * parent(w)
Base.copy(w::ModeWeights) = ModeWeights(copy(w.data), w.s, w.ℓₘᵢₙ, w.ℓₘₐₓ)

# Broadcasting preserves the `ModeWeights` wrapper when the shape is unchanged.  An
# elementwise operation cannot change which modes are held, so the spin weight and the range of
# ℓ that label them survive it; only a broadcast that changes the length drops back to a plain
# `Vector`.  This used a `Broadcast.ArrayStyle` before version 3, which is available only to an
# `AbstractArray`; a style of this type's own does the same job, given a `broadcastable` that
# hands back the container rather than `collect`ing it, and the `axes` and linear `getindex`
# defined above.  Writing *into* one with `.=` is handled with the other containers, in
# `array_view.jl`.
struct ModeWeightsStyle <: Broadcast.AbstractArrayStyle{1} end
ModeWeightsStyle(::Val{0}) = ModeWeightsStyle()
ModeWeightsStyle(::Val{1}) = ModeWeightsStyle()
ModeWeightsStyle(::Val{N}) where {N} = Broadcast.DefaultArrayStyle{N}()
Base.BroadcastStyle(::Type{<:ModeWeights}) = ModeWeightsStyle()
Base.Broadcast.broadcastable(w::ModeWeights) = w
function Base.similar(bc::Broadcast.Broadcasted{ModeWeightsStyle}, ::Type{S}) where {S}
    w = find_modeweights(bc)
    if axes(bc) == axes(w)
        ModeWeights(similar(Vector{S}, axes(bc)), w.s, w.ℓₘᵢₙ, w.ℓₘₐₓ)
    else
        similar(Vector{S}, axes(bc))
    end
end
find_modeweights(bc::Broadcast.Broadcasted) = find_modeweights(bc.args)
find_modeweights(args::Tuple) = find_modeweights(find_modeweights(args[1]), Base.tail(args))
find_modeweights(x) = x
find_modeweights(::Tuple{}) = nothing
find_modeweights(w::ModeWeights, rest) = w
find_modeweights(::Any, rest) = find_modeweights(rest)

# The non-mutating `copy(bc)` builds its destination with the `similar` above and then fills
# it, so a `ModeWeights` destination needs a `copyto!` of its own; `.=` into an existing one
# goes through `materialize!` in `array_view.jl` instead.
function Base.copyto!(w::ModeWeights, bc::Broadcast.Broadcasted)
    copyto!(w.data, bc)
    w
end

# `map` keeps the wrapper for the same reason broadcasting does.
Base.map(f, w::ModeWeights) = ModeWeights(map(f, w.data), w.s, w.ℓₘᵢₙ, w.ℓₘₐₓ)

# Arithmetic that cannot change which modes are held keeps the label; `similar` keeps it for
# the same length and falls back to a plain array for any other shape.  These were inherited
# from `AbstractVector` before version 3 and are written out now.
Base.:-(w::ModeWeights) = ModeWeights(-w.data, w.s, w.ℓₘᵢₙ, w.ℓₘₐₓ)
Base.:+(w::ModeWeights) = w
Base.similar(w::ModeWeights, n::Integer) = similar(w, eltype(w), n)
Base.similar(w::ModeWeights, ::Type{S}, n::Integer) where {S} =
    n == length(w) ? ModeWeights(similar(w.data, S), w.s, w.ℓₘᵢₙ, w.ℓₘₐₓ) : similar(w.data, S, n)
Base.similar(w::ModeWeights, dims::Dims) = similar(w, eltype(w), dims)
Base.similar(w::ModeWeights, ::Type{S}, dims::Dims) where {S} =
    dims == size(w) ? ModeWeights(similar(w.data, S), w.s, w.ℓₘᵢₙ, w.ℓₘₐₓ) : similar(w.data, S, dims)

# Comparison and reduction against plain storage.  Iteration gives `sum`, `maximum` and the
# rest for free; these are the ones that need the two representations to meet.
Base.:(==)(w::ModeWeights, v::AbstractVector) = w.data == v
Base.:(==)(v::AbstractVector, w::ModeWeights) = v == w.data
Base.isequal(w::ModeWeights, v::AbstractVector) = isequal(w.data, v)
Base.isequal(v::AbstractVector, w::ModeWeights) = isequal(v, w.data)
Base.isapprox(a::ModeWeights, b::ModeWeights; kwargs...) = isapprox(a.data, b.data; kwargs...)
Base.isapprox(a::ModeWeights, b::AbstractVector; kwargs...) = isapprox(a.data, b; kwargs...)
Base.isapprox(a::AbstractVector, b::ModeWeights; kwargs...) = isapprox(a, b.data; kwargs...)
Base.adjoint(w::ModeWeights) = adjoint(w.data)
Base.transpose(w::ModeWeights) = transpose(w.data)
LinearAlgebra.norm(w::ModeWeights, p::Real=2) = LinearAlgebra.norm(w.data, p)
LinearAlgebra.dot(a::ModeWeights, b::ModeWeights) = LinearAlgebra.dot(a.data, b.data)
LinearAlgebra.dot(a::ModeWeights, b::AbstractVector) = LinearAlgebra.dot(a.data, b)
LinearAlgebra.dot(a::AbstractVector, b::ModeWeights) = LinearAlgebra.dot(a, b.data)

# Natural indexing.
#
# Each of `w[ℓ, m]`, `w[ℓ, m] = v` and `w[ℓ, :]` has a method for each kind of index — with
# the indices of the same kind as `w`'s own — and a boundary method that accepts any other
# type, normalizes it, checks that it is of `w`'s kind, and re-dispatches.  The boundary
# is what admits `w[3//2, 1//2]`, and what turns an integer index applied to a half-integer
# `w` into an explanation rather than a `MethodError` deep inside `Yindex`.

# The comparisons here are defined for either kind of index, and between the two kinds, so
# this is one method; a mixed call never reaches it, because the boundary method refuses it.
@inline function check_mode(w::ModeWeights, ℓ, m)
    if !(w.ℓₘᵢₙ ≤ ℓ ≤ w.ℓₘₐₓ && -ℓ ≤ m ≤ ℓ)
        throw(BoundsError(w, (ℓ, m)))
    end
end

# Normalize the natural indices of `w` and require them to be of `w`'s kind.  `half_integers`
# has already refused a mixture of the two kinds among the indices themselves, so only the
# first need be compared with `w`.
function mode_kind_error(::Type{IT}, indices) where {IT<:IntegerHalf}
    kind, example = IT <: Integer ? ("integers", "3") : ("half-odd-integers", "7//2")
    ArgumentError(
        "The indices of this `ModeWeights` are $kind, like $example, so the indices used with "
        * "it must be too; got " * join(indices, ", ") * "."
    )
end
@inline function natural_indices(::ModeWeights{T, IT}, ℓ, m) where {T, IT}
    ℓ′, m′ = half_integers(ℓ, m)
    isindex(IT, ℓ′) || throw(mode_kind_error(IT, (ℓ, m)))
    ℓ′, m′
end
@inline function natural_index(::ModeWeights{T, IT}, ℓ) where {T, IT}
    ℓ′ = half_integer(ℓ)
    isindex(IT, ℓ′) || throw(mode_kind_error(IT, (ℓ,)))
    ℓ′
end

"""
    w[ℓ, m]

The mode weight of ``(ℓ, m)`` in the [`ModeWeights`](@ref) `w`.  For a `w` with half-integer
indices, `ℓ` and `m` may be passed as `Rational`s — `w[3//2, 1//2]` — or as
[`HalfOddInteger`](@ref)s; for a `w` with integer indices they must be integers.
"""
@propagate_inbounds function Base.getindex(w::ModeWeights{T, <:Integer}, ℓ::Integer, m::Integer) where {T}
    @boundscheck check_mode(w, ℓ, m)
    @inbounds w.data[Yindex(ℓ, m, w.ℓₘᵢₙ)]
end
@propagate_inbounds function Base.getindex(
    w::ModeWeights{T, HalfOddInteger}, ℓ::HalfOddInteger, m::HalfOddInteger
) where {T}
    @boundscheck check_mode(w, ℓ, m)
    @inbounds w.data[Yindex(ℓ, m, w.ℓₘᵢₙ)]
end
@propagate_inbounds function Base.getindex(w::ModeWeights, ℓ::IndexArgument, m::IndexArgument)
    w[natural_indices(w, ℓ, m)...]
end
@propagate_inbounds function Base.setindex!(w::ModeWeights{T, <:Integer}, v, ℓ::Integer, m::Integer) where {T}
    @boundscheck check_mode(w, ℓ, m)
    @inbounds w.data[Yindex(ℓ, m, w.ℓₘᵢₙ)] = v
end
@propagate_inbounds function Base.setindex!(
    w::ModeWeights{T, HalfOddInteger}, v, ℓ::HalfOddInteger, m::HalfOddInteger
) where {T}
    @boundscheck check_mode(w, ℓ, m)
    @inbounds w.data[Yindex(ℓ, m, w.ℓₘᵢₙ)] = v
end
@propagate_inbounds function Base.setindex!(w::ModeWeights, v, ℓ::IndexArgument, m::IndexArgument)
    ℓ′, m′ = natural_indices(w, ℓ, m)
    w[ℓ′, m′] = v
end

"""
    w[ℓ, :]

A view of the mode weights of the [`ModeWeights`](@ref) `w` for the given ``ℓ``, indexed by
`m ∈ -ℓ:ℓ`.  This is a [`DegreeBlock`](@ref) for either kind of index; where the indices are
half-odd-integers they may be passed either as `Rational`s or as [`HalfOddInteger`](@ref)s.
Writing through the view writes into `w`.
"""
function Base.getindex(w::ModeWeights{T, IT}, ℓ::IT, ::Colon) where {T, IT<:IntegerHalf}
    if !(w.ℓₘᵢₙ ≤ ℓ ≤ w.ℓₘₐₓ)
        throw(BoundsError(w, (ℓ, :)))
    end
    i₀ = Yindex(ℓ, -ℓ, w.ℓₘᵢₙ)
    DegreeBlock(view(w.data, i₀:i₀+2ℓ), ℓ)
end
Base.getindex(w::ModeWeights, ℓ::IndexArgument, ::Colon) = w[natural_index(w, ℓ), :]

function Base.show(io::IO, ::MIME"text/plain", w::ModeWeights{T}) where {T}
    println(io, "ModeWeights{$T} with s=$(w.s), ℓ ∈ $(w.ℓₘᵢₙ):$(w.ℓₘₐₓ):")
    Base.print_array(io, w.data)
end



### Operators on mode weights
#
# One method covers all twelve: the operator is a value, so it says its own effect on the spin
# weight through `Δspin`, and the container already holds three normalized indices of one kind,
# so `op(...)` reaches the worker directly rather than going through the `IndexArgument`
# boundary again.  The range of ℓ is unchanged even where the spin weight moves — entries that
# fall outside the new |s| are zeroed by the coefficients, not dropped.
function Base.:*(op::DifferentialOperator, w::ModeWeights{T}) where {T}
    Treal = real(float(T))
    out = similar(w.data, Base.promote_op(*, coefftype(op, Treal), T))
    apply_operator!(out, op, bandstructure(op), w.data, w.s, w.ℓₘᵢₙ, w.ℓₘₐₓ, Treal)
    ModeWeights(out, w.s + Δspin(op), w.ℓₘᵢₙ, w.ℓₘₐₓ)
end
(op::DifferentialOperator)(w::ModeWeights) = op * w

# The in-place form, for a loop over many sets of weights.  Aliasing is refused for the banded
# operators, whose kernels read a neighbor that an in-place write may already have clobbered;
# it would be safe for the diagonal ones, but allowing it there only would be a trap.
function LinearAlgebra.mul!(
    w′::ModeWeights, op::DifferentialOperator, w::ModeWeights{T}
) where {T}
    if spin(w′) != w.s + Δspin(op) || ℓₘᵢₙ(w′) != w.ℓₘᵢₙ || ℓₘₐₓ(w′) != w.ℓₘₐₓ
        error(
            "The output has s=$(spin(w′)) and ℓ ∈ $(ℓₘᵢₙ(w′)):$(ℓₘₐₓ(w′)), but $(nameof(op)) "
            * "applied to these weights gives s=$(w.s + Δspin(op)) and "
            * "ℓ ∈ $(w.ℓₘᵢₙ):$(w.ℓₘₐₓ)."
        )
    end
    if Base.mightalias(w′.data, w.data)
        error(
            "The output aliases the input.  $(nameof(op)) reads neighboring modes, so it "
            * "cannot be applied in place; pass a separate destination, such as `similar(w)`."
        )
    end
    apply_operator!(
        w′.data, op, bandstructure(op), w.data, w.s, w.ℓₘᵢₙ, w.ℓₘₐₓ, real(float(T))
    )
    w′
end
