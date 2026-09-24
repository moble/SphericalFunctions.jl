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
indexing works.  So does the arithmetic of mode weights as the weights of functions, which
keeps the labels: `a + b` and `a - b` (and their broadcast forms) when the labels agree,
`complex.(a, b)` from real and imaginary parts whose labels agree, and products and quotients
with numbers, or elementwise with a plain vector of factors (a diagonal operator, such as a
filter).  Anything else — adding a number to every weight, the elementwise product of two sets
of weights, `conj.(w)`, `abs2.(w)` — would label numbers that are not the mode weights of any
function, and is an error; arithmetic on the raw numbers goes through `array_view(w)`.  `≈` and
`dot` between two `ModeWeights` require their labels to agree, and `map` returns plain
numbers.  In addition
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
  correctly labelled destination allocates nothing at all (`w′` may also be a bare vector at
  least as long as the result, which then comes back as a `ModeWeights` over it),
- multiplying by an operator *matrix* instead — `ð(s, ℓₘᵢₙ, ℓₘₐₓ) * w` — gives a plain
  `Vector`, because a matrix of numbers cannot say what spin weight its result has; use the
  operators themselves when you want the answer labelled, and
- `w(R)` evaluates the function at the rotor `R` (see [`sYlm`](@ref)).

When constructed from `data` alone, `ℓₘₐₓ` is deduced from `length(data)` and `ℓₘᵢₙ`, which
must match exactly.  The `data` vector is used as storage, not copied, so it must keep its
length for as long as the `ModeWeights` is in use; the operators, `w[ℓ, :]`, and `w[ℓ, m]`
outside `@inbounds` throw a `DimensionMismatch` when they find that it has been resized.  The
`undef` forms allocate uninitialized storage of type `T` instead, and in the three-argument
form `ℓₘᵢₙ` defaults to `abs(s)`, as it does when `data` is given.

The weights are real or complex numbers, and a quaternion element type is refused.  A `Rotor`
is a quaternion, so `R * w` is an error; the rotation of the function by `R` is
`D(R, ℓₘₐₓ(w)) * w` (see [`D`](@ref)).

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

# The array-like interface, written out rather than inherited.  A `ModeWeights` is an
# [`AbstractModeContainer`](@ref) like the rest, not an `AbstractVector`, so `op * w` and
# `w .+ 1` do not come for free from the generic machinery; these are the methods that supply
# the useful part of that behavior.  Forgoing the subtyping costs less than it appears to: the
# transforms in `ssht/` reach for the raw storage before every `mul!` and `ldiv!` anyway,
# through [`array_view`](@ref).
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

# Broadcasting over mode weights is allowed only where the result is again the mode weights of
# a function, with labels that can be stated: sums and differences of weights with the same
# spin weight and range of ℓ, and products and quotients of weights with numbers — or with
# plain vectors, one factor per mode, which is a diagonal operator such as a filter.  The
# result then carries those labels.  Anything else that involves mode weights is refused,
# because it would put a label on numbers it does not describe: the mode weights of a product
# of two functions are not the product of their mode weights, those of the complex conjugate
# are not the conjugates, and `abs2.(w)` is not the mode weights of anything.  Such arithmetic
# on the raw numbers is still available through `array_view(w)`.  A `Broadcast.ArrayStyle`
# would be the usual way to get a wrapped result, but it is available only to an
# `AbstractArray`; a style of this type's own does the same job, given a `broadcastable` that
# hands back the container rather than `collect`ing it, and the `axes` and linear `getindex`
# defined above.  Where an array of another style takes part, such as a `StaticArray`, the two
# styles conflict and the result is a plain array, which has no labels to check.
struct ModeWeightsStyle <: Broadcast.AbstractArrayStyle{1} end
ModeWeightsStyle(::Val{0}) = ModeWeightsStyle()
ModeWeightsStyle(::Val{1}) = ModeWeightsStyle()
ModeWeightsStyle(::Val{N}) where {N} = Broadcast.DefaultArrayStyle{N}()
Base.BroadcastStyle(::Type{<:ModeWeights}) = ModeWeightsStyle()
Base.Broadcast.broadcastable(w::ModeWeights) = w
function Base.similar(bc::Broadcast.Broadcasted{ModeWeightsStyle}, ::Type{S}) where {S}
    label = broadcast_label(bc)
    if label !== nothing && length(axes(bc)) == 1 && length(bc) == Ysize(label[2], label[3])
        ModeWeights(similar(Vector{S}, axes(bc)), label...)
    else
        similar(Vector{S}, axes(bc))
    end
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
# function — sums, differences, and complex weights built from their real and imaginary parts.
# The mode weights among the arguments must have the same labels, and every other argument of
# at most one dimension must hold one value per mode.  That refuses a number, however it is
# written — a literal, a `Ref`, a zero-dimensional array, or a zero-dimensional broadcast such
# as the `a .* b` of `w .+ a .* b` — and an array or tuple too short to hold one value per mode,
# such as `[1.0]` or `(1,)`, which broadcasting would extend to every mode just as it does a
# number.  An argument of two or more dimensions makes the result a matrix, such as the outer
# sum `w .+ transpose(w)`, which is never labelled, so it is left alone.  The messages are
# built only on the branches that throw them, so that a broadcast that passes allocates
# nothing here.
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

# Writing into mode weights with `.=` checks the labels in the same way: whatever the right-hand
# side computes must be what the destination's labels say it holds.  (The other containers go
# through the methods in `array_view.jl`.)
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
# goes through `materialize!` in `array_view.jl` instead.
function Base.copyto!(w::ModeWeights, bc::Broadcast.Broadcasted)
    copyto!(w.data, bc)
    w
end

# `map` applies an arbitrary function, whose result there is no way to label, so it returns
# plain numbers; broadcasting is the way to keep the labels, where they still apply.
Base.map(f, w::ModeWeights) = map(f, w.data)

# The linear arithmetic of mode weights, which keeps the labels: sums and differences of
# weights whose labels agree, and products and quotients with numbers.  `similar` keeps the
# labels for the same length and falls back to a plain array for any other shape.
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
# Between two sets of mode weights the labels count too, as they do for `==`: the same numbers
# under different labels are the weights of different functions.
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
# Each of `w[ℓ, m]`, `w[ℓ, m] = v` and `w[ℓ, :]` has a method for each kind of index — with
# the indices of the same kind as `w`'s own — and a boundary method that accepts any other
# type, normalizes it, checks that it is of `w`'s kind, and re-dispatches.  The boundary
# is what admits `w[3//2, 1//2]`, and what turns an integer index applied to a half-integer
# `w` into an explanation rather than a `MethodError` deep inside `Yindex`.

# The comparisons here are defined for either kind of index, and between the two kinds, so
# this is one method; a mixed call never reaches it, because the boundary method refuses it.
# A mode within the labels is then looked up in the storage under `@inbounds`, at the position
# the labels give it, so the storage is compared with that position as well: it is the
# caller's vector, not a copy, and may have been resized since the constructor compared its
# length with the labels.
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
# The operator kernels index the storage up to the length the labels imply, under `@inbounds`,
# so they compare the whole length with the labels, once per call.
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

# As for `w[ℓ, m]`, a container with integer indices accepts an `Integer` of any type.  The
# worker cannot require the container's own type exactly: `natural_index` leaves an integer's
# type alone, so the boundary method below would then call itself forever for, say, an
# `Int32` ℓ.  The block is labelled in the container's own type, after the bounds check, so
# that an out-of-range index is a `BoundsError` rather than an `InexactError`.
"""
    w[ℓ, :]

A view of the mode weights of the [`ModeWeights`](@ref) `w` for the given ``ℓ``, indexed by
`m ∈ -ℓ:ℓ`.  This is a [`DegreeBlock`](@ref) for either kind of index; where the indices are
half-odd-integers they may be passed either as `Rational`s or as [`HalfOddInteger`](@ref)s.
Writing through the view writes into `w`.
"""
Base.getindex(w::ModeWeights{T, <:Integer}, ℓ::Integer, ::Colon) where {T} = degree_block(w, ℓ)
Base.getindex(w::ModeWeights{T, HalfOddInteger}, ℓ::HalfOddInteger, ::Colon) where {T} =
    degree_block(w, ℓ)
Base.getindex(w::ModeWeights, ℓ::IndexArgument, ::Colon) = w[natural_index(w, ℓ), :]
function degree_block(w::ModeWeights{T, IT}, ℓ) where {T, IT}
    if !(w.ℓₘᵢₙ ≤ ℓ ≤ w.ℓₘₐₓ)
        throw(BoundsError(w, (ℓ, :)))
    end
    ℓ = convert(IT, ℓ)
    i₀ = Yindex(ℓ, -ℓ, w.ℓₘᵢₙ)
    check_storage(w, i₀ + 2ℓ)
    DegreeBlock(view(w.data, i₀:i₀+2ℓ), ℓ)
end

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
    # The result is allocated at the length of the input's storage, so this one check covers
    # both of the vectors that the kernel indexes.
    check_storage_length(w)
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
    check_storage_length(w)
    check_storage_length(w′)
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
# Bare storage, at least as long as the result, is accepted as the output too, and the result
# comes back labelled, as a `ModeWeights` over it (see `mode_weights_view`).
function LinearAlgebra.mul!(w′::AbstractVector, op::DifferentialOperator, w::ModeWeights)
    mul!(mode_weights_view(w′, w.s + Δspin(op), w.ℓₘᵢₙ, w.ℓₘₐₓ), op, w)
end

# The in-place operations that write mode weights accept, as their output, a bare vector at
# least as long as the result, and return the result as a `ModeWeights` over its first entries
# — always a view, even when the length is exact, so that the storage is shared rather than
# copied and the type returned does not depend on the length.
function mode_weights_view(v::AbstractVector, s, ℓₘᵢₙ, ℓₘₐₓ)
    Base.require_one_based_indexing(v)
    n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
    if length(v) < n
        error(
            "The output has length $(length(v)); at least Ysize($ℓₘᵢₙ, $ℓₘₐₓ) = $n is needed."
        )
    end
    ModeWeights(view(v, 1:n), s, ℓₘᵢₙ, ℓₘₐₓ)
end
