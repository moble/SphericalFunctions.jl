### Products of the labelled containers.
#
#   𝔇 * w   rigidly rotates mode weights:  f′(𝐐) = f(𝐑⁻¹ 𝐐)
#   Y * w   evaluates the function at the rotor (or rotors) the harmonics were computed at
#
# Both are *bilinear* — neither conjugates anything — which is why both are `*` and neither is
# `⋅`; see the `dot` methods at the bottom of this file.
#
# The rotation law is derived in the conventions section, under "Rotation of mode weights":
#
#     f′_{ℓ,m′} = Σ_m 𝔇^{(ℓ)}_{m′,m}(𝐑) f_{ℓ,m}
#
# — no conjugate, summed over the *second* index of 𝔇, i.e. an ordinary matrix–vector product
# on each ℓ block.  Version 2's 𝔇 was the conjugate of this one, so code ported from it must
# *drop* a `conj`, not add one.


### Checks shared by the products below.
#
# The index *kind* is checked first, always.  `HalfOddInteger` and `Integer` deliberately do
# not promote — but `≤` between them is well defined and returns an ordinary `Bool`, so a kind
# mismatch sails straight through an ℓ-range test and surfaces much later as an `InexactError`
# from `convert`, or as a complaint from inside a loop.

index_kind_name(::Type{<:Integer}) = "integers"
index_kind_name(::Type{HalfOddInteger}) = "half-odd-integers"

function check_same_kind(::Type{IT}, w::ModeWeights{T, JT}, what) where {IT<:IntegerHalf, T, JT}
    isindex(IT, ℓₘᵢₙ(w)) && return nothing
    error(
        "These mode weights are indexed by $(index_kind_name(JT)) — "
        * "ℓ ∈ $(ℓₘᵢₙ(w)):$(ℓₘₐₓ(w)) — but $what is indexed by $(index_kind_name(IT)); "
        * "the two must be of one kind."
    )
end

# Containment, not equality: `D` has no `ℓₘᵢₙ` argument and always starts at 0 (or 1/2), while
# a `ModeWeights` usually starts at `abs(s)`.  Truncation the other way is never silent.
function check_ℓ_covers(lo, hi, w::ModeWeights, what, fix)
    if lo > ℓₘᵢₙ(w) || hi < ℓₘₐₓ(w)
        error(
            "The mode weights cover ℓ ∈ $(ℓₘᵢₙ(w)):$(ℓₘₐₓ(w)); the ℓ range of $what is only "
            * "$lo:$hi.  $fix"
        )
    end
    nothing
end

# The run of `w`'s modes inside the flat storage of a container whose ℓ range contains it.
@inline function shared_mode_range(container, w::ModeWeights)
    i₀ = Yindex(ℓₘᵢₙ(w), -ℓₘᵢₙ(w), ℓₘᵢₙ(container))
    i₀:(i₀ + length(w) - 1)
end


### Rotation.

# `eltype` of a `WignerSeries` is the *block* type — a series is a vector of blocks, not of
# numbers — so the number type is two `eltype`s deep.  Getting this wrong would quietly ask
# `similar` for a `Vector{WignerMatrix{…}}`.
number_type(𝔇::WignerSeries) = eltype(eltype(𝔇))
number_type(::WignerCalculator{IT, RT, NT}) where {IT, RT, NT} = NT

# A rotation mixes every m into every m′, so a block built with any of the m′/m restrictions
# cannot rotate anything: the missing rows and columns are not zero, they are absent.  The
# limits are `max(-ℓ, m′ₘᵢₙ):min(ℓ, m′ₘₐₓ)`, so a restriction is invisible below its own
# threshold — which is why a series is checked per ℓ rather than once.
function check_whole_block(B, ℓ)
    if m′ₘᵢₙ(B) != -ℓ || m′ₘₐₓ(B) != ℓ || mₘᵢₙ(B) != -ℓ || mₘₐₓ(B) != ℓ
        error(
            "The Wigner matrix for ℓ=$ℓ holds only m′ ∈ $(m′ₘᵢₙ(B)):$(m′ₘₐₓ(B)) and "
            * "m ∈ $(mₘᵢₙ(B)):$(mₘₐₓ(B)), but rotating mode weights needs the whole "
            * "(2ℓ+1)×(2ℓ+1) block, because every m mixes into every m′.  Rebuild the "
            * "matrices without the m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ or mₘᵢₙ restrictions."
        )
    end
    nothing
end

function check_rotation(𝔇::WignerSeries{IT}, w::ModeWeights) where {IT}
    check_same_kind(IT, w, "this WignerSeries")
    check_ℓ_covers(
        ℓₘᵢₙ(𝔇), ℓₘₐₓ(𝔇), w, "this WignerSeries",
        "Build the matrices with `D(R, $(ℓₘₐₓ(w)))`."
    )
end

function check_rotation(calc::WignerCalculator{IT}, w::ModeWeights) where {IT}
    check_same_kind(IT, w, "this calculator")
    check_ℓ_covers(
        ℓₘᵢₙ(calc), ℓₘₐₓ(calc), w, "this calculator", "Build it with ℓₘₐₓ=$(ℓₘₐₓ(w))."
    )
    if Nᵣ(calc) != 1
        error(
            "This calculator handles Nᵣ=$(Nᵣ(calc)) rotors at once, but a rotation of mode "
            * "weights is by one rotor; build the calculator with a single `Rotor`, or step "
            * "through the batch with `set_R!`."
        )
    end
    # The limits are fields of a calculator, so one check covers every ℓ.
    let L = ℓₘₐₓ(w)
        if m′ₘₐₓ(calc) < L || m′ₘᵢₙ(calc) > -L || mₘₐₓ(calc) < L || mₘᵢₙ(calc) > -L
            error(
                "This calculator computes only m′ ∈ $(m′ₘᵢₙ(calc)):$(m′ₘₐₓ(calc)) and "
                * "m ∈ $(mₘᵢₙ(calc)):$(mₘₐₓ(calc)), but rotating mode weights up to "
                * "ℓₘₐₓ=$L needs the whole block for every ℓ.  Build it without those "
                * "restrictions."
            )
        end
    end
    nothing
end

function check_rotation_output(w′::ModeWeights, w::ModeWeights)
    if spin(w′) != spin(w) || ℓₘᵢₙ(w′) != ℓₘᵢₙ(w) || ℓₘₐₓ(w′) != ℓₘₐₓ(w)
        error(
            "The output has s=$(spin(w′)) and ℓ ∈ $(ℓₘᵢₙ(w′)):$(ℓₘₐₓ(w′)), but a rotation "
            * "changes neither, so it must have s=$(spin(w)) and ℓ ∈ $(ℓₘᵢₙ(w)):$(ℓₘₐₓ(w))."
        )
    end
    if Base.mightalias(array_view(w′), array_view(w))
        error(
            "The output aliases the input.  A rotation mixes every m of a block into every "
            * "m′, so it cannot be done in place; pass a separate destination, such as "
            * "`similar(w)`."
        )
    end
    nothing
end


"""
    𝔇 * w
    calc * w

Rigidly rotate the mode weights `w`, giving the weights of ``f′(𝐐) = f(𝐑^{-1} 𝐐)`` — the
function actively rotated by the rotor whose Wigner matrices `𝔇` (a [`WignerSeries`](@ref)
from [`D`](@ref)) or `calc` (a [`DCalculator`](@ref)) hold.  Mode by mode,
```math
f′_{ℓ,m′} = \\sum_m 𝔇^{(ℓ)}_{m′,m}(𝐑)\\, f_{ℓ,m},
```
with **no** complex conjugate: version 2's ``𝔇`` was the conjugate of this one, so code ported
from it must drop a `conj` rather than add one.  The derivation is in the conventions section,
under [Rotation of mode weights](@ref conv_rotation_of_modes).

The result is a new `ModeWeights` with the same spin weight and the same range of ``ℓ`` as `w`
— a rotation changes neither.  `𝔇` must *cover* `w`'s range, which is not the same as matching
it: [`D`](@ref) has no `ℓₘᵢₙ` argument and always starts at ``0`` (or ``1/2``), while a
`ModeWeights` usually starts at ``|s|``.  The blocks must be whole: a ``𝔇`` built with any of
the `m′ₘₐₓ`, `m′ₘᵢₙ`, `mₘₐₓ` or `mₘᵢₙ` restrictions cannot rotate anything, because every ``m``
mixes into every ``m′``, and is refused.

Real matrices from [`d`](@ref) are accepted, and are exactly the rotation by
`from_euler_angles(0, β, 0)`, for which ``𝔇 = d``.

The calculator form holds one ``ℓ`` at a time instead of materializing every block, so it
allocates only the result; `mul!` into an existing container allocates nothing at all.  A
calculator holds workspace, so it must not be used from two threads at once.

```julia
w′ = D(R, ℓₘₐₓ(w)) * w
w′(Q) ≈ w(inv(R) * Q)          # the defining property
```

See also [`ModeWeights`](@ref) and [`HarmonicValues`](@ref).
"""
function Base.:*(𝔇::WignerSeries, w::ModeWeights)
    check_rotation(𝔇, w)
    w′ = similar(w, promote_type(number_type(𝔇), eltype(w)))
    rotate_modes!(array_view(w′), 𝔇, w)
    w′
end

function Base.:*(calc::WignerCalculator, w::ModeWeights)
    check_rotation(calc, w)
    w′ = similar(w, promote_type(number_type(calc), eltype(w)))
    rotate_modes!(array_view(w′), calc, w)
    w′
end

function LinearAlgebra.mul!(w′::ModeWeights, 𝔇::WignerSeries, w::ModeWeights)
    check_rotation(𝔇, w)
    check_rotation_output(w′, w)
    rotate_modes!(array_view(w′), 𝔇, w)
    w′
end

function LinearAlgebra.mul!(w′::ModeWeights, calc::WignerCalculator, w::ModeWeights)
    check_rotation(calc, w)
    check_rotation_output(w′, w)
    rotate_modes!(array_view(w′), calc, w)
    w′
end

# The workers, against the flat storage of both containers: one `mul!` per ℓ, on the
# contiguous run of modes that ℓ occupies.  Both ends of the loop come from `w`, because a
# range formed from two different objects' accessors is a `MethodError` waiting to happen when
# one side is a `HalfOddInteger` and the other is not.
function rotate_modes!(dst::AbstractVector, 𝔇::WignerSeries, w::ModeWeights)
    src = array_view(w)
    for ℓ ∈ ℓₘᵢₙ(w):ℓₘₐₓ(w)
        B = 𝔇[ℓ]
        check_whole_block(B, ℓ)
        r = mode_range(w, ℓ)
        mul!(view(dst, r), array_view(B), view(src, r))
    end
    dst
end

# Streaming: nothing is held but the block the calculator is standing on.  Beginning above the
# calculator's own ℓₘᵢₙ still runs the recurrence through the ℓ below, so the result is
# bit-for-bit what a sequential pass gives.  The m′/m limits were checked once, in
# `check_rotation`, because they are fields of the calculator rather than of each block.
function rotate_modes!(dst::AbstractVector, calc::WignerCalculator, w::ModeWeights)
    src = array_view(w)
    for ℓ ∈ ℓₘᵢₙ(w):ℓₘₐₓ(w)
        B = recurrence!(calc, ℓ)
        r = mode_range(w, ℓ)
        mul!(view(dst, r), array_view(B), view(src, r))
    end
    dst
end


### Evaluation.

# The contraction over modes, written out rather than reached through `dot`, which would
# conjugate.  `eachindex(y, x)` is also the length check, and `init` is what lets an empty
# container give zero rather than throw.  The association is `mapfoldl`'s, which is what the
# `w(R)` of earlier versions used, so the two agree bit for bit.
@inline function synthesize(y::AbstractVector, x::AbstractVector)
    T = promote_type(eltype(y), eltype(x))
    sum(y[i] * x[i] for i ∈ eachindex(y, x); init=zero(T))
end

# A difference of two `HalfOddInteger`s is an `Int` by construction, so this is the same
# expression for either kind of index.
@inline spin_row(Y, w::ModeWeights) = Int(spin(w) - first(spins(Y))) + 1

function check_spin_available(sr, w::ModeWeights, what)
    if !(first(sr) ≤ spin(w) ≤ last(sr))
        error(
            "These mode weights have spin weight s=$(spin(w)), but $what serves "
            * (length(sr) == 1 ? "s=$(first(sr))" : "s ∈ $(first(sr)):$(last(sr))")
            * ".  Evaluation pairs the weights of a function with the harmonics of its own "
            * "spin weight."
        )
    end
    nothing
end

function check_evaluation(Y::HarmonicValues{T, IT}, w::ModeWeights) where {T, IT}
    check_same_kind(IT, w, "these harmonic values")
    check_ℓ_covers(
        ℓₘᵢₙ(Y), ℓₘₐₓ(Y), w, "these harmonic values",
        "Compute them with `sYlm(R, $(ℓₘₐₓ(w)), $(spin(w)); ℓₘᵢₙ=$(ℓₘᵢₙ(w)))`."
    )
    check_spin_available(spins(Y), w, "these harmonic values")
end

function check_evaluation(calc::HarmonicCalculator{IT}, w::ModeWeights) where {IT}
    check_same_kind(IT, w, "this calculator")
    check_ℓ_covers(
        ℓₘᵢₙ(calc), ℓₘₐₓ(calc), w, "this calculator", "Build it with ℓₘₐₓ=$(ℓₘₐₓ(w))."
    )
    check_spin_available(spins(calc), w, "this calculator")
end


"""
    Y * w
    calc * w

Evaluate the function with mode weights `w` at the rotor (or rotors) whose harmonics `Y`
holds, ``f(𝐑) = \\sum_{ℓ,m} f_{ℓ,m}\\, {}_sY_{ℓ,m}(𝐑)``.  `Y` is a [`HarmonicValues`](@ref)
from [`sYlm`](@ref), or an [`sYlmCalculator`](@ref).

The result is a scalar when the harmonics were computed for a single rotor, and a `Vector` of
one value per rotor otherwise.  Where `Y` holds a *range* of spin weights, the row of `w`'s own
spin weight is the one used, and the result is as above; weights of one spin weight paired
against harmonics of another mean nothing, so a `w` whose spin weight is not among them is
refused.

`Y`'s range of ``ℓ`` must *cover* `w`'s, which is not the same as matching it: [`sYlm`](@ref)
takes an `ℓₘᵢₙ`, and the two need not agree.

This product does **not** conjugate anything, which is why it is `*` and not `⋅`: `⋅` is
`LinearAlgebra.dot`, which conjugates its first argument.  It is the same product that
[`sYlm_matrix`](@ref)'s docstring writes as `f = Y * f̃`.

The calculator form holds one ``ℓ`` at a time rather than every harmonic at once.  A
calculator holds workspace, so it must not be used from two threads at once.

```julia
Y = sYlm(R, ℓₘₐₓ(w), spin(w); ℓₘᵢₙ=ℓₘᵢₙ(w))
Y * w == w(R)
```
"""
function Base.:*(
    Y::HarmonicValues{T, IT, S, <:AbstractVector}, w::ModeWeights
) where {T, IT, S<:IntegerHalf}
    check_evaluation(Y, w)
    synthesize(view(array_view(Y), shared_mode_range(Y, w)), array_view(w))
end
function Base.:*(
    Y::HarmonicValues{T, IT, S, <:AbstractMatrix}, w::ModeWeights
) where {T, IT, S<:IntegerHalf}
    check_evaluation(Y, w)
    view(array_view(Y), :, shared_mode_range(Y, w)) * array_view(w)
end
function Base.:*(
    Y::HarmonicValues{T, IT, S, <:AbstractMatrix}, w::ModeWeights
) where {T, IT, S<:AbstractUnitRange}
    check_evaluation(Y, w)
    synthesize(view(array_view(Y), spin_row(Y, w), shared_mode_range(Y, w)), array_view(w))
end
function Base.:*(
    Y::HarmonicValues{T, IT, S, <:AbstractArray{T, 3}}, w::ModeWeights
) where {T, IT, S<:AbstractUnitRange}
    check_evaluation(Y, w)
    view(array_view(Y), :, spin_row(Y, w), shared_mode_range(Y, w)) * array_view(w)
end

# Streaming.  Batchedness is a type parameter, so these dispatch on it rather than branching
# inside; each then has one concrete return type — a scalar, or one value per rotor — instead
# of a union of the two.
function Base.:*(
    calc::HarmonicCalculator{IT, RT, YT, ST, S, false}, w::ModeWeights
) where {IT, RT, YT, ST, S}
    check_evaluation(calc, w)
    src = array_view(w)
    f = zero(promote_type(YT, eltype(w)))
    iₛ = spin_index(calc, convert(IT, spin(w)))
    for ℓ ∈ ℓₘᵢₙ(w):ℓₘₐₓ(w)
        recurrence!(calc, ℓ)
        Yˡ = spin_row(calc, ℓ, iₛ)
        f += synthesize(array_view(Yˡ), view(src, mode_range(w, ℓ)))
    end
    f
end
function Base.:*(
    calc::HarmonicCalculator{IT, RT, YT, ST, S, true}, w::ModeWeights
) where {IT, RT, YT, ST, S}
    check_evaluation(calc, w)
    src = array_view(w)
    NT = promote_type(YT, eltype(w))
    f = zeros(NT, Nᵣ(calc))
    iₛ = spin_index(calc, convert(IT, spin(w)))
    for ℓ ∈ ℓₘᵢₙ(w):ℓₘₐₓ(w)
        recurrence!(calc, ℓ)
        # β = 1 accumulates across ℓ, so nothing beyond `f` is ever held.
        mul!(f, array_view(spin_row(calc, ℓ, iₛ)), view(src, mode_range(w, ℓ)), one(NT), one(NT))
    end
    f
end


### Evaluation at a rotor, spelled as a call.

"""
    w(R)
    w(R⃗)

Evaluate the function with mode weights `w` at the rotor `R`, or at each of a vector of rotors
`R⃗`, ``f(𝐑) = \\sum_{ℓ,m} f_{ℓ,m}\\, {}_sY_{ℓ,m}(𝐑)``.  For spherical coordinates use
`R = from_spherical_coordinates(θ, ϕ)`.

This computes the harmonics afresh on every call.  For repeated evaluation build an
[`sYlmCalculator`](@ref) once and use `calc * w`; for many rotors at once see
[`sYlm_matrix`](@ref) and the transforms in the "Transformations" section.
"""
(w::ModeWeights)(R::Rotor) = sYlm(R, ℓₘₐₓ(w), spin(w); ℓₘᵢₙ=ℓₘᵢₙ(w)) * w
(w::ModeWeights)(R⃗::AbstractVector{<:Rotor}) = sYlm(R⃗, ℓₘₐₓ(w), spin(w); ℓₘᵢₙ=ℓₘᵢₙ(w)) * w


### `dot` is *not* evaluation.
#
# `⋅` is `LinearAlgebra.dot`, which conjugates its first argument; evaluation must not.  The
# package's own `dot(::ModeWeights, ::ModeWeights)` *is* the conjugating inner product, and
# `dot(array_view(Y), array_view(w))` conjugates too — so a non-conjugating `dot` here would make `⋅`
# mean two different things a few lines apart, with the disagreement showing up only as a wrong
# phase.  That is exactly the silent-wrong-answer failure these containers exist to prevent
# (see the header of `array_view.jl`).  These methods say so, rather than leaving a bare
# `MethodError` for someone to "fix" later by adding the harmful method.
const dot_is_not_evaluation = (
    "`dot` (`⋅`) conjugates its first argument, but evaluating a spin-weighted function — "
    * "f(𝐑) = Σ f_{ℓ,m} ₛYₗₘ(𝐑) — must not conjugate the harmonics.  Use `Y * w`, which is "
    * "the bilinear product that `sYlm_matrix` documents as `f = Y * f̃`.  For the conjugating "
    * "inner product of two sets of mode weights, `dot(w₁, w₂)` is still what you want."
)

LinearAlgebra.dot(::Union{HarmonicValues, HarmonicCalculator}, ::ModeWeights) =
    error(dot_is_not_evaluation)
LinearAlgebra.dot(::ModeWeights, ::Union{HarmonicValues, HarmonicCalculator}) =
    error(dot_is_not_evaluation)
