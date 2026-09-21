# The last parameter, `B`, is `Nᵣ > 1`.  It is redundant — `Nᵣ` is a field of the wedge
# inside `H`, and `isbatched(c)` could simply compare it to 1 — but it is the difference
# between `calc[ℓ]` having one return type and having a union of the batched and unbatched
# ones, because the branch in `block` below is then resolved at compile time.  Only the
# *predicate* is lifted into the type, not `Nᵣ` itself: `SSHTRS` sets `Nᵣ = Nθ`, which grows
# with ℓₘₐₓ, and parameterizing on the count would recompile the recurrence for every
# resolution.

"""
    WignerCalculator{IT, RT, NT}

Calculator producing Wigner's ``𝔇`` matrices (when `NT` is `Complex{RT}`) or ``d`` matrices
(when `NT` is `RT`) for `Nᵣ` rotors at a time, one ``ℓ`` at a time.  Use the constructors
[`DCalculator`](@ref) and [`dCalculator`](@ref).

Internally this wraps a [`HCalculator`](@ref), which does the actual recurrence, plus a
buffer into which the requested block of the matrix is written for the current ``ℓ``; that
block is returned by `calc[ℓ]` as an array indexed naturally by `[m′, m]` (or `[iᵣ, m′, m]`
when `Nᵣ > 1`).

Which of those two shapes `calc[ℓ]` returns is recorded in the type, as the `Bool` parameter
read by [`isbatched`](@ref), so that the return type is inferrable; see the comment on the
struct.
"""
struct WignerCalculator{IT, RT<:Real, NT<:Union{RT, Complex{RT}}, ST, B}
    H::HCalculator{IT, RT, ST}
    Wˡ::Array{NT, 3}  # [iᵣ, m′, m] block for the current ℓ, using the leading entries
    Z₊::Matrix{Complex{RT}}  # Z₊[k+1, iᵣ] = z₊^k for k ∈ 0:2ℓₘₐₓ; empty when NT is real
    Z₋::Matrix{Complex{RT}}  # Z₋[k+1, iᵣ] = z₋^k for k ∈ 0:2ℓₘₐₓ; empty when NT is real
    m′ₘₐₓ::IT
    m′ₘᵢₙ::IT
    mₘₐₓ::IT
    mₘᵢₙ::IT
    ℓ::Base.RefValue{IT}  # ℓ of the block currently in Wˡ; ℓₘᵢₙ-1 if none
end

# Allocate the buffers without touching them.  PRIVATE: see the note on `allocate_H`.  The
# rotor data here is the nested `H`'s phases *plus* this calculator's own power tables `Z₊`
# and `Z₋`, so a caller that copies rather than sets must copy all of them.
function allocate_W(
    ::Type{IT}, ::Type{RT}, ::Type{NT}, ℓₘₐₓ::IT,
    m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT, Nᵣ::Int
) where {IT<:IntegerHalf, RT<:Real, NT<:Union{RT, Complex{RT}}}
    validate_index_ranges(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    # The recurrence needs the wedge for |m′| up to the larger of the two m′ limits (and all
    # m); the four limits only select the block that is materialized.
    H = allocate_H(IT, RT, ℓₘₐₓ, max(m′ₘₐₓ, -m′ₘᵢₙ), Nᵣ)
    Wˡ = Array{NT, 3}(undef, Nᵣ, Int(m′ₘₐₓ - m′ₘᵢₙ) + 1, Int(mₘₐₓ - mₘᵢₙ) + 1)
    K = NT <: Complex ? Int(2ℓₘₐₓ) + 1 : 0
    Z₊ = Matrix{Complex{RT}}(undef, K, Nᵣ)
    Z₋ = Matrix{Complex{RT}}(undef, K, Nᵣ)
    # The reference is typed explicitly for the same reason as in `allocate_Y`: for a narrower
    # integer `IT`, `ℓₘᵢₙ(IT) - 1` is an `Int`, and the field is a `RefValue{IT}`.
    WignerCalculator{IT, RT, NT, typeof(parent(H.Hˡ)), Nᵣ > 1}(
        H, Wˡ, Z₊, Z₋, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, Ref{IT}(ℓₘᵢₙ(IT) - 1)
    )
end

function WignerCalculator{IT, RT, NT}(
    R, ℓₘₐₓ::IT;
    m′ₘₐₓ::IT=ℓₘₐₓ, m′ₘᵢₙ::IT=-m′ₘₐₓ, mₘₐₓ::IT=ℓₘₐₓ, mₘᵢₙ::IT=-mₘₐₓ
) where {IT<:IntegerHalf, RT<:Real, NT<:Union{RT, Complex{RT}}}
    set_rotors!(
        allocate_W(IT, RT, NT, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, nrotors(R)), R
    )
end

"""
    DCalculator(R, ℓₘₐₓ; m′ₘₐₓ=ℓₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ, mₘₐₓ=ℓₘₐₓ, mₘᵢₙ=-mₘₐₓ)

Calculator for Wigner's ``𝔇^{(ℓ)}_{m′,m}(R)`` matrices, for ``ℓ ≤ ℓₘₐₓ``, with elements of
type `Complex{RT}`.  The keyword arguments restrict the block of each matrix that is computed
and returned.

The rotor is given at construction, so that the calculator is usable the moment it exists,
and so that the element type is the rotor's own: a `Rotor{Float32}` gives a `Float32`
calculator, and `floattype(calc)` reports it.  There is no argument to override that — to
compute in another type, convert the rotor, which is the honest way to say that those are the
values to be treated as exact.  Give `R` as an `AbstractVector` of rotors to
get a calculator that handles all `Nᵣ = length(R)` of them at once — which is substantially
faster per rotor, for the reason described under
[Reusing the workspace](@ref interface_wigner_matrices).  Later rotors are supplied with
[`set_R!`](@ref).

The calculator is iterable, yielding one ``ℓ`` at a time:

```julia
calc = DCalculator(R, ℓₘₐₓ)
for (ℓ, 𝔇ˡ) ∈ calc
    # 𝔇ˡ[m′, m] with m′, m ∈ -ℓ:ℓ
end
```

With `Nᵣ > 1` each block is indexed as `[iᵣ, m′, m]` instead.  The block is a view into the
calculator's storage and is overwritten by the next step, so `copy` it if it must survive
(the copy keeps the natural indices), or `collect` it to get an ordinary 1-based `Matrix`.
`collect(calc)` copies every block, so it is safe; see [`eachℓ`](@ref) to restrict the range
of ``ℓ``, and [`recurrence!`](@ref) to step the calculator by hand.

The convention is ``𝔇^{(ℓ)}_{m′,m}(𝐑_{α,β,γ}) = e^{-im′α}\\, d^{(ℓ)}_{m′,m}(β)\\, e^{-imγ}``;
see the "Conventions" section of the documentation.

# Half-integer indices

`DCalculator(R, 7//2)` — a `Rational` `ℓₘₐₓ` with denominator 2 — gives a
calculator for half-integer ``ℓ, m′, m``.  All four keyword limits must then be
half-integers too, `recurrence!` accepts only half-integer `ℓ`, and `calc[ℓ]` returns a
[`WignerMatrix`](@ref) (or a [`WignerMatrixBatch`](@ref) when `Nᵣ > 1`) whose indices are
half-odd-integers; it is indexed the same way.  The double cover is respected exactly: ``𝔇(-R) = -𝔇(R)``.

See also [`D`](@ref) for a simpler interface when the matrices for only one rotor are needed,
[`dCalculator`](@ref) for the real ``d`` matrices, and [`recurrence!`](@ref).
"""
const DCalculator{IT, RT, ST, B} = WignerCalculator{IT, RT, Complex{RT}, ST, B} where {IT, RT<:Real, ST, B}

"""
    dCalculator(β, ℓₘₐₓ; m′ₘₐₓ=ℓₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ, mₘₐₓ=ℓₘₐₓ, mₘᵢₙ=-mₘₐₓ)

Calculator for Wigner's real ``d^{(ℓ)}_{m′,m}(β)`` matrices, for ``ℓ ≤ ℓₘₐₓ``, with elements
of the angle's own floating-point type.  The first argument may be the angle ``β``, the phase ``e^{iβ}``, or a `Rotor`
(of which only the ``β`` Euler angle is used), or an `AbstractVector` of `Nᵣ` of any one of
those; later values are supplied with [`set_β!`](@ref).  Otherwise this behaves exactly like
[`DCalculator`](@ref) — including iteration, and half-integer ``ℓ`` for a `Rational`
`ℓₘₐₓ`.

Half-integer ``d`` has period ``4π`` in ``β``, so an angle or a `Rotor` determines it
unambiguously, while a bare phase ``e^{iβ}`` determines it only up to the double-cover sign
``(-1)^{2ℓ}``; the branch ``β ∈ (-π, π]`` is used in that case.

See also [`d`](@ref).
"""
const dCalculator{IT, RT, ST, B} = WignerCalculator{IT, RT, RT, ST, B} where {IT, RT<:Real, ST, B}

# `ℓₘₐₓ` is constrained to `IntegerHalf` here (with the `Rational` methods at the bottom of
# this file taking the half-integer spelling) so that a call in the old argument order —
# `DCalculator(ℓₘₐₓ, Float64)` — is an immediate `MethodError` at the call site rather
# than something that dispatches with the element type in the rotor's place.
# The element type is derived here and passed on as a *type*, to the helpers below, rather
# than computed inside the body as a value: that is what lets the compiler settle the concrete
# return type, including the `B` parameter that `calc[ℓ]`'s type depends on.
#
# Those helpers are deliberately *not* methods of `DCalculator` and `dCalculator`.
# A three-argument method of either name would be a public way to override the element type,
# and version 3 has none by design — the type of the rotor data is the only thing that decides
# it.  `test/wigner/iteration.jl` asserts exactly that, with `@test_throws MethodError`.
function DCalculator(R, ℓₘₐₓ::IT; kwargs...) where {IT<:IntegerHalf}
    wigner_D_calculator(R, ℓₘₐₓ, rotor_basetype(R); kwargs...)
end
function dCalculator(β, ℓₘₐₓ::IT; kwargs...) where {IT<:IntegerHalf}
    wigner_d_calculator(β, ℓₘₐₓ, rotor_basetype(β); kwargs...)
end
function wigner_D_calculator(R, ℓₘₐₓ::IT, ::Type{RT}; kwargs...) where {IT<:IntegerHalf, RT<:Real}
    WignerCalculator{IT, RT, Complex{RT}}(R, ℓₘₐₓ; kwargs...)
end
function wigner_d_calculator(β, ℓₘₐₓ::IT, ::Type{RT}; kwargs...) where {IT<:IntegerHalf, RT<:Real}
    WignerCalculator{IT, RT, RT}(β, ℓₘₐₓ; kwargs...)
end

# A second workspace holding the same rotor data, with nothing computed — which is what a
# thread wants.  The assertion is what keeps this inferrable: `Nᵣ(c)` is a field lookup, so
# the constructor cannot know `B`, but the copy necessarily has the same parameters as the
# original.  The data is copied buffer-by-buffer rather than re-derived; see the comment on
# `similar(::HCalculator)` for why it cannot be re-derived at all.
function Base.similar(c::WignerCalculator{IT, RT, NT, ST, B}) where {IT, RT, NT, ST, B}
    c′ = allocate_W(
        IT, RT, NT, ℓₘₐₓ(c), c.m′ₘₐₓ, c.m′ₘᵢₙ, c.mₘₐₓ, c.mₘᵢₙ, Nᵣ(c)
    )::WignerCalculator{IT, RT, NT, ST, B}
    copyto!(c′.H.eⁱᵝ, c.H.eⁱᵝ)
    copyto!(c′.H.cβ½, c.H.cβ½)
    copyto!(c′.H.sβ½, c.H.sβ½)
    copyto!(c′.Z₊, c.Z₊)
    copyto!(c′.Z₋, c.Z₋)
    c′
end
function Base.similar(c::WignerCalculator{IT, RT, NT, ST, B}, R) where {IT, RT, NT, ST, B}
    if nrotors(R) != Nᵣ(c)
        error("This calculator handles Nᵣ=$(Nᵣ(c)) rotors, but got $(nrotors(R)).")
    end
    check_rotor_type(c, R)
    set_rotors!(
        allocate_W(
            IT, RT, NT, ℓₘₐₓ(c), c.m′ₘₐₓ, c.m′ₘᵢₙ, c.mₘₐₓ, c.mₘᵢₙ, Nᵣ(c)
        )::WignerCalculator{IT, RT, NT, ST, B},
        R
    )
end

ℓ(c::WignerCalculator) = c.ℓ[]
ℓₘᵢₙ(c::WignerCalculator{IT}) where {IT} = ℓₘᵢₙ(IT)
ℓₘₐₓ(c::WignerCalculator) = ℓₘₐₓ(c.H)
m′ₘₐₓ(c::WignerCalculator) = c.m′ₘₐₓ
m′ₘᵢₙ(c::WignerCalculator) = c.m′ₘᵢₙ
mₘₐₓ(c::WignerCalculator) = c.mₘₐₓ
mₘᵢₙ(c::WignerCalculator) = c.mₘᵢₙ
Nᵣ(c::WignerCalculator) = Nᵣ(c.H)
floattype(::WignerCalculator{IT, RT}) where {IT, RT} = RT
isbatched(::WignerCalculator{IT, RT, NT, ST, B}) where {IT, RT, NT, ST, B} = B

function Base.show(io::IO, c::WignerCalculator{IT, RT, NT}) where {IT, RT, NT}
    print(
        io,
        NT <: Complex ? "D" : "d", "Calculator{$IT, $RT} for ",
        "ℓₘₐₓ=$(ℓₘₐₓ(c)), m′=$(c.m′ₘᵢₙ):$(c.m′ₘₐₓ), m=$(c.mₘᵢₙ):$(c.mₘₐₓ), Nᵣ=$(Nᵣ(c))",
        c.ℓ[] < ℓₘᵢₙ(c) ? " (nothing computed yet)" : ", currently at ℓ=$(c.ℓ[])"
    )
end
function Base.show(io::IO, ::MIME"text/plain", c::WignerCalculator)
    show(io, c)
end

"""
    fill!(c::WignerCalculator, v)

Fill the axis, wedge and output buffers of `c` with the value `v` and mark the current
results as invalid.  The stored rotor data — `e^{iβ}`, the half angles, and the phase powers
`Z₊`, `Z₋` — is deliberately *not* touched, so `recurrence!(c, ℓ)` still has everything it
needs, exactly as for [`HCalculator`](@ref).  Useful for testing
that no uninitialized storage is ever read: everything the recurrence is responsible for
writing is poisoned, while everything `set_rotors!` is responsible for writing is left
alone.
"""
function Base.fill!(c::WignerCalculator{IT, RT, NT}, v::Number) where {IT, RT, NT}
    fill!(c.H, real(v))
    fill!(c.Wˡ, convert(NT, v))
    c.ℓ[] = ℓₘᵢₙ(IT) - 1
    c
end


### Rotor data

# For 𝔇 we need the full rotor: eⁱᵝ for the recurrence, and the powers of z₊ and z₋ for the
# phases.
function set_rotors!(
    c::WignerCalculator{IT, RT, Complex{RT}}, R::AbstractVector{<:Rotor}
) where {IT, RT<:Real}
    if length(R) != Nᵣ(c)
        error("Expected $(Nᵣ(c)) rotors (Nᵣ), but got $(length(R)).")
    end
    @inbounds for i ∈ eachindex(R)
        eⁱᵝ, z₊, z₋, cβ½, sβ½ = spinor_phases(R[i], RT)
        c.H.eⁱᵝ[i] = eⁱᵝ
        set_half_angles!(c.H, i, cβ½, sβ½)
        complex_powers!(view(c.Z₊, :, i), z₊)
        complex_powers!(view(c.Z₋, :, i), z₋)
    end
    c.H.axes_valid[] = false
    c.ℓ[] = ℓₘᵢₙ(IT) - 1
    c
end
function set_rotors!(c::WignerCalculator{IT, RT, Complex{RT}}, R::Rotor) where {IT, RT<:Real}
    if Nᵣ(c) != 1
        error("A single rotor was given, but this calculator expects Nᵣ=$(Nᵣ(c)) rotors.")
    end
    set_rotors!(c, @SVector [R])
end
function set_rotors!(c::WignerCalculator{IT, RT, Complex{RT}}, R) where {IT, RT<:Real}
    error(
        "A DCalculator needs rotors, given as `Rotor`s — one, or an AbstractVector "
        * "of $(Nᵣ(c)) of them — not $(typeof(R)); use a dCalculator if only β is "
        * "available."
    )
end

# For d we need only eⁱᵝ, and accept anything the H calculator accepts.
function set_rotors!(c::WignerCalculator{IT, RT, RT}, R) where {IT, RT<:Real}
    set_rotors!(c.H, R)
    c.ℓ[] = ℓₘᵢₙ(IT) - 1
    c
end


### Driver

function recurrence!(c::WignerCalculator, R, ℓ)
    check_ℓ(c.H, ℓ)
    set_rotors!(c, R)
    recurrence!(c, ℓ)
end
function recurrence!(c::WignerCalculator{IT}, ℓ) where {IT}
    let ℓ = convert(IT, ℓ)
        recurrence!(c.H, ℓ)
        materialize!(c, ℓ)
        current_block(c, ℓ)
    end
end

# The block for the ``ℓ`` just computed, restricted to the m′ and m limits the calculator was
# built with.
current_block(c::WignerCalculator{IT}, ℓ::IT) where {IT} =
    block(c, ℓ, m′range(c, ℓ), mrange(c, ℓ))

# Ranges of m′ and m in the block for a given ℓ
m′range(c::WignerCalculator, ℓ) = max(-ℓ, c.m′ₘᵢₙ):min(ℓ, c.m′ₘₐₓ)
mrange(c::WignerCalculator, ℓ) = max(-ℓ, c.mₘᵢₙ):min(ℓ, c.mₘₐₓ)

# Powers zᵏ for k of either sign, given Z[k+1, iᵣ] = zᵏ for k ≥ 0
@inline zpower(Z, iᵣ, k) = k ≥ 0 ? (@inbounds Z[k+1, iᵣ]) : conj(@inbounds Z[1-k, iᵣ])

# Write the block of the d matrix (real NT) or 𝔇 matrix (complex NT) for the current ℓ into
# c.Wˡ, reading the H wedge through `wedge_source` (the only place the symmetries of H are
# encoded) and applying the ϵ signs relating H to d, and for 𝔇 the phases e^{-i(m′α+mγ)}.
#
# Nothing here depends on whether the indices are integers or half-odd-integers: the
# exponents m′±m of z₊ and z₋ are `Integer`s either way (see the v3 design memo, §5.4), and
# ϵ and σ are already general.
function materialize!(c::WignerCalculator{IT, RT, NT}, ℓ::IT) where {IT, RT, NT}
    let H = c.H.Hˡ, Wˡ = c.Wˡ, Z₊ = c.Z₊, Z₋ = c.Z₋, Nᵣ = Nᵣ(c), Hp = parent(H)
        if H.ℓ != ℓ
            error("The H wedge holds ℓ=$(H.ℓ), but ℓ=$ℓ was requested.")
        end
        m′ₘₐₓw = m′ₘₐₓ(H)
        m′ₘᵢₙw = m′ₘᵢₙ(H)
        mlo = max(-ℓ, c.mₘᵢₙ)
        mhi = min(ℓ, c.mₘₐₓ)
        m′lo = max(-ℓ, c.m′ₘᵢₙ)
        m′hi = min(ℓ, c.m′ₘₐₓ)
        @inbounds for (j, m) ∈ enumerate(mlo:mhi)
            for (j′, m′) ∈ enumerate(m′lo:m′hi)
                a, b, σ = wedge_source(m′, m, m′ₘₐₓw)
                offset = wedge_offset(H, a, b, m′ₘᵢₙw)
                # d = ϵ(m′) ϵ(-m) H, with the transposition sign σ
                coefficient = convert(RT, σ * ϵ(m′) * ϵ(-m))
                if NT <: Real
                    @simd for iᵣ ∈ 1:Nᵣ
                        Wˡ[iᵣ, j′, j] = coefficient * Hp[offset + iᵣ]
                    end
                else
                    # 𝔇 = d e^{-i(m′α + mγ)} = d conj(z₊^(m′+m) z₋^(m′-m))
                    k₊ = m′ + m
                    k₋ = m′ - m
                    for iᵣ ∈ 1:Nᵣ
                        phase = conj(zpower(Z₊, iᵣ, k₊) * zpower(Z₋, iᵣ, k₋))
                        Wˡ[iᵣ, j′, j] = coefficient * Hp[offset + iᵣ] * phase
                    end
                end
            end
        end
    end
    c.ℓ[] = ℓ
    c
end

"""
    calc[ℓ]

The block of the Wigner matrix for the current ``ℓ`` of the calculator, as an array indexed
naturally: `calc[ℓ][m′, m]` for `Nᵣ == 1`, or `calc[ℓ][iᵣ, m′, m]` for `Nᵣ > 1`.  The
result is a view into the calculator's storage, valid until the next call to
[`recurrence!`](@ref).  `ℓ` must be the value passed to the most recent `recurrence!`.

The result is a [`WignerMatrix`](@ref) (for `Nᵣ == 1`) or a [`WignerMatrixBatch`](@ref), for
either kind of index, supporting `[m′, m]` and `[iᵣ, m′, m]` respectively.  `copy` keeps it,
with its natural indices, beyond the next `recurrence!`; `collect` gives an ordinary 1-based
`Array`; and [`strided`](@ref) gives a 1-based view of the same storage, which is what linear
algebra takes.
"""
function Base.getindex(c::WignerCalculator{IT}, ℓ) where {IT}
    ℓ = convert(IT, ℓ)
    if c.ℓ[] < ℓₘᵢₙ(c)
        error(
            "This calculator currently holds no result, because nothing has been computed "
            * "yet; iterate the calculator, or call `recurrence!(calc, ℓ)` first."
        )
    end
    if ℓ < ℓₘᵢₙ(c) || ℓ > ℓₘₐₓ(c)
        error(
            "ℓ=$ℓ is out of bounds [$(ℓₘᵢₙ(c)), $(ℓₘₐₓ(c))] for this calculator; "
            * "`recurrence!` accepts only ℓ in that range."
        )
    end
    if ℓ != c.ℓ[]
        error(
            "This calculator currently holds ℓ=$(c.ℓ[]), not ℓ=$ℓ; "
            * "call `recurrence!(calc, $ℓ)` first."
        )
    end
    m′r = m′range(c, ℓ)
    mr = mrange(c, ℓ)
    block(c, ℓ, m′r, mr)
end

# `isbatched(c)` reads the type parameter, so this branch is resolved at compile time and the
# method has a single concrete return type.  The same containers are returned for integer and
# half-odd-integer indices alike; see the note on `AbstractWignerMatrix` for why they are not
# `OffsetArray`s even where an `OffsetArray` could represent them.
function block(c::WignerCalculator{IT}, ℓ::IT, m′r, mr) where {IT<:IntegerHalf}
    if isbatched(c)
        WignerMatrixBatch(
            view(c.Wˡ, :, 1:length(m′r), 1:length(mr)), ℓ;
            m′ₘₐₓ=last(m′r), m′ₘᵢₙ=first(m′r), mₘₐₓ=last(mr), mₘᵢₙ=first(mr)
        )
    else
        WignerMatrix(
            view(c.Wˡ, 1, 1:length(m′r), 1:length(mr)), ℓ;
            m′ₘₐₓ=last(m′r), m′ₘᵢₙ=first(m′r), mₘₐₓ=last(mr), mₘᵢₙ=first(mr)
        )
    end
end


### Convenience functions

"""
    D(R, ℓₘₐₓ; m′ₘₐₓ=ℓₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ, mₘₐₓ=ℓₘₐₓ, mₘᵢₙ=-mₘₐₓ)

Wigner's ``𝔇^{(ℓ)}_{m′,m}(R)`` matrices for all ``ℓ ≤ ℓₘₐₓ``, for the single rotor `R`.

The result is a [`WignerSeries`](@ref), indexed by ``ℓ`` and then naturally by ``(m′, m)``:
`D(R, ℓₘₐₓ)[ℓ][m′, m]`.  Each block is a [`WignerMatrix`](@ref); `Matrix` (or `collect`) gives
the ordinary ``(2ℓ+1)×(2ℓ+1)`` `Matrix` with rows and columns in order of increasing ``m′``
and ``m``.
The keyword arguments restrict the block of each matrix that is computed.  The convention is
``𝔇^{(ℓ)}_{m′,m}(𝐑_{α,β,γ}) = e^{-im′α}\\, d^{(ℓ)}_{m′,m}(β)\\, e^{-imγ}``; see the
"Conventions" section of the documentation.

`ℓₘₐₓ` may be a half-integer `Rational` (denominator 2), as may the keyword limits; then
``ℓ`` runs over ``1/2, 3/2, …, ℓₘₐₓ`` rather than ``0, 1, …, ℓₘₐₓ``.  The container types are
the same either way, and `D(R, ℓₘₐₓ)[ℓ][m′, m]` reads the same.

This function allocates all of its output on every call.  To evaluate the matrices for many
rotors, or to avoid holding every ``ℓ`` at once, use a [`DCalculator`](@ref) instead,
which allocates once and computes one ``ℓ`` at a time.

See also [`d`](@ref) and [`sYlm`](@ref).
"""
function D(R::Rotor{T}, ℓₘₐₓ::IT; kwargs...) where {T<:Real, IT<:IntegerHalf}
    calc = DCalculator(R, ℓₘₐₓ; kwargs...)
    WignerSeries(
        [copy(recurrence!(calc, ℓ)) for ℓ ∈ ℓₘᵢₙ(IT):ℓₘₐₓ], ℓₘᵢₙ(IT), ℓₘₐₓ
    )
end

"""
    d(β, ℓₘₐₓ; m′ₘₐₓ=ℓₘₐₓ, m′ₘᵢₙ=-m′ₘₐₓ, mₘₐₓ=ℓₘₐₓ, mₘᵢₙ=-mₘₐₓ)

Wigner's real ``d^{(ℓ)}_{m′,m}(β)`` matrices for all ``ℓ ≤ ℓₘₐₓ``, for the single angle `β`,
which may also be given as the phase ``e^{iβ}`` or as a `Rotor` (of which only the ``β`` Euler
angle is used).  The result is indexed as `d(β, ℓₘₐₓ)[ℓ][m′, m]`.  See [`D`](@ref) for
details; this function is the real, ``β``-only analogue.

`ℓₘₐₓ` may be a half-integer `Rational` (denominator 2).  Half-integer ``d`` has period
``4π`` in ``β``, so an angle or a `Rotor` determines it unambiguously, but a bare phase
``e^{iβ}`` determines it only up to the sign ``(-1)^{2ℓ}`` (the branch ``β ∈ (-π, π]`` is
used).
"""
function d(β::Union{Real, Complex, Rotor}, ℓₘₐₓ::IT; kwargs...) where {IT<:IntegerHalf}
    calc = dCalculator(β, ℓₘₐₓ; kwargs...)
    WignerSeries(
        [copy(recurrence!(calc, ℓ)) for ℓ ∈ ℓₘᵢₙ(IT):ℓₘₐₓ], ℓₘᵢₙ(IT), ℓₘₐₓ
    )
end


# ASCII aliases for the public accessors (see also `ell`, `mpmax`, etc. in wigner_matrix.jl)
const ellmax = ℓₘₐₓ
const Nr = Nᵣ


### Boundary conversion.
#
# The recurrences want a `HalfOddInteger`, but `3//2` is what a caller naturally writes, and
# is what every previous version of this package accepted.  These methods convert and
# re-dispatch, so the `Rational` spelling never reaches the hot code.  See
# [`half_integer`](@ref).

function DCalculator(R, ℓₘₐₓ::Rational; kwargs...)
    DCalculator(R, half_integer(ℓₘₐₓ); half_integer_kwargs(kwargs)...)
end
function dCalculator(β, ℓₘₐₓ::Rational; kwargs...)
    dCalculator(β, half_integer(ℓₘₐₓ); half_integer_kwargs(kwargs)...)
end
function D(R::Rotor, ℓₘₐₓ::Rational; kwargs...)
    D(R, half_integer(ℓₘₐₓ); half_integer_kwargs(kwargs)...)
end
function d(β::Union{Real, Complex, Rotor}, ℓₘₐₓ::Rational; kwargs...)
    d(β, half_integer(ℓₘₐₓ); half_integer_kwargs(kwargs)...)
end
