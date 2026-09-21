"""
    WignerHCalculator(β, ℓₘₐₓ; m′ₘₐₓ=ℓₘₐₓ)

Engine for the Gumerov–Duraiswami recurrences, computing the ``H`` wedge (see [`HWedge`](@ref))
for one value of ``ℓ`` at a time, for `Nᵣ` rotors simultaneously.

The ``H`` matrix depends on the rotor only through ``β``, so the calculator stores one phase
``e^{iβ}`` per rotor.  The first argument supplies them: an angle ``β``, a phase ``e^{iβ}``, a
`Rotor` (of which only ``β`` is used), or an `AbstractVector` of any one of those — and it is
the vector case that makes the calculator handle `Nᵣ = length(β)` rotors at once.  Later
values are supplied by [`set_β!`](@ref).  The wedge is stored for ``|m′| ≤ m′ₘₐₓ`` and
``m ≥ |m′|``; every other element of the ``H`` matrix is obtained by symmetry through
[`wedge_value`](@ref).

The element type is the rotor data's own: an angle given as a `Float32` gives a `Float32`
calculator, and there is no argument to override that.  To compute in another type, convert
the data — `WignerHCalculator(Double64(β), ℓₘₐₓ)` — which says what is meant, that these are
the values to treat as exact.  `floattype(calc)` reports the type in use.

This is the low-level engine shared by [`WignerDCalculator`](@ref), [`WignerdCalculator`](@ref)
and the spin-weighted spherical harmonics; most users will want one of those instead.

# Usage

```julia
calc = WignerHCalculator(rotors, ℓₘₐₓ)
for ℓ ∈ 0:ℓₘₐₓ
    recurrence!(calc, ℓ)            # advance to the next ℓ (cheap when sequential)
    H = calc.Hˡ                     # HWedge for the current ℓ; H[iᵣ, m′, m] for m ≥ |m′|
end
```

Unlike [`WignerDCalculator`](@ref) and the other calculators built on it, this one is not
iterable: its wedge is one mutable object handed back by identity, rather than a view that
`copy` can preserve.

Requesting an ``ℓ`` smaller than the current one restarts the recurrence from ``ℓ_{min}``;
jumping forward advances through the intermediate values.  Nothing is allocated after
construction.

# Half-integer indices

Passing a `Rational` `ℓₘₐₓ` with denominator 2 — `WignerHCalculator(β, 7//2)` — gives a
calculator for half-integer ``ℓ, m′, m``.  Then `m′ₘₐₓ` must also be a half-integer, and
`recurrence!(calc, ℓ)` accepts only half-integer `ℓ`.  The recurrence is the same one:
the ``m'=0`` axis is run at the *integer* order ``j = ℓ - 1/2``, the rows ``m' = ±1/2`` are
seeded from it by a Clebsch–Gordan step, and the ``m'`` ladder then proceeds unchanged (see
the notes on the [``H`` recursion](@ref "Algorithm for computing ``H`` (redesigned)")).  Note that the
``H`` matrix is then symmetric only up to the sign
``σ = \\mathrm{sgn}(m)\\,\\mathrm{sgn}(m')`` (see [`transpose_sign`](@ref)), which
[`wedge_value`](@ref) applies for you.

Because half-integer ``d`` has period ``4π`` in ``β``, a rotor or an angle ``β`` determines
it unambiguously, but a bare phase ``e^{iβ}`` determines ``β`` only modulo ``2π`` and hence
``d`` only up to the double-cover sign ``(-1)^{2ℓ}``; the branch ``β ∈ (-π, π]`` is used.
"""
struct WignerHCalculator{IT, RT<:Real, ST}
    # The axes are always *integer*-indexed: for half-integer ℓ they encode the order
    # j = ℓ - 1/2, and `OffsetArray`-like half-integer labels would buy nothing.
    h⃗ᵃ::HAxis{Int, RT}
    h⃗ᵇ::HAxis{Int, RT}
    Hˡ::HWedge{IT, RT, ST}
    eⁱᵝ::FixedSizeVectorDefault{Complex{RT}}
    cβ½::FixedSizeVectorDefault{RT}  # cos(β/2) per rotor; length 0 unless IT <: HalfOddInteger
    sβ½::FixedSizeVectorDefault{RT}  # sin(β/2) per rotor; length 0 unless IT <: HalfOddInteger
    ℓₘₐₓ::IT
    m′ₘₐₓ::IT
    swapH::Base.RefValue{Bool}  # h⃗ˡ(w) returns h⃗ᵃ if `false`, otherwise h⃗ᵇ; and vice versa for h⃗ˡ⁺¹(w)
    axes_valid::Base.RefValue{Bool}  # h⃗ˡ and h⃗ˡ⁺¹ hold correct data for their ℓ labels
end

function WignerHCalculator(β, ℓₘₐₓ::Rational; kwargs...)
    WignerHCalculator(β, half_integer(ℓₘₐₓ); half_integer_kwargs(kwargs)...)
end
function WignerHCalculator(β, ℓₘₐₓ::IT; m′ₘₐₓ::IT=ℓₘₐₓ) where {IT<:HalfInteger}
    # `rotor_basetype` is called in argument position so that the element type reaches
    # `allocate_H` as a type rather than as a value, which is what keeps the result
    # inferrable.
    set_rotors!(allocate_H(IT, rotor_basetype(β), ℓₘₐₓ, m′ₘₐₓ, nrotors(β)), β)
end

# Allocate the buffers without touching them.  PRIVATE, and deliberately so: the returned
# calculator's rotor data (`eⁱᵝ`, and the half angles on the half-integer path) is
# uninitialized, and nothing in the object records that fact.  Every caller must therefore
# either `set_rotors!` it or copy those buffers in from another calculator before the object
# can escape — which is exactly what the public constructor above and `similar` below do.
function allocate_H(
    ::Type{IT}, ::Type{RT}, ℓₘₐₓ::IT, m′ₘₐₓ::IT, Nᵣ::Int
) where {IT<:HalfInteger, RT<:Real}
    if m′ₘₐₓ < 0 || m′ₘₐₓ > ℓₘₐₓ
        error("m′ₘₐₓ=$m′ₘₐₓ must satisfy 0 ≤ m′ₘₐₓ ≤ ℓₘₐₓ=$ℓₘₐₓ.")
    end
    validate_index_ranges(ℓₘₐₓ, m′ₘₐₓ, -m′ₘₐₓ)
    # The axis buffers hold the *integer* orders 0:axisℓₘₐₓ.  For integer ℓ that is
    # 0:ℓₘₐₓ+1, one extra because step 3 reads the ℓ+1 axis; for half-integer ℓ it is
    # 0:jₘₐₓ+1 with jₘₐₓ = ℓₘₐₓ - 1/2, one extra because the axis advance computes it anyway.
    axisℓₘₐₓ = Int(ℓₘₐₓ - ℓₘᵢₙ(IT)) + 1
    h⃗ᵃ = HAxis(RT, Nᵣ, axisℓₘₐₓ)
    h⃗ᵇ = HAxis(RT, Nᵣ, axisℓₘₐₓ)
    h⃗ᵇ.ℓ = 1
    Hˡ = HWedge(RT, Nᵣ, ℓₘₐₓ, m′ₘₐₓ, -m′ₘₐₓ)
    eⁱᵝ = FixedSizeVector{Complex{RT}}(undef, Nᵣ)
    # The half-angle pair (cos(β/2), sin(β/2)) is needed only by the half-integer seed; the
    # integer path gets zero-length buffers, so it allocates and touches nothing extra.
    Nₕ = IT <: HalfOddInteger ? Nᵣ : 0
    cβ½ = FixedSizeVector{RT}(undef, Nₕ)
    sβ½ = FixedSizeVector{RT}(undef, Nₕ)
    WignerHCalculator{IT, RT, typeof(parent(Hˡ))}(
        h⃗ᵃ, h⃗ᵇ, Hˡ, eⁱᵝ, cβ½, sβ½, ℓₘₐₓ, m′ₘₐₓ, Ref(false), Ref(false)
    )
end

# The rotor data is copied buffer-by-buffer rather than by re-running `set_rotors!`, because
# no calculator stores the rotor it was given: `eⁱᵝ` alone would fix `β` only modulo 2π, and
# so would flip the double-cover sign (-1)^{2ℓ} on the half-integer path.  The half-angle
# copies are no-ops on the integer path, where those buffers have length zero.
function Base.similar(w::WignerHCalculator{IT, RT}) where {IT, RT}
    w′ = allocate_H(IT, RT, w.ℓₘₐₓ, w.m′ₘₐₓ, Nᵣ(w))
    copyto!(w′.eⁱᵝ, w.eⁱᵝ)
    copyto!(w′.cβ½, w.cβ½)
    copyto!(w′.sβ½, w.sβ½)
    w′
end
function Base.similar(w::WignerHCalculator{IT, RT}, β) where {IT, RT}
    if nrotors(β) != Nᵣ(w)
        error("This calculator handles Nᵣ=$(Nᵣ(w)) rotors, but got $(nrotors(β)).")
    end
    check_rotor_type(w, β)
    set_rotors!(allocate_H(IT, RT, w.ℓₘₐₓ, w.m′ₘₐₓ, Nᵣ(w)), β)
end

ℓ(w::WignerHCalculator) = Hˡ(w).ℓ
ℓₘᵢₙ(w::WignerHCalculator{IT}) where {IT} = ℓₘᵢₙ(IT)
ℓₘₐₓ(w::WignerHCalculator) = w.ℓₘₐₓ
m′ₘₐₓ(w::WignerHCalculator) = w.m′ₘₐₓ
floattype(::WignerHCalculator{IT, RT}) where {IT, RT} = RT
m′ₘᵢₙ(w::WignerHCalculator) = -w.m′ₘₐₓ
Nᵣ(w::WignerHCalculator) = Nᵣ(Hˡ(w))

h⃗ˡ(w::WignerHCalculator) = w.swapH[] ? w.h⃗ᵇ : w.h⃗ᵃ
h⃗ˡ⁺¹(w::WignerHCalculator) = w.swapH[] ? w.h⃗ᵃ : w.h⃗ᵇ
Hˡ(w::WignerHCalculator) = w.Hˡ
eⁱᵝ(w::WignerHCalculator) = w.eⁱᵝ

"""
    fill!(w::WignerHCalculator, v)

Fill every internal buffer of `w` (both axes and the wedge) with the value `v`, and mark the
axis data as invalid.  The stored rotor data — the phases `e^{iβ}` and, on the half-integer
path, the half angles — is deliberately left alone, so that `recurrence!(w, ℓ)` still has
everything it needs.  Useful for testing that no uninitialized storage is ever read:
everything the recurrence is responsible for writing is poisoned, while everything the rotor
data consists of is preserved.
"""
function Base.fill!(w::WignerHCalculator{IT, RT}, v::Real) where {IT, RT}
    let v = convert(RT, v)
        fill!(parent(w.h⃗ᵃ), v)
        fill!(parent(w.h⃗ᵇ), v)
        fill!(parent(w.Hˡ), v)
    end
    w.axes_valid[] = false
    w
end

function increment_axes!(w::WignerHCalculator)
    # The data that is now stored as h⃗ˡ(w) will get swapped below so that it will be
    # returned by h⃗ˡ⁺¹(w), so we need to increment its ℓ value twice.
    let h⃗ˡ = h⃗ˡ(w)
        h⃗ˡ.ℓ = h⃗ˡ.ℓ + oftype(h⃗ˡ.ℓ, 2)
    end
    # The data that is now stored as h⃗ˡ⁺¹(w) will get swapped below so that it will be
    # returned by h⃗ˡ(w), which will already be correct for the next ℓ value.
    w.swapH[] = !w.swapH[]
    w
end

"""
    axis_ℓ(w::WignerHCalculator, ℓ)

Label of the integer axis that seeds the wedge of order `ℓ`: ``ℓ`` itself for integer
indices, and ``ℓ - 1/2`` for half-integer ones.  The axis buffers are always labelled by
this `Int`, never by `ℓ`.
"""
@inline axis_ℓ(::WignerHCalculator{IT}, ℓ) where {IT} = Int(ℓ - ℓₘᵢₙ(IT))

# Copy the m′ = 0 row of the wedge from the integer axis.  Integer indices only: for
# half-integer ℓ there is no m′ = 0 row, and `recurrence_seed!` writes the rows m′ = ±1/2
# instead.
function fillHˡ₀ₘ!(w::WignerHCalculator{IT}) where {IT<:Integer}
    let h⃗ˡ = h⃗ˡ(w), Hˡ = Hˡ(w)
        if h⃗ˡ.ℓ != axis_ℓ(w, Hˡ.ℓ)
            error("Cannot fill Hˡ₀ₘ for ℓ=$(Hˡ.ℓ) from h⃗ˡ for ℓ=$(h⃗ˡ.ℓ).")
        end
        # Get the index to the start of the central row, m′ = 0 or 1//2
        iˡ₀₀ = row_index(Hˡ, ℓₘᵢₙ(Hˡ))
        # Figure out how many entries to copy
        N = Nᵣ(Hˡ) * (Int(Hˡ.ℓ - ℓₘᵢₙ(Hˡ)) + 1)
        # Now just copy that many entries from h⃗ˡ₀ₘ into row Hˡ₀ₘ
        copyto!(parent(Hˡ), iˡ₀₀, parent(h⃗ˡ), 1, N)
    end
    w
end

function Base.show(io::IO, w::WignerHCalculator{IT, RT, ST}) where {IT, RT, ST}
    print(
        io,
        "WignerHCalculator{$IT, $RT} for ℓₘₐₓ=$(ℓₘₐₓ(w)), ",
        "m′ₘₐₓ=$(m′ₘₐₓ(w)), Nᵣ=$(Nᵣ(w))",
        w.axes_valid[] ? ", currently at ℓ=$(ℓ(w))" : " (nothing computed yet)"
    )
end
Base.show(io::IO, ::MIME"text/plain", w::WignerHCalculator) = show(io, w)


### Rotor data

"""
    spinor_phases(R::AbstractQuaternion, [F])

Return `(eⁱᵝ, z₊, z₋, cβ½, sβ½)` for the rotor `R`, where ``β`` is the Euler angle,
``z₊ = e^{i(α+γ)/2}``, ``z₋ = e^{i(α-γ)/2}``, ``cβ½ = \\cos(β/2)`` and
``sβ½ = \\sin(β/2)``.  These are the quantities the Wigner recurrences need:
``e^{i(m′α + mγ)} = z₊^{m′+m} z₋^{m′-m}``, with integer exponents even for half-integer
``m′, m``, while the half-angle pair seeds the half-integer recurrence.  At the poles
``β ∈ \\{0, π\\}`` the undefined phase is set to 1; the corresponding ``d`` elements vanish,
so the choice is immaterial.

The half-angles are taken as ``(\\sqrt{W²+Z²}, \\sqrt{X²+Y²})/\\|R\\|``, which is accurate
near both poles and is non-negative, so ``β ∈ [0, π]``; a rotor's double-cover sign is
encoded entirely in `z₊` and `z₋`, giving ``𝔇(-R) = -𝔇(R)`` for half-integer indices.

Callers that need only the first few outputs may drop the rest:
`eⁱᵝ, z₊, z₋ = spinor_phases(R, F)`.

`R` need not be normalized.

The optional second argument is the real type the phases are computed in; it defaults to
`float(eltype(R))`.  Pass the *calculator's* type whenever that is more precise than the
rotor's own — otherwise every later step inherits the rotor type's precision.
"""
spinor_phases(R::AbstractQuaternion{T}) where {T} = spinor_phases(R, float(T))
function spinor_phases(R::AbstractQuaternion, ::Type{F}) where {F<:Real}
    a = F(R[1])^2 + F(R[4])^2
    b = F(R[2])^2 + F(R[3])^2
    sqrta = √a
    sqrtb = √b
    z₊ = iszero(sqrta) ? one(Complex{F}) : Complex{F}(F(R[1]), F(R[4])) / sqrta  # exp[i(α+γ)/2]
    z₋ = iszero(sqrtb) ? one(Complex{F}) : Complex{F}(F(R[3]), -F(R[2])) / sqrtb  # exp[i(α-γ)/2]
    eⁱᵝ = Complex{F}(a - b, 2 * sqrta * sqrtb) / (a + b)
    nrm = √(a + b)
    cβ½ = sqrta / nrm
    sβ½ = sqrtb / nrm
    (eⁱᵝ, z₊, z₋, cβ½, sβ½)
end

"""
    half_angles(eⁱᵝ)

The pair ``(\\cos(β/2), \\sin(β/2))`` for the branch ``β ∈ (-π, π]`` of the phase
``e^{iβ}``, computed without cancellation at either pole.

A bare phase fixes ``β`` only modulo ``2π``, so for half-integer indices this fixes ``d``
only up to the double-cover sign ``(-1)^{2ℓ}``; pass the angle ``β`` itself or a `Rotor` if
that matters.
"""
@inline function half_angles(z::Complex{RT}) where {RT<:Real}
    cosβ, sinβ = reim(z)
    if cosβ ≥ 0
        c = √((1 + cosβ) / 2)
        s = sinβ / (2c)
    else
        s = copysign(√((1 - cosβ) / 2), sinβ)
        c = sinβ / (2s)
    end
    (c, s)
end

# Store the half-angle pair for rotor `i`.  Each of these is a no-op — emitting no code at
# all, and not even evaluating its source — for integer index types, whose buffers have
# length zero; that keeps the integer path bit-identical and allocation-free.
@inline set_half_angles!(::WignerHCalculator{IT}, ::Int, ::Real, ::Real) where {IT<:Integer} = nothing
@inline function set_half_angles!(
    w::WignerHCalculator{IT, RT}, i::Int, c::Real, s::Real
) where {IT<:HalfOddInteger, RT}
    @inbounds w.cβ½[i] = convert(RT, c)
    @inbounds w.sβ½[i] = convert(RT, s)
    nothing
end

@inline set_half_angles_from_phase!(::WignerHCalculator{IT}, ::Int, ::Complex) where {IT<:Integer} = nothing
@inline function set_half_angles_from_phase!(
    w::WignerHCalculator{IT}, i::Int, z::Complex
) where {IT<:HalfOddInteger}
    set_half_angles!(w, i, half_angles(z)...)
end

@inline set_half_angles_from_angle!(::WignerHCalculator{IT}, ::Int, ::Real) where {IT<:Integer} = nothing
@inline function set_half_angles_from_angle!(
    w::WignerHCalculator{IT, RT}, i::Int, β::Real
) where {IT<:HalfOddInteger, RT}
    # Straight from β, not from cis(β): this is the only input form that honours the true
    # 4π periodicity of half-integer d.
    let β½ = convert(RT, β) / 2
        set_half_angles!(w, i, cos(β½), sin(β½))
    end
end

function set_rotors!(w::WignerHCalculator{IT, RT}, eⁱᵝ::AbstractVector{<:Complex}) where {IT, RT<:Real}
    if length(eⁱᵝ) != Nᵣ(w)
        error("Expected $(Nᵣ(w)) rotors (Nᵣ), but got $(length(eⁱᵝ)).")
    end
    # The phase must be a unit complex number; anything else is almost certainly a mistake
    # (e.g., passing β itself as a complex number), so we refuse it rather than silently
    # producing garbage.  The tolerance is generous enough for `cis(β)` in any precision, and
    # the value is used as given (not renormalized), so that a phase and the corresponding
    # angle give bitwise-identical results.  Validate everything before mutating anything, so
    # that a rejected input leaves the calculator's state untouched.
    tolerance = 64 * max(eps(RT), eps(float(real(eltype(eⁱᵝ)))))
    for i ∈ eachindex(eⁱᵝ)
        z = convert(Complex{RT}, eⁱᵝ[i])
        if abs(abs2(z) - 1) > tolerance
            error(
                "The phase e^{iβ} must have unit modulus; got $(eⁱᵝ[i]) at index $i, "
                * "with |e^{iβ}|² = $(abs2(z))."
            )
        end
    end
    @inbounds for i ∈ eachindex(eⁱᵝ)
        z = convert(Complex{RT}, eⁱᵝ[i])
        w.eⁱᵝ[i] = z
        set_half_angles_from_phase!(w, i, z)
    end
    w.axes_valid[] = false
    w
end
function set_rotors!(w::WignerHCalculator, R)
    error(
        "Cannot set rotor data of type $(typeof(R)) for a calculator with Nᵣ=$(Nᵣ(w)); "
        * _rotor_input_forms * "."
    )
end
function set_rotors!(w::WignerHCalculator{IT, RT}, β::AbstractVector{<:Real}) where {IT, RT<:Real}
    if length(β) != Nᵣ(w)
        error("Expected $(Nᵣ(w)) rotors (Nᵣ), but got $(length(β)).")
    end
    @inbounds for i ∈ eachindex(β)
        w.eⁱᵝ[i] = cis(convert(RT, β[i]))
        set_half_angles_from_angle!(w, i, β[i])
    end
    w.axes_valid[] = false
    w
end
function set_rotors!(w::WignerHCalculator{IT, RT}, R::AbstractVector{<:Rotor}) where {IT, RT<:Real}
    if length(R) != Nᵣ(w)
        error("Expected $(Nᵣ(w)) rotors (Nᵣ), but got $(length(R)).")
    end
    @inbounds for i ∈ eachindex(R)
        eⁱᵝᵢ, _, _, cβ½, sβ½ = spinor_phases(R[i], RT)
        w.eⁱᵝ[i] = eⁱᵝᵢ
        set_half_angles!(w, i, cβ½, sβ½)
    end
    w.axes_valid[] = false
    w
end
function set_rotors!(w::WignerHCalculator{IT, RT}, R::Union{Real, Complex, Rotor}) where {IT, RT<:Real}
    if Nᵣ(w) != 1
        error("A single rotor was given, but this calculator expects Nᵣ=$(Nᵣ(w)) rotors.")
    end
    set_rotors!(w, @SVector [R])
end


### Driver

"""
    recurrence!(calc, R, ℓ)
    recurrence!(calc, ℓ)

Compute the Wigner quantities for index ``ℓ`` in the calculator `calc`.

In the first form, the rotor data `R` is stored in the calculator first.  For a calculator
with `Nᵣ` rotors, `R` is an `AbstractVector` of length `Nᵣ`; for `Nᵣ = 1` a single element is
also accepted.  The elements may be `Rotor`s, or — for
quantities that depend only on ``β``, such as ``d`` and ``H`` — the angle ``β`` or the phase
``e^{iβ}``.

In the second form, the rotor data from the previous call is reused.  Successive calls with
``ℓ, ℓ+1, ℓ+2, …`` are the cheap path: each costs ``O(N_r ℓ^2)``.  Requesting a smaller ``ℓ``
than the current one restarts the recurrence from ``ℓ_{min}``.

Returns `calc`; the result is then available as `calc[ℓ]` (for [`WignerDCalculator`](@ref) and
[`WignerdCalculator`](@ref)) or `calc.Hˡ` (for [`WignerHCalculator`](@ref)).

A phase given as a `Complex` number must have unit modulus (to within rounding), and is used
as given; an angle given as a `Real` is converted to the calculator's number type.
"""
function recurrence!(w::WignerHCalculator, R, ℓ)
    check_ℓ(w, ℓ)
    set_rotors!(w, R)
    recurrence!(w, ℓ)
end
# Every `recurrence!(calc, R, ℓ)` runs this *before* `set_rotors!`, so that a bad `ℓ` leaves
# the calculator's stored rotor data untouched ("validate everything before mutating
# anything"; see `set_rotors!` below).  It is also what rejects an `ℓ` of the wrong kind — a
# whole number for a half-integer calculator, or `5//3` for either.
function check_ℓ(w::WignerHCalculator{IT}, ℓ) where {IT}
    ℓ = convert(IT, ℓ)
    if ℓ < ℓₘᵢₙ(w) || ℓ > ℓₘₐₓ(w)
        error(
            "Requested ℓ=$(ℓ) is out of bounds [$(ℓₘᵢₙ(w)), $(ℓₘₐₓ(w))] for this calculator."
        )
    end
end

# Advance (or restart) the integer axis buffers so that h⃗ˡ holds order `j` and h⃗ˡ⁺¹ holds
# `j+1`.  Shared by the integer and half-integer drivers; `j` is `ℓ` in the former case and
# `ℓ - 1/2` in the latter.
function advance_axes!(w::WignerHCalculator, j::Int)
    if !w.axes_valid[] || h⃗ˡ(w).ℓ > j
        h⃗ˡ(w).ℓ = 0
        h⃗ˡ⁺¹(w).ℓ = 1
        recurrence_step1!(w)  # h⃗⁰₀₀ = 1
        recurrence_step2!(w)  # h⃗⁰₀ₘ -> h⃗¹₀ₘ
        w.axes_valid[] = true
    end
    while h⃗ˡ(w).ℓ < j
        increment_axes!(w)
        recurrence_step2!(w)  # h⃗ʲ₀ₘ -> h⃗ʲ⁺¹₀ₘ
    end
    w
end

function recurrence!(w::WignerHCalculator{IT, RT}, ℓ) where {IT<:Signed, RT}
    ℓ = convert(IT, ℓ)
    check_ℓ(w, ℓ)  # rejects an `ℓ` that is not of this calculator's index type

    # Steps 1 and 2 (the ℓ recurrence along the m′=0 axis).  The axis buffers h⃗ˡ and h⃗ˡ⁺¹
    # hold the axes for two successive ℓ values.  We restart from ℓₘᵢₙ if they are invalid or
    # ahead of the requested ℓ, and otherwise advance them one ℓ at a time.
    advance_axes!(w, axis_ℓ(w, ℓ))

    # Steps 3, 4, and 5 (the m′ recurrences at fixed ℓ), filling the wedge Hˡ.
    Hˡ(w).ℓ = ℓ
    fillHˡ₀ₘ!(w)  # Copy h⃗ˡ₀ₘ to Hˡ₀ₘ
    recurrence_step3!(w)  # Hˡ⁺¹₀ₘ -> Hˡ₁ₘ
    recurrence_step4!(w)  # Hˡₘ′ₘ₋₁, Hˡₘ′₋₁ₘ, Hˡₘ′ₘ₊₁ -> Hˡₘ′₊₁ₘ
    recurrence_step5!(w)  # Hˡₘ′ₘ₋₁, Hˡₘ′₊₁ₘ, Hˡₘ′ₘ₊₁ -> Hˡₘ′₋₁ₘ
    # Step 6 (the symmetries) is never applied to the wedge itself; elements outside the
    # wedge are read through `wedge_value`/`wedge_source` when the results are materialized.
    w
end

# The half-integer driver.  The only differences from the integer one are that the axis is
# run at the integer order j = ℓ - 1/2, and that steps "0" (`fillHˡ₀ₘ!`) and 3 — which
# together produce the rows m′ = 0 and m′ = 1 from the ℓ and ℓ+1 axes — are replaced by
# `recurrence_seed!`, which produces the rows m′ = ±1/2 from the j axis alone.  Steps 4 and
# 5 are then literally the same code (see the v3 design memo, §5).
function recurrence!(w::WignerHCalculator{IT, RT}, ℓ) where {IT<:HalfOddInteger, RT}
    ℓ = convert(IT, ℓ)
    check_ℓ(w, ℓ)  # rejects an `ℓ` that is not of this calculator's index type
    advance_axes!(w, axis_ℓ(w, ℓ))
    Hˡ(w).ℓ = ℓ
    recurrence_seed!(w)   # h⃗ʲ₀ₘ -> Hˡ₊₁⁄₂ₘ, Hˡ₋₁⁄₂ₘ
    recurrence_step4!(w)  # Hˡₘ′ₘ₋₁, Hˡₘ′₋₁ₘ, Hˡₘ′ₘ₊₁ -> Hˡₘ′₊₁ₘ
    recurrence_step5!(w)  # Hˡₘ′ₘ₋₁, Hˡₘ′₊₁ₘ, Hˡₘ′ₘ₊₁ -> Hˡₘ′₋₁ₘ
    w
end


### Recurrence steps.
#
# The comments give the equivalent expressions in terms of matrix indices [iᵣ, m′, m]; the
# code uses precomputed linear offsets so that the innermost loop over rotors vectorizes.
#
# Those innermost loops use `@simd ivdep`, which asserts two things the compiler cannot
# prove for itself and which hold at every site below: no iteration depends on a value
# written by an earlier one, and the row being written does not alias any row being read.
# Both follow from the wedge layout — each step writes a row (`m′±1`, or the seed's `±1/2`)
# that is structurally distinct from the rows it reads — so `@simd ivdep` would become
# undefined behaviour if that layout ever changed.  None of these loops is a reduction, so no
# arithmetic is reassociated and the results are unchanged; the only effect is vectorization.
# Measured on `recurrence_step5!`, the hottest of them, this is 1.5× at Nᵣ=8 with
# bit-identical output.  (It costs a few percent at Nᵣ=1, where there is nothing to vectorize.)

# The axis buffers are integer-indexed for every index type (for half-integer ℓ they encode
# the order j = ℓ - 1/2), so steps 1 and 2 are shared verbatim and never see a `Rational`.
function recurrence_step1!(w::WignerHCalculator{IT}) where {IT}
    let h⃗⁰ = h⃗ˡ(w)
        if h⃗⁰.ℓ ≠ ℓₘᵢₙ(h⃗⁰)
            error("recurrence_step1! can only be called for ℓ=$(ℓₘᵢₙ(h⃗⁰)); current ℓ=$(h⃗⁰.ℓ).")
        end
        @inbounds for i ∈ 1:Nᵣ(h⃗⁰)
            h⃗⁰[i] = 1  # h⃗⁰[i, 0, 0] = 1
        end
    end
    w
end

# Compute h⃗ⁿ₀ₘ = Hⁿ₀ₘ for m ∈ 0:n from h⃗ⁿ⁻¹₀ₘ = Hⁿ⁻¹₀ₘ, where n = ℓ+1 is the ℓ value of
# h⃗ˡ⁺¹.  This is the recurrence of Xing et al. (2020) for the normalized associated Legendre
# functions, in the notation of Gumerov and Duraiswami's step 2.
function recurrence_step2!(w::WignerHCalculator{IT, RT}) where {IT, RT}
    let h⃗ⁿ⁻¹ = h⃗ˡ(w), h⃗ⁿ = h⃗ˡ⁺¹(w), eⁱᵝ = eⁱᵝ(w)
        n = h⃗ⁿ⁻¹.ℓ + 1
        if h⃗ⁿ.ℓ ≠ n
            error("Inconsistent axes in recurrence_step2!: ℓ(h⃗ˡ)=$(h⃗ⁿ⁻¹.ℓ), ℓ(h⃗ˡ⁺¹)=$(h⃗ⁿ.ℓ).")
        end

        # Note that in this step only, we use notation derived from (but not the same as)
        # Xing et al., denoting the coefficients as b̄ₙ, c̄ₙₘ, d̄ₙₘ, ēₙₘ.  In the following
        # steps, we will use notation from Gumerov and Duraiswami, who denote their
        # different coefficients aₗᵐ, etc.
        @inbounds let √=sqrt∘RT, Nᵣ = Nᵣ(h⃗ⁿ⁻¹)
            if n == 1
                # We know that h⃗⁰₀₀ = 1, so this is just the general branch with the
                # (nonexistent) h⃗⁰₀₁ set to zero.
                invsqrt2 = inv(√2)
                @simd ivdep for i ∈ 1:Nᵣ
                    cosβ, sinβ = reim(eⁱᵝ[i])
                    h⃗ⁿ[i] = cosβ  # h⃗¹[i, 0, 0] = cosβ
                    h⃗ⁿ[Nᵣ + i] = invsqrt2 * sinβ  # h⃗¹[i, 0, 1] = sinβ / √2
                end
            else
                b̄ₙ = √(RT(n-1)/n)
                @simd ivdep for i ∈ 1:Nᵣ
                    cosβ, sinβ = reim(eⁱᵝ[i])
                    # h⃗ⁿ[i, 0, 0] = cosβ * h⃗ⁿ⁻¹[i, 0, 0] - b̄ₙ * sinβ * h⃗ⁿ⁻¹[i, 0, 1]
                    h⃗ⁿ[i] = cosβ * h⃗ⁿ⁻¹[i] - b̄ₙ * sinβ * h⃗ⁿ⁻¹[Nᵣ + i]
                end
                for m ∈ 1:n-2
                    c̄ₙₘ = √((n+m)*(n-m)) / n
                    d̄ₙₘ = √((n-m)*(n-m-1)) / 2n
                    ēₙₘ = √((n+m)*(n+m-1)) / 2n

                    i⁰ᵐ = Nᵣ * m
                    i⁰ᵐ⁺¹ = Nᵣ * (m + 1)
                    i⁰ᵐ⁻¹ = Nᵣ * (m - 1)

                    @simd ivdep for i ∈ 1:Nᵣ
                        cosβ, sinβ = reim(eⁱᵝ[i])
                        # h⃗ⁿ[i, 0, m] = (
                        #     c̄ₙₘ * cosβ * h⃗ⁿ⁻¹[i, 0, m]
                        #     - sinβ * (d̄ₙₘ * h⃗ⁿ⁻¹[i, 0, m+1] - ēₙₘ * h⃗ⁿ⁻¹[i, 0, m-1])
                        # )
                        h⃗ⁿ[i⁰ᵐ + i] = (
                            c̄ₙₘ * cosβ * h⃗ⁿ⁻¹[i⁰ᵐ + i]
                            - sinβ * (d̄ₙₘ * h⃗ⁿ⁻¹[i⁰ᵐ⁺¹ + i] - ēₙₘ * h⃗ⁿ⁻¹[i⁰ᵐ⁻¹ + i])
                        )
                    end
                end
                let m = n-1
                    # As above, but h⃗ⁿ⁻¹[i, 0, m+1] does not exist (it is zero).
                    c̄ₙₘ = √((n+m)*(n-m)) / n
                    ēₙₘ = √((n+m)*(n+m-1)) / 2n

                    i⁰ᵐ = Nᵣ * m
                    i⁰ᵐ⁻¹ = Nᵣ * (m - 1)

                    @simd ivdep for i ∈ 1:Nᵣ
                        cosβ, sinβ = reim(eⁱᵝ[i])
                        # h⃗ⁿ[i, 0, m] = c̄ₙₘ * cosβ * h⃗ⁿ⁻¹[i, 0, m] + sinβ * ēₙₘ * h⃗ⁿ⁻¹[i, 0, m-1]
                        h⃗ⁿ[i⁰ᵐ + i] = (
                            c̄ₙₘ * cosβ * h⃗ⁿ⁻¹[i⁰ᵐ + i]
                            + sinβ * ēₙₘ * h⃗ⁿ⁻¹[i⁰ᵐ⁻¹ + i]
                        )
                    end
                end
                let m = n
                    # As above, but now h⃗ⁿ⁻¹[i, 0, m] does not exist either.
                    ēₙₘ = √((n+m)*(n+m-1)) / 2n

                    i⁰ᵐ = Nᵣ * m
                    i⁰ᵐ⁻¹ = Nᵣ * (m - 1)

                    @simd ivdep for i ∈ 1:Nᵣ
                        cosβ, sinβ = reim(eⁱᵝ[i])
                        # h⃗ⁿ[i, 0, m] = sinβ * ēₙₘ * h⃗ⁿ⁻¹[i, 0, m-1]
                        h⃗ⁿ[i⁰ᵐ + i] = sinβ * ēₙₘ * h⃗ⁿ⁻¹[i⁰ᵐ⁻¹ + i]
                    end
                end
            end
        end
    end
    w
end

# Compute Hᴶ_{±1/2, m} for m ∈ 1/2:J from the integer axis h⃗ʲ₀ₖ = d^j_{0,k}(β), j = J - 1/2.
# This is the half-integer replacement for `fillHˡ₀ₘ!` together with step 3: it seeds the two
# rows the m′ ladder needs, at a cost of O(Nᵣ J), from Varshalovich Eqs. 4.8.2(14) and (15),
#
#     H^J_{+1/2, m} = [ √(J+m) c h_{m-1/2} - √(J-m) s h_{m+1/2} ] / √(J + 1/2),
#     H^J_{-1/2, m} = [ √(J+m) s h_{m-1/2} + √(J-m) c h_{m+1/2} ] / √(J + 1/2),
#
# with c = cos(β/2), s = sin(β/2).  Both rows are mandatory: the corner H_{-1/2,1/2} cannot
# be reached from the +1/2 row without leaving the wedge.  (See the v3 design memo, §5.3.)
function recurrence_seed!(w::WignerHCalculator{IT, RT}) where {IT<:HalfOddInteger, RT}
    let Hˡ = Hˡ(w), h⃗ʲ = h⃗ˡ(w), cβ½ = w.cβ½, sβ½ = w.sβ½
        @inbounds let √=sqrt∘RT, Nᵣ=Nᵣ(Hˡ), Hp=parent(Hˡ), hp=parent(h⃗ʲ)
            J = Hˡ.ℓ
            half = ℓₘᵢₙ(IT)  # the index 1/2, as a `HalfOddInteger`
            j = J - half     # an `Int`: the order of the integer axis
            if h⃗ʲ.ℓ ≠ j
                error("Inconsistent axis in recurrence_seed!: ℓ(h⃗ˡ)=$(h⃗ʲ.ℓ), j=$j.")
            end
            m′ₘᵢₙw = m′ₘᵢₙ(Hˡ)
            r₊ = row_index(Hˡ)[(half - m′ₘᵢₙw) + 1] - 1    # row m′ = +1/2
            r₋ = row_index(Hˡ)[(-half - m′ₘᵢₙw) + 1] - 1   # row m′ = -1/2
            invnrm = inv(√(J + half))                      # 1 / √(J + 1/2)

            for m ∈ half:(J-1)
                a = √(J + m)              # √(J + m)
                b = √(J - m)              # √(J - m)
                col = Nᵣ * (m - half)        # column of m in rows m′ = ±1/2
                i₊ = r₊ + col
                i₋ = r₋ + col
                iₗ = Nᵣ * (m - half)         # h_{m-1/2}, axis index k = m - 1/2
                iᵤ = Nᵣ * (m + half)         # h_{m+1/2}, axis index k = m + 1/2

                @simd ivdep for i ∈ 1:Nᵣ
                    c = cβ½[i]
                    s = sβ½[i]
                    hₗ = hp[iₗ + i]
                    hᵤ = hp[iᵤ + i]
                    # Hˡ[i, 1//2, m], Hˡ[i, -1//2, m]
                    Hp[i₊ + i] = (a * c * hₗ - b * s * hᵤ) * invnrm
                    Hp[i₋ + i] = (a * s * hₗ + b * c * hᵤ) * invnrm
                end
            end

            # The m = J case is peeled, because h_{j+1} does not exist: its slot in the axis
            # allocation holds another order's data — or, under `fill!(calc, NaN)`, a NaN
            # that `b == 0` would not annihilate.
            let m = J
                a = √(J + m)
                col = Nᵣ * (m - half)
                i₊ = r₊ + col
                i₋ = r₋ + col
                iₗ = Nᵣ * (m - half)

                @simd ivdep for i ∈ 1:Nᵣ
                    hₗ = hp[iₗ + i]
                    Hp[i₊ + i] = a * cβ½[i] * hₗ * invnrm
                    Hp[i₋ + i] = a * sβ½[i] * hₗ * invnrm
                end
            end
        end
    end
    w
end

# Compute Hˡ₁ₘ for m ∈ 1:ℓ from Hˡ⁺¹₀ₘ (Gumerov and Duraiswami's step 3).  Integer indices
# only; the half-integer route uses `recurrence_seed!` instead and never calls this.
function recurrence_step3!(w::WignerHCalculator{IT, RT}) where {IT<:Signed, RT}
    let Hˡ = Hˡ(w), h⃗ˡ⁺¹ = h⃗ˡ⁺¹(w), eⁱᵝ = eⁱᵝ(w)
        @inbounds let √=sqrt∘RT, ℓ=Hˡ.ℓ, Nᵣ = Nᵣ(Hˡ), m′ₘₐₓ=m′ₘₐₓ(Hˡ)
            if h⃗ˡ⁺¹.ℓ ≠ ℓ + 1
                error("Inconsistent axes in recurrence_step3!: ℓ(Hˡ)=$(ℓ), ℓ(h⃗ˡ⁺¹)=$(h⃗ˡ⁺¹.ℓ).")
            end
            if ℓ > 0 && m′ₘₐₓ ≥ 1
                c = 1 / √(ℓ*(ℓ+1))

                # Precompute base offset for m′=1 row in Hˡ
                r¹ = row_index(Hˡ, oneunit(IT)) - 1  # step 3 is integer-only, so `m′ = 1` exists

                for m ∈ 1:ℓ
                    āₗᵐ = √((ℓ+m+1)*(ℓ-m+1))
                    b̄ₗ₊₁ᵐ⁻¹ = √((ℓ-m+1)*(ℓ-m+2))
                    b̄ₗ₊₁⁻ᵐ⁻¹ = √((ℓ+m+1)*(ℓ+m+2))

                    # Column offsets in Hˡ row 1 and h⃗ˡ⁺¹ row 0
                    # `m` inherits the calculator's index type, which need not be `Int`;
                    # the axis labels are always `Int`, so convert here rather than handing
                    # an `Int128` or `BigInt` to `getindex(::HAxis, ...)`.
                    c¹ᵐ = Nᵣ * Int(m - 1)  # Hˡ[i, 1, m] has m′=1, so column is m-abs(1)=m-1
                    i¹ᵐ = r¹ + c¹ᵐ
                    i⁰ᵐ⁺¹ = Nᵣ * Int(m + 1)
                    i⁰ᵐ⁻¹ = Nᵣ * Int(m - 1)
                    i⁰ᵐ = Nᵣ * Int(m)

                    @simd ivdep for i ∈ 1:Nᵣ
                        cosβ, sinβ = reim(eⁱᵝ[i])
                        # Hˡ[i, 1, m] = -c * (
                        #     b̄ₗ₊₁⁻ᵐ⁻¹ * (1 - cosβ) / 2 * h⃗ˡ⁺¹[i, 0, m+1]
                        #     + b̄ₗ₊₁ᵐ⁻¹ * (1 + cosβ) / 2 * h⃗ˡ⁺¹[i, 0, m-1]
                        #     + āₗᵐ * sinβ * h⃗ˡ⁺¹[i, 0, m]
                        # )
                        Hˡ[i¹ᵐ + i] = -c * (
                            b̄ₗ₊₁⁻ᵐ⁻¹ * (1 - cosβ) / 2 * h⃗ˡ⁺¹[i⁰ᵐ⁺¹ + i]
                            + b̄ₗ₊₁ᵐ⁻¹ * (1 + cosβ) / 2 * h⃗ˡ⁺¹[i⁰ᵐ⁻¹ + i]
                            + āₗᵐ * sinβ * h⃗ˡ⁺¹[i⁰ᵐ + i]
                        )
                    end
                end
            end
        end
    end
    w
end

# Compute Hˡₘ′₊₁,ₘ for m′ ∈ (1-ℓₘᵢₙ):m′ₘₐₓ-1 and m ∈ m′+1:ℓ from the rows m′-1 and m′
# (Gumerov and Duraiswami's step 4).  The loop runs on twice-indices, so it is the same code
# for integer and half-integer ℓ; only the starting m′ differs (1 or 1/2), and it is a
# compile-time constant for each index type.
function recurrence_step4!(w::WignerHCalculator{IT, RT}) where {IT, RT}
    let Hˡ = Hˡ(w)
        @inbounds let √=sqrt∘RT, Nᵣ=Nᵣ(Hˡ), ri=row_index(Hˡ)
            ℓ = Hˡ.ℓ
            m′ₘₐₓw = m′ₘₐₓ(Hˡ)
            m′ₘᵢₙw = m′ₘᵢₙ(Hˡ)
            for m′ ∈ (1 - ℓₘᵢₙ(IT)):(m′ₘₐₓw - 1)
                # The m-side signs sgn(m) and sgn(m-1) are +1 throughout the range visited
                # here (m ≥ m′+1 ≥ 3/2 > 0), so they are left out.  The m′-side sign is
                # *not* always +1: at m′ = 1/2 the coefficient of Hˡ[m′-1, m] picks up
                # sgn(-1/2) = -1.  (See the v3 design memo, §5.2.)
                d̄ₗᵐ′ = √(δ²(ℓ, m′))
                d̄ₗᵐ′⁻¹ = sgn(m′ - 1) * √(δ²(ℓ, m′ - 1))
                inv_d̄ₗᵐ′ = inv(d̄ₗᵐ′)

                # Precompute base offsets for m′-1, m′, m′+1 rows.  Note that we subtract 1
                # here because row_index points to the beginning of the desired row, but we
                # just want the offset.
                rᵐ′⁻¹ = ri[((m′ - 1) - m′ₘᵢₙw) + 1] - 1
                rᵐ′ = ri[(m′ - m′ₘᵢₙw) + 1] - 1
                rᵐ′⁺¹ = ri[((m′ + 1) - m′ₘᵢₙw) + 1] - 1

                for m ∈ (m′+1):(ℓ-1)
                    d̄ₗᵐ⁻¹ = √(δ²(ℓ, m - 1))
                    d̄ₗᵐ = √(δ²(ℓ, m))

                    # Compute column offsets within each row.  For row m′, column m has
                    # offset Nᵣ * (m - abs(m′)).
                    cᵐ′⁻¹ᵐ = Nᵣ * (m - abs(m′ - 1))
                    cᵐ′ᵐ⁻¹ = Nᵣ * ((m - 1) - abs(m′))
                    cᵐ′ᵐ⁺¹ = Nᵣ * ((m + 1) - abs(m′))
                    cᵐ′⁺¹ᵐ = Nᵣ * (m - abs(m′ + 1))

                    # Final 1D index offsets (0-based for the loop)
                    iᵐ′⁻¹ᵐ = rᵐ′⁻¹ + cᵐ′⁻¹ᵐ
                    iᵐ′ᵐ⁻¹ = rᵐ′ + cᵐ′ᵐ⁻¹
                    iᵐ′ᵐ⁺¹ = rᵐ′ + cᵐ′ᵐ⁺¹
                    iᵐ′⁺¹ᵐ = rᵐ′⁺¹ + cᵐ′⁺¹ᵐ

                    @simd ivdep for i ∈ 1:Nᵣ
                        # Hˡ[i, m′+1, m] = (
                        #     d̄ₗᵐ′⁻¹ * Hˡ[i, m′-1, m]
                        #     - d̄ₗᵐ⁻¹ * Hˡ[i, m′, m-1]
                        #     + d̄ₗᵐ * Hˡ[i, m′, m+1]
                        # ) / d̄ₗᵐ′
                        Hˡ[iᵐ′⁺¹ᵐ + i] = (
                            d̄ₗᵐ′⁻¹ * Hˡ[iᵐ′⁻¹ᵐ + i]
                            - d̄ₗᵐ⁻¹ * Hˡ[iᵐ′ᵐ⁻¹ + i]
                            + d̄ₗᵐ * Hˡ[iᵐ′ᵐ⁺¹ + i]
                        ) * inv_d̄ₗᵐ′
                    end
                end

                # Now, we do the m=ℓ case separately, since there is no m+1, so we would get
                # out-of-bounds accesses; we just copy the body of the loop above, but
                # remove anything that involves m+1.
                let m = ℓ
                    d̄ₗᵐ⁻¹ = √(δ²(ℓ, m - 1))

                    cᵐ′⁻¹ᵐ = Nᵣ * (m - abs(m′ - 1))
                    cᵐ′ᵐ⁻¹ = Nᵣ * ((m - 1) - abs(m′))
                    cᵐ′⁺¹ᵐ = Nᵣ * (m - abs(m′ + 1))

                    iᵐ′⁻¹ᵐ = rᵐ′⁻¹ + cᵐ′⁻¹ᵐ
                    iᵐ′ᵐ⁻¹ = rᵐ′ + cᵐ′ᵐ⁻¹
                    iᵐ′⁺¹ᵐ = rᵐ′⁺¹ + cᵐ′⁺¹ᵐ

                    @simd ivdep for i ∈ 1:Nᵣ
                        # Hˡ[i, m′+1, m] = (
                        #     d̄ₗᵐ′⁻¹ * Hˡ[i, m′-1, m]
                        #     - d̄ₗᵐ⁻¹ * Hˡ[i, m′, m-1]
                        # ) / d̄ₗᵐ′
                        Hˡ[iᵐ′⁺¹ᵐ + i] = (
                            d̄ₗᵐ′⁻¹ * Hˡ[iᵐ′⁻¹ᵐ + i]
                            - d̄ₗᵐ⁻¹ * Hˡ[iᵐ′ᵐ⁻¹ + i]
                        ) * inv_d̄ₗᵐ′
                    end
                end
            end
        end
    end
    w
end

# Compute Hˡₘ′₋₁,ₘ for m′ ∈ -ℓₘᵢₙ:-1:m′ₘᵢₙ+1 and m ∈ -m′+1:ℓ from the rows m′ and m′+1
# (Gumerov and Duraiswami's step 5).  As in step 4, the loop runs on twice-indices and is
# shared between index types; only the starting m′ differs (0 or -1/2).
function recurrence_step5!(w::WignerHCalculator{IT, RT}) where {IT, RT}
    let Hˡ = Hˡ(w)
        @inbounds let √=sqrt∘RT, Nᵣ=Nᵣ(Hˡ), ri=row_index(Hˡ)
            ℓ = Hˡ.ℓ
            m′ₘᵢₙw = m′ₘᵢₙ(Hˡ)
            for m′ ∈ (-ℓₘᵢₙ(IT)):-1:(m′ₘᵢₙw + 1)
                d̄ₗᵐ′ = sgn(m′) * √(δ²(ℓ, m′))
                d̄ₗᵐ′⁻¹ = sgn(m′ - 1) * √(δ²(ℓ, m′ - 1))
                inv_d̄ₗᵐ′⁻¹ = inv(d̄ₗᵐ′⁻¹)

                # Precompute base offsets for m′-1, m′, m′+1 rows
                rᵐ′⁻¹ = ri[((m′ - 1) - m′ₘᵢₙw) + 1] - 1
                rᵐ′ = ri[(m′ - m′ₘᵢₙw) + 1] - 1
                rᵐ′⁺¹ = ri[((m′ + 1) - m′ₘᵢₙw) + 1] - 1

                for m ∈ (-m′+1):(ℓ-1)
                    d̄ₗᵐ = sgn(m) * √(δ²(ℓ, m))
                    d̄ₗᵐ⁻¹ = sgn(m - 1) * √(δ²(ℓ, m - 1))

                    # Compute column offsets within each row
                    cᵐ′⁺¹ᵐ = Nᵣ * (m - abs(m′ + 1))
                    cᵐ′ᵐ⁻¹ = Nᵣ * ((m - 1) - abs(m′))
                    cᵐ′ᵐ⁺¹ = Nᵣ * ((m + 1) - abs(m′))
                    cᵐ′⁻¹ᵐ = Nᵣ * (m - abs(m′ - 1))

                    iᵐ′⁺¹ᵐ = rᵐ′⁺¹ + cᵐ′⁺¹ᵐ
                    iᵐ′ᵐ⁻¹ = rᵐ′ + cᵐ′ᵐ⁻¹
                    iᵐ′ᵐ⁺¹ = rᵐ′ + cᵐ′ᵐ⁺¹
                    iᵐ′⁻¹ᵐ = rᵐ′⁻¹ + cᵐ′⁻¹ᵐ

                    @simd ivdep for i ∈ 1:Nᵣ
                        # Hˡ[i, m′-1, m] = (
                        #     d̄ₗᵐ′ * Hˡ[i, m′+1, m]
                        #     + d̄ₗᵐ⁻¹ * Hˡ[i, m′, m-1]
                        #     - d̄ₗᵐ * Hˡ[i, m′, m+1]
                        # ) / d̄ₗᵐ′⁻¹
                        Hˡ[iᵐ′⁻¹ᵐ + i] = (
                            d̄ₗᵐ′ * Hˡ[iᵐ′⁺¹ᵐ + i]
                            + d̄ₗᵐ⁻¹ * Hˡ[iᵐ′ᵐ⁻¹ + i]
                            - d̄ₗᵐ * Hˡ[iᵐ′ᵐ⁺¹ + i]
                        ) * inv_d̄ₗᵐ′⁻¹
                    end
                end
                let m = ℓ
                    d̄ₗᵐ⁻¹ = sgn(m - 1) * √(δ²(ℓ, m - 1))

                    cᵐ′⁺¹ᵐ = Nᵣ * (m - abs(m′ + 1))
                    cᵐ′ᵐ⁻¹ = Nᵣ * ((m - 1) - abs(m′))
                    cᵐ′⁻¹ᵐ = Nᵣ * (m - abs(m′ - 1))

                    iᵐ′⁺¹ᵐ = rᵐ′⁺¹ + cᵐ′⁺¹ᵐ
                    iᵐ′ᵐ⁻¹ = rᵐ′ + cᵐ′ᵐ⁻¹
                    iᵐ′⁻¹ᵐ = rᵐ′⁻¹ + cᵐ′⁻¹ᵐ

                    @simd ivdep for i ∈ 1:Nᵣ
                        # Hˡ[i, m′-1, m] = (
                        #     d̄ₗᵐ′ * Hˡ[i, m′+1, m]
                        #     + d̄ₗᵐ⁻¹ * Hˡ[i, m′, m-1]
                        # ) / d̄ₗᵐ′⁻¹
                        Hˡ[iᵐ′⁻¹ᵐ + i] = (
                            d̄ₗᵐ′ * Hˡ[iᵐ′⁺¹ᵐ + i]
                            + d̄ₗᵐ⁻¹ * Hˡ[iᵐ′ᵐ⁻¹ + i]
                        ) * inv_d̄ₗᵐ′⁻¹
                    end
                end
            end
        end
    end
    w
end
