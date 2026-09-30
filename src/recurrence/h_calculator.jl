# The supertype of the calculators: `HCalculator`, and `WignerCalculator` and
# `HarmonicCalculator`, which are built on one.  `IT` is the index type, `Int` or
# `HalfOddInteger`.
abstract type AbstractCalculator{IT} end

ℓₘᵢₙ(::AbstractCalculator{IT}) where {IT} = lowest_index(IT)

"""
    HCalculator(β, ℓₘₐₓ; m′ₘₐₓ=ℓₘₐₓ)

Engine for the Gumerov–Duraiswami recurrences, computing the ``H`` wedge (see [`HWedge`](@ref))
for one value of ``ℓ`` at a time, for `Nᵣ` rotors simultaneously.
    
The parameters of the type `HCalculator{IT, RT, ST}` are as follows:
- `IT` is the index type, `Int` or [`HalfOddInteger`](@ref).
- `RT` is the real type the calculator works in.
- `ST` is the storage type of the ``H`` wedge.

The ``H`` matrix depends on the rotor only through ``β``, so the calculator stores one phase
``e^{iβ}`` per rotor.  The first argument supplies them: an angle ``β``, a phase ``e^{iβ}``,
a `Rotor` or `Quaternion`, or an `AbstractVector` of any one of those — and it is the vector
case that makes the calculator handle `Nᵣ = length(β)` rotors at once.  A `Rotor`
contributes the ``β ∈ [0, π]`` of its canonical Euler decomposition; if the rotor was built
from a ``β`` outside that range, that ``β`` is folded into its ``α`` and ``γ``, which this
calculator does not see.  Later values are supplied by [`set_β!`](@ref).  The wedge is
stored for ``|m′| ≤ m′ₘₐₓ`` and ``m ≥ |m′|``; every other element of the ``H`` matrix is
obtained by symmetry through [`wedge_value`](@ref).  The keyword may also be spelled
`mp_max`.

The element type is the rotor data's own: an angle given as a `Float32` gives a `Float32`
calculator, and there is no argument to override that.  To compute in another type, convert
the data — `HCalculator(Double64(β), ℓₘₐₓ)` — which says what is meant, that these are
the values to treat as exact.  `floattype(calc)` reports the type in use.

This is the low-level calculator of ``H`` shared by [`DCalculator`](@ref),
[`dCalculator`](@ref), and the spin-weighted spherical harmonics; most users will want one
of those instead.

# Usage

```julia
calc = HCalculator(rotors, ℓₘₐₓ)
for ℓ ∈ 0:ℓₘₐₓ
    recurrence!(calc, ℓ)            # advance to the next ℓ (cheap when sequential)
    H = calc.Hˡ                     # HWedge for the current ℓ; H[iᵣ, m′, m] for m ≥ |m′|
end
```

Unlike [`DCalculator`](@ref) and the other calculators built on it, this one is not
iterable: its wedge is one mutable object handed back by identity, rather than a view that
`copy` can preserve.  The wedge belongs to the calculator, so its `ℓ` must not be
reassigned; `copy(calc.Hˡ)` gives an independent wedge that survives the next step.

Requesting an ``ℓ`` smaller than the current one restarts the recurrence from ``ℓ_{min}``;
jumping forward advances through the intermediate values.  Nothing is allocated after
construction.

A calculator is a mutable workspace, which every step and every setter overwrites, so one
calculator must not be used by two tasks at once; `similar(calc)` gives each task a
calculator of its own.

# Half-integer indices

Passing a `Rational` `ℓₘₐₓ` with denominator 2 — `HCalculator(β, 7//2)` — or a
[`HalfOddInteger`](@ref) gives a calculator for half-integer ``ℓ, m′, m``.  Then `m′ₘₐₓ` must
also be a half-integer, in either spelling, and `recurrence!(calc, ℓ)` accepts only
half-integer `ℓ`.  The recurrence is the same one:
the ``m'=0`` axis is run at the *integer* order ``j = ℓ - 1/2``, the rows ``m' = ±1/2`` are
seeded from it by a Clebsch–Gordan step, and the ``m'`` ladder then proceeds unchanged (see
the notes on the [``H`` recursion](@ref "Algorithm for computing ``H``")).  Note that the
``H`` matrix is then symmetric only up to the sign
``σ = \\mathrm{sgn}(m)\\,\\mathrm{sgn}(m')`` (see [`transpose_sign`](@ref)), which
[`wedge_value`](@ref) applies for you.

Because half-integer ``d`` has period ``4π`` in ``β``, a rotor or an angle ``β`` determines
it unambiguously, but a bare phase ``e^{iβ}`` determines ``β`` only modulo ``2π`` and hence
``d`` only up to the double-cover sign ``(-1)^{2ℓ}``; the branch ``β ∈ (-π, π]`` is used.
"""
struct HCalculator{IT, RT<:Real, ST} <: AbstractCalculator{IT}
    # The axes are always *integer*-indexed: for half-integer ℓ they encode the order
    # j = ℓ - 1/2, and `OffsetArray`-like half-integer labels would buy nothing.
    h⃗ᵃ::HAxis{RT}
    h⃗ᵇ::HAxis{RT}
    Hˡ::HWedge{IT, RT, ST}
    eⁱᵝ::FixedSizeVectorDefault{Complex{RT}}
    cβ½::FixedSizeVectorDefault{RT}  # cos(β/2) per rotor; length 0 unless IT <: HalfOddInteger
    sβ½::FixedSizeVectorDefault{RT}  # sin(β/2) per rotor; length 0 unless IT <: HalfOddInteger
    d̄ₗ::FixedSizeVectorDefault{RT}  # d̄ₗ[m-ℓₘᵢₙ+1] = √δ²(ℓ, m); see `recurrence_coefficients!`
    ℓₘₐₓ::IT
    m′ₘₐₓ::IT
    swapH::Base.RefValue{Bool}  # h⃗ˡ(w) returns h⃗ᵃ if `false`, otherwise h⃗ᵇ; and vice versa for h⃗ˡ⁺¹(w)
    axes_valid::Base.RefValue{Bool}  # h⃗ˡ and h⃗ˡ⁺¹ hold correct data for their ℓ labels
    # The recurrence steps index every buffer under `@inbounds`, for each of the wedge's `Nᵣ`
    # rotors, so each buffer must hold exactly that many, and the table of coefficients must
    # hold one entry for each m ∈ ℓₘᵢₙ:ℓ-1 at the largest ℓ of the wedge.  `allocate_H` always
    # builds them so; this checks the buffers wherever they come from, once, as they are
    # brought together.
    function HCalculator{IT, RT, ST}(
        h⃗ᵃ, h⃗ᵇ, Hˡ, eⁱᵝ, cβ½, sβ½, d̄ₗ, ℓₘₐₓ, m′ₘₐₓ, swapH, axes_valid
    ) where {IT, RT<:Real, ST}
        let n = Nᵣ(Hˡ), nₕ = IT <: HalfOddInteger ? Nᵣ(Hˡ) : 0, nₗ = Int(maxℓ(Hˡ) - lowest_index(IT))
            if !(
                Nᵣ(h⃗ᵃ) == n && Nᵣ(h⃗ᵇ) == n && length(eⁱᵝ) == n
                && length(cβ½) == nₕ && length(sβ½) == nₕ && length(d̄ₗ) ≥ nₗ
            )
                throw(DimensionMismatch(
                    "The buffers of an HCalculator must each hold one entry per rotor of its "
                    * "wedge, Nᵣ=$n (the half angles $nₕ), and the table of coefficients "
                    * "$nₗ entries; the axes hold $(Nᵣ(h⃗ᵃ)) and $(Nᵣ(h⃗ᵇ)), e^{iβ} has length "
                    * "$(length(eⁱᵝ)), the half angles $(length(cβ½)) and $(length(sβ½)), and "
                    * "the table $(length(d̄ₗ))."
                ))
            end
        end
        new{IT, RT, ST}(h⃗ᵃ, h⃗ᵇ, Hˡ, eⁱᵝ, cβ½, sβ½, d̄ₗ, ℓₘₐₓ, m′ₘₐₓ, swapH, axes_valid)
    end
end

@index_methods function HCalculator(
    β, ℓₘₐₓ::IT; mp_max::IndexType=ℓₘₐₓ, m′ₘₐₓ::IndexType=mp_max
) where {IT<:IndexType}
    # `floattype` is computed from the type of `β` alone, so the compiler settles it, and
    # with it the concrete type of the calculator.
    set_rotors!(allocate_H(IT, floattype(β), ℓₘₐₓ, m′ₘₐₓ, nrotors(β)), β)
end

# Allocate the buffers without touching them.  PRIVATE, and deliberately so: the returned
# calculator's rotor data (`eⁱᵝ`, and the half angles on the half-integer path) is
# uninitialized, and nothing in the object records that fact.  Every caller must therefore
# either `set_rotors!` it or copy those buffers in from another calculator before the object
# can escape — which is exactly what the public constructor above and `similar` below do.
function allocate_H(
    ::Type{IT}, ::Type{RT}, ℓₘₐₓ::IT, m′ₘₐₓ::IT, Nᵣ::Int
) where {IT<:IntegerHalf, RT<:Real}
    validate_degree(ℓₘₐₓ)
    validate_axis(ℓₘₐₓ, m′ₘₐₓ, -m′ₘₐₓ, "m′")
    # The axis buffers hold the *integer* orders 0:axisℓₘₐₓ.  For integer ℓ that is
    # 0:ℓₘₐₓ+1, one extra because step 3 reads the ℓ+1 axis; for half-integer ℓ it is
    # 0:jₘₐₓ+1 with jₘₐₓ = ℓₘₐₓ - 1/2, one extra because the axis advance computes it anyway.
    axisℓₘₐₓ = Int(ℓₘₐₓ - lowest_index(IT)) + 1
    h⃗ᵃ = HAxis(RT, Nᵣ, axisℓₘₐₓ)
    h⃗ᵇ = HAxis(RT, Nᵣ, axisℓₘₐₓ)
    h⃗ᵇ.ℓ = 1
    Hˡ = HWedge(RT, Nᵣ, ℓₘₐₓ, m′ₘₐₓ)
    eⁱᵝ = FixedSizeVector{Complex{RT}}(undef, Nᵣ)
    # The half-angle pair (cos(β/2), sin(β/2)) is needed only by the half-integer seed; the
    # integer path gets zero-length buffers, so it allocates and touches nothing extra.
    Nₕ = IT <: HalfOddInteger ? Nᵣ : 0
    cβ½ = FixedSizeVector{RT}(undef, Nₕ)
    sβ½ = FixedSizeVector{RT}(undef, Nₕ)
    # The table of the m-side coefficients of steps 4 and 5, refilled by every `recurrence!`,
    # so that it is not rotor data and is not copied by `copy_rotor_data!`.
    d̄ₗ = FixedSizeVector{RT}(undef, Int(ℓₘₐₓ - lowest_index(IT)))
    HCalculator{IT, RT, typeof(parent(Hˡ))}(
        h⃗ᵃ, h⃗ᵇ, Hˡ, eⁱᵝ, cβ½, sβ½, d̄ₗ, ℓₘₐₓ, m′ₘₐₓ, Ref(false), Ref(false)
    )
end

# The rotor data is copied buffer-by-buffer rather than by re-running `set_rotors!`, because
# no calculator stores the rotor it was given: `eⁱᵝ` alone would fix `β` only modulo 2π, and
# so would flip the double-cover sign (-1)^{2ℓ} on the half-integer path.  The half-angle
# copies are no-ops on the integer path, where those buffers have length zero.  Every
# calculator built on an `HCalculator` copies its rotor data through this, so that anything
# added to the rotor data of an `HCalculator` is copied everywhere.
function copy_rotor_data!(w′::HCalculator, w::HCalculator)
    copyto!(w′.eⁱᵝ, w.eⁱᵝ)
    copyto!(w′.cβ½, w.cβ½)
    copyto!(w′.sβ½, w.sβ½)
    w′
end
function Base.similar(w::HCalculator{IT, RT}) where {IT, RT}
    copy_rotor_data!(allocate_H(IT, RT, w.ℓₘₐₓ, w.m′ₘₐₓ, Nᵣ(w)), w)
end
function Base.similar(w::HCalculator{IT, RT}, β) where {IT, RT}
    check_rotor_count(w, β)
    check_rotor_type(w, β)
    set_rotors!(allocate_H(IT, RT, w.ℓₘₐₓ, w.m′ₘₐₓ, Nᵣ(w)), β)
end

# The wedge's own `ℓ` field says which order its storage was last laid out for, which is not
# the same thing: it starts at ℓₘᵢₙ before anything is computed, and it keeps its value when
# `set_β!` or `fill!` leaves the stored numbers stale.  `axes_valid` is what records whether the
# current rotor data have been taken through the recurrence at all (it is also what `show`
# consults), so it decides, as the other calculators' `ℓ` fields do.
ℓ(w::HCalculator) = w.axes_valid[] ? Hˡ(w).ℓ : ℓₘᵢₙ(w) - 1
ℓₘₐₓ(w::HCalculator) = w.ℓₘₐₓ
m′ₘₐₓ(w::HCalculator) = w.m′ₘₐₓ
floattype(::Type{<:HCalculator{IT, RT}}) where {IT, RT} = RT
m′ₘᵢₙ(w::HCalculator) = -w.m′ₘₐₓ
Nᵣ(w::HCalculator) = Nᵣ(Hˡ(w))

h⃗ˡ(w::HCalculator) = w.swapH[] ? w.h⃗ᵇ : w.h⃗ᵃ
h⃗ˡ⁺¹(w::HCalculator) = w.swapH[] ? w.h⃗ᵃ : w.h⃗ᵇ
Hˡ(w::HCalculator) = w.Hˡ
eⁱᵝ(w::HCalculator) = w.eⁱᵝ

"""
    fill!(w::HCalculator, v)

Fill every internal buffer of `w` (both axes, the wedge, and the table of recurrence
coefficients) with the value `v`, and mark the axis data as invalid.  The stored rotor data
— the phases `e^{iβ}` and, on the half-integer path, the half angles — is deliberately left
alone, so that `recurrence!(w, ℓ)` still has everything it needs.  Useful for testing that no
uninitialized storage is ever read: everything the recurrence is responsible for writing is
poisoned, while everything the rotor data consists of is preserved.
"""
function Base.fill!(w::HCalculator{IT, RT}, v::Real) where {IT, RT}
    let v = convert(RT, v)
        fill!(parent(w.h⃗ᵃ), v)
        fill!(parent(w.h⃗ᵇ), v)
        fill!(parent(w.Hˡ), v)
        fill!(w.d̄ₗ, v)
    end
    w.axes_valid[] = false
    w
end

function increment_axes!(w::HCalculator)
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
    axis_ℓ(w::HCalculator, ℓ)

Label of the integer axis that seeds the wedge of order `ℓ`: ``ℓ`` itself for integer
indices, and ``ℓ - 1/2`` for half-integer ones.  The axis buffers are always labelled by
this `Int`, never by `ℓ`.
"""
@inline axis_ℓ(::HCalculator{IT}, ℓ) where {IT} = Int(ℓ - lowest_index(IT))

# Copy the m′ = 0 row of the wedge from the integer axis.  Integer indices only: for
# half-integer ℓ there is no m′ = 0 row, and `recurrence_seed!` writes the rows m′ = ±1/2
# instead.
function fillHˡ₀ₘ!(w::HCalculator{IT}) where {IT<:Integer}
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

function Base.show(io::IO, w::HCalculator{IT, RT, ST}) where {IT, RT, ST}
    print(
        io,
        "HCalculator{$IT, $RT} for ℓₘₐₓ=$(ℓₘₐₓ(w)), ",
        "m′ₘₐₓ=$(m′ₘₐₓ(w)), Nᵣ=$(Nᵣ(w))",
        w.axes_valid[] ? ", currently at ℓ=$(ℓ(w))" : " (nothing computed yet)"
    )
end
Base.show(io::IO, ::MIME"text/plain", w::HCalculator) = show(io, w)


### Rotor data
#
# The rotor data of an `HCalculator` are the phase e^{iβ} of each rotor and, on the
# half-integer path, its half angles, as `spinor_phases` and `half_angles` compute them (see
# `src/calculators/rotors.jl`).

# Store the half-angle pair for rotor `i`.  Each of these is a no-op — emitting no code at
# all, and not even evaluating its source — for integer index types, whose buffers have
# length zero; that keeps the integer path bit-identical and allocation-free.
@inline set_half_angles!(::HCalculator{IT}, ::Int, ::Real, ::Real) where {IT<:Integer} = nothing
@inline function set_half_angles!(
    w::HCalculator{IT, RT}, i::Int, c::Real, s::Real
) where {IT<:HalfOddInteger, RT}
    @inbounds w.cβ½[i] = convert(RT, c)
    @inbounds w.sβ½[i] = convert(RT, s)
    nothing
end

@inline set_half_angles_from_phase!(::HCalculator{IT}, ::Int, ::Complex) where {IT<:Integer} = nothing
@inline function set_half_angles_from_phase!(
    w::HCalculator{IT}, i::Int, z::Complex
) where {IT<:HalfOddInteger}
    set_half_angles!(w, i, half_angles(z)...)
end

@inline set_half_angles_from_angle!(::HCalculator{IT}, ::Int, ::Real) where {IT<:Integer} = nothing
@inline function set_half_angles_from_angle!(
    w::HCalculator{IT, RT}, i::Int, β::Real
) where {IT<:HalfOddInteger, RT}
    # Straight from β, not from cis(β): this is the only input form that honours the true
    # 4π periodicity of half-integer d.
    let β½ = convert(RT, β) / 2
        set_half_angles!(w, i, cos(β½), sin(β½))
    end
end

# Store the rotor `R` as rotor `i` of the calculator — the phase e^{iβ} and, on the
# half-integer path, the half angles — and return the phases z₊ and z₋ that a calculator of
# 𝔇 or of the harmonics also needs.  Every calculator built on an `HCalculator` stores a
# rotor through this, so that anything added to its rotor data is stored everywhere.
@inline function store_rotor!(w::HCalculator{IT, RT}, i::Int, R) where {IT, RT}
    eⁱᵝ, z₊, z₋, cβ½, sβ½ = spinor_phases(R, RT)
    @inbounds w.eⁱᵝ[i] = eⁱᵝ
    set_half_angles!(w, i, cβ½, sβ½)
    (z₊, z₋)
end

function set_rotors!(w::HCalculator{IT, RT}, eⁱᵝ::AbstractVector{<:Complex}) where {IT, RT<:Real}
    # The loops below write the calculator's 1-based buffers at the input's own indices, under
    # `@inbounds`, so an offset vector would write outside them.
    Base.require_one_based_indexing(eⁱᵝ)
    check_rotor_count(w, eⁱᵝ)
    # The phase must be a unit complex number; anything else is almost certainly a mistake
    # (e.g., passing β itself as a complex number), so we refuse it rather than silently
    # producing garbage.  The comparison is written so that a NaN phase is refused too.  The
    # tolerance is generous enough for `cis(β)` in any precision, and the value is used as
    # given (not renormalized), so that for integer indices a phase and the corresponding
    # angle give bitwise-identical results; for half-integer ones the half angles are
    # reconstructed from the phase by `half_angles`, and agree with those of the angle to
    # within an ulp or so.  Validate everything before mutating anything, so that a rejected
    # input leaves the calculator's state untouched.
    tolerance = 64 * max(eps(RT), eps(float(real(eltype(eⁱᵝ)))))
    for i ∈ eachindex(eⁱᵝ)
        z = convert(Complex{RT}, eⁱᵝ[i])
        if !(abs(abs2(z) - 1) ≤ tolerance)
            throw(DomainError(
                eⁱᵝ[i],
                "The phase e^{iβ} must have unit modulus; got $(eⁱᵝ[i]) at index $i, "
                * "with |e^{iβ}|² = $(abs2(z))."
            ))
        end
    end
    # The axes are marked invalid before the first rotor is replaced, so that nothing
    # computed from the old rotors can be combined with the new ones, even if a replacement
    # were to fail part-way.
    w.axes_valid[] = false
    @inbounds for i ∈ eachindex(eⁱᵝ)
        z = convert(Complex{RT}, eⁱᵝ[i])
        w.eⁱᵝ[i] = z
        set_half_angles_from_phase!(w, i, z)
    end
    w
end
function set_rotors!(w::HCalculator, R)
    throw(ArgumentError(
        "Cannot set rotor data of type $(typeof(R)) for a calculator with Nᵣ=$(Nᵣ(w)); "
        * rotor_input_forms * "."
    ))
end
function set_rotors!(w::HCalculator{IT, RT}, β::AbstractVector{<:Real}) where {IT, RT<:Real}
    Base.require_one_based_indexing(β)  # as for eⁱᵝ above
    check_rotor_count(w, β)
    # As for the phases above, everything is validated before anything is replaced.  An
    # infinite angle has no phase — `cis`, and the `cos` and `sin` of the half angle, throw for
    # one in some types, such as `Float64` and `Double64`, and return NaN in others, such as
    # `BigFloat` — so it is refused here, for every type alike, rather than part-way through
    # the loop below with some rotors already replaced.  A NaN angle is refused with it, as a
    # NaN phase is above.
    for i ∈ eachindex(β)
        if !isfinite(β[i])
            throw(DomainError(β[i], "The angle of rotor $i is $(β[i]), so it has no phase."))
        end
    end
    w.axes_valid[] = false  # as for the phases above
    @inbounds for i ∈ eachindex(β)
        w.eⁱᵝ[i] = cis(convert(RT, β[i]))
        set_half_angles_from_angle!(w, i, β[i])
    end
    w
end
function set_rotors!(w::HCalculator{IT, RT}, R::AbstractVector{<:RotorLike}) where {IT, RT<:Real}
    Base.require_one_based_indexing(R)  # as for eⁱᵝ above
    check_rotor_count(w, R)
    # There is nothing further to validate: `spinor_phases` accepts every rotor, giving NaNs
    # for one with non-finite components rather than throwing.
    w.axes_valid[] = false  # as for the phases above
    @inbounds for i ∈ eachindex(R)
        store_rotor!(w, i, R[i])
    end
    w
end
function set_rotors!(w::HCalculator{IT, RT}, R::Union{Real, Complex, RotorLike}) where {IT, RT<:Real}
    check_rotor_count(w, R)
    set_rotors!(w, @SVector [R])
end


### Driver

"""
    recurrence!(calc, R, ℓ)
    recurrence!(calc, ℓ)

Compute the quantities for index ``ℓ`` in the calculator `calc`, and return the block holding
them.

In the first form, the rotor data `R` is stored in the calculator first.  For a calculator
with `Nᵣ` rotors, `R` is an `AbstractVector` of length `Nᵣ`; for `Nᵣ = 1` a single element
is also accepted.  The elements may be `Rotor`s or `Quaternion`s, or — for quantities that
depend only on ``β``, such as ``d`` and ``H`` — the angle ``β`` or the phase ``e^{iβ}``.

In the second form, the rotor data from the previous call is reused.  Successive calls with
``ℓ, ℓ+1, ℓ+2, …`` are the cheap path: each costs ``O(N_r ℓ^2)``.  Requesting a smaller ``ℓ``
than the current one restarts the recurrence from ``ℓ_{min}``, so a loop that reads two
neighboring ``ℓ`` together pays that restart at every step.

This is the only way to read a calculator by hand — the calculators are not indexed — and it
is what iterating one calls for each ``ℓ``:

```julia
for ℓ ∈ 2:ℓₘₐₓ
    𝔇ˡ = recurrence!(calc, ℓ)
    # 𝔇ˡ[m′, m] for m′, m ∈ -ℓ:ℓ
end
```

What comes back depends on the calculator, and in every case it is indexed by the natural
ranges rather than from 1:

| calculator | block |
|---|---|
| [`DCalculator`](@ref), [`dCalculator`](@ref) | [`WignerMatrix`](@ref), `[m′, m]` |
| [`sYlmCalculator`](@ref) for one spin weight | [`DegreeBlock`](@ref), `[m]` |
| [`sYlmCalculator`](@ref) for a range of them | [`SpinMatrix`](@ref), `[s, m]` |
| [`HCalculator`](@ref) | [`HWedge`](@ref), `[iᵣ, m′, m]` for ``m ≥ \\|m′\\|`` |

For a calculator built from a vector of rotor data — of any length, even one — each of the
others gains a leading rotor index, so that the first is `[iᵣ, m′, m]`; an `HWedge` has that
index even for a single rotor, where it is `1`.  One spin
weight of a block that holds several is `ₛYₗ[s, :]`.

!!! warning
    For every calculator but [`HCalculator`](@ref) the block is a *view* into storage that the
    next call overwrites, so `copy` it if it must outlive the step (the copy keeps the natural
    indices), or `collect` it for an ordinary 1-based array.  An `HCalculator` is the
    exception, and the sharper case: its wedge is one mutable object handed back by identity,
    so every call returns the same object, laid out afresh for the new ``ℓ``.  `copy(Hˡ)`
    gives an independent wedge that survives the next step; the wedge itself belongs to the
    calculator, and its `ℓ` must not be reassigned.

`ℓ` must be an index of the calculator's own kind — an `Int` for a calculator built with
integer indices, and a half-odd-integer, as a [`HalfOddInteger`](@ref) or a `Rational{Int}`
with denominator 2, for one built with half-integer indices — and must lie between
`ℓₘᵢₙ(calc)` and `ℓₘₐₓ(calc)`; anything else is refused with an `ArgumentError`.  An index
that is not of the calculator's kind, including an integer of another type such as an `Int8`
and a floating-point number such as `2.0`, is refused with a message that says how to write
it, as it is by every function that takes an index.

The rotor data given this way must be of the calculator's own floating-point type, exactly as
for [`set_R!`](@ref), [`set_β!`](@ref) and [`set_θ!`](@ref), and anything else is refused
rather than converted.  A phase given as a `Complex` number must have unit modulus (to within
rounding), and is used as given.
"""
function recurrence! end

function recurrence!(w::HCalculator, R, ℓ)
    check_ℓ(w, ℓ)
    check_rotor_type(w, R)  # as `set_β!` does, rather than silently converting
    set_rotors!(w, R)
    recurrence!(w, ℓ)
end
# Every `recurrence!(calc, R, ℓ)` runs this *before* `set_rotors!`, so that a bad `ℓ` leaves
# the calculator's stored rotor data untouched ("validate everything before mutating
# anything"; see `set_rotors!` below).  It is also what refuses an `ℓ` that is not an index
# of the calculator's kind, with the message of `checked_index`, which names `owner`: this
# calculator, or, when this is the `HCalculator` of a calculator of 𝔇, of `d`, or of the
# harmonics, that calculator, which is the one the caller holds.  It returns `ℓ` converted
# to the calculator's index type.
function check_ℓ(w::HCalculator{IT}, ℓ, owner=w) where {IT}
    ℓ = checked_index(IT, ℓ, owner, "ℓ")
    if ℓ < ℓₘᵢₙ(w) || ℓ > ℓₘₐₓ(w)
        throw(ArgumentError(
            "Requested ℓ=$(ℓ) is out of bounds [$(ℓₘᵢₙ(w)), $(ℓₘₐₓ(w))] for this calculator."
        ))
    end
    ℓ
end

# Advance (or restart) the integer axis buffers so that h⃗ˡ holds order `j` and h⃗ˡ⁺¹ holds
# `j+1`.  Shared by the integer and half-integer drivers; `j` is `ℓ` in the former case and
# `ℓ - 1/2` in the latter.
#
# Stepping on is valid only from axes of consecutive orders, which is what step 2 leaves and
# `increment_axes!` preserves.  Axes labelled otherwise — which a transform used from several
# tasks at once can leave behind — are rebuilt from the start rather than trusted, so that
# every later call succeeds instead of failing in the label check of step 2 or 3.
function advance_axes!(w::HCalculator, j::Int)
    if !w.axes_valid[] || h⃗ˡ(w).ℓ > j || h⃗ˡ⁺¹(w).ℓ != h⃗ˡ(w).ℓ + 1
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

function recurrence!(w::HCalculator{Int, RT}, ℓ) where {RT}
    ℓ = check_ℓ(w, ℓ)  # rejects an `ℓ` that is not an index of this calculator

    # Steps 1 and 2 (the ℓ recurrence along the m′=0 axis).  The axis buffers h⃗ˡ and h⃗ˡ⁺¹
    # hold the axes for two successive ℓ values.  We restart from ℓₘᵢₙ if they are invalid or
    # ahead of the requested ℓ, and otherwise advance them one ℓ at a time.
    advance_axes!(w, axis_ℓ(w, ℓ))

    # Steps 3, 4, and 5 (the m′ recurrences at fixed ℓ), filling the wedge Hˡ.
    Hˡ(w).ℓ = ℓ
    fillHˡ₀ₘ!(w)  # Copy h⃗ˡ₀ₘ to Hˡ₀ₘ
    recurrence_step3!(w)  # Hˡ⁺¹₀ₘ -> Hˡ₁ₘ
    recurrence_coefficients!(w)  # d̄ₗᵐ for steps 4 and 5
    recurrence_step4!(w)  # Hˡₘ′ₘ₋₁, Hˡₘ′₋₁ₘ, Hˡₘ′ₘ₊₁ -> Hˡₘ′₊₁ₘ
    recurrence_step5!(w)  # Hˡₘ′ₘ₋₁, Hˡₘ′₊₁ₘ, Hˡₘ′ₘ₊₁ -> Hˡₘ′₋₁ₘ
    # Step 6 (the symmetries) is never applied to the wedge itself; elements outside the
    # wedge are read through `wedge_value`/`wedge_source` when the results are materialized.
    Hˡ(w)
end

# The half-integer driver.  The only differences from the integer one are that the axis is
# run at the integer order j = ℓ - 1/2, and that steps "0" (`fillHˡ₀ₘ!`) and 3 — which
# together produce the rows m′ = 0 and m′ = 1 from the ℓ and ℓ+1 axes — are replaced by
# `recurrence_seed!`, which produces the rows m′ = ±1/2 from the j axis alone.  Steps 4 and
# 5 are then literally the same code (see the section "Steps to compute H" of
# `docs/src/50-notes/01-H_recurrence.md`).
function recurrence!(w::HCalculator{IT, RT}, ℓ) where {IT<:HalfOddInteger, RT}
    ℓ = check_ℓ(w, ℓ)  # rejects an `ℓ` that is not an index of this calculator
    advance_axes!(w, axis_ℓ(w, ℓ))
    Hˡ(w).ℓ = ℓ
    recurrence_seed!(w)   # h⃗ʲ₀ₘ -> Hˡ₊₁⁄₂ₘ, Hˡ₋₁⁄₂ₘ
    recurrence_coefficients!(w)  # d̄ₗᵐ for steps 4 and 5
    recurrence_step4!(w)  # Hˡₘ′ₘ₋₁, Hˡₘ′₋₁ₘ, Hˡₘ′ₘ₊₁ -> Hˡₘ′₊₁ₘ
    recurrence_step5!(w)  # Hˡₘ′ₘ₋₁, Hˡₘ′₊₁ₘ, Hˡₘ′ₘ₊₁ -> Hˡₘ′₋₁ₘ
    Hˡ(w)
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
# undefined behavior if that layout ever changed.  None of these loops is a reduction, so no
# arithmetic is reassociated and the results are unchanged; the only effect is vectorization.
# Measured on `recurrence_step5!`, the hottest of them, this is 1.5× at Nᵣ=8 with
# bit-identical output.
#
# At Nᵣ=1 there is nothing to vectorize over rotors, and the setup of a vector loop costs
# more than the one element it computes, so the inner loops of steps 4 and 5, which touch
# every element of the wedge, write that case out as a single statement.  The statement is
# the loop's own expression at i = 1, so the values are the same.  Since the coefficients on
# the m side are read from a table (see `recurrence_coefficients!`), the loop over m that
# holds the statement is then one the compiler vectorizes in turn, and measured on
# `recurrence_step5!` for one rotor at ℓ = 64 and ℓ = 200, the two together are about five
# times as fast as a loop over one rotor that takes its own square roots; either alone gains
# little.  The other loops run once per ℓ or once per row, and are left alone.

# The axis buffers are integer-indexed for every index type (for half-integer ℓ they encode
# the order j = ℓ - 1/2), so steps 1 and 2 are shared verbatim and never see a `Rational`.
"""
    recurrence_step1!(w::HCalculator)

Step 1 of the ``H`` recursion: set ``H^{0}_{0,0} = 1`` for every rotor, in the lower axis
buffer of `w`, `h⃗ˡ(w)`, whose order must be 0.
"""
function recurrence_step1!(w::HCalculator{IT}) where {IT}
    let h⃗⁰ = h⃗ˡ(w)
        if h⃗⁰.ℓ ≠ 0
            error("recurrence_step1! can only be called for ℓ=0; current ℓ=$(h⃗⁰.ℓ).")
        end
        @inbounds for i ∈ 1:Nᵣ(h⃗⁰)
            h⃗⁰[i] = 1  # h⃗⁰[i, 0, 0] = 1
        end
    end
    w
end

"""
    recurrence_step2!(w::HCalculator)

Step 2 of the ``H`` recursion: compute the axis ``H^{n}_{0,m}``, for ``0 ≤ m ≤ n``, from
``H^{n-1}_{0,m}``, for every rotor, where ``n`` is the order of the upper axis buffer of
`w`, `h⃗ˡ⁺¹(w)`, and ``n-1`` that of the lower one, `h⃗ˡ(w)`.  This is the recurrence of
Xing et al. (2020) for the normalized associated Legendre functions, in the notation of
Gumerov and Duraiswami's step 2.
"""
function recurrence_step2!(w::HCalculator{IT, RT}) where {IT, RT}
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

# The rows H^J_{±1/2, m} of the wedge of half-integer order J, for m ∈ 1/2:J, are computed
# at a cost of O(Nᵣ J) from the axis h⃗ʲ of the integer order j = J - 1/2, by Varshalovich
# Eqs. 4.8.2(14) and (15),
#
#     H^J_{+1/2, m} = [ √(J+m) c h_{m-1/2} - √(J-m) s h_{m+1/2} ] / √(J + 1/2),
#     H^J_{-1/2, m} = [ √(J+m) s h_{m-1/2} + √(J-m) c h_{m+1/2} ] / √(J + 1/2),
#
# with c = cos(β/2), s = sin(β/2), and h_k = h⃗ʲ₀ₖ = d^j_{0,k}(β).  Both rows are mandatory:
# the corner H_{-1/2,1/2} cannot be reached from the +1/2 row without leaving the wedge.
# (See "Step 3 for half-integer ℓ" in `docs/src/50-notes/01-H_recurrence.md`.)
"""
    recurrence_seed!(w::HCalculator{HalfOddInteger})

The half-integer form of step 3 of the ``H`` recursion: compute the rows
``H^{J}_{±1/2,m}`` of the wedge of `w`, for ``1/2 ≤ m ≤ J``, from the axis ``H^{j}_{0,k}``
of the integer order ``j = J - 1/2`` in the lower axis buffer, `h⃗ˡ(w)`, for every rotor.
These are the two rows from which steps 4 and 5 begin; for integer indices, those rows are
the row ``m' = 0``, copied from the axis, and the row ``m' = 1`` of step 3.
"""
function recurrence_seed!(w::HCalculator{IT, RT}) where {IT<:HalfOddInteger, RT}
    let Hˡ = Hˡ(w), h⃗ʲ = h⃗ˡ(w), cβ½ = w.cβ½, sβ½ = w.sβ½
        @inbounds let √=sqrt∘RT, Nᵣ=Nᵣ(Hˡ), Hp=parent(Hˡ), hp=parent(h⃗ʲ)
            J = Hˡ.ℓ
            half = lowest_index(IT)  # the index 1/2, as a `HalfOddInteger`
            j = J - half             # an `Int`: the order of the integer axis
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

"""
    recurrence_step3!(w::HCalculator{Int})

Step 3 of the ``H`` recursion, for integer indices: compute the row ``H^{ℓ}_{1,m}`` of the
wedge of `w`, for ``1 ≤ m ≤ ℓ``, from the axis ``H^{ℓ+1}_{0,m}`` in the upper axis buffer,
`h⃗ˡ⁺¹(w)`, for every rotor.  For half-integer indices, [`recurrence_seed!`](@ref) takes its
place.
"""
function recurrence_step3!(w::HCalculator{Int, RT}) where {RT}
    let Hˡ = Hˡ(w), h⃗ˡ⁺¹ = h⃗ˡ⁺¹(w), eⁱᵝ = eⁱᵝ(w)
        @inbounds let √=sqrt∘RT, ℓ=Hˡ.ℓ, Nᵣ = Nᵣ(Hˡ), m′ₘₐₓ=m′ₘₐₓ(Hˡ)
            if h⃗ˡ⁺¹.ℓ ≠ ℓ + 1
                error("Inconsistent axes in recurrence_step3!: ℓ(Hˡ)=$(ℓ), ℓ(h⃗ˡ⁺¹)=$(h⃗ˡ⁺¹.ℓ).")
            end
            if ℓ > 0 && m′ₘₐₓ ≥ 1
                c = 1 / √(ℓ*(ℓ+1))

                # Precompute base offset for m′=1 row in Hˡ
                r¹ = row_index(Hˡ, 1) - 1  # step 3 is integer-only, so `m′ = 1` exists

                for m ∈ 1:ℓ
                    āₗᵐ = √((ℓ+m+1)*(ℓ-m+1))
                    b̄ₗ₊₁ᵐ⁻¹ = √((ℓ-m+1)*(ℓ-m+2))
                    b̄ₗ₊₁⁻ᵐ⁻¹ = √((ℓ+m+1)*(ℓ+m+2))

                    # Column offsets in Hˡ row 1 and h⃗ˡ⁺¹ row 0
                    c¹ᵐ = Nᵣ * (m - 1)  # Hˡ[i, 1, m] has m′=1, so column is m-abs(1)=m-1
                    i¹ᵐ = r¹ + c¹ᵐ
                    i⁰ᵐ⁺¹ = Nᵣ * (m + 1)
                    i⁰ᵐ⁻¹ = Nᵣ * (m - 1)
                    i⁰ᵐ = Nᵣ * m

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

# Fill the table `w.d̄ₗ` with the coefficients d̄ₗᵐ = √δ²(ℓ, m) for m ∈ ℓₘᵢₙ:ℓ-1, at the ℓ of
# the wedge.  These are the coefficients on the m side of steps 4 and 5, which depend on ℓ
# and m but not on m′, so that the ladders would otherwise take the same square roots once
# for every row.  The table is refilled by every `recurrence!`, from the same expression the
# steps would evaluate, so it is never out of date, and the values are the same ones.  It is
# left alone when the wedge has no rows beyond those of the seed, since then neither step
# runs.
function recurrence_coefficients!(w::HCalculator{IT, RT}) where {IT, RT}
    let Hˡ = Hˡ(w), d̄ₗ = w.d̄ₗ
        if m′ₘₐₓ(Hˡ) > lowest_index(IT)
            @inbounds let √=sqrt∘RT, ℓ = Hˡ.ℓ
                for m ∈ lowest_index(IT):(ℓ - 1)
                    d̄ₗ[(m - lowest_index(IT)) + 1] = √(δ²(ℓ, m))
                end
            end
        end
    end
    w
end

# The index arithmetic is the same for integer and half-integer ℓ; only the starting m′
# differs (1 or 1/2), and it is a compile-time constant for each index type.  The
# coefficients d̄ₗᵐ on the m side are read from the table that `recurrence_coefficients!`
# fills, at position m - ℓₘᵢₙ + 1.
"""
    recurrence_step4!(w::HCalculator)

Step 4 of the ``H`` recursion: compute the rows ``H^{ℓ}_{m'+1,m}`` of the wedge of `w`, for
``m'`` from ``1`` (or ``1/2`` for half-integer indices) up to ``m'_{\\mathrm{max}} - 1`` and
``m'+1 ≤ m ≤ ℓ``, from ``H^{ℓ}_{m'-1,m}``, ``H^{ℓ}_{m',m-1}``, and ``H^{ℓ}_{m',m+1}``, for
every rotor.  The first rows read are ``m' = 0`` and ``1`` for integer indices, and
``m' = ±1/2`` for half-integer ones.
"""
function recurrence_step4!(w::HCalculator{IT, RT}) where {IT, RT}
    let Hˡ = Hˡ(w), d̄ₗ = w.d̄ₗ
        @inbounds let √=sqrt∘RT, Nᵣ=Nᵣ(Hˡ), ri=row_index(Hˡ), k₀ = lowest_index(IT)
            ℓ = Hˡ.ℓ
            m′ₘₐₓw = m′ₘₐₓ(Hˡ)
            m′ₘᵢₙw = m′ₘᵢₙ(Hˡ)
            for m′ ∈ (1 - lowest_index(IT)):(m′ₘₐₓw - 1)
                # The m-side signs sgn(m) and sgn(m-1) are +1 throughout the range visited
                # here (m ≥ m′+1 ≥ 3/2 > 0), so they are left out.  The m′-side sign is
                # *not* always +1: at m′ = 1/2 the coefficient of Hˡ[m′-1, m] picks up
                # sgn(-1/2) = -1.  (See step 4 in `docs/src/50-notes/01-H_recurrence.md`.)
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
                    d̄ₗᵐ⁻¹ = d̄ₗ[m - k₀]        # √(δ²(ℓ, m - 1))
                    d̄ₗᵐ = d̄ₗ[(m - k₀) + 1]    # √(δ²(ℓ, m))

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

                    # Hˡ[i, m′+1, m] = (
                    #     d̄ₗᵐ′⁻¹ * Hˡ[i, m′-1, m]
                    #     - d̄ₗᵐ⁻¹ * Hˡ[i, m′, m-1]
                    #     + d̄ₗᵐ * Hˡ[i, m′, m+1]
                    # ) / d̄ₗᵐ′
                    if Nᵣ == 1  # the same expression, for the one rotor
                        Hˡ[iᵐ′⁺¹ᵐ + 1] = (
                            d̄ₗᵐ′⁻¹ * Hˡ[iᵐ′⁻¹ᵐ + 1]
                            - d̄ₗᵐ⁻¹ * Hˡ[iᵐ′ᵐ⁻¹ + 1]
                            + d̄ₗᵐ * Hˡ[iᵐ′ᵐ⁺¹ + 1]
                        ) * inv_d̄ₗᵐ′
                    else
                        @simd ivdep for i ∈ 1:Nᵣ
                            Hˡ[iᵐ′⁺¹ᵐ + i] = (
                                d̄ₗᵐ′⁻¹ * Hˡ[iᵐ′⁻¹ᵐ + i]
                                - d̄ₗᵐ⁻¹ * Hˡ[iᵐ′ᵐ⁻¹ + i]
                                + d̄ₗᵐ * Hˡ[iᵐ′ᵐ⁺¹ + i]
                            ) * inv_d̄ₗᵐ′
                        end
                    end
                end

                # Now, we do the m=ℓ case separately, since there is no m+1, so we would get
                # out-of-bounds accesses; we just copy the body of the loop above, but
                # remove anything that involves m+1.
                let m = ℓ
                    d̄ₗᵐ⁻¹ = d̄ₗ[m - k₀]        # √(δ²(ℓ, m - 1))

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

# As in step 4, the code is shared between index types; only the starting m′ differs
# (0 or -1/2).  The m-side signs sgn(m) and sgn(m-1) are +1 throughout the range visited
# here (m ≥ -m′+1 ≥ 1, so m - 1 ≥ 0), and the coefficients themselves are read from the
# table, as in step 4.
"""
    recurrence_step5!(w::HCalculator)

Step 5 of the ``H`` recursion: compute the rows ``H^{ℓ}_{m'-1,m}`` of the wedge of `w`, for
``m'`` from ``0`` (or ``-1/2`` for half-integer indices) down to ``m'_{\\mathrm{min}} + 1``
and ``1 - m' ≤ m ≤ ℓ``, from ``H^{ℓ}_{m'+1,m}``, ``H^{ℓ}_{m',m-1}``, and
``H^{ℓ}_{m',m+1}``, for every rotor.
"""
function recurrence_step5!(w::HCalculator{IT, RT}) where {IT, RT}
    let Hˡ = Hˡ(w), d̄ₗ = w.d̄ₗ
        @inbounds let √=sqrt∘RT, Nᵣ=Nᵣ(Hˡ), ri=row_index(Hˡ), k₀ = lowest_index(IT)
            ℓ = Hˡ.ℓ
            m′ₘᵢₙw = m′ₘᵢₙ(Hˡ)
            for m′ ∈ (-lowest_index(IT)):-1:(m′ₘᵢₙw + 1)
                d̄ₗᵐ′ = sgn(m′) * √(δ²(ℓ, m′))
                d̄ₗᵐ′⁻¹ = sgn(m′ - 1) * √(δ²(ℓ, m′ - 1))
                inv_d̄ₗᵐ′⁻¹ = inv(d̄ₗᵐ′⁻¹)

                # Precompute base offsets for m′-1, m′, m′+1 rows
                rᵐ′⁻¹ = ri[((m′ - 1) - m′ₘᵢₙw) + 1] - 1
                rᵐ′ = ri[(m′ - m′ₘᵢₙw) + 1] - 1
                rᵐ′⁺¹ = ri[((m′ + 1) - m′ₘᵢₙw) + 1] - 1

                for m ∈ (-m′+1):(ℓ-1)
                    d̄ₗᵐ = d̄ₗ[(m - k₀) + 1]    # sgn(m) √(δ²(ℓ, m))
                    d̄ₗᵐ⁻¹ = d̄ₗ[m - k₀]        # sgn(m - 1) √(δ²(ℓ, m - 1))

                    # Compute column offsets within each row
                    cᵐ′⁺¹ᵐ = Nᵣ * (m - abs(m′ + 1))
                    cᵐ′ᵐ⁻¹ = Nᵣ * ((m - 1) - abs(m′))
                    cᵐ′ᵐ⁺¹ = Nᵣ * ((m + 1) - abs(m′))
                    cᵐ′⁻¹ᵐ = Nᵣ * (m - abs(m′ - 1))

                    iᵐ′⁺¹ᵐ = rᵐ′⁺¹ + cᵐ′⁺¹ᵐ
                    iᵐ′ᵐ⁻¹ = rᵐ′ + cᵐ′ᵐ⁻¹
                    iᵐ′ᵐ⁺¹ = rᵐ′ + cᵐ′ᵐ⁺¹
                    iᵐ′⁻¹ᵐ = rᵐ′⁻¹ + cᵐ′⁻¹ᵐ

                    # Hˡ[i, m′-1, m] = (
                    #     d̄ₗᵐ′ * Hˡ[i, m′+1, m]
                    #     + d̄ₗᵐ⁻¹ * Hˡ[i, m′, m-1]
                    #     - d̄ₗᵐ * Hˡ[i, m′, m+1]
                    # ) / d̄ₗᵐ′⁻¹
                    if Nᵣ == 1  # the same expression, for the one rotor
                        Hˡ[iᵐ′⁻¹ᵐ + 1] = (
                            d̄ₗᵐ′ * Hˡ[iᵐ′⁺¹ᵐ + 1]
                            + d̄ₗᵐ⁻¹ * Hˡ[iᵐ′ᵐ⁻¹ + 1]
                            - d̄ₗᵐ * Hˡ[iᵐ′ᵐ⁺¹ + 1]
                        ) * inv_d̄ₗᵐ′⁻¹
                    else
                        @simd ivdep for i ∈ 1:Nᵣ
                            Hˡ[iᵐ′⁻¹ᵐ + i] = (
                                d̄ₗᵐ′ * Hˡ[iᵐ′⁺¹ᵐ + i]
                                + d̄ₗᵐ⁻¹ * Hˡ[iᵐ′ᵐ⁻¹ + i]
                                - d̄ₗᵐ * Hˡ[iᵐ′ᵐ⁺¹ + i]
                            ) * inv_d̄ₗᵐ′⁻¹
                        end
                    end
                end
                let m = ℓ
                    d̄ₗᵐ⁻¹ = d̄ₗ[m - k₀]        # sgn(m - 1) √(δ²(ℓ, m - 1))

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
