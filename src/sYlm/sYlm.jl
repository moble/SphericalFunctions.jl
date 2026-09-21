"""
    HarmonicCalculator{IT, RT, NT}

Calculator producing the spin-weighted spherical harmonics ``{}_sY_{ℓ,m}`` (when `NT` is
`Complex{RT}`) or the real ``{}_sλ_{ℓ,m}`` (when `NT` is `RT`), for `Nᵣ` points at a time, one
``ℓ`` at a time.  Use the constructors [`sYlmCalculator`](@ref) and [`sλlmCalculator`](@ref).

Internally this wraps a [`WignerHCalculator`](@ref), which runs the recurrence that both
flavours share, plus a buffer holding the block for the current ``ℓ``.  The phase tables `Z₊`
and `Z₋` are empty for the real flavour, which is the whole of the saving: the ``H``
recurrence is real, and it is only the ``e^{-i(mα - sγ)}`` factor that ever made the result
complex.

Two further parameters are lifted into the type so that the return type of `calc[ℓ]` is
inferrable — `S`, which records whether the calculator serves one spin weight or a range, and
a `Bool` read by [`isbatched`](@ref); see the comment on the struct.
"""
struct HarmonicCalculator{IT, RT<:Real, NT<:Union{RT, Complex{RT}}, ST, S, B}
    # As for [`WignerCalculator`](@ref), the last parameter is `Nᵣ > 1`, lifted into the type so
    # that the branch in `spin_row` and `spin_block` — and hence the return type of `calc[ℓ]` —
    # is settled at compile time.  `S` does the same job for the spin weights: it is the index
    # type when the calculator was built for one of them and a `UnitRange` of it when it was
    # built for several, which is what decides whether a block has a spin axis at all.
    H::WignerHCalculator{IT, RT, ST}
    Yˡ::Array{NT, 3}  # [iᵣ, s, m] block for the current ℓ, using the leading m entries
    Z₊::Matrix{Complex{RT}}  # Z₊[k+1, iᵣ] = z₊^k for k ∈ 0:2ℓₘₐₓ
    Z₋::Matrix{Complex{RT}}  # Z₋[k+1, iᵣ] = z₋^k for k ∈ 0:2ℓₘₐₓ
    s::S
    ℓ::Base.RefValue{IT}  # ℓ of the block currently in Yˡ; ℓₘᵢₙ-1 if none
    phases::Base.RefValue{Bool}  # false when the rotor data are angles θ (ϕ = γ = 0)
end

"""
    sYlmCalculator(R, ℓₘₐₓ, s)

Calculator for the spin-weighted spherical harmonics ``{}_sY_{ℓ,m}`` for all ``ℓ ≤ ℓₘₐₓ``,
with elements of `Complex` of the rotor's own element type.  The spin weight `s` is either a
single value or an ascending range of them, such as `-2:2`; the calculator serves exactly
those, and works out for itself how much of the underlying recurrence they need.

The rotor is given at construction, so that the calculator is usable the moment it exists,
and so that the element type is the rotor's own: a `Rotor{Float32}` gives a `Float32`
calculator, and `floattype(calc)` reports it.  There is no argument to override that — to
compute in another type, convert the rotor.  Give `R` as an `AbstractVector` of rotors to
get a calculator that handles all `Nᵣ = length(R)` of them at once.  Later rotors are
supplied with [`set_R!`](@ref).

Only one ``ℓ`` is held at a time, and iteration is how to step through them:

```julia
calc = sYlmCalculator(R, ℓₘₐₓ, 2)
for (ℓ, ₛYₗ) ∈ calc
    # ₛYₗ[m] is available for m ∈ -ℓ:ℓ
end

calc = sYlmCalculator(R, ℓₘₐₓ, -2:2)
for (ℓ, ₛYₗ) ∈ calc
    # ₛYₗ[s, m] is available for s ∈ -2:2 and m ∈ -ℓ:ℓ
end
```

With `Nᵣ > 1` each block gains a leading rotor index, so it is read as `ₛYₗ[iᵣ, m]` or
`ₛYₗ[iᵣ, s, m]`.  The same block is what `calc[ℓ]` returns after an explicit
[`recurrence!`](@ref), and `calc[ℓ, s]` picks out the row of one spin weight of a calculator
built for several.  Each block is a view into the calculator's storage, overwritten by the
next step; `copy` it if it must survive (keeping the natural indices), or `collect` it to get
an ordinary 1-based array.  Wherever ``ℓ < |s|`` the values are zero.

[`spins`](@ref) reports the range of spin weights served, and [`spin`](@ref) the single value
when there is only one.  Lengthening the range costs storage and arithmetic in proportion to
its length, but it is the largest ``|s|`` alone that decides the size of the underlying Wigner
``H`` wedge: `sYlmCalculator(R, ℓₘₐₓ, 2)` and `sYlmCalculator(R, ℓₘₐₓ, -2:2)` run precisely
the same recurrence, and differ only in how much of it is read out.

Instead of rotors, the angle ``θ`` (or a vector of `Nᵣ` angles) may be given, at construction
or through [`set_θ!`](@ref).  This evaluates the harmonics at ``(θ, ϕ=0)``, which is what the
ring-based transforms need — though the values are still stored as complex numbers, and a
calculator that will only ever be used that way is better built as an
[`sλlmCalculator`](@ref), which stores them as reals.

The harmonics are defined by
```math
{}_sY_{ℓ,m}(𝐑) = (-1)^s \\sqrt{\\frac{2ℓ+1}{4π}}\\, \\overline{𝔇^{(ℓ)}_{m,-s}(𝐑)},
```
so that they are functions on the rotation group; for spherical coordinates use
`R = from_spherical_coordinates(θ, ϕ)`.  See the "Conventions" section of the documentation.

# Half-integer indices

`ℓₘₐₓ` and the spin weights may all be half-integers, spelled as `Rational`s with denominator
2 or as [`HalfOddInteger`](@ref)s, in which case ``ℓ``, ``m`` and ``s`` are all half-integers.
A range is written the same way, as `-3//2:3//2`.  The prefactor ``(-1)^s`` is then ``\\pm
i``; the principal branch ``(-1)^s ≡ e^{iπs} = i^{2s}`` settled in the "Conventions" section
is used, so the harmonics are not real multiples of ``\\overline{𝔇}``, and
``{}_sY_{ℓ,m}(θ, 0)`` is therefore not real (see [`sλlmCalculator`](@ref), which divides that
constant phase out).  The blocks are the same containers as for integer indices — a
[`DegreeBlock`](@ref) or [`DegreeBlockBatch`](@ref) for one spin weight, and a
[`SpinMatrix`](@ref) or [`SpinMatrixBatch`](@ref) for several — indexed exactly as described
above.  The flat
interfaces [`sYlm`](@ref), [`sYlm!`](@ref) and [`sYlm_matrix`](@ref) accept the same
spellings, and lay the values out in the canonical mode-weight ordering of [`Yindex`](@ref),
which holds for half-odd indices exactly as it does for integers.

See also [`sYlm`](@ref) and [`sYlm_matrix`](@ref) for simpler interfaces, and
[`WignerDCalculator`](@ref).
"""
const sYlmCalculator{IT, RT, ST, S, B} =
    HarmonicCalculator{IT, RT, Complex{RT}, ST, S, B} where {IT, RT<:Real, ST, S, B}

"""
    sλlmCalculator(θ, ℓₘₐₓ, s)

Calculator for the real functions

```math
{}_sλ_{ℓ,m}(θ) = {}_sY_{ℓ,m}(θ, 0) \\big/ i^{2s},
```

for all ``ℓ ≤ ℓₘₐₓ``, with elements of the angle's own real type.  The first argument is the
angle ``θ``, or an `AbstractVector` of `Nᵣ` of them; later values are supplied with
[`set_θ!`](@ref).  Otherwise this behaves exactly like [`sYlmCalculator`](@ref) — the same
blocks, the same iteration, the same half-integer spellings — but stores real numbers rather
than complex ones, and allocates no phase tables at all.

The division by ``i^{2s}`` is what makes the definition uniform in the kind of the indices.
For an integer spin weight ``i^{2s} = (-1)^s`` is already included in ``{}_sY_{ℓ,m}`` itself,
so ``{}_sλ_{ℓ,m}`` is just the harmonic evaluated at ``ϕ = 0``, exactly as the literature
writes it.  For a half-odd spin weight ``i^{2s}`` is ``\\pm i``, so ``{}_sY_{ℓ,m}(θ, 0)`` is
imaginary rather than real, and dividing that constant phase out is what leaves a real
function behind.  The result is bit-for-bit what the ring-based transforms used to extract
from a complex calculator by hand.

A `Rotor` is **not** accepted, here or through [`set_R!`](@ref): a rotor specifies the angles
``α`` and ``γ``, whose phases a real calculator cannot represent.  Use an
[`sYlmCalculator`](@ref) for that.

See also [`sλlm`](@ref) and [`sλlm_matrix`](@ref) for simpler interfaces, and
[`WignerdCalculator`](@ref), which stands in the same relation to [`WignerDCalculator`](@ref).
"""
const sλlmCalculator{IT, RT, ST, S, B} =
    HarmonicCalculator{IT, RT, RT, ST, S, B} where {IT, RT<:Real, ST, S, B}

# The spellings the spin-weight argument accepts: one index, however spelled, or a range of
# them.  The argument is typed rather than left open so that a call whose arguments are in the
# wrong order — `sYlmCalculator(3, 1, Float64)`, say — is still the `MethodError` it should be
# rather than a puzzling complaint from the index normalizer.
const SpinSpelling = Union{IndexSpelling, AbstractRange}

function sYlmCalculator(R, ℓₘₐₓ::IndexSpelling, s::SpinSpelling)
    sYlmCalculator_helper(R, spin_indices(ℓₘₐₓ, s)...)
end
function sYlmCalculator_helper(R, ℓₘₐₓ::IT, s) where {IT<:HalfInteger}
    # See the note on `WignerDCalculator`: the element type reaches `allocate_Y` as a type,
    # not as a value, so that the concrete result type is settled at compile time.
    RT = rotor_basetype(R)
    set_rotors!(allocate_Y(IT, RT, Complex{RT}, ℓₘₐₓ, s, nrotors(R)), R)
end

function sλlmCalculator(θ, ℓₘₐₓ::IndexSpelling, s::SpinSpelling)
    sλlmCalculator_helper(θ, spin_indices(ℓₘₐₓ, s)...)
end
function sλlmCalculator_helper(θ, ℓₘₐₓ::IT, s) where {IT<:HalfInteger}
    RT = rotor_basetype(θ)
    set_rotors!(allocate_Y(IT, RT, RT, ℓₘₐₓ, s, nrotors(θ)), θ)
end

# The largest and smallest ``|s|`` among the spin weights a calculator serves, and how many of
# them there are.  These are written from the endpoints rather than as `maximum(abs, s)` and
# `minimum(abs, s)` because a reduction over a range of `HalfOddInteger`s reaches for an
# identity element that the type deliberately lacks.  Note that `min_abs_spin` is not simply
# the smaller of the two magnitudes: a range that straddles zero contains the smallest index
# of its kind, which is 0 for integers and 1/2 for half-odd-integers.
max_abs_spin(s::HalfInteger) = abs(s)
max_abs_spin(s::AbstractUnitRange) = max(abs(first(s)), abs(last(s)))
min_abs_spin(s::HalfInteger) = abs(s)
function min_abs_spin(s::AbstractUnitRange{IT}) where {IT<:HalfInteger}
    first(s) ≤ 0 ≤ last(s) ? ℓₘᵢₙ(IT) : min(abs(first(s)), abs(last(s)))
end
nspins(::HalfInteger) = 1
nspins(s::AbstractUnitRange) = length(s)

# A single index of the same kind as the spin-weight argument, for unifying `ℓₘᵢₙ` against.
spin_representative(s::HalfInteger) = s
spin_representative(s::AbstractUnitRange) = first(s)

# Allocate the buffers without touching them.  PRIVATE: see the note on `allocate_H`.  Here
# the uninitialized state includes the `phases` flag, which decides whether `Z₊` and `Z₋` are
# ever read, so this must not escape without a `set_rotors!` or a full buffer copy.
function allocate_Y(
    ::Type{IT}, ::Type{RT}, ::Type{NT}, ℓₘₐₓ::IT, s::S, Nᵣ::Int
) where {IT<:HalfInteger, RT<:Real, NT<:Union{RT, Complex{RT}}, S}
    sₕ = max_abs_spin(s)
    if sₕ > ℓₘₐₓ
        error("The spin weights $s need |s| ≤ ℓₘₐₓ=$ℓₘₐₓ; the largest of them is $sₕ.")
    end
    # The wedge is sized by the largest |s| alone, and cannot be narrowed further for a single
    # spin weight, tempting though that looks.  `materialize!` asks for H[m, -s] with m over
    # the whole of -ℓ:ℓ.  Where |m| > |s| the stored representative is the row -s or +s
    # according to the sign of m; where |m| ≤ |s| the symmetries offer only the row ±m, and m
    # runs over both signs.  So a single spin weight touches every row of -|s|:|s| just as the
    # full range would, and |s| is the smallest wedge that can serve it.
    H = allocate_H(IT, RT, ℓₘₐₓ, sₕ, Nᵣ)
    Yˡ = Array{NT, 3}(undef, Nᵣ, nspins(s), 2ℓₘₐₓ + 1)
    # The phase tables are what a real calculator does not have: with `K = 0` they are empty
    # rather than merely unread, so the ``{}_sλ_{ℓ,m}`` flavour allocates nothing for them.
    # This is the same trick `allocate_W` uses to separate ``𝔇`` from ``d``.
    K = NT <: Complex ? Int(2ℓₘₐₓ) + 1 : 0
    Z₊ = Matrix{Complex{RT}}(undef, K, Nᵣ)
    Z₋ = Matrix{Complex{RT}}(undef, K, Nᵣ)
    # The reference is typed explicitly because `ℓₘᵢₙ(IT) - 1` is an `Int` whenever `IT` is a
    # narrower integer type, and a bare `Ref` of it would not fit the `RefValue{IT}` field.
    HarmonicCalculator{IT, RT, NT, typeof(parent(H.Hˡ)), S, Nᵣ > 1}(
        H, Yˡ, Z₊, Z₋, s, Ref{IT}(ℓₘᵢₙ(IT) - 1), Ref(NT <: Complex)
    )
end

# See the notes on `similar(::WignerCalculator)` for why the assertion is here and why the
# rotor data is copied rather than re-derived.  The `phases` flag is part of that data: it
# records whether the calculator was given rotors or bare angles, and hence whether `Z₊` and
# `Z₋` hold anything at all.
function Base.similar(c::HarmonicCalculator{IT, RT, NT, ST, S, B}) where {IT, RT, NT, ST, S, B}
    c′ = allocate_Y(IT, RT, NT, ℓₘₐₓ(c), c.s, Nᵣ(c))::HarmonicCalculator{IT, RT, NT, ST, S, B}
    copyto!(c′.H.eⁱᵝ, c.H.eⁱᵝ)
    copyto!(c′.H.cβ½, c.H.cβ½)
    copyto!(c′.H.sβ½, c.H.sβ½)
    copyto!(c′.Z₊, c.Z₊)
    copyto!(c′.Z₋, c.Z₋)
    c′.phases[] = c.phases[]
    c′
end
function Base.similar(c::HarmonicCalculator{IT, RT, NT, ST, S, B}, R) where {IT, RT, NT, ST, S, B}
    if nrotors(R) != Nᵣ(c)
        error("This calculator handles Nᵣ=$(Nᵣ(c)) rotors, but got $(nrotors(R)).")
    end
    check_rotor_type(c, R)
    set_rotors!(
        allocate_Y(IT, RT, NT, ℓₘₐₓ(c), c.s, Nᵣ(c))::HarmonicCalculator{IT, RT, NT, ST, S, B}, R
    )
end

ℓ(c::HarmonicCalculator) = c.ℓ[]
ℓₘᵢₙ(c::HarmonicCalculator{IT}) where {IT} = ℓₘᵢₙ(IT)
ℓₘₐₓ(c::HarmonicCalculator) = ℓₘₐₓ(c.H)
floattype(::HarmonicCalculator{IT, RT}) where {IT, RT} = RT
# The element type of the values themselves — `Complex{RT}` for ₛYₗₘ, `RT` for ₛλₗₘ.  The
# flat interfaces size their output buffers by this, so that one set of helpers serves both.
number_type(::HarmonicCalculator{IT, RT, NT}) where {IT, RT, NT} = NT
# Which of the two constructors produced a calculator, for messages that should name the
# type the caller actually wrote rather than the struct they share.
flavour_name(::Type{<:Complex}) = "sYlmCalculator"
flavour_name(::Type{<:Real}) = "sλlmCalculator"
Nᵣ(c::HarmonicCalculator) = Nᵣ(c.H)
isbatched(::HarmonicCalculator{IT, RT, NT, ST, S, B}) where {IT, RT, NT, ST, S, B} = B

# `spins` answers in the same currency however the calculator was built, so that code which
# does not care can loop over it either way; `spin` exists only where there is a single value
# to name, and a calculator built for several gives a `MethodError` rather than a value that
# would have to be wrong.
spins(c::HarmonicCalculator{IT, RT, NT, ST, S}) where {IT, RT, NT, ST, S<:HalfInteger} = c.s:c.s
spins(c::HarmonicCalculator{IT, RT, NT, ST, S}) where {IT, RT, NT, ST, S<:AbstractUnitRange} = c.s
spin(c::HarmonicCalculator{IT, RT, NT, ST, S}) where {IT, RT, NT, ST, S<:HalfInteger} = c.s

function Base.show(io::IO, c::HarmonicCalculator{IT, RT, NT}) where {IT, RT, NT}
    print(
        io,
        "$(flavour_name(NT)){$IT, $RT} for ℓₘₐₓ=$(ℓₘₐₓ(c)), s=$(c.s), Nᵣ=$(Nᵣ(c))",
        c.ℓ[] < ℓₘᵢₙ(c) ? " (nothing computed yet)" : ", currently at ℓ=$(c.ℓ[])"
    )
end
function Base.show(io::IO, ::MIME"text/plain", c::HarmonicCalculator)
    show(io, c)
end

"""
    fill!(c::sYlmCalculator, v)

Fill the axis, wedge and output buffers of `c` with the value `v` and mark the current
results as invalid.  The stored rotor data — `e^{iβ}`, the half angles, the phase powers
`Z₊`, `Z₋`, and the flag recording whether rotors or bare angles were given — is deliberately
*not* touched, so `recurrence!(c, ℓ)` still has everything it needs, exactly as for
[`WignerHCalculator`](@ref).  Useful for testing that no uninitialized storage is ever read.
"""
function Base.fill!(c::HarmonicCalculator{IT, RT, NT}, v::Number) where {IT, RT, NT}
    fill!(c.H, real(v))
    fill!(c.Yˡ, convert(NT, v))
    c.ℓ[] = ℓₘᵢₙ(IT) - 1
    c
end


### Rotor data: full rotors (as for 𝔇), or angles θ meaning (θ, ϕ=0)

function set_rotors!(c::HarmonicCalculator{IT, RT}, θ::AbstractVector{<:Real}) where {IT, RT<:Real}
    set_rotors!(c.H, θ)
    c.phases[] = false
    c.ℓ[] = ℓₘᵢₙ(IT) - 1
    c
end
function set_rotors!(c::HarmonicCalculator{IT, RT}, θ::Real) where {IT, RT<:Real}
    set_rotors!(c.H, θ)
    c.phases[] = false
    c.ℓ[] = ℓₘᵢₙ(IT) - 1
    c
end
function set_rotors!(c::sYlmCalculator{IT, RT}, R::AbstractVector{<:Rotor}) where {IT, RT<:Real}
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
    c.phases[] = true
    c.ℓ[] = ℓₘᵢₙ(IT) - 1
    c
end
function set_rotors!(c::sYlmCalculator{IT, RT}, R::Rotor) where {IT, RT<:Real}
    if Nᵣ(c) != 1
        error("A single rotor was given, but this calculator expects Nᵣ=$(Nᵣ(c)) rotors.")
    end
    set_rotors!(c, @SVector [R])
end
# A real calculator refuses rotors rather than silently dropping the phases they specify,
# exactly as `set_β!` refuses them for the real `d` matrices.  The check is a method rather
# than a branch so that the refusal happens at the outermost call, naming the type the caller
# actually has.
function set_rotors!(c::sλlmCalculator, R::Union{Rotor, AbstractVector{<:Rotor}})
    error(
        "An sλlmCalculator evaluates at (θ, ϕ=0) and stores real numbers, so it cannot take "
        * "a rotor, whose α and γ angles are phases it has nowhere to put.  Give the angle θ "
        * "instead, or use an sYlmCalculator."
    )
end
# One catch-all rather than two, so that it stays strictly less specific than every method
# above: a pair of them, keyed on the flavour, would be ambiguous with the angle methods,
# which name the argument type but not the calculator's.  The branch is on a type parameter,
# so only the message costs anything, and only on the way to an error.
function set_rotors!(c::HarmonicCalculator{IT, RT, NT}, R) where {IT, RT<:Real, NT}
    if NT <: Complex
        error(
            "An sYlmCalculator needs rotors (as `Rotor`s) or angles θ::Real — one, or an "
            * "AbstractVector of $(Nᵣ(c)) of either — not $(typeof(R))."
        )
    else
        error(
            "An sλlmCalculator needs angles θ::Real — one, or an AbstractVector of "
            * "$(Nᵣ(c)) of them — not $(typeof(R))."
        )
    end
end


### Driver

function recurrence!(c::HarmonicCalculator, R, ℓ)
    check_ℓ(c.H, ℓ)
    set_rotors!(c, R)
    recurrence!(c, ℓ)
end
function recurrence!(c::HarmonicCalculator{IT}, ℓ) where {IT}
    let ℓ = convert(IT, ℓ)
        recurrence!(c.H, ℓ)
        materialize!(c, ℓ)
    end
    c
end

# (-1)^s, as e^{iπs} = i^{2s}; for integer s this is ±1.
minus_one_to_the(s::Integer) = ifelse(iseven(s), 1, -1)

# i^k for an integer k, exactly.  With k = 2s this is the principal branch of (-1)^s chosen
# by the conventions (see `docs/src/30-conventions/02-details.md`, "Half-integer indices"): ±1 for
# integer s, and ±i for half-integer s.
@inline function im_power(::Type{RT}, k::Integer) where {RT<:Real}
    let o = one(RT), z = zero(RT)
        (Complex{RT}(o, z), Complex{RT}(z, o), Complex{RT}(-o, z), Complex{RT}(z, -o))[mod(k, 4) + 1]
    end
end

# Index of spin weight s along the second dimension of c.Yˡ.  For a calculator built for a
# single spin weight this is 1, whatever that weight is.
@inline spin_index(c::HarmonicCalculator, s) = Int(s - first(spins(c))) + 1

# The coefficient (-1)^s ϵ(m) ϵ(s) σ √((2ℓ+1)/4π) of Hₘ,₋ₛ in ₛYₗₘ, in whichever number type
# `NT` the calculator stores.  For integer indices it is real whatever `NT` is, and stays real
# so that the integer path is bit-for-bit what it always was.  For half-odd indices
# (-1)^s = i^{2s} is ±i, which a complex calculator includes and a real one divides out — that
# division being the whole of the difference between ₛYₗₘ(θ,0) and ₛλₗₘ(θ).
#
# The two half-odd coefficients differ by exactly the factor i^{2s}, and `Complex * Real` is
# computed componentwise, so the real flavour's value is bit-for-bit ±imag of the complex
# flavour's — which is what the transforms used to extract by hand, and what the regression
# test asserts with `==` rather than `≈`.
@inline function sYlm_coefficient(
    ::Type{NT}, ::Type{RT}, σ::Int, m::IT, s::IT, prefactor::RT
) where {NT, IT<:Integer, RT}
    convert(RT, σ * ϵ(m) * ϵ(s) * minus_one_to_the(s)) * prefactor
end
@inline function sYlm_coefficient(
    ::Type{NT}, ::Type{RT}, σ::Int, m::IT, s::IT, prefactor::RT
) where {NT<:Complex, IT<:HalfOddInteger, RT}
    (convert(RT, σ * ϵ(m) * ϵ(s)) * prefactor) * im_power(RT, 2s)
end
@inline function sYlm_coefficient(
    ::Type{NT}, ::Type{RT}, σ::Int, m::IT, s::IT, prefactor::RT
) where {NT<:Real, IT<:HalfOddInteger, RT}
    convert(RT, σ * ϵ(m) * ϵ(s)) * prefactor
end

# Write ₛYₗₘ for every spin weight the calculator serves and every m ∈ -ℓ:ℓ into
# c.Yˡ[iᵣ, s, m].  From the definition ₛYₗₘ = (-1)^s √((2ℓ+1)/4π) conj(𝔇ₘ,₋ₛ), with
# 𝔇ₘ,₋ₛ = ϵ(m) ϵ(s) Hₘ,₋ₛ e^{-i(mα - sγ)}, and e^{i(mα - sγ)} = z₊^(m-s) z₋^(m+s).  As in the
# Wigner `materialize!`, nothing here depends on whether the indices are integers or
# half-odd-integers.
function materialize!(c::HarmonicCalculator{IT, RT, NT}, ℓ::IT) where {IT, RT, NT}
    let H = c.H.Hˡ, Yˡ = c.Yˡ, Z₊ = c.Z₊, Z₋ = c.Z₋, Nᵣ = Nᵣ(c), Hp = parent(H)
        if H.ℓ != ℓ
            error("The H wedge holds ℓ=$(H.ℓ), but ℓ=$ℓ was requested.")
        end
        prefactor = √((2ℓ + 1) / (4 * RT(π)))
        srange = spins(c)
        m′ₘₐₓw = m′ₘₐₓ(H)
        m′ₘᵢₙw = m′ₘᵢₙ(H)
        @inbounds for (j, m) ∈ enumerate(-ℓ:ℓ)
            for (i, s) ∈ enumerate(srange)
                if abs(s) > ℓ
                    for iᵣ ∈ 1:Nᵣ
                        Yˡ[iᵣ, i, j] = 0
                    end
                    continue
                end
                a, b, σ = wedge_source(m, -s, m′ₘₐₓw)
                offset = wedge_offset(H, a, b, m′ₘᵢₙw)
                coefficient = sYlm_coefficient(NT, RT, σ, m, s, prefactor)
                # `NT <: Complex` is a compile-time constant, so a real calculator never
                # compiles the phase branch at all — which is what lets `Z₊` and `Z₋` be
                # empty rather than merely unread.
                if NT <: Complex && c.phases[]
                    k₊ = m - s
                    k₋ = m + s
                    for iᵣ ∈ 1:Nᵣ
                        phase = zpower(Z₊, iᵣ, k₊) * zpower(Z₋, iᵣ, k₋)
                        Yˡ[iᵣ, i, j] = coefficient * Hp[offset + iᵣ] * phase
                    end
                else
                    @simd for iᵣ ∈ 1:Nᵣ
                        Yˡ[iᵣ, i, j] = coefficient * Hp[offset + iᵣ]
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
    calc[ℓ, s]

The spin-weighted spherical harmonics ``{}_sY_{ℓ,m}`` for the current ``ℓ`` of the calculator,
as an array indexed naturally.

The first form gives everything the calculator serves.  For a calculator built for a single
spin weight that is `calc[ℓ][m]`, or `calc[ℓ][iᵣ, m]` when `Nᵣ > 1`; for one built for a range
of spin weights it is `calc[ℓ][s, m]`, or `calc[ℓ][iᵣ, s, m]`.  The second form picks the row
of one spin weight out of the range, and so always has the shape of the single-spin-weight
case.  In both, `m ∈ -ℓ:ℓ`.  One spin weight can equally be sliced out of a block that holds
several, as `calc[ℓ][s, :]`, which is spelled the same way whichever kind of index the
calculator has.

The result is a view into the calculator's storage, valid until the next call to
[`recurrence!`](@ref).  `ℓ` must be the value passed to the most recent `recurrence!`, and `s`
must be one of [`spins`](@ref)`(calc)`.  Wherever ``ℓ < |s|`` the elements are zero.

The result is a [`DegreeBlock`](@ref), [`DegreeBlockBatch`](@ref), [`SpinMatrix`](@ref) or
[`SpinMatrixBatch`](@ref) according to its shape, for either kind of index.  `copy` keeps it
with its natural indices, `collect` gives an ordinary 1-based `Array`, and [`strided`](@ref)
gives a 1-based view of the same storage for linear algebra.
"""
function Base.getindex(c::HarmonicCalculator{IT, RT, NT, ST, S}, ℓ) where {IT, RT, NT, ST, S<:HalfInteger}
    ℓ = convert(IT, ℓ)
    check_current_ℓ(c, ℓ)
    spin_row(c, ℓ, 1)
end
function Base.getindex(
    c::HarmonicCalculator{IT, RT, NT, ST, S}, ℓ
) where {IT, RT, NT, ST, S<:AbstractUnitRange}
    ℓ = convert(IT, ℓ)
    check_current_ℓ(c, ℓ)
    spin_block(c, ℓ)
end
function Base.getindex(c::HarmonicCalculator{IT}, ℓ, s) where {IT}
    ℓ = convert(IT, ℓ)
    s = convert(IT, s)
    check_current_ℓ(c, ℓ)
    check_spin(c, s)
    spin_row(c, ℓ, spin_index(c, s))
end

# The three ways a request for a block can be premature or out of range, each said separately
# because the three have quite different remedies.
function check_current_ℓ(c::HarmonicCalculator, ℓ)
    if c.ℓ[] < ℓₘᵢₙ(c)
        error(
            "This calculator currently holds no result, because nothing has been computed "
            * "yet; iterate it, or call `recurrence!(calc, ℓ)` first."
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
    nothing
end

function check_spin(c::HarmonicCalculator, s)
    let sr = spins(c)
        if !(first(sr) ≤ s ≤ last(sr))
            error("This calculator serves the spin weights $(c.s), so s=$s is not among them.")
        end
    end
    nothing
end

# One spin weight's row, selected by its position `i` in the calculator's storage, and the
# whole block of every spin weight.  `isbatched(c)` reads a type parameter, so every branch is
# resolved at compile time.  The same containers are returned for both kinds of index.
function spin_row(c::HarmonicCalculator{IT}, ℓ::IT, i::Int) where {IT<:HalfInteger}
    let mr = -ℓ:ℓ
        if isbatched(c)
            DegreeBlockBatch(view(c.Yˡ, :, i, 1:length(mr)), ℓ; mₘₐₓ=last(mr), mₘᵢₙ=first(mr))
        else
            DegreeBlock(view(c.Yˡ, 1, i, 1:length(mr)), ℓ; mₘₐₓ=last(mr), mₘᵢₙ=first(mr))
        end
    end
end

function spin_block(c::HarmonicCalculator{IT}, ℓ::IT) where {IT<:HalfInteger}
    let mr = -ℓ:ℓ, sr = spins(c), n = length(spins(c))
        if isbatched(c)
            SpinMatrixBatch(
                view(c.Yˡ, :, 1:n, 1:length(mr)), ℓ;
                sₘₐₓ=last(sr), sₘᵢₙ=first(sr), mₘₐₓ=last(mr), mₘᵢₙ=first(mr)
            )
        else
            SpinMatrix(
                view(c.Yˡ, 1, 1:n, 1:length(mr)), ℓ;
                sₘₐₓ=last(sr), sₘᵢₙ=first(sr), mₘₐₓ=last(mr), mₘᵢₙ=first(mr)
            )
        end
    end
end


### Convenience functions
#
# Each public function here is a boundary method: it accepts every spelling of an index —
# `Integer`, `HalfOddInteger`, or a `Rational` with denominator 2 — and of a range of them,
# normalizes them, and re-dispatches to a worker whose `where {IT<:HalfInteger}` signature is
# what the rest of the package sees.  The split is made at the function boundary rather than by
# a second method because `ℓₘᵢₙ` is a keyword argument, and keyword arguments take no part in
# dispatch: a method with `ℓₘᵢₙ::IT` in its signature would refuse `ℓₘᵢₙ=1//2` with a bare
# `TypeError` before anything could normalize it.  The normalization is `spin_indices`, from
# `half_odd_integer.jl`, which also refuses a mixture of the two kinds of index, and a range
# whose step is not 1, with an explanation.  For integer indices the normalization is the
# identity, inlined, so the integer path is unchanged.
#
# `ℓₘᵢₙ` defaults to the smallest ``|s|`` among the spin weights, which is `abs(s)` when there
# is only one and 0 (or 1/2) for a range that straddles zero.  That default is resolved here,
# after the two indices the caller wrote have been unified, rather than in the keyword's own
# default expression; `nothing` is what asks for it.  Doing it in that order is what keeps a
# message about a mixture of kinds citing the caller's own values rather than a default derived
# from one of them.
@inline function flat_indices(ℓₘₐₓ, s, ℓₘᵢₙ)
    ℓₘₐₓ, s = spin_indices(ℓₘₐₓ, s)
    if ℓₘᵢₙ === nothing
        (ℓₘₐₓ, s, min_abs_spin(s))
    else
        (ℓₘₐₓ, s, last(unify_indices(spin_representative(s), ℓₘᵢₙ)))
    end
end

"""
    sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s))
    sYlm(R⃗, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s))

The spin-weighted spherical harmonics ``{}_sY_{ℓ,m}(R)`` for all ``ℓₘᵢₙ ≤ ℓ ≤ ℓₘₐₓ`` and
``-ℓ ≤ m ≤ ℓ``, as a [`HarmonicValues`](@ref): indexed first by ``ℓ`` and then naturally, so
that `sYlm(R, ℓₘₐₓ, s)[ℓ][m]` is one value.  With `ℓₘᵢₙ` below `abs(s)` the entries for
``ℓ < |s|`` are zero.  The computation is done in the rotor's own floating-point type, and the
result is `Complex` of it; to compute in another type, convert the rotor.

The first argument may be a single rotor or a vector of them, and `s` may be a single spin
weight or an ascending range such as `-2:2`, which between them give a block four shapes:

| built for | `sY[ℓ]` is indexed |
|---|---|
| one rotor, one spin weight | `[m]` |
| many rotors, one spin weight | `[iᵣ, m]` |
| one rotor, a range of spin weights | `[s, m]` |
| many rotors, a range of spin weights | `[iᵣ, s, m]` |

For a range, `ℓₘᵢₙ` defaults to the smallest ``|s|`` in it, which is 0 (or 1/2) whenever the
range straddles zero, so rows below their own ``|s|`` are zero.

Underneath, the values are held in one array whose *last* axis is the modes in the canonical
ordering `[ₛYₗₘ for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ]` (see [`Yindex`](@ref)), and whose leading axes
are the rotors and spin weights.  [`strided`](@ref) hands that array back, which is the form a
product with mode weights takes; [`sYlm_matrix`](@ref) is the direct spelling of it for
callers who want the bare array.

For spherical coordinates use `R = from_spherical_coordinates(θ, ϕ)`.  For repeated
evaluation use an [`sYlmCalculator`](@ref) (with [`sYlm!`](@ref), or directly), which
allocates once.
The convention is ``{}_sY_{ℓ,m} = (-1)^s \\sqrt{(2ℓ+1)/4π}\\, \\overline{𝔇^{(ℓ)}_{m,-s}}``;
see the "Conventions" section of the documentation.

# Half-integer indices

`ℓₘₐₓ` and `s` may be half-integers, spelled as `Rational`s with denominator 2 — as in
`sYlm(R, 7//2, 1//2)`, or `sYlm(R, 7//2, -3//2:3//2)` — in which case every ``ℓ`` and ``m`` is
a half-odd-integer, and so is `ℓₘᵢₙ`, which may be as small as `1//2`.  The indices in one
call must all be of one kind, integers or half-odd-integers; a call that mixes them, such as
`sYlm(R, 7//2, 1)`, is an error.  For half-integer `s` the prefactor ``(-1)^s`` is
``i^{2s} = \\pm i``, so the values include that phase and are not real multiples of
``\\overline{𝔇}``; the choice of branch is explained under [`sYlmCalculator`](@ref).
"""
function sYlm(R::Rotor, ℓₘₐₓ::IndexSpelling, s::SpinSpelling; ℓₘᵢₙ=nothing)
    sYlm_helper(sYlmCalculator_helper, R, flat_indices(ℓₘₐₓ, s, ℓₘᵢₙ)...)
end
function sYlm_helper(make, R, ℓₘₐₓ::IT, s, ℓₘᵢₙ::IT) where {IT<:HalfInteger}
    check_sYlm_args(ℓₘₐₓ, s, ℓₘᵢₙ)
    # The calculator decides the element type, and the output buffer follows it, so that
    # there is exactly one place where that decision is made.  `make` is what chooses the
    # flavour: `sYlmCalculator_helper` for the complex harmonics, `sλlmCalculator_helper` for
    # the real ones.  Everything below is common to both.
    calc = make(R, ℓₘₐₓ, s)
    Y = allocate_sYlm(number_type(calc), s, ℓₘᵢₙ, ℓₘₐₓ)
    sYlm_helper!(Y, calc, R, s, ℓₘᵢₙ)
    HarmonicValues(Y, s, ℓₘᵢₙ, ℓₘₐₓ, 1)
end

# Many rotors at once.  The storage and the recursion are `sYlm_matrix`'s — that is the
# efficient path, and there is no reason to have two — so this labels the same array rather
# than computing it again.  `sYlm_matrix` remains the way to ask for the bare array.
function sYlm(R⃗::AbstractVector{<:Rotor}, ℓₘₐₓ::IndexSpelling, s::SpinSpelling; ℓₘᵢₙ=nothing)
    sYlm_batch_helper(sYlmCalculator_helper, R⃗, flat_indices(ℓₘₐₓ, s, ℓₘᵢₙ)...)
end
function sYlm_batch_helper(
    make, R⃗::AbstractVector, ℓₘₐₓ::IT, s, ℓₘᵢₙ::IT
) where {IT<:HalfInteger}
    check_sYlm_args(ℓₘₐₓ, s, ℓₘᵢₙ)
    HarmonicValues(sYlm_matrix_helper(make, R⃗, ℓₘₐₓ, s, ℓₘᵢₙ), s, ℓₘᵢₙ, ℓₘₐₓ, length(R⃗))
end

# The output of a flat call: a vector of modes for one spin weight, and a matrix of spin
# weights by modes for several.
function allocate_sYlm(::Type{T}, ::HalfInteger, ℓₘᵢₙ, ℓₘₐₓ) where {T}
    Vector{T}(undef, Ysize(ℓₘᵢₙ, ℓₘₐₓ))
end
function allocate_sYlm(::Type{T}, s::AbstractUnitRange, ℓₘᵢₙ, ℓₘₐₓ) where {T}
    Matrix{T}(undef, length(s), Ysize(ℓₘᵢₙ, ℓₘₐₓ))
end

# These checks hold for either kind of index, and for one spin weight or many.  The floor of
# ℓₘᵢₙ is 0 for integers and 1/2 for half-odd-integers, and `ℓₘᵢₙ < 0` is the right test for
# both, since no half-odd-integer lies between 0 and 1/2; the message names the floor of the
# kind at hand, which is what the `ℓₘᵢₙ` accessor gives for the index type.  (The accessor is
# called qualified, because the argument of the same name shadows it here.)
function check_sYlm_args(ℓₘₐₓ, s, ℓₘᵢₙ)
    sₕ = max_abs_spin(s)
    if sₕ > ℓₘₐₓ
        error("|s|=$sₕ exceeds ℓₘₐₓ=$ℓₘₐₓ; there are no such harmonics.")
    end
    if ℓₘᵢₙ < 0 || ℓₘᵢₙ > max(sₕ, ℓₘₐₓ)
        lowest = SphericalFunctions.ℓₘᵢₙ(typeof(ℓₘₐₓ))
        error("ℓₘᵢₙ=$ℓₘᵢₙ must satisfy $lowest ≤ ℓₘᵢₙ ≤ max(|s|, ℓₘₐₓ).")
    end
end

# The calculator forms of `sYlm!` take the spin weight and `ℓₘᵢₙ` apart from the calculator
# whose index type they must match, so a mismatch — an integer `s` handed to a half-integer
# calculator, say — is possible there in a way it is not for the other flat functions, where
# `where {IT}` unifies the indices.  `convert(IT, s)` would refuse it, but with a bare
# `InexactError` about the type; this says what the calculator's indices are instead.
function check_index_kind(::Type{IT}, x, name) where {IT<:HalfInteger}
    if !isindex(IT, x)
        kind, example = IT <: Integer ? ("integers", "3") : ("half-odd-integers", "7//2")
        error(
            "This calculator's indices are $kind, like $example, so $name must be one too; "
            * "got $x."
        )
    end
end

"""
    Ylm(R, ℓₘₐₓ; ℓₘᵢₙ=0)
    Ylm(R⃗, ℓₘₐₓ; ℓₘᵢₙ=0)

The ordinary scalar spherical harmonics ``Y_{ℓ,m}(R)`` for all ``ℓₘᵢₙ ≤ ℓ ≤ ℓₘₐₓ``, as a
[`HarmonicValues`](@ref) indexed by ``ℓ`` and then by ``m``: `Ylm(R, ℓₘₐₓ)[ℓ][m]`.

These are the spin-weight-zero case of the spin-weighted harmonics; this function is exactly
`sYlm(R, ℓₘₐₓ, 0; ℓₘᵢₙ)`, and everything [`sYlm`](@ref) says applies here too — including
that a vector of rotors gives blocks indexed `[iᵣ, m]`.  For spherical coordinates use
`R = from_spherical_coordinates(θ, ϕ)`.  See [`YlmCalculator`](@ref) for repeated
evaluation.

The ``ℓ`` arguments must be `Int`.  Half-integer ``ℓ`` goes with half-integer spin weight, so
there is no half-integer analogue of this function; see [`sYlm`](@ref) for half-integer spin
weights.
"""
function Ylm(R::Rotor, ℓₘₐₓ::Int; ℓₘᵢₙ::Int=0)
    sYlm(R, ℓₘₐₓ, 0; ℓₘᵢₙ)
end
function Ylm(R⃗::AbstractVector{<:Rotor}, ℓₘₐₓ::Int; ℓₘᵢₙ::Int=0)
    sYlm(R⃗, ℓₘₐₓ, 0; ℓₘᵢₙ)
end

"""
    YlmCalculator(R, ℓₘₐₓ)

Calculator for the ordinary scalar spherical harmonics ``Y_{ℓ,m}``, for ``ℓ ≤ ℓₘₐₓ``.

This is exactly `sYlmCalculator(R, ℓₘₐₓ, 0)`, and returns an [`sYlmCalculator`](@ref) rather
than a type of its own — there is nothing about spin weight zero to specialize.  Everything
`sYlmCalculator` documents applies, so blocks are indexed `Yₗ[m]` (or `Yₗ[iᵣ, m]` for a
collection of rotors), and it iterates as `for (ℓ, Yₗ) ∈ calc`.

``ℓₘₐₓ`` must be an `Int`: a half-integer ``ℓ`` goes with a half-integer spin weight, so
spin weight zero has no half-integer analogue.  See [`sYlmCalculator`](@ref) for those.
"""
YlmCalculator(R, ℓₘₐₓ::Int) = sYlmCalculator(R, ℓₘₐₓ, 0)

"""
    sYlm!(Y, R, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s))
    sYlm!(Y, calc::sYlmCalculator, R; ℓₘᵢₙ=abs(s))
    sYlm!(Y, calc::sYlmCalculator, R, s; ℓₘᵢₙ=abs(s))

In-place version of [`sYlm`](@ref): fills `Y` with ``{}_sY_{ℓ,m}(R)`` in the canonical
ordering, and returns `Y`.  For a single spin weight `Y` is a vector, and its first
`Ysize(ℓₘᵢₙ, ℓₘₐₓ)` elements are written; for a range of them `Y` is a matrix, and the first
`length(s)` rows and `Ysize(ℓₘᵢₙ, ℓₘₐₓ)` columns are.

The element type of `Y` must be `Complex` of the rotor's own floating-point type — or, in the
calculator forms, of the calculator's.  A mismatch is an error rather than a silent
conversion: the type of the rotor is what decides the arithmetic, and `Y` is where the answer
lands, so disagreement between them means one of the two is not what the caller thinks it is.

The calculator forms reuse `calc` (which must have `Nᵣ = 1`, and whose `ℓₘₐₓ` is used), so
that repeated evaluation at many rotors allocates nothing; the calculator's rotor data are
replaced by `R`.  Without a spin weight, everything the calculator serves is written, which
for a calculator built for several is the matrix form.  With one, that spin weight is picked
out of [`spins`](@ref)`(calc)` and the result is a vector; calling this once per spin weight is
the other way to evaluate several of them at one point.

The indices may be half-integers, spelled as `Rational`s with denominator 2, on the terms
described under [`sYlm`](@ref): all of one kind, with `ℓₘᵢₙ` defaulting to the smallest ``|s|``
and the values including the phase ``i^{2s}``.  In the calculator forms the kind is already
fixed by the calculator, whose ``ℓ`` are integers or half-odd-integers according to how it was
constructed, and a spin weight or `ℓₘᵢₙ` of the other kind is refused with a message saying so.
"""
function sYlm!(Y::HarmonicValues, args...; kwargs...)
    # A `HarmonicValues` is filled by writing through its storage; the labels do not change,
    # so the container itself comes back.
    sYlm!(strided(Y), args...; kwargs...)
    Y
end
function sYlm!(
    Y::AbstractVecOrMat{<:Complex}, R::Rotor, ℓₘₐₓ::IndexSpelling, s::SpinSpelling;
    ℓₘᵢₙ=nothing
)
    sYlm_flat_helper!(Y, R, flat_indices(ℓₘₐₓ, s, ℓₘᵢₙ)...)
end
function sYlm_flat_helper!(
    Y::AbstractVecOrMat{<:Complex}, R::Rotor, ℓₘₐₓ::IT, s, ℓₘᵢₙ::IT
) where {IT<:HalfInteger}
    check_sYlm_args(ℓₘₐₓ, s, ℓₘᵢₙ)
    sYlm_helper!(Y, sYlmCalculator_helper(R, ℓₘₐₓ, s), R, s, ℓₘᵢₙ)
end
# In the calculator forms the calculator fixes the kind of index, so the arguments are
# normalized one at a time and the worker checks each against it.
function sYlm!(
    Y::AbstractVecOrMat{<:Complex}, calc::sYlmCalculator, R::Rotor; ℓₘᵢₙ=nothing
)
    sYlm_helper!(
        Y, calc, R, calc.s,
        ℓₘᵢₙ === nothing ? min_abs_spin(calc.s) : half_integer(ℓₘᵢₙ)
    )
end
function sYlm!(
    Y::AbstractVector{<:Complex}, calc::sYlmCalculator, R::Rotor, s::IndexSpelling;
    ℓₘᵢₙ=nothing
)
    let s = half_integer(s)
        sYlm_helper!(Y, calc, R, s, ℓₘᵢₙ === nothing ? abs(s) : half_integer(ℓₘᵢₙ))
    end
end

function sYlm_helper!(
    Y::AbstractVector, calc::HarmonicCalculator{IT}, R, s::HalfInteger,
    ℓₘᵢₙ::HalfInteger
) where {IT<:HalfInteger}
    check_index_kind(IT, s, "the spin weight s")
    check_index_kind(IT, ℓₘᵢₙ, "ℓₘᵢₙ")
    ℓₘₐₓ = SphericalFunctions.ℓₘₐₓ(calc)
    s = convert(IT, s)
    ℓₘᵢₙ = convert(IT, ℓₘᵢₙ)
    check_sYlm_args(ℓₘₐₓ, s, ℓₘᵢₙ)
    check_spin(calc, s)
    check_sYlm_calculator(calc, R)
    check_sYlm_eltype(Y, calc)
    let needed = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
        if length(Y) < needed
            error("Output vector has length $(length(Y)); at least $needed is needed.")
        end
    end
    set_rotors!(calc, R)
    iₛ = spin_index(calc, s)
    Yˡ = calc.Yˡ
    @inbounds for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ
        recurrence!(calc, ℓ)
        i₀ = Yindex(ℓ, -ℓ, ℓₘᵢₙ) - 1
        for j ∈ 1:2ℓ+1
            Y[i₀ + j] = Yˡ[1, iₛ, j]
        end
    end
    Y
end
function sYlm_helper!(
    Y::AbstractMatrix, calc::HarmonicCalculator{IT}, R, s::AbstractUnitRange,
    ℓₘᵢₙ::HalfInteger
) where {IT<:HalfInteger}
    check_index_kind(IT, ℓₘᵢₙ, "ℓₘᵢₙ")
    ℓₘₐₓ = SphericalFunctions.ℓₘₐₓ(calc)
    ℓₘᵢₙ = convert(IT, ℓₘᵢₙ)
    check_sYlm_args(ℓₘₐₓ, s, ℓₘᵢₙ)
    check_spin(calc, first(s))
    check_spin(calc, last(s))
    check_sYlm_calculator(calc, R)
    check_sYlm_eltype(Y, calc)
    let needed = Ysize(ℓₘᵢₙ, ℓₘₐₓ), n = length(s)
        if size(Y, 1) < n || size(Y, 2) < needed
            error(
                "Output matrix has size $(size(Y)); at least ($n, $needed) is needed, for "
                * "$n spin weights and $needed modes."
            )
        end
    end
    set_rotors!(calc, R)
    n = length(s)
    i₁ = spin_index(calc, first(s))
    Yˡ = calc.Yˡ
    @inbounds for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ
        recurrence!(calc, ℓ)
        i₀ = Yindex(ℓ, -ℓ, ℓₘᵢₙ) - 1
        for j ∈ 1:2ℓ+1
            for i ∈ 1:n
                Y[i, i₀ + j] = Yˡ[1, i₁ + (i - 1), j]
            end
        end
    end
    Y
end

# `R` is untyped because the same check serves a rotor and a bare angle θ; `check_rotor_type`
# is what decides whether the argument suits the calculator at all.
function check_sYlm_calculator(calc::HarmonicCalculator, R)
    if Nᵣ(calc) != 1
        error("sYlm! needs a calculator with Nᵣ=1; this one has Nᵣ=$(Nᵣ(calc)).")
    end
    check_rotor_type(calc, R)
    nothing
end

function check_sYlm_eltype(Y, calc::HarmonicCalculator{IT, RT, NT}) where {IT, RT, NT}
    if eltype(Y) !== NT
        wanted = NT <: Complex ? "Complex{$RT}" : "$RT"
        error(
            "This calculator works in $RT, so the output array's element type must be "
            * "$wanted, not $(eltype(Y))."
        )
    end
    nothing
end

"""
    sYlm_matrix(R⃗, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s))

The dense matrix of spin-weighted spherical harmonics ``{}_sY_{ℓ,m}(R_i)``, with rows indexed
by the rotors in `R⃗` and columns by the modes ``(ℓ, m)`` in the canonical ordering
`[(ℓ, m) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ]` (see [`Yindex`](@ref)).  Row `i` equals
`sYlm(R⃗[i], ℓₘₐₓ, s; ℓₘᵢₙ)`.  The computation is done in the rotors' own floating-point type;
to compute in another type, convert them.

Multiplying this matrix by a vector of mode weights synthesizes the corresponding
spin-weighted function at the rotors (`f = Y * f̃`); its (pseudo)inverse performs the
analysis.  It is computed with all rotors batched, which is the efficient path; for
``ℓₘₐₓ ≳ 64`` the matrix becomes large and the transforms in the "Transformations" section
of the documentation should be preferred.

`s` may also be an ascending range of spin weights.  The result is then three-dimensional,
indexed `[rotor, spin, mode]`, and is a stack of the matrices above rather than one of them:
`Y[:, i, :]` is the synthesis matrix of spin weight `s[i]`, and is what a product with mode
weights takes.  The spin index is an ordinary 1-based position, and `ℓₘᵢₙ` defaults to the
smallest ``|s|`` in the range.

The indices may be half-integers, spelled as `Rational`s with denominator 2, on the terms
described under [`sYlm`](@ref): all of one kind, with `ℓₘᵢₙ` defaulting to the smallest
``|s|`` and the values including the phase ``i^{2s}``.  The columns are then indexed by half-odd
``(ℓ, m)`` in the same canonical ordering, and [`Yindex`](@ref) locates them as before.
"""
function sYlm_matrix(
    R⃗::AbstractVector{<:Rotor}, ℓₘₐₓ::IndexSpelling, s::SpinSpelling; ℓₘᵢₙ=nothing
)
    sYlm_matrix_helper(sYlmCalculator_helper, R⃗, flat_indices(ℓₘₐₓ, s, ℓₘᵢₙ)...)
end
function sYlm_matrix_helper(
    make, R⃗::AbstractVector, ℓₘₐₓ::IT, s, ℓₘᵢₙ::IT
) where {IT<:HalfInteger}
    check_sYlm_args(ℓₘₐₓ, s, ℓₘᵢₙ)
    calc = make(R⃗, ℓₘₐₓ, s)
    fill_sYlm_matrix!(
        allocate_sYlm_matrix(number_type(calc), s, length(R⃗), ℓₘᵢₙ, ℓₘₐₓ),
        calc, s, ℓₘᵢₙ, ℓₘₐₓ
    )
end

function allocate_sYlm_matrix(::Type{T}, ::HalfInteger, Nᵣ, ℓₘᵢₙ, ℓₘₐₓ) where {T}
    Matrix{T}(undef, Nᵣ, Ysize(ℓₘᵢₙ, ℓₘₐₓ))
end
function allocate_sYlm_matrix(::Type{T}, s::AbstractUnitRange, Nᵣ, ℓₘᵢₙ, ℓₘₐₓ) where {T}
    Array{T, 3}(undef, Nᵣ, length(s), Ysize(ℓₘᵢₙ, ℓₘₐₓ))
end

function fill_sYlm_matrix!(Y::AbstractMatrix, calc, s::HalfInteger, ℓₘᵢₙ, ℓₘₐₓ)
    Nᵣ = SphericalFunctions.Nᵣ(calc)
    iₛ = spin_index(calc, s)
    Yˡ = calc.Yˡ
    @inbounds for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ
        recurrence!(calc, ℓ)
        i₀ = Yindex(ℓ, -ℓ, ℓₘᵢₙ) - 1
        for j ∈ 1:2ℓ+1
            for iᵣ ∈ 1:Nᵣ
                Y[iᵣ, i₀ + j] = Yˡ[iᵣ, iₛ, j]
            end
        end
    end
    Y
end
function fill_sYlm_matrix!(Y::AbstractArray{<:Any, 3}, calc, s::AbstractUnitRange, ℓₘᵢₙ, ℓₘₐₓ)
    Nᵣ = SphericalFunctions.Nᵣ(calc)
    n = length(s)
    i₁ = spin_index(calc, first(s))
    Yˡ = calc.Yˡ
    @inbounds for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ
        recurrence!(calc, ℓ)
        i₀ = Yindex(ℓ, -ℓ, ℓₘᵢₙ) - 1
        for j ∈ 1:2ℓ+1
            for i ∈ 1:n
                for iᵣ ∈ 1:Nᵣ
                    Y[iᵣ, i, i₀ + j] = Yˡ[iᵣ, i₁ + (i - 1), j]
                end
            end
        end
    end
    Y
end


### The real flavour: ₛλₗₘ(θ) = ₛYₗₘ(θ, 0) / i^{2s}
#
# Every one of these is the corresponding ₛYₗₘ form with `sλlmCalculator_helper` in place of
# `sYlmCalculator_helper`, and an angle in place of a rotor.  There is no separate machinery:
# the helpers above take the factory as an argument precisely so that this file does not
# acquire a second copy of the same loops.

"""
    sλlm(θ, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s))

The real functions ``{}_sλ_{ℓ,m}(θ) = {}_sY_{ℓ,m}(θ, 0) / i^{2s}`` for all ``ℓ ≤ ℓₘₐₓ``,
returned in the same [`HarmonicValues`](@ref) container as [`sYlm`](@ref) and indexed exactly
as described there — but holding real numbers rather than complex ones.  `θ` is one angle or
an `AbstractVector` of them, and `s` one spin weight or an ascending range.

This is what the ring-based transforms want, and what they used to extract from a complex
calculator by hand; the values are bit-for-bit the same.  For an integer spin weight the
division by ``i^{2s} = (-1)^s`` is simply the factor that ``{}_sY_{ℓ,m}`` already includes, so
``{}_sλ_{ℓ,m}`` is the harmonic at ``ϕ = 0`` as the literature writes it; for a half-odd spin
weight it removes the ``\\pm i`` that would otherwise make the value imaginary.  See
[`sλlmCalculator`](@ref) for the details, and for repeated evaluation.

The indices may be half-integers on exactly the terms described under [`sYlm`](@ref).
"""
function sλlm(θ::Real, ℓₘₐₓ::IndexSpelling, s::SpinSpelling; ℓₘᵢₙ=nothing)
    sYlm_helper(sλlmCalculator_helper, θ, flat_indices(ℓₘₐₓ, s, ℓₘᵢₙ)...)
end
function sλlm(θ::AbstractVector{<:Real}, ℓₘₐₓ::IndexSpelling, s::SpinSpelling; ℓₘᵢₙ=nothing)
    sYlm_batch_helper(sλlmCalculator_helper, θ, flat_indices(ℓₘₐₓ, s, ℓₘᵢₙ)...)
end

"""
    sλlm!(Y, θ, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s))
    sλlm!(Y, calc::sλlmCalculator, θ; ℓₘᵢₙ=abs(s))
    sλlm!(Y, calc::sλlmCalculator, θ, s; ℓₘᵢₙ=abs(s))
    sλlm!(Y::HarmonicValues, args...; kwargs...)

Write ``{}_sλ_{ℓ,m}(θ)`` into the existing real array `Y`, which must have the layout
[`sλlm`](@ref) would return and an element type matching the calculator's.  This is
[`sYlm!`](@ref) for the real flavour, and behaves identically in every other respect.
"""
function sλlm!(Y::HarmonicValues, args...; kwargs...)
    sλlm!(strided(Y), args...; kwargs...)
    Y
end
function sλlm!(
    Y::AbstractVecOrMat{<:Real}, θ::Real, ℓₘₐₓ::IndexSpelling, s::SpinSpelling; ℓₘᵢₙ=nothing
)
    let (ℓₘₐₓ, s, ℓₘᵢₙ) = flat_indices(ℓₘₐₓ, s, ℓₘᵢₙ)
        check_sYlm_args(ℓₘₐₓ, s, ℓₘᵢₙ)
        sYlm_helper!(Y, sλlmCalculator_helper(θ, ℓₘₐₓ, s), θ, s, ℓₘᵢₙ)
    end
end
function sλlm!(Y::AbstractVecOrMat{<:Real}, calc::sλlmCalculator, θ::Real; ℓₘᵢₙ=nothing)
    sYlm_helper!(
        Y, calc, θ, calc.s, ℓₘᵢₙ === nothing ? min_abs_spin(calc.s) : half_integer(ℓₘᵢₙ)
    )
end
function sλlm!(
    Y::AbstractVector{<:Real}, calc::sλlmCalculator, θ::Real, s::IndexSpelling; ℓₘᵢₙ=nothing
)
    let s = half_integer(s)
        sYlm_helper!(Y, calc, θ, s, ℓₘᵢₙ === nothing ? abs(s) : half_integer(ℓₘᵢₙ))
    end
end

"""
    sλlm_matrix(θ⃗, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s))

The dense real matrix ``{}_sλ_{ℓ,m}(θ_i)``, with rows indexed by the angles in `θ⃗` and columns
by the modes in the canonical ordering of [`Yindex`](@ref).  This is [`sYlm_matrix`](@ref) for
the real flavour; `s` may likewise be a range, in which case the result is three-dimensional
and indexed `[angle, spin, mode]`.
"""
function sλlm_matrix(
    θ⃗::AbstractVector{<:Real}, ℓₘₐₓ::IndexSpelling, s::SpinSpelling; ℓₘᵢₙ=nothing
)
    sYlm_matrix_helper(sλlmCalculator_helper, θ⃗, flat_indices(ℓₘₐₓ, s, ℓₘᵢₙ)...)
end
