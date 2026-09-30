"""
    HarmonicCalculator{IT, RT, NT, ST, S, B, FT, L}

Calculator producing the spin-weighted spherical harmonics ``{}_sY_{ℓ,m}`` (when `NT` is
`Complex{RT}`) or the real ``{}_sλ_{ℓ,m}`` (when `NT` is `RT`), for `Nᵣ` points at a time,
one ``ℓ`` at a time.  Use the constructors [`sYlmCalculator`](@ref) and
[`sλlmCalculator`](@ref).

- `IT` is the index type, `Int` or [`HalfOddInteger`](@ref).
- `RT` is the real type the calculator works in.
- `NT` is the number type of the harmonics, `Complex{RT}` or `RT`.
- `ST` is the storage type of the ``H`` wedge.
- `S` is the type of the spin weights served: the index type for one spin weight, or a
  `UnitRange` of it for a range of them.
- `B` is `true` exactly when the calculator was built from a vector of rotor data, and is
  what [`isbatched`](@ref) reads.
- `FT` is the real type in which the recurrence runs, which is `RT` itself unless `RT`
  holds derivatives, as a dual number does.
- `L` is `Nothing`, unless the calculator lifts the blocks of a calculator of the values of
  its rotors into blocks that hold derivatives, when it is the type of the data for that.


Internally this wraps an [`HCalculator`](@ref), which runs the recurrence that both subtypes
use, plus a buffer backing the block for the current ``ℓ``.  The phase tables `Z₊` and `Z₋`
are empty for the real subtype, which is the whole of the saving: the ``H`` recurrence is
real, and it is only the ``e^{-i(mα - sγ)}`` factor that ever made the result complex.


Because `S` and `B` decide the shape of the block, the type of the block is known at compile
time, and a loop over the blocks is inferrable.
"""
struct HarmonicCalculator{IT, RT<:Real, NT<:Union{RT, Complex{RT}}, ST, S, B, FT<:Real, L}
    # As for [`WignerCalculator`](@ref), the parameter `B` says whether the calculator was
    # built from a vector of rotor data (of any length), lifted into the type so that the
    # branch in `spin_row` and `spin_block` — and hence the type of the block that
    # `recurrence!` returns — is settled at compile time.  `S` does the same job for the
    # spin weights: it is the index type when the calculator was built for one of them and a
    # `UnitRange` of it when it was built for several, which is what decides whether a block
    # has a spin axis at all.
    H::HCalculator{IT, FT, ST}
    Yˡ::Array{NT, 3}  # [iᵣ, s, m] block for the current ℓ, using the leading m entries
    Z₊::Matrix{Complex{FT}}  # Z₊[iᵣ, k+1] = z₊^k for k ∈ 0:2ℓₘₐₓ
    Z₋::Matrix{Complex{FT}}  # Z₋[iᵣ, k+1] = z₋^k for k ∈ 0:2ℓₘₐₓ
    rotors::Vector{Quaternion{RT}}  # the rotors, or those of the points (θ, 0); empty for ₛλₗₘ
    s::S
    ℓ::Base.RefValue{IT}  # ℓ of the block currently in Yˡ; ℓₘᵢₙ-1 if none
    phases::Base.RefValue{Bool}  # false when the rotor data are angles θ (ϕ = γ = 0)
    lift::L  # `nothing`, or the calculator of the rotors' values and their generators
    # `materialize!` writes `Yˡ` and reads the power tables under `@inbounds`, for every
    # rotor of `H`, every spin weight served and every m up to ℓₘₐₓ, so the buffers must be
    # large enough for those; as for `HCalculator`, this checks them once, as they are
    # brought together.
    function HarmonicCalculator{IT, RT, NT, ST, S, B, FT, L}(
        H, Yˡ, Z₊, Z₋, rotors, s, ℓ, phases, lift
    ) where {IT, RT<:Real, NT<:Union{RT, Complex{RT}}, ST, S, B, FT<:Real, L}
        let n = Nᵣ(H), M = 2ℓₘₐₓ(H) + 1, K = NT <: Complex ? 2ℓₘₐₓ(H) + 1 : 0
            if !(
                size(Yˡ, 1) ≥ n && size(Yˡ, 2) ≥ nspins(s) && size(Yˡ, 3) ≥ M
                && size(Z₊, 1) ≥ n && size(Z₊, 2) ≥ K && size(Z₋, 1) ≥ n && size(Z₋, 2) ≥ K
                && length(rotors) == (NT <: Complex ? n : 0)
            )
                throw(DimensionMismatch(
                    "The buffers of a $(flavor_name(NT)) for Nᵣ=$n rotors, the spin weights "
                    * "$s and ℓₘₐₓ=$(ℓₘₐₓ(H)) are too small: the block has size $(size(Yˡ)), "
                    * "the power tables $(size(Z₊)) and $(size(Z₋)), which need a row for "
                    * "each rotor and at least $K columns, and there are $(length(rotors)) "
                    * "rotors."
                ))
            end
        end
        new{IT, RT, NT, ST, S, B, FT, L}(H, Yˡ, Z₊, Z₋, rotors, s, ℓ, phases, lift)
    end
end

"""
    sYlmCalculator(R, ℓₘₐₓ, s)
    sYlmCalculator(θ, ϕ, ℓₘₐₓ, s)

Calculator for the spin-weighted spherical harmonics ``{}_sY_{ℓ,m}`` for all ``ℓ ≤ ℓₘₐₓ``,
with elements of `Complex` of the rotor's own element type.  The spin weight `s` is either a
single value or an ascending unit range of them, such as `-2:2`; the calculator serves
exactly those, and works out for itself how much of the underlying recurrence they need.

The rotor is given at construction, so that the calculator is usable the moment it exists,
and so that the element type is the rotor's own: a `Rotor{Float32}` gives a `Float32`
calculator, and `floattype(calc)` reports it.  There is no argument to override that — to
compute in another type, convert the rotor.  Give `R` as an `AbstractVector` of rotors to
get a calculator that handles all `Nᵣ = length(R)` of them at once.  Later rotors are
supplied with [`set_R!`](@ref).  A single rotor may instead be given by the spherical
coordinates of its point: `sYlmCalculator(θ, ϕ, ℓₘₐₓ, s)` is
`sYlmCalculator(from_spherical_coordinates(θ, ϕ), ℓₘₐₓ, s)`.

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

For a calculator built from a vector of rotors or angles — of any length, even one — each
block gains a leading rotor index, so it is read as `ₛYₗ[iᵣ, m]` or `ₛYₗ[iᵣ, s, m]`.  The
same block is what [`recurrence!`](@ref) returns when the calculator is stepped by hand, and
`ₛYₗ[s, :]` picks out the row of one spin weight of a calculator built for several.  Each
block is a view into the calculator's storage, overwritten by the next step; `copy` it if it
must survive (keeping the natural indices), or `collect` it to get an ordinary 1-based
array.  Wherever ``ℓ < |s|`` the values are zero.

[`spins`](@ref) reports the range of spin weights served, and [`spin`](@ref) the single
value when there is only one.  Lengthening the range costs storage and arithmetic in
proportion to its length, but it is the largest ``|s|`` alone that decides the size of the
underlying Wigner ``H`` wedge: `sYlmCalculator(R, ℓₘₐₓ, 2)` and `sYlmCalculator(R, ℓₘₐₓ,
-2:2)` run precisely the same recurrence, and differ only in how much of it is read out.

Instead of rotors, the angle ``θ`` (or a vector of `Nᵣ` angles) may be given, at
construction or through [`set_θ!`](@ref).  This evaluates the harmonics at ``(θ, ϕ=0)``,
which is what the ring-based transforms need — though the values are still stored as complex
numbers, and a calculator that will only ever be used that way is better built as an
[`sλlmCalculator`](@ref), which stores them as reals.

A calculator is a mutable workspace, which every step and every setter overwrites, so one
calculator must not be used by two tasks at once; `similar(calc)` gives each task a
calculator of its own, holding the same rotor data.

The harmonics are defined by
```math
{}_sY_{ℓ,m}(𝐑) = (-1)^s \\sqrt{\\frac{2ℓ+1}{4π}}\\, \\overline{𝔇^{(ℓ)}_{m,-s}(𝐑)},
```
so that they are functions on the rotation group, and the harmonics at the spherical
coordinates ``(θ, ϕ)`` are those at `R = from_spherical_coordinates(θ, ϕ)`.  See the
"Conventions" section of the documentation.

# Half-integer indices

`ℓₘₐₓ` and the spin weights may all be half-integers, passed as `Rational`s with denominator
2 or as [`HalfOddInteger`](@ref)s, in which case ``ℓ``, ``m`` and ``s`` are all
half-integers.  A range is written the same way, as `-3//2:3//2`.  The prefactor ``(-1)^s``
is then ``\\pm i``; the principal branch ``(-1)^s ≡ e^{iπs} = i^{2s}`` settled in the
"Conventions" section is used, so the harmonics are not real multiples of
``\\overline{𝔇}``, and ``{}_sY_{ℓ,m}(θ, 0)`` is therefore not real (see
[`sλlmCalculator`](@ref), which divides that constant phase out).  The blocks are the same
containers as for integer indices — a [`DegreeBlock`](@ref) or [`DegreeBlockBatch`](@ref)
for one spin weight, and a [`SpinMatrix`](@ref) or [`SpinMatrixBatch`](@ref) for several —
indexed exactly as described above.  The flat interfaces [`sYlm`](@ref), [`sYlm!`](@ref) and
[`sYlm_matrix`](@ref) accept the same types, and lay the values out in the canonical
mode-weight ordering of [`Yindex`](@ref), which holds for half-odd indices exactly as it
does for integers.

See also [`sYlm`](@ref) and [`sYlm_matrix`](@ref) for simpler interfaces, and
[`DCalculator`](@ref).
"""
const sYlmCalculator{IT, RT, ST, S, B} =
    HarmonicCalculator{IT, RT, Complex{RT}, ST, S, B} where {IT, RT<:Real, ST, S, B}

"""
    sλlmCalculator(θ, ℓₘₐₓ, s)
    slambdalmCalculator(θ, ℓₘₐₓ, s)

Calculator for the real functions

```math
{}_sλ_{ℓ,m}(θ) = \\begin{cases}
    {}_sY_{ℓ,m}(θ, 0), & s ∈ ℤ, \\\\
    {}_sY_{ℓ,m}(θ, 0) \\big/ i^{2s}, & s ∈ ℤ + \\tfrac{1}{2},
\\end{cases}
```

for all ``ℓ ≤ ℓₘₐₓ``, with elements of the angle's own real type.  The first argument is the
angle ``θ``, or an `AbstractVector` of `Nᵣ` of them; later values are supplied with
[`set_θ!`](@ref).  Otherwise this behaves exactly like [`sYlmCalculator`](@ref) — the same
blocks, the same iteration, the same half-integer types — but stores real numbers rather
than complex ones, and allocates no phase tables at all.  `slambdalmCalculator` is an ASCII
alias of the same constructor.

Both are real.  For an integer spin weight the prefactor ``(-1)^s`` in the definition of
``{}_sY_{ℓ,m}`` is ``\\pm 1``, so ``{}_sY_{ℓ,m}(θ, 0)`` is already real, and ``{}_sλ_{ℓ,m}``
is just the harmonic evaluated at ``ϕ = 0``, exactly as the literature writes it.  For a
half-odd spin weight that prefactor is ``i^{2s} = \\pm i``, so ``{}_sY_{ℓ,m}(θ, 0)`` is
imaginary rather than real, and dividing that constant phase out is what leaves a real
function behind.  (Dividing by ``i^{2s}`` in both cases would give the wrong sign for odd
integer ``s``, where ``i^{2s} = -1``.)  The values are bit-for-bit those of an
[`sYlmCalculator`](@ref) at the same angles: the real part for an integer spin weight, and
``±`` the imaginary part for a half-odd one.

A `Rotor` is **not** accepted, here or through [`set_R!`](@ref): a rotor specifies the
angles ``α`` and ``γ``, whose phases a real calculator cannot represent.  Use an
[`sYlmCalculator`](@ref) for that.

See also [`sλlm`](@ref) and [`sλlm_matrix`](@ref) for simpler interfaces, and
[`dCalculator`](@ref), which stands in the same relation to [`DCalculator`](@ref).
"""
const sλlmCalculator{IT, RT, ST, S, B} =
    HarmonicCalculator{IT, RT, RT, ST, S, B, RT, Nothing} where {IT, RT<:Real, ST, S, B}
# A calculator of the real harmonics never lifts the blocks of another, since it is given
# angles, whose recurrence is differentiated as it runs (see `src/derivatives/kernels.jl`),
# so the last two parameters above are fixed, and the type is concrete once the others are.

# The spin-weight argument is typed `IndexOrRange` rather than left open, so that a call
# whose arguments are in the wrong order — `sYlmCalculator(3, 1, Float64)`, say — is still
# the `MethodError` it should be rather than a puzzling complaint about the indices.
@index_methods function sYlmCalculator(R, ℓₘₐₓ::IndexType, s::IndexOrRange)
    sYlmCalculator_helper(R, ℓₘₐₓ, s)
end
@index_methods function sYlmCalculator(θ::Real, ϕ::Real, ℓₘₐₓ::IndexType, s::IndexOrRange)
    sYlmCalculator(from_spherical_coordinates(θ, ϕ), ℓₘₐₓ, s)
end
function sYlmCalculator_helper(R, ℓₘₐₓ::IT, s) where {IT<:IntegerHalf}
    check_spin_range(s)
    # See the note on `DCalculator`: the element type reaches `allocate_Y` as a type,
    # not as a value, so that the concrete result type is settled at compile time.
    RT = rotor_basetype(R)
    set_rotors!(allocate_Y(IT, RT, Complex{RT}, ℓₘₐₓ, s, nrotors(R), batched_data(R)), R)
end

@index_methods function sλlmCalculator(θ, ℓₘₐₓ::IndexType, s::IndexOrRange)
    sλlmCalculator_helper(θ, ℓₘₐₓ, s)
end
function sλlmCalculator_helper(θ, ℓₘₐₓ::IT, s) where {IT<:IntegerHalf}
    check_spin_range(s)
    RT = rotor_basetype(θ)
    set_rotors!(allocate_Y(IT, RT, RT, ℓₘₐₓ, s, nrotors(θ), batched_data(θ)), θ)
end

# A range of spin weights reaches the helpers above as a `UnitRange` (a range whose step is
# not 1 is refused before, with an explanation of how to write it), but a `UnitRange` may
# still be empty, which is also how a range written downward with a unit step, such as
# `2:-2`, comes out.  The two are a single mistake with a single remedy, and naming that
# remedy is more use than naming the emptiness.
check_spin_range(::IntegerHalf) = nothing
function check_spin_range(s::AbstractUnitRange)
    if isempty(s)
        throw(ArgumentError(
            "The range of spin weights $s runs downward or is empty.  A range runs from its "
            * "lower limit to its upper one, so 3//2:-1:-3//2 is written -3//2:3//2."
        ))
    end
    nothing
end

# The largest and smallest ``|s|`` among the spin weights a calculator serves, and how many of
# them there are.  These are written from the endpoints rather than as `maximum(abs, s)` and
# `minimum(abs, s)` because a reduction over a range of `HalfOddInteger`s reaches for an
# identity element that the type deliberately lacks.  Note that `min_abs_spin` is not simply
# the smaller of the two magnitudes: a range that straddles zero contains the smallest index
# of its kind, which is 0 for integers and 1/2 for half-odd-integers.
max_abs_spin(s::IntegerHalf) = abs(s)
max_abs_spin(s::AbstractUnitRange) = max(abs(first(s)), abs(last(s)))
min_abs_spin(s::IntegerHalf) = abs(s)
function min_abs_spin(s::AbstractUnitRange{IT}) where {IT<:IntegerHalf}
    first(s) ≤ 0 ≤ last(s) ? ℓₘᵢₙ(IT) : min(abs(first(s)), abs(last(s)))
end
nspins(::IntegerHalf) = 1
nspins(s::AbstractUnitRange) = length(s)

# Allocate the buffers without touching them.  PRIVATE: see the note on `allocate_H`.  Here
# the uninitialized state includes the `phases` flag, which decides whether `Z₊` and `Z₋`
# are ever read, so this must not escape without a `set_rotors!` or a full buffer copy.
function allocate_Y(
    ::Type{IT}, ::Type{RT}, ::Type{NT}, ℓₘₐₓ::IT, s::S, Nᵣ::Int, ::Val{B}
) where {IT<:IntegerHalf, RT<:Real, NT<:Union{RT, Complex{RT}}, S, B}
    sₕ = max_abs_spin(s)
    if sₕ > ℓₘₐₓ
        throw(ArgumentError(
            "The spin weights $s include |s|=$sₕ, which exceeds ℓₘₐₓ=$ℓₘₐₓ; there are no "
            * "such harmonics."
        ))
    end
    # The wedge is sized by the largest |s| alone, and cannot be narrowed further for a single
    # spin weight, tempting though that looks.  `materialize!` asks for H[m, -s] with m over
    # the whole of -ℓ:ℓ.  Where |m| > |s| the stored representative is the row -s or +s
    # according to the sign of m; where |m| ≤ |s| the symmetries offer only the row ±m, and
    # m runs over both signs.  So a single spin weight touches every row of -|s|:|s| just as
    # the full range would, and |s| is the smallest wedge that can serve it.
    Yˡ = Array{NT, 3}(undef, Nᵣ, nspins(s), 2ℓₘₐₓ + 1)
    rotors = Vector{Quaternion{RT}}(undef, NT <: Complex ? Nᵣ : 0)
    ℓ = Ref{IT}(ℓₘᵢₙ(IT) - 1)  # typed explicitly because the field is a `RefValue{IT}`
    # A calculator of ₛYₗₘ whose real type holds derivatives lifts the blocks of a
    # calculator of its rotors' values, whose `H`, power tables, and `phases` flag it
    # shares, as a `WignerCalculator` does; see `allocate_W`.
    if NT <: Complex && value_type(RT) !== RT
        let inner = allocate_Y(IT, value_type(RT), Complex{value_type(RT)}, ℓₘₐₓ, s, Nᵣ, Val(B))
            lift = allocate_lift(RT, inner, Nᵣ)
            HarmonicCalculator{IT, RT, NT, typeof(parent(inner.H.Hˡ)), S, B, float_type(RT), typeof(lift)}(
                inner.H, Yˡ, inner.Z₊, inner.Z₋, rotors, s, ℓ, inner.phases, lift
            )
        end
    else
        H = allocate_H(IT, RT, ℓₘₐₓ, sₕ, Nᵣ)
        # The phase tables are what a real calculator does not have: with `K = 0` they are
        # empty rather than merely unread, so the ``{}_sλ_{ℓ,m}`` flavor allocates nothing
        # for them.  This is the same trick `allocate_W` uses to separate ``𝔇`` from ``d``,
        # and the layout, with the rotor index first, is the same too.
        K = NT <: Complex ? 2ℓₘₐₓ + 1 : 0
        Z₊ = Matrix{Complex{RT}}(undef, Nᵣ, K)
        Z₋ = Matrix{Complex{RT}}(undef, Nᵣ, K)
        HarmonicCalculator{IT, RT, NT, typeof(parent(H.Hˡ)), S, B, RT, Nothing}(
            H, Yˡ, Z₊, Z₋, rotors, s, ℓ, Ref(NT <: Complex), nothing
        )
    end
end

# See the notes on `similar(::WignerCalculator)`, whose docstring covers these methods too,
# for why the assertion is here and why the rotor data is copied rather than re-derived.
# The `phases` flag is part of that data: it records whether the calculator was given rotors
# or bare angles, and hence whether `Z₊` and `Z₋` hold anything at all.
function Base.similar(
    c::HarmonicCalculator{IT, RT, NT, ST, S, B, FT, L}
) where {IT, RT, NT, ST, S, B, FT<:Real, L}
    c′ = allocate_Y(IT, RT, NT, ℓₘₐₓ(c), c.s, Nᵣ(c), Val(B))::HarmonicCalculator{IT, RT, NT, ST, S, B, FT, L}
    copy_rotor_state!(c′, c)
end
# A calculator that lifts the blocks of another shares that calculator's `H`, power tables,
# and `phases` flag, which are therefore copied once for both; its own
# rotors and generators are copied, and then the rest of the other calculator's state.
function copy_rotor_state!(c′::HarmonicCalculator, c::HarmonicCalculator)
    copy_rotor_data!(c′.H, c.H)
    copyto!(c′.Z₊, c.Z₊)
    copyto!(c′.Z₋, c.Z₋)
    c′.phases[] = c.phases[]
    copyto!(c′.rotors, c.rotors)
    if c.lift !== nothing
        copyto!(c′.lift.G, c.lift.G)
        copy_rotor_state!(c′.lift.inner, c.lift.inner)
    end
    c′
end
function Base.similar(
    c::HarmonicCalculator{IT, RT, NT, ST, S, B, FT, L}, R
) where {IT, RT, NT, ST, S, B, FT<:Real, L}
    check_rotor_count(c, R)
    check_rotor_type(c, R)
    set_rotors!(
        allocate_Y(IT, RT, NT, ℓₘₐₓ(c), c.s, Nᵣ(c), Val(B))::HarmonicCalculator{IT, RT, NT, ST, S, B, FT, L}, R
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
flavor_name(::Type{<:Complex}) = "sYlmCalculator"
flavor_name(::Type{<:Real}) = "sλlmCalculator"
Nᵣ(c::HarmonicCalculator) = Nᵣ(c.H)
isbatched(::HarmonicCalculator{IT, RT, NT, ST, S, B}) where {IT, RT, NT, ST, S, B} = B

# `spins` answers in the same currency however the calculator was built, so that code which
# does not care can loop over it either way; `spin` exists only where there is a single
# value to name, and a calculator built for several gives a `MethodError` rather than a
# value that would have to be wrong.
spins(c::HarmonicCalculator{IT, RT, NT, ST, S}) where {IT, RT, NT, ST, S<:IntegerHalf} = c.s:c.s
spins(c::HarmonicCalculator{IT, RT, NT, ST, S}) where {IT, RT, NT, ST, S<:AbstractUnitRange} = c.s
spin(c::HarmonicCalculator{IT, RT, NT, ST, S}) where {IT, RT, NT, ST, S<:IntegerHalf} = c.s

# As for the Wigner calculators, a batched calculator says so.
function Base.show(io::IO, c::HarmonicCalculator{IT, RT, NT}) where {IT, RT, NT}
    print(
        io,
        "$(flavor_name(NT)){$IT, $RT} for ℓₘₐₓ=$(ℓₘₐₓ(c)), s=$(c.s), Nᵣ=$(Nᵣ(c))",
        isbatched(c) ? ", batched" : "",
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
`Z₊`, `Z₋`, and the flag recording whether rotors or bare angles were given — is
deliberately *not* touched, so `recurrence!(c, ℓ)` still has everything it needs, exactly as
for [`HCalculator`](@ref).  Useful for testing that no uninitialized storage is ever read.
"""
function Base.fill!(c::HarmonicCalculator{IT, RT, NT}, v::Number) where {IT, RT, NT}
    fill!(c.H, real(v))
    fill!(c.Yˡ, convert(NT, v))
    c.ℓ[] = ℓₘᵢₙ(IT) - 1
    c
end


### Rotor data: full rotors (as for 𝔇), or angles θ meaning (θ, ϕ=0)

# The `HCalculator` validates everything before it replaces anything, so if it refuses the
# angles this calculator is left exactly as it was, and its own state is reset only once the
# new data are in place.  A calculator of ₛYₗₘ also keeps the rotors of the points (θ, 0),
# which are what the rules for automatic differentiation read (see
# `src/derivatives/kernels.jl`).
#
# As for a `WignerCalculator`, setting the rotor data is two steps: the rotors are copied by
# `store_rotors!` or `store_point_rotors!`, and everything else is computed by
# `set_rotor_data!`, which the extensions for Enzyme and Mooncake declare, for ₛYₗₘ, to have
# no derivatives.  A calculator that lifts the blocks of another gives it the values of its
# rotor data, and computes its own generators.
function set_rotors!(
    c::HarmonicCalculator{IT, RT, NT, ST, S, B, FT, Nothing}, θ::Union{Real, AbstractVector{<:Real}}
) where {IT, RT<:Real, NT, ST, S, B, FT<:Real}
    set_rotor_data!(c, θ)
    NT <: Complex && store_point_rotors!(c.rotors, θ)
    c
end
function set_rotor_data!(c::HarmonicCalculator{IT}, θ::Union{Real, AbstractVector{<:Real}}) where {IT}
    set_rotors!(c.H, θ)
    c.phases[] = false
    c.ℓ[] = ℓₘᵢₙ(IT) - 1
    nothing
end
function set_rotors!(
    c::HarmonicCalculator{IT, RT, NT, ST, S, B, FT, <:Lift}, θ::AbstractVector{<:Real}
) where {IT, RT<:Real, NT, ST, S, B, FT<:Real}
    set_rotors!(c.lift.inner, LiftedValues(θ))
    store_point_rotors!(c.rotors, θ)
    set_generators!(c.lift, true, c.rotors)
    c.ℓ[] = ℓₘᵢₙ(IT) - 1
    c
end
function set_rotors!(
    c::HarmonicCalculator{IT, RT, NT, ST, S, B, FT, <:Lift}, θ::Real
) where {IT, RT<:Real, NT, ST, S, B, FT<:Real}
    check_single_rotor(c)
    set_rotors!(c, @SVector [θ])
end

# The rotors of the points (θ, 0).
function store_point_rotors!(rotors::AbstractVector{Quaternion{T}}, θ) where {T}
    @inbounds for i ∈ eachindex(rotors)
        θᵢ = θ[i]
        rotors[i] = as_quaternion(T, from_spherical_coordinates(θᵢ, zero(θᵢ)))
    end
    rotors
end

function set_rotors!(
    c::HarmonicCalculator{IT, RT, Complex{RT}, ST, S, B, FT, Nothing}, R::AbstractVector{<:RotorLike}
) where {IT, RT<:Real, ST, S, B, FT<:Real}
    # The loops write the calculator's 1-based buffers at the input's own indices, under
    # `@inbounds`, so an offset vector would write outside them.
    Base.require_one_based_indexing(R)
    check_rotor_length(c, R)
    store_rotors!(c.rotors, R)
    set_rotor_data!(c, R)
    c
end
function set_rotor_data!(
    c::HarmonicCalculator{IT, RT, Complex{RT}}, R::AbstractVector{<:RotorLike}
) where {IT, RT<:Real}
    # As in `set_rotors!(::HCalculator, …)`: every rotor is acceptable, and the results are
    # marked invalid before the first one is replaced.
    c.H.axes_valid[] = false
    c.ℓ[] = ℓₘᵢₙ(IT) - 1
    @inbounds for i ∈ eachindex(R)
        z₊, z₋ = store_rotor!(c.H, i, R[i])
        complex_powers!(view(c.Z₊, i, :), z₊)
        complex_powers!(view(c.Z₋, i, :), z₋)
    end
    c.phases[] = true
    nothing
end
function set_rotors!(
    c::HarmonicCalculator{IT, RT, Complex{RT}, ST, S, B, FT, <:Lift}, R::AbstractVector{<:RotorLike}
) where {IT, RT<:Real, ST, S, B, FT<:Real}
    Base.require_one_based_indexing(R)
    check_rotor_length(c, R)
    store_rotors!(c.rotors, R)
    set_rotors!(c.lift.inner, LiftedValues(c.rotors))
    set_generators!(c.lift, true, c.rotors)
    c.ℓ[] = ℓₘᵢₙ(IT) - 1
    c
end
function set_rotors!(c::sYlmCalculator{IT, RT}, R::RotorLike) where {IT, RT<:Real}
    check_single_rotor(c)
    set_rotors!(c, @SVector [R])
end
# A real calculator refuses rotors rather than silently dropping the phases they specify,
# exactly as `set_β!` refuses them for the real `d` matrices.  The check is a method rather
# than a branch so that the refusal happens at the outermost call, naming the type the
# caller actually has.
function set_rotors!(c::sλlmCalculator, R::Union{RotorLike, AbstractVector{<:RotorLike}})
    throw(ArgumentError(
        "An sλlmCalculator evaluates at (θ, ϕ=0) and stores real numbers, so it cannot take "
        * "a rotor, whose α and γ angles are phases it has nowhere to put.  Give the angle θ "
        * "instead, or use an sYlmCalculator."
    ))
end
# One catch-all rather than two, so that it stays strictly less specific than every method
# above: a pair of them, keyed on the flavor, would be ambiguous with the angle methods,
# which name the argument type but not the calculator's.  The branch is on a type parameter,
# so only the message costs anything, and only on the way to an error.
function set_rotors!(c::HarmonicCalculator{IT, RT, NT}, R) where {IT, RT<:Real, NT}
    if NT <: Complex
        throw(ArgumentError(
            "An sYlmCalculator needs rotors (as `Rotor`s or `Quaternion`s) or angles θ::Real — one, or an "
            * "AbstractVector of $(Nᵣ(c)) of either — not $(typeof(R))."
        ))
    else
        throw(ArgumentError(
            "An sλlmCalculator needs angles θ::Real — one, or an AbstractVector of "
            * "$(Nᵣ(c)) of them — not $(typeof(R))."
        ))
    end
end


### Driver

function recurrence!(c::HarmonicCalculator, R, ℓ)
    check_ℓ(c.H, ℓ)
    check_rotor_type(c, R)  # as `set_R!` and `set_θ!` do, rather than silently converting
    set_rotors!(c, R)
    recurrence!(c, ℓ)
end
function recurrence!(c::HarmonicCalculator{IT}, ℓ) where {IT}
    let ℓ = calculator_index(IT, ℓ, "ℓ")
        compute_block!(c, ℓ)
        current_block(c, ℓ)
    end
end

# Compute the spin rows `is` of the block of degree ℓ into `Y[:, :, j₀ .+ (1:2ℓ+1)]`, where
# `Y` is the calculator's own block `Yˡ` with `j₀ = 0`, or another array of the same layout,
# as `sYlm_matrix` gives it; otherwise as for the `WignerCalculator` method, which describes
# the role of this function.  The destination is an array and an offset, rather than a view,
# so that the rules for it see an ordinary array.
compute_block!(c::HarmonicCalculator, ℓ) = compute_block!(c, ℓ, Base.OneTo(nspins(c.s)), c.Yˡ, 0)
function compute_block!(
    c::HarmonicCalculator{IT, RT, NT, ST, S, B, FT, Nothing}, ℓ::IT, is, Y, j₀::Int
) where {IT, RT, NT, ST, S, B, FT<:Real}
    recurrence!(c.H, ℓ)
    if Y === c.Yˡ && j₀ == 0
        materialize!(c, ℓ, is, c.Yˡ)
    else
        materialize!(c, ℓ, is, view(Y, :, :, (j₀ + 1):(j₀ + Int(2ℓ) + 1)))
    end
    c
end
function compute_block!(
    c::HarmonicCalculator{IT, RT, NT, ST, S, B, FT, <:Lift}, ℓ::IT, is, Y, j₀::Int
) where {IT, RT, NT, ST, S, B, FT<:Real}
    inner = c.lift.inner
    compute_block!(inner, ℓ, is, inner.Yˡ, 0)
    lift!(c, ℓ, is, Y, j₀)
    c.ℓ[] = Y === c.Yˡ && j₀ == 0 && length(is) == nspins(c.s) ? ℓ : ℓₘᵢₙ(IT) - 1
    c
end

# The block for the ``ℓ`` just computed: one spin weight's row when the calculator serves a
# single spin weight, and the whole spin axis when it serves a range.  Which of the two is a
# type parameter, so the choice is made at compile time and each method has one concrete
# return type.
function current_block(
    c::HarmonicCalculator{IT, RT, NT, ST, S}, ℓ::IT
) where {IT, RT, NT, ST, S<:IntegerHalf}
    spin_row(c, ℓ, 1)
end
function current_block(
    c::HarmonicCalculator{IT, RT, NT, ST, S}, ℓ::IT
) where {IT, RT, NT, ST, S<:AbstractUnitRange}
    spin_block(c, ℓ)
end

# (-1)^s, as e^{iπs} = i^{2s}; for integer s this is ±1.
minus_one_to_the(s::Integer) = ifelse(iseven(s), 1, -1)

# i^k for an integer k, exactly.  With k = 2s this is the principal branch of (-1)^s chosen
# by the conventions (see `docs/src/30-conventions/02-details.md`, "Half-integer indices"):
# ±1 for integer s, and ±i for half-integer s.
@inline function im_power(::Type{RT}, k::Integer) where {RT<:Real}
    let o = one(RT), z = zero(RT)
        (Complex{RT}(o, z), Complex{RT}(z, o), Complex{RT}(-o, z), Complex{RT}(z, -o))[mod(k, 4) + 1]
    end
end

# Index of spin weight s along the second dimension of c.Yˡ.  For a calculator built for a
# single spin weight this is 1, whatever that weight is.
@inline spin_index(c::HarmonicCalculator, s) = Int(s - first(spins(c))) + 1

# The coefficient (-1)^s ϵ(m) ϵ(s) σ √((2ℓ+1)/4π) of Hₘ,₋ₛ in ₛYₗₘ, in whichever number type
# `NT` the calculator stores.  For integer indices it is real whatever `NT` is, and stays
# real so that the integer path is bit-for-bit what it always was.  For half-odd indices
# (-1)^s = i^{2s} is ±i, which a complex calculator includes and a real one divides out —
# that division being the whole of the difference between ₛYₗₘ(θ,0) and ₛλₗₘ(θ).
#
# The prefactor is any `Real` rather than an `RT`: for ReverseDiff's number type, which
# records where a value came from as a type parameter, √((2ℓ+1)/4π) computed in `RT` is not
# itself an `RT`.
#
# The two half-odd coefficients differ by exactly the factor i^{2s}, and `Complex * Real` is
# computed componentwise, so the real flavor's value is bit-for-bit ±imag of the complex
# flavor's, which the regression test asserts with `==` rather than `≈`.
@inline function sYlm_coefficient(
    ::Type{NT}, ::Type{RT}, σ::Int, m::IT, s::IT, prefactor::Real
) where {NT, IT<:Integer, RT}
    convert(RT, σ * ϵ(m) * ϵ(s) * minus_one_to_the(s)) * prefactor
end
@inline function sYlm_coefficient(
    ::Type{NT}, ::Type{RT}, σ::Int, m::IT, s::IT, prefactor::Real
) where {NT<:Complex, IT<:HalfOddInteger, RT}
    (convert(RT, σ * ϵ(m) * ϵ(s)) * prefactor) * im_power(RT, 2s)
end
@inline function sYlm_coefficient(
    ::Type{NT}, ::Type{RT}, σ::Int, m::IT, s::IT, prefactor::Real
) where {NT<:Real, IT<:HalfOddInteger, RT}
    convert(RT, σ * ϵ(m) * ϵ(s)) * prefactor
end

# Write ₛYₗₘ for the spin weights at the positions `is` of `spins(c)` — by default all of
# them — and every m ∈ -ℓ:ℓ into c.Yˡ[iᵣ, s, m].  From the definition ₛYₗₘ = (-1)^s
# √((2ℓ+1)/4π) conj(𝔇ₘ,₋ₛ), with 𝔇ₘ,₋ₛ = ϵ(m) ϵ(s) Hₘ,₋ₛ e^{-i(mα - sγ)}, and e^{i(mα -
# sγ)} = z₊^(m-s) z₋^(m+s).  As in the Wigner `materialize!`, nothing here depends on
# whether the indices are integers or half-odd-integers.
#
# The element H[m, -s] is read, as `wedge_source` would read it, from one of three places:
#
#     m < -|s|    H[s, -m]                        along the row s
#     |m| ≤ |s|   through `wedge_source`
#     m > |s|     σ H[-s, m], σ = transpose_sign(m, -s)   along the row -s
#
# For one rotor, or one spin weight, each spin weight is written in turn, and the two outer
# parts are runs along a row of the wedge at a fixed stride, which is the cheapest order
# when the index arithmetic is what costs most.  With several of each the block is written
# in the order of its storage instead, a mode at a time, since the arithmetic is then spread
# over the rotors and it is the traffic to memory that costs most.  Both orders compute
# every element from the same source by the same expression.
#
# Only the rows `is` are written, and the block is marked as held only when every row was: a
# block with some rows left from another ℓ, or never written, must not be handed out as the
# block of this one.  `recurrence!` writes every row, and the flat interfaces that need
# fewer write only those (see `spin_row!` and `fill_sYlm!`).  The values may also be written
# into another array `Yˡ` of the same layout, [iᵣ, s, m] with m from -ℓ, as `sYlm_matrix`
# does; the calculator's own block is then not written at all.
function materialize!(
    c::HarmonicCalculator{IT, RT, NT}, ℓ::IT,
    is::AbstractUnitRange{Int}=Base.OneTo(nspins(c.s)), Yˡ::AbstractArray{NT, 3}=c.Yˡ
) where {IT, RT, NT}
    let H = c.H.Hˡ, Z₊ = c.Z₊, Z₋ = c.Z₋, Nᵣ = Nᵣ(c), Hp = parent(H)
        if H.ℓ != ℓ
            error("The H wedge holds ℓ=$(H.ℓ), but ℓ=$ℓ was requested.")
        end
        if !(
            1 ≤ first(is) && last(is) ≤ nspins(c.s) && size(Yˡ, 1) ≥ Nᵣ
            && size(Yˡ, 2) ≥ last(is) && size(Yˡ, 3) ≥ 2ℓ + 1
        )
            error(
                "Spin positions $is requested of a calculator holding $(nspins(c.s)), into "
                * "an array of size $(size(Yˡ)) for ℓ=$ℓ and Nᵣ=$Nᵣ."
            )
        end
        prefactor = √((2ℓ + 1) / (4 * RT(π)))
        srange = spins(c)
        W = m′ₘₐₓ(H)
        m′ₘᵢₙw = m′ₘᵢₙ(H)
        ri = row_index(H)
        # `NT <: Complex` is a compile-time constant, so a real calculator never compiles
        # the phase branch at all — which is what lets `Z₊` and `Z₋` be empty rather than
        # merely unread.
        phases = NT <: Complex && c.phases[]
        # Wherever ℓ < |s| the values are zero.
        @inbounds for i ∈ is
            if abs(srange[i]) > ℓ
                for j ∈ 1:2ℓ+1
                    for iᵣ ∈ 1:Nᵣ
                        Yˡ[iᵣ, i, j] = 0
                    end
                end
            end
        end
        if Nᵣ == 1 || length(is) == 1
            @inbounds for i ∈ is
                s = srange[i]
                a = abs(s)
                a > ℓ && continue
                # m < -|s|: H[m, -s] = H[s, -m]
                r = ri[(s - m′ₘᵢₙw) + 1] - 1
                for m ∈ -ℓ:(-a - 1)
                    coefficient = sYlm_coefficient(NT, RT, 1, m, s, prefactor)
                    materialize_element!(
                        Yˡ, Hp, Z₊, Z₋, Nᵣ, i, Int(m + ℓ) + 1, r + Nᵣ * Int(-m - a),
                        coefficient, m - s, m + s, phases
                    )
                end
                # |m| ≤ |s|
                for m ∈ -a:a
                    a′, b′, σ = wedge_source(m, -s, W)
                    coefficient = sYlm_coefficient(NT, RT, σ, m, s, prefactor)
                    materialize_element!(
                        Yˡ, Hp, Z₊, Z₋, Nᵣ, i, Int(m + ℓ) + 1, wedge_offset(H, a′, b′, m′ₘᵢₙw),
                        coefficient, m - s, m + s, phases
                    )
                end
                # m > |s|: H[m, -s] = σ H[-s, m]
                r = ri[(-s - m′ₘᵢₙw) + 1] - 1
                for m ∈ (a + 1):ℓ
                    σ = transpose_sign(m, -s)
                    coefficient = sYlm_coefficient(NT, RT, σ, m, s, prefactor)
                    materialize_element!(
                        Yˡ, Hp, Z₊, Z₋, Nᵣ, i, Int(m + ℓ) + 1, r + Nᵣ * Int(m - a),
                        coefficient, m - s, m + s, phases
                    )
                end
            end
        else
            @inbounds for m ∈ -ℓ:ℓ
                j = Int(m + ℓ) + 1
                for i ∈ is
                    s = srange[i]
                    a = abs(s)
                    a > ℓ && continue
                    if m < -a
                        σ = 1
                        offset = ri[(s - m′ₘᵢₙw) + 1] - 1 + Nᵣ * Int(-m - a)
                    elseif m > a
                        σ = transpose_sign(m, -s)
                        offset = ri[(-s - m′ₘᵢₙw) + 1] - 1 + Nᵣ * Int(m - a)
                    else
                        a′, b′, σ = wedge_source(m, -s, W)
                        offset = wedge_offset(H, a′, b′, m′ₘᵢₙw)
                    end
                    coefficient = sYlm_coefficient(NT, RT, σ, m, s, prefactor)
                    materialize_element!(
                        Yˡ, Hp, Z₊, Z₋, Nᵣ, i, j, offset, coefficient, m - s, m + s, phases
                    )
                end
            end
        end
    end
    c.ℓ[] = Yˡ === c.Yˡ && length(is) == nspins(c.s) ? ℓ : ℓₘᵢₙ(IT) - 1
    c
end

# One element (s, m) of the block, at position (i, j) of `Yˡ`, for every rotor, from the
# wedge element whose first rotor is at `offset + 1` of `Hp`.
@inline function materialize_element!(
    Yˡ::AbstractArray{NT}, Hp, Z₊, Z₋, Nᵣ, i, j, offset, coefficient, k₊, k₋, phases::Bool
) where {NT}
    # As in the Wigner `materialize_element!`, one rotor is written out, and otherwise the
    # signs of the powers are settled outside the loop; where both are negative the phase is
    # conj(a) * conj(b).
    if NT <: Complex && phases
        if Nᵣ == 1
            @inbounds Yˡ[1, i, j] = (
                coefficient * Hp[offset + 1] * (zpower(Z₊, 1, k₊) * zpower(Z₋, 1, k₋))
            )
        else
            with_power_signs(k₊, k₋) do n₊, n₋
                @inbounds @simd for iᵣ ∈ 1:Nᵣ
                    phase = zpower(Z₊, iᵣ, k₊, n₊) * zpower(Z₋, iᵣ, k₋, n₋)
                    Yˡ[iᵣ, i, j] = coefficient * Hp[offset + iᵣ] * phase
                end
            end
        end
    elseif Nᵣ == 1
        @inbounds Yˡ[1, i, j] = coefficient * Hp[offset + 1]
    else
        @inbounds @simd for iᵣ ∈ 1:Nᵣ
            Yˡ[iᵣ, i, j] = coefficient * Hp[offset + iᵣ]
        end
    end
    nothing
end

# Step the calculator to ℓ and return the row of the block for the spin weight at position
# `iₛ` of `spins(c)`, writing that row alone.  This is for callers that read one spin weight
# out of a calculator built for several, which would otherwise assemble every row at every ℓ
# only to use one of them.  Since the rest of the block is not written, the calculator does
# not count the block as held (see `materialize!`), and only the row returned may be read.
function spin_row!(c::HarmonicCalculator{IT}, ℓ, iₛ::Int) where {IT}
    let ℓ = calculator_index(IT, ℓ, "ℓ")
        compute_block!(c, ℓ, iₛ:iₛ, c.Yˡ, 0)
        spin_row(c, ℓ, iₛ)
    end
end

function check_spin(c::HarmonicCalculator, s)
    let sr = spins(c)
        if !(first(sr) ≤ s ≤ last(sr))
            throw(ArgumentError(
                "This calculator serves the spin weights $(c.s), so s=$s is not among them."
            ))
        end
    end
    nothing
end

# One spin weight's row, selected by its position `i` in the calculator's storage, and the
# whole block of every spin weight.  `isbatched(c)` reads a type parameter, so every branch
# is resolved at compile time.  The same containers are returned for both kinds of index.
function spin_row(c::HarmonicCalculator{IT}, ℓ::IT, i::Int) where {IT<:IntegerHalf}
    let mr = -ℓ:ℓ
        if isbatched(c)
            DegreeBlockBatch(view(c.Yˡ, :, i, 1:length(mr)), ℓ; mₘₐₓ=last(mr), mₘᵢₙ=first(mr))
        else
            DegreeBlock(view(c.Yˡ, 1, i, 1:length(mr)), ℓ; mₘₐₓ=last(mr), mₘᵢₙ=first(mr))
        end
    end
end

function spin_block(c::HarmonicCalculator{IT}, ℓ::IT) where {IT<:IntegerHalf}
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
# Each public function here is defined with `@index_methods`, so that it accepts an index of
# any type a caller may reasonably write — an `Int`, a `HalfOddInteger`, or a `Rational`
# with denominator 2 — and a unit range of them, while its body sees only `Int`s or only
# `HalfOddInteger`s.  A mixture of the two kinds, a narrow or unsigned integer, and a range
# whose step is not 1 are refused with an explanation before anything is computed.  The
# keyword `ℓₘᵢₙ`, and its ASCII alias `ell_min`, is normalized against the kind of the
# positional indices.
#
# `ℓₘᵢₙ` defaults to the smallest ``|s|`` among the spin weights, which is `abs(s)` when
# there is only one and 0 (or 1/2) for a range that straddles zero.  That default is
# evaluated in the work method, after the positional indices have been converted, so that it
# is of their kind.

"""
    sYlm(R, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s))
    sYlm(θ, ϕ, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s))
    sYlm(R⃗, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s))

The spin-weighted spherical harmonics ``{}_sY_{ℓ,m}(R)`` for all ``ℓₘᵢₙ ≤ ℓ ≤ ℓₘₐₓ`` and
``-ℓ ≤ m ≤ ℓ``, as a [`HarmonicValues`](@ref): indexed first by ``ℓ`` and then naturally, so
that `sYlm(R, ℓₘₐₓ, s)[ℓ][m]` is one value.  With `ℓₘᵢₙ` below `abs(s)` the entries for ``ℓ
< |s|`` are zero.  The keyword may also be given as `ell_min`.  The computation is done in
the rotor's own floating-point type, and the result is `Complex` of it; to compute in
another type, convert the rotor.

The first argument may be a single rotor or a vector of them, and `s` may be a single spin
weight or an ascending unit range such as `-2:2`, which between them give a block four
shapes:

| built for | `sY[ℓ]` is indexed |
|---|---|
| one rotor, one spin weight | `[m]` |
| many rotors, one spin weight | `[iᵣ, m]` |
| one rotor, a range of spin weights | `[s, m]` |
| many rotors, a range of spin weights | `[iᵣ, s, m]` |

A vector of rotors gives the batched shapes whatever its length, even one.  For a range,
`ℓₘᵢₙ` defaults to the smallest ``|s|`` in it, which is 0 (or 1/2) whenever the range
straddles zero, so rows below their own ``|s|`` are zero.

Underneath, the values are held in one array whose *last* axis is the modes in the canonical
ordering `[ₛYₗₘ for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ]` (see [`Yindex`](@ref)), and whose leading
axes are the rotors and spin weights.  [`array_view`](@ref) hands that array back, which is
the form a product with mode weights takes; [`sYlm_matrix`](@ref) is the direct name for it,
for callers who want the bare array.

A single rotor may instead be given by the spherical coordinates of its point: `sYlm(θ, ϕ,
ℓₘₐₓ, s)` is `sYlm(from_spherical_coordinates(θ, ϕ), ℓₘₐₓ, s)`.  For repeated evaluation use
an [`sYlmCalculator`](@ref) (with [`sYlm!`](@ref), or directly), which allocates once.  The
convention is ``{}_sY_{ℓ,m} = (-1)^s \\sqrt{(2ℓ+1)/4π}\\, \\overline{𝔇^{(ℓ)}_{m,-s}}``; see
the "Conventions" section of the documentation.

# Half-integer indices

`ℓₘₐₓ` and `s` may be half-integers, passed as `Rational`s with denominator 2 or as
[`HalfOddInteger`](@ref)s — as in `sYlm(R, 7//2, 1//2)`, or `sYlm(R, 7//2, -3//2:3//2)` — in
which case every ``ℓ`` and ``m`` is a half-odd-integer, and so is `ℓₘᵢₙ`, which may be as
small as `1//2`.  The indices in one call must all be of one kind, integers or
half-odd-integers; a call that mixes them, such as `sYlm(R, 7//2, 1)`, is refused with an
`ArgumentError`.  For half-integer `s` the prefactor ``(-1)^s`` is ``i^{2s} = \\pm i``, so
the values include that phase and are not real multiples of ``\\overline{𝔇}``; the choice
of branch is explained under [`sYlmCalculator`](@ref).
"""
@index_methods function sYlm(
    R::RotorLike, ℓₘₐₓ::IndexType, s::IndexOrRange;
    ell_min::IndexType=min_abs_spin(s), ℓₘᵢₙ::IndexType=ell_min
)
    HarmonicValues(sYlm_array(R, ℓₘₐₓ, s, ℓₘᵢₙ), s, ℓₘᵢₙ, ℓₘₐₓ, 1)
end
@index_methods function sYlm(
    θ::Real, ϕ::Real, ℓₘₐₓ::IndexType, s::IndexOrRange;
    ell_min::IndexType=min_abs_spin(s), ℓₘᵢₙ::IndexType=ell_min
)
    sYlm(from_spherical_coordinates(θ, ϕ), ℓₘₐₓ, s; ℓₘᵢₙ)
end
function sYlm_helper(make, R, ℓₘₐₓ::IT, s, ℓₘᵢₙ::IT) where {IT<:IntegerHalf}
    HarmonicValues(harmonic_array(make, R, ℓₘₐₓ, s, ℓₘᵢₙ), s, ℓₘᵢₙ, ℓₘₐₓ, 1)
end
function harmonic_array(make, R, ℓₘₐₓ::IT, s, ℓₘᵢₙ::IT) where {IT<:IntegerHalf}
    check_sYlm_args(ℓₘₐₓ, s, ℓₘᵢₙ)
    # The calculator decides the element type, and the output buffer follows it, so that
    # there is exactly one place where that decision is made.  `make` is what chooses the
    # flavor: `sYlmCalculator_helper` for the complex harmonics, `sλlmCalculator_helper` for
    # the real ones.  Everything below is common to both.
    calc = make(R, ℓₘₐₓ, s)
    Y = allocate_sYlm(number_type(calc), s, ℓₘᵢₙ, ℓₘₐₓ)
    fill_sYlm!(Y, calc, s, ℓₘᵢₙ)
    Y
end

# The values of `sYlm` for a single rotor, as the bare array that `HarmonicValues` labels: a
# vector of modes for one spin weight, or a matrix of spin weights by modes for a range of
# them.  Like `D_array`, this is the function to which the rules for automatic
# differentiation are attached (see `src/derivatives/kernels.jl`), because it takes the
# rotor and returns a plain array.
sYlm_array(R::RotorLike, ℓₘₐₓ::IT, s, ℓₘᵢₙ::IT) where {IT<:IntegerHalf} =
    harmonic_array(sYlmCalculator_helper, R, ℓₘₐₓ, s, ℓₘᵢₙ)

# Many rotors at once.  The storage and the recursion are `sYlm_matrix`'s — that is the
# efficient path, and there is no reason to have two — so this labels the same array rather
# than computing it again.  `sYlm_matrix` remains the way to ask for the bare array.
@index_methods function sYlm(
    R⃗::AbstractVector{<:RotorLike}, ℓₘₐₓ::IndexType, s::IndexOrRange;
    ell_min::IndexType=min_abs_spin(s), ℓₘᵢₙ::IndexType=ell_min
)
    HarmonicValues(sYlm_matrix_array(R⃗, ℓₘₐₓ, s, ℓₘᵢₙ), s, ℓₘᵢₙ, ℓₘₐₓ, length(R⃗))
end
function sYlm_batch_helper(
    make, R⃗::AbstractVector, ℓₘₐₓ::IT, s, ℓₘᵢₙ::IT
) where {IT<:IntegerHalf}
    HarmonicValues(sYlm_matrix_helper(make, R⃗, ℓₘₐₓ, s, ℓₘᵢₙ), s, ℓₘᵢₙ, ℓₘₐₓ, length(R⃗))
end

# The output of a flat call: a vector of modes for one spin weight, and a matrix of spin
# weights by modes for several.
function allocate_sYlm(::Type{T}, ::IntegerHalf, ℓₘᵢₙ, ℓₘₐₓ) where {T}
    Vector{T}(undef, Ysize(ℓₘᵢₙ, ℓₘₐₓ))
end
function allocate_sYlm(::Type{T}, s::AbstractUnitRange, ℓₘᵢₙ, ℓₘₐₓ) where {T}
    Matrix{T}(undef, length(s), Ysize(ℓₘᵢₙ, ℓₘₐₓ))
end

# These checks hold for either kind of index, and for one spin weight or many.  The floor of
# ℓₘₐₓ and ℓₘᵢₙ is 0 for integers and 1/2 for half-odd-integers, and `< 0` is the right test
# for both, since no half-odd-integer lies between 0 and 1/2; the message names the floor of
# the kind at hand, which is what the `ℓₘᵢₙ` accessor gives for the index type.  (The
# accessor is called qualified, because the argument of the same name shadows it here.)
function check_sYlm_args(ℓₘₐₓ, s, ℓₘᵢₙ)
    check_spin_range(s)
    lowest = SphericalFunctions.ℓₘᵢₙ(typeof(ℓₘₐₓ))
    if ℓₘₐₓ < 0
        throw(ArgumentError("ℓₘₐₓ=$ℓₘₐₓ must be at least $lowest."))
    end
    sₕ = max_abs_spin(s)
    if sₕ > ℓₘₐₓ
        throw(ArgumentError("|s|=$sₕ exceeds ℓₘₐₓ=$ℓₘₐₓ; there are no such harmonics."))
    end
    if ℓₘᵢₙ < 0 || ℓₘᵢₙ > ℓₘₐₓ
        throw(ArgumentError("ℓₘᵢₙ=$ℓₘᵢₙ must satisfy $lowest ≤ ℓₘᵢₙ ≤ ℓₘₐₓ=$ℓₘₐₓ."))
    end
end

"""
    Ylm(R, ℓₘₐₓ; ℓₘᵢₙ=0)
    Ylm(θ, ϕ, ℓₘₐₓ; ℓₘᵢₙ=0)
    Ylm(R⃗, ℓₘₐₓ; ℓₘᵢₙ=0)

The ordinary scalar spherical harmonics ``Y_{ℓ,m}(R)`` for all ``ℓₘᵢₙ ≤ ℓ ≤ ℓₘₐₓ``, as a
[`HarmonicValues`](@ref) indexed by ``ℓ`` and then by ``m``: `Ylm(R, ℓₘₐₓ)[ℓ][m]`.  The
keyword may also be given as `ell_min`.

These are the spin-weight-zero case of the spin-weighted harmonics; this function is exactly
`sYlm(R, ℓₘₐₓ, 0; ℓₘᵢₙ)`, and everything [`sYlm`](@ref) says applies here too — including
that a vector of rotors gives blocks indexed `[iᵣ, m]`, and that `Ylm(θ, ϕ, ℓₘₐₓ)` is
`Ylm(from_spherical_coordinates(θ, ϕ), ℓₘₐₓ)`.  See [`YlmCalculator`](@ref) for repeated
evaluation.

The ``ℓ`` arguments must be integers.  Half-integer ``ℓ`` goes with half-integer spin
weight, so there is no half-integer analogue of this function, and a half-integer index is
refused with an `ArgumentError` that says so; see [`sYlm`](@ref) for half-integer spin
weights.
"""
@index_methods integer_only (
    "Half-integer ℓ goes with half-integer spin weight, so spin weight zero has none; "
    * "`sYlm` and `sYlmCalculator` accept half-integer indices."
) function Ylm(R::RotorLike, ℓₘₐₓ::IndexType; ell_min::IndexType=0, ℓₘᵢₙ::IndexType=ell_min)
    sYlm(R, ℓₘₐₓ, 0; ℓₘᵢₙ)
end
@index_methods integer_only (
    "Half-integer ℓ goes with half-integer spin weight, so spin weight zero has none; "
    * "`sYlm` and `sYlmCalculator` accept half-integer indices."
) function Ylm(
    θ::Real, ϕ::Real, ℓₘₐₓ::IndexType; ell_min::IndexType=0, ℓₘᵢₙ::IndexType=ell_min
)
    sYlm(θ, ϕ, ℓₘₐₓ, 0; ℓₘᵢₙ)
end
@index_methods integer_only (
    "Half-integer ℓ goes with half-integer spin weight, so spin weight zero has none; "
    * "`sYlm` and `sYlmCalculator` accept half-integer indices."
) function Ylm(
    R⃗::AbstractVector{<:RotorLike}, ℓₘₐₓ::IndexType; ell_min::IndexType=0, ℓₘᵢₙ::IndexType=ell_min
)
    sYlm(R⃗, ℓₘₐₓ, 0; ℓₘᵢₙ)
end
sYlm(R::NonRotorData, ℓₘₐₓ, s; kwargs...) = throw(ArgumentError(not_a_rotor(R)))
Ylm(R::NonRotorData, ℓₘₐₓ; kwargs...) = throw(ArgumentError(not_a_rotor(R)))
sYlm_matrix(R⃗::NonRotorData, ℓₘₐₓ, s; kwargs...) = throw(ArgumentError(not_a_rotor(R⃗)))

"""
    YlmCalculator(R, ℓₘₐₓ)
    YlmCalculator(θ, ϕ, ℓₘₐₓ)

Calculator for the ordinary scalar spherical harmonics ``Y_{ℓ,m}``, for ``ℓ ≤ ℓₘₐₓ``.

This is exactly `sYlmCalculator(R, ℓₘₐₓ, 0)`, or `sYlmCalculator(θ, ϕ, ℓₘₐₓ, 0)`, and
returns an [`sYlmCalculator`](@ref) rather than a type of its own — there is nothing about
spin weight zero to specialize.  Everything `sYlmCalculator` documents applies, so blocks
are indexed `Yₗ[m]` (or `Yₗ[iᵣ, m]` for a collection of rotors), and it iterates as `for (ℓ,
Yₗ) ∈ calc`.

``ℓₘₐₓ`` must be an integer: a half-integer ``ℓ`` goes with a half-integer spin weight, so
spin weight zero has no half-integer analogue, and a half-integer `ℓₘₐₓ` is refused with an
`ArgumentError` that says so.  See [`sYlmCalculator`](@ref) for those.
"""
@index_methods integer_only (
    "Half-integer ℓ goes with half-integer spin weight, so spin weight zero has none; "
    * "`sYlm` and `sYlmCalculator` accept half-integer indices."
) function YlmCalculator(R, ℓₘₐₓ::IndexType)
    sYlmCalculator(R, ℓₘₐₓ, 0)
end
@index_methods integer_only (
    "Half-integer ℓ goes with half-integer spin weight, so spin weight zero has none; "
    * "`sYlm` and `sYlmCalculator` accept half-integer indices."
) function YlmCalculator(θ::Real, ϕ::Real, ℓₘₐₓ::IndexType)
    sYlmCalculator(θ, ϕ, ℓₘₐₓ, 0)
end

"""
    sYlm!(Y, R, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s))
    sYlm!(Y, calc::sYlmCalculator, R; ℓₘᵢₙ=abs(s))
    sYlm!(Y, calc::sYlmCalculator, R, s; ℓₘᵢₙ=abs(s))

In-place version of [`sYlm`](@ref): fills `Y` with ``{}_sY_{ℓ,m}(R)`` in the canonical
ordering.  For a single spin weight `Y` is a vector, and its first `Ysize(ℓₘᵢₙ, ℓₘₐₓ)`
elements are written; for a range of them `Y` is a matrix, and the first `length(s)` rows
and `Ysize(ℓₘᵢₙ, ℓₘₐₓ)` columns are.  An output of the other shape is refused with an
`ArgumentError`.  The returned object is a [`HarmonicValues`](@ref) over the part of `Y`
that was written — a view, sharing `Y`'s storage — so that the values can be indexed as
`Y[ℓ][m]`, as those of `sYlm` are, whatever the size of `Y`.  The keyword `ℓₘᵢₙ` may also be
given as `ell_min`.  Because a longer `Y` is accepted, the layout is fixed by `ℓₘᵢₙ`, which
defaults to ``|s|``; a buffer laid out from ``ℓ = 0``, such as the one version 2's
`sYlm_prep` allocated, needs `ℓₘᵢₙ=0` to be filled in the layout it was sized for.

`Y` may also be a [`HarmonicValues`](@ref) built for one rotor, such as an earlier result of
`sYlm`, which is then refilled in place.  Its labels do not change, so the arguments must
agree with them — the same ``ℓₘₐₓ`` and spin weights, and the same `ℓₘᵢₙ`, which defaults to
the container's own rather than to ``|s|`` — and anything else is an `ArgumentError`.  So is
a `HarmonicValues` computed for a vector of rotors, even a vector of one, whose storage has
a rotor axis.

The element type of `Y` must be `Complex` of the rotor's own floating-point type — or, in
the calculator forms, of the calculator's.  A mismatch is an error rather than a silent
conversion: the type of the rotor is what decides the arithmetic, and `Y` is where the
answer lands, so disagreement between them means one of the two is not what the caller
thinks it is.

The calculator forms reuse `calc`, which must have been built for a single rotor rather than
a vector of them, and whose `ℓₘₐₓ` is used, so that repeated evaluation at many rotors
allocates nothing; the calculator's rotor data are replaced by `R`.  Without a spin weight,
everything the calculator serves is written, which for a calculator built for several is the
matrix form.  With one, that spin weight is picked out of [`spins`](@ref)`(calc)`, only its
values are assembled, and the result is a vector.  Calling this once per spin weight is
another way to evaluate several of them at one point, but each call runs the recurrence
again, so the form without a spin weight, which runs it once for all of them, is the cheaper
one.  To evaluate at many rotors at once, build a calculator from the vector of rotors and
iterate it, or use [`sYlm_matrix`](@ref).

The indices may be half-integers, passed as `Rational`s with denominator 2 or as
[`HalfOddInteger`](@ref)s, on the terms described under [`sYlm`](@ref): all of one kind,
with `ℓₘᵢₙ` defaulting to the smallest ``|s|`` and the values including the phase
``i^{2s}``.  In the calculator forms the kind is already fixed by the calculator, whose
``ℓ`` are integers or half-odd-integers according to how it was constructed, and a spin
weight or `ℓₘᵢₙ` of the other kind is refused with a message saying so.  In every form an
integer index must be an `Int` (as opposed to `Int16`, for example).
"""
@index_methods function sYlm!(
    Y::HarmonicValues, R::RotorLike, ℓₘₐₓ::IndexType, s::IndexOrRange;
    ell_min::Union{Nothing, IndexType}=nothing, ℓₘᵢₙ::Union{Nothing, IndexType}=ell_min
)
    # A `HarmonicValues` is filled by writing through its storage, and its labels do not
    # change, so the arguments must describe exactly what those labels say it holds: one
    # rotor, the same ℓₘₐₓ and spin weights, and the same ℓₘᵢₙ.  An ℓₘᵢₙ that is not given
    # is therefore taken from the container, rather than defaulting to the smallest |s| and
    # writing a different layout under the old labels.  The container itself comes back.
    check_harmonic_labels(Y, ℓₘₐₓ, s, something(ℓₘᵢₙ, Y.ℓₘᵢₙ), ℓₘᵢₙ !== nothing)
    sYlm!(array_view(Y), R, ℓₘₐₓ, s; ℓₘᵢₙ=Y.ℓₘᵢₙ)
    Y
end
function sYlm!(
    Y::HarmonicValues, calc::sYlmCalculator, R::RotorLike;
    ell_min=nothing, ℓₘᵢₙ=calculator_ℓₘᵢₙ(calc, ell_min, "ell_min")
)
    check_harmonic_labels(Y, calc, calc.s, ℓₘᵢₙ)
    sYlm!(array_view(Y), calc, R; ℓₘᵢₙ=Y.ℓₘᵢₙ)
    Y
end
@index_methods function sYlm!(
    Y::HarmonicValues, calc::sYlmCalculator, R::RotorLike, s::IndexType;
    ell_min::Union{Nothing, IndexType}=nothing, ℓₘᵢₙ::Union{Nothing, IndexType}=ell_min
)
    check_harmonic_labels(Y, calc, s, ℓₘᵢₙ)
    sYlm!(array_view(Y), calc, R, s; ℓₘᵢₙ=Y.ℓₘᵢₙ)
    Y
end

# The calculator forms: the calculator fixes ℓₘₐₓ, and the spin weight is the calculator's
# or the one picked out of it.
function check_harmonic_labels(Y::HarmonicValues, calc::HarmonicCalculator, s, ℓₘᵢₙ)
    check_harmonic_labels(
        Y, SphericalFunctions.ℓₘₐₓ(calc), s,
        ℓₘᵢₙ === nothing ? Y.ℓₘᵢₙ : calculator_ℓₘᵢₙ(calc, ℓₘᵢₙ, "ℓₘᵢₙ"),
        ℓₘᵢₙ !== nothing
    )
end

# The `ℓₘᵢₙ` of a calculator form, whose kind is fixed by the calculator rather than by
# dispatch: `nothing` where it is not given, and otherwise checked and converted as the
# index of a container is, so that the calculator forms accept exactly what the forms
# defined with `@index_methods` accept, an `Int` for an integer calculator and a
# `HalfOddInteger` or a `Rational{Int}` with denominator 2 for a half-integer one.  The
# forms without a spin weight check each spelling of the keyword under its own name,
# `ell_min` as the default of `ℓₘᵢₙ`.
calculator_ℓₘᵢₙ(calc, ::Nothing, name) = nothing
calculator_ℓₘᵢₙ(calc, ℓₘᵢₙ, name) = container_index(indextype(calc), ℓₘᵢₙ, calc, name)
# A calculator is named in those messages by the constructor that built it.
container_name(::HarmonicCalculator{IT, RT, NT}) where {IT, RT, NT} = flavor_name(NT)

function check_harmonic_labels(Y::HarmonicValues, ℓₘₐₓ, s, ℓₘᵢₙ, ℓₘᵢₙ_given::Bool)
    # A batch is a batch however short it is: values computed for a vector of rotors have a
    # rotor axis in their storage, even for a vector of one, and cannot be refilled from a
    # single rotor.
    if isbatched(Y)
        throw(ArgumentError(
            "Only a HarmonicValues built for one rotor can be filled in place; this one was "
            * "built for a vector of rotors, of length $(Nᵣ(Y)), and its storage has a rotor "
            * "axis."
        ))
    end
    if ℓₘₐₓ != Y.ℓₘₐₓ || s != Y.s || (ℓₘᵢₙ_given && ℓₘᵢₙ != Y.ℓₘᵢₙ)
        asked = ℓₘᵢₙ_given ? "ℓ ∈ $ℓₘᵢₙ:$ℓₘₐₓ" : "ℓₘₐₓ=$ℓₘₐₓ"
        throw(ArgumentError(
            "The HarmonicValues hold s=$(Y.s) and ℓ ∈ $(Y.ℓₘᵢₙ):$(Y.ℓₘₐₓ), and cannot be "
            * "filled with s=$s and $asked; the labels would no longer describe the values."
        ))
    end
    nothing
end

# The flat form without a calculator.  It is `@inline`, like the calculator forms below, for
# the reason given above `harmonic_values_view`.
@index_methods @inline function sYlm!(
    Y::AbstractVecOrMat{<:Complex}, R::RotorLike, ℓₘₐₓ::IndexType, s::IndexOrRange;
    ell_min::IndexType=min_abs_spin(s), ℓₘᵢₙ::IndexType=ell_min
)
    check_sYlm_args(ℓₘₐₓ, s, ℓₘᵢₙ)
    let calc = sYlmCalculator_helper(R, ℓₘₐₓ, s)
        check_sYlm_output(Y, calc, s, ℓₘᵢₙ)
        harmonic_values_view(fill_sYlm!(Y, calc, s, ℓₘᵢₙ)...)
    end
end
# In the calculator forms the calculator fixes the kind of index, so `ℓₘᵢₙ` is checked
# against the calculator's kind and converted to it here, by `calculator_ℓₘᵢₙ`, and the
# worker checks it, and the spin weight, against the calculator.
@inline function sYlm!(
    Y::AbstractVecOrMat{<:Complex}, calc::sYlmCalculator, R::RotorLike;
    ell_min=nothing, ℓₘᵢₙ=calculator_ℓₘᵢₙ(calc, ell_min, "ell_min")
)
    harmonic_values_view(sYlm_helper!(
        Y, calc, R, calc.s,
        ℓₘᵢₙ === nothing ? min_abs_spin(calc.s) : calculator_ℓₘᵢₙ(calc, ℓₘᵢₙ, "ℓₘᵢₙ")
    )...)
end
@index_methods @inline function sYlm!(
    Y::AbstractVecOrMat{<:Complex}, calc::sYlmCalculator, R::RotorLike, s::IndexType;
    ell_min::IndexType=abs(s), ℓₘᵢₙ::IndexType=ell_min
)
    harmonic_values_view(sYlm_helper!(Y, calc, R, s, ℓₘᵢₙ)...)
end

# The name of the in-place function a caller used, for the messages of the workers the two
# flavors share.
inplace_name(::HarmonicCalculator{IT, RT, NT}) where {IT, RT, NT} =
    NT <: Complex ? "sYlm!" : "sλlm!"

# The calculator forms of `sYlm!` and `sλlm!` reach the calculator with nothing checked but
# the kinds of their own indices, so everything about the call is checked here, before the
# calculator's rotor data are replaced: the kind of the spin weight and of `ℓₘᵢₙ` against
# the calculator's, the limits, the spin weights against those the calculator serves, the
# calculator itself, and the output.  The loop is `fill_sYlm!`, which the forms that build a
# calculator of their own call directly, having checked what that calculator does not
# settle.
function sYlm_helper!(
    Y::AbstractVector, calc::HarmonicCalculator{IT}, R, s::IntegerHalf,
    ℓₘᵢₙ::IntegerHalf
) where {IT<:IntegerHalf}
    check_index_kind(IT, s, "the spin weight s")
    check_index_kind(IT, ℓₘᵢₙ, "ℓₘᵢₙ")
    s = convert(IT, s)
    ℓₘᵢₙ = convert(IT, ℓₘᵢₙ)
    check_sYlm_args(SphericalFunctions.ℓₘₐₓ(calc), s, ℓₘᵢₙ)
    check_spin(calc, s)
    check_sYlm_calculator(calc, R)
    check_sYlm_output(Y, calc, s, ℓₘᵢₙ)
    set_rotors!(calc, R)
    fill_sYlm!(Y, calc, s, ℓₘᵢₙ)
end
function sYlm_helper!(
    Y::AbstractMatrix, calc::HarmonicCalculator{IT}, R, s::AbstractUnitRange,
    ℓₘᵢₙ::IntegerHalf
) where {IT<:IntegerHalf}
    check_index_kind(IT, ℓₘᵢₙ, "ℓₘᵢₙ")
    ℓₘᵢₙ = convert(IT, ℓₘᵢₙ)
    check_sYlm_args(SphericalFunctions.ℓₘₐₓ(calc), s, ℓₘᵢₙ)
    check_spin(calc, first(s))
    check_spin(calc, last(s))
    check_sYlm_calculator(calc, R)
    check_sYlm_output(Y, calc, s, ℓₘᵢₙ)
    set_rotors!(calc, R)
    fill_sYlm!(Y, calc, s, ℓₘᵢₙ)
end
sYlm_helper!(Y, calc::HarmonicCalculator, R, s, ℓₘᵢₙ) = sYlm_shape_error(Y, calc, s)

# The output of a flat call must suit the spin argument — a vector for one spin weight, a
# matrix for a range of them — have the calculator's element type, and be large enough for
# the modes ℓₘᵢₙ:ℓₘₐₓ, all of which the loop below relies on under `@inbounds`.
function check_sYlm_output(Y::AbstractVector, calc::HarmonicCalculator, ::IntegerHalf, ℓₘᵢₙ)
    check_sYlm_eltype(Y, calc)
    Base.require_one_based_indexing(Y)  # Y is written at 1-based positions under `@inbounds`
    let needed = Ysize(ℓₘᵢₙ, SphericalFunctions.ℓₘₐₓ(calc))
        if length(Y) < needed
            throw(DimensionMismatch(
                "Output vector has length $(length(Y)); at least $needed is needed."
            ))
        end
    end
    nothing
end
function check_sYlm_output(
    Y::AbstractMatrix, calc::HarmonicCalculator, s::AbstractUnitRange, ℓₘᵢₙ
)
    check_sYlm_eltype(Y, calc)
    Base.require_one_based_indexing(Y)  # Y is written at 1-based positions under `@inbounds`
    let needed = Ysize(ℓₘᵢₙ, SphericalFunctions.ℓₘₐₓ(calc)), n = length(s)
        if size(Y, 1) < n || size(Y, 2) < needed
            throw(DimensionMismatch(
                "Output matrix has size $(size(Y)); at least ($n, $needed) is needed, for "
                * "$n spin weights and $needed modes."
            ))
        end
    end
    nothing
end
check_sYlm_output(Y, calc::HarmonicCalculator, s, ℓₘᵢₙ) = sYlm_shape_error(Y, calc, s)

# An output whose shape does not suit the spin argument: a matrix for one spin weight, a
# vector for a range of them, or storage with a rotor axis.
function sYlm_shape_error(Y, calc::HarmonicCalculator, s)
    throw(ArgumentError(
        "`$(inplace_name(calc))` fills a vector for one spin weight, and a matrix with "
        * "length(s) rows for a range of them; it was given an output of size $(size(Y)) "
        * "for the spin weights $s."
    ))
end

# Fill `Y` from the calculator, whose rotor data are already in place, in the canonical
# ordering from `ℓₘᵢₙ`, after everything above has been checked.  Only the spin weights `s`
# are assembled at each ℓ, even when the calculator serves more of them, so that reading one
# spin weight out of a calculator built for a range costs the recurrence and that one row.
function fill_sYlm!(
    Y::AbstractVector, calc::HarmonicCalculator, s::IntegerHalf, ℓₘᵢₙ::IntegerHalf
)
    ℓₘₐₓ = SphericalFunctions.ℓₘₐₓ(calc)
    iₛ = spin_index(calc, s)
    Yˡ = calc.Yˡ
    @inbounds for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ
        compute_block!(calc, ℓ, iₛ:iₛ, calc.Yˡ, 0)
        i₀ = Yindex(ℓ, -ℓ, ℓₘᵢₙ) - 1
        for j ∈ 1:2ℓ+1
            Y[i₀ + j] = Yˡ[1, iₛ, j]
        end
    end
    (Y, s, ℓₘᵢₙ, ℓₘₐₓ)  # the storage and its normalized labels; see `harmonic_values_view`
end
function fill_sYlm!(
    Y::AbstractMatrix, calc::HarmonicCalculator, s::AbstractUnitRange, ℓₘᵢₙ::IntegerHalf
)
    ℓₘₐₓ = SphericalFunctions.ℓₘₐₓ(calc)
    n = length(s)
    i₁ = spin_index(calc, first(s))
    Yˡ = calc.Yˡ
    @inbounds for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ
        compute_block!(calc, ℓ, i₁:(i₁ + n - 1), calc.Yˡ, 0)
        i₀ = Yindex(ℓ, -ℓ, ℓₘᵢₙ) - 1
        for j ∈ 1:2ℓ+1
            for i ∈ 1:n
                Y[i, i₀ + j] = Yˡ[1, i₁ + (i - 1), j]
            end
        end
    end
    (Y, s, ℓₘᵢₙ, ℓₘₐₓ)
end

# The in-place forms fill storage the caller supplies, which may be larger than needed, and
# return the values as a `HarmonicValues` over the part that was written.  It is always a
# view, even when the storage is exactly the right size, so that the storage is shared
# rather than copied and the type returned does not depend on the size.  (`sYlm` itself
# allocates exactly what it needs, and labels that array directly.)  The workers return the
# storage and its labels, and the wrapper is built in the small public methods, which are
# inlined: a caller that ignores the result — as a loop reusing a calculator does — then
# never allocates it, which it would if the wrapper were returned from the worker.
harmonic_values_view(Y::AbstractVector, s, ℓₘᵢₙ, ℓₘₐₓ) =
    HarmonicValues(view(Y, 1:Ysize(ℓₘᵢₙ, ℓₘₐₓ)), s, ℓₘᵢₙ, ℓₘₐₓ, 1)
harmonic_values_view(Y::AbstractMatrix, s, ℓₘᵢₙ, ℓₘₐₓ) =
    HarmonicValues(view(Y, 1:length(s), 1:Ysize(ℓₘᵢₙ, ℓₘₐₓ)), s, ℓₘᵢₙ, ℓₘₐₓ, 1)

# `R` is untyped because the same check serves a rotor and a bare angle θ;
# `check_rotor_type` is what decides whether the argument suits the calculator at all.  A
# calculator built from a vector is batched whatever its length, and its blocks have a rotor
# index that the output of these functions does not.
function check_sYlm_calculator(calc::HarmonicCalculator, R)
    if isbatched(calc)
        throw(ArgumentError(
            "`$(inplace_name(calc))` needs a calculator built for a single rotor; this one "
            * "was built from a vector of $(Nᵣ(calc)), and its blocks have a rotor index."
        ))
    end
    check_rotor_type(calc, R)
    nothing
end

function check_sYlm_eltype(Y, calc::HarmonicCalculator{IT, RT, NT}) where {IT, RT, NT}
    if eltype(Y) !== NT
        wanted = NT <: Complex ? "Complex{$RT}" : "$RT"
        throw(ArgumentError(
            "This calculator works in $RT, so the output array's element type must be "
            * "$wanted, not $(eltype(Y))."
        ))
    end
    nothing
end

"""
    sYlm_matrix(R⃗, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s))

The dense matrix of spin-weighted spherical harmonics ``{}_sY_{ℓ,m}(R_i)``, with rows
indexed by the rotors in `R⃗` and columns by the modes ``(ℓ, m)`` in the canonical ordering
`[(ℓ, m) for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ for m ∈ -ℓ:ℓ]` (see [`Yindex`](@ref)).  Row `i` equals
`sYlm(R⃗[i], ℓₘₐₓ, s; ℓₘᵢₙ)`.  The keyword may also be given as `ell_min`.  The computation
is done in the rotors' own floating-point type; to compute in another type, convert them.

Multiplying this matrix by the array of mode weights synthesizes the corresponding
spin-weighted function at the rotors, as `Y * array_view(f̃)`, provided the weights are
stored from the same `ℓₘᵢₙ`; the labelled `sYlm(R⃗, ℓₘₐₓ, s) * f̃` does the same, and checks
the labels.  The (pseudo)inverse of the matrix performs the analysis.  It is computed with
all rotors batched, which is the efficient path; for ``ℓₘₐₓ ≳ 64`` the matrix becomes large
and the transforms in the "Transformations" section of the documentation should be
preferred.

The argument `s` may also be a range of spin weights (ascending, and with step size 1).  The
result is then three-dimensional, indexed `[rotor, spin, mode]`, and is a stack of the
matrices above rather than one of them.  The spin index is an ordinary 1-based position, and
`ℓₘᵢₙ` defaults to the smallest ``|s|`` in the range, which is 0 (or 1/2) for a range that
straddles zero.  `Y[:, i, :]` is the synthesis matrix of spin weight `s[i]` for weights
stored from that same `ℓₘᵢₙ`.  Weights that start at their own ``|s[i]|``, as those the
transforms return do, need that `ℓₘᵢₙ` passed explicitly; the labelled `sYlm(R⃗, ℓₘₐₓ, s) *
w` picks the right spin row and ``ℓ`` range itself.

The indices may be half-integers, passed as `Rational`s with denominator 2 or as
[`HalfOddInteger`](@ref)s, on the terms described under [`sYlm`](@ref): all of one kind,
with `ℓₘᵢₙ` defaulting to the smallest ``|s|`` and the values including the phase
``i^{2s}``.  The columns are then indexed by half-odd ``(ℓ, m)`` in the same canonical
ordering, and [`Yindex`](@ref) locates them as before.
"""
@index_methods function sYlm_matrix(
    R⃗::AbstractVector{<:RotorLike}, ℓₘₐₓ::IndexType, s::IndexOrRange;
    ell_min::IndexType=min_abs_spin(s), ℓₘᵢₙ::IndexType=ell_min
)
    sYlm_matrix_array(R⃗, ℓₘₐₓ, s, ℓₘᵢₙ)
end
# The array of `sYlm_matrix`, and of `sYlm` of a vector of rotors, for which the extensions
# for ChainRulesCore and ReverseDiff define rules, as for `sYlm_array`.
sYlm_matrix_array(R⃗::AbstractVector, ℓₘₐₓ::IT, s, ℓₘᵢₙ::IT) where {IT<:IntegerHalf} =
    sYlm_matrix_helper(sYlmCalculator_helper, R⃗, ℓₘₐₓ, s, ℓₘᵢₙ)
function sYlm_matrix_helper(
    make, R⃗::AbstractVector, ℓₘₐₓ::IT, s, ℓₘᵢₙ::IT
) where {IT<:IntegerHalf}
    check_sYlm_args(ℓₘₐₓ, s, ℓₘᵢₙ)
    calc = make(R⃗, ℓₘₐₓ, s)
    # The values of each ℓ are written by `materialize!` straight into their columns of the
    # result, rather than into the calculator's block and then copied, which for many rotors
    # would be a second pass over an array of many gigabytes.  The result is allocated with
    # the calculator's own layout, [rotor, spin, mode], so that the columns of one ℓ have
    # the layout of the block; for one spin weight the spin axis has length 1, and the
    # `Matrix` returned shares its storage.
    Y = Array{number_type(calc), 3}(undef, length(R⃗), nspins(s), Ysize(ℓₘᵢₙ, ℓₘₐₓ))
    fill_sYlm_matrix!(Y, calc, ℓₘᵢₙ, ℓₘₐₓ)
    s isa AbstractUnitRange ? Y : reshape(Y, size(Y, 1), size(Y, 3))
end

# `calc` is a fresh calculator built for exactly the spin weights of `Y`, and `Y` has one
# row per rotor, so the columns of each ℓ are a block of the calculator's own shape.
function fill_sYlm_matrix!(Y::Array{<:Any, 3}, calc::HarmonicCalculator, ℓₘᵢₙ, ℓₘₐₓ)
    for ℓ ∈ ℓₘᵢₙ:ℓₘₐₓ
        i₀ = Yindex(ℓ, -ℓ, ℓₘᵢₙ) - 1
        compute_block!(calc, ℓ, Base.OneTo(nspins(calc.s)), Y, i₀)
    end
    Y
end


### The real flavor: ₛλₗₘ(θ), defined in the `sλlmCalculator` docstring
#
# Every one of these is the corresponding ₛYₗₘ form with `sλlmCalculator_helper` in place of
# `sYlmCalculator_helper`, and an angle in place of a rotor.  There is no separate
# machinery: the helpers above take the factory as an argument precisely so that this file
# does not acquire a second copy of the same loops.

"""
    sλlm(θ, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s))
    slambdalm(θ, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s))

The real functions ``{}_sλ_{ℓ,m}(θ)`` for all ``ℓₘᵢₙ ≤ ℓ ≤ ℓₘₐₓ``, returned in the same
[`HarmonicValues`](@ref) container as [`sYlm`](@ref) but holding real numbers.  These are
``{}_sY_{ℓ,m}(θ, 0)`` for integer ``s`` and ``{}_sY_{ℓ,m}(θ, 0) / i^{2s}`` for half-odd
``s`` (see [`sλlmCalculator`](@ref)), and are indexed exactly as described under
[`sYlm`](@ref).  `θ` is one angle or an `AbstractVector` of them, and `s` one spin weight or
an ascending unit range.  The keyword may also be given as `ell_min`, and `slambdalm` is an
ASCII alias of the function.

This is what the ring-based transforms use internally.  The values are bit-for-bit those of
[`sYlm`](@ref) at the rotor of ``(θ, 0)``, computed from the angle: the real part for an
integer spin weight, and ``±`` the imaginary part for a half-odd one.  For an integer spin
weight the division by ``i^{2s} = (-1)^s`` is simply the factor that ``{}_sY_{ℓ,m}`` already
includes, so ``{}_sλ_{ℓ,m}`` is the harmonic at ``ϕ = 0`` as the literature writes it; for a
half-odd spin weight it removes the ``\\pm i`` that would otherwise make the value
imaginary.  See [`sλlmCalculator`](@ref) for the details, and for repeated evaluation.

The indices may be half-integers on exactly the terms described under [`sYlm`](@ref).
"""
@index_methods function sλlm(
    θ::Real, ℓₘₐₓ::IndexType, s::IndexOrRange;
    ell_min::IndexType=min_abs_spin(s), ℓₘᵢₙ::IndexType=ell_min
)
    sYlm_helper(sλlmCalculator_helper, θ, ℓₘₐₓ, s, ℓₘᵢₙ)
end
@index_methods function sλlm(
    θ::AbstractVector{<:Real}, ℓₘₐₓ::IndexType, s::IndexOrRange;
    ell_min::IndexType=min_abs_spin(s), ℓₘᵢₙ::IndexType=ell_min
)
    sYlm_batch_helper(sλlmCalculator_helper, θ, ℓₘₐₓ, s, ℓₘᵢₙ)
end

"""
    sλlm!(Y, θ, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s))
    sλlm!(Y, calc::sλlmCalculator, θ; ℓₘᵢₙ=abs(s))
    sλlm!(Y, calc::sλlmCalculator, θ, s; ℓₘᵢₙ=abs(s))
    sλlm!(Y::HarmonicValues, args...; kwargs...)
    slambdalm!(Y, args...; kwargs...)

Write ``{}_sλ_{ℓ,m}(θ)`` into the existing real array `Y`, which must have the layout
[`sλlm`](@ref) would return and an element type matching the calculator's.  This is
[`sYlm!`](@ref) for the real flavor, and behaves identically in every other respect; in
particular, a `HarmonicValues` must be refilled with arguments that agree with its labels,
and `ℓₘᵢₙ` may also be given as `ell_min`.  `slambdalm!` is an ASCII alias of the function.
"""
@index_methods function sλlm!(
    Y::HarmonicValues, θ::Real, ℓₘₐₓ::IndexType, s::IndexOrRange;
    ell_min::Union{Nothing, IndexType}=nothing, ℓₘᵢₙ::Union{Nothing, IndexType}=ell_min
)
    # As for `sYlm!` into a `HarmonicValues`: the arguments must agree with the labels, and
    # ℓₘᵢₙ defaults to the container's own.
    check_harmonic_labels(Y, ℓₘₐₓ, s, something(ℓₘᵢₙ, Y.ℓₘᵢₙ), ℓₘᵢₙ !== nothing)
    sλlm!(array_view(Y), θ, ℓₘₐₓ, s; ℓₘᵢₙ=Y.ℓₘᵢₙ)
    Y
end
function sλlm!(
    Y::HarmonicValues, calc::sλlmCalculator, θ::Real;
    ell_min=nothing, ℓₘᵢₙ=calculator_ℓₘᵢₙ(calc, ell_min, "ell_min")
)
    check_harmonic_labels(Y, calc, calc.s, ℓₘᵢₙ)
    sλlm!(array_view(Y), calc, θ; ℓₘᵢₙ=Y.ℓₘᵢₙ)
    Y
end
@index_methods function sλlm!(
    Y::HarmonicValues, calc::sλlmCalculator, θ::Real, s::IndexType;
    ell_min::Union{Nothing, IndexType}=nothing, ℓₘᵢₙ::Union{Nothing, IndexType}=ell_min
)
    check_harmonic_labels(Y, calc, s, ℓₘᵢₙ)
    sλlm!(array_view(Y), calc, θ, s; ℓₘᵢₙ=Y.ℓₘᵢₙ)
    Y
end
@index_methods @inline function sλlm!(
    Y::AbstractVecOrMat{<:Real}, θ::Real, ℓₘₐₓ::IndexType, s::IndexOrRange;
    ell_min::IndexType=min_abs_spin(s), ℓₘᵢₙ::IndexType=ell_min
)
    check_sYlm_args(ℓₘₐₓ, s, ℓₘᵢₙ)
    let calc = sλlmCalculator_helper(θ, ℓₘₐₓ, s)
        check_sYlm_output(Y, calc, s, ℓₘᵢₙ)
        harmonic_values_view(fill_sYlm!(Y, calc, s, ℓₘᵢₙ)...)
    end
end
@inline function sλlm!(
    Y::AbstractVecOrMat{<:Real}, calc::sλlmCalculator, θ::Real;
    ell_min=nothing, ℓₘᵢₙ=calculator_ℓₘᵢₙ(calc, ell_min, "ell_min")
)
    harmonic_values_view(sYlm_helper!(
        Y, calc, θ, calc.s,
        ℓₘᵢₙ === nothing ? min_abs_spin(calc.s) : calculator_ℓₘᵢₙ(calc, ℓₘᵢₙ, "ℓₘᵢₙ")
    )...)
end
@index_methods @inline function sλlm!(
    Y::AbstractVecOrMat{<:Real}, calc::sλlmCalculator, θ::Real, s::IndexType;
    ell_min::IndexType=abs(s), ℓₘᵢₙ::IndexType=ell_min
)
    harmonic_values_view(sYlm_helper!(Y, calc, θ, s, ℓₘᵢₙ)...)
end

"""
    sλlm_matrix(θ⃗, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s))
    slambdalm_matrix(θ⃗, ℓₘₐₓ, s; ℓₘᵢₙ=abs(s))

The dense real matrix ``{}_sλ_{ℓ,m}(θ_i)``, with rows indexed by the angles in `θ⃗` and
columns by the modes in the canonical ordering of [`Yindex`](@ref).  This is
[`sYlm_matrix`](@ref) for the real flavor; `s` may likewise be a range, in which case the
result is three-dimensional and indexed `[angle, spin, mode]`, and `ℓₘᵢₙ` may also be given
as `ell_min`.  `slambdalm_matrix` is an ASCII alias of the function.
"""
@index_methods function sλlm_matrix(
    θ⃗::AbstractVector{<:Real}, ℓₘₐₓ::IndexType, s::IndexOrRange;
    ell_min::IndexType=min_abs_spin(s), ℓₘᵢₙ::IndexType=ell_min
)
    sYlm_matrix_helper(sλlmCalculator_helper, θ⃗, ℓₘₐₓ, s, ℓₘᵢₙ)
end

# ASCII aliases, for use where the Greek letter is inconvenient.
const slambdalm = sλlm
const slambdalm! = sλlm!
const slambdalm_matrix = sλlm_matrix
const slambdalmCalculator = sλlmCalculator
