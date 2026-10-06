### The state shared by the calculators of 𝔇 and of the harmonics
#
# A calculator of 𝔇 or of `d`, and a calculator of ₛYₗₘ or of ₛλₗₘ, each runs the ``H``
# recurrence and multiplies each element of the wedge by a sign and, for 𝔇 and ₛYₗₘ, by a
# phase from the powers of the rotor's phases z₊ and z₋.  What they share for that, in the
# real type in which the recurrence runs, is this engine: the `HCalculator` and the tables
# of the powers of z₊ and z₋.  The tables hold the powers for every calculator of 𝔇, and
# for a calculator of ₛYₗₘ exactly when it was given rotors rather than angles θ, which mean
# the points (θ, ϕ=0); a calculator of ₛYₗₘ records which in its own `phases` flag.  They
# are empty for every calculator of `d` or of ₛλₗₘ.  Each calculator keeps its own copy of
# its rotors, in its own real type, and its own block.  A calculator whose rotors hold
# derivatives lifts the blocks of a calculator of their values (see
# `src/derivatives/lifting.jl`), and holds that calculator's engine as its own.
#
# The engine is immutable, and stored inline in the calculator, so that building a
# calculator allocates nothing for it beyond its buffers.  The flag of ₛYₗₘ is a field of
# that calculator rather than of the engine, because a calculator of 𝔇 always holds the
# phases of its rotors and one of `d` never does, so neither of them needs a flag.
# - `IT` is the index type, `Int` or `HalfOddInteger`.
# - `FT` is the real type in which the recurrence runs.
# - `ST` is the storage type of the ``H`` wedge.
struct SphericalFunctionsEngine{IT, FT<:Real, ST}
    H::HCalculator{IT, FT, ST}
    # K is `power_extent(ℓₘₐₓ, m′ₘₐₓ)`, where m′ₘₐₓ is the width of the wedge in `H`.
    Z₊::Matrix{Complex{FT}}  # Z₊[iᵣ, k+K+1] = z₊^k for k ∈ -K:K, or empty
    Z₋::Matrix{Complex{FT}}  # Z₋[iᵣ, k+K+1] = z₋^k for k ∈ -K:K, or empty
end

# The largest |k| of a power zᵏ that a block reads, for a calculator whose wedge holds the
# rows |m′| ≤ W: every element of a block has min(|m′|, |m|) ≤ W, so |m′ ± m| ≤ ℓₘₐₓ + W,
# and never more than 2ℓₘₐₓ.  It is never less than 2 unless 2ℓₘₐₓ is, because the
# recurrence of `complex_powers!` computes the powers of a table of length 2 differently
# from those of a longer one, which would change the signs of zeros at z = ±1 and ±i.
power_extent(ℓₘₐₓ, W) = min(Int(2ℓₘₐₓ), max(Int(ℓₘₐₓ + W), 2))

# The column of the power zᵏ in a power table `Z`, whose middle column holds z⁰.  A size is
# never negative, so the shift halves it exactly, without the correction for a negative
# sign that `÷ 2` makes.
@inline power_column(Z, k) = Int(k) + (size(Z, 2) + 1) >> 1

# A row of a power table, as the vector of the powers z⁰, z¹, …, zᴷ that `complex_powers!`
# fills: each power is stored twice, zᵏ at the column of k and conj(zᵏ) at that of -k, which
# is exact for a phase.  The conjugate is stored first, so that at k = 0, where the two
# columns coincide, the power itself, 1 + 0i, remains rather than its conjugate, 1 - 0i.
struct SignedPowers{T, V<:AbstractVector{T}} <: AbstractVector{T}
    row::V
end
Base.size(p::SignedPowers) = ((length(p.row) + 1) ÷ 2,)
Base.@propagate_inbounds Base.getindex(p::SignedPowers, j::Int) = p.row[length(p) - 1 + j]
Base.@propagate_inbounds function Base.setindex!(p::SignedPowers, w, j::Int)
    K = length(p) - 1
    p.row[K + 2 - j] = conj(w)
    p.row[K + j] = w
    p
end

# Allocate the buffers without touching them.  PRIVATE: see the note on `allocate_H`.  The
# power tables have a column for each power zᵏ with |k| ≤ power_extent(ℓₘₐₓ, m′ₘₐₓ), where
# m′ₘₐₓ is the width of the wedge, when `tables` is true, as it is for a calculator of 𝔇 or
# of ₛYₗₘ, and none for one of `d` or of ₛλₗₘ, whose tables are then empty rather than
# merely unread, so that nothing is allocated for them.  The tables have the rotor index
# first, so that the innermost loop of `materialize!`, which runs over the rotors at a fixed
# power, reads them contiguously.
#
# This is never inlined.  Inlined into its caller, it builds the engine in a slot on the
# caller's stack and copies the `HCalculator` into it, and Enzyme's compiler (as of Enzyme
# 0.13.205 on Julia 1.13) misreads the layout of that slot and refuses to compile the
# caller, which is every allocation of a calculator of 𝔇 or of the harmonics.
@noinline function allocate_engine(
    ::Type{IT}, ::Type{FT}, ℓₘₐₓ::IT, m′ₘₐₓ::IT, Nᵣ::Int, tables::Bool
) where {IT<:IntegerHalf, FT<:Real}
    H = allocate_H(IT, FT, ℓₘₐₓ, m′ₘₐₓ, Nᵣ)
    K = tables ? 2power_extent(ℓₘₐₓ, m′ₘₐₓ) + 1 : 0
    Z₊ = Matrix{Complex{FT}}(undef, Nᵣ, K)
    Z₋ = Matrix{Complex{FT}}(undef, Nᵣ, K)
    SphericalFunctionsEngine(H, Z₊, Z₋)
end

ℓₘₐₓ(e::SphericalFunctionsEngine) = ℓₘₐₓ(e.H)
Nᵣ(e::SphericalFunctionsEngine) = Nᵣ(e.H)


### Rotor data

# The rotor data of every rotor in `R`: the phase e^{iβ} and, on the half-integer path, the
# half angles, in the `HCalculator`, and the powers of z₊ and z₋ in the tables.  Every rotor
# is acceptable, as in `set_rotors!(::HCalculator, …)`, and the results are marked invalid
# before the first one is replaced.  A calculator's own `set_rotor_data!` is a call to this,
# and so is what the extensions for Enzyme and Mooncake declare to have no derivatives.
#
# The rotors are taken two at a time, so that the recurrences for the powers of their four
# phases, each a chain of dependent operations, overlap; the values are those of one rotor
# at a time.  (Four at a time is no faster: the states of eight recurrences no longer fit in
# the registers.)
function set_rotor_data!(e::SphericalFunctionsEngine, R::AbstractVector{<:RotorLike})
    Base.require_one_based_indexing(R)
    check_rotor_count(e, R)
    e.H.axes_valid[] = false
    iᵣ = 1
    while iᵣ < length(R)
        set_rotor_data!(e, R, (iᵣ, iᵣ + 1))
        iᵣ += 2
    end
    if iᵣ == length(R)
        set_rotor_data!(e, R, (iᵣ,))
    end
    nothing
end
# The rotor data of the rotors `R[iᵣ]` for each `iᵣ` in the tuple `rotors`.
@inline function set_rotor_data!(e::SphericalFunctionsEngine, R, rotors::Tuple)
    z = map(iᵣ -> store_rotor!(e.H, iᵣ, @inbounds R[iᵣ]), rotors)
    phase_powers!(
        (
            map(iᵣ -> SignedPowers(view(e.Z₊, iᵣ, :)), rotors)...,
            map(iᵣ -> SignedPowers(view(e.Z₋, iᵣ, :)), rotors)...
        ),
        (map(first, z)..., map(last, z)...)
    )
    nothing
end
# Angles θ, meaning the points (θ, ϕ=0), at which the phases are all 1.  The `HCalculator`
# validates everything before it replaces anything, so if it refuses the angles the engine
# is left exactly as it was.
function set_rotor_data!(e::SphericalFunctionsEngine, θ::Union{Real, AbstractVector{<:Real}})
    set_rotors!(e.H, θ)
    nothing
end

# The rotor data is copied buffer by buffer rather than derived again; see the comment on
# `copy_rotor_data!(::HCalculator, ::HCalculator)` for why it cannot be derived again from
# the `HCalculator`.
function copy_rotor_data!(e′::SphericalFunctionsEngine, e::SphericalFunctionsEngine)
    copy_rotor_data!(e′.H, e.H)
    copyto!(e′.Z₊, e.Z₊)
    copyto!(e′.Z₋, e.Z₋)
    e′
end

# The whole of the rotor state of a calculator `c` of 𝔇 or of the harmonics, copied into a
# calculator `c′` of the same type that `similar` has just allocated: its rotors or angles,
# and either its rotor data (see `copy_rotor_data!`) or, for a calculator that lifts the
# blocks of another, the generators of its rotors' derivatives and the rotor state of that
# calculator, whose engine it shares.  So the engine is copied once, at the bottom, and
# serves every level above it.
function copy_rotor_state!(c′::C, c::C) where {C<:AbstractCalculator}
    copyto!(c′.rotors, c.rotors)
    copyto!(c′.angles, c.angles)
    if c.lift === nothing
        copy_rotor_data!(c′, c)
    else
        copyto!(c′.lift.G, c.lift.G)
        copy_rotor_state!(c′.lift.inner, c.lift.inner)
    end
    c′
end

# Fill the buffers of the `HCalculator` with `v`, and mark its axes invalid, leaving the
# rotor data alone, as `fill!(::HCalculator, v)` does; the power tables are rotor data.
function Base.fill!(e::SphericalFunctionsEngine, v::Real)
    fill!(e.H, v)
    e
end


### The elements of a block

# The phase z₊^k₊ z₋^k₋ of an element, conjugated for 𝔇 and not for ₛYₗₘ.
@inline conjugated(z, ::Val{true}) = conj(z)
@inline conjugated(z, ::Val{false}) = z

# The number of rotors that a method compiled for `Nᵣ` writes: one, given as `Val(1)`, or
# the `Int` it is.
@inline rotor_count(::Val{1}) = 1
@inline rotor_count(Nᵣ::Int) = Nᵣ

# The offset, before the first rotor, of the power zᵏ in a table `Z` of the powers of `n`
# rotors' phases (see `power_column`), as `materialize_run!` takes it.
@inline power_offset(Z, k, n) = n * (power_column(Z, k) - 1)

# The coefficient of element p of a run, from those of its first two elements: the signs ϵ
# along a run are constant or alternate, so these are the coefficients of the even and of
# the odd elements.  Each is computed as the element's own coefficient would be, rather than
# the one as the negative of the other, which for a number that holds derivatives would give
# the zeros among them the other sign.
@inline run_coefficient((c₀, c₁), p) = ifelse(isodd(p), c₁, c₀)

# The values of an element for one rotor, and of its partner (see `materialize_run!`), from
# their coefficients `c` and `c′`, the wedge element `h` they share, and the entries `i₊`
# and `i₋` of the power tables for the element's own powers.  The partner's powers are the
# opposites of the element's, so its phase is the conjugate of the element's, and one
# product serves both.  Without phases (a block of `d`) the tables are not read at all, so
# that they may be empty, and the partner's value is its coefficient times `h`.  Without a
# partner (`c′` is `nothing`) the second value is `nothing`.
@inline element_values(c, ::Nothing, h, Z₊, i₊, Z₋, i₋, ::Val{false}, conjugate) =
    (c * h, nothing)
@inline element_values(c, c′, h, Z₊, i₊, Z₋, i₋, ::Val{false}, conjugate) = (c * h, c′ * h)
@inline element_values(c, ::Nothing, h, Z₊, i₊, Z₋, i₋, ::Val{true}, conjugate) =
    (c * h * conjugated((@inbounds Z₊[i₊]) * (@inbounds Z₋[i₋]), conjugate), nothing)
@inline function element_values(c, c′, h, Z₊, i₊, Z₋, i₋, ::Val{true}, conjugate)
    q = (@inbounds Z₊[i₊]) * (@inbounds Z₋[i₋])
    (c * h * conjugated(q, conjugate), c′ * h * conjugated(conj(q), conjugate))
end

# The partners of a run are `nothing`, or `(b, Δb, coefficients)`: the offset of the partner
# of the run's first element in the block, the step from one partner to the next, and their
# coefficients, as for the run itself.
@inline partner_coefficient(::Nothing, p) = nothing
@inline partner_coefficient((b, Δb, coefficients), p) = run_coefficient(coefficients, p)
@inline store_partner!(A, ::Nothing, p, iᵣ, y) = nothing
@inline store_partner!(A, (b, Δb, coefficients), p, iᵣ, y) =
    (@inbounds A[b + Δb * p + iᵣ] = y; nothing)

# A run of `n` elements of a block, for every rotor, and, when `partners` is not `nothing`,
# of their partners: every element that a calculator's `materialize!` computes from the
# wedge is written by this.  Each of `A`, `Hp`, `Z₊`, and `Z₋` is read by linear index, and
# is given as the 0-based offset of the first rotor of the run's first element and the step
# to the next element; rotor iᵣ of an element at offset `o` is at `o + iᵣ`, which is the
# layout of the blocks ([iᵣ, …]), of the wedge (see `wedge_offset`), and of the power tables
# ([iᵣ, k]; see `power_offset`).  Element p ∈ 0:n-1 and its partner are written as
# `element_values` gives them, with the coefficients of `run_coefficient`, so the loops have
# no branches, and every load and store advances by a fixed step.  `ivdep` asserts that the
# block does not overlap the wedge or the tables, which spares each run the checks for
# overlap that would otherwise cost as much as a short run itself.
#
# With one rotor, `Nᵣ` is `Val(1)`, and the loop runs over the elements.  The steps are then
# constants wherever the caller's are, as they are for every step through the wedge and the
# tables and for a step down a column of a block, so the compiler sees the unit strides and
# vectorizes the loop.  Otherwise `Nᵣ` is the `Int` it is, and the loop runs over the rotors
# of each element, which are contiguous in every array.  Both compute each element by the
# same expression, so the values do not depend on `Nᵣ`.
@inline function materialize_run!(
    A, Hp, Z₊, Z₋, ::Val{1}, n::Int, (a, Δa), partners, (h, Δh), (t₊, Δ₊), (t₋, Δ₋),
    coefficients, phases::Val, conjugate::Val
)
    @inbounds @simd ivdep for p ∈ 0:n-1
        x, y = element_values(
            run_coefficient(coefficients, p), partner_coefficient(partners, p),
            Hp[h + Δh * p + 1], Z₊, t₊ + Δ₊ * p + 1, Z₋, t₋ + Δ₋ * p + 1, phases, conjugate
        )
        A[a + Δa * p + 1] = x
        store_partner!(A, partners, p, 1, y)
    end
    nothing
end
@inline function materialize_run!(
    A, Hp, Z₊, Z₋, Nᵣ::Int, n::Int, (a, Δa), partners, (h, Δh), (t₊, Δ₊), (t₋, Δ₋),
    coefficients, phases::Val, conjugate::Val
)
    @inbounds for p ∈ 0:n-1
        c, c′ = run_coefficient(coefficients, p), partner_coefficient(partners, p)
        aₚ, hₚ, t₊ₚ, t₋ₚ = a + Δa * p, h + Δh * p, t₊ + Δ₊ * p, t₋ + Δ₋ * p
        @simd ivdep for iᵣ ∈ 1:Nᵣ
            x, y = element_values(
                c, c′, Hp[hₚ + iᵣ], Z₊, t₊ₚ + iᵣ, Z₋, t₋ₚ + iᵣ, phases, conjugate
            )
            A[aₚ + iᵣ] = x
            store_partner!(A, partners, p, iᵣ, y)
        end
    end
    nothing
end
