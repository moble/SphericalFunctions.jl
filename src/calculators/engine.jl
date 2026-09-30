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
    Z₊::Matrix{Complex{FT}}  # Z₊[iᵣ, k+1] = z₊^k for k ∈ 0:2ℓₘₐₓ, or empty
    Z₋::Matrix{Complex{FT}}  # Z₋[iᵣ, k+1] = z₋^k for k ∈ 0:2ℓₘₐₓ, or empty
end

# Allocate the buffers without touching them.  PRIVATE: see the note on `allocate_H`.  The
# power tables have a column for each power in 0:2ℓₘₐₓ when `tables` is true, as it is for a
# calculator of 𝔇 or of ₛYₗₘ, and none for one of `d` or of ₛλₗₘ, whose tables are then
# empty rather than merely unread, so that nothing is allocated for them.  The tables have
# the rotor index first, so that the innermost loop of `materialize!`, which runs over the
# rotors at a fixed power, reads them contiguously.
#
# This is never inlined.  Inlined into its caller, it builds the engine in a slot on the
# caller's stack and copies the `HCalculator` into it, and Enzyme's compiler (as of Enzyme
# 0.13.205 on Julia 1.13) misreads the layout of that slot and refuses to compile the
# caller, which is every allocation of a calculator of 𝔇 or of the harmonics.
@noinline function allocate_engine(
    ::Type{IT}, ::Type{FT}, ℓₘₐₓ::IT, m′ₘₐₓ::IT, Nᵣ::Int, tables::Bool
) where {IT<:IntegerHalf, FT<:Real}
    H = allocate_H(IT, FT, ℓₘₐₓ, m′ₘₐₓ, Nᵣ)
    K = tables ? 2ℓₘₐₓ + 1 : 0
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
function set_rotor_data!(e::SphericalFunctionsEngine, R::AbstractVector{<:RotorLike})
    Base.require_one_based_indexing(R)
    check_rotor_count(e, R)
    e.H.axes_valid[] = false
    @inbounds for i ∈ eachindex(R)
        z₊, z₋ = store_rotor!(e.H, i, R[i])
        complex_powers!(view(e.Z₊, i, :), z₊)
        complex_powers!(view(e.Z₋, i, :), z₋)
    end
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
# calculator `c′` of the same type that `similar` has just allocated: its rotors, and either
# its rotor data (see `copy_rotor_data!`) or, for a calculator that lifts the blocks of
# another, the generators of its rotors' derivatives and the rotor state of that calculator,
# whose engine it shares.  So the engine is copied once, at the bottom, and serves every
# level above it.
function copy_rotor_state!(c′::C, c::C) where {C<:AbstractCalculator}
    copyto!(c′.rotors, c.rotors)
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

# Powers zᵏ for k of either sign, given Z[iᵣ, k+1] = zᵏ for k ≥ 0.  The three-argument form
# decides the sign itself; the four-argument form is told it, as `Val(k < 0)`, so that a
# loop over rotors at a fixed k has no branch in it (see `with_power_signs`).
@inline zpower(Z, iᵣ, k, ::Val{false}) = @inbounds Z[iᵣ, k+1]
@inline zpower(Z, iᵣ, k, ::Val{true}) = conj(@inbounds Z[iᵣ, 1-k])
@inline zpower(Z, iᵣ, k) = k ≥ 0 ? zpower(Z, iᵣ, k, Val(false)) : zpower(Z, iᵣ, k, Val(true))

# Call `f(Val(k₊ < 0), Val(k₋ < 0))`, so that the signs of two powers are compile-time
# constants inside `f`.  Each loop over rotors in `materialize_element!` is at fixed powers,
# and settling there, once, which of the two table entries are conjugated leaves a loop body
# with no branch at all.  The arithmetic is the same in each of the four cases, so the
# values are those the three-argument `zpower` gives.
@inline function with_power_signs(f, k₊, k₋)
    if k₊ ≥ 0
        k₋ ≥ 0 ? f(Val(false), Val(false)) : f(Val(false), Val(true))
    else
        k₋ ≥ 0 ? f(Val(true), Val(false)) : f(Val(true), Val(true))
    end
end

# The phase z₊^k₊ z₋^k₋ of an element, conjugated for 𝔇 and not for ₛYₗₘ.
@inline conjugated(z, ::Val{true}) = conj(z)
@inline conjugated(z, ::Val{false}) = z

# One element of a block, at position (i, j) of the last two axes of `A`, for every rotor,
# from the wedge element whose first rotor is at `offset + 1` of `Hp`: the product of
# `coefficient`, the wedge element, and, when `A` is complex and `phases` is true, the phase
# z₊^k₊ z₋^k₋, conjugated as `conjugate` says.  For 𝔇, (i, j) is the position of (m′, m),
# and the phase is conj(z₊^(m′+m) z₋^(m′-m)) = e^{-i(m′α + mγ)}; for ₛYₗₘ, (i, j) is that
# of (s, m), and the phase is z₊^(m-s) z₋^(m+s) = e^{i(mα - sγ)}.
@inline function materialize_element!(
    A::AbstractArray{NT}, Hp, Z₊, Z₋, Nᵣ, i, j, offset, coefficient, k₊, k₋,
    phases::Bool, conjugate::Val
) where {NT}
    # With one rotor there is nothing to vectorize, and the setup of a loop would cost more
    # than the element itself, so that case is written out.
    if NT <: Complex && phases
        if Nᵣ == 1
            @inbounds A[1, i, j] = (
                coefficient * Hp[offset + 1]
                * conjugated(zpower(Z₊, 1, k₊) * zpower(Z₋, 1, k₋), conjugate)
            )
        elseif Nᵣ < 8
            # For a few rotors the four loops of `with_power_signs`, and the setup of a
            # vectorized loop, cost more than they save, so the sign of each power is
            # tested on every rotor instead.  Each element is computed by the same
            # expression in every branch, so the values do not depend on `Nᵣ`.
            @inbounds for iᵣ ∈ 1:Nᵣ
                A[iᵣ, i, j] = (
                    coefficient * Hp[offset + iᵣ]
                    * conjugated(zpower(Z₊, iᵣ, k₊) * zpower(Z₋, iᵣ, k₋), conjugate)
                )
            end
        else
            with_power_signs(k₊, k₋) do n₊, n₋
                @inbounds @simd for iᵣ ∈ 1:Nᵣ
                    phase = conjugated(
                        zpower(Z₊, iᵣ, k₊, n₊) * zpower(Z₋, iᵣ, k₋, n₋), conjugate
                    )
                    A[iᵣ, i, j] = coefficient * Hp[offset + iᵣ] * phase
                end
            end
        end
    elseif Nᵣ == 1
        @inbounds A[1, i, j] = coefficient * Hp[offset + 1]
    else
        @inbounds @simd for iᵣ ∈ 1:Nᵣ
            A[iᵣ, i, j] = coefficient * Hp[offset + iᵣ]
        end
    end
    nothing
end
