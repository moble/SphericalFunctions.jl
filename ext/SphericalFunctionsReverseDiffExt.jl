module SphericalFunctionsReverseDiffExt

# ReverseDiff's rules for `D_array`, `sYlm_array`, and `sYlm_matrix_array`, from the
# generators, as described in `src/derivatives.jl`.  ReverseDiff tracks only real numbers
# and arrays of them, and a rule is an instruction recorded on its tape by `@grad`, whose
# inputs and output must be of those kinds.  So a rotor is passed to the instruction as its
# four tracked components, and the instruction's output is the real array of the values'
# real and imaginary parts, interleaved.  The complex values are then assembled from the
# elements of that tracked array, each of which sends its derivative back to the array, and
# so to the instruction.  The harmonics of a vector of rotors are recorded as one
# instruction for each rotor, since an instruction with every rotor's components as its
# inputs would have to be compiled anew for every number of rotors.
#
# ReverseDiff does not read ChainRules' `rrule` for a function unless it is told to with
# `@grad_from_chainrules`, which requires the arguments to be tracked arrays or reals, as a
# rotor is not; hence the separate rules.
#
# A calculator of tracked rotors runs the recurrence in a calculator of their values, as a
# calculator of ForwardDiff's dual numbers does (see `src/utilities/lifting.jl`), and each
# of its steps records one instruction for each rotor, whose output is that rotor's block
# and whose pullback is the kernel of `src/derivatives.jl`.  The instruction is given a copy
# of the block of values, which the next step overwrites.  The elements of a block are
# assembled from the instruction's output, as above, and so are of the type of an element of
# a tracked `Vector`; a calculator stores only that type, and its rotors' components are
# passed through an instruction of their own, so that they are of that type too, whatever
# the type of the tracked numbers they were given as.
#
# ReverseDiff replays a recorded tape by running each instruction's function again on the
# new inputs, but it keeps the pullback of the first run.  A rule defined by `@grad` then
# gives the derivatives at the point where the tape was recorded, so none of these rules is
# for a tape that is recorded once and replayed; each gradient should record its own.

import SphericalFunctions
import SphericalFunctions: D_array, sYlm_array, sYlm_matrix_array, D_array_with_stored,
    D_array_pullback, harmonic_array_pullback!, derivatives_from_left, rotor_cotangent,
    rotor_cotangents, IntegerHalf, Ysize, nspins, WignerCalculator, HarmonicCalculator, Lift,
    wigner_block_pullback!, harmonic_block_pullback!, stored_m′range, stored_mrange,
    m′range, mrange, rotor_value
using Quaternionic: Quaternionic, Rotor, Quaternion
import ReverseDiff
using ReverseDiff: TrackedReal, TrackedArray, @grad, value, track

const TrackedRotorLike = Union{Rotor{<:TrackedReal}, Quaternion{<:TrackedReal}}

# The element type of a tracked `Vector`, which is what a calculator stores.
const TrackedElement{V, D} = TrackedReal{V, D, TrackedArray{V, D, 1, Vector{V}, Vector{D}}}

SphericalFunctions.value_type(::Type{<:TrackedReal{V}}) where {V} = V
SphericalFunctions.real_value(x::TrackedReal) = value(x)
SphericalFunctions.ndirections(::Type{<:TrackedReal}) = 0
SphericalFunctions.working_type(::Type{<:TrackedReal{V}}) where {V} = TrackedElement{V, V}

# The rotor of the values of four tracked components, built without normalizing it.
function value_rotor(w, x, y, z)
    W, X, Y, Z = value(w), value(x), value(y), value(z)
    Quaternion{typeof(W)}(W, X, Y, Z)
end

# The real and imaginary parts of the complex array `Y`, interleaved in a real vector, and
# the complex array of the given shape assembled from such a vector, from the element
# after `offset` on.
interleaved(Y::AbstractArray{<:Complex}) = vec(permutedims(hcat(vec(real.(Y)), vec(imag.(Y)))))
function assembled(y::AbstractVector, dims, offset=0)
    reshape([Complex(y[offset + 2k - 1], y[offset + 2k]) for k ∈ 1:prod(dims)], dims)
end
deinterleaved(Δ::AbstractVector, dims, offset=0) =
    reshape(complex.(Δ[(offset + 1):2:(offset + 2prod(dims))], Δ[(offset + 2):2:(offset + 2prod(dims))]), dims)


## 𝔇

function SphericalFunctions.D_array(
    R::TrackedRotorLike, ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT
) where {IT<:IntegerHalf}
    y = tracked_D_array(R[1], R[2], R[3], R[4], ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    o = 0
    map(SphericalFunctions.ℓₘᵢₙ(IT):ℓₘₐₓ) do ℓ
        dims = (length(max(-ℓ, m′ₘᵢₙ):min(ℓ, m′ₘₐₓ)), length(max(-ℓ, mₘᵢₙ):min(ℓ, mₘₐₓ)))
        b = assembled(y, dims, o)
        o += 2prod(dims)
        b
    end
end

tracked_D_array(w::TrackedReal, x::TrackedReal, y::TrackedReal, z::TrackedReal, args...) =
    track(tracked_D_array, w, x, y, z, args...)

@grad function tracked_D_array(w, x, y, z, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    R = value_rotor(w, x, y, z)
    blocks, stored, calc = D_array_with_stored(R, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    function tracked_D_array_pullback(Δ)
        o = 0
        Ā = map(blocks) do b
            Āᵢ = deinterleaved(Δ, size(b), o)
            o += 2length(b)
            Āᵢ
        end
        g = D_array_pullback(calc, stored, Ā)
        (rotor_cotangent(derivatives_from_left(calc), R, g)..., ntuple(_ -> nothing, 5)...)
    end
    (reduce(vcat, [interleaved(b) for b ∈ blocks]), tracked_D_array_pullback)
end


## The harmonics

function SphericalFunctions.sYlm_array(
    R::TrackedRotorLike, ℓₘₐₓ::IT, s, ℓₘᵢₙ::IT
) where {IT<:IntegerHalf}
    y = tracked_sYlm_array(R[1], R[2], R[3], R[4], ℓₘₐₓ, s, ℓₘᵢₙ)
    n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
    assembled(y, s isa AbstractUnitRange ? (length(s), n) : (n,))
end

tracked_sYlm_array(w::TrackedReal, x::TrackedReal, y::TrackedReal, z::TrackedReal, args...) =
    track(tracked_sYlm_array, w, x, y, z, args...)

@grad function tracked_sYlm_array(w, x, y, z, ℓₘₐₓ, s, ℓₘᵢₙ)
    R = value_rotor(w, x, y, z)
    Y = sYlm_array(R, ℓₘₐₓ, s, ℓₘᵢₙ)
    function tracked_sYlm_array_pullback(Δ)
        Ȳ = deinterleaved(Δ, size(Y))
        Ḡ = harmonic_array_pullback!(zeros(real(eltype(Y)), 3, 1), Y, Ȳ, false, ℓₘᵢₙ, ℓₘₐₓ)
        (only(rotor_cotangents(true, [R], Ḡ))..., nothing, nothing, nothing)
    end
    (interleaved(Y), tracked_sYlm_array_pullback)
end

# One instruction for each rotor, whose values are the rotor's row of the result.
function SphericalFunctions.sYlm_matrix_array(
    R⃗::AbstractVector{<:TrackedRotorLike}, ℓₘₐₓ::IT, s, ℓₘᵢₙ::IT
) where {IT<:IntegerHalf}
    n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
    rows = [sYlm_array(R, ℓₘₐₓ, s, ℓₘᵢₙ) for R ∈ R⃗]
    Y = [rows[i][j] for i ∈ eachindex(R⃗), j ∈ 1:(nspins(s) * n)]
    s isa AbstractUnitRange ? reshape(Y, length(R⃗), nspins(s), n) : Y
end


## Calculators

# The rotors, each passed through an instruction that returns its four components as a
# tracked `Vector`, whose elements are of the type the calculator stores.
function SphericalFunctions.store_rotors!(
    rotors::AbstractVector{Quaternion{T}}, R::AbstractVector
) where {T<:TrackedReal}
    for i ∈ eachindex(rotors, R)
        q = tracked_components(R[i][1], R[i][2], R[i][3], R[i][4])
        rotors[i] = Quaternion{T}(q[1], q[2], q[3], q[4])
    end
    rotors
end

tracked_components(w::TrackedReal, x::TrackedReal, y::TrackedReal, z::TrackedReal) =
    track(tracked_components, w, x, y, z)

@grad function tracked_components(w, x, y, z)
    tracked_components_pullback(Δ) = (Δ[1], Δ[2], Δ[3], Δ[4])
    ([value(w), value(x), value(y), value(z)], tracked_components_pullback)
end

# The steps record the derivatives themselves, so there are no generators to compute.
SphericalFunctions.set_generators!(
    lift::Lift, ::Bool, ::AbstractVector{Quaternion{RT}}
) where {RT<:TrackedReal} = lift

# The labelled rows and columns of the block of degree ℓ, from the calculator of values, one
# instruction for each rotor.
function SphericalFunctions.lift!(
    c::WignerCalculator{IT, RT, NT, ST, B, FT, <:Lift}, ℓ::IT
) where {IT, RT<:TrackedReal, NT, ST, B, FT<:Real}
    inner = c.lift.inner
    rows, cols = stored_m′range(inner, ℓ), stored_mrange(inner, ℓ)
    outrows, outcols = m′range(c, ℓ), mrange(c, ℓ)
    o′ = Int(first(outrows) - first(stored_m′range(c, ℓ)))
    o = Int(first(outcols) - first(stored_mrange(c, ℓ)))
    left = derivatives_from_left(c)
    for iᵣ ∈ eachindex(c.rotors)
        q = c.rotors[iᵣ]
        A = inner.Wˡ[iᵣ:iᵣ, 1:length(rows), 1:length(cols)]
        y = tracked_wigner_block(
            q[1], q[2], q[3], q[4], A, rotor_value(q), ℓ, rows, cols, outrows, outcols, left
        )
        k = 0
        for j ∈ eachindex(outcols), j′ ∈ eachindex(outrows)
            c.Wˡ[iᵣ, o′ + j′, o + j] = Complex(y[2k + 1], y[2k + 2])
            k += 1
        end
    end
    c
end

tracked_wigner_block(w::TrackedReal, x::TrackedReal, y::TrackedReal, z::TrackedReal, args...) =
    track(tracked_wigner_block, w, x, y, z, args...)

@grad function tracked_wigner_block(w, x, y, z, A, R, ℓ, rows, cols, outrows, outcols, left)
    o′, o = Int(first(outrows) - first(rows)), Int(first(outcols) - first(cols))
    values = A[1, (o′ + 1):(o′ + length(outrows)), (o + 1):(o + length(outcols))]
    function tracked_wigner_block_pullback(Δ)
        Ā = reshape(deinterleaved(Δ, size(values)), 1, size(values)...)
        Ḡ = wigner_block_pullback!(zeros(real(eltype(A)), 3, 1), A, Ā, ℓ, rows, cols, outrows, outcols, left)
        (rotor_cotangent(left, R, (Ḡ[1], Ḡ[2], Ḡ[3]))..., ntuple(_ -> nothing, 8)...)
    end
    (interleaved(values), tracked_wigner_block_pullback)
end

# The spin rows `is` of the block of degree ℓ, written into `Y` after its first `j₀` modes,
# from the calculator of values, one instruction for each rotor.
function SphericalFunctions.lift!(
    c::HarmonicCalculator{IT, RT, NT, ST, S, B, FT, <:Lift}, ℓ::IT, is, Y, j₀::Int
) where {IT, RT<:TrackedReal, NT, ST, S, B, FT<:Real}
    inner = c.lift.inner
    n = Int(2ℓ) + 1
    for iᵣ ∈ eachindex(c.rotors)
        q = c.rotors[iᵣ]
        A = inner.Yˡ[iᵣ:iᵣ, :, 1:n]
        y = tracked_harmonic_block(q[1], q[2], q[3], q[4], A, rotor_value(q), ℓ, is)
        k = 0
        for j ∈ 1:n, i ∈ is
            Y[iᵣ, i, j₀ + j] = Complex(y[2k + 1], y[2k + 2])
            k += 1
        end
    end
    c
end

tracked_harmonic_block(w::TrackedReal, x::TrackedReal, y::TrackedReal, z::TrackedReal, args...) =
    track(tracked_harmonic_block, w, x, y, z, args...)

@grad function tracked_harmonic_block(w, x, y, z, A, R, ℓ, is)
    values = A[1, is, :]
    function tracked_harmonic_block_pullback(Δ)
        Ā = zeros(eltype(A), size(A))
        Ā[1, is, :] .= deinterleaved(Δ, size(values))
        Ḡ = harmonic_block_pullback!(zeros(real(eltype(A)), 3, 1), A, Ā, ℓ, is)
        (rotor_cotangent(true, R, (Ḡ[1], Ḡ[2], Ḡ[3]))..., ntuple(_ -> nothing, 4)...)
    end
    (interleaved(values), tracked_harmonic_block_pullback)
end

end # module SphericalFunctionsReverseDiffExt
