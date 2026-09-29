module SphericalFunctionsReverseDiffExt

# ReverseDiff's rules for `D_array` and `sYlm_array`, from the generators, as described in
# `src/derivatives.jl`.  ReverseDiff tracks only real numbers and arrays of them, and a rule
# is an instruction recorded on its tape by `@grad`, whose inputs and output must be of those
# kinds.  So the rotor is passed to the instruction as its four tracked components, and the
# instruction's output is the real array of the values' real and imaginary parts,
# interleaved.  The complex values are then assembled from the elements of that tracked
# array, each of which sends its derivative back to the array, and so to the instruction.
#
# ReverseDiff does not read ChainRules' `rrule` for a function unless it is told to with
# `@grad_from_chainrules`, which requires the arguments to be tracked arrays or reals, as a
# rotor is not; hence the separate rule.

import SphericalFunctions
import SphericalFunctions: D_array, sYlm_array, D_array_widened, D_is_widened, D_narrowed,
    D_pullback, sYlm_pullback, rotor_cotangent, D_offset, Ysize, IntegerHalf
using Quaternionic: Rotor
import ReverseDiff
using ReverseDiff: TrackedReal, @grad, value, track

# The rotor of the values of four tracked components, built without normalizing it.
function value_rotor(w, x, y, z)
    W, X, Y, Z = value(w), value(x), value(y), value(z)
    Rotor{typeof(W)}(W, X, Y, Z)
end

# The real and imaginary parts of the complex array `Y`, interleaved in a real vector, and
# the complex array of the given shape assembled from such a vector.
interleaved(Y::AbstractArray{<:Complex}) = vec(permutedims(hcat(vec(real.(Y)), vec(imag.(Y)))))
function assembled(y::AbstractVector, dims)
    reshape([Complex(y[2k - 1], y[2k]) for k ∈ 1:(length(y) ÷ 2)], dims)
end
deinterleaved(Δ::AbstractVector, dims) = reshape(complex.(Δ[1:2:end], Δ[2:2:end]), dims)


## The harmonics

function SphericalFunctions.sYlm_array(
    R::Rotor{<:TrackedReal}, ℓₘₐₓ::IT, s, ℓₘᵢₙ::IT
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
        (rotor_cotangent(R, sYlm_pullback(Y, Ȳ, ℓₘᵢₙ, ℓₘₐₓ))..., nothing, nothing, nothing)
    end
    (interleaved(Y), tracked_sYlm_array_pullback)
end


## Wigner's 𝔇

function SphericalFunctions.D_array(
    R::Rotor{<:TrackedReal}, ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT
) where {IT<:IntegerHalf}
    y = tracked_D_array(R[1], R[2], R[3], R[4], ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    assembled(y, (D_offset(ℓₘₐₓ + 1, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ),))
end

tracked_D_array(w::TrackedReal, x::TrackedReal, y::TrackedReal, z::TrackedReal, args...) =
    track(tracked_D_array, w, x, y, z, args...)

@grad function tracked_D_array(w, x, y, z, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    R = value_rotor(w, x, y, z)
    limits = (ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    Aʷ = D_array_widened(R, limits...)
    A = D_is_widened(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ) ? D_narrowed(Aʷ, limits...) : Aʷ
    function tracked_D_array_pullback(Δ)
        Ā = deinterleaved(Δ, size(A))
        (rotor_cotangent(R, D_pullback(Aʷ, Ā, limits...))..., ntuple(_ -> nothing, 5)...)
    end
    (interleaved(A), tracked_D_array_pullback)
end

end # module SphericalFunctionsReverseDiffExt
