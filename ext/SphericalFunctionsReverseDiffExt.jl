module SphericalFunctionsReverseDiffExt

# ReverseDiff's rules for `D_array`, `sYlm_array`, and `sYlm_matrix_array`, from the
# generators, as described in `src/derivatives/kernels.jl`.  ReverseDiff tracks only real
# numbers and arrays of them, and a rule is an instruction recorded on its tape by `@grad`,
# whose inputs and output must be of those kinds.  So a rotor is passed to the instruction
# as its four tracked components, and the instruction's output is the real array of the
# values' real and imaginary parts, interleaved.  The complex values are then assembled from
# the elements of that tracked array, each of which sends its derivative back to the array,
# and so to the instruction.  The harmonics of a vector of rotors are recorded as one
# instruction for each rotor, since an instruction with every rotor's components as its
# inputs would have to be compiled anew for every number of rotors.
#
# ReverseDiff does not read ChainRules' `rrule` for a function unless it is told to with
# `@grad_from_chainrules`, which requires the arguments to be tracked arrays or reals, as a
# rotor is not; hence the separate rules.
#
# A calculator of tracked rotor data runs the recurrence in a calculator of their values, as
# a calculator of ForwardDiff's dual numbers does (see `src/derivatives/lifting.jl`), and
# each of its steps records one instruction, whose output is the block of every rotor and
# whose pullback is the kernel of `src/derivatives/kernels.jl`.  The elements of a block are
# assembled from the instruction's output, as above, and so are of the type of an element of
# a tracked `Vector`; a calculator stores only that type, and its rotors' components, or its
# angles, are passed through an instruction of their own, so that they are of that type too,
# whatever the type of the tracked numbers they were given as.
#
# `d` of tracked rotor data is one instruction too, whose inputs are the data's tracked
# components, and whose pullback gives their cotangents from that of the data's angle (see
# `rotation_angle`).  A calculator of `d` or of ₛλₗₘ keeps a copy of its angles, as tracked
# numbers, and the instructions of its steps take those angles.
#
# ReverseDiff replays a recorded tape, as a compiled tape is replayed for every gradient, by
# calling each instruction's function again on the new values of its inputs, but it runs the
# pullback that the first call returned.  So every instruction here is given a `Ref` as an
# input of its own, into which every call writes what its pullback reads, and the pullback
# of the first call reads what the latest call wrote.  The instruction of a calculator's
# step computes its block from the values of its inputs, as described under "The steps"
# below, so that the replay of a loop over a calculator is right too.

import SphericalFunctions
import SphericalFunctions: D_array, sYlm_array, sYlm_matrix_array,
    wigner_arrays_with_derivative_values, wigner_arrays_pullback, harmonic_array_pullback!,
    derivatives_from_left, rotor_cotangent, rotor_cotangents, IntegerHalf, Ysize, nspins,
    WignerCalculator, HarmonicCalculator, Lift, wigner_block_pullback!,
    harmonic_block_pullback!, m′range, mrange, block_array, rotation_angle,
    rotation_angle_cotangent, zero_cotangents, floattype, lowest_index, angle_cotangent,
    AngleCotangents, value_type, compute_block!, holds_block, set_rotors!, rotor_data_components
using Quaternionic: Quaternionic, Rotor, Quaternion
import ReverseDiff
using ReverseDiff: TrackedReal, TrackedArray, @grad, value, track

const TrackedRotorLike = Union{Rotor{<:TrackedReal}, Quaternion{<:TrackedReal}}

# The element type of a tracked `Vector`, which is what a calculator stores.
const TrackedElement{V, D} = TrackedReal{V, D, TrackedArray{V, D, 1, Vector{V}, Vector{D}}}

SphericalFunctions.value_type(::Type{<:TrackedReal{V}}) where {V} = V
SphericalFunctions.real_value(x::TrackedReal) = value(x)
SphericalFunctions.ndirections(::Type{<:TrackedReal}) = 0
SphericalFunctions.floattype(::Type{<:TrackedReal{V}}) where {V} = TrackedElement{V, V}

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
    y = tracked_D_array(R[1], R[2], R[3], R[4], ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, Ref{Any}())
    o = 0
    map(SphericalFunctions.lowest_index(IT):ℓₘₐₓ) do ℓ
        dims = (length(max(-ℓ, m′ₘᵢₙ):min(ℓ, m′ₘₐₓ)), length(max(-ℓ, mₘᵢₙ):min(ℓ, mₘₐₓ)))
        b = assembled(y, dims, o)
        o += 2prod(dims)
        b
    end
end

tracked_D_array(w::TrackedReal, x::TrackedReal, y::TrackedReal, z::TrackedReal, args...) =
    track(tracked_D_array, w, x, y, z, args...)

@grad function tracked_D_array(w, x, y, z, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, state)
    R = value_rotor(w, x, y, z)
    blocks, values, calc = wigner_arrays_with_derivative_values(
        Complex{floattype(R)}, R, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ
    )
    state[] = (R, blocks, values, calc)
    tracked_D_array_pullback(Δ) = D_array_cotangents(state[], Δ)
    (reduce(vcat, [interleaved(b) for b ∈ blocks]), tracked_D_array_pullback)
end

# The cotangents of the inputs of `tracked_D_array`, from what its latest call wrote into its
# state.
function D_array_cotangents((R, blocks, values, calc), Δ)
    o = 0
    Ā = map(blocks) do b
        Āᵢ = deinterleaved(Δ, size(b), o)
        o += 2length(b)
        Āᵢ
    end
    Ḡ = wigner_arrays_pullback(calc, values, Ā)
    g = (Ḡ[1], Ḡ[2], Ḡ[3])
    (rotor_cotangent(derivatives_from_left(calc), R, g)..., ntuple(_ -> nothing, 6)...)
end


## d
#
# The instruction for `d` of tracked rotor data takes the data's components, one for an
# angle, two for a phase, and four for a rotor, and computes the values from the data itself,
# as `d` of the data's values would.

const TrackedRotorData = Union{TrackedReal, Complex{<:TrackedReal}, TrackedRotorLike}
# The rotor data of the values of its components, and the cotangents of its components
# from a cotangent of the data as `rotation_angle_cotangent` gives it.
data_value(β) = value(β)
data_value(x, y) = Complex(value(x), value(y))
data_value(w, x, y, z) = value_rotor(w, x, y, z)
component_cotangents(x̄::Real) = (x̄,)
component_cotangents(x̄::Complex) = (real(x̄), imag(x̄))
component_cotangents(x̄::Tuple) = x̄

function SphericalFunctions.d_array(
    x::TrackedRotorData, ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT
) where {IT<:IntegerHalf}
    y = tracked_d_array(rotor_data_components(x)..., ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, Ref{Any}())
    o = 0
    map(lowest_index(IT):ℓₘₐₓ) do ℓ
        dims = (length(max(-ℓ, m′ₘᵢₙ):min(ℓ, m′ₘₐₓ)), length(max(-ℓ, mₘᵢₙ):min(ℓ, mₘₐₓ)))
        b = reshape([y[o + k] for k ∈ 1:prod(dims)], dims)
        o += prod(dims)
        b
    end
end

tracked_d_array(β::TrackedReal, args...) = track(tracked_d_array, β, args...)
tracked_d_array(x::TrackedReal, y::TrackedReal, args...) = track(tracked_d_array, x, y, args...)
tracked_d_array(w::TrackedReal, x::TrackedReal, y::TrackedReal, z::TrackedReal, args...) =
    track(tracked_d_array, w, x, y, z, args...)

@grad tracked_d_array(β, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, state) =
    d_array_rule(data_value(β), (ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ), state)
@grad tracked_d_array(x, y, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, state) =
    d_array_rule(data_value(x, y), (ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ), state)
@grad tracked_d_array(w, x, y, z, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, state) =
    d_array_rule(data_value(w, x, y, z), (ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ), state)

# The output and the pullback of `tracked_d_array` of the rotor data `x`
function d_array_rule(x, limits, state)
    blocks, values, calc = wigner_arrays_with_derivative_values(floattype(x), x, limits...)
    state[] = (x, blocks, values, calc)
    tracked_d_array_pullback(Δ) = d_array_cotangents(state[], Δ)
    (reduce(vcat, [vec(b) for b ∈ blocks]), tracked_d_array_pullback)
end
function d_array_cotangents((x, blocks, values, calc), Δ)
    o = 0
    Ā = map(blocks) do b
        Āᵢ = reshape(Δ[(o + 1):(o + length(b))], size(b))
        o += length(b)
        Āᵢ
    end
    β̄ = angle_cotangent(wigner_arrays_pullback(calc, values, Ā), 1)
    (component_cotangents(rotation_angle_cotangent(x, β̄))..., ntuple(_ -> nothing, 6)...)
end


## The harmonics

function SphericalFunctions.sYlm_array(
    R::TrackedRotorLike, ℓₘₐₓ::IT, s, ℓₘᵢₙ::IT
) where {IT<:IntegerHalf}
    y = tracked_sYlm_array(R[1], R[2], R[3], R[4], ℓₘₐₓ, s, ℓₘᵢₙ, Ref{Any}())
    n = Ysize(ℓₘᵢₙ, ℓₘₐₓ)
    assembled(y, s isa AbstractUnitRange ? (length(s), n) : (n,))
end

tracked_sYlm_array(w::TrackedReal, x::TrackedReal, y::TrackedReal, z::TrackedReal, args...) =
    track(tracked_sYlm_array, w, x, y, z, args...)

@grad function tracked_sYlm_array(w, x, y, z, ℓₘₐₓ, s, ℓₘᵢₙ, state)
    R = value_rotor(w, x, y, z)
    Y = sYlm_array(R, ℓₘₐₓ, s, ℓₘᵢₙ)
    state[] = (R, Y, ℓₘᵢₙ, ℓₘₐₓ)
    tracked_sYlm_array_pullback(Δ) = sYlm_array_cotangents(state[], Δ)
    (interleaved(Y), tracked_sYlm_array_pullback)
end
function sYlm_array_cotangents((R, Y, ℓₘᵢₙ, ℓₘₐₓ), Δ)
    Ȳ = deinterleaved(Δ, size(Y))
    Ḡ = harmonic_array_pullback!(zeros(real(eltype(Y)), 3, 1), Y, Ȳ, false, ℓₘᵢₙ, ℓₘₐₓ)
    (only(rotor_cotangents(true, [R], Ḡ))..., nothing, nothing, nothing, nothing)
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

# The rotors of the points (θ, 0), (cos(θ/2), 0, sin(θ/2), 0), built from their components,
# since a compiled tape cannot replay the broadcast into an `SVector` of Quaternionic's
# `from_spherical_coordinates`.
function SphericalFunctions.store_point_rotors!(
    rotors::AbstractVector{Quaternion{T}}, θ
) where {T<:TrackedReal}
    for i ∈ eachindex(rotors)
        s, c = sincos(θ[i] / 2)
        q = tracked_components(c, zero(c), s, zero(c))
        rotors[i] = Quaternion{T}(q[1], q[2], q[3], q[4])
    end
    rotors
end

# The angles of a calculator of `d` or of ₛλₗₘ, each passed through an instruction that
# returns it as a tracked `Vector` of one element, which is of the type the calculator
# stores.
function SphericalFunctions.store_angles!(
    angles::AbstractVector{T}, R::AbstractVector
) where {T<:TrackedReal}
    for i ∈ eachindex(angles, R)
        angles[i] = tracked_angle(rotation_angle(R[i]))[1]
    end
    angles
end
function SphericalFunctions.store_angles!(angles::AbstractVector{T}, R) where {T<:TrackedReal}
    angles[1] = tracked_angle(rotation_angle(R))[1]
    angles
end

tracked_angle(β::TrackedReal) = track(tracked_angle, β)

@grad function tracked_angle(β)
    tracked_angle_pullback(Δ) = (Δ[1],)
    ([value(β)], tracked_angle_pullback)
end

# The steps
#
# A step of a calculator of tracked rotor data records one instruction, whose inputs are the
# components of all of the calculator's rotors, or its angles, as one tracked vector, and
# whose output is the block of every rotor.  The instruction brings the calculator of values
# to the degree ℓ at the values of those inputs itself: it gives that calculator the rotor
# data, unless it holds them already, and runs its recurrence to ℓ, which advances it by one
# degree when it is at ℓ - 1, as it is in a loop over the blocks.  So the instruction gives
# the right block when the tape is replayed at other values, whatever the order in which the
# blocks were computed as the tape was recorded.  For a calculator of tracked numbers, the
# matrix `G` of its `Lift`, which holds no generators, holds the values of the rotor data
# that the calculator of values was last given: the four components of each rotor, or each
# angle.  A calculator of `d` given phases or rotors gives its calculator of values their
# angles when a replay changes them.

SphericalFunctions.allocate_lift(::Type{RT}, ::Type{NT}, inner, Nᵣ::Int) where {RT<:TrackedReal, NT} =
    Lift(inner, fill(value_type(RT)(NaN), NT <: Complex ? 4 : 1, Nᵣ))

function SphericalFunctions.set_generators!(
    lift::Lift, ::Bool, rotors::AbstractVector{Quaternion{RT}}
) where {RT<:TrackedReal}
    for i ∈ eachindex(rotors), k ∈ 1:4
        lift.G[k, i] = value(rotors[i][k])
    end
    lift
end
function SphericalFunctions.set_generators!(
    lift::Lift, ::Bool, angles::AbstractVector{RT}
) where {RT<:TrackedReal}
    for i ∈ eachindex(angles)
        lift.G[1, i] = value(angles[i])
    end
    lift
end

# The components of a calculator's rotors, or its angles, as one vector, which is the input
# of each step's instruction.
data_vector(c) = isempty(c.rotors) ? copy(c.angles) : [q[k] for q ∈ c.rotors for k ∈ 1:4]

# Give the calculator of values of `lift` the rotor data whose values are `x`, as
# `data_vector` orders them, unless it holds them already.
function set_values!(lift::Lift, x::AbstractVector)
    G = lift.G
    if !isequal(vec(G), x)
        copyto!(G, x)
        set_rotors!(lift.inner, rotor_data_values(G))
    end
    lift
end

# The rotor data whose values `G` holds
function rotor_data_values(G::AbstractMatrix)
    size(G, 1) == 1 && return vec(copy(G))
    [Quaternion(G[1, i], G[2, i], G[3, i], G[4, i]) for i ∈ axes(G, 2)]
end

# The cotangents of the components in `x` from the cotangents of the generators, or of the
# angles, that a kernel gives.
data_cotangents(x, Ḡ::AbstractMatrix, left::Bool) = reduce(vcat, [
    collect(rotor_cotangent(left, view(x, (4i - 3):4i), (Ḡ[1, i], Ḡ[2, i], Ḡ[3, i])))
    for i ∈ axes(Ḡ, 2)
])
data_cotangents(x, Ḡ::AngleCotangents, ::Bool) = [angle_cotangent(Ḡ, i) for i ∈ eachindex(x)]

# The values as a real vector, and the values of the given size from such a vector
real_vector(A::AbstractArray{<:Complex}) = interleaved(A)
real_vector(A::AbstractArray{<:Real}) = vec(A)
from_real_vector(::Type{<:Complex}, Δ, dims) = deinterleaved(Δ, dims)
from_real_vector(::Type{<:Real}, Δ, dims) = reshape(Δ[1:prod(dims)], dims)

# Write the elements of the instruction's output `y` into the block `out`.
function assemble!(out::AbstractArray{<:Complex}, y)
    for k ∈ eachindex(out)
        out[k] = Complex(y[2k - 1], y[2k])
    end
    out
end
function assemble!(out::AbstractArray{<:Real}, y)
    for k ∈ eachindex(out)
        out[k] = y[k]
    end
    out
end

function SphericalFunctions.compute_block!(
    c::WignerCalculator{IT, RT, NT, ST, B, FT, <:Lift}, ℓ::IT, A, o::Int
) where {IT, RT<:TrackedReal, NT, ST, B, FT<:Real}
    inner = c.lift.inner
    y = tracked_wigner_step(
        data_vector(c), c.lift, ℓ, m′range(inner, ℓ), mrange(inner, ℓ), m′range(c, ℓ),
        mrange(c, ℓ), derivatives_from_left(c), Ref{Any}()
    )
    assemble!(block_array(c, A, ℓ, o), y)
    c.ℓ[] = holds_block(c, A, o) ? ℓ : lowest_index(IT) - 1
    c
end

tracked_wigner_step(x::AbstractVector{<:TrackedReal}, args...) = track(tracked_wigner_step, x, args...)

@grad function tracked_wigner_step(x, lift, ℓ, rows, cols, outrows, outcols, left, state)
    xᵥ = value(x)
    inner = set_values!(lift, xᵥ).inner
    compute_block!(inner, ℓ)
    values = copy(block_array(inner, inner.Wˡ, ℓ))
    o′, o = Int(first(outrows) - first(rows)), Int(first(outcols) - first(cols))
    block = values[:, (o′ + 1):(o′ + length(outrows)), (o + 1):(o + length(outcols))]
    state[] = (xᵥ, values)
    function tracked_wigner_step_pullback(Δ)
        xᵥ, values = state[]
        Ā = from_real_vector(eltype(values), Δ, size(block))
        Ḡ = zero_cotangents(eltype(values), size(values, 1))
        wigner_block_pullback!(Ḡ, values, Ā, ℓ, rows, cols, outrows, outcols, left)
        (data_cotangents(xᵥ, Ḡ, left), ntuple(_ -> nothing, 8)...)
    end
    (real_vector(block), tracked_wigner_step_pullback)
end

function SphericalFunctions.compute_block!(
    c::HarmonicCalculator{IT, RT, NT, ST, S, B, FT, <:Lift}, ℓ::IT, is, A, o::Int
) where {IT, RT<:TrackedReal, NT, ST, S, B, FT<:Real}
    y = tracked_harmonic_step(data_vector(c), c.lift, ℓ, is, Ref{Any}())
    assemble!(block_array(c, A, ℓ, is, o), y)
    c.ℓ[] = holds_block(c, is, A, o) ? ℓ : lowest_index(IT) - 1
    c
end

tracked_harmonic_step(x::AbstractVector{<:TrackedReal}, args...) =
    track(tracked_harmonic_step, x, args...)

@grad function tracked_harmonic_step(x, lift, ℓ, is, state)
    xᵥ = value(x)
    inner = set_values!(lift, xᵥ).inner
    compute_block!(inner, ℓ, is, inner.Yˡ, 0)
    values = copy(block_array(inner, inner.Yˡ, ℓ, is))
    state[] = (xᵥ, values)
    function tracked_harmonic_step_pullback(Δ)
        xᵥ, values = state[]
        Ȳ = from_real_vector(eltype(values), Δ, size(values))
        Ḡ = harmonic_block_pullback!(zero_cotangents(eltype(values), size(values, 1)), values, Ȳ, ℓ)
        (data_cotangents(xᵥ, Ḡ, true), nothing, nothing, nothing, nothing)
    end
    (real_vector(values), tracked_harmonic_step_pullback)
end

end # module SphericalFunctionsReverseDiffExt
