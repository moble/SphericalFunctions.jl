module SphericalFunctionsMooncakeExt

# Mooncake's rules for the calculators of 𝔇 and of the harmonics, from the generators, as
# described in `src/derivatives.jl`, in forward and in reverse mode.  They are the same
# rules as Enzyme's, and the comments at the top of that extension describe them: a
# calculator's copy of its rotors is differentiated as Mooncake finds it, `set_rotor_data!`
# has no derivatives, and each `compute_block!` gives its block's derivatives from those of
# the rotors, or adds its block's cotangents into the rotors' and zeroes them.  `D`, `sYlm`,
# and `sYlm_matrix` are computed by calculators, which Mooncake follows, so these rules
# serve them too.  The rules are primitives of Mooncake's `MinimalCtx`, the context for
# rules that are needed for correctness rather than for speed.
#
# Mooncake's pullbacks run with the primal state as it was when the corresponding call
# returned, so a rule that mutates memory must restore it: the pullback of an earlier read
# of the previous block may read the block's buffer again.  So the reverse rule for a step
# keeps the part of the buffer that the step overwrites, and puts it back in its pullback.
#
# The tangent of a calculator is structural, as Mooncake builds it: the tangent, or in
# reverse mode the forward data, of each field in turn.  That of a quaternion is the tuple
# of its four components, inside the tangent of the `SVector` that holds them.

import SphericalFunctions: WignerCalculator, HarmonicCalculator, compute_block!,
    set_rotor_data!, wigner_block_pushforward!, wigner_block_pullback!,
    harmonic_block_pushforward!, harmonic_block_pullback!, derivatives_from_left,
    stored_m′range, stored_mrange, m′range, mrange, rotor_generator, rotor_cotangent
import Mooncake
using Mooncake: MinimalCtx, @is_primitive, CoDual, NoRData, primal, tangent

const DCalc = WignerCalculator{IT, RT, Complex{RT}, ST, B, RT, Nothing} where {IT, RT, ST, B}
const YCalc = HarmonicCalculator{IT, RT, Complex{RT}, ST, S, B, RT, Nothing} where {IT, RT, ST, S, B}

# The components of the tangent of a quaternion, and the tangent with the given components.
quaternion_tangent_components(t::Mooncake.Tangent) = t.fields.components.fields.data
quaternion_tangent(q::NTuple{4}) = Mooncake.Tangent((components=Mooncake.Tangent((data=q,)),))

# The generators of the rotors' tangents in `ṙ`, and the addition of their cotangents from
# the vectors in the columns of `Ḡ` into `r̄`.
function generators(left::Bool, calc, ṙ::AbstractVector)
    G = Matrix{real(eltype(eltype(calc.rotors)))}(undef, 3, length(calc.rotors))
    for i ∈ eachindex(calc.rotors)
        G[1, i], G[2, i], G[3, i] =
            rotor_generator(left, calc.rotors[i], quaternion_tangent_components(ṙ[i]))
    end
    G
end
function add_rotor_cotangents!(r̄::AbstractVector, left::Bool, calc, Ḡ)
    for i ∈ eachindex(calc.rotors)
        R̄ = rotor_cotangent(left, calc.rotors[i], (Ḡ[1, i], Ḡ[2, i], Ḡ[3, i]))
        q = quaternion_tangent_components(r̄[i])
        r̄[i] = quaternion_tangent((q[1] + R̄[1], q[2] + R̄[2], q[3] + R̄[3], q[4] + R̄[4]))
    end
    nothing
end


## The rotor data

@is_primitive MinimalCtx Tuple{typeof(set_rotor_data!), DCalc, Any}
@is_primitive MinimalCtx Tuple{typeof(set_rotor_data!), YCalc, Any}

function Mooncake.frule!!(
    ::Mooncake.Dual{typeof(set_rotor_data!)}, c::Mooncake.Dual{<:Union{DCalc, YCalc}},
    R::Mooncake.Dual
)
    set_rotor_data!(primal(c), primal(R))
    Mooncake.zero_dual(nothing)
end
function Mooncake.rrule!!(
    ::CoDual{typeof(set_rotor_data!)}, c::CoDual{<:Union{DCalc, YCalc}}, R::CoDual
)
    set_rotor_data!(primal(c), primal(R))
    rdata = (Mooncake.zero_rdata(primal(c)), Mooncake.zero_rdata(primal(R)))
    set_rotor_data_pullback(::NoRData) = (NoRData(), rdata...)
    (Mooncake.zero_fcodual(nothing), set_rotor_data_pullback)
end


## 𝔇

function block_geometry(c::DCalc, ℓ)
    rows, cols = stored_m′range(c, ℓ), stored_mrange(c, ℓ)
    outrows, outcols = m′range(c, ℓ), mrange(c, ℓ)
    rows, cols, outrows, outcols
end
# The part of `Wˡ` that holds the block that `recurrence!` returns, within what is stored.
function labelled_view(Wˡ, rows, cols, outrows, outcols)
    o′, o = Int(first(outrows) - first(rows)), Int(first(outcols) - first(cols))
    view(Wˡ, :, (o′ + 1):(o′ + length(outrows)), (o + 1):(o + length(outcols)))
end
stored_view(Wˡ, rows, cols) = view(Wˡ, :, 1:length(rows), 1:length(cols))

@is_primitive MinimalCtx Tuple{typeof(compute_block!), DCalc, Any}

function Mooncake.frule!!(
    ::Mooncake.Dual{typeof(compute_block!)}, c::Mooncake.Dual{<:DCalc}, ℓ::Mooncake.Dual
)
    calc, l, ċ = primal(c), primal(ℓ), tangent(c)
    compute_block!(calc, l)
    left = derivatives_from_left(calc)
    rows, cols, outrows, outcols = block_geometry(calc, l)
    wigner_block_pushforward!(
        (x, ẋ) -> only(ẋ), labelled_view(ċ.fields.Wˡ, rows, cols, outrows, outcols),
        stored_view(calc.Wˡ, rows, cols), l, rows, cols, outrows, outcols, left,
        generators(left, calc, ċ.fields.rotors), Val(1)
    )
    c
end

function Mooncake.rrule!!(
    ::CoDual{typeof(compute_block!)}, c::CoDual{<:DCalc}, ℓ::CoDual
)
    calc, l, c̄ = primal(c), primal(ℓ), tangent(c)
    rows, cols, outrows, outcols = block_geometry(calc, l)
    previous = copy(stored_view(calc.Wˡ, rows, cols))
    compute_block!(calc, l)
    values = copy(stored_view(calc.Wˡ, rows, cols))
    function compute_block_pullback(::NoRData)
        left = derivatives_from_left(calc)
        Ā = labelled_view(c̄.data.Wˡ, rows, cols, outrows, outcols)
        Ḡ = zeros(real(eltype(values)), 3, length(calc.rotors))
        wigner_block_pullback!(Ḡ, values, Ā, l, rows, cols, outrows, outcols, left)
        add_rotor_cotangents!(c̄.data.rotors, left, calc, Ḡ)
        fill!(Ā, zero(eltype(Ā)))
        copyto!(stored_view(calc.Wˡ, rows, cols), previous)
        (NoRData(), NoRData(), NoRData())
    end
    (c, compute_block_pullback)
end


## The harmonics

block_view(Y, ℓ, j₀) = view(Y, :, :, (j₀ + 1):(j₀ + Int(2ℓ) + 1))

@is_primitive MinimalCtx Tuple{typeof(compute_block!), YCalc, Any, Any, Any, Int}

function Mooncake.frule!!(
    ::Mooncake.Dual{typeof(compute_block!)}, c::Mooncake.Dual{<:YCalc}, ℓ::Mooncake.Dual,
    is::Mooncake.Dual, Y::Mooncake.Dual, j₀::Mooncake.Dual
)
    calc, l, i, j = primal(c), primal(ℓ), primal(is), primal(j₀)
    compute_block!(calc, l, i, primal(Y), j)
    harmonic_block_pushforward!(
        (y, ẏ) -> only(ẏ), block_view(tangent(Y), l, j), block_view(primal(Y), l, j), l, i,
        generators(true, calc, tangent(c).fields.rotors), Val(1)
    )
    c
end

function Mooncake.rrule!!(
    ::CoDual{typeof(compute_block!)}, c::CoDual{<:YCalc}, ℓ::CoDual, is::CoDual,
    Y::CoDual, j₀::CoDual
)
    calc, l, i, j = primal(c), primal(ℓ), primal(is), primal(j₀)
    previous = copy(block_view(primal(Y), l, j))
    compute_block!(calc, l, i, primal(Y), j)
    values = copy(block_view(primal(Y), l, j))
    function compute_block_pullback(::NoRData)
        Ȳ = block_view(tangent(Y), l, j)
        Ḡ = zeros(real(eltype(values)), 3, length(calc.rotors))
        harmonic_block_pullback!(Ḡ, values, Ȳ, l, i)
        add_rotor_cotangents!(tangent(c).data.rotors, true, calc, Ḡ)
        fill!(view(Ȳ, :, i, :), zero(eltype(Ȳ)))
        copyto!(block_view(primal(Y), l, j), previous)
        (NoRData(), NoRData(), NoRData(), NoRData(), NoRData(), NoRData())
    end
    (c, compute_block_pullback)
end

end # module SphericalFunctionsMooncakeExt
