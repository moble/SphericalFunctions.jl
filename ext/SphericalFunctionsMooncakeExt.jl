module SphericalFunctionsMooncakeExt

# Mooncake's rules for the calculators of 𝔇, of `d`, and of the harmonics, from the
# generators, as described in `src/derivatives/kernels.jl`, in forward and in reverse mode.
# They are the same rules as Enzyme's, and the comments at the top of that extension
# describe them: a calculator's copy of its rotors or angles is differentiated as Mooncake
# finds it, `set_rotor_data!` has no derivatives, and each `compute_block!` gives its
# block's derivatives from those of the rotors or angles, or adds its block's cotangents
# into theirs and zeroes them.  `D`, `d`, `sYlm`, and `sYlm_matrix` are computed by
# calculators, which Mooncake follows, so these rules serve them too.  The rules are
# primitives of Mooncake's `MinimalCtx`, the context for rules that are needed for
# correctness rather than for speed.
#
# Mooncake's pullbacks run with the primal state as it was when the corresponding call
# returned, so a rule that mutates memory must restore it: the pullback of an earlier read
# of the previous block may read the block's buffer again.  So the reverse rule for a step
# keeps the part of the destination that the step overwrites, and puts it back in its
# pullback.
#
# The tangent of a calculator is structural, as Mooncake builds it: the tangent, or in
# reverse mode the forward data, of each field in turn.  That of a quaternion is the tuple
# of its four components, inside the tangent of the `SVector` that holds them.

import SphericalFunctions: WignerCalculator, HarmonicCalculator, compute_block!,
    set_rotor_data!, wigner_block_pushforward!, wigner_block_pullback!,
    harmonic_block_pushforward!, harmonic_block_pullback!, derivatives_from_left,
    derivative_m′range, derivative_mrange, m′range, mrange, rotor_generator, rotor_cotangent,
    block_array, derivative_values, angle_generators, add_angle_cotangents!, zero_cotangents,
    Nᵣ, store_angle!, angle_component_gradient
import Mooncake
using Mooncake: MinimalCtx, @is_primitive, CoDual, NoRData, primal, tangent

const DCalc = WignerCalculator{IT, RT, NT, ST, B, RT, Nothing} where {IT, RT, NT, ST, B}
const YCalc = HarmonicCalculator{IT, RT, NT, ST, S, B, RT, Nothing} where {IT, RT, NT, ST, S, B}
const RealCalc = Union{WignerCalculator{IT, RT, NT}, HarmonicCalculator{IT, RT, NT}} where {IT, RT, NT<:Real}

# The components of the tangent of a quaternion, and the tangent with the given components.
quaternion_tangent_components(t::Mooncake.Tangent) = t.fields.components.fields.data
quaternion_tangent(q::NTuple{4}) = Mooncake.Tangent((components=Mooncake.Tangent((data=q,)),))

# The generators of the tangents of the rotors in the tangent `ċ` of a calculator, and the
# addition of their cotangents from the vectors in the columns of `Ḡ` into the forward data
# `c̄` of the calculator; or, for a calculator of `d` or of ₛλₗₘ, the same for its angles.
function generators(left::Bool, calc, ċ)
    ṙ = ċ.fields.rotors
    G = Matrix{real(eltype(eltype(calc.rotors)))}(undef, 3, length(calc.rotors))
    for i ∈ eachindex(calc.rotors)
        G[1, i], G[2, i], G[3, i] =
            rotor_generator(left, calc.rotors[i], quaternion_tangent_components(ṙ[i]))
    end
    G
end
function add_rotor_cotangents!(c̄, left::Bool, calc, Ḡ)
    r̄ = c̄.data.rotors
    for i ∈ eachindex(calc.rotors)
        R̄ = rotor_cotangent(left, calc.rotors[i], (Ḡ[1, i], Ḡ[2, i], Ḡ[3, i]))
        q = quaternion_tangent_components(r̄[i])
        r̄[i] = quaternion_tangent((q[1] + R̄[1], q[2] + R̄[2], q[3] + R̄[3], q[4] + R̄[4]))
    end
    nothing
end
generators(left::Bool, calc::RealCalc, ċ) = angle_generators(ċ.fields.angles)
add_rotor_cotangents!(c̄, left::Bool, calc::RealCalc, Ḡ) =
    (add_angle_cotangents!(c̄.data.angles, Ḡ); nothing)


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


# `store_angle!` does not compute the angle in a calculator of floats, and its rules give
# the shadow of the angle from the components of the rotor data alone, as Enzyme's do (see
# the comment on them in that extension).
@is_primitive MinimalCtx Tuple{typeof(store_angle!), Vector{<:Base.IEEEFloat}, Int, Vararg{Base.IEEEFloat}}

function Mooncake.frule!!(
    ::Mooncake.Dual{typeof(store_angle!)}, angles::Mooncake.Dual, i::Mooncake.Dual,
    x::Mooncake.Dual...
)
    ẋ = map(tangent, x)
    g = angle_component_gradient(map(primal, x)...)
    dangles = tangent(angles)
    dangles[primal(i)] = all(iszero, ẋ) ? zero(eltype(dangles)) : sum(g .* ẋ)
    Mooncake.zero_dual(nothing)
end
function Mooncake.rrule!!(
    ::CoDual{typeof(store_angle!)}, angles::CoDual, i::CoDual, x::CoDual...
)
    dangles, k = tangent(angles), primal(i)
    xs = map(primal, x)
    function store_angle_pullback(::NoRData)
        g = angle_component_gradient(xs...)
        β̄ = dangles[k]
        dangles[k] = zero(β̄)
        x̄ = map((xⱼ, gⱼ) -> iszero(β̄) ? zero(xⱼ) : oftype(xⱼ, gⱼ * β̄), xs, g)
        (NoRData(), NoRData(), NoRData(), x̄...)
    end
    (Mooncake.zero_fcodual(nothing), store_angle_pullback)
end


## 𝔇

function block_geometry(c::DCalc, ℓ)
    rows, cols = derivative_m′range(c, ℓ), derivative_mrange(c, ℓ)
    outrows, outcols = m′range(c, ℓ), mrange(c, ℓ)
    rows, cols, outrows, outcols
end

@is_primitive MinimalCtx Tuple{typeof(compute_block!), DCalc, Any, Any, Int}

function Mooncake.frule!!(
    ::Mooncake.Dual{typeof(compute_block!)}, c::Mooncake.Dual{<:DCalc}, ℓ::Mooncake.Dual,
    A::Mooncake.Dual, o::Mooncake.Dual
)
    calc, l, j = primal(c), primal(ℓ), primal(o)
    compute_block!(calc, l, primal(A), j)
    left = derivatives_from_left(calc)
    rows, cols, outrows, outcols = block_geometry(calc, l)
    wigner_block_pushforward!(
        (x, ẋ) -> only(ẋ), block_array(calc, tangent(A), l, j),
        derivative_values(calc, l, primal(A), j), l, rows, cols, outrows, outcols, left,
        generators(left, calc, tangent(c)), Val(1)
    )
    c
end

function Mooncake.rrule!!(
    ::CoDual{typeof(compute_block!)}, c::CoDual{<:DCalc}, ℓ::CoDual, A::CoDual, o::CoDual
)
    calc, l, j, c̄ = primal(c), primal(ℓ), primal(o), tangent(c)
    rows, cols, outrows, outcols = block_geometry(calc, l)
    previous = copy(block_array(calc, primal(A), l, j))
    compute_block!(calc, l, primal(A), j)
    values = copy(derivative_values(calc, l, primal(A), j))
    function compute_block_pullback(::NoRData)
        left = derivatives_from_left(calc)
        Ā = block_array(calc, tangent(A), l, j)
        Ḡ = zero_cotangents(eltype(values), Nᵣ(calc))
        wigner_block_pullback!(Ḡ, values, Ā, l, rows, cols, outrows, outcols, left)
        add_rotor_cotangents!(c̄, left, calc, Ḡ)
        fill!(Ā, zero(eltype(Ā)))
        copyto!(block_array(calc, primal(A), l, j), previous)
        (NoRData(), NoRData(), NoRData(), NoRData(), NoRData())
    end
    (c, compute_block_pullback)
end


## The harmonics

@is_primitive MinimalCtx Tuple{typeof(compute_block!), YCalc, Any, Any, Any, Int}

function Mooncake.frule!!(
    ::Mooncake.Dual{typeof(compute_block!)}, c::Mooncake.Dual{<:YCalc}, ℓ::Mooncake.Dual,
    is::Mooncake.Dual, A::Mooncake.Dual, o::Mooncake.Dual
)
    calc, l, i, j = primal(c), primal(ℓ), primal(is), primal(o)
    compute_block!(calc, l, i, primal(A), j)
    Y = block_array(calc, primal(A), l, i, j)
    harmonic_block_pushforward!(
        (y, ẏ) -> only(ẏ), block_array(calc, tangent(A), l, i, j), Y, l,
        generators(true, calc, tangent(c)), Val(1)
    )
    c
end

function Mooncake.rrule!!(
    ::CoDual{typeof(compute_block!)}, c::CoDual{<:YCalc}, ℓ::CoDual, is::CoDual,
    A::CoDual, o::CoDual
)
    calc, l, i, j = primal(c), primal(ℓ), primal(is), primal(o)
    previous = copy(block_array(calc, primal(A), l, i, j))
    compute_block!(calc, l, i, primal(A), j)
    values = copy(block_array(calc, primal(A), l, i, j))
    function compute_block_pullback(::NoRData)
        Ȳ = block_array(calc, tangent(A), l, i, j)
        Ḡ = zero_cotangents(eltype(values), Nᵣ(calc))
        harmonic_block_pullback!(Ḡ, values, Ȳ, l)
        add_rotor_cotangents!(tangent(c), true, calc, Ḡ)
        fill!(Ȳ, zero(eltype(Ȳ)))
        copyto!(block_array(calc, primal(A), l, i, j), previous)
        (NoRData(), NoRData(), NoRData(), NoRData(), NoRData(), NoRData())
    end
    (c, compute_block_pullback)
end

end # module SphericalFunctionsMooncakeExt
