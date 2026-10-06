module SphericalFunctionsEnzymeCoreExt

# Enzyme's rules for the calculators of 𝔇, of `d`, and of the harmonics, from the
# generators, as described in `src/derivatives/kernels.jl`.  They are defined in EnzymeCore,
# which is all that a rule needs.  `D`, `d`, `sYlm`, and `sYlm_matrix` are computed by
# calculators, which Enzyme follows, so these rules serve them too.
#
# A calculator keeps a copy of its rotors, or, for `d` and ₛλₗₘ, of their angles, and
# everything else it derives from them is computed by `set_rotor_data!`, which is declared
# here to have no derivatives.  So Enzyme differentiates only the copy, which puts the
# rotors' or angles' tangents, or later their cotangents, into the shadow calculator's copy;
# the angles of a calculator of floats are given their shadows by the rules for
# `store_angle!` below.  Each step of a calculator is a `compute_block!`, whose rules give
# the block's derivatives from those of the rotors in forward mode, and add the cotangents
# of the block into those of the rotors in reverse mode, so that Enzyme never differentiates
# the recurrence.  A calculator that Enzyme is to differentiate must be `Duplicated`, as any
# mutable workspace must: if it is created within the function being differentiated, Enzyme
# makes it so.
#
# Each step writes its block into a destination, the calculator's own buffer or an array
# such as the result of `D` or `sYlm_matrix`, and the rules follow the shadow of that
# destination.  In reverse mode, every step of a calculator writes its block into the same
# part of its buffer, so the cotangents that the code after the step accumulates in the
# shadow of that part belong to that step alone.  The reverse rule reads them, and then
# zeroes them, before the code that read the block of the step before it adds its own.  The
# values from which the derivatives are computed are kept on the tape, since the next step
# overwrites them.

import SphericalFunctions: Nᵣ, WignerCalculator, HarmonicCalculator, allocate_W, allocate_Y,
    compute_block!, set_rotor_data!, wigner_block_pushforward!, wigner_block_pullback!,
    harmonic_block_pushforward!, harmonic_block_pullback!, derivatives_from_left,
    derivative_m′range, derivative_mrange, m′range, mrange, rotor_generator, rotor_cotangent,
    block_array, derivative_values, angle_generators, add_angle_cotangents!, zero_cotangents,
    AngleCotangents, check_storage, store_angle!, angle_component_gradient
using Quaternionic: Quaternion
using EnzymeCore: EnzymeCore, EnzymeRules, Annotation, Const, Active, Duplicated, BatchDuplicated
using EnzymeCore.EnzymeRules: FwdConfig, RevConfig, AugmentedReturn, needs_primal,
    needs_shadow, width

# The calculators whose rotor data are floats, which are those these rules serve: those of
# 𝔇 and of `d`, and those of ₛYₗₘ and of ₛλₗₘ.
const DCalc = WignerCalculator{IT, RT, NT, ST, B, RT, Nothing} where {IT, RT, NT, ST, B}
const YCalc = HarmonicCalculator{IT, RT, NT, ST, S, B, RT, Nothing} where {IT, RT, NT, ST, S, B}

EnzymeRules.inactive(::typeof(set_rotor_data!), ::DCalc, ::Any) = nothing
EnzymeRules.inactive(::typeof(set_rotor_data!), ::YCalc, ::Any) = nothing

# `array_view(b)` compares the length of a block's vector with the block, in
# `check_storage`, whose branch throws.  Inlined into a loop over the elements of a block,
# that branch leads Enzyme's optimizer to peel the loop's first iteration, whose copy of the
# first element then crashes Enzyme in batched reverse mode (EnzymeAD/Enzyme#3316).  The
# check has no derivative, and a function declared inactive is not inlined into the code
# that Enzyme compiles, although it is everywhere else, so this declaration keeps the branch
# out of the loop at no cost elsewhere.  The view itself is left to Enzyme: rules for it,
# which would also keep it out of line, crash Enzyme's reverse mode on `d` in the garbage
# collector (as of Enzyme 0.13.209 on Julia 1.13).  The declaration may go once the issue is
# fixed.
EnzymeRules.inactive(::typeof(check_storage), ::Any) = nothing

# The shadows of an annotated argument, as a tuple of `width` of them, or `nothing` for a
# constant.  Any other annotation is an error, rather than being taken to have no shadow.
shadows(::Const) = nothing
shadows(x::Duplicated) = (x.dval,)
shadows(x::BatchDuplicated) = x.dval

# What a forward rule returns, given the primal and the tuple of its shadows.
function forward_return(config::FwdConfig, x, ẋ)
    if needs_shadow(config) && ẋ !== nothing
        if needs_primal(config)
            width(config) == 1 ? Duplicated(x, only(ẋ)) : BatchDuplicated(x, ẋ)
        else
            width(config) == 1 ? only(ẋ) : ẋ
        end
    else
        needs_primal(config) ? x : nothing
    end
end

# The generators of the tangents of `c`'s rotors that are held in the shadow `dc`, or of
# those of its angles, for a calculator of `d` or of ₛλₗₘ.
function generators(left::Bool, c, dc)
    G = Matrix{real(eltype(eltype(c.rotors)))}(undef, 3, length(c.rotors))
    for i ∈ eachindex(c.rotors)
        G[1, i], G[2, i], G[3, i] = rotor_generator(left, c.rotors[i], dc.rotors[i])
    end
    G
end

# Add the cotangents of the rotors from the vectors in the columns of `Ḡ` into the shadow
# `dc`'s copy of them, or those of the angles from the cotangents of their generators.
function add_rotor_cotangents!(dc, left::Bool, c, Ḡ::AbstractMatrix)
    for i ∈ eachindex(c.rotors)
        R̄ = rotor_cotangent(left, c.rotors[i], (Ḡ[1, i], Ḡ[2, i], Ḡ[3, i]))
        q = dc.rotors[i]
        dc.rotors[i] = Quaternion(q[1] + R̄[1], q[2] + R̄[2], q[3] + R̄[3], q[4] + R̄[4])
    end
    nothing
end
generators(
    left::Bool, c::Union{WignerCalculator{IT, RT, NT}, HarmonicCalculator{IT, RT, NT}}, dc
) where {IT, RT, NT<:Real} = angle_generators(dc.angles)
add_rotor_cotangents!(dc, left::Bool, c, Ḡ::AngleCotangents) =
    (add_angle_cotangents!(dc.angles, Ḡ); nothing)

# The shadow of the returned calculator, as the configuration asks for it.
shadow_return(config, c) = needs_shadow(config) ? (width(config) == 1 ? only(shadows(c)) : shadows(c)) : nothing


## The angles of a calculator of floats
#
# `store_angle!` does not compute the angle in a calculator of floats (see `store_angles!`),
# so its rules give the shadow of the angle from the components of the rotor data alone: the
# angle's tangent from theirs in forward mode, and in reverse mode their cotangents from the
# angle's, which the shadow then no longer holds, since the angle is overwritten.  A zero
# tangent or cotangent gives zero, even at a pole, where the gradient of the angle of a
# rotor is not finite.

function EnzymeRules.forward(
    config::FwdConfig, ::Const{typeof(store_angle!)}, ::Type{<:Annotation},
    angles::Annotation{<:AbstractVector{<:Base.IEEEFloat}}, i::Annotation,
    x::Annotation{<:Base.IEEEFloat}...
)
    if shadows(angles) !== nothing
        g = angle_component_gradient(map(a -> a.val, x)...)
        for (k, dangles) ∈ enumerate(shadows(angles))
            ẋ = map(a -> shadows(a) === nothing ? zero(a.val) : shadows(a)[k], x)
            dangles[i.val] = all(iszero, ẋ) ? zero(eltype(dangles)) : sum(g .* ẋ)
        end
    end
    nothing
end

function EnzymeRules.augmented_primal(
    ::RevConfig, ::Const{typeof(store_angle!)}, ::Type{<:Annotation},
    angles::Annotation{<:AbstractVector{<:Base.IEEEFloat}}, i::Annotation,
    x::Annotation{<:Base.IEEEFloat}...
)
    AugmentedReturn(nothing, nothing, nothing)
end

function EnzymeRules.reverse(
    config::RevConfig, ::Const{typeof(store_angle!)}, ::Type{<:Annotation}, tape,
    angles::Annotation{<:AbstractVector{<:Base.IEEEFloat}}, i::Annotation,
    x::Annotation{<:Base.IEEEFloat}...
)
    g = angle_component_gradient(map(a -> a.val, x)...)
    β̄ = if shadows(angles) === nothing
        ntuple(_ -> zero(eltype(angles.val)), Val(width(config)))
    else
        map(shadows(angles)) do dangles
            b = dangles[i.val]
            dangles[i.val] = zero(b)
            b
        end
    end
    x̄ = ntuple(Val(length(x))) do j
        if x[j] isa Active
            c = map(b -> iszero(b) ? zero(x[j].val) : oftype(x[j].val, g[j] * b), β̄)
            width(config) == 1 ? only(c) : c
        else
            nothing
        end
    end
    (nothing, nothing, x̄...)
end


## Allocation
#
# A new calculator is given a zeroed copy of itself as its shadow, rather than having Enzyme
# differentiate its allocation.  That allocation stores empty buffers — the half angles of
# an integer calculator, for example — in the fields of a struct, and Julia gives every empty
# buffer of a type the same constant object, which Enzyme cannot prove is never written
# through, and so would refuse without runtime activity.  Nothing about the allocation has a
# derivative: every argument is an index, a type or a size.

for allocate ∈ (allocate_W, allocate_Y)
    @eval begin
        function EnzymeRules.forward(
            config::FwdConfig, ::Const{typeof($allocate)}, ::Type{<:Annotation}, args::Const...
        )
            c = $allocate(map(a -> a.val, args)...)
            needs_shadow(config) || return needs_primal(config) ? c : nothing
            forward_return(config, c, ntuple(_ -> EnzymeCore.make_zero(c), Val(width(config))))
        end
        function EnzymeRules.augmented_primal(
            config::RevConfig, ::Const{typeof($allocate)}, ::Type{<:Annotation}, args::Const...
        )
            c = $allocate(map(a -> a.val, args)...)
            shadow = if needs_shadow(config)
                width(config) == 1 ? EnzymeCore.make_zero(c) :
                    ntuple(_ -> EnzymeCore.make_zero(c), Val(width(config)))
            else
                nothing
            end
            AugmentedReturn(needs_primal(config) ? c : nothing, shadow, nothing)
        end
        function EnzymeRules.reverse(
            ::RevConfig, ::Const{typeof($allocate)}, ::Type{<:Annotation}, tape, args::Const...
        )
            # One `nothing` for each argument, in a tuple whose length the compiler knows,
            # as Enzyme requires of a reverse rule.
            map(_ -> nothing, args)
        end
    end
end


## 𝔇

# The rows and columns from which the derivatives of the block of degree ℓ are computed, and
# those of the block itself.
function block_geometry(c::DCalc, ℓ)
    rows, cols = derivative_m′range(c, ℓ), derivative_mrange(c, ℓ)
    outrows, outcols = m′range(c, ℓ), mrange(c, ℓ)
    rows, cols, outrows, outcols
end

function EnzymeRules.forward(
    config::FwdConfig, ::Const{typeof(compute_block!)}, ::Type{<:Annotation},
    c::Annotation{<:DCalc}, ℓ::Annotation, A::Annotation, o::Annotation
)
    calc, l, j = c.val, ℓ.val, o.val
    compute_block!(calc, l, A.val, j)
    # The shadow block is written whether or not the returned calculator's shadow is asked
    # for, since the code after the step reads the block from the destination's shadow.
    if shadows(c) !== nothing && shadows(A) !== nothing
        left = derivatives_from_left(calc)
        rows, cols, outrows, outcols = block_geometry(calc, l)
        values = derivative_values(calc, l, A.val, j)
        for (dc, dA) ∈ zip(shadows(c), shadows(A))
            wigner_block_pushforward!(
                (x, ẋ) -> only(ẋ), block_array(calc, dA, l, j), values, l,
                rows, cols, outrows, outcols, left, generators(left, calc, dc), Val(1)
            )
        end
    end
    forward_return(config, calc, shadows(c))
end

function EnzymeRules.augmented_primal(
    config::RevConfig, ::Const{typeof(compute_block!)}, ::Type{<:Annotation},
    c::Annotation{<:DCalc}, ℓ::Annotation, A::Annotation, o::Annotation
)
    calc, l, j = c.val, ℓ.val, o.val
    compute_block!(calc, l, A.val, j)
    tape = copy(derivative_values(calc, l, A.val, j))
    AugmentedReturn(needs_primal(config) ? calc : nothing, shadow_return(config, c), tape)
end

function EnzymeRules.reverse(
    ::RevConfig, ::Const{typeof(compute_block!)}, ::Type{<:Annotation}, tape,
    c::Annotation{<:DCalc}, ℓ::Annotation, A::Annotation, o::Annotation
)
    calc, l, j = c.val, ℓ.val, o.val
    if shadows(c) !== nothing && shadows(A) !== nothing
        left = derivatives_from_left(calc)
        rows, cols, outrows, outcols = block_geometry(calc, l)
        for (dc, dA) ∈ zip(shadows(c), shadows(A))
            Ā = block_array(calc, dA, l, j)
            Ḡ = zero_cotangents(eltype(tape), Nᵣ(calc))
            wigner_block_pullback!(Ḡ, tape, Ā, l, rows, cols, outrows, outcols, left)
            add_rotor_cotangents!(dc, left, calc, Ḡ)
            fill!(Ā, zero(eltype(Ā)))
        end
    end
    (nothing, nothing, nothing, nothing)
end


## The harmonics

function EnzymeRules.forward(
    config::FwdConfig, ::Const{typeof(compute_block!)}, ::Type{<:Annotation},
    c::Annotation{<:YCalc}, ℓ::Annotation, is::Annotation, A::Annotation, o::Annotation
)
    calc, l, i, j = c.val, ℓ.val, is.val, o.val
    compute_block!(calc, l, i, A.val, j)
    # As for 𝔇, the shadow block is written whether or not the returned calculator's shadow
    # is asked for.
    if shadows(c) !== nothing && shadows(A) !== nothing
        Y = block_array(calc, A.val, l, i, j)
        for (dc, dA) ∈ zip(shadows(c), shadows(A))
            harmonic_block_pushforward!(
                (y, ẏ) -> only(ẏ), block_array(calc, dA, l, i, j), Y, l,
                generators(true, calc, dc), Val(1)
            )
        end
    end
    forward_return(config, calc, shadows(c))
end

function EnzymeRules.augmented_primal(
    config::RevConfig, ::Const{typeof(compute_block!)}, ::Type{<:Annotation},
    c::Annotation{<:YCalc}, ℓ::Annotation, is::Annotation, A::Annotation, o::Annotation
)
    calc, l, i, j = c.val, ℓ.val, is.val, o.val
    compute_block!(calc, l, i, A.val, j)
    tape = copy(block_array(calc, A.val, l, i, j))
    AugmentedReturn(needs_primal(config) ? calc : nothing, shadow_return(config, c), tape)
end

function EnzymeRules.reverse(
    ::RevConfig, ::Const{typeof(compute_block!)}, ::Type{<:Annotation}, tape,
    c::Annotation{<:YCalc}, ℓ::Annotation, is::Annotation, A::Annotation, o::Annotation
)
    calc, l, i, j = c.val, ℓ.val, is.val, o.val
    if shadows(c) !== nothing && shadows(A) !== nothing
        for (dc, dA) ∈ zip(shadows(c), shadows(A))
            Ȳ = block_array(calc, dA, l, i, j)
            Ḡ = zero_cotangents(eltype(tape), Nᵣ(calc))
            harmonic_block_pullback!(Ḡ, tape, Ȳ, l)
            add_rotor_cotangents!(dc, true, calc, Ḡ)
            fill!(Ȳ, zero(eltype(Ȳ)))
        end
    end
    (nothing, nothing, nothing, nothing, nothing)
end

end # module SphericalFunctionsEnzymeCoreExt
