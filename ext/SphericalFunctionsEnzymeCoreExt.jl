module SphericalFunctionsEnzymeCoreExt

# Enzyme's rules for the calculators of 𝔇 and of the harmonics, from the generators, as
# described in `src/derivatives.jl`.  They are defined in EnzymeCore, which is all that a
# rule needs.  `D`, `sYlm`, and `sYlm_matrix` are computed by calculators, which Enzyme
# follows, so these rules serve them too.
#
# A calculator keeps a copy of its rotors, and everything else it derives from them is
# computed by `set_rotor_data!`, which is declared here to have no derivatives.  So Enzyme
# differentiates only the copy, which puts the rotors' tangents, or later their cotangents,
# into the shadow calculator's copy.  Each step of a calculator is a `compute_block!`, whose
# rules give the block's derivatives from those of the rotors in forward mode, and add the
# cotangents of the block into those of the rotors in reverse mode, so that Enzyme never
# differentiates the recurrence.  A calculator that Enzyme is to differentiate must be
# `Duplicated`, as any mutable workspace must: if it is created within the function being
# differentiated, Enzyme makes it so.
#
# In reverse mode, every step of a calculator writes its block into the same buffer, so the
# cotangents that the code after the step accumulates in the shadow of that buffer belong to
# that step alone.  The reverse rule reads them, and then zeroes them, before the code that
# read the block of the step before it adds its own.  The values of the block are kept on
# the tape, since the next step overwrites them.

import SphericalFunctions: WignerCalculator, HarmonicCalculator, allocate_W, allocate_Y,
    compute_block!, set_rotor_data!, wigner_block_pushforward!, wigner_block_pullback!,
    harmonic_block_pushforward!, harmonic_block_pullback!, derivatives_from_left,
    stored_m′range, stored_mrange, m′range, mrange, rotor_generator, rotor_cotangent
using Quaternionic: Quaternion
using EnzymeCore: EnzymeCore, EnzymeRules, Annotation, Const, Duplicated, BatchDuplicated
using EnzymeCore.EnzymeRules: FwdConfig, RevConfig, AugmentedReturn, needs_primal,
    needs_shadow, width

# The complex calculators whose rotor data are floats, which are those these rules serve.
const DCalc = WignerCalculator{IT, RT, Complex{RT}, ST, B, RT, Nothing} where {IT, RT, ST, B}
const YCalc = HarmonicCalculator{IT, RT, Complex{RT}, ST, S, B, RT, Nothing} where {IT, RT, ST, S, B}

EnzymeRules.inactive(::typeof(set_rotor_data!), ::DCalc, ::Any) = nothing
EnzymeRules.inactive(::typeof(set_rotor_data!), ::YCalc, ::Any) = nothing

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

# The generators of the tangents of `c`'s rotors that are held in the shadow `dc`.
function generators(left::Bool, c, dc)
    G = Matrix{real(eltype(eltype(c.rotors)))}(undef, 3, length(c.rotors))
    for i ∈ eachindex(c.rotors)
        G[1, i], G[2, i], G[3, i] = rotor_generator(left, c.rotors[i], dc.rotors[i])
    end
    G
end

# Add the cotangents of the rotors from the vectors in the columns of `Ḡ` into the shadow
# `dc`'s copy of them.
function add_rotor_cotangents!(dc, left::Bool, c, Ḡ)
    for i ∈ eachindex(c.rotors)
        R̄ = rotor_cotangent(left, c.rotors[i], (Ḡ[1, i], Ḡ[2, i], Ḡ[3, i]))
        q = dc.rotors[i]
        dc.rotors[i] = Quaternion(q[1] + R̄[1], q[2] + R̄[2], q[3] + R̄[3], q[4] + R̄[4])
    end
    nothing
end

# The shadow of the returned calculator, as the configuration asks for it.
shadow_return(config, c) = needs_shadow(config) ? (width(config) == 1 ? only(shadows(c)) : shadows(c)) : nothing


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
            ntuple(_ -> nothing, length(args))
        end
    end
end


## 𝔇

# The block of degree ℓ: the stored rows of `c.Wˡ`, its columns, the rows that the block
# returns, and their offset among those stored.
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

function EnzymeRules.forward(
    config::FwdConfig, ::Const{typeof(compute_block!)}, ::Type{<:Annotation},
    c::Annotation{<:DCalc}, ℓ::Annotation
)
    calc, l = c.val, ℓ.val
    compute_block!(calc, l)
    # The shadow block is written whether or not the returned calculator's shadow is asked
    # for, since the code after the step reads the block from the shadow's buffer.
    left = derivatives_from_left(calc)
    rows, cols, outrows, outcols = block_geometry(calc, l)
    A = view(calc.Wˡ, :, 1:length(rows), 1:length(cols))
    for dc ∈ something(shadows(c), ())
        wigner_block_pushforward!(
            (x, ẋ) -> only(ẋ), labelled_view(dc.Wˡ, rows, cols, outrows, outcols), A, l,
            rows, cols, outrows, outcols, left, generators(left, calc, dc), Val(1)
        )
    end
    forward_return(config, calc, shadows(c))
end

function EnzymeRules.augmented_primal(
    config::RevConfig, ::Const{typeof(compute_block!)}, ::Type{<:Annotation},
    c::Annotation{<:DCalc}, ℓ::Annotation
)
    calc, l = c.val, ℓ.val
    compute_block!(calc, l)
    rows, cols, _, _ = block_geometry(calc, l)
    tape = calc.Wˡ[:, 1:length(rows), 1:length(cols)]
    AugmentedReturn(needs_primal(config) ? calc : nothing, shadow_return(config, c), tape)
end

function EnzymeRules.reverse(
    ::RevConfig, ::Const{typeof(compute_block!)}, ::Type{<:Annotation}, tape,
    c::Annotation{<:DCalc}, ℓ::Annotation
)
    calc, l = c.val, ℓ.val
    left = derivatives_from_left(calc)
    rows, cols, outrows, outcols = block_geometry(calc, l)
    for dc ∈ something(shadows(c), ())
        Ā = labelled_view(dc.Wˡ, rows, cols, outrows, outcols)
        Ḡ = zeros(real(eltype(tape)), 3, length(calc.rotors))
        wigner_block_pullback!(Ḡ, tape, Ā, l, rows, cols, outrows, outcols, left)
        add_rotor_cotangents!(dc, left, calc, Ḡ)
        fill!(Ā, zero(eltype(Ā)))
    end
    (nothing, nothing)
end


## The harmonics

# The part of the destination `Y` that holds the block of degree ℓ; see `compute_block!`.
block_view(Y, ℓ, j₀) = view(Y, :, :, (j₀ + 1):(j₀ + Int(2ℓ) + 1))

function EnzymeRules.forward(
    config::FwdConfig, ::Const{typeof(compute_block!)}, ::Type{<:Annotation},
    c::Annotation{<:YCalc}, ℓ::Annotation, is::Annotation, Y::Annotation, j₀::Annotation
)
    calc, l, i, j = c.val, ℓ.val, is.val, j₀.val
    compute_block!(calc, l, i, Y.val, j)
    # As for 𝔇, the shadow block is written whether or not the returned calculator's shadow
    # is asked for.
    if shadows(c) !== nothing && shadows(Y) !== nothing
        for (dc, dY) ∈ zip(shadows(c), shadows(Y))
            harmonic_block_pushforward!(
                (y, ẏ) -> only(ẏ), block_view(dY, l, j), block_view(Y.val, l, j), l, i,
                generators(true, calc, dc), Val(1)
            )
        end
    end
    forward_return(config, calc, shadows(c))
end

function EnzymeRules.augmented_primal(
    config::RevConfig, ::Const{typeof(compute_block!)}, ::Type{<:Annotation},
    c::Annotation{<:YCalc}, ℓ::Annotation, is::Annotation, Y::Annotation, j₀::Annotation
)
    calc, l, i, j = c.val, ℓ.val, is.val, j₀.val
    compute_block!(calc, l, i, Y.val, j)
    tape = copy(block_view(Y.val, l, j))
    AugmentedReturn(needs_primal(config) ? calc : nothing, shadow_return(config, c), tape)
end

function EnzymeRules.reverse(
    ::RevConfig, ::Const{typeof(compute_block!)}, ::Type{<:Annotation}, tape,
    c::Annotation{<:YCalc}, ℓ::Annotation, is::Annotation, Y::Annotation, j₀::Annotation
)
    calc, l, i, j = c.val, ℓ.val, is.val, j₀.val
    if shadows(c) !== nothing && shadows(Y) !== nothing
        for (dc, dY) ∈ zip(shadows(c), shadows(Y))
            Ȳ = block_view(dY, l, j)
            Ḡ = zeros(real(eltype(tape)), 3, length(calc.rotors))
            harmonic_block_pullback!(Ḡ, tape, Ȳ, l, i)
            add_rotor_cotangents!(dc, true, calc, Ḡ)
            fill!(view(Ȳ, :, i, :), zero(eltype(Ȳ)))
        end
    end
    (nothing, nothing, nothing, nothing, nothing)
end

end # module SphericalFunctionsEnzymeCoreExt
