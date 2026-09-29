module SphericalFunctionsEnzymeCoreExt

# Enzyme's rules for `D_array` and `sYlm_array`, from the generators, as described in
# `src/derivatives.jl`.  They are defined in EnzymeCore, which is all that a rule needs.
#
# Enzyme gives a shadow the type of its primal, so the tangent of a rotor arrives as a
# `Rotor`, and its cotangent must be returned as one, although neither has unit norm.  Both
# are read and built by their components alone.
#
# The index arguments are integers or ranges of them, which Enzyme treats as constants, but
# the rules accept any annotation of them rather than only `Const`: a rule whose signature
# does not match the call is not an error, but is silently passed over, and Enzyme would
# then differentiate the recurrence instead.

import SphericalFunctions: D_array, sYlm_array, D_array_widened, D_is_widened, D_narrowed,
    D_pushforward, D_pullback, sYlm_pushforward, sYlm_pullback, rotor_generator,
    rotor_cotangent
using Quaternionic: Rotor
using EnzymeCore: EnzymeRules, Annotation, Const, Active, Duplicated, BatchDuplicated
using EnzymeCore.EnzymeRules: FwdConfig, RevConfig, AugmentedReturn, needs_primal,
    needs_shadow, width

# The tangents of the rotor in each of the `width` directions, or `nothing` for a constant.
# Any other annotation is an error, rather than being taken to have no tangent.
rotor_tangents(::Const) = nothing
rotor_tangents(R::Duplicated) = (R.dval,)
rotor_tangents(R::BatchDuplicated) = R.dval

# What a forward rule returns, given the primal and a tuple of `width` shadows.
function forward_return(config::FwdConfig, Y, Ẏ::Tuple)
    if needs_primal(config)
        width(config) == 1 ? Duplicated(Y, only(Ẏ)) : BatchDuplicated(Y, Ẏ)
    else
        width(config) == 1 ? only(Ẏ) : Ẏ
    end
end

zero_shadows(Y, ::Val{1}) = zero(Y)
zero_shadows(Y, ::Val{N}) where {N} = ntuple(_ -> zero(Y), Val(N))

# The derivative of the rotor as a `Rotor` of its type, from the components of its cotangent,
# and the zero derivative, in each of the `width` directions, of an active rotor on which
# nothing that was differentiated depends.
as_rotor(::Rotor{T}, R̄) where {T} = Rotor{T}(R̄...)
zero_rotors(R::Rotor{T}, ::Val{1}) where {T} = as_rotor(R, ntuple(_ -> zero(T), 4))
zero_rotors(R::Rotor, ::Val{N}) where {N} = ntuple(_ -> zero_rotors(R, Val(1)), Val(N))


## The harmonics

function EnzymeRules.forward(
    config::FwdConfig, ::Const{typeof(sYlm_array)}, ::Type{<:Annotation},
    R::Annotation{<:Rotor}, ℓₘₐₓ::Annotation, s::Annotation, ℓₘᵢₙ::Annotation
)
    Y = sYlm_array(R.val, ℓₘₐₓ.val, s.val, ℓₘᵢₙ.val)
    needs_shadow(config) || return needs_primal(config) ? Y : nothing
    Ṙ = rotor_tangents(R)
    Ẏ = if Ṙ === nothing
        ntuple(_ -> zero(Y), Val(width(config)))
    else
        map(Ṙ) do Ṙₖ
            sYlm_pushforward(Y, rotor_generator(R.val, Ṙₖ), ℓₘᵢₙ.val, ℓₘₐₓ.val)
        end
    end
    forward_return(config, Y, Ẏ)
end

# The tape holds the values, copied if the caller is also given them and so might change
# them, and the shadows into which the cotangents are accumulated.
function EnzymeRules.augmented_primal(
    config::RevConfig, ::Const{typeof(sYlm_array)}, ::Type{<:Annotation},
    R::Annotation{<:Rotor}, ℓₘₐₓ::Annotation, s::Annotation, ℓₘᵢₙ::Annotation
)
    Y = sYlm_array(R.val, ℓₘₐₓ.val, s.val, ℓₘᵢₙ.val)
    Ȳ = needs_shadow(config) ? zero_shadows(Y, Val(width(config))) : nothing
    AugmentedReturn(needs_primal(config) ? Y : nothing, Ȳ, (needs_primal(config) ? copy(Y) : Y, Ȳ))
end

function EnzymeRules.reverse(
    config::RevConfig, ::Const{typeof(sYlm_array)}, ::Type{<:Annotation}, tape,
    R::Annotation{<:Rotor}, ℓₘₐₓ::Annotation, ::Annotation, ℓₘᵢₙ::Annotation
)
    Y, Ȳ = tape
    dR = if !(R isa Active)
        nothing
    elseif Ȳ === nothing
        zero_rotors(R.val, Val(width(config)))
    else
        cotangent(Ȳₖ) = as_rotor(R.val, rotor_cotangent(
            R.val, sYlm_pullback(Y, Ȳₖ, ℓₘᵢₙ.val, ℓₘₐₓ.val)
        ))
        width(config) == 1 ? cotangent(Ȳ) : map(cotangent, Ȳ)
    end
    (dR, nothing, nothing, nothing)
end


## Wigner's 𝔇

function EnzymeRules.forward(
    config::FwdConfig, ::Const{typeof(D_array)}, ::Type{<:Annotation},
    R::Annotation{<:Rotor}, ℓₘₐₓ::Annotation, m′ₘₐₓ::Annotation, m′ₘᵢₙ::Annotation,
    mₘₐₓ::Annotation, mₘᵢₙ::Annotation
)
    limits = (ℓₘₐₓ.val, m′ₘₐₓ.val, m′ₘᵢₙ.val, mₘₐₓ.val, mₘᵢₙ.val)
    needs_shadow(config) || return needs_primal(config) ? D_array(R.val, limits...) : nothing
    Aʷ = D_array_widened(R.val, limits...)
    A = D_is_widened(limits[1:3]...) ? D_narrowed(Aʷ, limits...) : Aʷ
    Ṙ = rotor_tangents(R)
    Ȧ = if Ṙ === nothing
        ntuple(_ -> zero(A), Val(width(config)))
    else
        map(Ṙₖ -> D_pushforward(Aʷ, rotor_generator(R.val, Ṙₖ), limits...), Ṙ)
    end
    forward_return(config, A, Ȧ)
end

function EnzymeRules.augmented_primal(
    config::RevConfig, ::Const{typeof(D_array)}, ::Type{<:Annotation},
    R::Annotation{<:Rotor}, ℓₘₐₓ::Annotation, m′ₘₐₓ::Annotation, m′ₘᵢₙ::Annotation,
    mₘₐₓ::Annotation, mₘᵢₙ::Annotation
)
    limits = (ℓₘₐₓ.val, m′ₘₐₓ.val, m′ₘᵢₙ.val, mₘₐₓ.val, mₘᵢₙ.val)
    Aʷ = D_array_widened(R.val, limits...)
    widened = D_is_widened(limits[1:3]...)
    A = widened ? D_narrowed(Aʷ, limits...) : Aʷ
    Ā = needs_shadow(config) ? zero_shadows(A, Val(width(config))) : nothing
    # `Aʷ` is `A` itself when nothing was widened, and must then be copied as `A` would be.
    kept = needs_primal(config) && !widened ? copy(Aʷ) : Aʷ
    AugmentedReturn(needs_primal(config) ? A : nothing, Ā, (kept, Ā))
end

function EnzymeRules.reverse(
    config::RevConfig, ::Const{typeof(D_array)}, ::Type{<:Annotation}, tape,
    R::Annotation{<:Rotor}, ℓₘₐₓ::Annotation, m′ₘₐₓ::Annotation, m′ₘᵢₙ::Annotation,
    mₘₐₓ::Annotation, mₘᵢₙ::Annotation
)
    limits = (ℓₘₐₓ.val, m′ₘₐₓ.val, m′ₘᵢₙ.val, mₘₐₓ.val, mₘᵢₙ.val)
    Aʷ, Ā = tape
    dR = if !(R isa Active)
        nothing
    elseif Ā === nothing
        zero_rotors(R.val, Val(width(config)))
    else
        cotangent(Āₖ) = as_rotor(R.val, rotor_cotangent(R.val, D_pullback(Aʷ, Āₖ, limits...)))
        width(config) == 1 ? cotangent(Ā) : map(cotangent, Ā)
    end
    (dR, nothing, nothing, nothing, nothing, nothing)
end

end # module SphericalFunctionsEnzymeCoreExt
