module SphericalFunctionsChainRulesCoreExt

# `frule` and `rrule` for `D_array` and `sYlm_array`, from the generators, as described in
# `src/derivatives.jl`.  These serve the tools that read ChainRules directly, such as
# Zygote and Diffractor.
#
# A tangent of the rotor may arrive as any quaternion, as a structural `Tangent` of the
# `Rotor`, or as a vector of four components; it is read by its components, whatever its
# type.  The cotangent returned is a `Quaternion`, not a `Rotor`: a cotangent is not a unit
# quaternion, and a `Rotor` would be taken to be one by any method that reached it.

import SphericalFunctions: D_array, sYlm_array, D_array_widened, D_is_widened, D_narrowed,
    D_pushforward, D_pullback, sYlm_pushforward, sYlm_pullback, rotor_generator,
    rotor_cotangent, IntegerHalf
using Quaternionic: AbstractQuaternion, Rotor, Quaternion
import ChainRulesCore
using ChainRulesCore: AbstractZero, AbstractThunk, NoTangent, ZeroTangent, Tangent, unthunk,
    backing

# The four components of a tangent of a rotor, or `nothing` for a zero tangent.
tangent_components(::AbstractZero) = nothing
tangent_components(Ṙ::AbstractThunk) = tangent_components(unthunk(Ṙ))
tangent_components(Ṙ::AbstractQuaternion) = (Ṙ[1], Ṙ[2], Ṙ[3], Ṙ[4])
tangent_components(Ṙ::AbstractVector) = (Ṙ[1], Ṙ[2], Ṙ[3], Ṙ[4])
tangent_components(Ṙ::Tuple) = (Ṙ[1], Ṙ[2], Ṙ[3], Ṙ[4])
# A structural tangent with no fields is a zero tangent.
function tangent_components(Ṙ::Tangent{<:AbstractQuaternion})
    haskey(backing(Ṙ), :components) ? tangent_components(backing(Ṙ).components) : nothing
end
function tangent_components(Ṙ::Tangent)  # an `SVector`'s
    haskey(backing(Ṙ), :data) ? tangent_components(backing(Ṙ).data) : nothing
end

function ChainRulesCore.frule(
    (_, Ṙ, _, _, _), ::typeof(sYlm_array), R::Rotor, ℓₘₐₓ::IT, s, ℓₘᵢₙ::IT
) where {IT<:IntegerHalf}
    Y = sYlm_array(R, ℓₘₐₓ, s, ℓₘᵢₙ)
    ṙ = tangent_components(Ṙ)
    ṙ === nothing && return (Y, ZeroTangent())
    (Y, sYlm_pushforward(Y, rotor_generator(R, ṙ), ℓₘᵢₙ, ℓₘₐₓ))
end

function ChainRulesCore.rrule(
    ::typeof(sYlm_array), R::Rotor, ℓₘₐₓ::IT, s, ℓₘᵢₙ::IT
) where {IT<:IntegerHalf}
    Y = sYlm_array(R, ℓₘₐₓ, s, ℓₘᵢₙ)
    function sYlm_array_pullback(ΔY)
        Ȳ = unthunk(ΔY)
        R̄ = if Ȳ isa AbstractZero
            ZeroTangent()
        else
            Quaternion(rotor_cotangent(R, sYlm_pullback(Y, Ȳ, ℓₘᵢₙ, ℓₘₐₓ))...)
        end
        (NoTangent(), R̄, NoTangent(), NoTangent(), NoTangent())
    end
    (Y, sYlm_array_pullback)
end

function ChainRulesCore.frule(
    (_, Ṙ, _, _, _, _, _), ::typeof(D_array),
    R::Rotor, ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT
) where {IT<:IntegerHalf}
    Aʷ = D_array_widened(R, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    A = D_is_widened(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ) ? D_narrowed(Aʷ, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ) : Aʷ
    ṙ = tangent_components(Ṙ)
    ṙ === nothing && return (A, ZeroTangent())
    (A, D_pushforward(Aʷ, rotor_generator(R, ṙ), ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ))
end

function ChainRulesCore.rrule(
    ::typeof(D_array), R::Rotor, ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT
) where {IT<:IntegerHalf}
    Aʷ = D_array_widened(R, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    A = D_is_widened(ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ) ? D_narrowed(Aʷ, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ) : Aʷ
    function D_array_pullback(ΔA)
        Ā = unthunk(ΔA)
        R̄ = if Ā isa AbstractZero
            ZeroTangent()
        else
            g = D_pullback(Aʷ, Ā, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
            Quaternion(rotor_cotangent(R, g)...)
        end
        (NoTangent(), R̄, NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent())
    end
    (A, D_array_pullback)
end

end # module SphericalFunctionsChainRulesCoreExt
