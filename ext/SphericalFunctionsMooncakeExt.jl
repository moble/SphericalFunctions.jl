module SphericalFunctionsMooncakeExt

# Mooncake's rules for `D_array` and `sYlm_array`, from the generators, as described in
# `src/derivatives.jl`, in forward and in reverse mode.  They are primitives of Mooncake's
# `MinimalCtx`, the context for rules that are needed for correctness rather than for speed,
# since without them Mooncake would differentiate the recurrence, whose derivatives are
# wrong near the poles.  Mooncake's `@from_chainrules` could instead adopt the ChainRules
# rules, but it converts tangents between the two systems generically, which its
# documentation recommends only for floats and arrays of them, and the rotor is neither.
#
# Mooncake's tangent of a `Rotor` is structural: the tuple of its four components, inside
# the tangent of the `SVector` that holds them, inside the tangent of the `Rotor`.  The
# tangent of an array of complex numbers is an array of complex numbers, with the
# convention for cotangents that the reverse rules assume.

import SphericalFunctions: D_array, sYlm_array, D_array_widened, D_is_widened, D_narrowed,
    D_pushforward, D_pullback, sYlm_pushforward, sYlm_pullback, rotor_generator,
    rotor_cotangent, IEEEFloat
using Quaternionic: Rotor
import Mooncake
using Mooncake: MinimalCtx, @is_primitive, CoDual, NoRData, RData, primal, tangent

# The four components of the tangent of a rotor, and the reverse data of a rotor from the
# four components of its cotangent.
rotor_tangent_components(Ṙ::Mooncake.Tangent) = Ṙ.fields.components.fields.data
rotor_rdata(R̄::NTuple{4}) = RData((components=RData((data=R̄,)),))

@is_primitive MinimalCtx Tuple{typeof(sYlm_array), Rotor{<:IEEEFloat}, Any, Any, Any}
@is_primitive MinimalCtx Tuple{typeof(D_array), Rotor{<:IEEEFloat}, Any, Any, Any, Any, Any}


## The harmonics

function Mooncake.frule!!(
    ::Mooncake.Dual{typeof(sYlm_array)}, R::Mooncake.Dual{<:Rotor},
    ℓₘₐₓ::Mooncake.Dual, s::Mooncake.Dual, ℓₘᵢₙ::Mooncake.Dual
)
    R₀ = primal(R)
    Y = sYlm_array(R₀, primal(ℓₘₐₓ), primal(s), primal(ℓₘᵢₙ))
    v = rotor_generator(R₀, rotor_tangent_components(tangent(R)))
    Mooncake.Dual(Y, sYlm_pushforward(Y, v, primal(ℓₘᵢₙ), primal(ℓₘₐₓ)))
end

function Mooncake.rrule!!(
    ::CoDual{typeof(sYlm_array)}, R::CoDual{<:Rotor},
    ℓₘₐₓ::CoDual, s::CoDual, ℓₘᵢₙ::CoDual
)
    R₀ = primal(R)
    Y = sYlm_array(R₀, primal(ℓₘₐₓ), primal(s), primal(ℓₘᵢₙ))
    # The caller is given `Y`, and may change it, so the values that the pullback reads are
    # kept apart from it.
    Y₀ = copy(Y)
    Ȳ = zero(Y)
    function sYlm_array_pullback(::NoRData)
        R̄ = rotor_cotangent(R₀, sYlm_pullback(Y₀, Ȳ, primal(ℓₘᵢₙ), primal(ℓₘₐₓ)))
        (NoRData(), rotor_rdata(R̄), NoRData(), NoRData(), NoRData())
    end
    (CoDual(Y, Ȳ), sYlm_array_pullback)
end


## Wigner's 𝔇

function Mooncake.frule!!(
    ::Mooncake.Dual{typeof(D_array)}, R::Mooncake.Dual{<:Rotor},
    ℓₘₐₓ::Mooncake.Dual, m′ₘₐₓ::Mooncake.Dual, m′ₘᵢₙ::Mooncake.Dual,
    mₘₐₓ::Mooncake.Dual, mₘᵢₙ::Mooncake.Dual
)
    R₀ = primal(R)
    limits = (primal(ℓₘₐₓ), primal(m′ₘₐₓ), primal(m′ₘᵢₙ), primal(mₘₐₓ), primal(mₘᵢₙ))
    Aʷ = D_array_widened(R₀, limits...)
    A = D_is_widened(limits[1:3]...) ? D_narrowed(Aʷ, limits...) : Aʷ
    v = rotor_generator(R₀, rotor_tangent_components(tangent(R)))
    Mooncake.Dual(A, D_pushforward(Aʷ, v, limits...))
end

function Mooncake.rrule!!(
    ::CoDual{typeof(D_array)}, R::CoDual{<:Rotor},
    ℓₘₐₓ::CoDual, m′ₘₐₓ::CoDual, m′ₘᵢₙ::CoDual, mₘₐₓ::CoDual, mₘᵢₙ::CoDual
)
    R₀ = primal(R)
    limits = (primal(ℓₘₐₓ), primal(m′ₘₐₓ), primal(m′ₘᵢₙ), primal(mₘₐₓ), primal(mₘᵢₙ))
    Aʷ = D_array_widened(R₀, limits...)
    widened = D_is_widened(limits[1:3]...)
    # As for the harmonics, the values that the pullback reads are kept apart from those the
    # caller is given.
    A = widened ? D_narrowed(Aʷ, limits...) : copy(Aʷ)
    Ā = zero(A)
    function D_array_pullback(::NoRData)
        R̄ = rotor_cotangent(R₀, D_pullback(Aʷ, Ā, limits...))
        (NoRData(), rotor_rdata(R̄), NoRData(), NoRData(), NoRData(), NoRData(), NoRData())
    end
    (CoDual(A, Ā), D_array_pullback)
end

end # module SphericalFunctionsMooncakeExt
