module SphericalFunctionsForwardDiffExt

# Rules for ForwardDiff, which has no rule system of its own: a rule is a method for dual
# numbers.  `D_array` and `sYlm_array` of a rotor of duals evaluate the values at the rotor
# of their values, and combine those values into each partial derivative by the generators,
# as described in `src/derivatives.jl`.  The values are computed by calling the same
# function again, so a rotor of nested duals, as in a Hessian, comes back here once for
# every level of nesting, and each level of derivatives is exact.
#
# The recurrence is never differentiated here, so the rotor's partials may point in any
# direction, including out of the unit sphere, and the values are differentiated as
# functions of R/‖R‖, which is how they are computed.

import SphericalFunctions
import SphericalFunctions: D_array, sYlm_array, D_array_widened, D_pushforward,
    sYlm_pushforward, rotor_generator, IntegerHalf
using Quaternionic: Rotor
using ForwardDiff: Dual, Partials, value, partials

# The rotor of the values of the components, built without normalizing it, and the tangent
# given by the k-th partial of each component.
primal_rotor(R::Rotor{Dual{T, V, N}}) where {T, V, N} =
    Rotor{V}(value(R[1]), value(R[2]), value(R[3]), value(R[4]))
tangent_components(R::Rotor{<:Dual}, k) =
    (partials(R[1], k), partials(R[2], k), partials(R[3], k), partials(R[4], k))

# A complex dual number with tag `T`, from its value `x` and the tuple `ẋ` of its partial
# derivatives, as the `combine` argument of the pushforwards builds each element.
dual(::Type{T}, x::Complex, ẋ::Tuple) where {T} = Complex(
    Dual{T}(real(x), Partials(map(real, ẋ))), Dual{T}(imag(x), Partials(map(imag, ẋ)))
)

function SphericalFunctions.sYlm_array(
    R::Rotor{Dual{T, V, N}}, ℓₘₐₓ::IT, s, ℓₘᵢₙ::IT
) where {T, V, N, IT<:IntegerHalf}
    R₀ = primal_rotor(R)
    vs = ntuple(k -> rotor_generator(R₀, tangent_components(R, k)), Val(N))
    sYlm_pushforward((y, ẏ) -> dual(T, y, ẏ), sYlm_array(R₀, ℓₘₐₓ, s, ℓₘᵢₙ), vs, ℓₘᵢₙ, ℓₘₐₓ)
end

function SphericalFunctions.D_array(
    R::Rotor{Dual{T, V, N}}, ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT
) where {T, V, N, IT<:IntegerHalf}
    R₀ = primal_rotor(R)
    vs = ntuple(k -> rotor_generator(R₀, tangent_components(R, k)), Val(N))
    Aʷ = D_array_widened(R₀, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
    D_pushforward((x, ẋ) -> dual(T, x, ẋ), Aʷ, vs, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ)
end

end # module SphericalFunctionsForwardDiffExt
