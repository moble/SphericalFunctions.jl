module SphericalFunctionsForwardDiffExt

# ForwardDiff has no rule system of its own: a rule is a method for dual numbers.  Here the
# methods are those that tell the calculators how to read and build dual numbers (see
# `src/derivatives/lifting.jl`).  A calculator whose rotor data are dual numbers then runs
# the recurrence on the values of those data, and gives each block its partial derivatives
# from the angular-momentum operators, as described in `src/derivatives/kernels.jl`; the
# recurrence is never differentiated.  `D`, `d`, `sYlm`, and `sYlm_matrix` are computed by
# calculators, so they are differentiated the same way.
#
# The values are computed by a calculator of the dual numbers' values, which are themselves
# dual numbers under nested differentiation, as for a Hessian; that calculator lifts the
# blocks of another in turn, so that every level of derivatives is exact.  The rotor's
# partials may point in any direction, including out of the unit sphere, and the values are
# differentiated as functions of R/‖R‖, which is how they are computed.

import SphericalFunctions
using Quaternionic: Quaternion
using ForwardDiff: Dual, Partials, value, partials

SphericalFunctions.value_type(::Type{Dual{T, V, N}}) where {T, V, N} = V
SphericalFunctions.real_value(x::Dual) = value(x)
SphericalFunctions.ndirections(::Type{Dual{T, V, N}}) where {T, V, N} = N

function SphericalFunctions.rotor_tangents(q::Quaternion{Dual{T, V, N}}) where {T, V, N}
    ntuple(k -> (partials(q[1], k), partials(q[2], k), partials(q[3], k), partials(q[4], k)), Val(N))
end
SphericalFunctions.angle_tangents(β::Dual{T, V, N}) where {T, V, N} =
    ntuple(k -> partials(β, k), Val(N))

# A dual number with tag `T`, complex or real, from its value `x` and the tuple `ẋ` of its
# partial derivatives.
struct DualCombine{T} end
(::DualCombine{T})(x::Complex, ẋ) where {T} = Complex(
    Dual{T}(real(x), Partials(map(real, ẋ))), Dual{T}(imag(x), Partials(map(imag, ẋ)))
)
(::DualCombine{T})(x::Real, ẋ) where {T} = Dual{T}(x, Partials(ẋ))
SphericalFunctions.lift_combine(::Type{Dual{T, V, N}}) where {T, V, N} = DualCombine{T}()

end # module SphericalFunctionsForwardDiffExt
