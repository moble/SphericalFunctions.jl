# Explicit versions of the left and right angular-momentum operators.
#
# These implement the settled definitions from the conventions summary,
#
#     L_𝐮 f(𝐑) =  i d/dϵ f(e^{-ϵ𝐮/2} 𝐑),        R_𝐮 f(𝐑) = -i d/dϵ f(𝐑 e^{-ϵ𝐮/2}),
#
# by automatic differentiation, for testing purposes.  Writing e^{-ϵ𝐮/2} = e^{θ𝐮} with θ = -ϵ/2
# gives L_𝐮 f = -(i/2) d/dθ f(e^{θ𝐮} 𝐑) and R_𝐮 f = +(i/2) d/dθ f(𝐑 e^{θ𝐮}), which is what is
# coded below.
#
# This used to spell `exp(θ*g)` out as `cos(θ) + sin(θ)*g`, because `exp` once had a special
# case at zero that returned a constant and so flattened the ForwardDiff derivative there.
# Quaternionic's `exp(::QuatVec)` now expands that branch as a series instead — its comment says
# "to obtain accurate ForwardDiff derivative" — so the workaround is no longer needed, and the
# two agree to 6e-17 where it was used.  Writing `exp` is worth the change beyond tidiness: it
# returns a `Rotor`, whereas `cos(θ) + sin(θ)*g` is a `Quaternion` that is not of unit magnitude
# unless `g` is, which made these tests the one thing in the package that needed the harmonics
# to accept a quaternion that does not denote a rotation on its face.
@testmodule ExplicitOperators begin
    using Quaternionic
    import ForwardDiff

    function L(g::QuatVec{T}, f) where T
        function L_g(Q)
            -im * ForwardDiff.derivative(θ -> f(exp(θ*g) * Q), zero(T)) / 2
        end
    end
    function R(g::QuatVec{T}, f) where T
        function R_g(Q)
            im * ForwardDiff.derivative(θ -> f(Q * exp(θ*g)), zero(T)) / 2
        end
    end

end  # @testmodule ExplicitOperators
