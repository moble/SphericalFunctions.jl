# Explicit versions of the left and right angular-momentum operators.
#
# These implement the settled definitions from the conventions summary,
#
#     L_𝐮 f(𝐑) =  i d/dϵ f(e^{-ϵ𝐮/2} 𝐑),        R_𝐮 f(𝐑) = -i d/dϵ f(𝐑 e^{-ϵ𝐮/2}),
#
# by automatic differentiation, for testing purposes.  Writing e^{-ϵ𝐮/2} = e^{θ𝐮} with θ =
# -ϵ/2 gives L_𝐮 f = -(i/2) d/dθ f(e^{θ𝐮} 𝐑) and R_𝐮 f = +(i/2) d/dθ f(𝐑 e^{θ𝐮}),
# which is what is coded below.
#
# The derivative of `exp(θ*g)` at θ = 0 is taken through Quaternionic's `exp(::QuatVec)`,
# which expands its small-argument branch as a series rather than returning a constant
# there, so that ForwardDiff sees the true derivative.  `exp` is used rather than `cos(θ) +
# sin(θ)*g` because it returns a `Rotor`, whereas the other form is a `Quaternion` that is
# not of unit magnitude unless `g` is, and the package's harmonics are defined only on
# quaternions that denote rotations.
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
