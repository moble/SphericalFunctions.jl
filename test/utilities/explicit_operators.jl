"""Explicit versions of the raising and lowering operators.

These are just simple versions of the operators defined in
notes/operators/explicit_definition.jl, for testing purposes.  Note that we explicitly use
`cos(θ) + sin(θ)*g` instead of simply `exp(θ*g)`, because the `exp` implementation currently
has a special case at zero, which messes with the derivative at that point.  But also note
that these are incorrect for `g=0` because we oversimplify.

"""
@testmodule ExplicitOperators begin
    using Quaternionic
    import ForwardDiff

    function L(g::QuatVec{T}, f) where T
        function L_g(Q)
            -im * ForwardDiff.derivative(θ -> f((cos(θ) + sin(θ)*g) * Q), zero(T)) / 2
        end
    end
    function R(g::QuatVec{T}, f) where T
        function R_g(Q)
            -im * ForwardDiff.derivative(θ -> f(Q * (cos(θ) + sin(θ)*g)), zero(T)) / 2
        end
    end

end  # @testmodule ExplicitOperators
