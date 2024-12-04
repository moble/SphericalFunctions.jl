# NOTE: This file is for utilities used in tests — not for tests of utilities

"""
This module is copied from [this post by Brian Guenter](
https://discourse.julialang.org/t/treating-nan-as-error-helping-debugging/36933/9).
Essentially, you can wrap your `Real` type in `NaNCheck`, and pass it (or an
array) of such objects to a function you want to check to track down where
NaNs are coming into a calculation.  As soon as such an object is used in a
calculation, it raises an exception — essentially making this a
[signaling NaN](https://en.wikipedia.org/wiki/NaN#Signaling_NaN).

This can also be helpful for ensuring that no un-initialized memory is used in
a computation (which might silently give bad results).  Instead of creating an
`undef` array, just create the array from `NaNCheck`s.  The elements can be
set to some non-NaN value before use, and the calculation will proceed as
usual.  If they are not set to non-NaN, they will raise an exception whenever
they are used.

Such an array can be created for type `T` with something like

    arr = fill(NaNCheck{T}(NaN), arr_size)

"""
@testmodule NaNChecker begin

struct NaNCheck{T<:Real} <: Real
    val::T
    function NaNCheck{T}(a::S) where {T<:Real, S<:Real}
        @assert !(T <: NaNCheck)
        new{T}(T(a))
    end
end
export NaNCheck
Base.isnan(a::NaNCheck{T}) where{T} = isnan(a.val)
Base.isinf(a::NaNCheck{T}) where{T} = isinf(a.val)
Base.typemin(::Type{NaNCheck{T}}) where{T} = NaNCheck{T}(typemin(T))
Base.typemax(::Type{NaNCheck{T}}) where{T} = NaNCheck{T}(typemax(T))
Base.eps(::Type{NaNCheck{T}}) where {T} = NaNCheck{T}(eps(T))
Base.decompose(a::NaNCheck{T}) where {T} = Base.decompose(a.val)
Base.round(a::NaNCheck{T}, m::RoundingMode) where {T} = NaNCheck{T}(round(a.val, m))

struct NaNException <: Exception end

# (::Type{Float64})(a::NaNCheck{S}) where {S<:Real} = NaNCheck{Float64}(Float64(a.val))
(::Type{T})(a::NaNCheck{S}) where {T<:Integer,S<:Real} = T(a.val)
(::Type{NaNCheck{T}})(a::NaNCheck{S}) where {T<:Real,S<:Real} = NaNCheck{T}(T(a.val))
Base.promote_rule(::Type{NaNCheck{T}}, ::Type{T}) where {T<:Number} = NaNCheck{T}
Base.promote_rule(::Type{T}, ::Type{NaNCheck{T}}) where {T<:Number} = NaNCheck{T}
Base.promote_rule(::Type{S}, ::Type{NaNCheck{T}}) where {T<:Number, S<:Number} = NaNCheck{promote_type(T,S)}
Base.promote_rule(::Type{NaNCheck{T}}, ::Type{S}) where {T<:Number, S<:Number} = NaNCheck{promote_type(T,S)}
Base.promote_rule(::Type{NaNCheck{S}}, ::Type{NaNCheck{T}}) where {T<:Number, S<:Number} = NaNCheck{promote_type(T,S)}

for op = (:sin, :cos, :tan, :log, :exp, :sqrt, :abs, :-, :atan, :acos, :asin, :log1p, :floor, :ceil, :float)
    eval(quote
        function Base.$op(a::NaNCheck{T}) where{T}
            temp = NaNCheck{T}(Base.$op(a.val))
            if isnan(temp)
                throw(NaNException())
            end
            return temp
        end
    end)
end

for op = (:+, :-, :/, :*, :^, :atan)
    eval(quote
        function Base.$op(a::NaNCheck{T}, b::NaNCheck{T}) where{T}
            temp = NaNCheck{T}(Base.$op(a.val, b.val))
            if isnan(temp)
                throw(NaNException())
            end
            return temp
        end
    end)
end

for op =  (:<, :>, :<=, :>=, :(==), :isequal)
    eval(quote
        function Base.$op(a::NaNCheck{T}, b::NaNCheck{T}) where{T}
            temp = Base.$op(a.val, b.val)
            return temp
        end
    end)
end

end  # module NaNChecker
NaNCheck = NaNChecker.NaNCheck;



@testsnippet Utilities begin

αrange(::Type{T}, n=15) where T = T[
    0; nextfloat(T(0)); rand(T(0):eps(T(π)):T(π), n÷2); prevfloat(T(π)); T(π);
    nextfloat(T(π)); rand(T(π):eps(2T(π)):2T(π), n÷2); prevfloat(T(π)); 2T(π)
]
βrange(::Type{T}, n=15) where T = T[
    0; nextfloat(T(0)); rand(T(0):eps(T(π)):T(π), n); prevfloat(T(π)); T(π)
]
γrange(::Type{T}, n=15) where T = αrange(T, n)
v̂range(::Type{T}, n=15) where T = QuatVec{T}[
    𝐢; 𝐣; 𝐤;
    -𝐢; -𝐣; -𝐤;
    Quaternionic.normalize.(randn(QuatVec{T}, n))
]
function Rrange(::Type{T}, n=15) where T
    invsqrt2 = inv(√T(2))
    [
        [
            sign*R
            for R in [
                Rotor{T}(1);
                [Rotor{T}(𝐯) for 𝐯 in (𝐢,𝐣,𝐤)];
                [Rotor{T}(invsqrt2 + invsqrt2*𝐯) for 𝐯 in (𝐢,𝐣,𝐤)];
                [Rotor{T}(invsqrt2 - invsqrt2*𝐯) for 𝐯 in (𝐢,𝐣,𝐤)]
        ]
        for sign in (1,-1)
        ];
        randn(Rotor{T}, n)
    ]
end
epsilon(k) = ifelse(k>0 && isodd(k), -1, 1)

"""
    array_equal(a1, a2, equal_nan=false)

Ensure that arrays have same types and shapes, and all elements are the same.
If `equal_nan` is `true`, NaNs in the same place in each array will be
considered to be equal.

Note that this is slightly stricter than the numpy version of this function,
because arrays of different type will not be considered equal.

"""
function array_equal(a1::T1, a2::T2, equal_nan=false) where {T1, T2}
    if T1 !== T2 || size(a1) != size(a2)
        return false
    end
    all(e->e[1]==e[2] || (equal_nan && isnan(e1) && isnan(e2)), zip(a1, a2))
end

function sYlm(s::Int, ell::Int, m::Int, theta::T, phi::T) where {T<:Real}
    # Eqs. (II.7) and (II.8) of https://arxiv.org/abs/0709.0093v3 [Ajith_2007](@cite)
    # Note their weird definition w.r.t. `-s`
    k_min = max(0, m + s)
    k_max = min(ell + m, ell + s)
    sin_half_theta, cos_half_theta = sincos(theta / 2)
    return (-1)^(-s) * sqrt((2 * ell + 1) / (4 * T(π))) *
        T(sum(
            (-1) ^ (k)
            * sqrt(factorial(big(ell + m)) * factorial(big(ell - m)) * factorial(big(ell - s)) * factorial(big(ell + s)))
            * (cos_half_theta ^ (2 * ell + m + s - 2 * k))
            * (sin_half_theta ^ (2 * k - s - m))
            / (factorial(big(ell + m - k)) * factorial(big(ell + s - k)) * factorial(big(k)) * factorial(big(k - s - m)))
            for k in k_min:k_max
        )) *
        cis(m * phi)
end

ε(j,k,l) = ifelse(
    (j,k,l)∈((1,2,3),(2,3,1),(3,1,2)),
    1,
    ifelse(
        (j,k,l)∈((2,1,3),(1,3,2),(3,2,1)),
        -1,
        0
    )
)

end  # Utilities snippet


@testmodule ExplicitOperators begin
    using Quaternionic
    import ForwardDiff

    # These are just simple versions of the operators defined in
    # notes/operators/explicit_definition.jl, for testing purposes.  Note that we explicitly
    # use `cos(θ) + sin(θ)*g` instead of simply `exp(θ*g)`, because the `exp` implementation
    # currently has a special case at zero, which messes with the derivative at that point.
    # But also note that these are incorrect for `g=0` because we oversimplify.
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

end


"""Write factorials in the obvious way, for testing purposes.

This snippet lets us write things like `5❗` to get `120`, for example.  This should
probably only be used in tests, because factorials are not very efficient.  Therefore,
to ensure accuracy and avoid overflow, the argument is first converted to a `BigInt`,
making it even more inefficient.  This "big" conversion should be preserved as long as
possible to ensure accurate cancellations in the final result.
"""
@testsnippet NaiveFactorials begin
    struct Factorial end
    Base.:*(n::Integer, ::Factorial) = factorial(big(n))
    const ❗ = Factorial()
end
