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
module NaNChecker

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


module ExplicitWignerMatrices

function d_explicit(n, m′, m, expiβ::Complex{T}) where T
    if abs(m′) < abs(m)
        return (-1)^(m-m′) * d_explicit(n, m, m′, expiβ)
    end
    if m′ < 0
        return (-1)^(m-m′) * d_explicit(n, -m′, -m, expiβ)
    end
    cosβ = expiβ.re
    sinβ = expiβ.im
    if (n,m′,m) == (0,0,0)
        T(1)
    elseif (n,m′,m) == (1,0,0)
        cosβ
    elseif (n,m′,m) == (1,1,-1)
        (1-cosβ) / 2
    elseif (n,m′,m) == (1,1,0)
        -sinβ / √T(2)
    elseif (n,m′,m) == (1,1,1)
        (1+cosβ) / 2
    elseif (n,m′,m) == (2,0,0)
        (3cosβ^2-1) / 2
    elseif (n,m′,m) == (2,1,-1)
        (1+cosβ-2cosβ^2) / 2
    elseif (n,m′,m) == (2,1,0)
        -√(T(3)/8) * 2 * sinβ * cosβ
    elseif (n,m′,m) == (2,1,1)
        (-1+cosβ+2cosβ^2) / 2
    elseif (n,m′,m) == (2,2,-2)
        (1-cosβ)^2 / 4
    elseif (n,m′,m) == (2,2,-1)
        -sinβ * (1-cosβ) / 2
    elseif (n,m′,m) == (2,2,0)
        √(T(3)/8) * sinβ^2
    elseif (n,m′,m) == (2,2,1)
        -sinβ * (1+cosβ) / 2
    elseif (n,m′,m) == (2,2,2)
        (1+cosβ)^2/4
    else
        T(NaN)
    end
end

function D_explicit(n, m′, m, expiα::Complex{T}, expiβ::Complex{T}, expiγ::Complex{T}) where T
    return expiα^(-m′) * d_explicit(n, m′, m, expiβ) * expiγ^(-m)
end


function d_formula(n, m′, m, expiβ::Complex{T}) where T
    # https://en.wikipedia.org/wiki/Wigner_D-matrix#Wigner_.28small.29_d-matrix
    cosβ = expiβ.re
    sin½β = √((1-cosβ)/2)
    cos½β = √((1+cosβ)/2)
    prefactor =√T(
        factorial(big(n + m′))
        * factorial(big(n - m′))
        * factorial(big(n + m))
        * factorial(big(n - m))
    )
    prefactor * sum(
        ifelse(iseven(m′ - m + s), 1, -1)
        * cos½β ^ (2n + m - m′ - 2s)
        * sin½β ^ (m′ - m + 2s)
        / (
            factorial(big(n + m - s))
            * factorial(big(s))
            * factorial(big(m′ - m + s))
            * factorial(big(n - m′ - s))
        )
        for s in max(0, m - m′):min(n + m, n - m′)
    )
end

function D_formula(n, m′, m, expiα::Complex{T}, expiβ::Complex{T}, expiγ::Complex{T}) where T
    # https://en.wikipedia.org/wiki/Wigner_D-matrix#Definition_of_the_Wigner_D-matrix
    # Note that the convention in this package is conjugated relative to the convention
    # used by Wikipedia, so we include that conjugation here.
    return expiα^(m′) * d_formula(n, m′, m, expiβ) * expiγ^(m)
end

end  # module ExplicitWignerMatrices


module NINJA

function Wigner_d(ι::T, ell, m, s) where {T<:Real}
    # Eq. II.8 of Ajith et al. (2007) 'Data formats...'
    k_min = max(0, m - s)
    k_max = min(ell + m, ell - s)
    prefactor = √T(
        factorial(big(ell + m))
        * factorial(big(ell - m))
        * factorial(big(ell + s))
        * factorial(big(ell - s))
    )
    sum(
        ifelse(iseven(k), 1, -1)
         * cos(ι / 2) ^ (2 * ell + m - s - 2 * k)
         * sin(ι / 2) ^ (2 * k + s - m)
         * prefactor
         / T(
             factorial(big(ell + m - k))
             * factorial(big(ell - s - k))
             * factorial(big(k))
             * factorial(big(k + s - m))
         )
        for k in k_min:k_max
    )
end

function sYlm(s, ell, m, ι::T, ϕ::T) where {T<:Real}
    # Eq. II.7 of Ajith et al. (2007) 'Data formats...'
    # Note the weird definition w.r.t. `-s`
    if abs(s) > ell || abs(m) > ell
        return zero(complex(T))
    end
    (
        ifelse(iseven(s), 1, -1)
        * √((2ell + 1) / (4T(π)))
        * Wigner_d(ι, ell, m, -s)
        * cis(m * ϕ)
    )
end

# Eqs. (II.9) through (II.13) of https://arxiv.org/abs/0709.0093v3 [Ajith_2007](@cite)
m2Y22(ι::T, ϕ::T) where {T<:Real} = √(5 / (64T(π))) * (1 + cos(ι))^2 * cis(2ϕ)
m2Y21(ι::T, ϕ::T) where {T<:Real} = √(5 / (16T(π))) * sin(ι) * (1 + cos(ι)) * cis(ϕ)
m2Y20(ι::T, ϕ::T) where {T<:Real} = √(15 / (32T(π))) * sin(ι)^2
m2Y2m1(ι::T, ϕ::T) where {T<:Real} = √(5 / (16T(π))) * sin(ι) * (1 - cos(ι)) * cis(-1ϕ)
m2Y2m2(ι::T, ϕ::T) where {T<:Real} = √(5 / (64T(π))) * (1 - cos(ι))^2 * cis(-2ϕ)

end  # module NINJA


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
