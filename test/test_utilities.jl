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

function d(n, m′, m, expiβ::Complex{T}) where T
    if abs(m′) < abs(m)
        return (-1)^(m-m′) * d(n, m, m′, expiβ)
    end
    if m′ < 0
        return (-1)^(m-m′) * d(n, -m′, -m, expiβ)
    end
    cosβ = expiβ.re
    sinβ = expiβ.im
    if (n,m′,m) == (0,0,0)
        1
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
    end
end

function D(n, m′, m, expiα::Complex{T}, expiβ::Complex{T}, expiγ::Complex{T}) where T
    return expiα^(-m′) * d(n, m′, m, expiβ) * expiγ^(-m)
end

end
