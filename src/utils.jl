using Base: @propagate_inbounds


mutable struct OffsetVec{T,VV<:AbstractVector{T}} <: AbstractVector{T}
    parent::VV
    offset::Int
end
offset(V::OffsetVec) = V.offset
Base.parent(V::OffsetVec) = V.parent
@propagate_inbounds Base.getindex(V::OffsetVec, i::Int) = parent(V)[i-offset(V)]
@propagate_inbounds function Base.setindex!(V::OffsetVec, val, i::Int)
    parent(V)[i-offset(V)] = val
    V
end


mutable struct OffsetMat{T,MM<:AbstractMatrix{T}} <: AbstractMatrix{T}
    parent::MM
    offset1::Int
    offset2::Int
end
offset1(M::OffsetMat) = M.offset1
offset2(M::OffsetMat) = M.offset2
Base.parent(M::OffsetMat) = M.parent
@propagate_inbounds Base.getindex(M::OffsetMat, i::Int, j::Int) = parent(M)[i-offset1(M), j-offset2(M)]
#@propagate_inbounds Base.getindex(M::OffsetMat, ::Colon, j::Int) = parent(M)[:, j-offset2(M)]


loggamma(a, T::Type{<:DoubleFloats.MultipartFloat}) = DoubleFloats.loggamma(T(a))
loggamma(a, T::Type) = SpecialFunctions.loggamma(T(a))
function logbinomial(n::T, k::T, S=float(T)) where {T<:Integer}
    if k == 0 || k == n
        return zero(S)
    end
    if k > (n>>1)
        k = n - k
    end
    if k == 1
        return log(S(n))
    else
        return -log1p(S(n)) - loggamma(n - k + one(T), S) - loggamma(k + one(T), S) + loggamma(n + 2one(T), S)
    end
end

"""
    sqrtbinomial(n, k, [T])

Evaluate the square-root of the binomial coefficient `binomial(n,k)` for large coefficients.

Ordinarily, when `n` and `k` are standard `Int` arguments, the built-in `binomial` function
will overflow around `n=66`, because it results in `Int`s.  We need much larger values.
This function, which is based on [`a related one in
SpecialFunctions.jl`](https://specialfunctions.juliamath.org/latest/functions_list/#SpecialFunctions.logabsbinomial),
returns reasonably accurate results up to `n ≈ 1026` when `k ≈ n/2` (which is the case of
interest in many applications in this package).

Computations are carried out (and returned) in type `T`, which defaults to `Float64`.
"""
function sqrtbinomial(n, k, T=Float64)
    exp(logbinomial(n, k, T)/2)
end
