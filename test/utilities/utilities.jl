@testsnippet Utilities begin

ℓmrange(ℓₘᵢₙ, ℓₘₐₓ) = eachrow(SphericalFunctions.Yrange(ℓₘᵢₙ, ℓₘₐₓ))
ℓmrange(ℓₘₐₓ) = ℓmrange(0, ℓₘₐₓ)
function sℓmrange(ℓₘₐₓ, sₘₐₓ)
    sₘₐₓ = min(abs(sₘₐₓ), ℓₘₐₓ)
    [
        (s, ℓ, m)
        for s in -sₘₐₓ:sₘₐₓ
        for ℓ in abs(s):ℓₘₐₓ
        for m in -ℓ:ℓ
    ]
end
function ℓm′mrange(ℓₘₐₓ)
    [
        (ℓ, m′, m)
        for ℓ in 0:ℓₘₐₓ
        for m′ in -ℓ:ℓ
        for m in -ℓ:ℓ
    ]
end

αrange(::Type{T}, n=15) where T = T[
    0; nextfloat(T(0)); rand(T(0):eps(T(π)):T(π), n÷2); prevfloat(T(π)); T(π);
    nextfloat(T(π)); rand(T(π):eps(2T(π)):2T(π), n÷2); prevfloat(T(π)); 2T(π)
]
βrange(::Type{T}=Float64, n=15; avoid_poles=0) where T = T[
    avoid_poles; nextfloat(T(avoid_poles));
    rand(T(0):eps(T(π)):T(π), n);
    prevfloat(T(π)-avoid_poles); T(π)-avoid_poles
]
γrange(::Type{T}, n=15) where T = αrange(T, n)
αβγrange(::Type{T}=Float64, n=15; avoid_poles=0) where T = vec(collect(
    Iterators.product(αrange(T, n), βrange(T, n; avoid_poles), γrange(T, n))
))

const θrange = βrange
const φrange = αrange
θϕrange(::Type{T}=Float64, n=15; avoid_poles=0) where T = vec(collect(
    Iterators.product(θrange(T, n; avoid_poles), φrange(T, n))
))

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
    # Eqs. (II.7) and (II.8) of https://arxiv.org/abs/0709.0093v3 [AjithEtAl_2011](@cite)
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
