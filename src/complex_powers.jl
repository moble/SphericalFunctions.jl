"""
    complex_powers!(zpowers, z)

Compute integer powers of `z` from `z^0` through `z^m`, recursively, where `m` is
one less than the length of the input `zpowers` vector.

Note that `z` is assumed to be normalized, with complex amplitude approximately 1.

See also: [`complex_powers`](@ref)
"""
function complex_powers!(zpowers, z)
    Base.require_one_based_indexing(zpowers)
    @fastmath @inbounds begin
        M = length(zpowers)
        if M == 0
            return zpowers
        end
        M -= 1
        θ = one(z)
        zpowers[1] = θ
        if M == 0
            return zpowers
        end
        if M == 1
            zpowers[2] = z
            return zpowers
        end
        while z.re<0 || z.im<0
            θ *= 1im
            z /= 1im
        end
        zpowers[2] = z
        clock = θ
        dc = -2 * sqrt(z).im^2
        t = 2 * dc
        dz = dc * (1 + 2 * z) + 1im * sqrt(-dc * (2 + dc))
        for m in 2:M
            zpowers[m+1] = zpowers[m] + dz
            zpowers[m] *= clock
            clock *= θ
            dz += t * zpowers[m+1]
        end
        zpowers[M+1] *= clock
    end
    zpowers
end


"""
    complex_powers(z, m)

Compute integer powers of `z` from `z^0` through `z^m`, recursively.

Note that `z` is assumed to be normalized, with complex amplitude approximately 1.

This algorithm is mostly due to Stoer and Bulirsch in "Introduction to
Numerical Analysis" (page 24) — with a little help from de Moivre's formula,
which is essentially exp(iθ)ⁿ = exp(inθ), as well as my own alterations to deal
with different behaviors in different quadrants.

There isn't usually a huge advantage to using this specialized function.  If
you just need a particular power, it will generally be far more efficient and
just as accurate to compute either exp(iθ)ⁿ or exp(inθ) explicitly.  However,
if you need all powers from 0 to m, this function is about 10 or 5 times faster
than those options, respectively, for large m.  Like those options, this function
is numerically stable, in the sense that its errors are usually smaller than `m`
times the error from machine-precision errors in the input argument — or at worst
about 50% larger, which occurs as the phase approaches multiples of π/2.

See also: [`complex_powers!`](@ref)
"""
function complex_powers(z, m::Int)
    if abs(z) ≉ one(z)
        throw(DomainError("z = $z",
            "This function is only valid for `z` with complex amplitude approximately 1; "
            * "abs(z) = $(abs(z))"
        ))
    end
    zpowers = zeros(typeof(z), m+1)
    complex_powers!(zpowers, z)
end


struct ComplexPowers{T}
    z¹::T
    ϕ::Complex{Int}
    τ::T
    function ComplexPowers{T}(z::T, factor_phase=true) where {T}
        z¹ = z
        ϕ = 1 + 0im
        if factor_phase
            for _ ∈ 1:3  # We need at most 3 iterations to get the phase right
                if z¹.re < 0 || abs(z¹.im) > z¹.re
                    ϕ *= im
                    z¹ *= -im
                end
            end
        end
        τ = -4 * sqrt(z¹).im^2
        new(z¹, ϕ, τ)
    end
end

@doc raw"""
    ComplexPowers(z)

Construct an iterator to compute powers of the complex phase factor ``z`` (assumed to have
magnitude 1).  The iterator will return the complex number ``zᵐ`` for each integer ``m = 0,
1, 2, \ldots``.

# Example
```julia-repl
julia> cp = ComplexPowers(cis(0.1));

julia> first(cp, 5)  # Get the first 5 values from the iterator
5-element Vector{ComplexF64}:
                1.0 + 0.0im
 0.9950041652780258 + 0.09983341664682815im
 0.9800665778412417 + 0.19866933079506122im
 0.9553364891256061 + 0.2955202066613396im
 0.9210609940028851 + 0.3894183423086505im

julia> cis(0.1).^(0:4)
5-element Vector{ComplexF64}:
                1.0 + 0.0im
 0.9950041652780258 + 0.09983341664682815im
 0.9800665778412417 + 0.19866933079506124im
 0.9553364891256062 + 0.2955202066613396im
 0.9210609940028853 + 0.3894183423086506im
```

# Notes

[StoerBulirsch_2002](@citet) described the basic algorithm on page 24 (Example 4), though
there is a dramatic improvement to be made.  The basic idea is a recurrence relation, where
``zᵐ`` is computed from ``zᵐ⁻¹`` by adding a small increment ``δz``, which itself is updated
by adding a small increment given by ``zᵐ`` times a constant ``τ``.

Although this algorithm is numerically stable, we can improve its accuracy by factoring out
``ϕ``, the smallest power of ``i`` that minimizes the phase of ``z``.  This power of ``i``
can be separately exponentiated exactly and efficiently because it is exactly representable
as a complex integer, while the error in the computation of ``zᵐ`` is reduced significantly
for certain values of ``z``.

As implemented here, this algorithm achieves a worst-case accuracy of roughly ``m ϵ`` for
``zᵐ`` — where ``ϵ`` is the precision of the type of the input argument — across the range
of inputs ``z = \exp^{iθ}`` for ``θ ∈ [0, 2π]``.  The original algorithm can be far worse
for values of ``θ`` close to ``π`` — often orders of magnitude worse.

"""
ComplexPowers(cisθ::T, factor_phase=true) where {T} = ComplexPowers{T}(cisθ, factor_phase)


function Base.iterate(cp::ComplexPowers{T}) where {T}
    z⁰ = one(T)
    δc = cp.τ / 2
    δz = δc + im * √(-δc * (2 + δc)) * (cp.z¹.im ≥ 0 ? 1 : -1)
    Φ = 1 + 0im
    z⁰, (z⁰, δz, Φ)
end

function Base.iterate(cp::ComplexPowers{T}, state) where {T}
    (zᵐ⁻¹, δz, Φ) = state
    zᵐ = zᵐ⁻¹ + δz
    δz = δz + zᵐ * cp.τ
    Φ = Φ * cp.ϕ
    zᵐ*Φ, (zᵐ, δz, Φ)
end

Base.IteratorSize(::Type{<:ComplexPowers}) = Base.IsInfinite()

Base.eltype(::Type{ComplexPowers{T}}) where {T} = T

Base.isdone(iterator::ComplexPowers{T}) where {T} = false
Base.isdone(iterator::ComplexPowers{T}, state) where {T} = false


@testitem "ComplexPowers" begin
    mₘₐₓ = 10_000
    for θ ∈ BigFloat(0):big(1//10):2big(π)
        z¹ = cis(θ)
        p = ComplexPowers(ComplexF64(z¹))
        for (i, zᵐ) in enumerate(p)
            m = i-1
            err = Float64(abs(zᵐ - z¹^m))
            @test err < mₘₐₓ * eps(Float64)
            if m == mₘₐₓ
                break
            end
        end
    end
end
