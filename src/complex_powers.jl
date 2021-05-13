"""
    complex_powers!(zpowers, z)

Compute integer powers of `z` from `z^0` through `z^m`, recursively, where `m` is
one less than the length of the input `zpowers` vector.

Note that `z` is assumed to be normalized, with complex amplitude approximately 1.

See also: [`complex_powers`](@ref)
"""
function complex_powers!(zpowers::Vector{Complex{T}}, z::Complex{T}) where {T<:Real}
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
function complex_powers(z::Complex{T}, m::Int) where {T<:Real}
    if abs(z) ≉ one(z)
        throw(DomainError("z = $z",
            "This function is only valid for `z` with complex amplitude approximately 1; "
            * "abs(z) = $(abs(z))"
        ))
    end
    zpowers = Array{typeof(z)}(undef, m+1)
    complex_powers!(zpowers, z)
end

