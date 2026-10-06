# The powers of a complex phase z, computed by the recurrence of Stoer and Bulirsch, in
# which each power is the previous one plus an increment, and each increment the previous
# one plus the new power times a real constant:
#
#     zᵐ⁺¹ = zᵐ + δzᵐ,    δzᵐ⁺¹ = δzᵐ + t zᵐ⁺¹,    t = 2(Re z - 1),
#
# which holds because δzᵐ = zᵐ (z - 1) and zᵐ⁺¹ (z - 1)² = zᵐ⁺¹ (z - 2 + z⁻¹) for |z| = 1.  The
# in-place function `complex_powers!` and the iterator `ComplexPowers` are two ways of
# running this one recurrence, and they share its start and its step, below, so that they
# compute the same values.

# The sum x² + y², rounded once.  `fma` rounds once on every machine, in hardware where
# there is an FMA instruction and in software where there is not; `muladd` only permits
# that, and on a machine where it is not fused it rounds twice.  Number types without an
# `fma` method, such as ReverseDiff's tracked reals, take `muladd`.
@inline fused_abs2(x::T, y::T) where {T<:AbstractFloat} = fma(x, x, y * y)
@inline fused_abs2(x, y) = muladd(x, x, y * y)

# A quarter turn, z / i, of the phase `z`.  For a z near the unit circle whose components
# are both nonzero, that is exactly (Im z, -Re z) for every division in Base — the robust
# one of `ComplexF64`, the widening one of `ComplexF32` and `ComplexF16`, and the generic
# one that `BigFloat` and `Double64` use — so the components are exchanged rather than
# divided, which saves a division that the recurrence below would otherwise wait for.  Where
# a component is zero or not finite, or `z` is far from the circle, where the division of
# `ComplexF64` rescales its arguments, the result, and the signs of its zeros, are left to
# the division.
@inline quarter_turn(z::Complex) = z / 1im
@inline function quarter_turn(z::Complex{T}) where {T<:AbstractFloat}
    x, y = reim(z)
    if !iszero(x) & !iszero(y) & (T(1)/2 ≤ max(abs(x), abs(y)) ≤ 2)
        Complex(y, -x)
    else
        z / 1im
    end
end

# The start of the recurrence: `z` rotated by a power of i into the sector -π/4 < arg z ≤
# π/4, that power `θ`, which is factored out and restored exactly at each power, the first
# increment δz¹ = z² - z, and the constant `t`.  The recurrence is most accurate for z near
# 1, where t is small, and least accurate as Re z approaches 0, where t approaches -2; the
# sector is as close to 1 as a power of i can bring every z, and within it |t| ≤ 2 - √2.
# The rotation is exact, since it only exchanges and negates the components, and z = 1 is
# not rotated.
#
# The sector is half open, so that it and its three rotations by powers of i divide the
# plane (without 0) between them, and z, i z, -z and -i z are all rotated to the same phase.
# The powers of i z and of -z are then exactly iᵐ and (-1)ᵐ times those of z, which is what
# makes identities such as 𝔇(-R) = (-1)^{2ℓ} 𝔇(R) hold exactly.  The price is paid on the
# diagonals, |Re z| = |Im z|, where z and conj(z) are rotated to the same phase rather than
# to conjugate phases, so that only there are the powers of conj(z) not exactly the
# conjugates of the powers of z.  Every nonzero z reaches the sector in at most three
# rotations; the loop is bounded, rather than run until z is in the sector, so that it ends
# for z = 0 too, which only the unchecked `complex_powers!` accepts.
@inline function complex_powers_start(z)
    θ = one(z)
    for _ in 1:3
        -z.re < z.im ≤ z.re && break
        θ *= 1im
        z = quarter_turn(z)
    end
    # dc = -2 (Im √z)² = Re z - |z| = -(Im z)² / (Re z + |z|).  The last form is the one used
    # here: it avoids `Base.sqrt(::Complex)`, whose `nextfloat` rules out element types such
    # as `ForwardDiff.Dual`, and it is free of the cancellation that `Re z - |z|` suffers
    # when `z` is near 1 (which is exactly the small-angle case).  `Re z ≥ |Im z|` and `Re z
    # > 0` here, so the denominator cannot cancel.
    #
    # `modulus` must be computed with a single rounding, by `fused_abs2`.  A one-ulp error
    # here feeds `dc`, which the form above has just gone to some trouble to keep free of
    # cancellation, and the recurrence then amplifies it linearly in `m`.  Measured at m =
    # 4096, ϕ = 0.3: `√(abs2(z))` gives 9.1e-14, this gives 7.7e-15.  (`hypot` does not
    # help; nor does a `muladd` in the recurrence itself.)
    modulus = √(fused_abs2(z.re, z.im))
    dc = -z.im^2 / (z.re + modulus)
    t = 2 * dc
    # The first increment δz¹ = z² - z is (z - 1) + 2dc z for |z| = 1, since then
    # z² = 2 Re z · z - 1, and so it is formed as dc (1 + 2z) + i Im(z - 1).  In the
    # original form of the algorithm the imaginary part of z - 1 is written √(-dc (2 + dc)),
    # which is |Im z|, and so is right only where Im z ≥ 0; `Im z · √((2 + dc) / (Re z +
    # |z|))` equals it there, because `-dc = (Im z)² / (Re z + |z|)`, and it has the sign of
    # Im z, which is negative in half of the sector.  The second form is also the better one
    # to differentiate, because the square root's argument is near 1 rather than near 0: `√`
    # has an infinite derivative at 0, so the first form gives a `NaN` derivative under
    # automatic differentiation at `z = 1`, which is exactly the phase the ring-based
    # transforms use, at `ϕ = 0`.
    dz = dc * (1 + 2 * z) + 1im * (z.im * sqrt((2 + dc) / (z.re + modulus)))
    z, θ, dz, t
end

# One step of the recurrence: from zᵐ and δzᵐ to zᵐ⁺¹ and δzᵐ⁺¹.
@inline function complex_powers_step(zᵐ, δzᵐ, t)
    zᵐ⁺¹ = zᵐ + δzᵐ
    zᵐ⁺¹, δzᵐ + t * zᵐ⁺¹
end

# The state of the recurrence for one phase, for `complex_powers!` below: the power zᵐ of
# the rotated phase, stored as an element of the output would be; the power of i that was
# factored out; the increment δzᵐ; the constant `t`; and θᵐ, by which zᵐ is multiplied on
# its way out.
@inline function complex_powers_state(zpowers, z)
    z¹, θ, dz, t = complex_powers_start(z)
    (convert(eltype(zpowers), z¹), θ, dz, t, θ)
end
# Store zᵐ⁻¹ as element `m` of `zpowers`, and advance the state by one step.
@inline function complex_powers_step!(zpowers, m, (zᵐ, θ, dz, t, θᵐ))
    zᵐ⁺¹, dz = complex_powers_step(zᵐ, dz, t)
    @inbounds zpowers[m] = zᵐ * θᵐ
    (convert(eltype(zpowers), zᵐ⁺¹), θ, dz, t, θᵐ * θ)
end
@inline function complex_powers_store!(zpowers, m, (zᵐ, θ, dz, t, θᵐ))
    @inbounds zpowers[m] = zᵐ * θᵐ
    nothing
end

# The functions that allocate the powers, `complex_powers` and `ComplexPowers`, are refused
# for a `z` whose modulus is not close to 1, for which the recurrence would compute the
# powers of a different number.  The in-place kernel, which the calculators call with phases
# they have just computed from normalized rotors, does not check.
function check_unit_modulus(z)
    if abs(z) ≉ one(z)
        throw(DomainError(z,
            "The powers are computed only for `z` with complex amplitude approximately 1; "
            * "abs(z) = $(abs(z))."
        ))
    end
    z
end


"""
    complex_powers!(zpowers, z)

Compute integer powers of `z` from `z^0` through `z^m`, recursively, where `m` is one less
than the length of the input `zpowers` vector.

Note that `z` is assumed to be normalized, with complex amplitude approximately 1; this is
not checked.  The algorithm, and its accuracy, are described under [`complex_powers`](@ref).

See also: [`complex_powers`](@ref), [`ComplexPowers`](@ref)
"""
function complex_powers!(zpowers, z)
    phase_powers!((zpowers,), (z,))
    zpowers
end

# The powers of several phases at once: each vector of `zpowers` is filled with the powers
# of the phase at the same position of `z`, exactly as `complex_powers!` fills one, and the
# vectors must have the same length.  The recurrence for each phase is a chain of dependent
# operations, so computing several together lets the processor overlap them.  It uses `map`
# rather than `foreach` throughout, because `map` over tuples is unrolled, and the states of
# the recurrences then stay in registers; `foreach` over several tuples goes through `zip`
# and a call that is not inlined.
@inline function phase_powers!(zpowers::Tuple, z::Tuple)
    map(Base.require_one_based_indexing, zpowers)
    M = length(first(zpowers))
    all(v -> length(v) == M, zpowers) || throw(DimensionMismatch(
        "The vectors of powers have lengths $(map(length, zpowers)), which must be equal."
    ))
    M == 0 && return zpowers
    map((v, z) -> (@inbounds v[1] = one(z)), zpowers, z)
    M == 1 && return zpowers
    if M == 2
        map((v, z) -> (@inbounds v[2] = z), zpowers, z)
        return zpowers
    end
    states = map(complex_powers_state, zpowers, z)
    for m ∈ 2:M-1
        states = map(@inline((v, s) -> complex_powers_step!(v, m, s)), zpowers, states)
    end
    map(@inline((v, s) -> complex_powers_store!(v, M, s)), zpowers, states)
    zpowers
end


"""
    complex_powers(z, m)

Compute integer powers of `z` from `z^0` through `z^m`, recursively, and return them as a
vector of length `m+1`.

The number `z` must have complex amplitude approximately 1, or a `DomainError` is thrown.  A
real `z` or one with integer components is converted to a complex floating-point number
first, so that the result is, for example, a `Vector{ComplexF64}` for `z = im`.  The largest
power `m` may be any non-negative integer.

This algorithm is mostly due to Stoer and Bulirsch in "Introduction to Numerical Analysis"
(page 24) — with a little help from de Moivre's formula, which is essentially exp(iθ)ⁿ =
exp(inθ), as well as my own alterations to deal with different behaviors in different
quadrants.

There isn't usually a huge advantage to using this specialized function.  If you just need a
particular power, it will generally be far more efficient and just as accurate to compute
either exp(iθ)ⁿ or exp(inθ) explicitly.  However, if you need all powers from 0 to m, this
function is several times faster than the first of those options, and about twice as fast as
the second, for large m.  Like those options, this function is numerically stable: measured
against the exact powers of the number `z` it was given, its error in ``zᵐ`` grows linearly
in ``m``, and beyond the first few powers, whose errors are a few roundings, it is at most
about ``0.6 m ϵ``, where ``ϵ`` is the machine precision of `z`, and typically about a third
of that, with the largest errors occurring for phases near the odd multiples of π/4.  That
supposes that ``|z|`` is 1 to within about ``ϵ/3``, as it is for a correctly rounded phase
such as `cis(θ)`; a larger departure ``δ = |z| - 1`` adds up to about ``m |δ|`` to the
error.

The powers of `-z` and of `im * z` are exactly ``(-1)ᵐ`` and ``iᵐ`` times those of `z`, and,
except where `abs(real(z)) == abs(imag(z))`, the powers of `conj(z)` are exactly the
conjugates of those of `z`.  In each case the values compare equal with `==`, though the
signs of zeros may differ.

See also: [`complex_powers!`](@ref), [`ComplexPowers`](@ref)
"""
function complex_powers(z::Number, m::Integer)
    if m < 0
        throw(ArgumentError("The largest power `m` must be non-negative; got m=$m."))
    end
    z = check_unit_modulus(complex(float(z)))
    zpowers = zeros(typeof(z), m+1)
    complex_powers!(zpowers, z)
end


struct ComplexPowers{T<:Complex, RT<:Real}
    z¹::T   # `z` rotated into the sector -π/4 < arg z ≤ π/4 by a power of i
    θ::T    # that power of i, so that z = θ z¹
    δz¹::T  # the first increment of the recurrence, (z¹)² - z¹
    t::RT   # its constant, 2(Re z¹ - 1)
end

@doc raw"""
    ComplexPowers(z)

Construct an iterator to compute powers of the complex phase factor ``z``, which must have
magnitude approximately 1.  The iterator will return the complex number ``zᵐ`` for each
integer ``m = 0, 1, 2, \ldots``.

The parameters of the iterator's type, `ComplexPowers{T, RT}`, are as follows:
- `T` is the complex type of the powers.
- `RT` is the real type of their components.

A real `z` or one with integer components is converted to a complex floating-point number
first, and a `z` whose magnitude is not approximately 1 is refused with a `DomainError`.

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
``ϕ``, the power of ``i`` that rotates ``z`` into the sector ``-π/4 < \arg z ≤ π/4``, which
is as close to 1 as such a rotation can bring it.  This power of ``i`` can be separately
exponentiated exactly and efficiently because it is exactly representable as a complex
integer, while the error in the computation of ``zᵐ`` is reduced significantly for certain
values of ``z``.

This is the recurrence of [`complex_powers!`](@ref), which gives the same values; beyond the
first few powers, it achieves a worst-case accuracy of about ``0.6 m ϵ`` for ``zᵐ`` — where
``ϵ`` is the precision of the type of the input argument — across the range of inputs ``z =
e^{iθ}`` for ``θ ∈ [0, 2π]``, as described under [`complex_powers`](@ref).  The original
algorithm can be far worse for values of ``θ`` close to ``π`` — often orders of magnitude
worse.
"""
function ComplexPowers(z::Number)
    z = check_unit_modulus(complex(float(z)))
    z¹, θ, δz¹, t = complex_powers_start(z)
    ComplexPowers(z¹, θ, δz¹, t)
end

# The state is the unrotated power zᵐ, its increment δzᵐ, and θᵐ, by which zᵐ is multiplied
# on its way out.  These are the steps of `complex_powers!`, in the same order, so that the
# two agree bit for bit from the second power on.
Base.iterate(cp::ComplexPowers) = one(cp.z¹), (cp.z¹, cp.δz¹, cp.θ)
function Base.iterate(cp::ComplexPowers, (zᵐ, δzᵐ, θᵐ))
    zᵐ⁺¹, δzᵐ⁺¹ = complex_powers_step(zᵐ, δzᵐ, cp.t)
    zᵐ * θᵐ, (zᵐ⁺¹, δzᵐ⁺¹, θᵐ * cp.θ)
end

Base.IteratorSize(::Type{<:ComplexPowers}) = Base.IsInfinite()

Base.eltype(::Type{<:ComplexPowers{T}}) where {T} = T

Base.isdone(::ComplexPowers) = false
Base.isdone(::ComplexPowers, ::Any) = false
