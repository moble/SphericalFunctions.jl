### Replacing a calculator's rotor data
#
# Every calculator is constructed with its rotor data, so these functions are for the second
# and subsequent values: build one workspace, then walk it over many rotors.  They are thin
# wrappers on the internal `set_rotors!`, which does the work and the validation; what they
# add is a name for each kind of data, so that no name promises something it does not do —
# `set_β!` on a `dCalculator` really does keep only the β angle of a rotor handed to it, and
# `set_θ!` really does mean the (θ, ϕ=0) evaluation rather than a full rotation.  A call
# that names the wrong kind of data for a calculator is refused with an `ArgumentError` that
# names the setter that does serve it.

"""
    set_R!(calc, R)

Replace the rotor data of `calc` with the rotor `R`, and return `calc`.  Any results the
calculator was holding are discarded, so the next step of the recurrence starts from
``ℓ_{min}`` again.

`R` must be a `Rotor`, or an `AbstractVector` of exactly as many rotors as the calculator
handles, which is required when the calculator was built from a vector of rotors; a vector
of one rotor is accepted by a calculator built from a single one.  A `Quaternion` may also
stand in for a rotor, as the rotation of its normalization; a `QuatVec` is refused, as
[`not_a_rotor`](@ref) explains.  The reason we *must* accept `Quaternion`s is that some
automatic-differentiation packages require the tangent of a type to have the same type; the
tangent to the space of `Rotor`s is not a `Rotor`, but a `Quaternion`, so we must accept
`Quaternion`s to begin with.  The number of rotors is fixed at construction and cannot be
changed here; for a different number of rotors, build another calculator.

This applies to [`DCalculator`](@ref) and [`sYlmCalculator`](@ref), both of which need the
whole rotor.  The real ``d`` matrices and the ``H`` wedge depend on ``β`` alone, so their
calculators take [`set_β!`](@ref) instead, and the ``(θ, ϕ=0)`` evaluation of the harmonics
takes [`set_θ!`](@ref); an angle given to `set_R!` is refused with a message naming it.

The new rotor's floating-point type must be the one the calculator already works in, which
was fixed by the rotor it was built from.  A mismatch is an error rather than a silent
conversion — narrowing a `BigFloat` rotor into a `Float64` calculator would throw away
precision that nobody chose to throw away — so convert whichever side you meant, or build a
calculator from data of the type you want.  `floattype(calc)` reports the type in use.

Everything is checked before anything is replaced, so data that are refused leave the
calculator exactly as it was, including whatever it had computed.

```julia
calc = DCalculator(first(rotors), ℓₘₐₓ)
for R ∈ rotors
    set_R!(calc, R)
    for (ℓ, 𝔇ˡ) ∈ calc
        # ...
    end
end
```
"""
function set_R!(c::WignerCalculator{IT, RT, Complex{RT}}, R) where {IT, RT<:Real}
    check_rotor_type(c, R)
    set_rotors!(c, R)
end
function set_R!(c::sYlmCalculator, R)
    check_rotor_type(c, R)
    set_rotors!(c, R)
end
# An sYlmCalculator also evaluates at (θ, ϕ=0), but that is what `set_θ!` is named for, and
# `set_R!` does not quietly switch the calculator to it.
function set_R!(::sYlmCalculator, ::Union{Real, AbstractVector{<:Real}})
    throw(ArgumentError(
        "`set_R!` takes rotors; an angle θ means the point (θ, ϕ=0), for which the setter "
        * "is `set_θ!`.  For a general point use `set_R!` with "
        * "`from_spherical_coordinates(θ, ϕ)`."
    ))
end
# An `sλlmCalculator` stores real numbers, so there is nowhere to put the α and γ phases a
# rotor specifies.  This mirrors the refusal below for a `dCalculator`, and names the setter
# that does serve it.
function set_R!(::sλlmCalculator, ::Any)
    throw(ArgumentError(
        "An sλlmCalculator evaluates at (θ, ϕ=0) and stores real numbers, so it has nowhere "
        * "to put the α and γ phases a rotor specifies; its setter is `set_θ!`.  Use an "
        * "sYlmCalculator if the full harmonics are wanted."
    ))
end

function set_R!(::WignerCalculator{IT, RT, RT}, ::Any) where {IT, RT<:Real}
    throw(ArgumentError(
        "A dCalculator holds only the angle β, not a whole rotor, so its setter is "
        * "`set_β!` — which accepts a rotor and keeps just its β.  Use a DCalculator "
        * "if the full 𝔇 matrices are wanted."
    ))
end
function set_R!(::HCalculator, ::Any)
    throw(ArgumentError(
        "An HCalculator holds only the angle β, not a whole rotor, so its setter is "
        * "`set_β!` — which accepts a rotor and keeps just its β."
    ))
end

"""
    set_β!(calc, β)
    set_beta!(calc, β)

Replace the rotor data of `calc` with the angle `β`, and return `calc`.  Any results the
calculator was holding are discarded.  `set_beta!` is an ASCII alias of the same function.

`β` may be given as the angle itself, as the phase ``e^{iβ}``, or as a `Rotor` or other
`Quaternion`, of which only the ``β`` Euler angle is kept — the ``d`` matrices and the
``H`` wedge depend on nothing else.  A `Rotor` contributes the ``β ∈ [0, π]`` of its
canonical Euler decomposition; if the rotor was built from a ``β`` outside that range, that
``β`` is folded into its ``α`` and ``γ``, which are discarded here.  For a calculator built
from a vector, pass an `AbstractVector` of exactly as many of any one of those forms; a
vector of one is accepted by a calculator built from a single value.

This applies to [`dCalculator`](@ref) and [`HCalculator`](@ref).  Wigner's ``𝔇`` matrices
and the spin-weighted harmonics need the whole rotor, so their calculators take
[`set_R!`](@ref) instead.

As for [`set_R!`](@ref), the new value's floating-point type must match the calculator's,
and data that are refused leave the calculator as it was.  An angle that is infinite or NaN
has no phase, and a phase whose modulus is not 1 is not one; both are refused with a
`DomainError`.

Half-integer ``d`` has period ``4π`` in ``β``, so an angle or a rotor determines it
unambiguously, while a bare phase determines it only up to the double-cover sign
``(-1)^{2ℓ}``; the branch ``β ∈ (-π, π]`` is used in that case.
"""
function set_β!(c::WignerCalculator{IT, RT, RT}, β) where {IT, RT<:Real}
    check_rotor_type(c, β)
    set_rotors!(c, β)
end
function set_β!(w::HCalculator, β)
    check_rotor_type(w, β)
    set_rotors!(w, β)
end

function set_β!(::WignerCalculator{IT, RT, Complex{RT}}, ::Any) where {IT, RT<:Real}
    throw(ArgumentError(
        "A DCalculator needs the whole rotor — 𝔇 depends on all three Euler angles — "
        * "so its setter is `set_R!`.  Use a dCalculator if only β is available."
    ))
end

"""
    set_θ!(calc, θ)
    set_theta!(calc, θ)

Replace the rotor data of the [`sYlmCalculator`](@ref) or [`sλlmCalculator`](@ref) `calc`
with the angle `θ`, meaning the point ``(θ, ϕ = 0)``, and return `calc`.  Any results the
calculator was holding are discarded.  `set_theta!` is an ASCII alias of the same
function.

For an `sYlmCalculator` the values are then ``{}_sY_{ℓ,m}(θ, 0)``, stored as complex
numbers.  For an integer spin weight these are the real functions ``{}_sλ_{ℓ,m}(θ)`` that
the ring-based transforms are built on, with zero imaginary part; for a half-odd one they
are ``i^{2s}`` times ``{}_sλ_{ℓ,m}(θ)``, and so imaginary.  An `sλlmCalculator` gives
``{}_sλ_{ℓ,m}(θ)`` directly, as real numbers.  For a calculator built from a vector, pass an
`AbstractVector` of exactly as many angles; a vector of one is accepted by a calculator
built from a single angle.

As for [`set_R!`](@ref), the angle's floating-point type must match the calculator's, and
data that are refused leave the calculator as it was.  An angle that is infinite or NaN is
refused with a `DomainError`, and a rotor with an `ArgumentError` that names the setter for
rotors.

For a general point on the sphere use [`set_R!`](@ref) with `from_spherical_coordinates(θ,
ϕ)`.
"""
function set_θ!(c::HarmonicCalculator, θ::Union{Real, AbstractVector{<:Real}})
    check_rotor_type(c, θ)
    set_rotors!(c, θ)
end
function set_θ!(::sYlmCalculator, ::Union{RotorLike, AbstractVector{<:RotorLike}})
    throw(ArgumentError(
        "`set_θ!` takes angles θ, meaning the points (θ, ϕ=0); the setter for rotors is "
        * "`set_R!`."
    ))
end
function set_θ!(::sλlmCalculator, ::Union{RotorLike, AbstractVector{<:RotorLike}})
    throw(ArgumentError(
        "`set_θ!` takes angles θ, meaning the points (θ, ϕ=0).  An sλlmCalculator stores "
        * "real numbers, so it has nowhere to put the α and γ phases a rotor specifies; use "
        * "an sYlmCalculator and `set_R!` if the full harmonics are wanted."
    ))
end

function set_θ!(::WignerCalculator{IT, RT, Complex{RT}}, ::Any) where {IT, RT<:Real}
    throw(ArgumentError(
        "Only an sYlmCalculator or an sλlmCalculator evaluates at (θ, ϕ=0); a DCalculator's "
        * "setter is `set_R!`."
    ))
end
function set_θ!(::WignerCalculator{IT, RT, RT}, ::Any) where {IT, RT<:Real}
    throw(ArgumentError(
        "Only an sYlmCalculator or an sλlmCalculator evaluates at (θ, ϕ=0); a dCalculator's "
        * "setter is `set_β!`."
    ))
end
function set_θ!(::HCalculator, ::Any)
    throw(ArgumentError(
        "Only an sYlmCalculator or an sλlmCalculator evaluates at (θ, ϕ=0); an "
        * "HCalculator's setter is `set_β!`."
    ))
end
