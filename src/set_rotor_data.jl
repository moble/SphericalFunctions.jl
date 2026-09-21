### Replacing a calculator's rotor data
#
# Every calculator is constructed with its rotor data, so these functions are for the second
# and subsequent values: build one workspace, then walk it over many rotors.  They are thin
# wrappers on the internal `set_rotors!`, which does the work and the validation; what they
# add is a name for each kind of data, so that no name promises something it does not do —
# `set_β!` on a `WignerdCalculator` really does keep only the β angle of a rotor handed to it,
# and `set_θ!` really does mean the (θ, ϕ=0) evaluation rather than a full rotation.

"""
    set_R!(calc, R)

Replace the rotor data of `calc` with the rotor `R`, and return `calc`.  Any results the
calculator was holding are discarded, so the next step of the recurrence starts from
``ℓ_{min}`` again.

`R` must be a `Rotor`, or — for a calculator built for `Nᵣ > 1` rotors — an `AbstractVector`
of exactly that many.  A quaternion that is not a `Rotor` is refused rather than normalized;
see [`not_a_rotor`](@ref).  The number of rotors is fixed at
construction and cannot be changed here; for a different number, build another calculator.

This applies to [`WignerDCalculator`](@ref) and [`sYlmCalculator`](@ref), both of which need
the whole rotor.  The real ``d`` matrices and the ``H`` wedge depend on ``β`` alone, so their
calculators take [`set_β!`](@ref) instead, and the ``(θ, ϕ=0)`` evaluation of the harmonics
takes [`set_θ!`](@ref).

The new rotor's floating-point type must be the one the calculator already works in, which
was fixed by the rotor it was built from.  A mismatch is an error rather than a silent
conversion — narrowing a `BigFloat` rotor into a `Float64` calculator would throw away
precision that nobody chose to throw away — so convert whichever side you meant, or build a
calculator from data of the type you want.  `floattype(calc)` reports the type in use.

```julia
calc = WignerDCalculator(first(rotors), ℓₘₐₓ)
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
function set_R!(c::sYlmCalculator, R::Union{Rotor, AbstractVector})
    check_rotor_type(c, R)
    set_rotors!(c, R)
end
# An `sλlmCalculator` stores real numbers, so there is nowhere to put the α and γ phases a
# rotor specifies.  This mirrors the refusal above for a `WignerdCalculator`, and names the
# setter that does serve it.
function set_R!(::sλlmCalculator, ::Any)
    error(
        "An sλlmCalculator evaluates at (θ, ϕ=0) and stores real numbers, so it has nowhere "
        * "to put the α and γ phases a rotor specifies; its setter is `set_θ!`.  Use an "
        * "sYlmCalculator if the full harmonics are wanted."
    )
end

function set_R!(::WignerCalculator{IT, RT, RT}, ::Any) where {IT, RT<:Real}
    error(
        "A WignerdCalculator holds only the angle β, not a whole rotor, so its setter is "
        * "`set_β!` — which accepts a rotor and keeps just its β.  Use a WignerDCalculator "
        * "if the full 𝔇 matrices are wanted."
    )
end
function set_R!(::WignerHCalculator, ::Any)
    error(
        "A WignerHCalculator holds only the angle β, not a whole rotor, so its setter is "
        * "`set_β!` — which accepts a rotor and keeps just its β."
    )
end

"""
    set_β!(calc, β)

Replace the rotor data of `calc` with the angle `β`, and return `calc`.  Any results the
calculator was holding are discarded.

`β` may be given as the angle itself, as the phase ``e^{iβ}``, or as a `Rotor`, of which only
the ``β`` Euler angle is kept — the ``d`` matrices and the ``H`` wedge depend on nothing else.
For a calculator built for `Nᵣ > 1`, pass an `AbstractVector` of exactly that many of any one
of those forms.

This applies to [`WignerdCalculator`](@ref) and [`WignerHCalculator`](@ref).  Wigner's ``𝔇``
matrices and the spin-weighted harmonics need the whole rotor, so their calculators take
[`set_R!`](@ref) instead.

As for [`set_R!`](@ref), the new value's floating-point type must match the calculator's.

Half-integer ``d`` has period ``4π`` in ``β``, so an angle or a rotor determines it
unambiguously, while a bare phase determines it only up to the double-cover sign
``(-1)^{2ℓ}``; the branch ``β ∈ (-π, π]`` is used in that case.
"""
function set_β!(c::WignerCalculator{IT, RT, RT}, β) where {IT, RT<:Real}
    check_rotor_type(c, β)
    set_rotors!(c, β)
end
function set_β!(w::WignerHCalculator, β)
    check_rotor_type(w, β)
    set_rotors!(w, β)
end

function set_β!(::WignerCalculator{IT, RT, Complex{RT}}, ::Any) where {IT, RT<:Real}
    error(
        "A WignerDCalculator needs the whole rotor — 𝔇 depends on all three Euler angles — "
        * "so its setter is `set_R!`.  Use a WignerdCalculator if only β is available."
    )
end

"""
    set_θ!(calc, θ)

Replace the rotor data of the [`sYlmCalculator`](@ref) `calc` with the angle `θ`, meaning the
point ``(θ, ϕ = 0)``, and return `calc`.  Any results the calculator was holding are
discarded.

This is what gives the real functions ``{}_sλ_{ℓ,m}(θ) = {}_sY_{ℓ,m}(θ, 0)`` that the
ring-based transforms are built on; the values are still stored as complex numbers, with zero
imaginary part for integer spin weight.  For a calculator built for `Nᵣ > 1`, pass an
`AbstractVector` of exactly that many angles.

As for [`set_R!`](@ref), the angle's floating-point type must match the calculator's.

For a general point on the sphere use [`set_R!`](@ref) with
`from_spherical_coordinates(θ, ϕ)`.
"""
function set_θ!(c::HarmonicCalculator, θ::Union{Real, AbstractVector{<:Real}})
    check_rotor_type(c, θ)
    set_rotors!(c, θ)
end

function set_θ!(::WignerCalculator{IT, RT, Complex{RT}}, ::Any) where {IT, RT<:Real}
    error("Only an sYlmCalculator evaluates at (θ, ϕ=0); a WignerDCalculator's setter is `set_R!`.")
end
function set_θ!(::WignerCalculator{IT, RT, RT}, ::Any) where {IT, RT<:Real}
    error("Only an sYlmCalculator evaluates at (θ, ϕ=0); a WignerdCalculator's setter is `set_β!`.")
end
function set_θ!(::WignerHCalculator, ::Any)
    error("Only an sYlmCalculator evaluates at (θ, ϕ=0); a WignerHCalculator's setter is `set_β!`.")
end
