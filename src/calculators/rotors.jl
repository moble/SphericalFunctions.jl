### Rotor data supplied to a calculator
#
# Every calculator takes its rotor data as a constructor argument, which fixes both the
# floating-point type it works in and the number of rotors it handles at once.  The two
# functions below, `floattype` and `nrotors`, derive those from the argument.  The first is
# defined on the argument's *type*: a runtime reduction over a vector would infer only as
# `Type`, which would make the element type a runtime value and cost the calculators their
# type stability.

# The rotor data that denote rotations: a `Rotor`, or any other `Quaternion`, which denotes
# the rotation of its normalization.  The recurrence divides out the norm in any case.
const RotorLike = Union{Rotor, Quaternionic.Quaternion}

# A rotation as a `Quaternion` with components of type `T`, read by its components alone.
@inline as_quaternion(::Type{T}, R) where {T} = Quaternion{T}(R[1], R[2], R[3], R[4])

const rotor_input_forms = (
    "the accepted forms are a Rotor or Quaternion, the angle β::Real, or the phase "
    * "e^{iβ}::Complex — or, "
    * "for a batch, a non-empty AbstractVector of any one of those, which makes a batch "
    * "however short, and whose element type must say which (a `Vector{Any}`, or one with a "
    * "`Union` element type, does not, and should be converted before it is passed)"
)

# Quaternions that do not denote rotations — `QuatVec`s — singly or in a vector.  The
# functions that take rotations — `D`, `d`, `sYlm`, `Ylm`, `sYlm_matrix`, `w(R)` and the
# transforms — never reach `floattype` with these, so each has a method on this type that
# refuses them with the message of `not_a_rotor` rather than a bare `MethodError`.
const NonRotorData = Union{QuatVec, AbstractVector{<:QuatVec}}

"""
    not_a_rotor(R)

The message for rotor data `R`, or for rotor data of type `R`, that is a quaternion but does
not denote a rotation, which is to say a `QuatVec`, or a vector of them.

These functions are defined on the rotation group, so a rotation is what they take: a
`Rotor`, or any other `Quaternion`, which denotes the rotation of its normalization, since
the recurrence divides out its magnitude.  A `QuatVec` is a different thing — it represents
a vector, and reading one as a rotation by ``π`` about its own direction would be a category
error rather than a convenience; `exp(v/2)` gives the rotation that a vector generates.
"""
not_a_rotor(R) = not_a_rotor(typeof(R))
function not_a_rotor(::Type{D}) where {D}
    T = D <: AbstractVector ? eltype(D) : D
    (
        (
            D <: AbstractVector ? "These are `$T`s, which are vectors, not rotations" :
                "A `$T` is a vector, not a rotation"
        )
        * ".  Rotations are taken as `Rotor`s or `Quaternion`s, and `exp(v/2)` gives the "
        * "rotation that a `QuatVec` generates."
    )
end

floattype(x) = floattype(typeof(x))
function floattype(::Type{D}) where {D}
    T = component_type(D)
    T === nothing && throw(not_rotor_data_error(D))
    isconcretetype(T) || throw(abstract_components_error(T, D))
    floattype(T)
end
# `float(Real)` is `Float64`, so without the check a type of rotor data whose components are
# abstract would be answered with a guess — the one thing `floattype` is not allowed to do.
# A concrete `T` is fine even when it is not itself a float: `float(Int)` is a derivation,
# exactly as it is for a single `Int` angle.  The check is on a type known when the method
# is compiled, so it folds away.
function floattype(::Type{T}) where {T<:Real}
    isconcretetype(T) || throw(abstract_components_error(T, T))
    float(T)
end

# The type of the components of rotor data of type `D`, or `nothing` if `D` is not a type
# of rotor data.  A vector is rotor data only if its element type says which kind of rotor
# data it holds, and whether that type is concrete is left for `floattype` to check, so that
# it can refuse an abstract one with a message of its own.  `Vector{Rotor}` and
# `Vector{Rotor{<:Real}}` hold rotations, but their element type does not say in what
# precision, so they are not rotor data here; neither are the vectors of `Quaternion`s like
# them.
component_type(::Type) = nothing
component_type(::Type{T}) where {T<:Real} = T
component_type(::Type{Complex{T}}) where {T<:Real} = T
component_type(::Type{Q}) where {T<:Real, Q<:Union{Rotor{T}, Quaternion{T}}} =
    Quaternionic.basetype(Q)
component_type(::Type{<:AbstractVector{E}}) where {E<:Union{Real, Complex, RotorLike}} =
    component_type(E)

# A `QuatVec` is refused here, which is where the calculators catch it: the constructors
# call `floattype` directly, and the setters through `check_rotor_type`.
@noinline function not_rotor_data_error(::Type{D}) where {D}
    ArgumentError(
        D <: NonRotorData ? not_a_rotor(D) :
            "Cannot build a calculator from rotor data of type $D; $rotor_input_forms."
    )
end
@noinline function abstract_components_error(T, D)
    ArgumentError(
        "The rotor data, of type $D, has components of type $T, which does not say what "
        * "floating-point type to work in; convert the data to a concrete type first."
    )
end

"""
    nrotors(R)

The number of rotors `Nᵣ` that the rotor data `R` describes: one for a single rotor, angle
or phase, and `length(R)` for a vector of them.  This is how a calculator learns its batch
size, which is why there is no `Nᵣ` keyword argument on any constructor.
"""
nrotors(::Number) = 1  # `AbstractQuaternion <: Number`, so this covers a single rotor too
function nrotors(R::AbstractVector)
    if isempty(R)
        throw(ArgumentError(
            "A calculator needs at least one rotor, but got an empty $(typeof(R))."
        ))
    end
    length(R)
end
function nrotors(R)
    throw(ArgumentError(
        "Cannot build a calculator from rotor data of type $(typeof(R)); $rotor_input_forms."
    ))
end

# Whether a calculator built from the rotor data `R` is batched, with blocks that have a
# leading rotor index: exactly when `R` is a vector, however long.  A `Val`, so that the
# calculator's type — and with it the type of its blocks — is known at compile time.
batched_data(::AbstractVector) = Val(true)
batched_data(::Any) = Val(false)

# The rotor data that replace a calculator's own, or that `similar(calc, R)` is given, must
# describe as many rotors as the calculator handles.  A single rotor, angle, or phase fills
# a calculator that handles exactly one, as a calculator built from a vector of one does,
# and the message names it as a single rotor.
function check_rotor_count(c, R)
    if nrotors(R) != Nᵣ(c)
        throw(DimensionMismatch(
            "This calculator handles Nᵣ=$(Nᵣ(c)) rotors, but "
            * (R isa AbstractVector ? "got $(length(R))." : "a single rotor was given.")
        ))
    end
    nothing
end

"""
    check_rotor_type(calc, data)

Throw an `ArgumentError` unless `data` would give the element type `calc` already works in.

A calculator's element type is fixed by the data it was constructed from, so replacing that
data later — through [`set_R!`](@ref), [`set_β!`](@ref), [`set_θ!`](@ref) or `similar(calc,
data)` — cannot change it.  Rather than convert silently, which is how precision gets lost
without anyone choosing to lose it, the mismatch is an error and the caller converts
whichever side they meant.  Both sides are given by [`floattype`](@ref), whose method for
each calculator is defined alongside it.
"""
function check_rotor_type(calc, data)
    RT = floattype(data)
    if RT !== floattype(calc)
        throw(ArgumentError(
            "This calculator works in $(floattype(calc)), but the given data would give "
            * "$RT.  A calculator's element type is fixed by the data it was built from; "
            * "convert the data to $(floattype(calc)), or build a calculator from data of "
            * "the type you want."
        ))
    end
    nothing
end


### Rotor phases
#
# What the recurrence and the phases of 𝔇 and of the harmonics need from each rotor, which
# every calculator stores when its rotors are set.

"""
    spinor_phases(R::AbstractQuaternion, [F])

Return `(eⁱᵝ, z₊, z₋, cβ½, sβ½)` for the rotor `R`, where ``β`` is the Euler angle, ``z₊ =
e^{i(α+γ)/2}``, ``z₋ = e^{i(α-γ)/2}``, ``cβ½ = \\cos(β/2)`` and ``sβ½ = \\sin(β/2)``.  These
are the quantities the Wigner recurrences need: ``e^{i(m′α + mγ)} = z₊^{m′+m} z₋^{m′-m}``,
with integer exponents even for half-integer ``m′, m``, while the half-angle pair seeds the
half-integer recurrence.  At the poles ``β ∈ \\{0, π\\}`` the undefined phase is set to 1;
the corresponding ``d`` elements vanish, so the choice is immaterial.

The half-angles are taken as ``(\\sqrt{W²+Z²}, \\sqrt{X²+Y²})/\\|R\\|``, which is accurate
near both poles and is non-negative, so ``β ∈ [0, π]``; a rotor's double-cover sign is
encoded entirely in `z₊` and `z₋`, giving ``𝔇(-R) = -𝔇(R)`` for half-integer indices.

Callers that need only the first few outputs may drop the rest: `eⁱᵝ, z₊, z₋ =
spinor_phases(R, F)`.

`R` need not be normalized.

These phases are not differentiable at ``β = 0`` or ``β = π``, where `√b` or `√a` is taken
of an exact zero, and a derivative taken through them there by automatic differentiation is
`NaN`.  No local rule can repair this — `sβ½ * z₋` is smooth in the rotor, but `sβ½` and
`z₋` separately are not, so treating the zero as exact would silently drop the first-order
term.  Near a pole the derivatives are finite but inaccurate, the ``k``-th by about ``ε
r^{-k}`` relative to their size at a distance ``r`` from it.  The values are accurate at
every rotor, the poles included, so the calculators of ``𝔇`` and of ``{}_sY_{ℓ,m}``, which
are smooth there, are never differentiated through this: the rules for automatic
differentiation give their derivatives in terms of their values (see
`src/derivatives/kernels.jl`).  ``H`` of a rotor has no such rules, and keeps the `NaN`.
``d`` of a rotor is differentiated through ``β``, by the same rules as ``d`` of an angle.
Since ``β`` has a cone-shaped singularity at each pole, the derivatives of those elements of
``d`` that actually have no derivative there are `NaN`, and those of the others are zero.

The optional second argument is the real type the phases are computed in; it defaults to
`float(eltype(R))`.  Pass the *calculator's* type whenever that is more precise than the
rotor's own — otherwise every later step inherits the rotor type's precision.
"""
function spinor_phases end

spinor_phases(R::AbstractQuaternion{T}) where {T} = spinor_phases(R, float(T))
function spinor_phases(R::AbstractQuaternion, ::Type{F}) where {F<:Real}
    a = F(R[1])^2 + F(R[4])^2
    b = F(R[2])^2 + F(R[3])^2
    sqrta = √a
    sqrtb = √b
    z₊ = iszero(sqrta) ? one(Complex{F}) : Complex{F}(F(R[1]), F(R[4])) / sqrta  # exp[i(α+γ)/2]
    z₋ = iszero(sqrtb) ? one(Complex{F}) : Complex{F}(F(R[3]), -F(R[2])) / sqrtb  # exp[i(α-γ)/2]
    eⁱᵝ = Complex{F}(a - b, 2 * sqrta * sqrtb) / (a + b)
    nrm = √(a + b)
    cβ½ = sqrta / nrm
    sβ½ = sqrtb / nrm
    (eⁱᵝ, z₊, z₋, cβ½, sβ½)
end

# The angle β of rotor data: the angle itself, the argument of a phase e^{iβ}, and the β ∈
# [0, π] of a rotor's Euler decomposition, 2 atan(√(X²+Y²), √(W²+Z²)), as the half angles of
# `spinor_phases` give it.  A calculator of `d` or of ₛλₗₘ keeps a copy of these angles,
# whose tangents and cotangents are what the rules for automatic differentiation read (see
# `src/derivatives/kernels.jl`); its values are computed from the rotor data, not from
# these.
rotation_angle(β::Real) = β
rotation_angle(z::Complex) = angle(z)
rotation_angle(R::RotorLike) = 2atan(sqrt(R[2]^2 + R[3]^2), sqrt(R[1]^2 + R[4]^2))

# The components of rotor data, as a tuple of one, two, or four reals, and the rotor data
# with those components.
rotor_data_components(β::Real) = (β,)
rotor_data_components(z::Complex) = (real(z), imag(z))
rotor_data_components(R::RotorLike) = (R[1], R[2], R[3], R[4])
rotor_data(β::Real) = β
rotor_data(x::Real, y::Real) = Complex(x, y)
rotor_data(w::Real, x::Real, y::Real, z::Real) = Quaternion(w, x, y, z)

"""
    half_angles(eⁱᵝ)

The pair ``(\\cos(β/2), \\sin(β/2))`` for the branch ``β ∈ (-π, π]`` of the phase
``e^{iβ}``, computed without cancellation at either pole.

A bare phase fixes ``β`` only modulo ``2π``, so for half-integer indices this fixes ``d``
only up to the double-cover sign ``(-1)^{2ℓ}``; pass the angle ``β`` itself or a `Rotor` if
that matters.
"""
@inline function half_angles(z::Complex{RT}) where {RT<:Real}
    cosβ, sinβ = reim(z)
    if cosβ ≥ 0
        c = √((1 + cosβ) / 2)
        s = sinβ / (2c)
    else
        s = copysign(√((1 - cosβ) / 2), sinβ)
        c = sinβ / (2s)
    end
    (c, s)
end
