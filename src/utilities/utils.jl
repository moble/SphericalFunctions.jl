# The logarithm of `binomial(n, k)`, computed in the float type `S`, which is what
# `sqrtbinomial` exponentiates.  The general case is `-log(n+1) - log B(n-k+1, k+1)`, where
# B is the beta function, as in `SpecialFunctions.logabsbinomial`.  In `Float64` this is
# several times more accurate at large `n` than the difference of three `loggamma`s, which
# are large and nearly cancel, and unlike `logabsbinomial` it accepts any float type.  The
# coefficient is symmetric in `k ↔ n - k`, so the smaller of the two is used.
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
        return -log1p(S(n)) - SpecialFunctions.logbeta(S(n - k + 1), S(k + 1))
    end
end

"""
    sqrtbinomial(n, k, [T=Float64])

The square root of the binomial coefficient `binomial(n, k)`, computed in the float type `T`
from its logarithm, so that it is finite where the coefficient itself would overflow.

For `Int` arguments, `binomial` overflows at about `n = 66` when `k ≈ n/2`, but the square
root of the coefficient, which is what many normalization constants need, is representable
in `Float64` up to about `n = 2050`.  The logarithm is computed through the beta function,
as it is by [`logabsbinomial` in
SpecialFunctions.jl](https://specialfunctions.juliamath.org/latest/functions_list/#SpecialFunctions.logabsbinomial),
but in any float type `T` that `SpecialFunctions.logbeta` accepts, including `BigFloat` and
`Double64`.  Exponentiating the logarithm magnifies its rounding error in proportion to its
size, so the relative error is a few ulp for small coefficients and grows with their
logarithm, to about a thousand ulp near `n = 2050` in `Float64`.  `n` and `k` are integers
of one type, and the result is zero when `k` is negative or greater than `n`.
"""
function sqrtbinomial(n, k, ::Type{T}=Float64) where T
    exp(logbinomial(n, k, T)/2)
end


### Rotor data supplied to a calculator
#
# Every calculator takes its rotor data as a constructor argument, which fixes both the
# floating-point type it works in and the number of rotors it handles at once.  The two
# helpers below derive those from the argument.  Both dispatch on the argument's *type*: a
# runtime reduction over a vector would infer only as `Type`, which would make the element
# type a runtime value and cost the calculators their type stability.

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

"""
    rotor_basetype(R)

The floating-point type in which to work, given the rotor data `R`: the `float` of its
component type.  This is the *only* thing that decides the element type a calculator works
in, so that a `Rotor{Float32}` gives a `Float32` calculator and a `Rotor{Double64}` a
`Double64` one.  A `Quaternion` counts as a rotor here; we always divide out by the norm,
and accept `Quaternion` for compatibility with various automatic-differentiation packages.
To compute in some other type, convert the rotor data — which is also the honest way to say
it, since the type of the data is the claim being made about the points.

The rotor data must therefore commit to a type to go on: data whose component type is
abstract, such as a `Complex{Real}` or a `Vector{Rotor{Real}}`, or a vector whose element
type is abstract or a `Union`, is refused with an `ArgumentError` rather than guessed at.
"""
rotor_basetype(R::RotorLike) = concrete_float(Quaternionic.basetype(R), R)
rotor_basetype(β::Real) = float(typeof(β))
rotor_basetype(z::Complex{T}) where {T<:Real} = concrete_float(T, z)
rotor_basetype(R::AbstractVector{<:Rotor{T}}) where {T<:Real} = concrete_float(T, R)
rotor_basetype(R::AbstractVector{<:Quaternionic.Quaternion{T}}) where {T<:Real} = concrete_float(T, R)
rotor_basetype(β::AbstractVector{T}) where {T<:Real} = concrete_float(T, β)
rotor_basetype(z::AbstractVector{<:Complex{T}}) where {T<:Real} = concrete_float(T, z)
function rotor_basetype(R)
    throw(ArgumentError(
        "Cannot build a calculator from rotor data of type $(typeof(R)); $rotor_input_forms."
    ))
end
# `Vector{Rotor}` and `Vector{Rotor{<:Real}}` hold rotations, but their element type does not
# say in what precision, so they get the message above rather than the one below; so do the
# vectors of `Quaternion`s like them.
function rotor_basetype(R::AbstractVector{<:RotorLike})
    throw(ArgumentError(
        "Cannot build a calculator from rotor data of type $(typeof(R)); $rotor_input_forms."
    ))
end
# A `QuatVec` is refused here, which is where the calculators catch it: the constructors
# call this directly, and the setters through `check_rotor_type`.
rotor_basetype(R::Union{QuatVec, AbstractVector{<:QuatVec}}) = throw(ArgumentError(not_a_rotor(R)))

# The float type of rotor data `R` whose components are of type `T`.  `float(Real)` is
# `Float64`, so without the check data of an abstract component type would be answered with
# a guess — the one thing `rotor_basetype` is not allowed to do.  A concrete `T` is fine
# even when it is not itself a float: `float(Int)` is a derivation, exactly as it is for a
# single `Int` angle.  The check is on a type known when the method is compiled, so it folds
# away.
@inline function concrete_float(::Type{T}, R) where {T}
    isconcretetype(T) || throw(abstract_components_error(T, R))
    working_type(float(T))
end
@noinline function abstract_components_error(T, R)
    ArgumentError(
        "The rotor data, of type $(typeof(R)), has components of type $T, which does not say "
        * "what floating-point type to work in; convert the data to a concrete type first."
    )
end

# Quaternions that do not denote rotations — `QuatVec`s — singly or in a vector.  The
# functions that take rotations — `D`, `d`, `sYlm`, `Ylm`, `sYlm_matrix`, `w(R)` and the
# transforms — never reach `rotor_basetype` with these, so each has a method on this type
# that refuses them with the message of `not_a_rotor` rather than a bare `MethodError`.
const NonRotorData = Union{QuatVec, AbstractVector{<:QuatVec}}

"""
    not_a_rotor(R)

The message for rotor data that is a quaternion but does not denote a rotation, which is to
say a `QuatVec`.

These functions are defined on the rotation group, so a rotation is what they take: a
`Rotor`, or any other `Quaternion`, which denotes the rotation of its normalization, since
the recurrence divides out its magnitude.  A `QuatVec` is a different thing — it represents
a vector, and reading one as a rotation by ``π`` about its own direction would be a category
error rather than a convenience; `exp(v/2)` gives the rotation that a vector generates.
"""
function not_a_rotor(R)
    T = R isa AbstractVector ? eltype(R) : typeof(R)
    (
        (
            R isa AbstractVector ? "These are `$T`s, which are vectors, not rotations" :
                "A `$T` is a vector, not a rotation"
        )
        * ".  Rotations are taken as `Rotor`s or `Quaternion`s, and `exp(v/2)` gives the "
        * "rotation that a `QuatVec` generates."
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


"""
    floattype(calc)

The floating-point type a calculator works in, which is the type of the rotor data it was
constructed from.  Its results are that type, or `Complex` of it.

There is no way to set this independently of the data: converting the data is how one asks
for a different type, and the methods that replace a calculator's data — [`set_R!`](@ref),
[`set_β!`](@ref), [`set_θ!`](@ref) — require the new data to agree with what this reports.
"""
function floattype end

"""
    check_rotor_type(calc, data)

Throw an `ArgumentError` unless `data` would give the element type `calc` already works in.

A calculator's element type is fixed by the data it was constructed from, so replacing that
data later — through [`set_R!`](@ref), [`set_β!`](@ref), [`set_θ!`](@ref) or `similar(calc,
data)` — cannot change it.  Rather than convert silently, which is how precision gets lost
without anyone choosing to lose it, the mismatch is an error and the caller converts
whichever side they meant.  `floattype` is defined alongside each calculator.
"""
function check_rotor_type(calc, data)
    RT = rotor_basetype(data)
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
