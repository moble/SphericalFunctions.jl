# `loggamma` and `logbinomial` exist to support `sqrtbinomial`.  As of 3.0 nothing inside
# the package calls any of the three — `Deprecated`, deleted in 3.0, was the last caller —
# but `sqrtbinomial` is documented on `docs/src/20-interface/05-utilities.md` and §9 of the v3
# design memo directs callers to it in preference to `binomial`, so they are kept and
# tested (the "Combinatorics: sqrtbinomial and logbinomial" test item).
loggamma(a, ::Type{T}) where T = SpecialFunctions.loggamma(T(a))
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
        return (
            -log1p(S(n)) - loggamma(n - k + one(T), S)
            - loggamma(k + one(T), S) + loggamma(n + 2one(T), S)
        )
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

const _rotor_input_forms = (
    "the accepted forms are a Rotor, the angle β::Real, or the "
    * "phase e^{iβ}::Complex — or, for Nᵣ > 1, an AbstractVector of length Nᵣ of any one of "
    * "those, whose element type must say which (a `Vector{Any}`, or one with a `Union` "
    * "element type, does not, and should be converted before it is passed)"
)

"""
    rotor_basetype(R)

The floating-point type in which to work, given the rotor data `R`: the `float` of its
component type.  This is the *only* thing that decides the element type a calculator works
in, so that a `Rotor{Float32}` gives a `Float32` calculator and a `Rotor{Double64}` a
`Double64` one.  To compute in some other type, convert the rotor data — which is also the
honest way to say it, since the type of the data is the claim being made about the points.

The rotor data must therefore commit to a type to go on: a vector whose element type is abstract
or a `Union` is rejected rather than guessed at.
"""
rotor_basetype(R::Rotor) = float(Quaternionic.basetype(R))
rotor_basetype(β::Real) = float(typeof(β))
rotor_basetype(::Complex{T}) where {T<:Real} = float(T)
rotor_basetype(::AbstractVector{<:Rotor{T}}) where {T<:Real} = float(T)
function rotor_basetype(R::AbstractVector{T}) where {T<:Real}
    # `float(Real)` is `Float64`, so without this an abstractly-typed vector of angles would
    # be answered with a guess — the one thing this function is not allowed to do.  A
    # concrete element type is fine even when it is not itself a float: `float(Int)` is a
    # derivation, exactly as it is for a single `Int` angle.
    if !isconcretetype(T)
        error(
            "The element type of the given angles is $T, which does not say what "
            * "floating-point type to work in; convert them to a concrete type first."
        )
    end
    float(T)
end
rotor_basetype(::AbstractVector{<:Complex{T}}) where {T<:Real} = float(T)
function rotor_basetype(R)
    error("Cannot build a calculator from rotor data of type $(typeof(R)); $_rotor_input_forms.")
end
# `Vector{Rotor}` and `Vector{Rotor{<:Real}}` hold rotations, but their element type does not
# say in what precision, so they get the message above rather than the one below.
function rotor_basetype(R::AbstractVector{<:Rotor})
    error("Cannot build a calculator from rotor data of type $(typeof(R)); $_rotor_input_forms.")
end
# A quaternion that is not a `Rotor` is refused here, which is where every entry point
# catches it: the constructors call this directly, the setters through `check_rotor_type`, and
# `D`, `d`, `sYlm`, `Ylm` and `sYlm_matrix` because each of them builds a calculator.  Adding
# per-function refusals would only introduce dispatch ambiguities.
rotor_basetype(R::Union{AbstractQuaternion, AbstractVector{<:AbstractQuaternion}}) = error(not_a_rotor(R))

"""
    not_a_rotor(R)

The message for rotor data that is a quaternion but not a `Rotor`.

These functions are defined on the rotation group, so a rotation is what they take, and
`Rotor` is the type that says a quaternion is one.  A general `Quaternion` has a magnitude
that the recurrence would simply divide out; `rotor(q)` normalizes it into the rotation it
denotes, and says so at the call site.  A `QuatVec` is further still from a rotation — it
represents a vector, and reading one as a rotation by ``π`` about its own direction would be a
category error rather than a convenience; `exp(v/2)` gives the rotation a vector generates.
"""
function not_a_rotor(R)
    T = R isa AbstractVector ? eltype(R) : typeof(R)
    (
        (R isa AbstractVector ? "These are `$T`s, which are not `Rotor`s" : "A `$T` is not a `Rotor`")
        * ".  Rotations are taken as `Rotor`s, which is what says a quaternion denotes one: "
        * "`rotor(q)` normalizes a `Quaternion` into the rotation it points at, and `exp(v/2)` "
        * "gives the rotation a `QuatVec` generates."
    )
end

"""
    nrotors(R)

The number of rotors `Nᵣ` that the rotor data `R` describes: one for a single rotor, angle or
phase, and `length(R)` for a vector of them.  This is how a calculator learns its batch size,
which is why there is no `Nᵣ` keyword argument on any constructor.
"""
nrotors(::Number) = 1  # `AbstractQuaternion <: Number`, so this covers a single rotor too
function nrotors(R::AbstractVector)
    if isempty(R)
        error("A calculator needs at least one rotor, but got an empty $(typeof(R)).")
    end
    length(R)
end
function nrotors(R)
    error("Cannot build a calculator from rotor data of type $(typeof(R)); $_rotor_input_forms.")
end


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

Throw unless `data` would give the element type `calc` already works in.

A calculator's element type is fixed by the data it was constructed from, so replacing that
data later — through [`set_R!`](@ref), [`set_β!`](@ref), [`set_θ!`](@ref) or
`similar(calc, data)` — cannot change it.  Rather than convert silently, which is how
precision gets lost without anyone choosing to lose it, the mismatch is an error and the
caller converts whichever side they meant.  `floattype` is defined alongside each calculator.
"""
function check_rotor_type(calc, data)
    RT = rotor_basetype(data)
    if RT !== floattype(calc)
        error(
            "This calculator works in $(floattype(calc)), but the given data would give "
            * "$RT.  A calculator's element type is fixed by the data it was built from; "
            * "convert the data to $(floattype(calc)), or build a calculator from data of "
            * "the type you want."
        )
    end
    nothing
end
