### Rotor data, and rotors that hold derivatives
#
# A calculator of 𝔇 or of the harmonics whose rotors hold derivatives — dual numbers, say —
# runs the recurrence on the values of those rotors, and lifts each block of values into a
# block that holds the derivatives, which the angular-momentum operators give in terms of
# the values (see `src/derivatives/kernels.jl`).  This file holds, in order, what that needs
# to know about a number type, as functions that an extension for a tool such as ForwardDiff
# extends; the data that a lifting calculator keeps; the copying of rotor data into its
# storage; and the view of a vector of rotor data that holds derivatives as a vector of
# their values.  `lift!` itself is in `src/derivatives/kernels.jl`.

# The type of the values of a real type that holds derivatives, and the type itself for one
# that does not.  A calculator whose real type differs from its `value_type` lifts the
# blocks of a calculator of that type.  `recurrence_type` is the type at the bottom of that
# chain, in which the recurrence finally runs.
value_type(::Type{T}) where {T<:Real} = T
recurrence_type(::Type{T}) where {T<:Real} = value_type(T) === T ? T : recurrence_type(value_type(T))

# The value of a number that holds derivatives, and the number itself otherwise; the
# rotation of those values; the number of directions in which a number of type `T` holds
# derivatives; the tangents of a quaternion in each of those directions, as a tuple of
# tuples of four components; and the function that assembles an element of type `Complex{T}`
# from a value and the tuple of its derivatives in those directions.  An extension defines
# all of these for its number type.
real_value(x::Real) = x
rotor_value(q::Quaternion{T}) where {T} =
    Quaternion{value_type(T)}(real_value(q[1]), real_value(q[2]), real_value(q[3]), real_value(q[4]))
function ndirections end
function rotor_tangents end
function lift_combine end

# The data of a calculator that lifts the blocks of another: that calculator, of the values
# of its rotors, and the generators of its rotors' derivatives, with the generator of rotor
# iᵣ in direction d in `G[3d-2:3d, iᵣ]` (see `set_generators!`).
# - `C` is the type of the calculator of the values.
# - `M` is the type of the matrix of generators.
struct Lift{C, M<:AbstractMatrix}
    inner::C
    G::M
end
allocate_lift(::Type{RT}, inner, Nᵣ::Int) where {RT} =
    Lift(inner, Matrix{value_type(RT)}(undef, 3ndirections(RT), Nᵣ))

# Copy the rotation data `R` into `rotors`, as quaternions.
function store_rotors!(rotors::AbstractVector{Quaternion{T}}, R::AbstractVector) where {T}
    @inbounds for i ∈ eachindex(rotors, R)
        rotors[i] = as_quaternion(T, R[i])
    end
    rotors
end

# The values of a vector of quaternions or angles that hold derivatives, without copying
# them.  This is what a lifting calculator gives the calculator of its rotors' values.
# - `T` is the element type of the values.
# - `V` is the type of the vector of quaternions or angles.
struct LiftedValues{T, V<:AbstractVector} <: AbstractVector{T}
    data::V
end
LiftedValues(q::AbstractVector{Quaternion{T}}) where {T} =
    LiftedValues{Quaternion{value_type(T)}, typeof(q)}(q)
LiftedValues(θ::AbstractVector{T}) where {T<:Real} = LiftedValues{value_type(T), typeof(θ)}(θ)
Base.size(v::LiftedValues) = size(v.data)
Base.IndexStyle(::Type{<:LiftedValues}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(v::LiftedValues{<:Quaternion}, i::Int) =
    rotor_value(v.data[i])
Base.@propagate_inbounds Base.getindex(v::LiftedValues{<:Real}, i::Int) = real_value(v.data[i])
