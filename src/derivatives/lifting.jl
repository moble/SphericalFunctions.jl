### Rotor data, and rotors that hold derivatives
#
# A calculator whose rotor data — rotors, angles, or phases — hold derivatives, dual numbers
# say, runs the recurrence on the values of those data, and lifts each block of values into
# a block that holds the derivatives, which the angular-momentum operators give in terms of
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
# tuples of four components, and those of an angle, as a tuple of numbers; and the function
# that assembles an element of type `Complex{T}` or `T` from a value and the tuple of its
# derivatives in those directions.  An extension defines all of these for its number type.
real_value(x::Real) = x
rotor_value(q::Quaternion{T}) where {T} =
    Quaternion{value_type(T)}(real_value(q[1]), real_value(q[2]), real_value(q[3]), real_value(q[4]))
function ndirections end
function rotor_tangents end
function angle_tangents end
function lift_combine end

# The data of a calculator that lifts the blocks of another: that calculator, of the values
# of its rotors, and the generators of its rotors' derivatives, with the generator of rotor
# iᵣ in direction d in `G[3d-2:3d, iᵣ]` (see `set_generators!`).  The rotors of a calculator
# of `d` or of ₛλₗₘ are rotations about the y axis, whose generators have only a y
# component, so that `G[d, iᵣ]` holds it alone (see `AngleGenerators`).
# - `C` is the type of the calculator of the values.
# - `M` is the type of the matrix of generators.
struct Lift{C, M<:AbstractMatrix}
    inner::C
    G::M
end
allocate_lift(::Type{RT}, ::Type{NT}, inner, Nᵣ::Int) where {RT, NT} =
    Lift(inner, Matrix{value_type(RT)}(undef, (NT <: Complex ? 3 : 1) * ndirections(RT), Nᵣ))

# Copy the rotation data `R` into `rotors`, as quaternions.
function store_rotors!(rotors::AbstractVector{Quaternion{T}}, R::AbstractVector) where {T}
    @inbounds for i ∈ eachindex(rotors, R)
        rotors[i] = as_quaternion(T, R[i])
    end
    rotors
end

# Copy the angles of the rotor data `R`, one datum or a vector of them, into `angles` (see
# `rotation_angle`), each by `store_angle!` from the datum's components.  A calculator of
# floats has no use for the angles, whose tangents and cotangents alone are what the rules
# of Enzyme and Mooncake read; so for floats `store_angle!` does not compute the angle, and
# those tools have rules for it that put the angle's tangent, or the cotangents of the
# components, into their shadows (see `rotation_angle_gradient`).  This keeps the angles'
# cost, chiefly that of `atan`, off the calculators of floats.
function store_angles!(angles::AbstractVector, R::AbstractVector)
    @inbounds for i ∈ eachindex(angles, R)
        store_angle!(angles, i, rotor_data_components(R[i])...)
    end
    angles
end
function store_angles!(angles::AbstractVector, R)
    store_angle!(angles, 1, rotor_data_components(R)...)
    angles
end
Base.@propagate_inbounds function store_angle!(angles::AbstractVector{T}, i::Int, x::Real...) where {T}
    angles[i] = convert(T, rotation_angle(rotor_data(x...)))
    nothing
end
# For floats, the first component of the datum is stored in the angle's place, and is never
# read.  That store costs a fraction of a nanosecond, but a method that did nothing at all
# would have its call deleted by the compiler in the code that Enzyme and Mooncake
# differentiate, and their rules for it would never run.  (A store of the angle's own value
# back into it is no better: Enzyme then takes the datum to be constant.)
Base.@propagate_inbounds function store_angle!(
    angles::AbstractVector{T}, i::Int, x::Real...
) where {T<:Base.IEEEFloat}
    angles[i] = x[1]
    nothing
end

# The values of a vector of quaternions, angles, or phases that hold derivatives, without
# copying them.  This is what a lifting calculator gives the calculator of its rotors'
# values.
# - `T` is the element type of the values.
# - `V` is the type of the vector of quaternions, angles, or phases.
struct LiftedValues{T, V<:AbstractVector} <: AbstractVector{T}
    data::V
end
LiftedValues(q::AbstractVector{<:Union{Rotor{T}, Quaternion{T}}}) where {T} =
    LiftedValues{Quaternion{value_type(T)}, typeof(q)}(q)
LiftedValues(θ::AbstractVector{T}) where {T<:Real} = LiftedValues{value_type(T), typeof(θ)}(θ)
LiftedValues(z::AbstractVector{Complex{T}}) where {T<:Real} =
    LiftedValues{Complex{value_type(T)}, typeof(z)}(z)
# Rotor data that is not a vector of one of those is passed on as it is, for the calculator
# of the values to refuse with its own message.
LiftedValues(R) = R
Base.size(v::LiftedValues) = size(v.data)
Base.IndexStyle(::Type{<:LiftedValues}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(v::LiftedValues{<:Quaternion}, i::Int) =
    rotor_value(v.data[i])
Base.@propagate_inbounds Base.getindex(v::LiftedValues{<:Real}, i::Int) = real_value(v.data[i])
Base.@propagate_inbounds function Base.getindex(v::LiftedValues{<:Complex}, i::Int)
    z = v.data[i]
    Complex(real_value(real(z)), real_value(imag(z)))
end
