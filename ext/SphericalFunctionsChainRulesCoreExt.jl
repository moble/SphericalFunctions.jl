module SphericalFunctionsChainRulesCoreExt

# `frule` and `rrule` for `D_array`, `sYlm_array`, and `sYlm_matrix_array`, from the
# generators, as described in `src/derivatives/kernels.jl`.  These serve the tools that read
# ChainRules directly, such as Zygote and Diffractor, which cannot follow the mutation in a
# calculator, and so differentiate `D`, `sYlm`, and `sYlm_matrix` of rotors through these
# functions instead.
#
# A tangent of a rotor may arrive as any quaternion, as a structural `Tangent` of the
# `Rotor`, or as a vector of four components; it is read by its components, whatever its
# type.  The cotangent returned is a `Quaternion`, not a `Rotor`: a cotangent is not a unit
# quaternion, and a `Rotor` would be taken to be one by any method that reached it.

import SphericalFunctions: D_array, d_array, D_series, sYlm_array, sYlm_matrix_array,
    wigner_arrays_with_derivative_values, wigner_arrays_pushforward, wigner_arrays_pullback,
    harmonic_array_pushforward, harmonic_array_pullback!, derivatives_from_left,
    rotor_generator, rotor_generators, rotor_cotangent, rotor_cotangents,
    rotation_angle_gradient, rotation_angle_cotangent, angle_generators, angle_cotangent,
    floattype, RotorLike, IntegerHalf, leading_view
using Quaternionic: AbstractQuaternion, Quaternion
import ChainRulesCore
using ChainRulesCore: AbstractZero, AbstractThunk, NoTangent, ZeroTangent, Tangent, unthunk,
    backing, RuleConfig, HasReverseMode, rrule_via_ad

# The four components of a tangent of a rotor, or `nothing` for a zero tangent.
tangent_components(::AbstractZero) = nothing
tangent_components(Ṙ::AbstractThunk) = tangent_components(unthunk(Ṙ))
tangent_components(Ṙ::AbstractQuaternion) = (Ṙ[1], Ṙ[2], Ṙ[3], Ṙ[4])
tangent_components(Ṙ::AbstractVector) = (Ṙ[1], Ṙ[2], Ṙ[3], Ṙ[4])
tangent_components(Ṙ::Tuple) = (Ṙ[1], Ṙ[2], Ṙ[3], Ṙ[4])
# A structural tangent with no fields is a zero tangent.
function tangent_components(Ṙ::Tangent{<:AbstractQuaternion})
    haskey(backing(Ṙ), :components) ? tangent_components(backing(Ṙ).components) : nothing
end
function tangent_components(Ṙ::Tangent)  # an `SVector`'s
    haskey(backing(Ṙ), :data) ? tangent_components(backing(Ṙ).data) : nothing
end

# The tangents of a vector of rotors, each as its four components, with a zero tangent as
# four zeros; or `nothing` when every one of them is zero.
function vector_tangent_components(Ṙ⃗, R⃗::AbstractVector)
    Ṙ⃗ = unthunk(Ṙ⃗)
    Ṙ⃗ isa AbstractZero && return nothing
    ṙ = [tangent_components(Ṙ) for Ṙ ∈ Ṙ⃗]
    all(isnothing, ṙ) && return nothing
    z = zero(real(eltype(eltype(R⃗))))
    [x === nothing ? (z, z, z, z) : x for x ∈ ṙ]
end

# The cotangents of the blocks of `D_array`, with the zero ones as `nothing`.
cotangent_blocks(Ā::AbstractVector) =
    [unthunk(Āᵢ) isa Union{AbstractZero, Nothing} ? nothing : unthunk(Āᵢ) for Āᵢ ∈ Ā]


## 𝔇

function ChainRulesCore.frule(
    (_, Ṙ, _, _, _, _, _), ::typeof(D_array),
    R::RotorLike, ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT
) where {IT<:IntegerHalf}
    blocks, values, calc = wigner_arrays_with_derivative_values(
        Complex{floattype(R)}, R, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ
    )
    ṙ = tangent_components(Ṙ)
    ṙ === nothing && return (blocks, ZeroTangent())
    v = rotor_generator(derivatives_from_left(calc), R, ṙ)
    (blocks, wigner_arrays_pushforward(calc, values, reshape([v[1], v[2], v[3]], 3, 1)))
end

function ChainRulesCore.rrule(
    ::typeof(D_array), R::RotorLike, ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT
) where {IT<:IntegerHalf}
    blocks, values, calc = wigner_arrays_with_derivative_values(
        Complex{floattype(R)}, R, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ
    )
    function D_array_rrule_pullback(ΔA)
        Ā = unthunk(ΔA)
        R̄ = if Ā isa AbstractZero
            ZeroTangent()
        else
            Ḡ = wigner_arrays_pullback(calc, values, cotangent_blocks(Ā))
            Quaternion(rotor_cotangent(derivatives_from_left(calc), R, (Ḡ[1], Ḡ[2], Ḡ[3]))...)
        end
        (NoTangent(), R̄, NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent())
    end
    (blocks, D_array_rrule_pullback)
end


## d
#
# The derivatives of `d` are taken with respect to the angle β of its rotor data (see
# `rotation_angle`), whose tangent and cotangent are formed here from those of the rotor
# data by the gradient of β, so that a phase or a rotor is differentiated through its angle.

# The tangent β̇ of the angle of `x` from a tangent `ẋ` of `x`, or `nothing` for a zero one.
angle_tangent(x::Real, ẋ) = (ẋ = unthunk(ẋ); ẋ isa AbstractZero ? nothing : ẋ)
angle_tangent(z::Complex, ż) =
    (ż = unthunk(ż); ż isa AbstractZero ? nothing : real(conj(rotation_angle_gradient(z)) * ż))
function angle_tangent(R::RotorLike, Ṙ)
    ṙ = tangent_components(Ṙ)
    ṙ === nothing && return nothing
    sum(rotation_angle_gradient(R) .* ṙ)
end
# The cotangent of `x` from that of its angle.
data_cotangent(x::Union{Real, Complex}, β̄) = rotation_angle_cotangent(x, β̄)
data_cotangent(R::RotorLike, β̄) = Quaternion(rotation_angle_cotangent(R, β̄)...)

function ChainRulesCore.frule(
    (_, ẋ, _, _, _, _, _), ::typeof(d_array),
    x::Union{Real, Complex, RotorLike}, ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT, mₘₐₓ::IT, mₘᵢₙ::IT
) where {IT<:IntegerHalf}
    blocks, values, calc = wigner_arrays_with_derivative_values(
        floattype(x), x, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ
    )
    β̇ = angle_tangent(x, ẋ)
    β̇ === nothing && return (blocks, ZeroTangent())
    (blocks, wigner_arrays_pushforward(calc, values, angle_generators([β̇])))
end

function ChainRulesCore.rrule(
    ::typeof(d_array), x::Union{Real, Complex, RotorLike}, ℓₘₐₓ::IT, m′ₘₐₓ::IT, m′ₘᵢₙ::IT,
    mₘₐₓ::IT, mₘᵢₙ::IT
) where {IT<:IntegerHalf}
    blocks, values, calc = wigner_arrays_with_derivative_values(
        floattype(x), x, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ
    )
    function d_array_rrule_pullback(ΔA)
        Ā = unthunk(ΔA)
        x̄ = if Ā isa AbstractZero
            ZeroTangent()
        else
            Ḡ = wigner_arrays_pullback(calc, values, cotangent_blocks(Ā))
            data_cotangent(x, angle_cotangent(Ḡ, 1))
        end
        (NoTangent(), x̄, NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent())
    end
    (blocks, d_array_rrule_pullback)
end


# `D_series` labels the blocks of `D_array` without copying them, and its pullback takes the
# cotangent of the series apart into those of the blocks' matrices.  The cotangent may come
# as a structural `Tangent` of the series, or as a `NamedTuple` of its fields, and its
# blocks likewise; any part of it may be a zero.
field_cotangent(Δ, name::Symbol) = ZeroTangent()
field_cotangent(Δ::Tangent, name::Symbol) =
    haskey(backing(Δ), name) ? unthunk(getproperty(Δ, name)) : ZeroTangent()
field_cotangent(Δ::NamedTuple, name::Symbol) = haskey(Δ, name) ? unthunk(Δ[name]) : ZeroTangent()
# The storage of each block is `vec` of its matrix, so the cotangent of that storage is
# given the matrix's shape.
shaped_like(Δ, b) = Δ
shaped_like(Δ::AbstractArray, b::AbstractArray) = reshape(Δ, size(b))
block_cotangents(Δblocks, blocks) = fill(ZeroTangent(), length(blocks))
block_cotangents(Δblocks::AbstractVector, blocks) =
    [shaped_like(field_cotangent(unthunk(Δb), :parent), b) for (Δb, b) ∈ zip(Δblocks, blocks)]

function ChainRulesCore.rrule(::typeof(D_series), blocks::AbstractVector, ℓₘₐₓ, limits...)
    S = D_series(blocks, ℓₘₐₓ, limits...)
    function D_series_pullback(ΔS)
        Δblocks = field_cotangent(unthunk(ΔS), :blocks)
        B̄ = Δblocks isa AbstractZero ? ZeroTangent() : block_cotangents(Δblocks, blocks)
        (NoTangent(), B̄, NoTangent(), map(_ -> NoTangent(), limits)...)
    end
    (S, D_series_pullback)
end


## The harmonics

function ChainRulesCore.frule(
    (_, Ṙ, _, _, _), ::typeof(sYlm_array), R::RotorLike, ℓₘₐₓ::IT, s, ℓₘᵢₙ::IT
) where {IT<:IntegerHalf}
    Y = sYlm_array(R, ℓₘₐₓ, s, ℓₘᵢₙ)
    ṙ = tangent_components(Ṙ)
    ṙ === nothing && return (Y, ZeroTangent())
    G = rotor_generators(true, [R], [ṙ])
    (Y, harmonic_array_pushforward((y, ẏ) -> only(ẏ), Y, false, ℓₘᵢₙ, ℓₘₐₓ, G, Val(1)))
end

function ChainRulesCore.rrule(
    ::typeof(sYlm_array), R::RotorLike, ℓₘₐₓ::IT, s, ℓₘᵢₙ::IT
) where {IT<:IntegerHalf}
    Y = sYlm_array(R, ℓₘₐₓ, s, ℓₘᵢₙ)
    function sYlm_array_rrule_pullback(ΔY)
        Ȳ = unthunk(ΔY)
        R̄ = if Ȳ isa AbstractZero
            ZeroTangent()
        else
            Ḡ = harmonic_array_pullback!(zeros(real(eltype(Y)), 3, 1), Y, Ȳ, false, ℓₘᵢₙ, ℓₘₐₓ)
            Quaternion(only(rotor_cotangents(true, [R], Ḡ))...)
        end
        (NoTangent(), R̄, NoTangent(), NoTangent(), NoTangent())
    end
    (Y, sYlm_array_rrule_pullback)
end

function ChainRulesCore.frule(
    (_, Ṙ⃗, _, _, _), ::typeof(sYlm_matrix_array),
    R⃗::AbstractVector{<:RotorLike}, ℓₘₐₓ::IT, s, ℓₘᵢₙ::IT
) where {IT<:IntegerHalf}
    Y = sYlm_matrix_array(R⃗, ℓₘₐₓ, s, ℓₘᵢₙ)
    ṙ = vector_tangent_components(Ṙ⃗, R⃗)
    ṙ === nothing && return (Y, ZeroTangent())
    G = rotor_generators(true, R⃗, ṙ)
    (Y, harmonic_array_pushforward((y, ẏ) -> only(ẏ), Y, true, ℓₘᵢₙ, ℓₘₐₓ, G, Val(1)))
end

function ChainRulesCore.rrule(
    ::typeof(sYlm_matrix_array), R⃗::AbstractVector{<:RotorLike}, ℓₘₐₓ::IT, s, ℓₘᵢₙ::IT
) where {IT<:IntegerHalf}
    Y = sYlm_matrix_array(R⃗, ℓₘₐₓ, s, ℓₘᵢₙ)
    function sYlm_matrix_array_rrule_pullback(ΔY)
        Ȳ = unthunk(ΔY)
        R̄⃗ = if Ȳ isa AbstractZero
            ZeroTangent()
        else
            Ḡ = harmonic_array_pullback!(
                zeros(real(eltype(Y)), 3, length(R⃗)), Y, Ȳ, true, ℓₘᵢₙ, ℓₘₐₓ
            )
            [Quaternion(R̄...) for R̄ ∈ rotor_cotangents(true, R⃗, Ḡ)]
        end
        (NoTangent(), R̄⃗, NoTangent(), NoTangent(), NoTangent())
    end
    (Y, sYlm_matrix_array_rrule_pullback)
end


## The storage of a block

# The elements of a block are a view of the leading entries of its storage, which
# `leading_view` builds: it is `view(p, Base.OneTo(n))`, constructed by `Base.unsafe_view`
# (see `src/containers/blocks.jl`).  Zygote cannot differentiate that constructor, so `leading_view` is differentiated as the
# `view` that it is, by the tool's own rule for `view`.
function ChainRulesCore.rrule(
    config::RuleConfig{>:HasReverseMode}, ::typeof(leading_view), p::AbstractVector, n::Int
)
    x, view_pullback = rrule_via_ad(config, view, p, Base.OneTo(n))
    leading_view_pullback(Δ) = (NoTangent(), view_pullback(Δ)[2], NoTangent())
    (x, leading_view_pullback)
end

end # module SphericalFunctionsChainRulesCoreExt
