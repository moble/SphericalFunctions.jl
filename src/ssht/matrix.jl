"""
    SSHTMatrix(s, ℓₘₐₓ; decomposition=LinearAlgebra.lu, T=Float64, Rθϕ=golden_ratio_spiral_rotors(s, ℓₘₐₓ, T), inplace=…)

Construct an ``s``-SHT object that uses the "Matrix" method: the dense matrix of harmonics
[`sYlm_matrix`](@ref) evaluated at the rotors `Rθϕ` is stored, so that synthesis is a
matrix-vector product and analysis is a linear solve using the given `decomposition` of that
matrix.  Also see [`SSHT`](@ref) for general information about how to use these objects.

By default, this uses precisely optimal sampling — meaning that the number of points on which
the function is evaluated, represented by `Rθϕ`, is *equal to* the number of modes.  However,
it is equally possible to evaluate on *more* points than there are modes, in which case the
analysis is a least-squares solve.  This can be useful, for example, when processing multiple
fields with different spin weights; the function could be evaluated on points appropriate
for the lowest value of ``|s|``, and therefore could also be used to solve for fields of all
other spin weights.

LU decomposition (`LinearAlgebra.lu`) is used by default for a square matrix; any function
that decomposes the matrix into something capable of solving the linear problem may be passed
as `decomposition`.  In particular, QR decomposition (`LinearAlgebra.qr`) will typically be
around 3 times slower, but have an error level roughly 10 times lower; it is also the default
when there are more points than modes.

In-place operation is possible for this type when the length of the input `Rθϕ` is equal to
the number of modes given `s` and `ℓₘₐₓ` — and is the default behavior when possible.  See
[`SSHT`](@ref) for a description of in-place operation.

This method is typically better than other current implementations for ``ℓₘₐₓ ≲ 24``, both
in terms of speed and accuracy.  However, this advantage quickly falls away.  A warning is
issued if `ℓₘₐₓ` is greater than about 64, because this method is not likely to be the most
efficient or most accurate choice.
"""
struct SSHTMatrix{T<:Real, Inplace, Tdecomp, IT<:HalfInteger} <: SSHT{T}
    s::IT
    ℓₘₐₓ::IT
    Rθϕ::Vector{Rotor{T}}
    Y::Matrix{Complex{T}}  # Spin-weighted spherical harmonic values, [pixel, mode]
    Ydecomposition::Tdecomp
end

inplaceable(s, ℓₘₐₓ, Rθϕ) = Ysize(abs(s), ℓₘₐₓ) == length(Rθϕ)

# The public constructor is the boundary: it normalizes the two indices and re-dispatches
# to the worker, whose keyword defaults are then computed from indices of one kind.
function SSHTMatrix(s::IndexSpelling, ℓₘₐₓ::IndexSpelling; T::Type{TT}=Float64, kwargs...) where {TT}
    SSHTMatrix(transform_indices(s, ℓₘₐₓ)..., TT; kwargs...)
end
function SSHTMatrix(
    s::IT, ℓₘₐₓ::IT, ::Type{TT};
    Rθϕ=golden_ratio_spiral_rotors(s, ℓₘₐₓ, TT),
    decomposition=(inplaceable(s, ℓₘₐₓ, Rθϕ) ? LinearAlgebra.lu : LinearAlgebra.qr),
    inplace=inplaceable(s, ℓₘₐₓ, Rθϕ)
) where {IT<:HalfInteger, TT}
    if abs(s) > ℓₘₐₓ
        error("|s|=$(abs(s)) exceeds ℓₘₐₓ=$ℓₘₐₓ; there are no such modes.")
    end
    if Ysize(abs(s), ℓₘₐₓ)^2 > 65^4
        @warn """
        The "Matrix" method for s-SHT is only recommended for fairly small ℓ values (or comparably large s values).
        Using it with ℓₘₐₓ=$ℓₘₐₓ and s=$s will be slow due to large memory requirements, and may be inaccurate.
        You will likely benefit from trying other methods for these parameters.
        """
    end
    if length(Rθϕ) < Ysize(abs(s), ℓₘₐₓ)
        error(
            "There are $(length(Rθϕ)) sample points but $(Ysize(abs(s), ℓₘₐₓ)) modes; "
            * "the analysis would be underdetermined."
        )
    end
    if inplace && !inplaceable(s, ℓₘₐₓ, Rθϕ)
        error("In-place operation requires exactly as many sample points as modes.")
    end
    Rs = Vector{Rotor{TT}}(Rθϕ)
    Y = sYlm_matrix(Rs, ℓₘₐₓ, s)
    Ydecomp = decomposition(Y)
    SSHTMatrix{TT, inplace, typeof(Ydecomp), IT}(s, ℓₘₐₓ, Rs, Y, Ydecomp)
end

pixels(𝒯::SSHTMatrix) = Quaternionic.to_spherical_coordinates.(rotors(𝒯))
rotors(𝒯::SSHTMatrix) = 𝒯.Rθϕ
npixels(𝒯::SSHTMatrix) = length(𝒯.Rθϕ)

# Dense linear algebra wants a vector or a matrix, so any dimensions beyond the first are
# folded into the columns (and unfolded again on the way out).  `reshape` of an `Array`
# shares its memory, so the in-place methods stay in place.
@inline flatten_trailing(a) = ndims(a) ≤ 2 ? a : reshape(a, size(a, 1), :)
@inline function unflatten_trailing(flat, dims)
    length(dims) ≤ 2 ? flat : reshape(flat, size(flat, 1), Base.tail(dims)...)
end
function check_trailing(f, f̃)
    if size(f)[2:end] != size(f̃)[2:end]
        error("Trailing dimensions of f $(size(f)[2:end]) and f̃ $(size(f̃)[2:end]) differ.")
    end
end

function Base.:*(𝒯::SSHTMatrix, f̃)
    check_modes(𝒯, f̃)
    d = strided(f̃)
    unflatten_trailing(𝒯.Y * flatten_trailing(d), size(d))
end
function LinearAlgebra.mul!(f, 𝒯::SSHTMatrix, f̃)
    check_modes(𝒯, f̃)
    check_pixels(𝒯, f)
    check_trailing(f, f̃)
    mul!(flatten_trailing(f), 𝒯.Y, flatten_trailing(strided(f̃)))
    f
end

function Base.:\(𝒯::SSHTMatrix, f)
    check_pixels(𝒯, f)
    f̃ = 𝒯.Ydecomposition \ flatten_trailing(f)
    ndims(f) == 1 ? ModeWeights(f̃, 𝒯.s, abs(𝒯.s), 𝒯.ℓₘₐₓ) : unflatten_trailing(f̃, size(f))
end
function Base.:\(𝒯::SSHTMatrix{T, true}, ff̃) where {T}
    check_pixels(𝒯, ff̃)
    ldiv!(𝒯.Ydecomposition, flatten_trailing(strided(ff̃)))
    ff̃
end
function LinearAlgebra.ldiv!(f̃, 𝒯::SSHTMatrix, f)
    check_modes(𝒯, f̃)
    check_pixels(𝒯, f)
    check_trailing(f, f̃)
    ldiv!(flatten_trailing(strided(f̃)), 𝒯.Ydecomposition, flatten_trailing(f))
    f̃
end
function LinearAlgebra.ldiv!(𝒯::SSHTMatrix, ff̃)
    check_pixels(𝒯, ff̃)
    ldiv!(𝒯.Ydecomposition, flatten_trailing(strided(ff̃)))
    ff̃
end

