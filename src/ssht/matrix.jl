"""
    SSHTMatrix(s, ℓₘₐₓ, [T=Float64]; Rθϕ=leja_rotors(s, ℓₘₐₓ, T), decomposition, inplace)

Construct an ``s``-SHT object that uses the "Matrix" method.  The dense matrix of harmonics
[`sYlm_matrix`](@ref) evaluated at the rotors `Rθϕ` is stored, so that synthesis is a
matrix-vector product and analysis is a linear solve using the given `decomposition` of that
matrix.  Also see [`SSHT`](@ref) for general information about how to use these objects.

By default, this uses precisely optimal sampling — meaning that *the number of points* on
which the function is evaluated (the length of `Rθϕ`), *is equal to the number of modes*
``(ℓₘₐₓ+1)²``.  However, it is equally possible to evaluate on *more* points than there are
modes, in which case the analysis is a least-squares solve.  This can be useful, for
example, when processing multiple fields with different spin weights; the function could be
evaluated on points appropriate for the lowest value of ``|s|``, and therefore could also be
used to solve for fields of all other spin weights.  The rotors must be `Rotor{T}`s, in the
type `T` that the transform works in.

With precisely optimal sampling, the accuracy is limited by the conditioning of the points.
The default points, those of [`leja_rotors`](@ref), are chosen to keep the matrix well
conditioned: with the default LU decomposition, a round trip of random mode weights loses
about 2.5 digits (in the norm of the error relative to that of the weights) at ℓₘₐₓ = 32, and
3 at ℓₘₐₓ = 64.  Choosing them costs about 2–3 times as much as the decomposition itself.
Points that merely spread evenly over the sphere are not enough: on the golden-ratio spiral of
[`golden_ratio_spiral_rotors`](@ref), a round trip loses about 4 digits at ℓₘₐₓ = 16, 5 at
32, and all of them at 64.  The constructor measures the error of one round trip, and warns
when fewer than half the digits of `T` would survive.

LU decomposition (`LinearAlgebra.lu`) is used by default for a square matrix; any function
that decomposes the matrix into something capable of solving the linear problem may be
passed as `decomposition`.  In particular, QR decomposition (`LinearAlgebra.qr`) will
typically be about twice as slow, with errors a few times smaller; it is also the default
when there are more points than modes, where an LU decomposition, which cannot solve a
least-squares problem, is rejected.

In-place operation is possible for this type when the length of the input `Rθϕ` is equal to
the number of modes given `s` and `ℓₘₐₓ` — and is the default behavior when possible.  It
applies to analysis only: with `inplace=true`, `𝒯 \\ f` overwrites the storage of `f` with
the mode weights, and returns them as a `ModeWeights` wrapping it, while `𝒯 * f̃` always
allocates a new array for the function values, since a matrix product cannot write over its
own input.  The two-argument `LinearAlgebra.ldiv!(𝒯, f)` analyzes in place whatever the
option, including for more points than modes, when the mode weights are left in the first
`nmodes(𝒯)` entries of `f`; there is no two-argument `mul!`.  Real function values are
analyzed with `inplace=false`, since in place the complex mode weights would have to be
written into real storage.  See [`SSHT`](@ref) for a description of in-place operation.

The object holds no workspace: its transforms only read the matrix and its decomposition,
through BLAS and LAPACK, so one object may be used by any number of tasks at the same time,
and `copy(𝒯)` is `𝒯` itself.

This method is fastest for ``ℓₘₐₓ ≲ 24``, where its round-trip errors are about
``10^{-14}``, somewhat larger than those of the `"RS"` method.  However, this advantage
quickly falls away.  A warning is issued if `ℓₘₐₓ` is greater than about 64, because this
method is not likely to be the most efficient or most accurate choice.
"""
struct SSHTMatrix{T<:Real, Inplace, Tdecomp, IT<:IntegerHalf} <: SSHT{T}
    s::IT
    ℓₘₐₓ::IT
    Rθϕ::Vector{Rotor{T}}
    Y::Matrix{Complex{T}}  # Spin-weighted spherical harmonic values, [pixel, mode]
    Ydecomposition::Tdecomp
end

# The keyword defaults are `nothing` and are resolved in the body, after the indices and the
# type have been checked: the default points of `leja_rotors` would otherwise be computed,
# at some cost, for a call that is then refused, and the defaults of `decomposition` and
# `inplace` depend on the points.
@index_methods function SSHTMatrix(
    s::IT, ℓₘₐₓ::IT, ::Type{TT}=Float64;
    Rθϕ=nothing, decomposition=nothing, inplace=nothing
) where {IT<:IndexType, TT}
    check_transform_type(TT)
    check_band_limit(s, ℓₘₐₓ)
    n = Ysize(abs(s), ℓₘₐₓ)
    Rθϕ = Rθϕ === nothing ? leja_rotors(s, ℓₘₐₓ, TT) : Rθϕ
    check_sample_rotors(TT, Rθϕ)
    if n^2 > 65^4
        @warn """
        The "Matrix" method for s-SHT is only recommended for fairly small ℓ values (or comparably large s values).
        Using it with ℓₘₐₓ=$ℓₘₐₓ and s=$s will be slow due to large memory requirements, and may be inaccurate.
        You will likely benefit from trying other methods for these parameters.
        """
    end
    if length(Rθϕ) < n
        throw(DimensionMismatch(
            "There are $(length(Rθϕ)) sample points but $n modes; the analysis would be "
            * "underdetermined."
        ))
    end
    square = length(Rθϕ) == n
    inplace = inplace === nothing ? square : inplace
    check_inplace_option(inplace)
    if inplace && !square
        throw(ArgumentError(
            "In-place operation requires exactly as many sample points as modes, but there "
            * "are $(length(Rθϕ)) points and $n modes."
        ))
    end
    decomposition = decomposition === nothing ?
        (square ? LinearAlgebra.lu : LinearAlgebra.qr) : decomposition
    Rs = Vector{Rotor{TT}}(Rθϕ)  # a copy, never a conversion: `check_sample_rotors` saw to that
    Y = sYlm_matrix(Rs, ℓₘₐₓ, s)
    Ydecomp = decomposition(Y)
    if !square && Ydecomp isa LinearAlgebra.LU
        throw(ArgumentError(
            "An LU decomposition cannot solve the least-squares problem of $(length(Rs)) "
            * "points and $n modes; use a QR decomposition, `LinearAlgebra.qr`, which is the "
            * "default."
        ))
    end
    𝒯 = SSHTMatrix{TT, inplace, typeof(Ydecomp), IT}(s, ℓₘₐₓ, Rs, Y, Ydecomp)
    warn_if_inaccurate(𝒯, "Matrix")
    𝒯
end

pixels(𝒯::SSHTMatrix) = Quaternionic.to_spherical_coordinates.(𝒯.Rθϕ)
rotors(𝒯::SSHTMatrix) = copy(𝒯.Rθϕ)
npixels(𝒯::SSHTMatrix) = length(𝒯.Rθϕ)

# The transforms only read the fields, and `rotors` returns a copy of the rotors, so nothing
# that a copy would make independent is ever modified: the copy is the object itself.
Base.copy(𝒯::SSHTMatrix) = 𝒯

# Dense linear algebra uses a vector or a matrix, so any dimensions beyond the first are
# folded into the columns (and unfolded again on the way out).  `reshape` of an `Array`
# shares its memory, so the in-place methods stay in place.
@inline flatten_trailing(a) = ndims(a) ≤ 2 ? a : reshape(a, size(a, 1), :)
@inline function unflatten_trailing(flat, dims)
    length(dims) ≤ 2 ? flat : reshape(flat, size(flat, 1), Base.tail(dims)...)
end

# LAPACK solves in place only in storage whose columns are contiguous, so any other storage
# — a strided view, or a row of a matrix — is copied into an array, solved there, and copied
# back.
function solve_in_place!(F, x)
    if x isa StridedArray && stride(x, 1) == 1
        ldiv!(F, x)
    else
        y = Array(x)
        ldiv!(F, y)
        copyto!(x, y)
    end
    x
end

function Base.:*(𝒯::SSHTMatrix, f̃::SSHTData)
    d = synthesis_modes(𝒯, f̃)
    unflatten_trailing(𝒯.Y * flatten_trailing(d), size(d))
end
# A matrix product cannot write its output over its own input — BLAS reads the input while
# writing the output, so that `mul!(x, 𝒯, x)` would return zeros — and an input that shares
# memory with the output is therefore copied first.  The same holds for the solve below.
function LinearAlgebra.mul!(f, 𝒯::SSHTMatrix, f̃)
    check_modes(𝒯, f̃)
    check_pixels(𝒯, f)
    check_trailing(f, f̃)
    check_complex_output(f, "f")
    d, F = array_view(f̃), array_view(f)
    mul!(flatten_trailing(F), 𝒯.Y, flatten_trailing(Base.mightalias(F, d) ? copy(d) : d))
    f
end

function Base.:\(𝒯::SSHTMatrix, f::SSHTData)
    check_pixels(𝒯, f)
    d = array_view(f)
    f̃ = 𝒯.Ydecomposition \ flatten_trailing(d)
    ndims(d) == 1 ? ModeWeights(f̃, 𝒯.s, abs(𝒯.s), 𝒯.ℓₘₐₓ) : unflatten_trailing(f̃, size(d))
end
function Base.:\(𝒯::SSHTMatrix{T, true}, ff̃::SSHTData) where {T}
    check_pixels(𝒯, ff̃)
    check_complex_storage(ff̃, "𝒯 \\ f")
    solve_in_place!(𝒯.Ydecomposition, flatten_trailing(array_view(ff̃)))
    in_place_modes(𝒯, ff̃)
end
function LinearAlgebra.ldiv!(f̃, 𝒯::SSHTMatrix, f)
    f̃ = analysis_output(𝒯, f̃, f)
    check_modes(𝒯, f̃)
    check_pixels(𝒯, f)
    check_trailing(f, f̃)
    check_complex_output(f̃, "f̃")
    d, F = array_view(f̃), array_view(f)
    ldiv!(
        flatten_trailing(d), 𝒯.Ydecomposition,
        flatten_trailing(Base.mightalias(d, F) ? copy(F) : F)
    )
    f̃
end
function LinearAlgebra.ldiv!(𝒯::SSHTMatrix, ff̃::SSHTData)
    check_pixels(𝒯, ff̃)
    check_complex_storage(ff̃, "ldiv!(𝒯, f)", "𝒯 \\ f")
    solve_in_place!(𝒯.Ydecomposition, flatten_trailing(array_view(ff̃)))
    solution_modes(𝒯, ff̃)
end
# An in-place transform is square, so its solution fills the whole argument.
solution_modes(𝒯::SSHTMatrix{T, true}, ff̃) where {T} = in_place_modes(𝒯, ff̃)
