"""
    SSHTDirect(s, ℓₘₐₓ; decomposition=LinearAlgebra.qr, T=Float64, Rθϕ=golden_ratio_spiral_rotors(s, ℓₘₐₓ, T), inplace=true)

Construct an ``s``-SHT object that uses the "Direct" method; see [`ₛ𝐘`](@ref) for details about the
method and optional arguments.  Also see [`SSHT`](@ref) for general information about how to use
these objects.

By default, this uses precisely optimal sampling — meaning that the number of points on which the
function is evaluated, represented by `Rθϕ`, is *equal to* the number of modes.  However, it is
equally possible to evaluate on *more* points than there are modes.  This can be useful, for
example, when processing multiple fields with different spin weights; the function could be
evaluated on points appropriate for the lowest value of ``|s|``, and therefore could also be used to
solve for fields of all other spin weights.

Note that in-place operation is possible for this type when the length of the input `Rθϕ` is equal
to the number of modes given `s` and `ℓₘₐₓ` — and is the default behavior when possible.  See
[`SSHT`](@ref) for description of in-place operation.

This method is typically better than other current implementations for ``ℓₘₐₓ ≲ 24``, both in terms
of speed and accuracy.  However, this advantage quickly falls away.  A warning will be issued if
`ℓₘₐₓ` is greater than about 64, because this method is not likely to be the most efficient or most
accurate choice.

"""
struct SSHTDirect{T<:Real, Inplace, Tdecomp} <: SSHT{T}
    """Spin weight"""
    s::Integer

    """Highest ℓ value present in the data"""
    ℓₘₐₓ::Integer

    """Rotors at which to evalue the ``s``-SH"""
    Rθϕ::Vector{Rotor{T}}

    """Spin-weighted spherical harmonic values"""
    ₛ𝐘::Matrix{Complex{T}}

    """Decomposed ₛ𝐘 matrix used in inversion

    LU decomposition (`LinearAlgebra.lu`) is used by default, but any function that decomposes ₛ𝐘
    into something capable of solving the linear problem may be passed to the constructor.  In
    particular, using QR decomposition (`LinearAlgebra.qr`) will typically be around 3 times slower,
    but have an error level roughly 10 times lower.  Also note that QR decomposition does not
    currently scale effectively with multiple threads.
    """
    ₛ𝐘decomposition::Tdecomp
end

function SSHTDirect(
    s, ℓₘₐₓ;
    decomposition=LinearAlgebra.lu,
    T=Float64, Rθϕ=golden_ratio_spiral_rotors(s, ℓₘₐₓ, T),
    inplace=inplaceable(s, ℓₘₐₓ, Rθϕ)
)
    if ((ℓₘₐₓ+1)^2-s^2)^2 > 65^4
        @warn """
        The "Direct" method for s-SHT is only recommended for fairly small ℓ values (or comparably large s values).
        Using it with ℓₘₐₓ=$ℓₘₐₓ and s=$s will be slow due to large memory requirements, and may be inaccurate.
        You will likely benefit from trying other methods for these parameters.
        """
    end
    if (message = check_blas_threads()) != ""
        @warn """$message
        Computations with SSHTDirect can benefit greatly from using all available threads if many
        functions are to be transformed.
        """
    end
    let ₛ𝐘 = ₛ𝐘(s, ℓₘₐₓ, T, Rθϕ)
        ₛ𝐘decomp = decomposition(ₛ𝐘)
        SSHTDirect{T, inplace, typeof(ₛ𝐘decomp)}(s, ℓₘₐₓ, Rθϕ, ₛ𝐘, ₛ𝐘decomp)
    end
end

function pixels(𝒯::SSHTDirect)
    to_spherical_coordinates.(rotors(𝒯))
end

function rotors(𝒯::SSHTDirect)
    𝒯.Rθϕ
end

function Base.:*(𝒯::SSHTDirect, f̃)
    𝒯.ₛ𝐘 * f̃
end

function LinearAlgebra.mul!(f, 𝒯::SSHTDirect, f̃)
    mul!(f, 𝒯.ₛ𝐘, f̃)
end

function Base.:\(𝒯::SSHTDirect, f)
    𝒯.ₛ𝐘decomposition \ f
end

function Base.:\(𝒯::SSHTDirect{T, true}, ff̃) where {T}
    ldiv!(𝒯.ₛ𝐘decomposition, ff̃)
end

function LinearAlgebra.ldiv!(f̃, 𝒯::SSHTDirect, f)
    ldiv!(f̃, 𝒯.ₛ𝐘decomposition, f)
end

function LinearAlgebra.ldiv!(𝒯::SSHTDirect, ff̃)
    ldiv!(𝒯.ₛ𝐘decomposition, ff̃)
end
