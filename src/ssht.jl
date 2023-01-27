"""Supertype of storage for spin-spherical-harmonic transforms"""
abstract type SSHT{T<:Real} end


@doc raw"""
    SSHT(s, ℓₘₐₓ; [method="Direct"], [T=Float64], [kwargs...])

Construct an `SSHT` object to transform between spin-weighted spherical-harmonic mode weights and
function values — performing an ``s``-SHT.

This object behaves similarly to an `AbstractFFTs.Plan` object — specfically in the ability to use
the semantics of algebra to perform transforms.  For example, if the function values are stored as
a vector `f`, the mode weights as `f̃`, and the `SSHT` as `𝒯`, then we can compute the function
values from the mode weights as

    f = 𝒯 * f̃

or solve for the mode weights from the function values as

    f̃ = 𝒯 \ f

The first dimensions of `f̃` must index the mode weights (as usual, for `ℓ∈abs(m):ℓₘₐₓ` and `m∈-ℓ:ℓ`)
and the first index of `f` must index the locations at which the function is evaluated.  Any
following dimensions will be broadcast over.  Note that certain types will broadcast using Julia
threads, while others will broadcast using BLAS threads.  The relevant number of threads must be set
appropriately.

Certain `SSHT` types (currently, only `Minimal` and `Direct`) also have an option to *always* act in
place — meaning that they simply re-use the input storage, even when used in an expression like
`𝒯 \ f`.  The option must be passed as the `inplace` argument to the constructors, and part of the
type of the resulting object.  Regardless of the value of that option, for those types where the
option exists, it is also possible to use `mul!` and `ldiv!` from the `LinearAlgebra` package to
force operation in place.

"""
function SSHT(s, ℓₘₐₓ; method="Direct", kwargs...)
    if method == "Direct"
        return SSHTDirect(s, ℓₘₐₓ; kwargs...)
    elseif method == "Minimal"
        return SSHTMinimal(s, ℓₘₐₓ; kwargs...)
    elseif method == "RS"
        return SSHTRS(s, ℓₘₐₓ; kwargs...)
    else
        error("""Unrecognized s-SHT method "$method".""")
    end
end


"""
    pixels(𝒯)

Return the spherical coordinates (θ, ϕ) on which the spin-weighted spherical harmonics are
evaluated.  See also [`rotors`](@ref), which provides the actual `Rotor`s on which they are
evaluated.
"""
function pixels(𝒯::SSHT)
    error("Not implemented for abstract supertype")
end

"""
    rotors(𝒯)

Return the `Rotor`s on which the spin-weighted spherical harmonics are evaluated.  See also
[`pixels`](@ref), which provides the corresponding spherical coordinates.
"""
function rotors(𝒯::SSHT)
    error("Not implemented for abstract supertype")
end


function Base.show(io::IO, ssht::SSHT)
    print(io, typeof(ssht), "($(ssht.s), $(ssht.ℓₘₐₓ))")
end


inplaceable(s, ℓₘₐₓ, Rθϕ) = (ℓₘₐₓ + 1)^2 - s^2 == length(Rθϕ)


function check_threads()
    cores = num_physical_cores()
    threads = nthreads()
    if threads < cores
        """
        You are using $threads threads on a machine with $cores physical cores.  While this
        will leave cores free to work on other tasks, you may wish to increase the
        number of threads available to Julia by restarting with the `--threads` flag.
        """
    else
        ""
    end
end


function check_blas_threads()
    cores = num_physical_cores()
    blas_threads = LinearAlgebra.BLAS.get_num_threads()
    if blas_threads < cores
        """
        You are using $blas_threads BLAS threads on a machine with $cores physical cores.  While
        this will leave cores free to work on other tasks, you may wish to increase the
        number of threads with something like

            LinearAlgebra.BLAS.set_num_threads($cores)
        """
    else
        ""
    end
end


include("ssht/direct.jl")
include("ssht/minimal.jl")
include("ssht/rs.jl")
