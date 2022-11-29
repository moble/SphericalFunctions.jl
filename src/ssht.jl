"""Supertype of storage for spin-spherical-harmonic transforms"""
abstract type SSHT end

@doc raw"""
    SSHT(s, ℓₘₐₓ; [method="Direct"], [T=Float64], [kwargs...])

Construct an `SSHT` object to transform between spin-weighted spherical-harmonic mode weights and
function values — performing an ``s``-SHT.

This object behaves similarly to an `AbstractFFTs.Plan` object — specfically in the ability to use
the semantics of algebra to perform transforms.  For example, if the mode weights are stored as a
vector `𝐦`, the function values as `𝐟`, and the `SSHT` as `𝒯`, then we can compute the function
values from the mode weights as

    𝐟 = 𝒯 * 𝐦

or solve for the mode weights from the function values as

    𝐦 = 𝒯 \ 𝐟

The first dimensions of `𝐦` must index the mode weights (as usual, for `ℓ∈abs(m):ℓₘₐₓ` and `m∈-ℓ:ℓ`)
and the first index of `𝐟` must index the locations at which the function is evaluated.  Any
following dimensions will be broadcast over.  Note that certain types will broadcast using Julia
threads, while others will broadcast using BLAS threads.  The relevant number of threads must be set
appropriately.

Certain `SSHT` types also have an option to *always* act in place — meaning that they simply re-use
the input storage, even when used in an expression like `𝒯 \ 𝐟`.  Regardless of the value of that
option, for those types where the option exists, it is also possible to use `mul!` and `ldiv!` from
the `LinearAlgebra` package to force operation in place.

"""
function SSHT(s, ℓₘₐₓ; method="Direct", kwargs...)
    if method == "Direct"
        return SSHTDirect(s, ℓₘₐₓ; kwargs...)
    elseif method == "EKKM"
        return SSHTEKKM(s, ℓₘₐₓ; kwargs...)
    elseif method == "RS"
        return SSHTRS(s, ℓₘₐₓ; kwargs...)
    else
        @error """Unrecognized s-SHT method "$method"."""
    end
end


function Base.show(io::IO, ssht::SSHT)
    print(io, typeof(ssht), "($(ssht.s), $(ssht.ℓₘₐₓ))")
end


inplaceable(s, ℓₘₐₓ, Rθϕ) = (ℓₘₐₓ+1)^2 - s^2 == length(Rθϕ)


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
include("ssht/ekkm.jl")
include("ssht/eg.jl")
include("ssht/rs.jl")
