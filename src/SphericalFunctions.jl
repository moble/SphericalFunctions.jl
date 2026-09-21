module SphericalFunctions

using TestItems: @testitem, @testsnippet
using FastTransforms: FastTransforms, FFTW, ifft, irfft, plan_fft!, plan_bfft!, fftshift!, ifftshift!
using LinearAlgebra: LinearAlgebra, mul!, ldiv!
using Base.Threads: @threads
using Quaternionic: Quaternionic, AbstractQuaternion, Rotor, QuatVec, from_spherical_coordinates
using StaticArrays: @SVector
using SpecialFunctions
using LinearAlgebra: Diagonal, Bidiagonal, Tridiagonal
using FixedSizeArrays: FixedSizeVectorDefault, FixedSizeVector
using OffsetArrays: OffsetArray, OffsetVector, OffsetMatrix
import Base: @propagate_inbounds


# Base.IEEEFloat is not public, so we just define our own
const IEEEFloat = Union{Float16, Float32, Float64}

include("utilities/utils.jl")

include("utilities/half_odd_integer.jl")

include("utilities/pixelizations.jl")
export golden_ratio_spiral_pixels, golden_ratio_spiral_rotors
export sorted_rings, sorted_ring_pixels, sorted_ring_rotors
export fejer1_rings, fejer2_rings, clenshaw_curtis_rings

include("utilities/complex_powers.jl")
export complex_powers, complex_powers!, ComplexPowers

include("utilities/weights.jl")
export fejer1, fejer2, clenshaw_curtis

include("utilities/operators.jl")
export L², Lz, L₊, L₋, Lx, Ly, R², Rz, R₊, R₋, ð, ð̄

include("wigner/wigner.jl")
export AbstractWignerMatrix, WignerMatrix, WignerDMatrix, WignerdMatrix
export WignerMatrixBatch, DegreeBlock, DegreeBlockBatch, WignerSeries
export SpinMatrix, SpinMatrixBatch
export WignerCalculator, DCalculator, dCalculator, HCalculator
export recurrence!, D, d

include("mode_weights/indexing.jl")
export Ysize, Yindex, Yrange

include("mode_weights/containers.jl")
export HarmonicValues

include("sYlm/sYlm.jl")
export sYlmCalculator, sYlm, sYlm!, sYlm_matrix, Ylm, YlmCalculator

include("set_rotor_data.jl")
export set_R!, set_β!, set_θ!

include("iteration.jl")

include("mode_weights/mode_weights.jl")
export ModeWeights, modes, spin

include("array_view.jl")
export array_view, relabel

include("mode_weights/operations.jl")

include("ssht/ssht.jl")
export SSHT, SSHTMatrix, SSHTRS, SSHTMinimal, pixels, rotors, map2salm, salm2map

# Names that are part of the documented interface but are not exported, either because they
# are accessors whose names are too generic to export, or because they are storage types that
# most users never name.  `public` is a keyword only from Julia 1.11 on; the package supports
# Julia 1.10, where this is simply skipped.
VERSION ≥ v"1.11.0-DEV.469" && eval(Meta.parse(
    "public ℓ, ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, sₘₐₓ, sₘᵢₙ, spins, Nᵣ, "
    * "AbstractModeContainer, DifferentialOperator, Δspin, "
    * "HarmonicCalculator, sλlmCalculator, sλlm, sλlm!, sλlm_matrix, "
    * "ell, ellmin, ellmax, mpmax, mpmin, mmax, mmin, smax, smin, Nr, ishalfinteger, isbatched, "
    * "HalfOddInteger, IntegerHalf, "
    * "nmodes, npixels, HWedge, HAxis, rotor_basetype, nrotors, floattype, "
    * "driscoll_healy_pixels, driscoll_healy_rotors, mcewen_wiaux_pixels, mcewen_wiaux_rotors"
))

end # module SphericalFunctions
