module SphericalFunctions

using TestItems: @testitem
using FFTW: FFTW, ifft, irfft, plan_fft!, plan_bfft!
# GenericFFT is loaded for its methods of the FFT functions above, which serve the element
# types that FFTW does not, such as `Float16`, `Double64` and `BigFloat`.
import GenericFFT
using LinearAlgebra: LinearAlgebra, mul!, ldiv!
using Quaternionic: Quaternionic, AbstractQuaternion, Quaternion, Rotor, QuatVec,
    from_spherical_coordinates
using StaticArrays: @SVector
import SpecialFunctions
import Serialization
using LinearAlgebra: Diagonal, Bidiagonal, Tridiagonal
using FixedSizeArrays: FixedSizeVectorDefault, FixedSizeVector
import Base: @propagate_inbounds
import PrecompileTools


# Base.IEEEFloat is not public, so we just define our own
const IEEEFloat = Union{Float16, Float32, Float64}

include("utilities/utils.jl")
include("utilities/lifting.jl")

include("utilities/half_odd_integer.jl")
include("utilities/index_methods.jl")

include("utilities/pixelizations.jl")
export golden_ratio_spiral_pixels, golden_ratio_spiral_rotors
export leja_pixels, leja_rotors
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

include("derivatives.jl")

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

# Names that are part of the documented interface but are not exported: accessors whose names
# are too generic to export, storage types and tools that most users never name, and the ASCII
# spellings of names written in Unicode.  `public` is a keyword only from Julia 1.11 on; the
# package supports Julia 1.10, where this is simply skipped.
VERSION ≥ v"1.11.0-DEV.469" && eval(Meta.parse(
    "public ℓ, ℓₘᵢₙ, ℓₘₐₓ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ, sₘₐₓ, sₘᵢₙ, spins, Nᵣ, "
    * "ell, ell_min, ell_max, mp_max, mp_min, m_max, m_min, s_max, s_min, Nr, "
    * "ishalfinteger, isbatched, "
    * "AbstractModeContainer, DifferentialOperator, Δspin, Deltaspin, "
    * "L2, Lplus, Lminus, R2, Rplus, Rminus, eth, ethbar, "
    * "HarmonicCalculator, sλlmCalculator, sλlm, sλlm!, sλlm_matrix, "
    * "slambdalmCalculator, slambdalm, slambdalm!, slambdalm_matrix, set_beta!, set_theta!, "
    * "HalfOddInteger, IntegerHalf, IndexType, IndexRange, IndexOrRange, @index_methods, "
    * "nmodes, npixels, HWedge, wedge_value, rotor_basetype, nrotors, floattype, sqrtbinomial, "
    * "driscoll_healy_pixels, driscoll_healy_rotors, mcewen_wiaux_pixels, mcewen_wiaux_rotors, "
    * "minimal_rings, map2salm_plan"
))

include("precompile.jl")

end # module SphericalFunctions
