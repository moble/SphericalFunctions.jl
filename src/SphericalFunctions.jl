module SphericalFunctions

using TestItems: @testitem, @testsnippet
using FastTransforms: ifft, irfft
using Quaternionic: Quaternionic, from_spherical_coordinates
using StaticArrays: @SVector
using SpecialFunctions
using DoubleFloats
using FixedSizeArrays: FixedSizeArrayDefault, FixedSizeVectorDefault, FixedSizeArray,
    FixedSizeVector


const MachineFloat = Union{Float16, Float32, Float64}

# Utilities (kept top-level; code lives in `src/utilities/`)
include("utilities/utils.jl")

include("utilities/pixelizations.jl")
export golden_ratio_spiral_pixels, golden_ratio_spiral_rotors
export sorted_rings, sorted_ring_pixels, sorted_ring_rotors
export fejer1_rings, fejer2_rings, clenshaw_curtis_rings

include("utilities/complex_powers.jl")
export complex_powers, complex_powers!, ComplexPowers

include("utilities/weights.jl")
export fejer1, fejer2, clenshaw_curtis

include("wigner/wigner.jl")
export AbstractWignerMatrix, WignerMatrix, WignerDMatrix, WignerdMatrix
export WignerCalculator, WignerDCalculator, WignerdCalculator, WignerHCalculator
export recurrence!
public ℓ, ℓₘᵢₙ, m′ₘₐₓ, m′ₘᵢₙ, mₘₐₓ, mₘᵢₙ
public ell, ellmin, mpmax, mpmin, mmax, mmin

# Legacy API (no legacy names exported from the top-level module)
include("deprecated/Deprecated.jl")
export Deprecated

end # module SphericalFunctions
