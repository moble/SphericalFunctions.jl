module SphericalFunctions

using TestItems: @testitem, @testsnippet
using FastTransforms: ifft, irfft
using Quaternionic: from_spherical_coordinates
using StaticArrays: @SVector
using SpecialFunctions, DoubleFloats

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

# New public API: promote redesign to the top-level interface
include("redesign/SphericalFunctions.jl")

# Legacy API (no legacy names exported from the top-level module)
include("deprecated/Deprecated.jl")
export Deprecated

end # module SphericalFunctions
