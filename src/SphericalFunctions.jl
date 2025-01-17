module SphericalFunctions

using FastTransforms: FFTW, fft, fftshift!, ifft, ifftshift!, irfft,
                      plan_bfft!, plan_fft, plan_fft!
using LinearAlgebra: LinearAlgebra, Bidiagonal, Diagonal, convert, ldiv!, mul!
using OffsetArrays: OffsetArray, OffsetVector
using ProgressMeter: Progress, next!
using Quaternionic: Quaternionic, Rotor, from_spherical_coordinates,
                    to_euler_phases, to_spherical_coordinates
using StaticArrays: @SVector
using SpecialFunctions, DoubleFloats
using LoopVectorization: @turbo
using Base.Threads: @threads, nthreads
using TestItems: @testitem

const MachineFloat = Union{Float16, Float32, Float64}


include("utils.jl")

include("pixelizations.jl")
export golden_ratio_spiral_pixels, golden_ratio_spiral_rotors
export sorted_rings, sorted_ring_pixels, sorted_ring_rotors
export fejer1_rings, fejer2_rings, clenshaw_curtis_rings

include("complex_powers.jl")
export complex_powers, complex_powers!

include("indexing.jl")
export Ysize, Yrange, Yindex, deduce_limits, theta_phi, phi_theta
export WignerHsize, WignerHindex, _WignerHindex, WignerHrange
export WignerDsize, WignerDindex, WignerDrange

include("iterators.jl")
export D_iterator, d_iterator, sYlm_iterator, λ_iterator
# Legacy API:
export Diterator, diterator, Yiterator, λiterator

include("associated_legendre.jl")
export ALFRecursionCoefficients, ALFrecurse!, ALFcompute!, ALFcompute

include("Hrecursion.jl")
export H!, H_recursion_coefficients

include("evaluate.jl")
export d_prep, d_matrices, d_matrices!
export D_prep, D_matrices, D_matrices!
export sYlm_prep, sYlm_values, sYlm_values!
# Legacy API:
export d!, d, D!, Y!, dprep, Dprep, Yprep, ₛ𝐘

include("weights.jl")
export fejer1, fejer2, clenshaw_curtis

include("ssht.jl")
export SSHT, pixels, rotors

include("map2salm.jl")
export map2salm, map2salm!, plan_map2salm

include("operators.jl")
export L², Lz, L₊, L₋, R², Rz, R₊, R₋, ð, ð̄

#include("rotate.jl")
#export rotate!



end # module
