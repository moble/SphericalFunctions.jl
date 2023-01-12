module SphericalFunctions

using FastTransforms, LinearAlgebra, ProgressMeter, Quaternionic, OffsetArrays, AbstractFFTs
import LoopVectorization: @turbo
import Hwloc: num_physical_cores
import Base.Threads: @threads, nthreads, threadid

const MachineFloat = Union{Float16, Float32, Float64}


include("utils.jl")

include("pixelizations.jl")
export golden_ratio_spiral, sorted_rings
export fejer1_rings, fejer2_rings, clenshaw_curtis_rings

include("complex_powers.jl")
export complex_powers, complex_powers!

include("indexing.jl")
export Ysize, Yrange, Yindex, deduce_limits, theta_phi, phi_theta
export WignerHsize, WignerHindex, _WignerHindex, WignerHrange
export WignerDsize, WignerDindex, WignerDrange

include("iterators.jl")
export Diterator, diterator, Yiterator

include("associated_legendre.jl")
export ALFRecursionCoefficients, ALFrecurse!, ALFcompute!, ALFcompute

include("Hrecursion.jl")
export H!, H_recursion_coefficients

include("evaluate.jl")
export d!, d, D!, Y!
export dprep, dstorage, Dprep, Dstorage, Yprep, Ystorage
export ₛ𝐘

include("weights.jl")
export fejer1, fejer2, clenshaw_curtis

include("ssht.jl")
export SSHT, pixels

include("map2salm.jl")
export map2salm, map2salm!, plan_map2salm

#include("rotate.jl")
#export rotate!


end # module
