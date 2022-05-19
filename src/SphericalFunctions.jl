module SphericalFunctions

using FastTransforms, LinearAlgebra, ProgressMeter, Quaternionic
import ChainRulesCore: rrule, frule
import Base.Threads: @threads, nthreads

const MachineFloat = Union{Float16, Float32, Float64}


include("utils.jl")

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
export dprep, dstorage, Dprep, Dstorage, Yprep, Ystorage#, ϵ

include("weights.jl")
export fejer1, fejer2, clenshaw_curtis

include("map2salm.jl")
export map2salm, map2salm!, plan_map2salm

#include("rotate.jl")
#export rotate!


end # module
