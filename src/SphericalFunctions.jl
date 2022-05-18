module SphericalFunctions

using FastTransforms, LinearAlgebra, ProgressMeter, Quaternionic
import Base.Threads: @threads, nthreads

export complex_powers, complex_powers!
export Ysize, Yrange, Yindex, deduce_limits, theta_phi, phi_theta
export WignerDsize, WignerHsize, WignerDindex, WignerHindex, _WignerHindex
export H!, H_recursion_coefficients
export d!, d, D!, Y!
export dprep, dstorage, Dprep, Dstorage, Yprep, Ystorage#, ϵ
export ALFRecursionCoefficients, ALFrecurse!, ALFcompute!, ALFcompute
export fejer1, fejer2, clenshaw_curtis
export map2salm, map2salm!, plan_map2salm

const MachineFloat = Union{Float16, Float32, Float64}


include("utils.jl")

include("complex_powers.jl")

include("indexing.jl")

include("Hrecursion.jl")

include("evaluate.jl")

include("associated_legendre.jl")

include("weights.jl")

include("map2salm.jl")


end # module
