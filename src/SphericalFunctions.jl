module SphericalFunctions

using FastTransforms, LinearAlgebra, ProgressMeter, Quaternionic
import Base.Threads: @threads, nthreads

export complex_powers, complex_powers!
export Ysize, Yrange, Yindex, deduce_limits, theta_phi, phi_theta
export ALFRecursionCoefficients, ALFrecurse!, ALFcompute!, ALFcompute
export fejer1, fejer2, clenshaw_curtis
export map2salm, map2salm!, plan_map2salm

const MachineFloat = Union{Float16, Float32, Float64}


include("utils.jl")

include("complex_powers.jl")

include("indexing.jl")

include("wigner_matrices/evaluate.jl")

include("associated_legendre.jl")

include("weights.jl")

include("map2salm.jl")

include("calculators.jl")

include("Hrecursion.jl")


end # module
