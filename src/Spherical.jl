module Spherical

using FastTransforms, LinearAlgebra

export complex_powers, complex_powers!
export theta_phi, phi_theta
export WignerMatrixCalculator, H!, d!, D!, Y!
export ALFRecursionCoefficients, ALFrecurse!, ALFcompute!
export fejer1, fejer2, clenshaw_curtis
export map2salm, map2salm!, plan_map2salm

const MachineFloat = Union{Float16, Float32, Float64}


include("utils.jl")

include("complex_powers.jl")

include("indexing.jl")

include("wigner_matrices/evaluate.jl")

include("associated_legendre/calculator.jl")
using .AssociatedLegendreFunction: ALFRecursionCoefficients, ALFrecurse!, ALFcompute!

include("weights.jl")

include("map2salm.jl")


end # module
