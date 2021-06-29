module Spherical

using FastTransforms

export complex_powers, complex_powers!
export WignerMatrixCalculator, H!, d!, D!
export ALFRecursionCoefficients, ALFrecurse!, ALFcompute!
export fejer1, fejer2, clenshaw_curtis


include("utils.jl")
include("complex_powers.jl")
include("indexing.jl")
include("wigner_matrices/evaluate.jl")
include("associated_legendre/calculator.jl")
include("weights.jl")


using .WignerMatrices
using .AssociatedLegendreFunction: ALFRecursionCoefficients, ALFrecurse!, ALFcompute!


end # module
