@testset verbose=true "wigner_matrices" begin

    # Some helpful definitions
    include("test_utilities.jl")
    using SphericalFunctions: WignerHsize, WignerDsize, WignerHindex, WignerDindex, abd

    # Now the actual tests themselves
    include("wigner_matrices/small_d.jl")
    include("wigner_matrices/H.jl")
    include("wigner_matrices/big_D.jl")

end
