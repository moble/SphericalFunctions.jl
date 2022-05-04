@testset verbose=true "wigner_matrices" begin
    include("wigner_matrices/small_d.jl")
    include("wigner_matrices/H.jl")
    include("wigner_matrices/big_D.jl")
end
