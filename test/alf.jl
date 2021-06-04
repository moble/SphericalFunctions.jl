using OffsetArrays

@testset "Associated Legendre Functions" begin
    nmax = 1_000
    min_length = (nmax * (nmax + 3)) ÷ 2 + 1

    for T in [BigFloat, Float64, Float32]
        ϵ = eps(T)
        P̄ = zeros(T, min_length)
        recursion_coefficients = ALFRecursionCoefficients(nmax, T)

        θs = collect(zero(T):one(T)/10:T(π)/2)
        
        for (i, θ) in enumerate(θs)
            expiβ = exp(θ * 1im)
            ALFcompute!(P̄, expiβ, nmax, recursion_coefficients)
            offset = -1
            for n in 0:nmax
                offset -= n
                P̄ₙ = OffsetVector(P̄, offset)
                total = sum(P̄ₙ[m]^2 for m in 0:n)
                @test total≈T(2n+1) rtol=2ϵ*n
            end
        end
    end
end
