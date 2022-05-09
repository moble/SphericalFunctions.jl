using OffsetArrays

@testset "Associated Legendre Functions" begin
    nmax = 1_000
    min_length = (nmax * (nmax + 3)) ÷ 2 + 1

    for T in [BigFloat, Float64, Float32]
        ϵ = eps(T)
        P̄ = zeros(T, min_length)
        recursion_coefficients = ALFRecursionCoefficients(nmax, T)

        for β in βrange(T)
            expiβ = cis(β)
            ALFcompute!(P̄, expiβ, nmax, recursion_coefficients)
            P̄′ = ALFcompute(expiβ, nmax)
            @test eltype(P̄′) === eltype(P̄)

            @test_throws ArgumentError ALFcompute!(P̄, expiβ, nmax+1, recursion_coefficients)

            offset = -1
            for n in 0:nmax
                offset -= n
                P̄ₙ = OffsetVector(P̄, offset)
                # Test sum of squares
                total = sum(P̄ₙ[m]^2 for m in 0:n)
                @test total ≈ T(2n+1) rtol=2n*ϵ
                # Test computation without recursion_coefficients
                P̄ₙ′ = OffsetVector(P̄′, offset)
                @test P̄ₙ′[0:n] ≈ P̄ₙ[0:n] rtol=2n*ϵ
            end
        end
    end
end
