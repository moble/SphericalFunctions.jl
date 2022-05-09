using OffsetArrays

@testset "Associated Legendre Functions" begin
    nmax = 1_000
    min_length = (nmax * (nmax + 3)) ÷ 2 + 1

    for T in [BigFloat, Float64, Float32]
        ϵ = eps(T)
        P̄ = fill(T(NaN), min_length)
        P̄′ = Vector{T}(undef, min_length)
        recursion_coefficients = ALFRecursionCoefficients(nmax, T)

        for β in βrange(T)
            expiβ = cis(β)

            # Make sure we're checking the size of P̄ against nmax
            @test_throws ArgumentError ALFcompute!(P̄, expiβ, nmax+1, recursion_coefficients)

            ALFcompute!(P̄, expiβ, nmax, recursion_coefficients)
            ALFcompute!(P̄′, expiβ, nmax)
            P̄″ = ALFcompute(expiβ, nmax)
            @test eltype(P̄″) === eltype(P̄)

            offset = -1
            for n in 0:nmax
                offset -= n
                P̄ₙ = OffsetVector(P̄, offset)
                # Test sum of squares
                total = sum(P̄ₙ[m]^2 for m in 0:n)
                @test total ≈ T(2n+1) rtol=2n*ϵ
                # Test in-place computation without recursion_coefficients
                P̄ₙ′ = OffsetVector(P̄′, offset)
                @test P̄ₙ′[0:n] ≈ P̄ₙ[0:n] rtol=2n*ϵ
                # Test allocating computation without recursion_coefficients
                P̄ₙ″ = OffsetVector(P̄″, offset)
                @test P̄ₙ″[0:n] ≈ P̄ₙ[0:n] rtol=2n*ϵ
            end
        end
    end
end
