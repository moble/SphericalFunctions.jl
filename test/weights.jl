@testset "weights" begin
    @testset "$T" for T in [BigFloat, Float64, Float32]
        ϵ = 50eps(T)

        ϑ(k, n) = k * (π / T(n))  # k ∈ 0:n
        x(n, N) = cospi(n / T(N))  # = cos(ϑ(n, N))

        b(j, n) = j==n/2 ? 1 : 2
        c(k, n) = k%n==0 ? 1 : 2

        wᶠ¹(k, n) = (2/T(n)) * (1 - 2sum(cos(j*ϑ(2k+1, n))/T(4j^2-1) for j ∈ 1:(n÷2)))  # Eq. (2.3a) of Waldvogel; k ∈ 0:n-1
        wᶠ²(k, n) = (4/T(n)) * sin(ϑ(k, n)) * sum(sin((2j-1)*ϑ(k, n))/T(2j-1) for j ∈ 1:(n÷2))  # Eq. (2.3b) of Waldvogel; k ∈ 0:n
        wᶜᶜ(k, n) = (c(k, n)/T(n)) * (1 - sum(b(j, n) * cos(2j*ϑ(k, n))/T(4j^2-1) for j ∈ 1:(n÷2)))  # Eq. (4) of Waldvogel; k ∈ 0:n

        for n in [3, 10, 11, 12, 13, 170, 171, 1070, 1071]
            @test size(fejer1(n)) == (n,)
            @test size(fejer1(n, T)) == (n,)
            @test eltype(fejer1(n)) === Float64
            @test eltype(fejer1(n, T)) === T
            @test fejer1(n) ≈ fejer1(n, T) rtol=2max(eps(Float64), eps(T))
            # v1, v2 = fejer1(n, T), wᶠ¹.(0:n, n)
            # if ≉(v1, v2, rtol=ϵ, atol=ϵ)
            #     println("atol: ", maximum(abs, v1 .- v2) / ϵ)
            #     println("rtol: ", maximum(abs.(v1 .- v2) ./ abs.(v1)) / ϵ)
            #     println()
            # end
            @test fejer1(n, T) ≈ wᶠ¹.(0:n-1, n) rtol=ϵ atol=ϵ

            @test size(fejer2(n)) == (n,)
            @test size(fejer2(n, T)) == (n,)
            @test eltype(fejer2(n)) === Float64
            @test eltype(fejer2(n, T)) === T
            @test fejer2(n) ≈ fejer2(n, T) rtol=2max(eps(Float64), eps(T))
            # v1, v2 = fejer2(n, T), wᶠ².(1:n-1, n)
            # if ≉(v1, v2, rtol=ϵ, atol=ϵ)
            #     println("atol: ", maximum(abs, v1 .- v2) / ϵ)
            #     println("rtol: ", maximum(abs.(v1 .- v2) ./ abs.(v1)) / ϵ)
            #     println()
            # end
            @test fejer2(n, T) ≈ wᶠ².(1:n, n+1) rtol=ϵ atol=ϵ

            @test size(clenshaw_curtis(n)) == (n,)
            @test size(clenshaw_curtis(n, T)) == (n,)
            @test eltype(clenshaw_curtis(n)) === Float64
            @test eltype(clenshaw_curtis(n, T)) === T
            @test clenshaw_curtis(n) ≈ clenshaw_curtis(n, T) rtol=2max(eps(Float64), eps(T))
            # v1, v2 = clenshaw_curtis(n, T), wᶜᶜ.(0:n, n)
            # if ≉(v1, v2, rtol=ϵ, atol=ϵ)
            #     println("atol: ", maximum(abs, v1 .- v2) / ϵ)
            #     println("rtol: ", maximum(abs.(v1 .- v2) ./ abs.(v1)) / ϵ)
            #     println()
            # end
            @test clenshaw_curtis(n, T) ≈ wᶜᶜ.(0:n-1, n-1) rtol=ϵ atol=ϵ
        end
    end

    @testset "FastTransforms" begin
        ϵ = eps()

        for N in [3, 10, 11, 12, 13, 170, 171, 1070, 1071]
            μ1 = FastTransforms.chebyshevmoments1(Float64, N)
            μ2 = FastTransforms.chebyshevmoments2(Float64, N)
            w1 = fejerweights1(μ1)
            w2 = fejerweights2(μ2)
            wc = clenshawcurtisweights(μ1)

            @test w1 ≈ fejer1(N) rtol=ϵ atol=ϵ
            @test w2 ≈ fejer2(N) rtol=ϵ atol=ϵ
            @test wc ≈ clenshaw_curtis(N) rtol=ϵ atol=ϵ
        end
    end
end
