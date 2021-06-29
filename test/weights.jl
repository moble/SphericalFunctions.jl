@testset "weights" begin
    @testset "$T" for T in [Float64, Float32, BigFloat]
        ϑ(k, n) = k * π / T(n)  # k ∈ 0:n
        x(n, N) = cospi(n / T(N))  # = cos(ϑ(n, N))

        b(j, n) = j==n/2 ? 1 : 2
        c(k, n) = k%n==0 ? 1 : 2

        wᶠ¹(k, n) = (2/T(n)) * (1 - 2sum(cos(j*ϑ(2k+1, n))/T(4j^2-1) for j ∈ 1:(n÷2)))  # Eq. (2.3a) of Waldvogel; k ∈ 0:n-1
        wᶠ²(k, n) = (4/T(n)) * sin(ϑ(k, n)) * sum(sin((2j-1)*ϑ(k, n))/T(2j-1) for j ∈ 1:(n÷2))  # Eq. (2.3b) of Waldvogel; k ∈ 0:n
        wᶜᶜ(k, n) = (c(k, n)/T(n)) * (1 - sum(b(j, n) * cos(2j*ϑ(k, n))/T(4j^2-1) for j ∈ 1:(n÷2)))  # Eq. (4) of Waldvogel; k ∈ 0:n

        for n in [3, 10, 11, 12, 13, 170, 171, 1070, 1071]
            @test size(fejer1(n)) == (n+1,)
            @test size(fejer1(n, T)) == (n+1,)
            @test eltype(fejer1(n)) === complex(Float64)
            @test eltype(fejer1(n, T)) === complex(T)
            @test fejer1(n) ≈ fejer1(n, T) rtol=2max(eps(Float64), eps(T))
            @test fejer1(n, T) ≈ wᶠ¹.(0:n, n)) < 2eps(T)

            @test size(fejer2(n)) == (n-1,)
            @test size(fejer2(n, T)) == (n-1,)
            @test eltype(fejer2(n)) === complex(Float64)
            @test eltype(fejer2(n, T)) === complex(T)
            @test fejer2(n) ≈ fejer2(n, T) rtol=2max(eps(Float64), eps(T))
            @test fejer2(n, T) ≈ wᶠ².(1:n-1, n)) < 2eps(T)

            @test size(curtis_clenshaw(n)) == (n+1,)
            @test size(curtis_clenshaw(n, T)) == (n+1,)
            @test eltype(curtis_clenshaw(n)) === complex(Float64)
            @test eltype(curtis_clenshaw(n, T)) === complex(T)
            @test curtis_clenshaw(n) ≈ curtis_clenshaw(n, T) rtol=2max(eps(Float64), eps(T))
            @test curtis_clenshaw(n, T) ≈ wᶜᶜ.(0:n, n)) < 2eps(T)
        end
    end
end
