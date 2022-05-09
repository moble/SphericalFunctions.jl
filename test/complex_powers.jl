@testset "complex_powers" begin
    complex_powers_comparison(z, m, T=Float64) = (
        complex_powers(Complex{T}(z), m),
        z.^collect(0:m),
        eps(T)*m
    )

    for T in [Float64, Float32, Float16]
        nozpowers = Vector{Complex{T}}(undef, 0)
        fudge = one(T) + 2 * sqrt(eps(T))
        for k in 0:25
            z = cis(k*big(π)/10)
            for m in [0, 1, 2, 3, 4, 1_000]
                mine, theirs, ϵ = complex_powers_comparison(z, m, T)
                @test mine ≈ theirs rtol=2ϵ
                @test_throws DomainError complex_powers(Complex{T}(z*fudge), m)
                @test_throws DomainError complex_powers(Complex{T}(z/fudge), m)
                inplace = zeros(Complex{T}, size(mine))
                complex_powers!(inplace, Complex{T}(z))
                @test array_equal(mine, inplace)
            end
            @test length(complex_powers!(nozpowers, Complex{T}(z))) == 0
        end
    end

end
