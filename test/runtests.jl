using Spherical
using Test

@testset "complex_powers" begin
    complex_powers_comparison(z, m, T=Float64) = (
        complex_powers(Complex{T}(z), m),
        z.^collect(0:m),
        eps(T)*m
    )

    for T in [Float64, Float32, Float16]
        for k in 0:25
            z = exp(k*big(π)im/10)
            for m in [0, 1, 2, 3, 4, 1_000]
                mine, theirs, ϵ = complex_powers_comparison(z, m, T)
                @test mine ≈ theirs rtol=2ϵ
            end
        end
    end
end
